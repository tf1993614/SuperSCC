import re
import os 
from os.path import basename, splitext

import numpy as np
import pandas as pd
from uuid import uuid4, UUID
import textwrap
from dotenv import load_dotenv
from typing import List

from SuperSCC import list_files, jaccard_score


from langchain_qdrant import QdrantVectorStore
from langchain_qdrant.qdrant import QdrantVectorStoreError
from langchain.prompts import ChatPromptTemplate, PromptTemplate
from langchain_openai import ChatOpenAI
from langchain_community.document_loaders import PyMuPDFLoader
from langchain_community.document_loaders.csv_loader import CSVLoader
from langchain_community.cross_encoders import HuggingFaceCrossEncoder
from langchain.retrievers.document_compressors.cross_encoder_rerank import CrossEncoderReranker
from langchain.retrievers.document_compressors.cross_encoder import BaseCrossEncoder
from langchain.retrievers.contextual_compression import ContextualCompressionRetriever
from langchain.retrievers import BM25Retriever, EnsembleRetriever
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_huggingface import HuggingFaceEmbeddings
from langchain_core.output_parsers import StrOutputParser, PydanticOutputParser
from langchain_core.runnables import RunnablePassthrough
from langchain.globals import set_llm_cache
from langchain.cache import InMemoryCache
from qdrant_client.http.models import Filter, FieldCondition, MatchValue
from pydantic import BaseModel, Field
from qdrant_client import QdrantClient
from qdrant_client.http.models import Distance, VectorParams, VectorParamsDiff

import streamlit as st



class HighlightDocuments(BaseModel):
    id: List[str] = Field(..., description="List of id of docs used for answering the question")
    title: List[str] = Field(..., description="List of titles used to answers the question")
    source: List[str] = Field( ..., description="List of sources used to answers the question")
    segment: List[str] = Field( ..., description="List of direct segements from used documents that answers the question")


class ScoreDocuments(BaseModel):
    binary_score: str = Field(..., description = "Documents are relevant to the question. 'Yes' or 'No'")


class TranslateQuery(BaseModel):
     binary_score: str = Field(..., description = "Should do translation. 'Yes' or 'No'")
     translation: str = Field(..., description = "translated or not translated content")


class SimpleRAG:
    def __init__(self, file_path: str, file_type: str):
        self.file_path = file_path
        self.file_type = file_type
        self.recursive_search_cond = False if re.search(self.file_type, os.path.splitext(basename(self.file_path))[1]) else True
        self.file_pool = list()
        self.split_document = list()
        self.retrieve_docs = list()
        self.relevant_segments = list()
        self.answer = list()
        self.dense_retriver = list()
        self.gene_list = None
        self.doc_filter = None
        self.retrieve_docs_by_BM25 = None
        self.retrieve_docs_by_dense = None
        self.retrieve_docs_after_rerank = None
        self.vector_store = None
        self.rag_chain = None
        self.retriever = None
        self.hybrid_retriever = None
        self.prompt = None
        self.qdrant_client = None
        self.qdrant_collection_name = None
        self.qdrant_host = None 
        self.__qdrant_api_key = None
        self.__llm_api_key = None
        self.llm_api_base = None
        self.llm_model = None
        self.llm = None
        self.text_embedding_model = None
        self.model_kwargs = None
        self.encode_kwargs = None
        self.rerank_model = None
        self.query = None
        self.hybrid_search_cond = None

    def __repr__(self):
       info = f"""
        SimpleRAG object with {len(self.split_document)} documents x {len(np.concatenate(self.split_document))} text chunk

        Text embedding model: {self.text_embedding_model}
        LLM model: {self.llm_model}
        Rerank model: {self.rerank_model}
        Answer prompt: {self.prompt}
        Qdrant collection name: {self.qdrant_collection_name}
        Current query: {self.query}
        No.relevant segments for the current query: {len(self.relevant_segments.segment) if self.relevant_segments != None else None} 
       """
       return info

    def recursive_search(self, path, type = None):
        current_path = path 
        current_path_content = os.listdir(current_path)
        for i in current_path_content:
            if os.path.isdir(i):
                new_path = os.path.join(current_path, i)
                self.recursive_search(path = new_path)
            else:
                if re.search(self.file_type if type == None else type, i):
                    self.file_pool.append(os.path.join(current_path, i))
                else:
                    continue

    def data_loader(self, file_path, mode = "single", metadata_columns = None):
        if self.file_type == "pdf" or self.file_type[-1] == "pdf":
            loader = PyMuPDFLoader(file_path = file_path, mode = mode)
        elif self.file_type == "csv" or self.file_type[-1] == "csv":
            loader = CSVLoader(file_path = file_path, metadata_columns = metadata_columns)
        
        docs = loader.load()

        if self.file_type == "csv" or self.file_type[-1] == "csv":
            for i in docs:
                i.page_content = i.page_content.replace("Markers:", "").strip()
        return docs


    def text_split(self, docs, chunk_size = 1000, chunk_overlap = 200, length_function = len):
        text_splitter = RecursiveCharacterTextSplitter(
                            chunk_size = chunk_size,
                            chunk_overlap = chunk_overlap,
                            length_function = length_function
                        )
        texts = text_splitter.split_documents(docs)
        for doc in texts:
            doc.page_content = doc.page_content.replace("\t", " ")
        return texts

    def text_encode(self, 
                    text,
                    model_name,
                    location, 
                    host = None,
                    api_key = None,
                    url = None,
                    collection_name = "SimpleRAG",
                    model_kwargs = {"device": "cpu"}, 
                    encode_kwargs = {"normalize_embeddings": True},
                    vectors_config = VectorParams(size = 1024, distance=Distance.COSINE),
                    timeout = 1000):
    
  
        embedding_model = HuggingFaceEmbeddings(model_name = model_name, model_kwargs= model_kwargs, encode_kwargs=encode_kwargs)
        self.model_kwargs = model_kwargs
        self.encode_kwargs = encode_kwargs

        client = QdrantClient(location = location, host = host, api_key = api_key, timeout = timeout, url = url)
        self.qdrant_collection_name = collection_name
        self.text_embedding_model = model_name

        if self.qdrant_host == None and self.__qdrant_api_key == None:
            self.qdrant_host = host
            self.__qdrant_api_key = api_key

        if not client.collection_exists(collection_name = collection_name):
            client.create_collection(collection_name = collection_name, vectors_config = vectors_config)
        self.qdrant_client = client

        vector_store = QdrantVectorStore(client = client, collection_name = collection_name, embedding = embedding_model)
        ids = [str(uuid4()) for _ in range(len(text))]
        vector_store.add_documents(documents=text, ids = ids)
        self.vector_store = vector_store
        return vector_store
    
    
    def create_rag_chain(self, 
                         vector_store, 
                         model,
                         api_key,
                         base_url,
                         prompt = None,
                         search_kwargs = {"k": 10}
                         ):
        if self.retriever == None:
            retriever = vector_store.as_retriever(search_kwargs = search_kwargs)
            self.retriever = retriever
        
        self.search_kwargs = search_kwargs
        
        template = """
                 Answer the question with using following the context and relevant information : 
                  {context}

                   Question: {question}
                 """
        if prompt is None:
            prompt = ChatPromptTemplate.from_template(template)
        else:
            prompt = prompt
        self.prompt = prompt
        self.llm_model = model
        
        llm = ChatOpenAI(model = model, 
                 temperature = 0, 
                 api_key = api_key,
                 base_url= base_url)
        self.llm = llm
        
        if self.__llm_api_key is None and self.llm_api_base is None:
            self.__llm_api_key = api_key
            self.llm_api_base = base_url
        
        rag_chain = prompt | llm | StrOutputParser()
        

        self.rag_chain = rag_chain
        return rag_chain
    
    def update_rag_chain(self, model = None, api_key = None, base_url = None, prompt = None, update_similarity_search = False, search_kwargs = {"k": 10}):
        if update_similarity_search is False:
            
            if prompt == None:
                prompt = self.prompt
            else:
                prompt = prompt

            self.prompt = prompt

            if api_key == None and base_url == None:
                api_key = self.__llm_api_key
                base_url = self.llm_api_base
            else:
                assert api_key != None and base_url != None, "Please provide api key and base url simultaneously"
                api_key = api_key
                base_url = base_url

            self.api_key = api_key
            self.base_url = base_url
                
            if model == None:
                model = self.llm_model
            else:
                model = model

            self.llm_model = model

            llm = ChatOpenAI(
                model = model, 
                temperature = 0, 
                api_key = api_key,
                base_url= base_url
            )
            self.llm = llm
            
            # rag_chain = (
            #     {"context": self.retrieve_docs_after_rerank, "question": RunnablePassthrough()}
            #     | prompt
            #     | llm
            #     | StrOutputParser()
            # )
            self.rag_chain = prompt | llm | StrOutputParser()
        else:
            self.rag_chain = self.create_rag_chain(
                vector_store = self.vector_store,
                model = model,
                api_key = api_key,
                base_url = base_url,
                prompt = prompt,
                search_kwargs = search_kwargs
            )

    def add_documents(self, file_path, file_type, chunk_overlap = 200, chunk_size = 1024, recursive_search = True):
        if recursive_search:
            files = self.recursive_search(file_path, type = file_type)
        else:
            files = [file_path]
        self.file_path = [self.file_path, file_path]
        self.file_type = [self.file_type, file_type]
        
        for file in files:
            print(f"Start adding document {file}")
            docs = self.data_loader(file_path = file)
            if file_type == "pdf":
                split_text = self.text_split(docs = docs, chunk_overlap = chunk_overlap, chunk_size = chunk_size)
            else:
                split_text = docs 
            self.split_document.append(split_text)   
            ids = [str(uuid4()) for _ in range(len(split_text))]
            self.vector_store.add_documents(documents=split_text, ids = ids)
            print(f"Finish adding document {file}")
    
    def change_text_embedding(self, model_name, model_kwargs = None, encode_kwargs = None, vector_config = VectorParams(size = 1024, distance=Distance.COSINE)):

        if model_kwargs == None:
            model_kwargs = self.model_kwargs
        else:
            model_kwargs = model_kwargs

        if encode_kwargs == None:
            encode_kwargs = self.encode_kwargs
        else:
            encode_kwargs = encode_kwargs

        embedding_model = HuggingFaceEmbeddings(model_name = model_name, model_kwargs= model_kwargs, encode_kwargs = encode_kwargs)
        self.text_embedding_model = model_name
        try:
            self.vector_store = QdrantVectorStore(client = self.qdrant_client, collection_name = self.qdrant_collection_name, embedding = embedding_model)
        except QdrantVectorStoreError as e:
            print(e)
            print("The dimension of the new text embedding model is not consistent with the settings of initating the qdrant collecction.\nProvide the correct dimesnison of the new embedding model via 'vector_config' argument")
           
            print(f"Delet the old collection {self.qdrant_collection_name}")
            self.qdrant_client.delete_collection(collection_name = self.qdrant_collection_name)

            print(f"Build the new collection {self.qdrant_collection_name}")
            self.qdrant_client.create_collection(collection_name = self.qdrant_collection_name, vectors_config = vector_config)
            
            print("Add vectors")
            self.vector_store = QdrantVectorStore(client = self.qdrant_client, collection_name = self.qdrant_collection_name, embedding = embedding_model)
            if self.split_document != None:
                for idx, i in enumerate(self.split_document):
                    print(f"Processing Splitted document {idx}")
                    ids = [str(uuid4()) for _ in range(len(i))]
                    self.vector_store.add_documents(documents = i, ids = ids)
            self.retriever = self.vector_store.as_retriever(search_kwargs = self.search_kwargs)
            self.update_rag_chain()

    def format_docs(self, docs):
        return "\n".join(
            f"<doc{i + 1}-{doc.metadata['_id']}>:\ngene_list_in_retrieval_context:{doc.page_content}\nData:{doc.metadata['Data']}\nCluster:{doc.metadata['Cluster']}\nOrganization:{doc.metadata['Organization']}\ncell_type:{doc.metadata['cell_type']}\nlevel:{doc.metadata['level']}\n</doc{i + 1}>\n" for i, doc in enumerate(docs)
        )


    def run_rag(self,
                qdrant_location, 
                text_embedding_model,
                llm_model,
                llm_api_key,
                llm_base_url,
                qdrant_host = None,
                qdrant_api_key = None,
                qdrant_url = None,
                metadata_columns = None,
                chunk_size = 1000,
                chunk_overlap = 200,
                text_model_kwargs = {"device": "cpu"}, 
                text_encode_kwargs = {"normalize_embeddings": True},
                qdrant_collection_name = "SimpleRAG",
                qdrant_search_kwargs = {"k": 10},
                vectors_config = VectorParams(size = 1024, distance=Distance.COSINE)
                ):
        
        if self.recursive_search_cond:
            self.recursive_search(self.file_path)

        for file in self.file_pool:
            print(f"Processing with file {file}")
            
            docs = self.data_loader(file_path = file, metadata_columns = metadata_columns)
            print(f"Finishing loading file")
            
            # only do text splitting when file type is PDF
            if self.file_type == "pdf":
                split_text = self.text_split(docs = docs, chunk_overlap = chunk_overlap, chunk_size = chunk_size)
                print(f"Finishing document splitting")
                self.split_document.append(split_text)

            vector_store = self.text_encode(text = split_text if self.file_type == "pdf" else docs, # if file type is not "PDF", A Document object list without text splitting as the input
                                            host = qdrant_host, 
                                            url = qdrant_url,
                                            location = qdrant_location,
                                            api_key = qdrant_api_key, 
                                            model_name = text_embedding_model, 
                                            model_kwargs = text_model_kwargs, 
                                            encode_kwargs = text_encode_kwargs, 
                                            collection_name = qdrant_collection_name,
                                            vectors_config = vectors_config)
            
            print("Finishing text embedding")
            print("=" *20)
        
        print(f"Similarity search to find the top {qdrant_search_kwargs['k']} relevant text chunk")
        self.create_rag_chain(vector_store = vector_store, model = llm_model, api_key = llm_api_key, base_url = llm_base_url, search_kwargs = qdrant_search_kwargs)
        print("Finishing constructing rag chain")

    
    def score_documents(self, docs = None):
        parser = PydanticOutputParser(pydantic_object = ScoreDocuments)
        system = """
                 You are a scorer assessing the relevance of a retrieved document to a user question.
                 if the document contains keyword(s) or semantic meaning related to the user question,
                 score it as relevant. It does not need to be a stringent test. The goal is to filter out 
                 erroneous retrivals. Give a binary score 'Yes' or 'No' score to indicate whether the document
                 is relevant to the question.

                 Retrieve document: <docs>{document}</docs> \n\n User question: <question>{question}</question>

                 <format_instruction>
                 {format_instructions}
                 </format_instruction>
                 """

        score_prompt = PromptTemplate(
            template = system,
            input_variables = ["question", "document"],
            partial_variables = {"format_instructions": parser.get_format_instructions()}
        )

        
        scorer = score_prompt | self.llm | parser
        
        ls = list()
        if docs == None:
            docs = self.retrieve_docs_after_rerank
        else:
            docs = docs
        for idx, doc in enumerate(docs):
            res = scorer.invoke({"document": doc.page_content, "question": self.gene_list})
            if res.binary_score == "Yes":
                ls.append(idx)
        
        print(f"Only {len(ls)} members out of  {len(docs)} retrieved text chunks are correlated with the current query")
        self.retrieve_docs = [docs[i] for i in ls]


    def refine_query(self):
        print(f"{'='*10} Rewrite query {'='*10}")
        print(f"Original query: {self.query}")
        if(len(self.query) < 20):
            refine_query_template = "You are an AI assistant with reformulating user queries to improve retrieval in a RAG system. Given the original query, rewrite to be more specific, detailed, and likely to retrieve relevant information. (Important) only return rewritten query.  \n\nOrigiinal query: {query}\n\nRewritten query: "
            refine_prompt = PromptTemplate(template = refine_query_template, input_variables = ["query"] )
            refiner = refine_prompt | self.llm
            self.query = refiner.invoke(self.query).content
        else:
            refine_query_template = """
                                    You are an AI assistant with breaking down complex queries into simpler sub-queries for a RAG system.
                                    Given the originary query, decompose it into 2-4 simpler sub-queries that, when answered together, would 
                                    provide a comprehensive response to the original query. Only return sub-querires after decomposition. 

                                    Original_query: {query}

                                    example: What are the impact of climate change on the environment?

                                    Sub-queries:
                                    1. What are the impacts of climate change on biodiversity?
                                    2. How does climate change affect the oceans?
                                    3. What are the effects of climate change on agriculture?
                                    4. What are the impacts of climate change on human health?
                                    """
            refine_prompt = PromptTemplate(template = refine_query_template, input_variables = ["query"] )
            refiner = refine_prompt | self.llm
            self.query = refiner.invoke(self.query).content.replace("/n", "")
        print(f"Refine query: {self.query}")

    def translator(self, query = None):
        print(f"{'='*10} translate  query in English {'='*10}")
        parser = PydanticOutputParser(pydantic_object = TranslateQuery)
        system = """
                 You are a Chinese-English translation expert, translating the Chinese input by the user into English.
                 The user provides the content to be translated to the assistant. 
                 The assistant first determines whether the content provided by the user contains Chinese characters and gives a binary score 
                 'Yes' or 'No' score to indicate whether the content has Chinese characters. 
                 If with Chinese characters, it performs the translation, otherwise it does not perform the translation and returns the content provided by the user. 
                 When performing the translation, the following requirements must be met: the original text must be translated into a translation that meets the standards of faithfulness,
                 expressiveness and elegance. "Faithfulness" means being faithful to the content and intention of the original text; 
                 "expressiveness" means that the translation should be fluent, easy to understand and clearly expressed; 
                 "elegance" pursues the cultural aesthetics of the translation and the beauty of the language. 
                 The goal is to create a translation that is both faithful to the spirit of the original work and in line with the target language culture and the aesthetics of the readers.

                 content: {query}

                 <format_instruction>
                 {format_instructions}
                 </format_instruction>
                 """
        prompt = PromptTemplate(template=system, input_variables=["query"], partial_variables= {"format_instructions": parser.get_format_instructions()})
        translator = prompt | self.llm | parser
        if query != None:
            self.query = query
        res = translator.invoke({"query": self.query})
        print(f"Do translation: {res.binary_score}")
        self.query = res.translation.replace("\n", "")
        print(f"Current query: {self.query}")

    def get_all_ids(self):
        ids = list()
        scroll_response = self.qdrant_client.scroll(collection_name = self.qdrant_collection_name, limit = 100, offset = None, with_payload = ["id"], with_vectors = False)
        records = scroll_response[0]
        offset = scroll_response[1]
        ids.extend([records[i].id for i in range(len(records))])
        while offset is not None:
            scroll_response = self.qdrant_client.scroll(collection_name = self.qdrant_collection_name, limit = 100, offset = offset, with_payload = ["id"], with_vectors = False)
            records = scroll_response[0]
            offset = scroll_response[1]
            ids.extend([records[i].id for i in range(len(records))])
        return ids 
    
    def rerank(self, model = "/home/fengtang/hugging_face_model/bge-reranker-v2-m3", top_n = None):
        
        print(f"{'='*10} Rerank retrieved documents {'='*10}")
        self.rerank_model = model
        rerank_model = HuggingFaceCrossEncoder(model_name = model)
        if top_n == None:
            top_n = int(np.around(0.5 * len(self.retrieve_docs)))
        else:
            top_n = top_n
        print(f"Only top {top_n} members after reranking out of {len(self.retrieve_docs)} retrieved documents remained")
        reranker = CrossEncoderReranker(model = rerank_model, top_n =  top_n)
        
        if self.hybrid_search_cond:
            reranker_retriever = ContextualCompressionRetriever(base_compressor = reranker, base_retriever = self.hybrid_retriever)  
        else:
            reranker_retriever = ContextualCompressionRetriever(base_compressor = reranker, base_retriever = self.retriever)  
        self.retriever = reranker_retriever
        self.retrieve_docs_after_rerank = reranker_retriever.invoke(self.gene_list)

    def hybrid_search(self, hierarchy_search = False, key = None, value = None):

        print(f"{'='*10} Hybrid search to find relevant documents {'='*10}")
        
        ids = self.get_all_ids()
        docs = self.vector_store.get_by_ids(ids) # [self.vector_store.get_by_ids([id])[0] for id in ids]

        if hierarchy_search:
            docs_filter = Filter(
                                must=[
                                    FieldCondition(
                                        key=f"metadata.{key}",
                                        match=MatchValue(value=value)
                                    )
                                ]
                            )
            
            self.doc_filter = docs_filter
            dense_retriever = self.vector_store.as_retriever(search_kwargs = {"k": 10, "filter": docs_filter})
            self.dense_retriver.append(dense_retriever)
            self.retrieve_docs_by_dense = self.vector_store.as_retriever(search_kwargs = {"filter": docs_filter, "k": 10}).invoke(self.gene_list)
            
            filter_docs = [doc for doc in docs if doc.metadata.get(key) == value]
            bm25_retriever = BM25Retriever.from_documents(documents = filter_docs)
            bm25_retriever.k = 10

            hybrid_retriever = EnsembleRetriever(
                    retrievers = [bm25_retriever, dense_retriever],
                    weights = [0.5, 0.5]
                                )
        
        else:
            bm25_retriever = BM25Retriever.from_documents(documents = docs)
            bm25_retriever.k = 10
            hybrid_retriever = EnsembleRetriever(
                retrievers = [bm25_retriever, self.retriever],
                weights = [0.5, 0.5]
            )
        
        self.retrieve_docs_by_BM25 = bm25_retriever.invoke(self.gene_list)


        
        self.hybrid_retriever = hybrid_retriever
        self.retrieve_docs = hybrid_retriever.invoke(self.gene_list)
        print(f"{len(self.retrieve_docs)} relevant documents found by hybrid search")

    def summary_res(self, res):
         system = """
                    You are an advanced AI summarizer. 

                    Your task is to summarize cell type annotation evidences from two levels provided by users to give a comprehensive answer.

                    input requirement:
                        - level 1: {level1} 
                        - level 2: {level2}

                    output requirement:
                    - Format your response with a clear structure: 
                        - an introduction stating the inferred cell type at two levels (e.g. B cell in levele 1, stress B cell in level 2)
                        - a detailed evidence section with subheadings (e.g., "Key Evidence from the Gene List," "Key Datasets Supporting This Conclusion"). In each subsheading,  another subheading to distinguish the information from two levels. 
                        - Use bold formatting for emphasis (e.g., **cell type**, **gene names**) and maintain a scientific tone.
                    """
         prompt = PromptTemplate(
                    template = system,
                    input_variables = ["level1", "level2"]
                )
        
         chain = prompt | self.llm 
            
         summary_res = chain.invoke({"level1": res[len(res)-2], "level2": res[len(res)-1]})
         summary_res = summary_res.content
         self.answer = summary_res
         return summary_res
    
    def get_answer(self, gene_list, query = None, hierarchy_search = True, hybrid_search = True, rerank_model = None, auto_translate = True, auto_refine_query = True, highlight_docs = True, score_docs = True):
        set_llm_cache(InMemoryCache())
        if query == None:
            self.query = f""" What cell type does the provided gene list suggest? \n\n gene list: {gene_list}:
          
                            Output Requirement:
                                - Break down the evidence into categories (e.g., key gene groups like surfactant proteins, epithelial markers, or supporting genes) and explain their relevance to potential cell types.
                                - Incorporating inline citations in the format such as [Glasner_2023 (F0)] in the context that might support your inference.
                                - Rule out alternative cell types by noting the absence of conflicting markers (e.g., SCGB1A1 for Clara cells).
                                - Provide a clear conclusion identifying the most likely cell type, supported by the combined evidence.
                                - Format your response with a clear structure: an introduction stating the inferred cell type, a detailed evidence section with subheadings 
                                (e.g., "Key Evidence from the Gene List," "Key Datasets Supporting This Conclusion"), and a concise final answer. 
                                Use bold formatting for emphasis (e.g., **cell type**, **gene names**) and maintain a scientific tone.
                          """
        else:
            self.query = query

        self.gene_list = gene_list
        self.hybrid_search_cond = hybrid_search
        
        if auto_translate: 
            self.translator()
        if auto_refine_query:
            self.refine_query()
        
        if hybrid_search:
            if hierarchy_search:
                for i in ["M", "F"]:
                    self.hybrid_search(hierarchy_search = True, key = "level", value = i)
                    self.rerank()
                    if score_docs:
                        self.score_documents()
                        self.answer.append(self.rag_chain.invoke({"question": self.query, "context": self.format_docs(self.retrieve_docs)}))
                    else:
                        self.answer.append(self.rag_chain.invoke({"question": self.query, "context": self.retrieve_docs_after_rerank}))
                self.answer = self.summary_res(self.answer)
            else:
                self.hybrid_search(hierarchy_search = False, key = None, value = None)
                self.rerank()
                if score_docs:
                    self.score_documents()
                    self.answer = self.rag_chain.invoke({"question": self.query, "context": self.format_docs(self.retrieve_docs)})
                else:
                    self.answer = self.rag_chain.invoke({"question": self.query, "context": self.format_docs(self.retrieve_docs_after_rerank)})
        else:
            self.retrieve_docs.extend(self.retriever.invoke(self.query))
            self.rerank()
            if score_docs:
                self.score_documents()
                self.answer = self.rag_chain.invoke({"question": self.query, "context": self.format_docs(self.retrieve_docs)})
            else:
                self.answer = self.rag_chain.invoke({"question": self.query, "context": self.retrieve_docs_after_rerank})
        
        if highlight_docs:
            self.highlight_docs()
        
        return self.answer
    
    def highlight_docs(self):
        parser = PydanticOutputParser(pydantic_object=HighlightDocuments)
        system = """You are an advanced assistant for document search and retrieval. You are provided with the following:
                    1. A question.
                    2. A generated answer based on the question.
                    3. A set of documents that were referenced in generating the answer.

                    Your task is to identify and extract the exact inline segments from the provided documents that directly correspond to the content used to 
                    generate the given answer. The extracted segments must be verbatim snippets from the documents, ensuring a word-for-word match with the text 
                    in the provided documents.

                    Ensure that:
                    - (Important) Each segment is an exact match to a part of the document and is fully contained within the document text.
                    - The relevance of each segment to the generated answer is clear and directly supports the answer provided.
                    - (Important) If you didn't used the specific document don't mention it.

                    Used documents: <docs>{documents}</docs> \n\n User question: <question>{question}</question> \n\n Generated answer: <answer>{generation}</answer>

                    <format_instruction>
                    {format_instructions}
                    </format_instruction>
                    """
        prompt = PromptTemplate(
                template = system,
                input_variables = ["documents", "question", "generation"],
                partial_variables = {"format_instructions": parser.get_format_instructions()},
            )
        chain = prompt | self.llm | parser
        relevant_segments =  chain.invoke(
            {"documents": self.format_docs(self.retrieve_docs),
             "question": self.query,
             "generation": self.answer}
        )
        self.relevant_segments = relevant_segments
        return relevant_segments

    def get_relevant_segments(self):
        ls = list()
        for id, title, source, segment in zip(
                self.relevant_segments.id,
                self.relevant_segments.title,
                self.relevant_segments.source,
                self.relevant_segments.segment,
            ):
            ls.append(f"ID: {id}\nTitle: {title}\nSource: {source}\nText Segment: {segment}\n")
        return ls
    

class ConnectRAG(SimpleRAG):

    def __init__(self, host, api_key, location, url, collection_name, embedding_model):
        super().__init__(file_path="", file_type="")
        self.qdrant_host = host
        self.__qdrant_api_key = api_key
        self.qdrant_collection_name = collection_name
        self.text_embedding_model  = embedding_model
        self.qdrant_location = location
        self.qdrant_url = url

    def __repr__(self):
        info = f"""
        ConnectRAG object

        Qdrant client: {self.qdrant_host}
        Collection name: {self.qdrant_collection_name}
        Text embedding model: {self.text_embedding_model}
        LLM model: {self.llm_model}
        """
        return info

    def connect_client(self, model_kwargs = {"device": "cpu"}, encode_kwargs = {"normalize_embeddings": True}, timeout = 1000):
        client = QdrantClient(host = self.qdrant_host, api_key = self.__qdrant_api_key, timeout = timeout, url = self.qdrant_url, location = self.qdrant_location)
        model = HuggingFaceEmbeddings(model_name =  self.text_embedding_model, model_kwargs= model_kwargs, encode_kwargs=encode_kwargs)
        vector_store = QdrantVectorStore(client = client, collection_name = self.qdrant_collection_name, embedding = model)
        self.vector_store = vector_store
        self.qdrant_client = client

    def run_rag(self,
                llm_model,
                llm_api_key,
                llm_base_url,
                qdrant_kwgars = {"k": 10}):
        self.create_rag_chain(vector_store = self.vector_store, model = llm_model, api_key = llm_api_key, base_url = llm_base_url, search_kwargs = qdrant_kwgars)
