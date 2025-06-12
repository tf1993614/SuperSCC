# %%
from rag import *
import streamlit as st

# --- Page Configuration ---
# st.set_page_config(layout="wide", page_title="SuperSCC-RAG Interface")
st.title("SuperSCC-RAG Interactive Interface")

# --- Initialize Session State ---
if "rag_instance" not in st.session_state:
    st.session_state.rag_instance = None
if "rag_initialized" not in st.session_state:
    st.session_state.rag_initialized = False
if "answer" not in st.session_state:
    st.session_state.answer = None
if "relevant_segments" not in st.session_state:
    st.session_state.relevant_segments = None

# --- Sidebar for Configuration ---
st.sidebar.header("RAG System Configuration")
setup_mode = st.sidebar.radio("Choose RAG Setup Mode:", ("Create New RAG from Local Files", "Connect to Existing RAG"))

# --- LLM Configuration ---
st.sidebar.subheader("LLM Configuration")
llm_model_name = st.sidebar.text_input("LLM Model Name (e.g., Qwen/Qwen2-7B-Instruct, gpt-3.5-turbo)", "Qwen/Qwen2-7B-Instruct")
llm_api_key = st.sidebar.text_input("LLM API Key", type="password", help="Required for OpenAI or commercial LLMs")
llm_api_base = st.sidebar.text_input("LLM API Base URL (e.g., for local LLMs like Ollama/LMStudio)", "http://localhost:11434/v1", help="For OpenAI, leave as default or set to official endpoint if needed.")


if setup_mode == "Create New RAG from Local Files":
    st.sidebar.subheader("Local Data Source")
    file_path_input = st.sidebar.text_input("File/Directory Path for RAG Data", "/path/to/your/data_directory_or_file.csv")
    file_type_input = st.sidebar.selectbox("File Type", ("csv", "pdf"), help="Select the type of files in the directory or the type of the single file.")

    st.sidebar.subheader("Embedding & Vector Store (Qdrant)")
    embedding_model_name = st.sidebar.text_input("Text Embedding Model (HuggingFace)", "BAAI/bge-m3")
    qdrant_url_local = st.sidebar.text_input("Qdrant URL (local/cloud)", "http://localhost:6333")
    qdrant_collection_name_local = st.sidebar.text_input("Qdrant Collection Name", "superscc_rag_collection")
    # qdrant_api_key_local = st.sidebar.text_input("Qdrant API Key (if cloud/secured)", type="password", key="qdrant_api_local") # Optional

    # Advanced options for SimpleRAG
    with st.sidebar.expander("Advanced Data Processing Options"):
        metadata_columns_str = st.text_input("CSV Metadata Columns (comma-separated)", "Data,Cluster,Organization,cell_type", help="For CSV files, columns to include in metadata.")
        chunk_size = st.number_input("PDF Chunk Size", min_value=100, max_value=4000, value=1000)
        chunk_overlap = st.number_input("PDF Chunk Overlap", min_value=0, max_value=1000, value=200)
        embedding_vector_size = st.number_input("Embedding Vector Size", min_value=128, max_value=4096, value=1024, help="Ensure this matches the chosen embedding model's output dimension.")


    if st.sidebar.button("Initialize New RAG System", key="init_simple_rag"):
        if not file_path_input or not os.path.exists(file_path_input):
            st.sidebar.error("Please provide a valid File/Directory Path.")
        elif not llm_model_name or not llm_api_base: # API key might be optional for local LLMs
            st.sidebar.error("Please configure LLM Model Name and API Base URL.")
        else:
            with st.spinner("Initializing SimpleRAG... This may take a while depending on data size."):
                try:
                    rag = SimpleRAG(file_path=file_path_input, file_type=file_type_input)
                    
                    metadata_cols_list = [col.strip() for col in metadata_columns_str.split(',') if col.strip()] if metadata_columns_str else None

                    rag.run_rag(
                        qdrant_location=None, # For local Qdrant, URL is primary
                        qdrant_url=qdrant_url_local,
                        # qdrant_api_key=qdrant_api_key_local if qdrant_api_key_local else None,
                        text_embedding_model=embedding_model_name,
                        llm_model=llm_model_name,
                        llm_api_key=llm_api_key if llm_api_key else "ollama", # ollama doesn't need key
                        llm_base_url=llm_api_base,
                        qdrant_collection_name=qdrant_collection_name_local,
                        metadata_columns=metadata_cols_list,
                        chunk_size=chunk_size,
                        chunk_overlap=chunk_overlap,
                        vectors_config=VectorParams(size=embedding_vector_size, distance=Distance.COSINE),
                        # text_model_kwargs = {"device": "cuda" if torch.cuda.is_available() else "cpu"} # Example
                    )
                    st.session_state.rag_instance = rag
                    st.session_state.rag_initialized = True
                    st.sidebar.success("SimpleRAG Initialized!")
                    st.toast(f"RAG system initialized with collection: {qdrant_collection_name_local}")
                    print(rag) # For debugging, prints to console
                except Exception as e:
                    st.sidebar.error(f"Error initializing SimpleRAG: {e}")
                    st.error(f"Detailed Error: {str(e)}") # Show error in main panel too


elif setup_mode == "Connect to Existing RAG":
    st.sidebar.subheader("Existing Qdrant Vector Store")
    qdrant_url_connect = st.sidebar.text_input("Qdrant URL", "http://localhost:6333", key="q_url_conn")
    qdrant_collection_name_connect = st.sidebar.text_input("Qdrant Collection Name", "superscc_rag_collection", key="q_coll_conn")
    embedding_model_connect = st.sidebar.text_input("Text Embedding Model (used for the collection)", "BAAI/bge-m3", key="emb_conn")
    qdrant_api_key_connect = st.sidebar.text_input("Qdrant API Key (if secured)", type="password", key="q_api_conn")
    qdrant_host_connect = st.sidebar.text_input("Qdrant Host (optional, if not part of URL)", key="q_host_conn")
    qdrant_location_connect = st.sidebar.text_input("Qdrant Location (e.g., :memory: or path for local)", key="q_loc_conn", help="Usually for on-disk local, otherwise URL is preferred.")


    if st.sidebar.button("Connect to Existing RAG System", key="connect_rag"):
        if not qdrant_collection_name_connect or not embedding_model_connect or not llm_model_name or not llm_api_base:
            st.sidebar.error("Please fill all required Qdrant and LLM fields.")
        else:
            with st.spinner("Connecting to RAG system..."):
                try:
                    rag = ConnectRAG(
                        host=qdrant_host_connect if qdrant_host_connect else None,
                        api_key=qdrant_api_key_connect if qdrant_api_key_connect else None,
                        location=qdrant_location_connect if qdrant_location_connect else None,
                        url=qdrant_url_connect,
                        collection_name=qdrant_collection_name_connect,
                        embedding_model=embedding_model_connect
                    )
                    rag.connect_client() # Initialize Qdrant client and vector store
                    # For ConnectRAG, run_rag mainly sets up the RAG chain
                    rag.run_rag(
                        llm_model=llm_model_name,
                        llm_api_key=llm_api_key if llm_api_key else "ollama",
                        llm_base_url=llm_api_base
                        # qdrant_kwgars can be added if needed for retriever settings
                    )
                    st.session_state.rag_instance = rag
                    st.session_state.rag_initialized = True
                    st.sidebar.success("Connected to RAG System!")
                    st.toast(f"Connected to RAG collection: {qdrant_collection_name_connect}")
                    print(rag) # For debugging
                except Exception as e:
                    st.sidebar.error(f"Error connecting to RAG: {e}")
                    st.error(f"Detailed Error: {str(e)}")


# --- Main Panel for Querying and Results ---
if st.session_state.rag_initialized and st.session_state.rag_instance:
    st.header("Query the RAG System")

    # Custom Prompt (Optional)
    st.subheader("Customize RAG Prompt (Optional)")
    default_prompt_template_text = """
     Answer the question with using following the context and relevant information : 
     {context}

    Question: {question}
    """
    custom_prompt_text = st.text_area("RAG Prompt Template (uses {context} and {question})", 
                                      value=default_prompt_template_text, height=250)

    if st.button("Update RAG Prompt", key="update_prompt_btn"):
        if custom_prompt_text:
            try:
                custom_prompt = ChatPromptTemplate.from_template(custom_prompt_text)
                st.session_state.rag_instance.update_rag_chain(prompt=custom_prompt)
                st.success("RAG prompt updated successfully!")
            except Exception as e:
                st.error(f"Error updating prompt: {e}")
        else:
            # Optionally revert to default prompt if text area is cleared
            # This would require storing/recreating the original default prompt
            st.warning("Prompt text is empty. No changes made or revert to default if implemented.")


    st.subheader("Enter Your Query")
    query_text = st.text_area("Main Question / Task Description:", 
                              """What cell type does the provided gene list suggest?
                                 Output Requirement:
                                    - Break down the evidence into categories (e.g., key gene groups like surfactant proteins, epithelial markers, or supporting genes) and explain their relevance to potential cell types.
                                    - Incorporating inline citations in the format [Glasner_2023 (F0)] in the context that might support your inference.
                                    - Rule out alternative cell types by noting the absence of conflicting markers (e.g., SCGB1A1 for Clara cells).
                                    - Provide a clear conclusion identifying the most likely cell type, supported by the combined evidence.
                                    - Format your response with a clear structure: an introduction stating the inferred cell type, a detailed evidence section with subheadings 
                                    (e.g., "Key Evidence from the Gene List," "Key Datasets Supporting This Conclusion"), and a concise final answer. 
                                    Use bold formatting for emphasis (e.g., **cell type**, **gene names**) and maintain a scientific tone""",
                              height=100)
    gene_list_text = st.text_input("Gene List for Retrieval (comma-separated, e.g., CD74,CD4,NKG7):", "SFTPA1,SFTPB,SLPI,EPCAM,CLDN7,IFT57,LPCAT1")

    full_query_for_llm = f"{query_text}\n\nGene List: {gene_list_text}"
    st.text_area("Full query to be sent to LLM (after potential refinement):", value=full_query_for_llm, height=100, disabled=True)
    
    col1, col2 = st.columns(2)
    with col1:
        st.write("Query Options:")
        opt_hybrid_search = st.checkbox("Enable Hybrid Search (BM25 + Dense)", value=True)
        opt_hierarchy_search = st.checkbox("Enable Hierarchy Search to get hierarchical cell labels", value=True)
        opt_score_docs = st.checkbox("Enable LLM to evaulate whether the retirval segement is relevant with the query", value=True)
        opt_auto_translate = st.checkbox("Enable Auto Query Translation (Chinese -> English)", value=False) # Default False as it's specific
        opt_auto_refine = st.checkbox("Enable Auto Query Refinement", value=False)
        opt_highlight_docs = st.checkbox("Highlight Relevant Document Segments", value=False)
    
    with col2:
        st.write("Reranking:")
        # Default from your notebook was a local path, using a common HF model name as a placeholder
        rerank_model_name_input = st.text_input("Reranker Model (HuggingFace name or path)", "BAAI/bge-reranker-large")

        st.write("Change LLM model:")
        custom_model = st.text_input("LLM model (model name)", llm_model_name)
        if custom_model != llm_model_name:
            try:
                st.session_state.rag_instance.update_rag_chain(model=custom_model)
                st.success("LLM model updated successfully!")
            except Exception as e:
                    st.error(f"Error updating prompt: {e}")


    if st.button("Get Answer", key="get_answer_btn"):
        if not query_text: # gene_list_text can be optional if query_text already has genes
            st.error("Please enter a main question.")
        else:
            with st.spinner("Retrieving and generating answer..."):
                try:
                    rag = st.session_state.rag_instance
                    # The 'query' param for get_answer should be the full text for the LLM
                    # The 'gene_list' param is specifically for retrieval/reranking stages
                    
                    final_query_to_llm = f"{query_text}\n\nGene List: {gene_list_text}" # Construct as per example
                    if rag.query != final_query_to_llm : # If user changed input after last run
                        rag.query = final_query_to_llm # Update internal query state if needed before refine/translate

                    answer = rag.get_answer(
                        query=final_query_to_llm, # This becomes self.query, used for refine/translate and final LLM
                        gene_list=gene_list_text, # This is self.gene_list, used for hybrid search & rerank
                        hybrid_search = opt_hybrid_search,
                        hierarchy_search = opt_hierarchy_search,
                        rerank_model = rerank_model_name_input,
                        auto_translate=opt_auto_translate,
                        auto_refine_query=opt_auto_refine,
                        highlight_docs=opt_highlight_docs,
                        score_docs=opt_score_docs 
                    )
                    st.session_state.answer = answer
                    if opt_highlight_docs and hasattr(rag, 'relevant_segments') and rag.relevant_segments:
                        st.session_state.relevant_segments = rag.get_relevant_segments()
                    else:
                        st.session_state.relevant_segments = None
                    
                    st.success("Answer generated!")

                except Exception as e:
                    st.error(f"Error during RAG execution: {e}")
                    st.error(f"Detailed Error: {str(e)}") # More detailed error

    if st.session_state.answer:
        st.subheader("Generated Answer")
        st.markdown(st.session_state.answer)

    if st.session_state.relevant_segments:
        st.subheader("Relevant Document Segments")
        for i, segment_info in enumerate(st.session_state.relevant_segments):
            with st.expander(f"Segment {i+1} - Source: {segment_info.splitlines()[2].replace('Source: ','')} | Title: {segment_info.splitlines()[1].replace('Title: ','')}"):
                st.markdown(segment_info)

elif not st.session_state.rag_initialized:
    st.info("Please configure and initialize or connect to a RAG system using the sidebar.")

# --- Display RAG object info for debugging ---
if st.session_state.rag_instance:
    with st.expander("Current RAG Instance Details (for debugging)"):
        st.json(st.session_state.rag_instance.__repr__()) # Needs __repr__ to be JSON serializable or use st.text


