import re
import os 
from os.path import basename, splitext

import numpy as np
import pandas as pd
from uuid import uuid4, UUID
import textwrap
from dotenv import load_dotenv
from typing import List

from ..SuperSCC import list_files

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

from .rag import *
