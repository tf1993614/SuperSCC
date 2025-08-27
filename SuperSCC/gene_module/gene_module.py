"""
Gene module for functions specific to gene data analysis.
"""

from langchain_community.chat_models import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.schema import StrOutputParser
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


def get_gene_module(data, intersect_size = 10, intersect_group_size = 5, parallel_num = 8, lib_loc = None):
    """
    A function to find conserved gene modules across multiple gene sets.

    Parameters
    ----------
    data: 
        A data frame with each column representing a gene set list.
    intersection_size: 
        A int control the minimal intersection cutoff for adding a new gene set to the forming gene module.
    intersection_group_size: 
        A int control the minimal group size to consider for defining the first gene module. 
    parallel_num: 
        A int to control the number of thread used for the program.
    lib_loc:
        A string to specify the path where geneModule R package is installed. 
    """
    # convert pandas dataframe into R dataframe
    pandas2ri.activate()
    r_df = pandas2ri.py2rpy(data)
    
    geneModule = importr("geneModule", lib_loc = lib_loc)
    res = geneModule.get_gene_module(data = r_df, intersect_size = intersect_size, intersect_group_size = intersect_group_size, parallel_num = parallel_num)
    res = dict(zip(res.names, map(lambda x: x[0] if len(x) == 1 else x, res)))
    return res


def compare_gene_modules(module1, module2, api_key, model = "deepseek-chat", base_url = "https://api.deepseek.com/v1"):
    """
    Compare two gene modules and analyze their similarities and differences.
    
    Parameters
    ----------
    module1:
        A gene set list representing the first gene module.
    module2:
        A gene set list representing the second gene module.
    api_key:
        A sting for the api key of the LLM provider.
    model:
        A string for the LLM model name.
    base_url:
        A string to base URL for API requests.
    """
    # Genes are already converted to names
    module1_genes = module1
    module2_genes = module2

    # Find common and unique genes
    common_genes = set(module1_genes).intersection(set(module2_genes))
    unique_to_module1 = set(module1_genes) - set(module2_genes)
    unique_to_module2 = set(module2_genes) - set(module1_genes)

    # Create comparison template
    comparison_template = """You are a bioinformatics expert. Compare these two gene modules:
    
    Module 1: {module1_genes}
    Module 2: {module2_genes}
    
    Analyze:
    1. Common biological pathways between modules
    2. Unique pathways in each module
    3. Potential functional relationships between modules
    4. Disease associations shared between modules
    5. Tissue/cell type specificity differences
    
    Provide your analysis in clear, structured paragraphs."""

    # Create comparison prompt
    comparison_prompt = ChatPromptTemplate.from_template(comparison_template)

    # set model
    model = ChatOpenAI(
        model=model,
        temperature=0.7,
        openai_api_key=api_key,
        openai_api_base=base_url
    )

    # Run comparison analysis
    comparison_chain = comparison_prompt | model | StrOutputParser()
    comparison_result = comparison_chain.invoke(
        {
            "module1_genes": ", ".join(module1_genes),
            "module2_genes": ", ".join(module2_genes),
        }
    )

    return {
        "common_genes": list(common_genes),
        "unique_to_module1": list(unique_to_module1),
        "unique_to_module2": list(unique_to_module2),
        "comparison_analysis": comparison_result,
    }



def analyse_one_gene_module(module_genes, api_key, model = "deepseek-chat", base_url = "https://api.deepseek.com/v1"):
    """
    Analyze a single gene module using the DeepSeek model
    
    Parameters
    ----------
    module_genes:
        A gene set list representing the gene module.
    api_key:
        A sting for the api key of the LLM provider.
    model:
        A string for the LLM model name.
    base_url:
        A string to base URL for API requests.
    """
    # Create prompt template
    template = """You are a bioinformatics expert. Analyze this list of genes and provide a detailed functional interpretation of the gene module:
    {module_genes}
    
    Consider:
    1. Common biological pathways
    2. Cellular processes involved
    3. Potential tissue/cell type specificity
    4. Disease associations
    5. Functional relationships between genes
    
    Provide your analysis in clear, structured paragraphs."""
    prompt = ChatPromptTemplate.from_template(template)

    # set model
    model = ChatOpenAI(
        model=model,
        temperature=0.7,
        openai_api_key=api_key,
        openai_api_base=base_url
    )

    chain = prompt | model | StrOutputParser()

    # Run analysis for gene module
    analysis = chain.invoke({"module_genes": ", ".join(module_genes)})

    return analysis
