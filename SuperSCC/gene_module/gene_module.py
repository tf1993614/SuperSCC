"""
Gene module for functions specific to gene data analysis.
"""

from langchain_community.chat_models import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.schema import StrOutputParser


import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import warnings

# Suppress pandas warnings about fragmented dataframes, which can occur during column additions.
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def find_overlap_singlet(data, full_search=False):
    """
    Calculates the overlap (intersection size) between all pairs of columns in a DataFrame.
    This is a single-threaded version.

    Args:
        data (pd.DataFrame): DataFrame where each column is a set of items.
        full_search (bool): If True, compares every column with every other column.
                            If False, only compares each column with subsequent columns to avoid duplicates.

    Returns:
        pd.DataFrame: A DataFrame with columns 'number', 'base_name', 'compare_name'.
    """
    num_cols = data.shape[1]
    results = []
    
    with tqdm(total=num_cols, desc="Finding overlaps (single-threaded)") as pbar:
        for i in range(num_cols):
            base_name = data.columns[i]
            base_set = set(data[base_name].dropna())
            
            if full_search:
                compare_cols = data.columns
            else:
                compare_cols = data.columns[i + 1:]

            for compare_name in compare_cols:
                compare_set = set(data[compare_name].dropna())
                intersection_size = len(base_set.intersection(compare_set))
                
                results.append({
                    'number': intersection_size,
                    'base_name': base_name,
                    'compare_name': compare_name
                })
            pbar.update(1)
            
    return pd.DataFrame(results)

def _find_overlap_base_worker(index, data, full_search=False):
    """
    Worker function to calculate pairwise intersection size for one base column.
    
    Args:
        index (int): The column index to use as the base group.
        data (pd.DataFrame): The full DataFrame of gene sets.
        full_search (bool): If False, compares the base column only with subsequent columns.

    Returns:
        pd.DataFrame: A DataFrame of intersection results for the given base column.
    """
    base_name = data.columns[index]
    base_set = set(data[base_name].dropna())
    
    if full_search:
        compare_cols = data.columns
    else:
        # To match R's `seq(index)`, which includes the index itself, we use `index:`
        # But the original R code `! colnames(data) %in% colnames(data)[seq(num)]`
        # actually means comparing with columns *after* the current one.
        compare_cols = data.columns[index:]

    results = []
    for compare_name in compare_cols:
        compare_set = set(data[compare_name].dropna())
        intersection_size = len(base_set.intersection(compare_set))
        results.append({
            'number': intersection_size,
            'base_name': base_name,
            'compare_name': compare_name
        })
    return pd.DataFrame(results)


def _find_overlap_founder_worker(index, data):
    """
    Worker function to calculate intersection size between a 'founder' group and other groups.

    Args:
        index (int): The column index of the comparison group (in data without 'founder').
        data (pd.DataFrame): The DataFrame containing a 'founder' column.

    Returns:
        pd.DataFrame: A one-row DataFrame with the intersection result.
    """
    founder_set = set(data['founder'].dropna())
    base_name = "founder"
    
    compare_cols = data.columns.drop("founder")
    compare_name = compare_cols[index]
    compare_set = set(data[compare_name].dropna())
    
    intersection_size = len(founder_set.intersection(compare_set))
    
    return pd.DataFrame([{
        'number': intersection_size,
        'base_name': base_name,
        'compare_name': compare_name
    }])


# ==============================================================================
# Parallel Intersection Functions
# ==============================================================================

def do_intersection_base(data, parallel_num=8):
    """
    A function to run find_overlap_base_worker in parallel.

    Args:
        data (pd.DataFrame): DataFrame where each column is a set of items.
        parallel_num (int): Number of parallel processes to use. If None, runs single-threaded.

    Returns:
        pd.DataFrame: A concatenated DataFrame of all pairwise intersection results.
    """
    if not parallel_num:
        # R code calls find_overlap_singlet without full_search, but its logic
        # is equivalent to the parallel version which avoids duplicate pairs.
        # The parallel version compares i with j>=i, then we filter i==j later.
        # This is slightly different than R, but achieves same goal of finding best pair.
        return find_overlap_singlet(data, full_search=False)

    worker_func = partial(_find_overlap_base_worker, data=data, full_search=False)
    
    with Pool(processes=parallel_num) as pool:
        results = list(tqdm(pool.imap(worker_func, range(data.shape[1])), 
                            total=data.shape[1], 
                            desc="Finding overlaps (parallel)"))
        
    intersect_data = pd.concat(results, ignore_index=True)
    # Filter out self-comparisons
    intersect_data = intersect_data[intersect_data['base_name'] != intersect_data['compare_name']].reset_index(drop=True)
    return intersect_data

def do_intersection_founder(data, parallel_num=8):
    """
    A function to run find_overlap_founder_worker in parallel.

    Args:
        data (pd.DataFrame): DataFrame containing a 'founder' column and other gene sets.
        parallel_num (int): Number of parallel processes.

    Returns:
        pd.DataFrame: A concatenated DataFrame of intersection results with the founder.
    """
    worker_func = partial(_find_overlap_founder_worker, data=data)
    num_compare_cols = data.shape[1] - 1
    
    with Pool(processes=parallel_num) as pool:
        results = list(tqdm(pool.imap(worker_func, range(num_compare_cols)), 
                            total=num_compare_cols, 
                            desc="Finding founder overlaps (parallel)"))

    return pd.concat(results, ignore_index=True)

# ==============================================================================
# Main Gene Module Logic
# ==============================================================================

def core_get_gene_module(data, intersect_size=10, intersect_group_size=5,
                         parallel_num=8, init_signal=True):
    """
    A recursive underlying function to identify and merge gene modules.
    """
    # In pandas, it's more efficient to work with split columns from the start if needed
    # The R code `str_remove(.x, "/.+")` keeps the part before the slash
    pre_data_for_intersection = data.apply(lambda col: col.str.split('/').str[0], axis=0)
    
    if init_signal:
        intersect_data = do_intersection_base(data=pre_data_for_intersection, parallel_num=parallel_num)
    else:
        intersect_data = do_intersection_founder(data=pre_data_for_intersection, parallel_num=parallel_num)
        
    intersect_data = intersect_data[intersect_data['number'] > intersect_size]
    intersect_data = intersect_data.sort_values('number', ascending=False, ignore_index=True, kind = "stable")

    # if not intersect_data.empty and intersect_data.shape[0] >= intersect_group_size:
    if intersect_data.shape[0] >= intersect_group_size:

        # 1. Identify the best pair to merge
        max_overlap_row = intersect_data.iloc[0]
        member1_name, member2_name = max_overlap_row['base_name'], max_overlap_row['compare_name']
        max_overlap_members = [member1_name, member2_name]

        # 2. Get gene sets (intersection, union, difference)
        genes1 = set(pre_data_for_intersection[member1_name].dropna())
        genes2 = set(pre_data_for_intersection[member2_name].dropna())
        
        intersection_genes = genes1.intersection(genes2)
        union_genes = genes1.union(genes2)
        diff_genes = union_genes.difference(intersection_genes)


        # 3. Create ranked gene pools from original data (with scores)
        def create_ranking_pool(data, col_name):
            pool = data[col_name].dropna().str.split('/', n=1, expand=True)
            pool.columns = ['genes', 'scores']
            # pool['scores'] = pd.to_numeric(pool['scores'])
            return pool

        genes_ranking_pool_1 = create_ranking_pool(data, member1_name)
        genes_ranking_pool_2 = create_ranking_pool(data, member2_name)
        
        # 4. Rank intersection genes by highest score from either pool
        intersection_genes_ranking_1 = genes_ranking_pool_1[genes_ranking_pool_1['genes'].isin(intersection_genes)]
        intersection_genes_ranking_2 = genes_ranking_pool_2[genes_ranking_pool_2['genes'].isin(intersection_genes)]
        
        intersection_genes_ranking_final = pd.concat([intersection_genes_ranking_1, intersection_genes_ranking_2])
        intersection_genes_ranking_final = intersection_genes_ranking_final.sort_values('scores', ascending=False, kind = "stable").drop_duplicates('genes', keep='first')


        # 5. Build the new 50-gene "founder" module
        # The R code uses a fixed size of 50.
        TARGET_SIZE = 50
        final_top50_ranking = None

        if len(intersection_genes) != TARGET_SIZE: 
            # Need to fill remaining slots from difference genes
            robust_ranking = intersection_genes_ranking_final
            
            diff_genes_ranking_1 = genes_ranking_pool_1[genes_ranking_pool_1['genes'].isin(diff_genes)]
            diff_genes_ranking_2 = genes_ranking_pool_2[genes_ranking_pool_2['genes'].isin(diff_genes)]
            
            final_ranking = pd.concat([intersection_genes_ranking_final, diff_genes_ranking_1, diff_genes_ranking_2])
            final_ranking = final_ranking.sort_values("scores", ascending=False, kind = "stable")

            border_ranking = final_ranking.loc[final_ranking.genes.isin(diff_genes), :]
            num_needed = TARGET_SIZE - len(intersection_genes)
            border_ranking = border_ranking.head(num_needed)

            robust_ranking = final_ranking.loc[final_ranking.genes.isin(intersection_genes), :]
            
            final_top50_ranking = pd.concat([border_ranking, robust_ranking])
        else:
            final_top50_ranking = intersection_genes_ranking_final.head(TARGET_SIZE)
           
        
        # Format back to "gene/score" string
        founder = final_top50_ranking.apply(lambda row: f"{row['genes']}/{row['scores']}", axis=1)

        # 6. Prepare data for recursive call
        update_data = data.drop(columns=max_overlap_members)
        
        # Add the new founder column, padding with NaN if necessary
        founder_series = pd.Series(founder.values, name='founder')
        update_data = update_data.reset_index(drop=True) # Ensure index alignment
        update_data['founder'] = founder_series

        return core_get_gene_module(data = update_data, parallel_num = parallel_num, init_signal = False)
    else:
        # Base case: no more sufficiently large intersections found
        return data


def get_gene_module(data, intersect_size=10, intersect_group_size=5, parallel_num=8):
    """
    A function to iteratively find gene modules from a collection of gene sets.

    Args:
        data (pd.DataFrame): A DataFrame where each column contains a gene set as a
                             list of strings (e.g., "GENE/SCORE"). Columns should be padded
                             with NaN for unequal lengths.
        intersect_size (int): The minimum intersection size to consider merging two sets.
        intersect_group_size (int): The minimum number of high-intersection pairs required
                                    to proceed with a merge.
        parallel_num (int): The number of parallel processes to use.

    Returns:
        dict: A dictionary containing:
              - 'gene_module': A list of the identified gene modules (each a list of strings).
              - 'module_members': A list of the original column names that formed each module.
              - 'remained_gene_sets': The final DataFrame of gene sets that were not merged.
    """
    meta_program_list = []
    meta_program_members = []
    remained_gene_sets = []
    
    current_data = data.copy()
    
    iteration = 1
    while True:
        # num_cols_before = current_data.shape[1]
        # if num_cols_before <= 1:
        #     break

        output = core_get_gene_module(
            data=current_data,
            intersect_size=intersect_size,
            intersect_group_size=intersect_group_size,
            parallel_num=parallel_num,
            init_signal=True # Always start with a base search in the main loop
        )
        


        # If 'founder' is in the output, a module was created
        if "founder" not in output.columns:
                # No module was found, terminate the loop
                print("No more modules found satisfying the criteria.")
                break
        else:
            print(f"Finding the gene module {iteration}")
            
            # Store the identified module
            module = output['founder'].dropna().tolist()
            meta_program_list.append(module)
            
            # Identify which original members were consumed
            # update_data = output.drop(columns='founder')
            update_data = output.drop(columns="founder")
            consumed_members = list(set(current_data.columns) - set(update_data.columns))
            meta_program_members.append(consumed_members)
            
            remained_gene_sets.append(current_data)

            # Update data for the next iteration
            current_data = update_data
            iteration += 1
            
    return {
        'gene_module': meta_program_list,
        'module_members': meta_program_members,
        'remained_gene_sets': remained_gene_sets
    }


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
