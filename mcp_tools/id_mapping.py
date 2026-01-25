# %%
from mcp.server.fastmcp import FastMCP
import mygene

# %%
mcp = FastMCP("id_mapping")

# %%
@mcp.tool()
async def id_mapping(gene_id_list: list, species:str = "human"):
    """
    Convert Entrez Gene IDs or Ensembl IDs to gene symbols using MyGene.info API
    
    Parameters:
    -----------
    gene_id_list : list or str
        Single gene ID or list of gene IDs (Entrez or Ensembl)
    species : str
        Species (default: human)
    
    Returns:
    --------
    list : Mapping of gene_id to gene_symbol.
    """
    mg = mygene.MyGeneInfo()
    result = mg.querymany(gene_id_list, scopes='entrezgene,ensembl.gene', 
                      fields='symbol,name', species=species)
    
    return [i["symbol"] for i in result]


# %%

if __name__ == "__main__":
    # Initialize and run the server
    mcp.run(transport='stdio')


