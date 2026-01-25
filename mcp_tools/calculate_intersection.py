# %%
from mcp.server.fastmcp import FastMCP
import re

# %%
mcp = FastMCP("calculate_intersection")


# %%
@mcp.tool()
async def calculate_intersection(evaluation_set: list, refernece_set: list):
    """
    Calculate intersection statistics between an evaluation set and a reference set.

    Computes the total size of the evaluation set, the number of shared elements,
    and the intersection ratio (proportion of evaluation set elements found in
    the reference set)

    Parameters
    ----------
    evaluation_set : list
        A list of strings containing gene symbols in the evaluation set.
        Example: ["TP53", "BRCA1", "EGFR", "MYC"]
    refernece_set : list
        A list of strings containing gene symbols in the reference set.
        Example: ["EGFR", "MYC", "KRAS", "BRAF"]

    Returns
    -------
    tuple
        A 3-tuple containing:
        - int: Total number of elements in evaluation_set
        - int: Number of elements in the intersection (shared between both sets)
        - float: Intersection ratio (intersection_size / len(evaluation_set)).
          Returns 0.0 if evaluation_set is empty. Range: [0.0, 1.0]
    """
    intersection_size = len(set(evaluation_set).intersection(refernece_set))
    return len(evaluation_set), intersection_size, intersection_size / len(evaluation_set)


# %%
if __name__ == "__main__":
    # Initialize and run the server
    mcp.run(transport='stdio')


