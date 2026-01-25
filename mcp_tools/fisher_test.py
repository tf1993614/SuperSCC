# %%
from scipy.stats import fisher_exact
from mcp.server.fastmcp import FastMCP

# %%
mcp = FastMCP("fisher_test")

# %%
@mcp.tool()
async def fisher_test(data, alternative = "two.sided"):
    """
    Perform Fisher's exact test on a 2x2 contingency table.

    Parameters
    ----------
    data : array_like, shape (2, 2)
        A 2x2 contingency table in the format:
        
                        Condition Positive    Condition Negative
        Group 1             a                    b
        Group 2             c                    d
        
        For example:
        >>> data = [[5, 15],    # Group 1: 5 positive, 15 negative
        ...         [15, 5]]    # Group 2: 15 positive, 5 negative
        
        Where:
        - a (data[0][0]): Count of Group 1 with condition present
        - b (data[0][1]): Count of Group 1 with condition absent  
        - c (data[1][0]): Count of Group 2 with condition present
        - d (data[1][1]): Count of Group 2 with condition absent
        
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis (default: 'two-sided'):
        
        - 'two-sided': The odds ratio of the two groups is not equal to 1
        - 'less': The odds ratio is less than 1 (Group 1 has lower odds)
        - 'greater': The odds ratio is greater than 1 (Group 1 has higher odds)

    Returns
    -------
    p_value : float
        The p-value under the null hypothesis that the odds ratio is 1.
    """
    return fisher_exact(data, alternative = alternative)[1]


# %%

if __name__ == "__main__":
    # Initialize and run the server
    mcp.run(transport='stdio')


