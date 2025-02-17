# SuperSCC

SuperSCC (Super single cell clustering), a novel algorithm that utilizes machine-learning models to discover cell identities and gene modules from multiple datasets without the need of data integration. Of note, SuperSCC can be implemented both in cell lineage and cell state level, thereby building the hierarchy of cell programs with specific cell identity and gene modules. Such information has the great potential to identify the shared rare populations across datasets regardless of batch effect and benefits label transfer for mapping cell labels from reference to query. 

![img](https://github.com/tf1993614/SuperSCC/blob/main/img/workflow.png)

## Installation

To install SuperSCC, you should run the following code if you already get [SuperSCC's tar ball](https://github.com/tf1993614/SuperSCC/tree/main/dist)

        pip install SuperSCC.tar.gz

Usually, all dependencies should be downloaded and installed automatically. 

To install SuperSCC via GitHub, you can do::

        git clone https://github.com/tf1993614/SuperSCC/
        cd SuperSCC
        python setup.py install

Since the `get_gene_module` function was written in R, you also need install [geneModule](https://github.com/tf1993614/SuperSCC/tree/main/geneModule) R package.

        R
        install.packages("geneModule_0.1.tar.gz", repos = .libPaths()[1], type = "source")
        .libPaths()[1] # get the location where packages are installed
        q()

Note: `get_gene_module` function can be also implemented in Python environment via rpy2-wrappered function.

After installing correctly, you can load SuperSCC::

        import SuperSCC as scc

## Documentation

For how to use SuperSCC, please read the [SuperSCC's documentation](https://superscc.readthedocs.io/en/latest/index.html).

## Citation
