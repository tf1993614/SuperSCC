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


After installing correctly, you can load SuperSCC

        import SuperSCC as scc


# Bonus for R person
Since the gene module-relevant functions was also written in R, you can also do the gene module analysis in R environment aftering installing [geneModule](https://github.com/tf1993614/SuperSCC/tree/main/geneModule) R package.
	
        R
        install.packages("geneModule_0.1.tar.gz", repos = .libPaths()[1], type = "source")
        .libPaths()[1] # get the location where packages are installed
        q()


## Documentation

For how to use SuperSCC, please read the [SuperSCC's documentation](https://superscc.readthedocs.io/en/latest/index.html).

## Citation

	@article {Tang2025.03.12.642774,
	author = {Tang, Feng and Zhang, Zhongmin and Zhou, Weige and Li, Guangpeng and Tian, Luyi},
	title = {Unveiling Gene Modules at Atlas Scale through Hierarchical Clustering of Single-Cell Data},
	elocation-id = {2025.03.12.642774},
	year = {2025},
	doi = {10.1101/2025.03.12.642774},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {A major challenge in scRNAseq analysis is how to recover the biologically meaningful cell ontology tree and conserved gene modules across datasets. Data integration and batch-effect correction have been the key to effectively analyze multiple datasets, but often fail to disentangle cell states in heterogeneous samples, such as in cancer and the immune system. Here we present super single cell clustering (SuperSCC), a novel algorithm that utilizes machine-learning models to discover cell identities and gene modules from multiple datasets without the need of data integration. Of note, SuperSCC can be implemented both in cell lineage and cell state level, thereby building the hierarchy of cell programs with specific cell identity and gene modules. Such information has the great potential to identify the shared rare populations across datasets regardless of batch effect and benefits label transfer for mapping cell labels from reference to query. We used SuperSCC to perform atlas level data analysis on more than 90 datasets and build a cell state map of complex tissue in healthy and diseased stages, such as human lung. We show that SuperSCC outperforms existing approaches in identifying cellular context, has better annotation accuracy, and outlines gene modules that indicate conserved immune cell status in lung microenvironments.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2025/03/14/2025.03.12.642774},
	eprint = {https://www.biorxiv.org/content/early/2025/03/14/2025.03.12.642774.full.pdf},
	journal = {bioRxiv}
}

