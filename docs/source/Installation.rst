Installation 
======================

We encourage to create a new virtual environment before installing SuperSCC::
        
        conda create -n superscc_env
        conda activate superscc_env

To install SuperSCC, you should run the following code if you already get `SuperSCC's tar ball <https://github.com/tf1993614/SuperSCC/tree/main/dist>`_::

        pip install SuperSCC.tar.gz

Usually, all dependencies should be downloaded and installed automatically. 

To install SuperSCC via GitHub, you can do::

        git clone https://github.com/tf1993614/SuperSCC/
        cd SuperSCC
        python setup.py install

Since the `get_gene_module` function was written in R, you also need install `geneModule <https://github.com/tf1993614/SuperSCC/tree/main/geneModule>`_ R package::

        R
        install.packages("geneModule_0.1.tar.gz", repos = .libPaths()[1], type = "source")
        .libPaths()[1] # get the location where packages are installed
        q()

Note: `get_gene_module` function can be also implemented in Python environment via rpy2-wrappered function.

After installing correctly, you can load SuperSCC::

        import SuperSCC as scc
 

