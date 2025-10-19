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

After installing correctly, you can load SuperSCC::

        import SuperSCC as scc

Bonus for R person
=======================

Since the gene module-relevant functions was also written in R, you can also do the gene module analysis in R environment aftering installing [geneModule](https://github.com/tf1993614/SuperSCC/tree/main/geneModule) R package::

        R
        install.packages("geneModule_0.1.tar.gz", repos = .libPaths()[1], type = "source")
        .libPaths()[1] # get the location where packages are installed
        q()
 

