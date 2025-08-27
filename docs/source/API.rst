.. automodule:: SuperSCC

API
==============

Import SuperSCC as::
        
        import SuperSCC

Feature selection
~~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
  :toctree: generated

  SuperSCC.feature_selection.feature_selection
  SuperSCC.feature_selection.find_signature_genes
  SuperSCC.feature_selection.find_markers_ovr

Clustering
~~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.clustering.global_consensus_cluster
   SuperSCC.clustering.sub_consensus_cluster

Label transfer
~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.label_transfer.model_training
   SuperSCC.label_transfer.predict_label

Gene module
~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.gene_module.get_gene_module
   SuperSCC.gene_module.compare_gene_modules
   SuperSCC.gene_module.analyse_one_gene_module

Retrival augment generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
  :toctree: generated

  SuperSCC.rag.SimpleRAG
  SuperSCC.rag.SimpleRAG.run_rag
  SuperSCC.rag.SimpleRAG.get_answer
  SuperSCC.rag.SimpleRAG.update_rag_chain
  
  SuperSCC.rag.ConnectRAG
  SuperSCC.rag.ConnectRAG.connect_client
  SuperSCC.rag.ConnectRAG.run_rag 


