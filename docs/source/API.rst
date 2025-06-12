.. automodule:: SuperSCC

API
==============

Import SuperSCC as::
        
        import SuperSCC as scc

Feature selection
~~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
  :toctree: generated

  SuperSCC.feature_selection
  SuperSCC.find_signature_genes
  SuperSCC.find_markers_ovr

Clustering
~~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.global_consensus_cluster
   SuperSCC.sub_consensus_cluster

Label transfer
~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.model_training
   SuperSCC.predict_label

Gene module
~~~~~~~~~~~~~~~~~
.. currentmodule:: SuperSCC

.. autosummary::
   :toctree: generated

   SuperSCC.get_gene_module
   SuperSCC.compare_gene_modules
   SuperSCC.analyse_one_gene_module

Retrival augment generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. module:: SuperSCC.rag
.. currentmodule:: SuperSCC

.. autosummary::
  :toctree: generated

  rag.SimpleRAG
  rag.SimpleRAG.run_rag
  rag.SimpleRAG.get_answer
  rag.SimpleRAG.update_rag_chain
  
  rag.ConnectRAG
  rag.ConnectRAG.connect_client
  rag.ConnectRAG.run_rag 


