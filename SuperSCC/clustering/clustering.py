"""
Clustering module for functions related to data clustering.
"""

import os
import pickle
import warnings
import copy
from datetime import datetime
from collections import Counter
from scipy.spatial.distance import pdist, squareform
from functools import reduce
import re

from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.graph_objects as go


from ..utils import record_time, key_search, retrieve_positive_markers
from ..feature_selection import pairwise_compare, find_signature_genes, find_markers_ovr


def core_consensus_cluster(data, 
                           n_components = 30,
                           resolution = 1,
                           return_hvgs = False):
    """
    A function to do KNN clustering via using Scanpy API.

    Parameters
    ----------
    data: 
        A pandas data frame. Rows are cells; columns are features.
    n_components: 
        A int to decide how many principle components will be used for KNN clustering.
    resolution: 
        A int to control the coarseness of the clustering. Higher values lead to more clusters.
    return_hvgs: 
        A Bool value to decide whether highly variable genes detected by Scanpy will be returned. Default is False. 
    """
    adata = sc.AnnData(data.copy())
    sc.pp.highly_variable_genes(adata, n_top_genes = 3000)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs = n_components)
    sc.tl.leiden(adata, n_iterations = 2, resolution = resolution)


    if return_hvgs:
        return list(adata.obs.leiden), list(adata.var.loc[list(adata.var.highly_variable), :].index)
    else:    
        return list(adata.obs.leiden)


def global_consensus_cluster(       
                      data, 
                      n_components = 30,
                      resolution = 1,
                      only_positive = True,
                      ratio_of_none_zero_counts = 0.1,
                      n_features_to_select = 0.15,
                      cut_off = 0.1,
                      focus = "expression1",
                      class_weight = "balanced",
                      ep_cut_off = 1,
                      pct_cut_off = 0.5, 
                      robust = True,
                      save = True,
                      file_name = None,
                      n_jobs = -1,
                      logger = None,
                      **kwargs
                     ):
    
    """
    A function to merge global clusters and find markers for each global cluster.

    Parameters
    ----------
    data: 
        A log normalized counts matrix. Rows are cells; Columns are features. 
    n_components: 
        A int to decide how many principle components will be used for KNN clustering. Default is 30.
    resolution: 
        A int to control the coarseness of the clustering. Higher values lead to more clusters. Default is 1.
    only_positive: 
        A Bool value to control whether only positive features are used. Default is True.
    ratio_of_none_zero_counts: 
        A float to determine the cutoff in which a feature will be omited when below specified value. A higher value leads to less features will be kept for calculation. Default is 0.1. 
    n_features_to_select: 
        A int or float to control the number of features to select. If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default is 0.15.
    cut_off: 
        A int or float to control the cutoff of combined score (jaccard index * distance correlation) that will used to decide whether a pair of cluster should be merged or not. Default is 0.1.
    focus: 
        A string to decide the value used for distance correlation calculation. "expression1" - the mean expression of feauture will be used. "rank1" - the rank of feature based on expression will be used. Default is 'expression1'.
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    ep_cut_off: 
        A int or float value to decide the expression cutoff to select positive informative features. Default is 1.
    pct_cut_off: 
        A int or float value to decide the expression percentage cutoff to select positive informative features. Default is 0.5. 
    robust: 
        A Bool value to decider whether re calculate the markers of clusters after merging. Default is True. 
    save: 
        A Bool value to decide whether write the output into the disk. Default is true.
    filename: 
        A string to control the name of output. Default is None. 
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    logger: 
        A log_file object. Default is None.
    **kwargs: 
        Other paremeters passed to feature_selection function. 
    """
    
    global_cluster = core_consensus_cluster(data = data.copy(), n_components = n_components, resolution = resolution)
    logger.write("info", "{} clusters are found under resolution {}".format(len(set(global_cluster)), resolution))
    
    # add the 'cluster' column in the data
    global_cluster = [str(i) for i in global_cluster]
    data.loc[:, "cluster"] = global_cluster
    data.loc[:, "cluster"] = data.loc[:, "cluster"].astype(str)
                              
    # find sinataure genes for whole global cluster
    logger.write("info", "finding overall signature genes")

    global_features = find_signature_genes(data = data.copy(), label_column = "cluster", n_features_to_select = n_features_to_select, ratio_of_none_zero_counts = ratio_of_none_zero_counts,logger = logger, save = False, filename = "global", class_weight = class_weight, n_jobs = n_jobs,**kwargs)
    global_features_intersection = global_features["features"]["final_feature_selection_by_intersection"]
    global_features_ensemble = [i[0] for i in global_features["features"]["final_feature_selection_by_ensemble"]]
    # global_features_combined = set(global_features_intersection).intersection(global_features_ensemble)
    if len(global_features_ensemble) > len(global_features_intersection):
        global_features_combined = global_features_ensemble
    else:
        global_features_combined = global_features_intersection
    target = ["cluster"]
    target.extend(global_features_combined)
    
    # remove "sub_high_expression_genes" key
    # global_features.pop("sub_high_expression_genes")

    logger.write("info", "{} overall signature genes are found".format(len(global_features_combined)))

    copy_data = data.copy()
    copy_data = copy_data.loc[:, copy_data.columns.isin(target)]
   

    global_sign_features_pool = find_markers_ovr(copy_data, label_column = "cluster", n_features_to_select = n_features_to_select, ratio_of_none_zero_counts = ratio_of_none_zero_counts, logger = logger, save = False, class_weight = class_weight, n_jobs = n_jobs, **kwargs) 
    global_sign_features_pool_copy = copy.deepcopy(global_sign_features_pool)
    
    # merge similar cluster
    pre_merge_dict = {"global_sign_features": global_sign_features_pool_copy}
    merge_cluster = between_group_merge(pre_merge_dict, data2 = data.copy(), cut_off = cut_off, focus = focus, ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, only_positive = only_positive, logger = logger)
    merge_keys = list(merge_cluster["global_sign_features"].keys())
    
    data.loc[:, "cluster"] = data.loc[:, "cluster"].astype(str)
    for i in merge_keys:
        if re.search("&", i):
            merge_candidate = re.split("&", i)
            for j in merge_candidate:
                data.loc[data.cluster == j, "cluster"] = str(i)

    # relabel the cluster label afetr merging
    encoder = LabelEncoder().fit(data.cluster)
    relabel_classes = list(encoder.classes_)
    relabel_cluster_labels = encoder.transform(data.cluster)

    if robust:
        data.loc[:, "cluster"] = [str(i) for i in relabel_cluster_labels]
        data.loc[:, "cluster"] = data.loc[:, "cluster"].astype(str)
        global_sign_features_after_merging = find_markers_ovr(data.copy(), label_column = "cluster", n_features_to_select = n_features_to_select, ratio_of_none_zero_counts = ratio_of_none_zero_counts, logger = logger, save = False, class_weight = class_weight, n_jobs = n_jobs, **kwargs) 

        res = {
                "global_features_before_merging": global_features,
                "global_sign_features_before_merging": global_sign_features_pool,
                "global_cluster_before_merging": global_cluster,
                "global_sign_features_after_merging": global_sign_features_after_merging,
                "global_cluster_after_merging": {"label_classes": relabel_classes, "labels": relabel_cluster_labels, "merge_label": [relabel_classes[i] for i in relabel_cluster_labels]}
              } 
    
    else:
        res = {
                "global_features_before_merging": global_features,
                "global_sign_features_before_merging": global_sign_features_pool,
                "global_cluster_before_merging": global_cluster,
                "global_sign_features_after_merging": merge_cluster["global_sign_features"],
                "global_cluster_after_merging": {"label_classes": relabel_classes, "labels": relabel_cluster_labels, "merge_label": [relabel_classes[i] for i in relabel_cluster_labels]}
              } 
    
    if save == True:
        if file_name != None:
            file_name = f"consensus_cluster_result_{file_name}_{record_time()}.pkl"
        else:
            file_name = f"consensus_cluster_result_{record_time()}.pkl"
        
        with open(file_name, "wb") as file:
            pickle.dump(res, file)
    
    logger.write("info", "finish global consensus clustering")
    return res


def within_group_merge(data, ep_cut_off = 1, pct_cut_off = 0.5, num_pos_genes = 10):
    """
    A function to find the mergeing candidates in sub clustering level. 

    Parameters
    ----------
    data: 
        A dict returned by find_signature_genes function. 
    ep_cut_off: 
        A int or float value to decide the expression cutoff to select positive informative features. Default is 1.
    pct_cut_off: 
        A int or float value to decide the expression percentage cutoff to select positive informative features. Default is 0.5. 
    num_pos_genes: 
        A int to decide the cutoff for positive markers. If the number of positive markers in a cluster is above this number, it will not be merged. Otherwise, it would be merged. Default is 10. 
    """

    keep = dict()
    merge = dict()
    
    key1 = data["sub_sign_genes"].keys()
    
    for i in key1:

        keep_ls = list()
        merge_ls = list()

        key2 = data["sub_sign_genes"][i].keys()
        
        for j in key2:
            candidates = retrieve_positive_markers({"sub_sign_genes_before_merging": data["sub_sign_genes"]}, i, j, pct_cut_off = pct_cut_off, ep_cut_off = ep_cut_off, level = "before")

            score_pool = [z[1] for z in data["sub_sign_genes"][i][j]["features"]["final_feature_selection_by_ensemble"]]
            names = [z[0] for z in data["sub_sign_genes"][i][j]["features"]["final_feature_selection_by_ensemble"]]
            score_df = pd.DataFrame({"name": names, "score": score_pool})
            score_p80 = np.percentile(score_df.score, 80)

            candidates = candidates.loc[candidates.score > score_p80]
 
            if candidates.shape[0] >= num_pos_genes:
                keep_ls.append(j)
            else:
                merge_ls.append(j)
        
        keep[i] = keep_ls
        merge[i] = merge_ls

    return keep, merge


def sub_consensus_cluster(data, 
                          data2, 
                          n_components = 30, 
                          resolution = 0.2, 
                          class_weight = "balanced", 
                          n_features_to_select = 0.15,
                          ratio_of_none_zero_counts = 0.1,
                          ep_cut_off = 1,
                          pct_cut_off = 0.5,
                          num_pos_genes = 10,
                          n_neighbors = 100,
                          robust = True,
                          n_jobs = -1,
                          file_name = None,
                          save = True, 
                          logger = None, 
                          **kwargs):

    """
    A function to merge sub clusters in different global clusters and find markers for corresponding sub cluster. 
    
    Parameters
    ----------
    data: 
        A log-normalized expression matrix. Rows are cells; Columns are features.
    data2: 
        A dict returned by running global_consensus_cluster function.
    n_components: 
        A int to decide how many principle components will be used for KNN clustering. Default is 30.
    resolution: 
        A int to control the coarseness of the clustering. Higher values lead to more clusters. Default is 0.2.
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    n_features_to_select:
        A int or float to control the number of features to select. If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default is 0.15. 
    ep_cut_off: 
        A float value to control the selection criteria of positive markers. Default is 1.
    pct_cut_off: 
        A float value to control the selection criteria of positibe markers. Default is 0.5.
    num_pos_genes: 
        A int to decide the cutoff for positive markers. If the number of positive markers in a cluster is above this number, it will not be merged. Default is 10. 
    n_neighbors:
        A int to control the number of neighborhoods being considered for merging. Default is 100. 
    robust: 
        A Bool value to decide whether re find the markers of each cluster after cluster merging. Default is true.
    filename: 
        A string to control the name of output. Default is None. 
    save: 
        A Bool value to decide whether write the output into the disk. Default is true.
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    logger: 
        A log_file object. Default is None.
    **kwargs: 
        Other paremeters passed to feature_selection function. 
    """
    
    # add a cluster column in the data
    data.loc[:, "cluster"] = data2["global_cluster_after_merging"]["labels"]
    
    # group data by global cluster
    group_data = data.groupby("cluster")

    # make containers
    sub_sign_features = dict()
    sub_clusters = dict()
    compare = dict()
    discard = list()
    
    # redo the clustering in each sub population
    for key1 in list(set(data.cluster)):
       
        # get the subset data 
        sub_data = group_data.get_group(key1)
        
        # remove cluster column before clustering
        sub_data = sub_data.loc[:, sub_data.columns != "cluster"]
        
        logger.write("info", "doing sub clustering in global cluster {}".format(key1))  

        # when speciifed n_components > min(n_samples, n_features)
        n_components_update = None
        if n_components > np.min(sub_data.shape):
            n_components_update = np.min(sub_data.shape) - 1
    
        if n_components_update == None:
            sub_cluster = core_consensus_cluster(data = sub_data.copy(), n_components = n_components, resolution = resolution, return_hvgs = True)
        else:
            sub_cluster = core_consensus_cluster(data = sub_data.copy(), n_components = n_components_update, resolution = resolution, return_hvgs = True)
            
        sub_clusters["sub_clusters_in_global_cluster_{}".format(key1)] = sub_cluster[0]
        
        # add the 'cluster' column in the sub data
        sub_data.loc[:, "cluster"] = sub_cluster[0]

        # do PCA 
        filter_data = sub_data.loc[:, sub_data.columns.isin(sub_cluster[1])]

        if n_components_update == None:
            pca = PCA(n_components = n_components).fit_transform(filter_data)
        else:
            pca = PCA(n_components = n_components_update).fit_transform(filter_data)
        
        pca = pd.DataFrame(pca, columns = ["PC_{}".format(i) for i in range(pca.shape[1])])
        pca.loc[:, "cluster"] = sub_cluster[0]

        pca_group = pca.groupby("cluster")

        pca_centroids = list()
        for i in pca_group.groups.keys():
            individual = pca_group.get_group(i)
            individual = pd.DataFrame(individual.select_dtypes("number").mean(), columns=[i]).T
            pca_centroids.append(individual)

        # use PCA centroids to calculate the distance between clusters
        pca_centroids = pd.concat(pca_centroids)
        distance = pdist(pca_centroids, metric = 'euclidean')
        square_distances = pd.DataFrame(squareform(distance), columns = pca_centroids.index)
        square_distances.index = pca_centroids.index

        
        for i in range(square_distances.shape[0]):
            index = square_distances.iloc[i, :].sort_values()[1:(n_neighbors+1)].index
            key = square_distances.index[i]
            if i == 0:
                compare[str(key1)] = {str(key): list(index)}
            else:
                compare[str(key1)].update({str(key): list(index)})

        # find sinataure genes for each sub cluster within corresponding global cluster
        sub_sign_features_pool = dict()
        
        for key2 in list(set([str(i) for i in sub_cluster[0]])):
            logger.write("info", "finding signature genes in sub cluster {} of global cluster {}".format(key2, key1))
            intermediate = sub_data.copy()
            # intermediate = intermediate.loc[:, intermediate.columns.isin(target)]
            intermediate.loc[:, "cluster"] = intermediate.loc[:, "cluster"].astype(str)
            selection = compare[str(key1)][key2]
            selection.append(key2)
            intermediate = intermediate.loc[intermediate.cluster.isin(selection), :] # only do One vs Rest for the nearby clusters, aka One vs Neighbors
            intermediate.loc[intermediate.cluster != str(key2), "cluster"] = "OnevsRest" 
            if len(set(intermediate.cluster)) >= 2:
                intermediate_features = find_signature_genes(data = intermediate, label_column = "cluster", HVG = False, logger = logger, save = False, filename = "sub_cluster_{}".format(key), n_features_to_select = n_features_to_select, ratio_of_none_zero_counts = ratio_of_none_zero_counts, class_weight = class_weight, n_jobs = n_jobs, **kwargs)
                sub_sign_features_pool[str(key2)] = intermediate_features
            else:
                sub_sign_features_pool[str(key2)] = f"no more clusters are found in sub cluster {key2} global cluster {key1} before cluster merging"
                discard.append(f"{key1}_{key2}")
        
        sub_sign_features[str(key1)] = sub_sign_features_pool
    
    robust_d = dict()  
    if robust:
        # generate new label 
        logger.write("info", "re find signature genes of clusters after cluster merging")
        copy_data = data.copy()
        copy_data.loc[:, "row"] = range(copy_data.shape[0])
        
        # remove sub cluster with only one member 
        sub_sign_features_copy = copy.deepcopy(sub_sign_features)
        for i in discard:
            key1 = re.split("_", i)[0]
            key2 = re.split("_", i)[1]
            sub_sign_features_copy[key1].pop(key2)

        merge_sub = within_group_merge({"sub_sign_genes": sub_sign_features_copy}, num_pos_genes = num_pos_genes, ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off)
   

        ls = list()
        for i in sub_clusters.keys():
            key = re.findall("\d+", i)[0]
            sub_data = copy_data.loc[copy_data.cluster == int(key), :]
            sub_data.loc[:, "sub_cluster"] = [f"{key}_{z}" for z in sub_clusters[i]]
            condition1 = list()
            condition2 = list()
            if len(merge_sub[1][key]) >= 2:
                for j in merge_sub[1][key]:
                    condition1.append(f"{key}_{j}")
                    condition2.append(j)
            sub_data.loc[sub_data.sub_cluster.isin(condition1), "sub_cluster"] = f"{key}_merge_{'-'.join(condition2)}"
            ls.append(sub_data)
        
        # label encoding
        new_label = pd.concat(ls).sort_values("row")
        le = LabelEncoder().fit(new_label.loc[:, "sub_cluster"])
        new_label = le.transform(new_label.loc[:, "sub_cluster"])
        
        copy_data.loc[:, "final_cluster"] = new_label

    
        # find signature genes of clusters after cluster merging
        for i in copy_data.cluster.unique():
            wk_data = copy_data.loc[copy_data.cluster == i]
            wk_data = wk_data.loc[:, wk_data.columns != "cluster"]
            wk_data = wk_data.loc[:, wk_data.columns != "row"]
            if (wk_data.shape[1] - 1) != (data.shape[1] - 1):
                raise Exception("Extra none-useful columns in the dataset. Please double check the input data only contains feature columns representating expression in each cell")
            
            wk_data.loc[:, "final_cluster"] = le.inverse_transform(wk_data.loc[:, "final_cluster"])
            if len(set(wk_data.loc[:, "final_cluster"])) >= 2:
                robust_d[str(i)] = find_markers_ovr(wk_data, label_column = "final_cluster", n_features_to_select = n_features_to_select, ratio_of_none_zero_counts = ratio_of_none_zero_counts, logger = logger, save = False, class_weight = class_weight, n_jobs = n_jobs, **kwargs) 
            else:
                robust_d[str(i)] = f"no more sub clusters in global cluster {i} after cluster merging"

    if robust:
        res = {"sub_cluster": sub_clusters, "sub_sign_genes_after_merging": robust_d, "sub_sign_genes_before_merging": sub_sign_features, "neighbors": compare, "merge_labels": new_label, "merge_candidates": merge_sub, "label_encoder": le}
    else:    
        res = {"sub_cluster": sub_clusters, "sub_sign_genes": sub_sign_features, "neighbors": compare}

    if save == True:
        if file_name != None:
            file_name = f"sub_consensus_cluster_result_{file_name}_{record_time()}.pkl"
        else:
            file_name = f"sub_consensus_cluster_result_{record_time()}.pkl"
        
        with open(file_name, "wb") as file:
            pickle.dump(res, file)
    
    logger.write("info", "finish sub consensus cluster")
    
    return res


def between_group_merge(data, data2, cut_off = 0.1, increase = 0, focus = "expression1", n_merge = 100, ep_cut_off = 1, pct_cut_off = 0.5, only_positive = True, logger = None):

   """
   A function to merge clusters based on features' similarity evulated by jaccard index and distance correlation.
   
   Parameters
   ----------
   data: 
        A dict return by running find_signature_genes function.
   data2: 
        A count expression matrix. Rows are cells; Columns are features. 
   cut_off: 
        A int or float to control the cutoff of combined score (jaccard index * distance correlation) that will used to decide whether a pair of cluster should be merged or not. Default is 0.1.
   penalty_increase: 
        A float to control the augment of cut_off in the next round when the function is recursively called. 
   focus: 
        A string to decide the value used for distance correlation calculation. "expression1" - the mean expression of feauture in the interest group will be used. "rank1" - the rank of feature based on expression in the interest group will be used. Default is 'expression1'.
   n_merge:
        A int value to control the number of candidates which should be merged in each round. Default is 100.
   ep_cut_off: 
        A int or float value to decide the expression cutoff to select positive informative features. Default is 1.
   pct_cut_off: 
        A int or float value to decide the expression percentage cutoff to select positive informative features. Default is 0.5. 
   only_positive: 
        A Bool value to control whether only positive features are used. Default is True. 
   logger: 
        A log_file object. Default is None. a
   """
      
   compare_data = pairwise_compare(data, focus = focus, only_positive = only_positive, ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off)
   try:
      distance_cor = [*zip(compare_data["distance_cor"].keys(), compare_data["distance_cor"].values())]
      jaccard_score = [*zip(compare_data["jaccard_index"].keys(), compare_data["jaccard_index"].values())]
   except:
      return compare_data
      
   # times the distance correlation and jaccard score and only keep the score above the cutoff 
   combined_score = list()
   for i in range(len(distance_cor)):
      score = distance_cor[i][1] * jaccard_score[i][1]
      if score > cut_off:
         print(score)
         if distance_cor[i][0] != jaccard_score[i][0]:
            raise ValueError("The key between distance_cor and jaccard_index is not consistent")
         else:
            combined_score.append((distance_cor[i][0], score))
      else:
         continue

   if len(combined_score) != 0:
      # detect the label tie
      # e.g 0vs1, 0vs10, 1vs10
      combined_score_label = [ i[0] for i in combined_score ]
      
      former_num = set(sorted([*map(detect_former_num, combined_score_label)]))
      label_tie = dict()
      for i in former_num:
         ls = []
         for j in combined_score_label:
            if re.search("global cluster {} vs global cluster \d+".format(i), j) and len(re.findall("\d+", j)[0]) == len(i):
                  ls.append(j)
            elif re.search("global cluster \d+ vs global cluster {}".format(i), j) and len(re.findall("\d+", j)[1]) == len(i):
               ls.append(j)
            else:
               target = find_keys(j)
               target1 = target[0]
               target2 = target[1]
               if target1 == i or target2 == i:
                  ls.append(j)
            
         label_tie[i] = ls
      
      # combine the two groups with label tie together  
      new_d = key_search(label_tie)
      
      # remove individual keys that have been combined in the new_d
      # from label_tie
      discard_keys = list()
      condition = re.compile("((\d+&?){1,100}|(\d+))")
      for i in new_d.keys():
         target = [*condition.finditer(i)]
         for j in range(len(target)):
            target1 = i[target[j].start():target[j].end()]
            discard_keys.append(target1)
      for i in discard_keys:
         label_tie.pop(i)
      
      # combine label_tie and new_d
      label_tie.update(new_d)
      
      # merge the group with the highest score in each pair group 
      merge_candidate = list()
      for key in label_tie.keys():
         label = label_tie[key]
         merge_candidate.append(sorted([*filter(lambda x: x[0] in label, combined_score)], key = lambda x: -x[1]))
      
      pool = list()
      final_candidate = list()
      num = 0
      for i in merge_candidate:
         for j in i:
            labels = find_keys(j[0])
            if labels[0] not in pool and labels[1] not in pool:
               pool.append(labels[0])
               pool.append(labels[1])
               if num < n_merge:
                  final_candidate.append(j[0])
                  num += 1
               else:
                  break
         num = 0
         pool = list()

   
      high_expression = data["global_sign_features"]
      sign_genes = data["global_sign_features"]

      # for i in merge_candidate:
      for i in final_candidate:
         
         cluster_label = find_keys(i)
         compare1 = cluster_label[0]
         compare2 = cluster_label[1]

         print("compare1: {}; compare2: {}".format(compare1, compare2))
         
         if isinstance(sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"][0], tuple) and isinstance(sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"][0], tuple):
            geneset1 = [ j[0] for j in sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"] ]
            geneset2 = [ j[0] for j in sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"] ]
            union_geneset = list(set(geneset1).union(geneset2))

         elif isinstance(sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"][0], tuple) and isinstance(sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"][0], tuple) == False:
            geneset1 = [ j[0] for j in sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"] ]
            geneset2 = sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"] 
            union_geneset = list(set(geneset1).union(geneset2))

         elif isinstance(sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"][0], tuple) == False and isinstance(sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"][0], tuple):
            geneset1 = sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"]
            geneset2 = [ j[0] for j in sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"] ]
            union_geneset = list(set(geneset1).union(geneset2))

         else:
            geneset1 = sign_genes[compare1]["features"]["final_feature_selection_by_ensemble"]
            geneset2 = sign_genes[compare2]["features"]["final_feature_selection_by_ensemble"] 
            union_geneset = list(set(geneset1).union(geneset2))
            
         high_expression_dataset1 = high_expression[compare1]["sub_high_expression_genes"][compare1]
         high_expression_dataset2 = high_expression[compare2]["sub_high_expression_genes"][compare2]
         union_high_exp_dataset = pd.concat([high_expression_dataset1, high_expression_dataset2])
         union_high_exp_dataset = union_high_exp_dataset.sort_values("feature")
         union_high_exp_dataset.index = range(union_high_exp_dataset.shape[0])
         
         duplicate_feature_rows = np.ravel(np.array(np.where(union_high_exp_dataset.duplicated("feature", keep=False) == True)))
         
         # when shared features found, calucate the mean expression and rank between two datasets
         if len(duplicate_feature_rows) != 0:
            uniq_high_exp_dataset = union_high_exp_dataset.loc[union_high_exp_dataset.index.isin(duplicate_feature_rows) == False, :]
            dup_high_exp_dataset = union_high_exp_dataset.loc[union_high_exp_dataset.index.isin(duplicate_feature_rows) == True, :].sort_values("feature")
            dup_features = list(dup_high_exp_dataset.drop_duplicates("feature").feature)
            dup_high_exp_dataset = dup_high_exp_dataset.select_dtypes("number").rolling(window = 2).mean()[1::2]
            dup_high_exp_dataset.loc[:, "feature"] = dup_features # add feature column back
            try:
               dup_high_exp_dataset = dup_high_exp_dataset[["feature", "expression1", "pct1", "rank1", "expression2", "pct2", "rank2"]] # reorder the column
            except:
               dup_high_exp_dataset = dup_high_exp_dataset[["feature", "expression1", "pct1", "rank1"]]
               
            final_high_exp_dataset = pd.concat([uniq_high_exp_dataset, dup_high_exp_dataset])
         else:
            final_high_exp_dataset = union_high_exp_dataset
         
         # add combined keys
         if logger != None: 
            logger.write("info", "{} merge {} and {}".format(record_time(), compare1, compare2))
         print("{} merge {} and {}".format(record_time(), compare1, compare2))
         
         data["global_sign_features"].update({"{}&{}".format(compare1, compare2): {"features": {"final_feature_selection_by_ensemble": union_geneset}}})
         data["global_sign_features"]["{}&{}".format(compare1, compare2)].update({"sub_high_expression_genes": {"{}&{}".format(compare1, compare2): final_high_exp_dataset}})

         # remove keys that has been merged    
         data["global_sign_features"].pop(compare1)
         data["global_sign_features"].pop(compare2)

      # re calculate the markers of each cluster after merging 
      merge_keys = list(data["global_sign_features"].keys())
      data2_copy = data2.copy()
      data2_copy.loc[:, "cluster"] = data2_copy.loc[:, "cluster"].astype(str)
      for i in merge_keys:
         if re.search("&", i):
            merge_candidate = re.split("&", i)
            for j in merge_candidate:
                  data2_copy.loc[data2_copy.cluster == j, "cluster"] = str(i)
         
      new_markers = find_markers_ovr(data2_copy, 
                                       label_column = "cluster", 
                                       n_features_to_select = 0.15, 
                                       ratio_of_none_zero_counts = 0.1, 
                                       logger = logger, 
                                       save = False, 
                                       class_weight = "balanced", 
                                       n_jobs = -1, 
                                       variance_threshold = "mean",
                                       model = "svm",
                                       mutual_info = False,
                                       F_test = True,
                                       normalization_method = "Min-Max")
      
      data = {"global_sign_features": new_markers}

      
      return between_group_merge(data, data2 = data2, cut_off = cut_off + increase, increase = increase, n_merge = n_merge, focus = focus, ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, logger = logger)
   
   else:
      return data


def find_keys(data, split_by = "vs"):
    """
    A function to split string by specified symbol.

    Parameters
    ----------
    data: 
        A string object.
    split_by: 
        A string to decide the regular expression to split the data.
    """
    
    res = re.split(split_by, data)
    output = list()
    for i in res:
        try:
            reg = [*re.compile("((\d+&?){1,100}|(\d+))").finditer(i)][0]
        except:
            raise Exception("Can't find numeric key")
        output.append(i[reg.start():reg.end()])
    return output


def detect_former_num(data):
    """
    A function to find the first numeric elements in a string.
    """
    res = re.split("vs", data)
    output = [*re.compile("((\d+&?){1,100}|(\d+))").finditer(res[0])][0]
    return res[0][output.start():output.end()]


def get_sankey_dataframe(data, plot = False):
    """
    A function to generate the dataframe for sankey plotting.
    
    Parameters
    ----------
    data: 
        A pandas dataframe. From first column to last column, it must conver the original source and final target in order. 
    plot: 
        A Bool value to decide whether plotting basic sankey diagram. 
    """

    # encode source and target labels
    colnames = list(data.columns)
    labels_pool = list()
    for i in range(data.shape[1]):
        data.iloc[:, i] = [f"{colnames[i]}_{j}" for j in data.iloc[:, i]]
        labels_pool.extend(list(set(data.iloc[:, i])))
    
    le = LabelEncoder().fit(labels_pool)

    for i in range(data.shape[1]):
        data.iloc[:, i] = le.transform(data.iloc[:, i])


    # for each source, calculate the numer of occasions to different target  
    source_target_num = data.shape[1]
    target = dict()
    for i in range(source_target_num-1):
        start = set(data.iloc[:, i])
        column = colnames[i]
        for j in start:
            sub_data = data.query("{} == @j".format(column))
            target_num = Counter(sub_data.iloc[:, i+1])
            for z in target_num.keys():
                target["{}_{}_{}".format(i, j, z)] = target_num.get(z)
    

    # detect whether one source has different targets
    # if yes, record the number of targets and target names
    replicate = list()
    for j in target.keys():
        new_source_key = re.split("_", j)[1]
        column = re.split("_", j)[0]
        key = f'{column}_{new_source_key}'
        replicate.append(key)
    
    replicate_count = dict()
    for i in replicate:
        count = 0
        member = list()
        for j in target.keys():
            if re.search(f"{i}_\d+", j):
              count += 1
              member.append(re.split("_", j)[2])
        replicate_count[f"{i}"] = count, member

    # generate the output sankey dataframe
    new_source = list()
    new_target = list()
    value = list()

    for i in replicate_count.keys():
        replicate_num = replicate_count[i][0]
        replicate_member = replicate_count[i][1]

        if replicate_num == 1:
            new_target.append(replicate_member[0])
            new_source.append(re.split("_", i)[1])
            value.append(target[f"{i}_{replicate_member[0]}"])
        else:
            new_target.extend(replicate_member)
            new_source.extend([re.split("_", i)[1]] * replicate_num)
            for j in replicate_member:
                value.append(target[f'{i}_{j}'])

    sankey_df = pd.DataFrame({"source": new_source, "target": new_target, "value": value})
 
    sankey_df = sankey_df.astype(int)
    labels = le.inverse_transform(np.unique(sankey_df[["source", "target"]].values))

    # plotting
    if plot:

        fig = go.Figure(data=[go.Sankey(
            node = dict(
            pad = 15,
            thickness = 20,
            line = dict(color = "black", width = 0.5),
            label = labels
            ),
            link = dict(
            source = sankey_df.source, 
            target = sankey_df.target,
            value =  sankey_df.value
            ))])

        fig.update_layout(title_text="Basic Sankey Diagram", font_size=10,width=1200, height=800)
        fig.show()
    
    return sankey_df, labels
