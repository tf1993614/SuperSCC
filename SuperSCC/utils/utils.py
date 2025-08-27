"""
Utils module for generic helper and utility functions used by other modules.
"""

import re
import os
import pickle
import warnings
import logging
from datetime import datetime
from dcor import distance_correlation
import scanpy as sc
from scanpy import AnnData
import pandas as pd
import numpy as np


def jaccard_score(x, y):
    """
    A function to calculate jaccard score.
    
    Parameters
    ----------
    x: 
        A list object.
    y: 
        A list object. 
    """
    intersection = len(set(x).intersection(set(y)))
    union = len(set(x).union(set(y)))
    return intersection / union


def list_files(path, pattern, recursive = True, full_name = True):
    """
    A function to produce a list of the names of files in the named directory.

    Parameters
    -----------
    path: 
        A string to indicate the working directory where the files will be recursively found.
    pattern: 
        A regular expression to match the name of wanted files. 
    recursive: 
        A Bool value to decide whether to recursive into the directories. Default is True.
    full_name: 
        A Bool value to decide whether the directory path is prepended to the file names. Default is True.
    """
    output = []
    def list_files_core(current_path = path, current_pattern = pattern, current_recursive = recursive, current_full_name = full_name):
        nonlocal output
        files = os.listdir(current_path)
        for file in files:
            file_path = os.path.join(current_path, file)
            
            if os.path.isdir(file_path) and current_recursive:
                list_files_core(file_path, current_pattern, current_recursive, current_full_name)
            
            else:
                if re.search(current_pattern, file):
                    if full_name == True:
                        file = os.path.join(current_path, file)
                        output.append(file)
                    else:
                        output.append(file)
    list_files_core()
    return output


def convert_data(data = None, path = None, assay = ".X", label_column = None):
    """
    A function to convert sparse matrix into pandas data frame object. 
    
    Parameters
    ----------
    data: 
        a Anno Data object. If speciefied, path should be specified in 'None'.
    path: 
        path to h5ad file. If speciefied, data should be specified in 'None'.
    assay: 
        ".X" or "raw" assay in Annot Data object could be specified
    label_column: 
        the name of cell type column in the h5ad file. If specified, the cell type column will be added into output.
    """
    if path != None:
        data = sc.read(path)
    
    if assay == ".X":
        try:
            counts = pd.DataFrame.sparse.from_spmatrix(data.X)
        except:
            counts = pd.DataFrame(np.array(data.X))
        finally:
            counts = pd.DataFrame(np.array(counts))
        features = data.var_names.tolist()
        index  = data.obs_names.tolist()
    else:
        try:
            counts = pd.DataFrame.sparse.from_spmatrix(data.raw.X)
        except:
            counts = pd.DataFrame(np.array(data.raw.X))
        finally:
            counts = pd.DataFrame(np.array(counts))
        
        features = data.var_names.tolist()
        index  = data.obs_names.tolist()
   
    counts.columns = features
    counts.index = index
    
    if label_column != None:
        try:
            labels = data.obs[label_column].tolist()
            counts["cell_type"] = labels
        except:
            raise ValueError("The length of cell type column is not consistent with matrix")

    return counts


def quantile_normalize(data):
    """
    A function to do quantile normalization.  
    """
    data = data.loc[:, data.columns != "cell_type"]
    ranks = (data.rank(method = "first").stack())
    rank_mean = (data.stack().groupby(ranks).mean())
    # add interproblated values in between ranks
    finer_ranks = ((rank_mean.index + 0.5).tolist() + rank_mean.index.tolist())
    rank_mean = rank_mean.reindex(finer_ranks).sort_index().interpolate()
    data = data.rank(method = "average").stack().map(rank_mean).unstack()
    
    return data


def record_time():
    """
    A function to call out current time.
    """
    current_second_time = datetime.now()
    return current_second_time.strftime("%Y-%m-%d %H:%M:%S")


def remove_time(string, count = 1):
    """
    A function to remove time information in the stdout message"
    """
    pattern = re.compile("\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")
    res = re.sub(pattern, "", string, count = count)
    return res.rstrip().lstrip()


def load_pickle_file(file_name, wk_dir = os.getcwd(), recursive = True):
    """
    A function to find pickle file and then load it. 

    Parameters
    ----------
    file_name: 
        A regular expression string matching the target files.
    wk_dir:
         A string to specify the directory where target files should be searched.
    recursive:
        A Bool value to control whether recursive search. Default is True.
    """
    current_dir = os.getcwd()
    res = dict()
    try:
        os.chdir(wk_dir)
        file_paths = list_files(os.getcwd(), pattern = file_name, recursive = recursive)
        for path in file_paths:
            with open(path, "rb") as f:
                pickle_file = pickle.load(f)
                res[os.path.basename(path)] = pickle_file
    except:
        raise Exception("can't find the target file or the target file is broken.")
    finally:
        os.chdir(current_dir)
        
    return res


class log_file:
    """
    A class to easily write log information into different log files.
    
    Parameters
    ----------
    filename: 
        A string to indicate the log file name.
    mode: 
        A string to choose the log file mode. "a" means 'append log information'; "w" means 'clean up old log information and then write new information'.
    """
    n_class = 0
    def __init__(self, filename, mode):
        self.time = record_time()
        self.filename = filename
        self.mode = mode
        self.logger = None
        self.n_object = 0
        log_file.n_class += 1
        
    def set(self):
        logger = logging.getLogger("my_logger")
        logger.setLevel(logging.INFO)
        log_file_name = "{}_{}.log".format(self.filename, self.time)
        file_handler = logging.FileHandler(log_file_name, mode = self.mode)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        self.logger = logger
        
    def start(self):
        if log_file.n_class == 1 and self.n_object == 0:
            self.set()
        if log_file.n_class > 1 and self.n_object == 0:
            self.set()
        
    def write(self, category, content):
        self.start()
        self.n_object += 1
        if category == "info":
            self.logger.info(content)
        elif category == "error":
            self.logger.error(content)
        elif category == "warning":
            self.logger.warning(content)
        elif category == "debug":
            self.logger.debug(content)
        elif category == "critical":
            self.logger.critical(content)
    
    def mode_reset(self, mode):
        self.mode = mode
        
    def clear():
        log_file.n_class = 0


def pre_processing(path, 
                   label_column = None,
                   assay = ".X",
                   convert = True,
                   log_normalize = True, 
                   scale_data = False, 
                   quantile_normalize = False,
                   save = False,
                   logger = None):
    """
    A function to do preliminary normalization, encoding label and convert data into pandas data frame for downstream sklearn-based
    machine learning workflow.
    
    Parameters
    ----------
    path: 
        Path to the h5ad file or a AnnData object.
    label_column: 
        The name of cell type column in the h5ad file. If specified, the cell type column will be added into output.
    assay: 
        ".X" or "raw" assay in Annot Data object could be specified. Default is ".X".
    convert: 
        A Bool value to decide whether convert Annot Data object into pandas data frame object. Default is True.
    log_normalize: 
        A Bool value to decide whether standard log normalization to be done. Default is True.
    scale_data: 
        A Bool value to decide whether standadlize data or not. Default is False.
    quantile_normalize: 
        A Bool value to decide whether quantile normalize data or not. Default is False.
    save: 
        A Bool value to decide whether write the pre-processed data into the disk. Default is False.
    logger: 
        A log_file object to write log information into disk. Default is None. Default is None.
    """
    if logger != None:
        logger.write("info", "start pre processing")
        logger.write("critical", "Parameters used for pre_processing")
        parameters = {"path": path,
                      "label_column": label_column,
                      "assay": assay,
                      "convert": convert,
                      "log_normalize": log_normalize,
                      "scale_data": scale_data,
                      "quantile_normalize": quantile_normalize,
                      "save": save
                     }
        for key,value in parameters.items():
            logger.write("critical", "{}: {}".format(key, value))
    
    if isinstance(path, str):
        data = sc.read(path)
    elif isinstance(path, AnnData):
        data = path.copy()
    else:
        data = sc.AnnData(path.copy())
        
    if log_normalize:
        sc.pp.normalize_total(data, target_sum = 1e4)
        sc.pp.log1p(data)
    
    if scale_data:
        sc.pp.scale(data)
        
    if convert:
        counts = convert_data(data = data, assay = assay, label_column = label_column)
    else:
        counts = data.X
        
    if quantile_normalize:
        quant_norm_data = convert_data(data, assay = assay)
        counts = quantile_normalize(quant_norm_data)        
    
    if save:
        file_name = "Preprocessing_data_{}.pkl".format(record_time())
        with open(file_name, "wb") as output:
            pickle.dump(counts, output)
    
    if logger != None: 
        logger.write("info", "finish pre processing")
    
    return counts


def check_duplicated_elements(x):
    """
    A function to fix duplicated string elements in the list by adding ordinal suffix.

    Paremeters
    ----------
    x: 
        A list object containing string elements.
    """
    
    assert isinstance(x, list), "x should be list."
    assert all(map(lambda x: isinstance(x, str), x)), "elements in the list should be string."
    
    if len(x) != len(set(x)):
        warnings.warn("Duplicated elements exist. To fix this, duplicated elements will be added suffix such as _1, _2..._n")
        duplicate = []
        unique = []
        index = []
        
        for i in x:
            if i not in unique:
                unique.append(i)
            else:
                duplicate.append(i)
        
        d = dict()
        for i in set(duplicate):
            for j in range(len(x)):
                if x[j] == i:
                    index.append(j)
            d[i] = index
            index = []
            
        for key, value in d.items():
            for i in range(len(value)):
                replacement = key + "_" + str(i+1)
                x[value[i]] = replacement
        
    return x


def flatten_list(data):
    """
    A function to flatten list.
    """
    res = []
    
    def core_flatten_list(data):
      nonlocal res
      for i in data:
        if isinstance(i, list):
          core_flatten_list(i)
        else:
          res.append(i)
    
    core_flatten_list(data)
    
    return res


def key_search(data):

    """
    A function to find the tie between different keys and values.

    Parameters
    ----------
    data: 
        A dict containing the label tie information. 
    """
    
    new_d = dict()
    condition = re.compile("((\d+&?){1,100}|(\d+))")
   
    # check whether there are shared values betweeen different keys 
    # if yes, combined keys and corresponding values
    keys = list(data.keys())
    for i in range(len(keys)-1):
        for j in range(i+1, len(keys)):
            set1 = set(data[keys[i]])
            set2 = set(data[keys[j]])
            if len(set1.intersection(set2)) != 0:
                new_d["{}_{}".format(keys[i], keys[j])] = set(flatten_list([data[keys[i]], data[keys[j]]]))
            else:
                group1 = [condition.finditer(x) for x in list(set1)]
                group2 = [condition.finditer(x) for x in list(set2)]

                groupA = list()
                for index in range(len(group1)):
                    intermediate = [*group1[index]]
                    for index2 in range(len(intermediate)):
                        match = intermediate[index2]
                        groupA.append(list(set1)[index][match.start():match.end()])

                groupB = list()
                for index in range(len(group2)):
                    intermediate = [*group2[index]]
                    for index2 in range(len(intermediate)):
                        match = intermediate[index2]
                        groupB.append(list(set2)[index][match.start():match.end()])
                
                if len(set(groupA).intersection(groupB)) != 0:
                    new_d["{}_{}".format(keys[i], keys[j])] = set(flatten_list([data[keys[i]], data[keys[j]]]))
    
    signal = False

    # check whether there are shared keys
    # if yes, combined keys and corresponding values

    for _ in range(10000): # first 'for' loop to re initialize second 'for' loop when length of new_d has been changed
        for i in range(len(new_d.keys()) - 1): # second 'for' loop
            new_d_keys = list(new_d.keys())
            if not signal:
                for j in range(i+1, len(new_d_keys)):
                    res1 = [*condition.finditer(new_d_keys[i])]
                    res2 = [*condition.finditer(new_d_keys[j])]

                    res3 = list()    
                    for z in range(len(res1)):
                        res3.append(new_d_keys[i][res1[z].start():res1[z].end()])

                    res4 = list()    
                    for z in range(len(res2)):
                        res4.append(new_d_keys[j][res2[z].start():res2[z].end()])

                    
                    if len(set(res3).intersection(res4)) != 0:
                        new_d.update({"_".join(list(set(res3).union(res4))) : new_d[new_d_keys[i]].union(new_d[new_d_keys[j]])})
                        # when combined key is not equivalent to either res3 or res4.
                        if len(set(res3).union(res4)) != len(res3) and len(set(res3).union(res4)) != len(res4):
                            new_d.pop(new_d_keys[i])
                            new_d.pop(new_d_keys[j])
                        # when combined key is equivalent to res3.
                        elif len(set(res3).union(res4)) == len(res3):
                             new_d.pop(new_d_keys[j])
                        # when combined key is equivalent to res4. 
                        elif len(set(res3).union(res4)) == len(res4):
                             new_d.pop(new_d_keys[i])
                        
                        signal = True # terminate the second 'for' loop since the length of new_d is changing
                        break
            else:
                signal = False
                break
            
    return new_d


def retrieve_positive_markers(data, m_key, f_key = None, pct_cut_off = 0.5, ep_cut_off = 1, level = "after"):
    """
    A function to retrieve positive markers of clusters after global clustering or sub clustering.

    Parameters
    ----------
    data: 
        A dict return by global_consensus_cluster or sub_consensus_cluster function.
    m_key: 
        A string to control selecting the global cluster. 
    f_key: 
        A string to control selecting the sub cluster. Default is None. When f_key is set, m_key should be set as the origin gloabl cluster of the selected sub cluster.
    ep_cut_off: 
        A int or float value to decide the expression cutoff to select positive informative features.
    pct_cut_off: 
        A int or float value to decide the expression percentage cutoff to select positive informative features.
    level: 
        A string to control selecting positive markers before or after clustering. Default is 'after'. Other available value is 'before'. 
    """
    
    if f"global_sign_features_{level}_merging" in data.keys():
        important_features = data[f"global_sign_features_{level}_merging"][m_key]["features"]["final_feature_selection_by_ensemble"]
        if isinstance(important_features[0], tuple):
            important_features = [i[0] for i in data[f"global_sign_features_{level}_merging"][m_key]["features"]["final_feature_selection_by_ensemble"]]
        exp_df = data[f"global_sign_features_{level}_merging"][m_key]["sub_high_expression_genes"][m_key]
        exp_df = exp_df.loc[exp_df.feature.isin(important_features)]
        exp_df.loc[:, "log2_fold_change"] = np.log2(exp_df.expression1 / exp_df.expression2)

        # retrieve the feature importance
        score = list()
        name = list()
        try:
            score_pool = dict(data[f"global_sign_features_{level}_merging"][m_key]["features"]["final_feature_selection_by_ensemble"])
            for i in exp_df.feature:
                name.append(i)
                score.append(score_pool[i])

            score_df = pd.DataFrame({"name": name, "score": score}) 

            exp_df = exp_df.join(score_df.set_index("name"), on = "feature")
        except:
            pass
        
    else:
        if f_key == None:
            raise ValueError("f_key can't be None and should be specified when data is the result of running sub_consensus_cluster function")
        else:
            if not re.search("_", f_key):
                try:
                    f_key = data["label_encoder"].inverse_transform([int(f_key)])[0]
                except:
                    f_key = f_key
            else:
                f_key = f_key
            important_features = data[f"sub_sign_genes_{level}_merging"][m_key][f_key]["features"]["final_feature_selection_by_ensemble"]
            if isinstance(important_features[0], tuple):
                important_features = [i[0] for i in data[f"sub_sign_genes_{level}_merging"][m_key][f_key]["features"]["final_feature_selection_by_ensemble"]]
            exp_df = data[f"sub_sign_genes_{level}_merging"][m_key][f_key]["sub_high_expression_genes"][f_key]
            exp_df = exp_df.loc[exp_df.feature.isin(important_features)]
            exp_df.loc[:, "log2_fold_change"] = np.log2(exp_df.expression1 / exp_df.expression2)

        score = list()
        name = list()
        try:
            score_pool = dict(data[f"sub_sign_genes_{level}_merging"][m_key][f_key]["features"]["final_feature_selection_by_ensemble"])
            for i in exp_df.feature:
                name.append(i)
                score.append(score_pool[i])

            score_df = pd.DataFrame({"name": name, "score": score}) 
            exp_df = exp_df.join(score_df.set_index("name"), on = "feature")
        except:
            pass
    
    try:
        exp_df = exp_df.loc[((exp_df.expression1 > ep_cut_off) & (1 - (exp_df.pct2/exp_df.pct1) > pct_cut_off)) | (np.isnan(exp_df.pct2)) ].sort_values("score", ascending = False)
    except:
        exp_df = exp_df.loc[((exp_df.expression1 > ep_cut_off) & (1 - (exp_df.pct2/exp_df.pct1) > pct_cut_off)) | (np.isnan(exp_df.pct2)) ].sort_values("log2_fold_change", ascending = False)

   
    return exp_df