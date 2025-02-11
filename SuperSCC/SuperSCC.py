# %%
import re
import os
import sys
import pickle
import time
import random
import warnings
from datetime import datetime
import logging
from functools import reduce
from multiprocessing import Process
from collections import Counter
import copy

from langchain_community.chat_models import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain.schema import StrOutputParser
import plotly.graph_objects as go
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from dcor import distance_correlation
import scanpy as sc
from scanpy import AnnData
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, gmean
from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef, balanced_accuracy_score, make_scorer
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.feature_selection import RFECV,VarianceThreshold,SelectKBest,chi2, SelectFromModel, f_classif, mutual_info_classif
from sklearn.model_selection import cross_val_score, GridSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler
import matplotlib.pyplot as plt

# %%
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
            

# %%
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

# %%
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
    

# %%
def record_time():
    """
    A function to call out current time.
    """
    current_second_time = datetime.now()
    return current_second_time.strftime("%Y-%m-%d %H:%M:%S")

# %%
def remove_time(string, count = 1):
    """
    A function to remove time information in the stdout message"
    """
    pattern = re.compile("\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")
    res = re.sub(pattern, "", string, count = count)
    return res.rstrip().lstrip()

#%%
def load_pickle_file(file_name, wk_dir = os.getcwd(), recursive = True):
    """
    A function to find pickle file and then load it. 

    Parameters:
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

# %%
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

# %%
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
        counts = quantile_normalize(data.X)        
    
    if save:
        file_name = "Preprocessing_data_{}.pkl".format(record_time())
        with open(file_name, "wb") as output:
            pickle.dump(counts, output)
    
    if logger != None: 
        logger.write("info", "finish pre processing")
    
    return counts

# %%
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

# %%
def model_accuracy(x, y, model, cv = 5, linear_svm_multi_class = "ovr", logistic_multi_class = "ovr", class_weight = "balanced", n_estimators = 100, random_state = 10):
    """
    A function to evulate the model accuracy on cross-validation splits.

    Parameters
    ----------
    data: 
        A pandas data frame object. Rows are cells. Columns are features.
    y: 
        A list of labels for each row in x.
    model: 
        A string to determine the evulation model. Three available values including 'svm', 'logistic' and 'random_forest'.
    cv: 
        A int to control the number of cross validation. 
    linear_svm_multi_class: 
        A string to decide which mode to deal with multiclassification in the linear svm model. Default is "ovr". This parameter only takes effect when model is set in 'svm'. Other available words include "ovo".
    logistic_multi_class: 
        A string to decide which mode to deal with multiclassification in the logistic model. Default is "ovr". Other available words include "multinomial" and "auto".
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    random_state: 
        A int to control the randomness of the bootstrapping of the samples. It takes effect when model is set in 'random_foreast' or in "logistic". Default is 10.
    """
    # prepare containers
    score1_ls = list()
    score2_ls = list()
    score3_ls = list()
    score4_ls = list()

    for i in range(cv):
        X_train, X_test, Y_train, Y_test = train_test_split(x, y, test_size = 0.2)
        if model == "svm":
            classifier = LinearSVC(multi_class =  linear_svm_multi_class, class_weight = class_weight).fit(X_train, Y_train)
        elif model == "random_foreast":
            classifier = RandomForestClassifier(n_estimators = n_estimators, random_state = random_state).fit(X_train, Y_train)
        elif model == "logistic":
            classifier = LogisticRegression(multi_class = logistic_multi_class, random_state = random_state, max_iter=200, class_weight = class_weight).fit(X_train, Y_train)

        y_pred = classifier.predict(X_test)
        score1 = cohen_kappa_score(Y_test, y_pred)
        score1_ls.append(score1)
        score2 = matthews_corrcoef(Y_test, y_pred)
        score2_ls.append(score2)
        score3 = balanced_accuracy_score(Y_test, y_pred)
        score3_ls.append(score3)
        score4 = accuracy_score(Y_test, y_pred)
        score4_ls.append(score4)

    accuracy = "cohen_kappa_score: {}; matthews_corrcoef: {}; balanced_accuracy_score: {}; accuracy_score: {}".format(np.mean(score1), np.mean(score2), np.mean(score3), np.mean(score4))

    return accuracy

# %%
def reverse_rank_weight(x):
    """
    A function to reverse rank weight to makre sure the greatest ranking value representing the most importance.
    """
    assert isinstance(x, dict), "x should be dictionary"
    unique_values = list(set(x.values()))
    unique_values.sort()
    half_length = int(len(unique_values) / 2)
    
    d = dict()
    
    for key, value in x.items():
        if value > unique_values[half_length]:
            for i in range(-1, -half_length-1, -1):
                if unique_values[i] == value:
                    j = i + (i * -2) - 1
                    d[key] = unique_values[j]
        
        else:
            for i in range(half_length+1):
                if unique_values[i] == value:
                    j = i + 1
                    d[key] = unique_values[-j]
    return d    

# %%
def feature_selection(data, 
                      label_column,
                      filename = None,
                      logger = None,
                      rank_method = "dense",
                      merge_rank_method = "geom.mean",
                      variance_threshold = "mean",
                      mutual_info = False, 
                      chi_square_test = False,
                      F_test = True,
                      model = "svm",
                      random_foreast_threshold = None,
                      n_estimators = 100,
                      random_state = 10,
                      normalization_method = "Min-Max",
                      logistic_multi_class = "ovr",
                      linear_svm_multi_class = "ovr",
                      class_weight = "balanced",
                      n_features_to_select = 0.15,
                      step = 100,
                      cv = 5,
                      n_jobs = -1,
                      save = True):
    """
    A function to do feature seletion based on filtering, embedding and wrapping method respectively or combing those methods together.
    
    Parameters
    ----------
    data: 
        A log normalized expression matrix (rows are cells; columns are features) with an extra column containing clustering or cell type labels. 
    label_column: 
        The name of cell type column in the data.
    filename: 
        A string to name the output file. Default is None.
    logger: 
        A log_file object to write log information into disk. Default is None. 
    rank_method: 
        A string to decide which rank method will be used to rank the coefficient values returned by different estimators. Default is "dense". Other available words including "average", "min", "max" and "ordinal".
    merge_rank_method: 
        A string to decide which method will be used to combine the rankings from different estimators. Default is "geom.mean". Other available words including "mean", "min", and "max".
    variance_threshold: 
        A string to decide which variance cutoff is used to filter out features. "zero" or "mean" could be selected. Default is 'mean'.
    mutual_info: 
        A Bool value decide whether a mutual information method is employed to filtering out features further. When it's true, F_test and chi_sqaure_test should be specified in false. Default is False.
    chi_sqaure_test: 
        A Bool value decide whether a chi square test method is employed to filtering out features further. When it's true, F_test and mutual_info should be specified in false. Default is False.
    F_test: 
        A Bool value decide whether a F test method is employed to filtering out features further. When it's true, chi_sqaure_test and mutual_info should be specified in false. Default is True.
    model: 
        A string to decide which model is used by embedding- and wrapper- based feature selection. "random_foreast", "logistic" and "svm" could be selected. Default is 'svm'.
    random_foreast_threshold:
        A float or int value to set the cutoff (feature_importance_) by random foreast model-basedd embedding feature selection. It only takes effect when model is set in 'random_foreast'. Default is None. When it is None, `1 / the number of all features` will be automaticcally used for this value.
    n_estimators: 
        A int to indicate the number of trees in the forest. It only takes effect when model is set in 'random_foreast'. Default is 100.
    random_state: 
        A int to control the randomness of the bootstrapping of the samples. It takes effect when model is set in 'random_foreast' or in "logistic". Default is 10.
    normalization_method: 
        A string to decide how to normalize data. Default is "Min-Max". Other available words including "Standardization".
    logistic_multi_class: 
        A string to decide which mode to deal with multiclassification in the logistic model. Default is "ovr". Other available words include "multinomial" and "auto".
    linear_svm_multi_class: 
        A string to decide which mode to deal with multiclassification in the linear svm model. Default is "ovr". This parameter only takes effect when model is set in 'svm'. Other available words include "ovo".
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    n_featurs_to_selct: 
        A int or float to control the number of features to select. If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default is 0.15. 
    step: 
        A int or float to control the number of features be removed in each round of RFECV algorithm. If greater than or equal to 1, then ``step`` corresponds to the (integer) number of features to remove at each iteration. If within (0.0, 1.0), then ``step`` corresponds to the percentage (rounded down) of features to remove at each iteration. Default is 100.
    cv: 
        A int to decide the number of cross validation in RFECV algorithm. Default is 5.
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    save: 
        A Bool value to decide whether write the output into the disk.
    """
    
    message = "{} {}".format(record_time(), "start feature selection")
    print(message)
    
    params = {
                "label_column": label_column,
                 "file_name": filename,
                 "rank_method": rank_method,
                 "merge_rank_method": merge_rank_method,
                 "variance_threshold": variance_threshold,
                 "mutual_info": mutual_info,
                 "chi_square_test": chi_square_test,
                 "F_test": F_test,
                 "model" : model,
                 "random_foreast_threshold": random_foreast_threshold,
                 "n_estimators": n_estimators,
                 "random_state": random_state,
                 "normalization_method": normalization_method,
                 "logistic_multi_class": logistic_multi_class,
                 "linear_svm_multi_class": linear_svm_multi_class,
                 "class_weight": class_weight,
                 "n_features_to_select": n_features_to_select,
                 "step" : step,
                 "cv": cv, 
                 "n_jobs": n_jobs,
                 "save" : save
            }
    
    if logger != None:
        logger.write("info", remove_time(message))

        for key, value in params.items():
            logger.write("critical", "{}: {}".format(key, value))
    
    # step 1 - convert category label into numeric label
    message = "{} step 1 - converting categoric label into numeric label".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    # if label is string, converse to numeric label
    if any(map(lambda x: isinstance(x, str), data[label_column])):
        le = LabelEncoder().fit(data[label_column])
        label = le.transform(data[label_column])
        data[label_column] = label
    
    # fix duplicated column names by adding ordinal suffix
    data.columns = check_duplicated_elements(data.columns.tolist())
    
    # get data and label
    X = data.iloc[:, 0:-1]
    y = data.iloc[:, -1]
    
    # # get all feature names
    # all_features = list(data.iloc[:, 0:-1].columns)

    # step 2 - do feature selection
    message = "{} step 2 - do feature selection".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    message = "{} ======== filtering based selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
        
    # filtering-based feature selection
    
    # filter out by variance 
    if variance_threshold == "zero":
        var_selector = VarianceThreshold()
        X_var = var_selector.fit_transform(X)
        
        # features selection by ensemble mode
        all_features = list(var_selector.feature_names_in_)
        remained_features = var_selector.get_feature_names_out()
        remained_index = [all_features.index(i) for i in remained_features]
        remained_features_var = var_selector.variances_[remained_index]

        retained_features_ranking_by_variance = dict([*zip(remained_features, rankdata(remained_features_var, method = rank_method))])
        message = "* {} {}".format(record_time(), "feature ranking based on variance of each feature")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        # features selection by intersection mode
        retained_features_by_filter = X.columns[var_selector.get_support()]
        message = "* {} {} features remained after filter out features with 0 variance".format(record_time(), X_var.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

        y_var = y

        try:
            X_var = pd.DataFrame(X_var.todense())
        except:
            X_var = pd.DataFrame(np.array(X_var))
        
        X_var.columns = remained_features
        

    elif variance_threshold == "mean":
        var_selector = VarianceThreshold(np.mean(np.var(np.array(X), axis = 0)))
        X_var = var_selector.fit_transform(X)
        
        # features selection by ensemble mode
        all_features = list(var_selector.feature_names_in_)
        remained_features = var_selector.get_feature_names_out()
        remained_index = [all_features.index(i) for i in remained_features]
        remained_features_var = var_selector.variances_[remained_index]
        
        retained_features_ranking_by_variance = dict([*zip(remained_features, rankdata(remained_features_var, method = rank_method))])
        message = "* {} {}".format(record_time(), "feature ranking based on variance of each feature")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
                
        # features selection by intersection mode
        retained_features_by_filter = X.columns[var_selector.get_support()]
        message = "* {} {} features remained after filter out features below mean variance of all features".format(record_time(), X_var.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

        y_var = y

        try:
            X_var = pd.DataFrame(X_var.todense())
        except:
            X_var = pd.DataFrame(np.array(X_var))
        
        X_var.columns = remained_features

    # filter by chi sqaure test
    if chi_square_test and (F_test == False and mutual_info == False):
        
        # features selection by ensemble mode
        chivalue, pvalues_chi = chi2(X_var, y_var)
        # when pvalue is NaN, it will be replaced with 1
        pvalues_chi = [1 if np.isnan( _ ) else _ for _ in pvalues_chi]
        retained_features_ranking_by_correlation = dict([*zip(remained_features, rankdata(pvalues_chi, method = rank_method))])
        # reverse the ranking weight to make sure the greatest weight representing the most importance
        retained_features_ranking_by_correlation = reverse_rank_weight(retained_features_ranking_by_correlation)
        message = "* {} {}".format(record_time(), "feature ranking based on p value calculated by chi square test to check the correlation between feature and lable")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        # features selection by intersection mode
        chivalue, pvalues_chi = chi2(X_var, y_var)
        k = chivalue.shape[0] - (pvalues_chi > 0.05).sum()
        selector = SelectKBest(chi2, k = k)
        X_fschi = selector.fit_transform(X_var, y_var)
        retained_features_by_filter = retained_features_by_filter[selector.get_support()]
        message = "** {} {} features remained after further chi sqaure test filtering".format(record_time(), X_fschi.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    # filter by F test
    if F_test and (chi_square_test == False and mutual_info == False):
        
        # features selection by ensemble mode
        F, pvalues_f = f_classif(np.array(X_var), y_var)

        # when pvalue is NaN, it will be replaced with 1
        pvalues_f = [1 if np.isnan( _ ) else _ for _ in pvalues_f]
        retained_features_ranking_by_correlation = dict([*zip(remained_features, rankdata(pvalues_f, method =  rank_method))])
        # reverse the ranking weight to make sure the greatest weight representing the most importance
        retained_features_ranking_by_correlation = reverse_rank_weight(retained_features_ranking_by_correlation)
        message = "* {} {}".format(record_time(), "feature ranking based on p value calculated by F test to check the correlation between feature and lable")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        # features selection by intersection mode
        F, pvalues_f = f_classif(np.array(X_var), y_var)
        k = F.shape[0] - (pvalues_f > 0.05).sum()
        selector = SelectKBest(f_classif, k = k)
        X_fsF = selector.fit_transform(X_var, y)
        retained_features_by_filter = retained_features_by_filter[selector.get_support()]
        message = "** {} {} features remained after further F test filtering".format(record_time(), X_fsF.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    # filter by mutual infomation
    if mutual_info and (F_test == False and chi_square_test == False):
        
        # features selection by ensemble mode
        res = mutual_info_classif(X_var, y_var)
        retained_features_ranking_by_correlation = dict([*zip(remained_features, rankdata(res, method = rank_method))])
        message = "* {} {}".format(record_time(), "feature ranking based on p value calculated by mutual information to check the correlation between feature and lable")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        # features selection by intersection mode
        # res = mutual_info_classif(X_var, y_var)
        k = res.shape[0] - sum(res <= 0)
        selector = SelectKBest(mutual_info_classif, k = k)
        X_fsmic = selector.fit_transform(X_var, y)
        retained_features_by_filter = retained_features_by_filter[selector.get_support()]
        message = "** {} {} features remained after further mutual information filtering".format(record_time(), X_fsmic.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    message = "{} ======== embedding based selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))

    # embedding-based on feature selection
    # select by random foreast model
    if model == "random_foreast":
        RFC_ = RandomForestClassifier(n_estimators = n_estimators, random_state = random_state)
        # when random_foreast_threshold is None, 
        #  `1 / number of features` will be used as threshold
        if random_foreast_threshold == None:
            random_foreast_threshold = 1 / X_var.shape[1]
        RFC_embedding_selector = SelectFromModel(RFC_, threshold = random_foreast_threshold)
      
        X_RFC_embedding = RFC_embedding_selector.fit_transform(X_var, y_var)

        # evulate the accuracy of classifier
        accuracy = model_accuracy(X_var, y_var, model = model, cv = cv, linear_svm_multi_class = linear_svm_multi_class, class_weight = class_weight, n_estimators = n_estimators, random_state = random_state, logistic_multi_class = logistic_multi_class)

        # features selection by ensemble mode
        retained_features_ranking_by_embedding = dict([*zip(RFC_embedding_selector.feature_names_in_ , rankdata(RFC_embedding_selector.estimator_.feature_importances_, method = rank_method))])
        message = "* {} {}".format(record_time(), "feature ranking based on feature importance calcualted by random foreast model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

        # when the number of unique labels is 2
        if len(set(y_var)) == 2:
            embedding_positive_marker_selection = dict([*zip(X_var.columns, RFC_embedding_selector.estimator_.feature_importances_)])
        
        # features selection by intersection mode
        retained_features_by_embedding = RFC_embedding_selector.get_feature_names_out()
        message = "* {} {} features remained after random foreast based embedding filtering".format(record_time(), X_RFC_embedding.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    # select by logistic regression model
    elif model == "logistic":
        if normalization_method == "Min-Max":
           X_stand = MinMaxScaler().fit_transform(np.array(X_var))
        elif normalization_method == "Standardization":
           X_stand = StandardScaler().fit_transform(np.array(X_var))
        
        logistic_ = LogisticRegression(multi_class = logistic_multi_class, random_state = random_state, max_iter=200, class_weight = class_weight)
        log_embedding_selector = SelectFromModel(logistic_, norm_order = 1)
        X_log_embedding = log_embedding_selector.fit_transform(np.array(X_stand), y_var)

        # evulate the accuracy of classifier
        accuracy = model_accuracy(X_stand, y_var, model = model, cv = cv, linear_svm_multi_class = linear_svm_multi_class, class_weight = class_weight, n_estimators = n_estimators, random_state = random_state, logistic_multi_class = logistic_multi_class)

        # features selection by ensemble mode
        retained_features_ranking_by_embedding = dict([*zip(X_var.columns, rankdata(np.sum(np.abs(log_embedding_selector.estimator_.coef_), axis = 0), method = rank_method))]) # sum all absolute coefficients of each feature 
        
        # when the number of unique labels is 2
        if len(set(y_var)) == 2:
            embedding_positive_marker_selection = dict([*zip(X_var.columns, log_embedding_selector.estimator_.coef_)])
        
        message = "* {} {}".format(record_time(), "feature ranking based on coefficient calcualted by logistic model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        # features selection by intersection mode
        retained_features_by_embedding = log_embedding_selector.get_feature_names_out()

        # test1 = retained_features_by_embedding
        # test2 = X_var.columns[log_embedding_selector.get_support()]
        # print(len(test1), np.sum(np.array(test1) == np.array(test2)))

        message = "* {} {} features remained after logistic regression based embedding filtering".format(record_time(), X_log_embedding.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    # select by SVM model
    elif model  == "svm":
        # Standardization or Min-Max scaling is required for linear SVM 
        if normalization_method == "Min-Max":
            X_stand = MinMaxScaler().fit_transform(np.array(X_var))
        elif normalization_method == "Standardization":
            X_stand = StandardScaler().fit_transform(np.array(X_var))
        
        SVC_ = LinearSVC(multi_class = linear_svm_multi_class, class_weight = class_weight)
        SVC_ = SVC_.fit(X_stand, y_var)
        SVC_embedding_selector = SelectFromModel(SVC_, norm_order = 1, prefit=True)
        X_SVC_embedding = SVC_embedding_selector.transform(np.array(X_stand))

        # evulate the accuracy of classifier
        accuracy = model_accuracy(X_stand, y_var, model = model, cv = cv, linear_svm_multi_class = linear_svm_multi_class, class_weight = class_weight, n_estimators = n_estimators, random_state = random_state, logistic_multi_class = logistic_multi_class)
        
        # features selection by ensemble mode
        feature_importance = rankdata(np.sum(np.abs(SVC_.coef_), axis = 0), method = rank_method) # sum all absolute coefficients of each feature
        retained_features_ranking_by_embedding = dict([*zip(X_var.columns, feature_importance)])
        message = "* {} {}".format(record_time(), "feature ranking based on coefficient calcualted by linear SVM model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

        # when the number of unique labels is 2
        if len(set(y_var)) == 2:
            embedding_positive_marker_selection = dict([*zip(X_var.columns, SVC_.coef_[0])])
       
        # features selection by intersection mode
        retained_features_by_embedding = X_var.columns[SVC_embedding_selector.get_support()]
        message = "* {} {} features remained after svm based embedding filtering".format(record_time(), X_SVC_embedding.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    # selected by wrapping method
    message = "{} ======== wrapping based selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    if n_features_to_select == None:
        # when features to select is None,
        # 50% of all features will be used as threshold
        n_features_to_select = int(X_var.shape[1] * 0.5)
    elif isinstance(n_features_to_select, float):
        n_features_to_select = int(n_features_to_select * X_var.shape[1])

    if model == "random_foreast":
        RFC_ = RandomForestClassifier(n_estimators = n_estimators, random_state = random_state)
        RFC_wrapping_selector = RFECV(RFC_, min_features_to_select = n_features_to_select, step = step, cv = cv, n_jobs= n_jobs)
        X_RFC_wrapping = RFC_wrapping_selector.fit_transform(X_var, y_var)
        
        # features selection by ensemble mode
        retained_features_ranking_by_wrapping = dict([*zip(RFC_wrapping_selector.feature_names_in_ , RFC_wrapping_selector.ranking_)])
        retained_features_ranking_by_wrapping = reverse_rank_weight(retained_features_ranking_by_wrapping)

        message = "* {} {}".format(record_time(), "feature ranking based on ranking provided RFE - random foreast model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

        
        # features selection by intersection mode
        retained_features_by_wrapping = X_var.columns[RFC_wrapping_selector.support_]

        message = "* {} {} features remained after RFE - random foreast based wrapping filtering".format(record_time(), X_RFC_wrapping.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    elif model == "logistic":
        logistic_ = LogisticRegression(multi_class = logistic_multi_class, random_state = random_state, max_iter=200, class_weight = class_weight)
        log_wrapping_selector = RFECV(logistic_, min_features_to_select = n_features_to_select, step = step, cv = cv, n_jobs = n_jobs)
        X_log_wrapping = log_wrapping_selector.fit_transform(np.array(X_stand), y_var)
        
        # features selection by ensemble mode
        retained_features_ranking_by_wrapping = dict([*zip(X_var.columns, log_wrapping_selector.ranking_)])
        retained_features_ranking_by_wrapping = reverse_rank_weight(retained_features_ranking_by_wrapping)
        message = "* {} {}".format(record_time(), "feature ranking based on ranking provided RFE - logistic model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
                
        # features selection by intersection mode
        retained_features_by_wrapping = X_var.columns[log_wrapping_selector.support_]
        message = "* {} {} features remained after RFE - logistic regression based wrapping filtering".format(record_time(), X_log_wrapping.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    elif model == "svm":
        SVC_ = LinearSVC(multi_class = linear_svm_multi_class, class_weight = class_weight)
        SVC_wrapping_selector = RFECV(SVC_, min_features_to_select = n_features_to_select, step = step, cv = cv, n_jobs = n_jobs)
        X_SVC_wrapping = SVC_wrapping_selector.fit_transform(X_stand, y_var)
        
        # features selection by ensemble mode
        retained_features_ranking_by_wrapping = dict([*zip(X_var.columns, SVC_wrapping_selector.ranking_)])
        retained_features_ranking_by_wrapping = reverse_rank_weight(retained_features_ranking_by_wrapping)
        message = "* {} {}".format(record_time(), "feature ranking based on ranking provided RFE - SVM model")
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
                
       # features selection by intersection mode
        retained_features_by_wrapping = X_var.columns[SVC_wrapping_selector.support_]
        message = "* {} {} features remained after RFE - svm based wrapping filtering".format(record_time(), X_SVC_wrapping.shape[1])
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))

    message = "{} ======== final feature selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
                
        
    #  features selection by ensemble mode
    rank_ls = list()

    for feature in X_var.columns:
        ranks = [retained_features_ranking_by_variance[feature], retained_features_ranking_by_correlation[feature], retained_features_ranking_by_embedding[feature], retained_features_ranking_by_wrapping[feature]]
        if merge_rank_method == "max":
            rank_ls.append((feature, max(ranks)))
        elif merge_rank_method == "median":
            rank_ls.append((feature, np.median(np.array(ranks))))
        elif merge_rank_method == "mean":
            rank_ls.append((feature, np.mean(np.array(ranks))))
        elif merge_rank_method == "geom.mean":
            rank_ls.append((feature, gmean(np.array(ranks))))

    rank_ls = sorted(rank_ls, key = lambda x: -x[1])
    rank_ls = [ _ for _ in rank_ls[0:n_features_to_select] ]

    #  features selection by intersection mode
    retained_features_by_filter = set(retained_features_by_filter)
    retained_features_by_embedding = set(retained_features_by_embedding)
    retained_features_by_wrapping = set(retained_features_by_wrapping)

    final_feture_selection = reduce(lambda x,y: x.intersection(y), [retained_features_by_embedding, retained_features_by_filter, retained_features_by_wrapping])

    output1 = {"retained_features_ranking_by_variance": retained_features_ranking_by_variance,
              "retained_features_ranking_by_correlation": retained_features_ranking_by_correlation,
              "retained_features_ranking_by_embedding": retained_features_ranking_by_embedding,
              "retained_features_ranking_by_wrapping": retained_features_ranking_by_wrapping,
              "final_feature_selection_by_ensemble": rank_ls,
              "model_accuracy": accuracy,
              "params_used_for_feature_selection": params}
    
    # only run when the number of unique labels is 2
    if len(set(y_var)) == 2:
        output1.update({"positive_marker_selection": embedding_positive_marker_selection})
              
    output2 = {"retained_features_by_filtering": retained_features_by_filter,
              "retained_features_by_embedding": retained_features_by_embedding,
              "retained_features_by_wrapping": retained_features_by_wrapping,
              "final_feature_selection_by_intersection": final_feture_selection,
              "model_accuracy": accuracy,
              "params_used_for_feature_selection": params}
    
    message = "* {} {} features remained after intersecting the key features found by filtering, embedding and wrapping-based feature selection methods".format(record_time(), len(final_feture_selection))
    message2 = "* {} select top {} features ranked by using {} method to combine feature rankings obtained by different estimators".format(record_time(), len(rank_ls), merge_rank_method)
    print(message2)
    print(message)
    if logger != None:
        logger.write("info", remove_time(message2))
        logger.write("info", remove_time(message))
    
    if save == True:
        if filename != None:
            filename1 = filename + "_" + model + "_" + "ensemble" + "_" + "feature_selection" + "_" + record_time() + ".pkl"
            with open(filename1, "wb") as file:
                pickle.dump(output1, file)

            filename2 = filename + "_" + model + "_" + "intersection" + "_" + "feature_selection" + "_" + record_time() + ".pkl"
            with open(filename2, "wb") as file:
                pickle.dump(output2, file)
        else:
            filename1 = model + "_" + "ensemble" + "_" + "feature_selection" + "_" + record_time() + ".pkl"
            with open(filename1, "wb") as file:
                pickle.dump(output1, file)

            filename2 = model + "_" + "intersection" + "_" + "feature_selection" + "_" + record_time() + ".pkl"
            with open(filename2, "wb") as file:
                pickle.dump(output2, file)
    
    # merge two feature selection results
    output1.update(output2)

    try:
        output1.update({"label_classes": le.classes_})
    except:
        message = "{} finish feature selection".format(record_time())
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
    
    return output1

# %%
def model_training(data,
                   label_column,
                   features,
                   model,
                   normalization_method = "Min-Max",
                   feature_selection_params = None,
                   parameters = None,
                   random_state = 10,
                   cv = 5,
                   n_jobs = -1,
                   filename = None,
                   save = True ,
                   logger = None
                  ):
    """
    A function to do model training based on selected features and model.
    
    Parameters
    ----------
    data: 
        A log normalized expression matrix. Rows are cells; Columns are features. 
    label_column: 
        A string to specify the name of cell type column in the data.
    features: 
        A list to decide the features remained for model training.
    model: 
        A string to decide the training model. "random_foreast", "svm" or "logistic" could be selected.
    normalization_method: 
        A string to decide how to normalize data. Default is "Min-Max". Other available words including "Standardization". If features are returned by feature_selection function, this parameter should be consistent with the counterpart used for running feature_selection function.
    parameters: 
        A dict to elucidate the parameter name and value pair that should be searched in the grid search algorithm for corresponding model. Default is None.
    random_state: 
        A int to control the randomness of the bootstrapping of the samples. It takes effect when model is set in 'random_foreast' or in "logistic". Default is 10.
    cv: 
        A int to decide the number of cross validation for grid serach. Default is 5.
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    filename: 
        A string to indicte the name of output file.
    save: 
        A Bool value to decide whether the result will be written into disk. Default is True.
    logger: 
        A log_file object to write log information into disk. Default is None. 
    """
    
    message = "{} start model training".format(record_time())
    print(message)
    print("{} model traning based on {} algorithm".format(record_time(), model))
    
    if logger != None:
        logger.write("info", remove_time(message))
        
        arguments = {"label_column": label_column,
                     "normalization_method": normalization_method,
                     "model": model,
                     "random_state": random_state,
                     "cv": cv,
                     "n_jobs": n_jobs,
                     "save": save}
        
        for key,value in arguments.items():
            logger.write("critical", "{}: {}".format(key, value))
        
        logger.write("info", "model traning based on {} algorithm".format(model))
    
    # get data and target
    X = data.loc[:, data.columns != label_column]
    y = data.loc[:, data.columns == label_column][label_column].values
    
    # only keep informative features
    X = data.loc[:, data.columns.isin(features)] 
    X = X.loc[:, list(features)] # reorder the columns
    
    # normalize data if necessary
    if normalization_method == "Min-Max":
        message = "{} doing {} normalization".format(record_time(), normalization_method)
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
            
        X = MinMaxScaler().fit_transform(np.array(X))
        X = pd.DataFrame(np.array(X))
        X.columns = features
        
    elif normalization_method == "Standardization":
        message = "{} doing {} normalization".format(record_time(), normalization_method)
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
            
        X = StandardScaler().fit_transform(np.array(X))
        X = pd.DataFrame(np.array(X))
        X.columns = features
        
    # check whethere there are string labels, if so, convert them to numeric labels
    if any(map(lambda x: isinstance(x, str), y)):
        message = "{} doing label encoding".format(record_time())
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
            
        le = LabelEncoder().fit(y)
        y_trans = le.transform(y)
    else:
        y_trans = y
    
    
    # model training on random foreast model
    if model == "random_foreast":
        if parameters == None:
            parameters = {
                         "n_estimators" : np.arange(10, 101, 10),
                         "criterion" : ["gini", "entropy"],
                         # "max_depth": np.linspace(10, 50, 5),
                         "max_features": np.linspace(0.2, 1, 5)
                         }
                         #"min_samples_leaf": np.arange(10, 300, 20),
                         #"min_samples_split": np.arange(2, 100, 10)}
        
        message = "{} grid search below paramters getting the best model".format(record_time())
        print(message)
        
        if logger != None:
            logger.write("info", remove_time(message))
        
        for key,value in parameters.items():
            print("* {}: {}".format(key, value))
            if logger != None:
                logger.write("critical", "{}: {}".format(key, value))
           
        RFC_ = RandomForestClassifier(random_state = random_state)
        GS = GridSearchCV(RFC_, parameters, cv = cv, n_jobs = n_jobs)
        GS.fit(X, y_trans)
        best_parameters = GS.best_params_
        best_core = GS.best_score_
    
    # model training on logistic regression model
    elif model == "logistic":
        if parameters == None:
            parameters = {
                          "penalty": ["l1", "l2"],
                          "C": np.linspace(0, 1, 5),
                          "multi_class": ["ovr", "multinomial"]
                         }

        message = "{} grid search below paramters getting the best model".format(record_time())
        print(message)
        
        if logger != None:
            logger.write("info", remove_time(message))
        
        for key,value in parameters.items():
            print("* {}: {}".format(key, value))
            if logger != None:
                logger.write("critical", "{}: {}".format(key, value))
        
        logistic_ = LogisticRegression(random_state = random_state)
        GS = GridSearchCV(logistic_, parameters, cv = cv, n_jobs = n_jobs)
        GS.fit(X, y_trans)
        best_parameters = GS.best_params_
        best_core = GS.best_score_
    
    # model training on SVM model
    elif model == "svm":
        if parameters == None:
            parameters = {"C": np.linspace(0.01, 1, 10),
                         "kernel": ["rbf", "poly", "sigmoid", "linear"]}
                         #"coef0": np.linspace(0,5,10)}
        
        message = "{} grid search below paramters getting the best model".format(record_time())
        print(message)
        
        if logger != None:
            logger.write("info", remove_time(message))
        
        for key,value in parameters.items():
            print("* {}: {}".format(key, value))
            if logger != None:
                logger.write("critical", "{}: {}".format(key, value))
        
        SVM_ = SVC(decision_function_shape = "ovr")
        GS = GridSearchCV(SVM_, parameters, cv = cv, n_jobs = n_jobs)
        GS.fit(X, y_trans)
        best_parameters = GS.best_params_
        best_core = GS.best_score_

    output = {"model": GS,
              "best_score": best_core,
              "best_parameters": best_parameters,
              "features_used_for_training": features,
              "params_used_for_feature_selection": feature_selection_params}
    try:
        output.update({"original_label_classes": le.classes_})
    except:
        pass
    
    if logger != None:
        for key, value in output.items():
            logger.write("info", "{}: {}".format(key, value))
    
    if save == True:
        if filename == None:
            filename = model + "_" + "training_model" + "_" + record_time() + ".pkl"
        else:
            filename = filename + "_" + model + "_"  + "training_model" + "_" + record_time() + ".pkl"
        
        with open(filename, "wb") as file:
            pickle.dump(output, file)
    
    message = "{} finish model training".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    return(output)

# %%
def predict_label(query_data, 
                  models, 
                  wk_dir = os.getcwd(), 
                  normalization_method = "Min-Max",
                  filname = None,
                  save = True,
                  logger = None
                  ):
    """
    A function to predict label and score the prediction by using pre-trained model.
    
    Parameters
    ----------
    query_data: 
       A log normalized expression matrix. Rows are cells; Columns are features. 
    models: 
        A regular expression to match the name of pre-trained models. Otherwise, a dict returned by load_pick_file function. 
    wk_dir: 
        A string to specify the directory where pre-trained model files should be searched.
    normalization_method: 
        A string to decide how to normalize query data. Default is "Min-Max". Other available words including "Standardization". This parameter should be same as the one when running the model training, thereby guaranteeing same normalization on query and reference. 
    filename: 
        A string decide the name of output. Default is None. 
    save: 
        A Bool value to decide whether the result will be written into disk. Default is True.
    logger: 
        A log_file object to write log information into disk. Default is None. 
    """
    
    # load pre-trained models in
    if isinstance(models, dict):
        models = models
    else:
        models = load_pickle_file(models, wk_dir = wk_dir)
   
    prediction_res = dict()
    
    for model_name in models:
        message = "{} start label prediction based on {} model".format(record_time(), model_name)
        print(message)
        if logger != None:
            logger.write("info", remove_time(message))
        
        
        # get the one of pre-trained model
        model = models[model_name]["model"]
        model_class = models[model_name]["original_label_classes"]

        # get the features used in the pre-trained model
        selected_features = model.feature_names_in_
        
        # make query data only with features in the pre-trained model
        query_data = query_data.loc[:, query_data.columns.isin(selected_features)]
        
        # check whether query data includes all features used for pre-trained model
        feature_consistence = int(query_data.shape[1]) == len(selected_features)
        
        # when query data doesn't include all features used for pre-trained model
        if not feature_consistence:
            without_features = []
            for i in selected_features:
                if i not in query_data.columns:
                    without_features.append(i)
        
            message = "{} {} features used for pre-trained model are not in the query dataset".format(record_time(), len(without_features))
            print(message)
            if logger != None:
                logger.write("info", message)

            message = "{} insert lost features back to query dataset and make those features' expression value with 0"
            print(message)
            if logger != None:
                logger.write("info", message)

            for lost in without_features:
                query_data.loc[:, lost] = 0
        
        # reorder columns to make sure
        # the order of features being consistent with model
        query_data = query_data.loc[:, list(selected_features)]
                
        # do normalization
        if normalization_method == "Min-Max":
            query_data = MinMaxScaler().fit_transform(np.array(query_data))
        elif normalization_method == "Standilization":
            query_data = StandardScaler().fit_transform(np.array(query_data))
        
        # convert to dataframe and name the column
        query_data = pd.DataFrame(np.array(query_data))
        query_data.columns = list(selected_features)
        
        # do prediction
        try:
            prediction = model.predict(query_data)
            try:
                prediction_prob = model.predict_proba(query_data)
            except:
                pass
        except:
            prediction = model.predict(np.array(query_data))
            try:
                prediction_prob = model.predict_proba(np.array(query_data))
            except:
                pass
        
        try:
            prediction_res[model_name] = {"prediction": [model_class[i] for i in prediction],
                                          "prediction_prob": prediction_prob}
        except:
            prediction_res[model_name] = {"prediction": [model_class[i] for i in prediction]}
        
        message = "{} finish label prediction based on {}".format(record_time(), model_name)
        print(message)
        if logger != None:
                logger.write("info", remove_time(message))
        
        if save:
            model_name = re.sub("\.pkl", "", remove_time(model_name)).rstrip()
            output_name = model_name + "prediction_result" + "_" + record_time() + ".pkl"

            with open(output_name, "wb") as file:
                pickle.dump(prediction_res, file)
    
    return prediction_res       

# %%
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

# %%
def find_signature_genes(data, 
                         label_column = "cluster", 
                         n_features_to_select = 0.15, 
                         ratio_of_none_zero_counts = 0.1,
                         class_weight = "balanced",
                         HVG = True,
                         save = False,
                         filename = None,
                         n_jobs = -1,
                         logger = None,
                         **kwargs):
    
    """
    A simple wrapper of feature_selection function for finding highly variable genes for overall clusters or for markers of the corresponding cluster.

    Parameters
    ----------
    data: 
        A log normalized expression matrix (rows are cells; columns are features) with an extra column containing clustering or cell type labels. 
    label_column: 
        A string to specify the name of cell type column in the data.
    n_features_to_select: 
        A int or float to control the number of features to select. If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default is 0.15.
    ratio_of_none_zero_counts: 
        A float to determine the cutoff in which a feature will be omited when below specified value. A higher value leads to less features will be kept for calculation. 
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    HVG: 
        A Bool value to decider whether only detected highly variable features will be returned. 
    save: 
        A Bool value to decide whether the result will be written into disk. Default is True.
    filename: 
        A string to name the output file. Default is None.
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    logger: 
        A log_file object to write log information into disk. Default is None. 
    **kwargs: 
        Other paremeters passed to feature_selection function. 
    """
    
    # step1 - do feature selection
    if isinstance(n_features_to_select, float):
        n_features_to_select = int(n_features_to_select * data.shape[1])
    
    features = feature_selection(data = data.copy(), 
                                 label_column = label_column, 
                                 n_features_to_select = n_features_to_select, 
                                 save = save,
                                 filename = filename,
                                 class_weight = class_weight, 
                                 logger = logger, 
                                 n_jobs = n_jobs,
                                 **kwargs)
    if not HVG:
        # group data
        group_data = data.groupby(label_column)
        
        # prepare containers
        sub_high_expression_genes = dict()
        
        for key in list(set(data.loc[:, label_column])):
                sub_data = group_data.get_group(key)
                sub_data = sub_data.iloc[:, 0:-1]
                try:
                    key = features["label_classes"][key]
                except:
                    key = key
                
                sub_data_none_zero = map(lambda x: sub_data.loc[:, x][sub_data.loc[:, x] > 0], sub_data.columns)
                sub_data_none_zero = [*map(lambda x: x.shape[0], sub_data_none_zero)]
                sub_data_none_zero_ratio  = pd.DataFrame(np.array(sub_data_none_zero) / sub_data.shape[0], columns = ["ratio_of_none_zero_counts"])
                sub_data_none_zero_ratio.loc[:, "feature"] = sub_data.columns
                
                # get the features with at least 10% (default) cells having none-zero expression values
                remained_features = sub_data_none_zero_ratio.loc[sub_data_none_zero_ratio.ratio_of_none_zero_counts > ratio_of_none_zero_counts, "feature"]
                
                # get the mean expression of each feature
                each_feature_mean = sub_data.mean()
                intermediate = pd.DataFrame(each_feature_mean, columns=["expression"])
                intermediate.loc[:, "feature"] = intermediate.index
                intermediate = intermediate.join(sub_data_none_zero_ratio.set_index("feature"), on = "feature")

                # only keep features with at least 10% (default) cells having none-zero expression values
                intermediate = intermediate.loc[intermediate.feature.isin(remained_features), :] 
                intermediate = intermediate.sort_values(by = "expression")
                intermediate.loc[:, "rank"] = rankdata(intermediate.loc[:, "expression"], method="dense")
                sub_high_expression_genes[str(key)] = intermediate
        
        condition = np.array(list(sub_high_expression_genes.keys()))
        data1 = sub_high_expression_genes[condition[condition != "OnevsRest"][0]]
        data1.rename(columns={'ratio_of_none_zero_counts': 'pct1', "expression": "expression1", "rank": "rank1"}, inplace=True)
        data2 = sub_high_expression_genes["OnevsRest"] 
        data2.rename(columns={'ratio_of_none_zero_counts': 'pct2', "expression": "expression2", "rank": "rank2"}, inplace=True)

        data3 = data1.join(data2.set_index("feature"), on = "feature")
        data3 = data3[["feature", "expression1", "pct1", "rank1", "expression2", "pct2", "rank2"]]

        scores = dict(features["final_feature_selection_by_ensemble"])
        score_keys = list(scores.keys())
        for j in score_keys:
            if j not in data3.feature.tolist():
                scores.pop(j)

        new_score_key = [i for i in scores.keys()]
        new_score_value = [i for i in scores.values()]
        new_scores = pd.DataFrame({"feature": new_score_key, "score": new_score_value})

        data4 = data3.join(new_scores.set_index("feature"), on = "feature")

        output = {"features": features, "sub_high_expression_genes": {condition[condition != "OnevsRest"][0]: data4}}
    else:
        output = {"features": features}

    if save:
        with open("markers_of_clusters_after_1st_merging_{}.pkl".format(record_time()), "wb") as file:
            pickle.dump(output, file)

    return output

# %%
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

# %%
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

# %%
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

# %%
def jaccard_score(x, y):
    """
    A function to calculate jaccard index.
    """
    x = set(x)
    y = set(y)
    n1 = len(x.intersection(y))
    n2 = len(x.union(y))
    
    return n1/n2

# %%
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

# %%
def pairwise_compare(data, level = "between", focus = "expression1", ep_cut_off = 1, pct_cut_off = 0.5, only_positive = True):

    """
    A function to calculate the jaccard index and distance correlation between features of paired group.

    Parameters
    ----------
    data: 
        A dict returned by global_consensus_cluster or by sub_consensus_cluster function.
    level: 
        A string. "between" should be used when data is the result of global_consensus_cluster. "within" should be used when data is the result of sub_consensus_cluster. Default is "between".
    focus: 
        A string to decide the value used for distance correlation calculation. "expression1" - the mean expression of feauture in the interest cluster will be used. "rank1" - the rank of feature based on expression in the interest cluster will be used. Default is 'expression1'.
    ep_cut_off: 
        A int or float value to decide the expression cutoff to select positive informative features.
    pct_cut_off: 
        A int or float value to decide the expression percentage cutoff to select positive informative features.
    only_positive: 
        A Bool value to decide whether only positive informative features will be used to calculate the feature similarities between clusters. Default is True. 
    """
    
    assert isinstance(data, dict), "data should be the result after running consensus_cluster function"
    
    jaccard_index = dict()
    distance_cor = dict()
   
    # if level == "between":
      
    keys = sorted(list(data["global_sign_features"].keys()))
    
    for key1 in range(len(keys)-1):
        for key2 in range(key1+1, len(keys)):
            high_expression_1 = data["global_sign_features"][keys[key1]]["sub_high_expression_genes"][keys[key1]]
            sign_genes_1 = data["global_sign_features"][keys[key1]]["features"]["final_feature_selection_by_ensemble"]

            high_expression_2 = data["global_sign_features"][keys[key2]]["sub_high_expression_genes"][keys[key2]]
            sign_genes_2 = data["global_sign_features"][keys[key2]]["features"]["final_feature_selection_by_ensemble"]    

            if isinstance(sign_genes_1[0], tuple) and isinstance(sign_genes_2[0], tuple):  
                sub_data1_feature_names = [i[0] for i in sign_genes_1]
                sub_data2_feature_names = [i[0] for i in sign_genes_2]
            elif isinstance(sign_genes_1[0], tuple) and isinstance(sign_genes_2[0], tuple) == False:
                sub_data1_feature_names = [i[0] for i in sign_genes_1]
                sub_data2_feature_names = sign_genes_2
            elif isinstance(sign_genes_1[0], tuple) == False and isinstance(sign_genes_2[0], tuple):
                sub_data1_feature_names = sign_genes_1
                sub_data2_feature_names = [i[0] for i in sign_genes_2]
            else:
                sub_data1_feature_names = sign_genes_1
                sub_data2_feature_names = sign_genes_2
            
            if only_positive:
                new_data = {"global_sign_features_before_merging": data["global_sign_features"]}
                positive_genes1 = retrieve_positive_markers(new_data, keys[key1], ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, level = "before")
                positive_genes2 = retrieve_positive_markers(new_data, keys[key2], ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, level = "before")

            
                # calculate jaccard index based on all positive genes between comparisons
                jaccard_index["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] = jaccard_score(positive_genes1.feature, positive_genes2.feature)

                # calculate distance correlation based on shared positive genes between comparisons
                shared_positive_genes = set(positive_genes1.feature).intersection(positive_genes2.feature)
                sub_data1_feature_ranking = positive_genes1.loc[positive_genes1.feature.isin(shared_positive_genes)].sort_values("feature").loc[:, focus].values
                sub_data2_feature_ranking = positive_genes2.loc[positive_genes2.feature.isin(shared_positive_genes)].sort_values("feature").loc[:, focus].values

                distance_cor["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] = distance_correlation(sub_data1_feature_ranking, sub_data2_feature_ranking)

            else:
                # calculate jaccard index based on all positive genes between comparisons
                jaccard_index["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] =  jaccard_score(sub_data1_feature_names, sub_data2_feature_names)    

                # calculate the distance correlation between the all shared genes consisting of positive/negative markers between comparisons
                shared_sign_features = set(sub_data1_feature_names).intersection(sub_data2_feature_names)
                sub_data1_feature_ranking = high_expression_1.loc[high_expression_1.feature.isin(shared_sign_features), :]
                sub_data2_feature_ranking = high_expression_2.loc[high_expression_2.feature.isin(shared_sign_features), :]
                
                shared_high_exp_features = set(sub_data1_feature_ranking.feature).intersection(sub_data2_feature_ranking.feature)
                sub_data1_feature_ranking = sub_data1_feature_ranking.loc[sub_data1_feature_ranking.feature.isin(shared_high_exp_features), :].sort_values("feature")
                sub_data2_feature_ranking = sub_data2_feature_ranking.loc[sub_data2_feature_ranking.feature.isin(shared_high_exp_features), :].sort_values("feature")
                
                distance_cor["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] = distance_correlation(sub_data1_feature_ranking.loc[:, focus].values, sub_data2_feature_ranking.loc[:, focus].values)
                

    # elif level == "within":
    #     wk_data = data["sub_sign_genes"]
    #     global_keys = wk_data.keys()
    #     group_d = dict()

    #     for j in global_keys:
    #         group_d[j] = list(wk_data[j].keys())
        
    #     for i in group_d.keys():

    #         keys = group_d[i]
    #         sign_genes = dict()
    #         high_expression = dict()
    #         if len(keys) < 2:
    #             distance_cor["global cluster {}".format(i)] = "no positive markers are found in global cluster {}".format(i)
    #             jaccard_index["global cluster {}".format(i)] = "no positive markers are found in global cluster {}".format(i)
    #             continue

    #         for j in keys:
    #             sign_genes[j] = wk_data[i][j]
    #             high_expression[j] = wk_data[i][j]["sub_high_expression_genes"][j]

    #         for key1 in range(len(keys)-1):

    #             for key2 in range(key1+1, len(keys)):

    #                 sub_data1 = sign_genes[keys[key1]]["features"]["final_feature_selection_by_ensemble"]
    #                 sub_data2 = sign_genes[keys[key2]]["features"]["final_feature_selection_by_ensemble"]

    #                 if len(sub_data1) == 0:
    #                     distance_cor["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = "no positive markers are found in cluster {}".format(keys[key1])
    #                     jaccard_index["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = "no positive markers are found in cluster {}".format(keys[key1])
    #                     continue

    #                 elif len(sub_data2) == 0:
    #                     distance_cor["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = "no positive markers are found in cluster {}".format(keys[key2])
    #                     jaccard_index["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = "no positive markers are found in cluster {}".format(keys[key2])
    #                     continue
                    
    #                 if isinstance(sub_data1[0], tuple) and isinstance(sub_data2[0], tuple):  
    #                     sub_data1_feature_names = [i[0] for i in sub_data1]
    #                     sub_data2_feature_names = [i[0] for i in sub_data2]
    #                 elif isinstance(sub_data1[0], tuple) and isinstance(sub_data2[0], tuple) == False:
    #                     sub_data1_feature_names = [i[0] for i in sub_data1]
    #                     sub_data2_feature_names = sub_data2
    #                 elif isinstance(sub_data1[0], tuple) == False and isinstance(sub_data2[0], tuple):
    #                     sub_data1_feature_names = sub_data1
    #                     sub_data2_feature_names = [i[0] for i in sub_data2]
    #                 else:
    #                     sub_data1_feature_names = sub_data1
    #                     sub_data2_feature_names = sub_data2

    #                 if only_positive:
    #                     new_data = {"sub_sign_genes_before_merging": data["sub_sign_genes"]}
    #                     positive_genes1 = retrieve_positive_markers(new_data, i, keys[key1], ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, level = "before")
    #                     positive_genes2 = retrieve_positive_markers(new_data, i, keys[key2], ep_cut_off = ep_cut_off, pct_cut_off = pct_cut_off, level = "before")

    #                     # calculate jaccard index based on all positive genes between comparisons
    #                     jaccard_index["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] =  jaccard_score(positive_genes1.feature.tolist(), positive_genes2.feature.tolist())

    #                     # calculate distance correlation based on shared positive genes between comparisons
    #                     shared_positive_genes = set(positive_genes1.feature.tolist()).intersection(positive_genes2.feature)
    #                     sub_data1_feature_ranking = positive_genes1.loc[positive_genes1.feature.isin(shared_positive_genes), :].sort_values("feature").loc[:, focus].values
    #                     sub_data2_feature_ranking = positive_genes2.loc[positive_genes2.feature.isin(shared_positive_genes), :].sort_values("feature").loc[:, focus].values

    #                     distance_cor["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = distance_correlation(sub_data1_feature_ranking, sub_data2_feature_ranking)
                    
    #                 else:
    #                     # calculate the jaccard score
    #                     jaccardScore = jaccard_score(sub_data1_feature_names, sub_data2_feature_names)
    #                     jaccard_index["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = jaccardScore
                        
    #                     # calculate the pearson coefficient between the shared features of two datasets
    #                     shared_sign_features = set(sub_data1_feature_names).intersection(sub_data2_feature_names)
    #                     sub_data1_feature_ranking = high_expression[keys[key1]].loc[high_expression[keys[key1]].feature.isin(shared_sign_features), :]
    #                     sub_data2_feature_ranking = high_expression[keys[key2]].loc[high_expression[keys[key2]].feature.isin(shared_sign_features), :]
                        
    #                     shared_high_exp_features = set(sub_data1_feature_ranking.feature).intersection(sub_data2_feature_ranking.feature)
    #                     sub_data1_feature_ranking = sub_data1_feature_ranking.loc[sub_data1_feature_ranking.feature.isin(shared_high_exp_features), :].sort_values("feature")
    #                     sub_data2_feature_ranking = sub_data2_feature_ranking.loc[sub_data2_feature_ranking.feature.isin(shared_high_exp_features), :].sort_values("feature")
                        
    #                     d_cor = distance_correlation(sub_data1_feature_ranking.loc[:, focus].values, sub_data2_feature_ranking.loc[:, focus].values)
    #                     distance_cor["{} vs {} in global cluster {}".format(keys[key1], keys[key2], i)] = d_cor # 0 if np.isnan(d_cor) else d_cor
                    
    res = {"jaccard_index": jaccard_index, "distance_cor": distance_cor}
    
    return res

# %%
def detect_former_num(data):
    """
    A function to find the first numeric elements in a string.
    """
    res = re.split("vs", data)
    output = [*re.compile("((\d+&?){1,100}|(\d+))").finditer(res[0])][0]
    return res[0][output.start():output.end()]

# %%
def find_markers_ovr(
                    data,
                    label_column = "cluster", 
                    n_features_to_select = 0.15, 
                    ratio_of_none_zero_counts = 0.1,
                    class_weight = "balanced",
                    save = False,
                    logger = None,
                    filename = None,
                    n_jobs = -1,
                    **kwargs
                    ):
    
    """
    A simple wrapper of fing_signature_genes function to run in one vs rest mode.

    Parameters
    ----------
    data: 
        A log normalized expression matrix (rows are cells; columns are features) with an extra column containing clustering or cell type labels.  
    label_column: 
        A string to specify the name of cell type column in the data.
    n_features_to_select: 
        A int or float to control the number of features to select. If integer, the parameter is the absolute number of features to select. If float between 0 and 1, it is the fraction of features to select. Default is 0.15.
    ratio_of_none_zero_counts: 
        A float to determine the cutoff in which a feature will be omited when below specified value. A higher value leads to less features will be kept for calculation. 
    class_weight: 
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
    save: 
        A Bool value to decide whether the result will be written into disk. Default is True.
    logger: 
        A log_file object to write log information into disk. Default is None. 
    filename: 
        A string to name the output file. Default is None.
    n_jobs: 
        A int to decide the number of thread used for the program. Default is -1, meaning using all available threads.
    **kwargs: 
        Other paremeters passed to feature_selection function. 
    """
    res = dict()
    keys = data.loc[:, label_column].unique()
    for key in keys:
        wk_data = data.copy()
        wk_data.loc[:, label_column] = wk_data.loc[:, label_column].astype(str)
        wk_data.loc[wk_data[label_column] != key, label_column] = "OnevsRest"
        res[key] = find_signature_genes(data = wk_data, 
                                        label_column = label_column, 
                                        logger = logger, 
                                        save = False, 
                                        filename = "sub_cluster_{}".format(key), 
                                        class_weight = class_weight,
                                        ratio_of_none_zero_counts = ratio_of_none_zero_counts,
                                        n_features_to_select = n_features_to_select, 
                                        n_jobs = n_jobs,
                                        HVG = False,
                                        **kwargs)
        
    if save:
        if filename != None:
            with open("{}_{}.pkl".format(filename, record_time()), "wb") as file:
                pickle.dump(res, file)
        else:
            with open("markers_{}.pkl".format(filename, record_time()), "wb") as file:
                pickle.dump(res, file)

    return res

# %%
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
        A log_file object. Default is None. 
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

# %%
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

# %%
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

# %%
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

# %%
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

# %%
def get_gene_module(data, intersect_size = 10, intersect_group_size = 5, parallel_num = 8, lib_loc = None):
    """
    A function to find conserved gene modules across multiple gene sets.

    Parameters
    ----------
    data: 
        A data frame with each column representing a gene set list.
    intersection_size: 
        A int control the minimal intersection cutoff for adding a new gene set to the forming gene module.
    intersection_group_size: 
        A int control the minimal group size to consider for defining the first gene module. 
    parallel_num: 
        A int to control the number of thread used for the program.
    lib_loc:
        A string to specify the path where geneModule R package is installed. 
    """
    # convert pandas dataframe into R dataframe
    pandas2ri.activate()
    r_df = pandas2ri.py2rpy(data)
    
    geneModule = importr("geneModule", lib_loc = lib_loc)
    res = geneModule.get_gene_module(data = r_df, intersect_size = intersect_size, intersect_group_size = intersect_group_size, parallel_num = parallel_num)
    res = dict(zip(res.names, map(lambda x: x[0] if len(x) == 1 else x, res)))
    return res

# %%
def compare_gene_modules(module1, module2, api_key):
    """
    Compare two gene modules and analyze their similarities and differences.
    
    Parameters
    ----------
    module1:
        A gene set list representing the first gene module.
    module2:
        A gene set list representing the second gene module.
    api_key:
        A sting for the api key of the DeepSeek model.
    """
    # Genes are already converted to names
    module1_genes = module1
    module2_genes = module2

    # Find common and unique genes
    common_genes = set(module1_genes).intersection(set(module2_genes))
    unique_to_module1 = set(module1_genes) - set(module2_genes)
    unique_to_module2 = set(module2_genes) - set(module1_genes)

    # Create comparison template
    comparison_template = """You are a bioinformatics expert. Compare these two gene modules:
    
    Module 1: {module1_genes}
    Module 2: {module2_genes}
    
    Analyze:
    1. Common biological pathways between modules
    2. Unique pathways in each module
    3. Potential functional relationships between modules
    4. Disease associations shared between modules
    5. Tissue/cell type specificity differences
    
    Provide your analysis in clear, structured paragraphs."""

    # Create comparison prompt
    comparison_prompt = ChatPromptTemplate.from_template(comparison_template)

    # set model
    model = ChatOpenAI(
        model="deepseek-chat",
        temperature=0.7,
        openai_api_key=api_key,
        openai_api_base="https://api.deepseek.com/v1"
    )

    # Run comparison analysis
    comparison_chain = comparison_prompt | model | StrOutputParser()
    comparison_result = comparison_chain.invoke(
        {
            "module1_genes": ", ".join(module1_genes),
            "module2_genes": ", ".join(module2_genes),
        }
    )

    return {
        "common_genes": list(common_genes),
        "unique_to_module1": list(unique_to_module1),
        "unique_to_module2": list(unique_to_module2),
        "comparison_analysis": comparison_result,
    }

#%%
def analyse_one_gene_module(module_genes, api_key):
    """
    Analyze a single gene module using the DeepSeek model
    
    Parameters
    ----------
    module_genes:
        A gene set list representing the gene module.
    api_key:
        A sting for the api key of the DeepSeek model.
    """
    # Create prompt template
    template = """You are a bioinformatics expert. Analyze this list of genes and provide a detailed functional interpretation of the gene module:
    {module_genes}
    
    Consider:
    1. Common biological pathways
    2. Cellular processes involved
    3. Potential tissue/cell type specificity
    4. Disease associations
    5. Functional relationships between genes
    
    Provide your analysis in clear, structured paragraphs."""
    prompt = ChatPromptTemplate.from_template(template)

    # set model
    model = ChatOpenAI(
        model="deepseek-chat",
        temperature=0.7,
        openai_api_key=api_key,
        openai_api_base="https://api.deepseek.com/v1"
    )

    chain = prompt | model | StrOutputParser()

    # Run analysis for gene module
    analysis = chain.invoke({"module_genes": ", ".join(module_genes)})

    return analysis


