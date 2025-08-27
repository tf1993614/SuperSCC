"""
Feature selection module for selecting data features.
"""

import re
import os
import pickle
import warnings
from datetime import datetime
from functools import reduce

from dcor import distance_correlation
import numpy as np
import pandas as pd
from scipy.stats import rankdata, gmean
from sklearn.metrics import accuracy_score, cohen_kappa_score, matthews_corrcoef, balanced_accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFECV, VarianceThreshold, SelectKBest, chi2, SelectFromModel, f_classif, mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder, MinMaxScaler

from ..utils import record_time, remove_time, check_duplicated_elements, retrieve_positive_markers, jaccard_score


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
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The "balanced" mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
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


def filtering_based_selection(X, y, variance_threshold="mean", mutual_info=False, chi_square_test=False, F_test=True, rank_method="dense", logger=None):
    """
    Perform filtering-based feature selection using variance threshold, chi-square test, F-test, and mutual information.
    
    Parameters
    ----------
    X : pandas.DataFrame
        Feature matrix (rows are samples, columns are features)
    y : pandas.Series
        Target labels
    variance_threshold : str
        Variance cutoff method - "zero" or "mean"
    mutual_info : bool
        Whether to use mutual information filtering
    chi_square_test : bool
        Whether to use chi-square test filtering
    F_test : bool
        Whether to use F-test filtering
    rank_method : str
        Ranking method for features
    logger : object
        Logger object for writing logs
        
    Returns
    -------
    dict
        Dictionary containing filtered features and rankings
    """
    message = "{} ======== filtering based selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    # Filter out by variance 
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

    # filter by chi square test
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

    # filter by mutual information
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
    
    return {
        'X_filtered': X_var,
        'y_filtered': y_var,
        'retained_features_ranking_by_variance': retained_features_ranking_by_variance,
        'retained_features_ranking_by_correlation': retained_features_ranking_by_correlation,
        'retained_features_by_filter': retained_features_by_filter
    }


def embedding_based_selection(X_var, y_var, model="svm", random_foreast_threshold=None, n_estimators=100, 
                             random_state=10, normalization_method="Min-Max", logistic_multi_class="ovr", 
                             linear_svm_multi_class="ovr", class_weight="balanced", rank_method="dense", 
                             cv=5, logger=None):
    """
    Perform embedding-based feature selection using random forest, logistic regression, and SVM.
    
    Parameters
    ----------
    X_var : pandas.DataFrame
        Filtered feature matrix from filtering step
    y_var : pandas.Series
        Target labels
    model : str
        Model type - "random_foreast", "logistic", or "svm"
    random_foreast_threshold : float
        Threshold for random forest feature importance
    n_estimators : int
        Number of trees for random forest
    random_state : int
        Random state for reproducibility
    normalization_method : str
        Normalization method - "Min-Max" or "Standardization"
    logistic_multi_class : str
        Multi-class strategy for logistic regression
    linear_svm_multi_class : str
        Multi-class strategy for SVM
    class_weight : str
        Class weight strategy
    rank_method : str
        Ranking method for features
    cv : int
        Cross-validation folds
    logger : object
        Logger object for writing logs
        
    Returns
    -------
    dict
        Dictionary containing embedding results and rankings
    """
    message = "{} ======== embedding based selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))

    embedding_positive_marker_selection = None
    
    # embedding-based on feature selection
    # select by random forest model
    if model == "random_foreast":
        RFC_ = RandomForestClassifier(n_estimators = n_estimators, random_state = random_state)
        # when random_foreast_threshold is None, 
        #  `1 / number of features` will be used as threshold
        if random_foreast_threshold == None:
            random_foreast_threshold = 1 / X_var.shape[1]
        RFC_embedding_selector = SelectFromModel(RFC_, threshold = random_foreast_threshold)
      
        X_RFC_embedding = RFC_embedding_selector.fit_transform(X_var, y_var)

        # evaluate the accuracy of classifier
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

        # evaluate the accuracy of classifier
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

        # evaluate the accuracy of classifier
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
    
    return {
        'retained_features_ranking_by_embedding': retained_features_ranking_by_embedding,
        'retained_features_by_embedding': retained_features_by_embedding,
        'embedding_positive_marker_selection': embedding_positive_marker_selection,
        'accuracy': accuracy,
        'X_stand': X_stand if model in ["logistic", "svm"] else X_var
    }


def wrapping_based_selection(X_var, y_var, X_stand, model="svm", n_features_to_select=0.15, 
                            step=100, cv=5, n_jobs=-1, n_estimators=100, random_state=10, 
                            logistic_multi_class="ovr", linear_svm_multi_class="ovr", 
                            class_weight="balanced", logger=None):
    """
    Perform wrapping-based feature selection using RFECV with specified model.
    
    Parameters
    ----------
    X_var : pandas.DataFrame
        Filtered feature matrix from filtering step
    y_var : pandas.Series
        Target labels
    X_stand : numpy.ndarray
        Standardized feature matrix (for logistic/svm models)
    model : str
        Model type - "random_foreast", "logistic", or "svm"
    n_features_to_select : int or float
        Number of features to select
    step : int or float
        Step size for RFECV
    cv : int
        Cross-validation folds
    n_jobs : int
        Number of parallel jobs
    n_estimators : int
        Number of trees for random forest
    random_state : int
        Random state for reproducibility
    logistic_multi_class : str
        Multi-class strategy for logistic regression
    linear_svm_multi_class : str
        Multi-class strategy for SVM
    class_weight : str
        Class weight strategy
    logger : object
        Logger object for writing logs
        
    Returns
    -------
    dict
        Dictionary containing wrapping results and rankings
    """
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
    
    return {
        'retained_features_ranking_by_wrapping': retained_features_ranking_by_wrapping,
        'retained_features_by_wrapping': retained_features_by_wrapping
    }


def aggregate_features(X_var, y_var, filtering_results, embedding_results, wrapping_results, 
                      merge_rank_method="geom.mean", n_features_to_select=0.15, logger=None):
    """
    Combine results from all feature selection methods using ensemble or intersection approach.
    
    Parameters
    ----------
    X_var : pandas.DataFrame
        Filtered feature matrix
    y_var : pandas.Series
        Target labels
    filtering_results : dict
        Results from filtering-based selection
    embedding_results : dict
        Results from embedding-based selection
    wrapping_results : dict
        Results from wrapping-based selection
    merge_rank_method : str
        Method to combine rankings - "geom.mean", "mean", "median", "max"
    n_features_to_select : int or float
        Number of features to select for ensemble method
    logger : object
        Logger object for writing logs
        
    Returns
    -------
    dict
        Dictionary containing final aggregated results
    """
    message = "{} ======== final feature selection ========".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    if isinstance(n_features_to_select, float):
        n_features_to_select = int(n_features_to_select * X_var.shape[1])
                
    # Features selection by ensemble mode
    rank_ls = list()

    for feature in X_var.columns:
        ranks = [filtering_results['retained_features_ranking_by_variance'][feature], 
                filtering_results['retained_features_ranking_by_correlation'][feature], 
                embedding_results['retained_features_ranking_by_embedding'][feature], 
                wrapping_results['retained_features_ranking_by_wrapping'][feature]]
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

    # Features selection by intersection mode
    retained_features_by_filter = set(filtering_results['retained_features_by_filter'])
    retained_features_by_embedding = set(embedding_results['retained_features_by_embedding'])
    retained_features_by_wrapping = set(wrapping_results['retained_features_by_wrapping'])

    final_feture_selection = reduce(lambda x,y: x.intersection(y), [retained_features_by_embedding, retained_features_by_filter, retained_features_by_wrapping])

    message2 = "* {} select top {} features ranked by using {} method to combine feature rankings obtained by different estimators".format(record_time(), len(rank_ls), merge_rank_method)
    message = "* {} {} features remained after intersecting the key features found by filtering, embedding and wrapping-based feature selection methods".format(record_time(), len(final_feture_selection))
    print(message2)
    print(message)
    if logger != None:
        logger.write("info", remove_time(message2))
        logger.write("info", remove_time(message))
    
    return {
        'final_feature_selection_by_ensemble': rank_ls,
        'final_feature_selection_by_intersection': final_feture_selection,
        'retained_features_by_filtering': retained_features_by_filter,
        'retained_features_by_embedding': retained_features_by_embedding,
        'retained_features_by_wrapping': retained_features_by_wrapping
    }


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
        A string to decide whether class weights will be considered. If None, all classes are supposed to have weight one. The "balanced" mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)). Default is 'balanced'.
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
    le = None
    if any(map(lambda x: isinstance(x, str), data[label_column])):
        le = LabelEncoder().fit(data[label_column])
        label = le.transform(data[label_column])
        data[label_column] = label
    
    # fix duplicated column names by adding ordinal suffix
    data.columns = check_duplicated_elements(data.columns.tolist())
    
    # get data and label
    X = data.iloc[:, 0:-1]
    y = data.iloc[:, -1]
    
    # step 2 - do feature selection
    message = "{} step 2 - do feature selection".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    # Step 2.1: Filtering-based selection
    filtering_results = filtering_based_selection(
        X, y, 
        variance_threshold=variance_threshold, 
        mutual_info=mutual_info, 
        chi_square_test=chi_square_test, 
        F_test=F_test, 
        rank_method=rank_method, 
        logger=logger
    )
    
    # Step 2.2: Embedding-based selection
    embedding_results = embedding_based_selection(
        filtering_results['X_filtered'], filtering_results['y_filtered'],
        model=model,
        random_foreast_threshold=random_foreast_threshold,
        n_estimators=n_estimators,
        random_state=random_state,
        normalization_method=normalization_method,
        logistic_multi_class=logistic_multi_class,
        linear_svm_multi_class=linear_svm_multi_class,
        class_weight=class_weight,
        rank_method=rank_method,
        cv=cv,
        logger=logger
    )
    
    # Step 2.3: Wrapping-based selection
    wrapping_results = wrapping_based_selection(
        filtering_results['X_filtered'], filtering_results['y_filtered'], 
        embedding_results['X_stand'],
        model=model,
        n_features_to_select=n_features_to_select,
        step=step,
        cv=cv,
        n_jobs=n_jobs,
        n_estimators=n_estimators,
        random_state=random_state,
        logistic_multi_class=logistic_multi_class,
        linear_svm_multi_class=linear_svm_multi_class,
        class_weight=class_weight,
        logger=logger
    )
    
    # Step 2.4: Aggregate features
    aggregation_results = aggregate_features(
        filtering_results['X_filtered'], filtering_results['y_filtered'],
        filtering_results, embedding_results, wrapping_results,
        merge_rank_method=merge_rank_method,
        n_features_to_select=n_features_to_select,
        logger=logger
    )
    
    # Build output dictionaries
    output1 = {
        "retained_features_ranking_by_variance": filtering_results['retained_features_ranking_by_variance'],
        "retained_features_ranking_by_correlation": filtering_results['retained_features_ranking_by_correlation'],
        "retained_features_ranking_by_embedding": embedding_results['retained_features_ranking_by_embedding'],
        "retained_features_ranking_by_wrapping": wrapping_results['retained_features_ranking_by_wrapping'],
        "final_feature_selection_by_ensemble": aggregation_results['final_feature_selection_by_ensemble'],
        "model_accuracy": embedding_results['accuracy'],
        "params_used_for_feature_selection": params
    }
    
    # only run when the number of unique labels is 2
    if len(set(filtering_results['y_filtered'])) == 2 and embedding_results['embedding_positive_marker_selection'] is not None:
        output1.update({"positive_marker_selection": embedding_results['embedding_positive_marker_selection']})
              
    output2 = {
        "retained_features_by_filtering": aggregation_results['retained_features_by_filtering'],
        "retained_features_by_embedding": aggregation_results['retained_features_by_embedding'],
        "retained_features_by_wrapping": aggregation_results['retained_features_by_wrapping'],
        "final_feature_selection_by_intersection": aggregation_results['final_feature_selection_by_intersection'],
        "model_accuracy": embedding_results['accuracy'],
        "params_used_for_feature_selection": params
    }
    
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

    if le is not None:
        output1.update({"label_classes": le.classes_})
    
    message = "{} finish feature selection".format(record_time())
    print(message)
    if logger != None:
        logger.write("info", remove_time(message))
    
    return output1


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


def pairwise_compare(data, focus = "expression1", ep_cut_off = 1, pct_cut_off = 0.5, only_positive = True):

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

                try:
                    # calculate jaccard index based on all positive genes between comparisons
                    jaccard_index["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] = jaccard_score(positive_genes1.feature, positive_genes2.feature)

                    # calculate distance correlation based on shared positive genes between comparisons
                    shared_positive_genes = set(positive_genes1.feature).intersection(positive_genes2.feature)
                    sub_data1_feature_ranking = positive_genes1.loc[positive_genes1.feature.isin(shared_positive_genes)].sort_values("feature").loc[:, focus].values
                    sub_data2_feature_ranking = positive_genes2.loc[positive_genes2.feature.isin(shared_positive_genes)].sort_values("feature").loc[:, focus].values

                    distance_cor["global cluster {} vs global cluster {}".format(keys[key1], keys[key2])] = distance_correlation(sub_data1_feature_ranking, sub_data2_feature_ranking)

                # when positive markers can't be extracted under the current cutoff
                # it would use all informative genes (containing both negative and positive genes)
                # for calculation
                except:
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
                    
    res = {"jaccard_index": jaccard_index, "distance_cor": distance_cor}
    
    return res
