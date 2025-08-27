"""
Label transfer module for functions that transfer labels between datasets.
"""

import os
import re
import pickle
import numpy as np
import pandas as pd
import magic
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV

from ..utils import record_time, remove_time, load_pickle_file


def model_training(data,
                   label_column,
                   features,
                   model,
                   normalization_method = "Min-Max",
                   feature_selection_params = None,
                   parameters = None,
                   probability = False,
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
    feature_selection_params: 
        A dict to show params used for extracting informaive features used for model training.
    parameters: 
        A dict to elucidate the parameter name and value pair that should be searched in the grid search algorithm for corresponding model. Default is None.
    probability:
        A Bool value to decide whether the logits should be return when the model is "svm". Default is False. 
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
        
        SVM_ = SVC(decision_function_shape = "ovr", probability = probability)
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


def predict_label(query_data, 
                  models, 
                  wk_dir = os.getcwd(), 
                  normalization_method = "Min-Max",
                  magic_based_imputation = False,
                  pred_confidence_cutoff = None, 
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
    magic_based_imputation:
        A Bool value to decide whether MAGIC-based imputation should be employed. When False, the reference features lost in the query would be padded with 0. Default is False. 
    pred_confidence_cutoff:
        A float to decide the prediction accuracy cutoff. If the prediction score is lower than the cut off, the predicted cell label will return 'uncertain'. Default is None. This argument  
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


            if magic_based_imputation:
                magic_operator = magic.MAGIC()
                X_magic = magic_operator.fit_transform(query_data, genes = list(without_features))
                query_data = pd.concat([query_data, X_magic], axis = 1)
            else:
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
                if pred_confidence_cutoff != None:
                    prediction_prob = pd.DataFrame(prediction_prob)
                    prediction_prob.columns = model_class
                    pred_score = list()
                    label = list()
                    for i in range(prediction_prob.shape[0]):
                        max_value_index = prediction_prob.iloc[i, :].argmax()
                        pred_score.append(prediction_prob.iloc[i, :][max_value_index])
                        label.append(prediction_prob.iloc[i, :].index[max_value_index])
                    pred_df = pd.DataFrame({"label": label, "pred_score": pred_score}, index = range(len(label)))
                    for idx, df in pred_df.iterrows():
                        if df.pred_score < pred_confidence_cutoff:
                            pred_df.loc[idx, "label"] = "uncertain"
            except:
                pass
        except:
            prediction = model.predict(np.array(query_data))
            try:
                if pred_confidence_cutoff != None:
                    prediction_prob = pd.DataFrame(prediction_prob)
                    prediction_prob.columns = model_class
                    pred_score = list()
                    label = list()
                    for i in range(prediction_prob.shape[0]):
                        max_value_index = prediction_prob.iloc[i, :].argmax()
                        pred_score.append(prediction_prob.iloc[i, :][max_value_index])
                        label.append(prediction_prob.iloc[i, :].index[max_value_index])
                    pred_df = pd.DataFrame({"label": label, "pred_score": pred_score}, index = range(len(label)))
                    for idx, df in pred_df.iterrows():
                        if df.pred_score < pred_confidence_cutoff:
                            pred_df.loc[idx, "label"] = "uncertain"
            except:
                pass
        
        try:
            prediction_res[model_name] = {"prediction": [model_class[i] for i in prediction],
                                          "prediction_prob": prediction_prob}
        except:
            prediction_res[model_name] = {"prediction": [model_class[i] for i in prediction]}

        if pred_confidence_cutoff != None:
            prediction_res[model_name]["prediction"] = pred_df.label.tolist()
        
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
