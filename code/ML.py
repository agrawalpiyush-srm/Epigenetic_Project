##### Python code for Implementing Machine Learning using Sklearn Library ########

from numpy import arange
from pandas import read_csv
from sklearn.model_selection import RepeatedKFold
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier

###### Loading datasets ########

# Load train dataset
train_url = '/Users/agrawalp4/Documents/Varun_work/ml/MAIN_files/Revised_ML_files_RH4/all/train.csv'
dataframe_train = pd.read_csv(train_url, header=None)
data_train = dataframe_train.values
X_train, y_train = data_train[:, :-1], data_train[:, -1]

# Load test dataset

test_url = '/Users/agrawalp4/Documents/Varun_work/ml/MAIN_files/Revised_ML_files_RH4/all/test.csv'
dataframe_test = pd.read_csv(test_url, header=None)
data_test = dataframe_test.values
X_test, y_test = data_test[:, :-1], data_test[:, -1]


####### Support Vector Classifier ####

clf = SVC(kernel='rbf', C=1, gamma=0.001, probability=True)
clf.fit(X_train, y_train)
np.random.seed(786)

# make a prediction
y_pred1 = clf.predict_proba(X_test)

score_res = roc_auc_score(y_test,y_pred1[:,1])          #### ROC using Probability scores
print("roc_auc on test data : %f" % score_res)


####### Random Forest Classifier ####
import random
from sklearn.ensemble import RandomForestClassifier
np.random.seed(786)

clf = RandomForestClassifier(n_estimators=500)
clf.fit(X_train, y_train)

# make a prediction
y_pred1 = clf.predict_proba(X_test)

score_res = roc_auc_score(y_test,y_pred1[:,1])          #### ROC using Probability scores
print("roc_auc on test data : %f" % score_res)


####### Gradient Boosting Classifier ####
from sklearn.ensemble import GradientBoostingClassifier

clf = GradientBoostingClassifier(n_estimators=500, learning_rate=0.1, max_depth=1, random_state=0)
clf.fit(X_train, y_train)

# make a prediction
#y_pred = clf.predict(X_test)
y_pred1 = clf.predict_proba(X_test)
# summarize prediction

#print(y_pred1)

score_res = roc_auc_score(y_test,y_pred1[:,1])          #### ROC using Probability scores
print("roc_auc on test data : %f" % score_res)
