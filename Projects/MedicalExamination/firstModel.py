from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score,precision_score, recall_score, f1_score, make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import label_binarize
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# Import File
# survey_file_terminal = "./Projects/MedicalExamination/data/survey.csv"
survey_file = "./data/survey.csv"
survey_data = pd.read_csv(survey_file, delimiter=",")
data = survey_data.drop(survey_data.columns[0], axis=1)

# Set the file into x and y 
data_extracted = data.iloc[:,-2].str.split('_').str[0].rename('Difficulty')
data.iloc[:,-2] = data_extracted
x_data = data.iloc[:,:-1]
y_data = data.iloc[:,-1].rename("Result")
# print(data.iloc[:, -2].unique())

# Split train and test data
X_train, X_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2, random_state=77)

# Set the model and train 
firstModel = RandomForestClassifier(n_estimators=100, random_state=77, n_jobs=-1)
firstModel.fit(X_train, y_train)
y_pred = firstModel.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)

cm = confusion_matrix(y_test, y_pred)
print("混同行列:\n", cm)

precision = precision_score(y_test, y_pred, average='macro')
recall = recall_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')

print("Precision:", precision)
print("Recall:", recall)
print("F1スコア:", f1)

y_test.unique()
# ラベルのバイナライズ（マルチクラス対応）
label_mapping = {'Cluster1': 0, 'Cluster2': 1, 'Cluster3': 2, 'Cluster4': 3, 'Cluster5': 4, 'Cluster6': 5, 'Cluster7': 6}
y_test_numeric = y_test.map(label_mapping)
y_test_numeric.unique()
y_test_binarized = label_binarize(y_test_numeric, classes=[0,1,2,3,4,5,6])
y_score = firstModel.predict_proba(X_test)

# マクロ平均のAUC
auc = roc_auc_score(y_test_binarized, y_score, average='macro', multi_class='ovr')
print("AUC (Macro):", auc)