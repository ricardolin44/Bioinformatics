
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

# GridSearchで最もいい候補を見つける
model = RandomForestClassifier(random_state=77)
param_grid = {
    'n_estimators': [200, 250, 300, 350],        # Number of trees in the forest
    'max_depth': [None, 2, 5, 10],        # Maximum depth of the tree
    'min_samples_split': [2, 3, 4],        # Minimum samples required to split a node
    'min_samples_leaf': [1, 2, 3],          # Minimum samples required at each leaf node
    'max_features': ['sqrt', 'log2', None]  # Number of features to consider for best split
}
gridSearchModel = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='roc_auc_ovr', n_jobs=-1, verbose=2)
gridSearchModel.fit(X_train, y_train)
print("Best Parameters:", gridSearchModel.best_params_)
print("Best Score:", gridSearchModel.best_score_)
best_rf_model = gridSearchModel.best_estimator_
y_pred = best_rf_model.predict_proba(X_test)
# accuracy = accuracy_score(y_test, y_pred)
# print("Test Set Accuracy:", accuracy)
auc = roc_auc_score(y_test_binarized, y_pred, average='macro', multi_class='ovr')
print("AUC (Macro):", auc)  

#Sort and select the feature with the most importance rate
features = X_test.columns
importance_df = pd.DataFrame({
    'features' : features,
    'importance' : best_rf_model.feature_importances_
})
importance_df = importance_df.sort_values(by='importance', ascending=False, ignore_index=True)
importance_df.loc[:50,'importance'].sum()
importance_df['cumulative_importance'] = \
    importance_df['importance'].cumsum()

# Set the threshold for cumulative importance and extract the top importance until threshold
# threshold = 0.35
# threshold_index = importance_df[importance_df['cumulative_importance'] > threshold].index[0]
# print("Current question number: "+len(selected_features))

question_range = range(5, 92)
auc_scores = []
# Re set our data and train it into model
for num_features in question_range:
    selected_features = importance_df['features'].iloc[:num_features].values

    X_data_filtered = x_data[selected_features]
    X_train_filtered, X_test_filtered, y_train, y_test = train_test_split(X_data_filtered, y_data, test_size=0.2, random_state=77)

    best_rf_model.fit(X_train_filtered, y_train)
    y_pred_best = best_rf_model.predict_proba(X_test_filtered)
    auc = roc_auc_score(y_test_binarized, y_pred_best, average='macro', multi_class='ovr')
    auc_scores.append(auc)

# show the importance into plot
plt.figure(figsize=(12, 6))
plt.bar(range(len(importance_df)), importance_df["cumulative_importance"]*100)
plt.xlabel("Data Point")
plt.ylabel("cumulative_importance")
plt.title("cumulative_importance of Each Data Point")
plt.show()

plt.figure(figsize=(10,6))
plt.plot(question_range, auc_scores, marker='o')
plt.title("AUC Score vs. Number of Questions")
plt.xlabel("Number of Questions")
plt.ylabel("AUC Score (Macro)")
plt.grid(True)
plt.show()

# Convert question_range and auc_scores to numpy arrays
x = np.array(question_range[7:30])
y = np.array(auc_scores[7:30])

m, b = np.polyfit(x, y, 1)
y_1streg = m * x + b

# Fit a polynomial regression (2nd degree is usually sufficient for such curves)
coefficients = np.polyfit(x, y, 2)
polynomial = np.poly1d(coefficients)

# Generate values for plotting the regression line
x_reg = np.linspace(x.min(), x.max(), 500)
y_reg = polynomial(x_reg)

# Find the x value (number of questions) that maximizes the fitted polynomial curve
optimal_questions = x_reg[np.argmax(y_reg)]
best_auc = max(y_reg)

# Plot the original data and the fitted curve
plt.figure(figsize=(10, 6))
plt.plot(x, y, 'o', label="Original Data")
plt.plot(x, y_1streg, '-', label="Regression Line (Polynomial)")
plt.axvline(optimal_questions, color='r', linestyle='--', label=f"Optimal Questions: {int(optimal_questions)}")
plt.title("AUC Score vs. Number of Questions")
plt.xlabel("Number of Questions")
plt.ylabel("AUC Score (Macro)")
plt.legend()
plt.grid(True)
plt.show()

print(f"Optimal number of questions for best AUC: {int(optimal_questions)}, Best AUC: {best_auc:.4f}")


# Using cross validation
scoring = make_scorer(roc_auc_score, needs_proba = True, multi_class = "ovr", average ="macro")
cv_mean_scores = []

question_range = range(5, 92)

cv_scores = cross_val_score(best_rf_model, x_data, y_data, cv=5, scoring=scoring)

print(f"Cross-validated AUC scores: {cv_scores}")
print(f"Mean AUC: {cv_scores.mean():.3f}")
print(f"Standard deviation of AUC: {cv_scores.std():.3f}")

for num_features in question_range:
    selected_features = importance_df['features'].iloc[:num_features].values

    x_data_selected = x_data[selected_features]
    cv_scores = cross_val_score(best_rf_model, x_data_selected, y_data, cv=5, scoring=scoring)
    cv_mean_scores.append(cv_scores.mean())

x = np.array(question_range[:30])
y = np.array(cv_mean_scores[:30])

m, b = np.polyfit(x, y, 1)
y_1streg = m * x + b

# Fit a polynomial regression (2nd degree is usually sufficient for such curves)
coefficients = np.polyfit(x, y, 2)
polynomial = np.poly1d(coefficients)

# Generate values for plotting the regression line
x_reg = np.linspace(x.min(), x.max(), 500)
y_reg = polynomial(x_reg)

# Find the x value (number of questions) that maximizes the fitted polynomial curve
optimal_questions = x_reg[np.argmax(y_reg)]
best_auc = max(y_reg)

plt.figure(figsize=(10, 6))
plt.plot(x, y, 'o', label="Original Data")
plt.plot(x, y_1streg, '-', label="Regression Line (Polynomial)")
plt.axvline(optimal_questions, color='r', linestyle='--', label=f"Optimal Questions: {int(optimal_questions)}")
plt.title("AUC Score vs. Number of Questions")
plt.xlabel("Number of Questions")
plt.ylabel("AUC Score (Macro)")
plt.legend()
plt.grid(True)
plt.show()