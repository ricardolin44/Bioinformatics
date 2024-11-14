
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
x = data.iloc[:,:-1]
y = data.iloc[:,-1].rename("Result")
# print(data.iloc[:, -2].unique())

# GridSearchで最もいい候補を見つける
model = RandomForestClassifier(random_state=77)
param_grid = {
    'classifier__n_estimators': [200, 250, 300, 350],
    'classifier__max_depth': [None, 5, 10, 15],
    'classifier__min_samples_split': [2, 5, 10],
    'classifier__min_samples_leaf': [1, 2, 4],
    'classifier__max_features': ['sqrt', 'log2', None],
    'classifier__bootstrap': [True, False],
    'classifier__max_samples': [0.5, 0.75, 1.0],       # Only when bootstrap=True
    'classifier__criterion': ['gini', 'entropy'],
    'classifier__min_impurity_decrease': [0.0, 0.01, 0.1],
    'classifier__class_weight': [None, 'balanced', 'balanced_subsample'],
    'classifier__oob_score': [True, False]              # Only when bootstrap=True
}
gridSearchModel = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='roc_auc_ovr', n_jobs=-1, verbose=2)
gridSearchModel.fit(x, y)
print("Best Parameters:", gridSearchModel.best_params_)
print("Best Score:", gridSearchModel.best_score_)

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


# show the importance into plot
plt.figure(figsize=(12, 6))
plt.bar(range(len(importance_df)), importance_df["cumulative_importance"]*100)
plt.xlabel("Data Point")
plt.ylabel("cumulative_importance")
plt.title("cumulative_importance of Each Data Point")
plt.show()


# Optimized using best random forest model
cv_mean_scores = []
question_range = range(5, 92)

scoring = make_scorer(roc_auc_score, needs_proba = True, multi_class = "ovr", average ="macro")

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