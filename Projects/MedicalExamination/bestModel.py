from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer, roc_auc_score
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
import pandas as pd

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

# Define your models
# Define your models with pipelines to handle missing values
models = {
    "Random Forest": Pipeline([
        ("imputer", SimpleImputer(strategy="mean")),
        ("classifier", RandomForestClassifier(n_estimators=200, max_depth=None, random_state=42, max_features='sqrt'))
    ]),
    "K-Nearest Neighbors": Pipeline([
        ("imputer", SimpleImputer(strategy="mean")),
        ("classifier", KNeighborsClassifier(n_neighbors=7))
    ]),
    "Support Vector Machine": Pipeline([
        ("imputer", SimpleImputer(strategy="mean")),
        ("classifier", SVC(kernel='linear', probability=True, random_state=42))
    ])
}

# Define the scoring metric as AUC
scoring = make_scorer(roc_auc_score, needs_proba=True, multi_class='ovr', average='macro')

# Dictionary to store cross-validation results
cv_results = {}

# Perform cross-validation for each model
for model_name, model in models.items():
    # Use 5-fold cross-validation for each model
    scores = cross_val_score(model, x, y, cv=5, scoring=scoring)
    cv_results[model_name] = scores
    print(f"{model_name} AUC scores: {scores}")
    print(f"{model_name} Mean AUC: {scores.mean():.3f} | Std Dev: {scores.std():.3f}\n")

# Summarize results
print("Summary of Model Performance (AUC):")
for model_name, scores in cv_results.items():
    print(f"{model_name}: Mean AUC = {scores.mean():.3f}, Std Dev = {scores.std():.3f}")
