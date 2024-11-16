
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
    'n_estimators': [200, 250, 300, 350],
    'max_depth': [None, 5, 10, 15],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['sqrt', 'log2', None],
    'bootstrap': [True, False],
    'max_samples': [0.5, 0.75, 1.0],       # Only when bootstrap=True
    'criterion': ['gini', 'entropy'],
    'min_impurity_decrease': [0.0, 0.01, 0.1],
    'class_weight': [None, 'balanced', 'balanced_subsample'],
    'oob_score': [True, False]              # Only when bootstrap=True
}
gridSearchModel = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='roc_auc_ovr', n_jobs=-1, verbose=2)
gridSearchModel.fit(x, y)
print("Best Parameters:", gridSearchModel.best_params_)
print("Best Score:", gridSearchModel.best_score_)

bestModel = RandomForestClassifier(
    random_state=77, 
    n_estimators=200, 
    max_depth=None, 
    bootstrap=True,
    class_weight='balanced',
    criterion='entropy',
    max_features=None,
    max_samples=0.75,
    min_impurity_decrease=0.0,
    min_samples_leaf=1,
    min_samples_split=2,
    oob_score=True
)
bestModel.fit(x, y)
scoring = make_scorer(roc_auc_score, needs_proba = True, multi_class = "ovr", average ="macro")
cv_scores = cross_val_score(bestModel, x, y, cv=5, scoring=scoring)
print("Mean Score of Best Model: ",round(cv_scores.mean(),3))

#Sort and select the feature with the most importance rate
features = x.columns
importance_df = pd.DataFrame({
    'features' : features,
    'importance' : bestModel.feature_importances_
})
importance_df = importance_df.sort_values(by='importance', ascending=False, ignore_index=True)
# importance_df.loc[:50,'importance'].sum()
importance_df['cumulative_importance'] = \
    importance_df['importance'].cumsum()


# show the importance into plot
plt.figure(figsize=(12, 6))
plt.bar(range(len(importance_df)), importance_df["cumulative_importance"]*100)
plt.xlabel("Data Point")
plt.ylabel("Cumulative Importance")
plt.title("Cumulative Importance of Sorted Data Point")
plt.show()

plt.figure(figsize=(12, 6))
plt.bar(importance_df["features"], importance_df["importance"])
plt.xticks(np.arange(0, 93, 10))  # Sets ticks at intervals of 10 between 0 and 50
plt.xlabel("Data Point")
plt.ylabel("Importance")
plt.title("Importance of Each Data Point")
plt.show()


# Optimized using best random forest model
cv_mean_scores = []
question_range = range(1,94)

for num_features in question_range:
    print(num_features)
    selected_features = importance_df['features'].iloc[:num_features].values

    x_selected = x[selected_features]
    cv_scores = cross_val_score(bestModel, x_selected, y, cv=5, scoring=scoring)
    cv_mean_scores.append(cv_scores.mean())

# x = np.array(question_range)
# y = np.array(cv_mean_scores)

# m, b = np.polyfit(x, y, 1)
# y_1streg = m * x + b

# # Fit a polynomial regression (2nd degree is usually sufficient for such curves)
# coefficients = np.polyfit(x, y, 2)
# polynomial = np.poly1d(coefficients)
# # Generate values for plotting the regression line
# x_reg = np.linspace(x.min(), x.max(), 500)
# y_reg = polynomial(x_reg)
# # Find the x value (number of questions) that maximizes the fitted polynomial curve
# optimal_questions = x_reg[np.argmax(y_reg)]
# best_auc = max(y_reg)


# plt.figure(figsize=(10, 6))
# plt.plot(x, y, 'o')
# plt.plot(x_reg, y_reg, '-', label="Regression Line (Polynomial)")
# plt.axvline(optimal_questions, color='r', linestyle='--', label=f"Optimal Questions: {int(optimal_questions)}")
# plt.title("AUC Score vs. Number of Questions (Sorted Based on Importance)")
# plt.xlabel("Number of Questions")
# plt.ylabel("AUC Score (Macro)")
# plt.legend()
# plt.grid(True)
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

def plot_auc_vs_questions(question_range, auc_scores, type='linear', degree=2, start=0, end=95):
    # Convert question_range and auc_scores to numpy arrays and slice based on input range
    x = np.array(question_range[start:end])
    y = np.array(auc_scores[start:end])

    # Fit a higher-degree polynomial regression
    coefficients = np.polyfit(x, y, degree)
    polynomial = np.poly1d(coefficients)

    # Generate values for plotting the polynomial regression line
    x_reg = np.linspace(x.min(), x.max(), 500)
    y_poly_reg = polynomial(x_reg)

    # Find the x value (number of questions) that maximizes the fitted polynomial curve
    optimal_questions = x_reg[np.argmax(y_poly_reg)]
    best_auc = max(y_poly_reg)

    # Linear Regression
    m, b = np.polyfit(x, y, 1)
    y_1streg = m * x + b

    # Define the exponential plateau function
    def plateau_func(x, a, b):
        return a * (1 - np.exp(-b * x))

    # Fit the data to the plateau function
    popt, _ = curve_fit(plateau_func, x, y, maxfev=10000)
    c, d = popt

    x_fit = np.linspace(x.min(), x.max(), 500)
    y_fit = plateau_func(x_fit, c, d)

    # Plot the original data and the polynomial regression curve
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, 'o', label="Original Data")
    if (type == "linear"):
        plt.plot(x, y_1streg, "-", label="Linear Regression")
        equation_text = f"y = {m:.2f}x + {b:.2f}"
        plt.text(0.05, 0.95, equation_text, transform=plt.gca().transAxes, 
             fontsize=12, color='orange', verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    elif (type == 'polynomial') :
        plt.plot(x_reg, y_poly_reg, '-', label=f"Polynomial Regression (Degree {degree})")
        plt.axvline(optimal_questions, color='r', linestyle='--', label=f"Optimal Questions: {int(optimal_questions)}")
    else:
        plt.plot(x_fit, y_fit, '-', label=f"Exponential Plateau Model: y = {c:.2f} * (1 - e^(-{d:.2f} * x))", color='orange')
    plt.title("AUC Score vs. Number of Questions (Sorted Based on Importance)")
    plt.xlabel("Number of Questions")
    plt.ylabel("AUC Score (Macro)")
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"Optimal number of questions for best AUC: {int(optimal_questions)}, Best AUC: {best_auc:.4f}")

plot_auc_vs_questions(question_range, cv_mean_scores, type='pol', start=4)