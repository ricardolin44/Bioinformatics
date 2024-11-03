
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, auc
from sklearn.model_selection import GridSearchCV
import pandas as pd
import matplotlib.pyplot as plt

# Import File
survey_file_terminal = "./Projects/MedicalExamination/data/survey.csv"
# survey_file = "./data/survey.csv"
survey_data = pd.read_csv(survey_file_terminal, delimiter=",")
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


# model = RandomForestClassifier(random_state=77)
# param_grid = {
#     'n_estimators': [200, 250, 300, 350],        # Number of trees in the forest
#     'max_depth': [None, 2, 5, 10],        # Maximum depth of the tree
#     'min_samples_split': [2, 3, 4],        # Minimum samples required to split a node
#     'min_samples_leaf': [1, 2, 3],          # Minimum samples required at each leaf node
#     'max_features': ['sqrt', 'log2', None]  # Number of features to consider for best split
# }
# model = RandomForestClassifier(random_state=77)
# gridSearchModel = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='accuracy', n_jobs=-1, verbose=2)
# gridSearchModel.fit(X_train, y_train)
# print("Best Parameters:", gridSearchModel.best_params_)
# print("Best Score:", gridSearchModel.best_score_)
# best_rf_model = gridSearchModel.best_estimator_
# y_pred = best_rf_model.predict(X_test)
# accuracy = accuracy_score(y_test, y_pred)
# print("Test Set Accuracy:", accuracy)

#Sort and select the feature with the most importance rate
features = X_test.columns
importance_df = pd.DataFrame({
    'features' : features,
    'importance' : firstModel.feature_importances_
})
importance_df = importance_df.sort_values(by='importance', ascending=False, ignore_index=True)
importance_df.loc[:50,'importance'].sum()
importance_df['cumulative_importance'] = \
    importance_df['importance'].cumsum()

# Set the threshold for cumulative importance and extract the top importance until threshold
threshold = 0.8
threshold_index = importance_df[importance_df['cumulative_importance'] > threshold].index[0]
selected_features = importance_df['features'].iloc[:threshold_index].values
# print("Current question number: "+len(selected_features))

# Re set our data and train it into model
X_data_filtered = x_data[selected_features]
X_train_filtered, X_test_filtered, y_train, y_test = train_test_split(X_data_filtered, y_data, test_size=0.2, random_state=77)

model = RandomForestClassifier(n_estimators=100, 
    random_state=77)
model.fit(X_train_filtered, y_train)
y_pred = model.predict(X_test_filtered)
accuracy = accuracy_score(y_test, y_pred)

# show the importance into plot
plt.figure(figsize=(12, 6))
plt.bar(range(len(importance_df)), importance_df["cumulative_importance"]*100)
plt.xlabel("Data Point")
plt.ylabel("cumulative_importance")
plt.title("cumulative_importance of Each Data Point")
plt.show()