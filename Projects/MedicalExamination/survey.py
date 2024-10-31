#%%
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, auc
import pandas as pd
import matplotlib.pyplot as plt

survey_file = "./data/survey.csv"
survey_data = pd.read_csv(survey_file, delimiter=",")
data = survey_data.drop(survey_data.columns[0], axis=1)

x_data = data.iloc[:,:-2]
y_data = data.iloc[:,-1].rename("Result")

#%%
X_train, X_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2, random_state=42)

model = RandomForestClassifier(n_estimators=100, random_state=77)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)

#%%
features = X_test.columns
importance_df = pd.DataFrame({
    'features' : features,
    'importance' : model.feature_importances_
})

importance_df = importance_df.sort_values(by='importance', ascending=False, ignore_index=True)
importance_df.loc[:50,'importance'].sum()

importance_df['cumulative_importance'] = \
    importance_df['importance'].cumsum()

threshold = 0.8
threshold_index = importance_df[importance_df['cumulative_importance'] > threshold].index[0]

selected_features = importance_df['features'].iloc[:threshold_index].values
print(len(selected_features))

plt.figure(figsize=(12, 6))
plt.bar(range(len(importance_df)), importance_df["cumulative_importance"]*100)
plt.xlabel("Data Point")
plt.ylabel("cumulative_importance")
plt.title("cumulative_importance of Each Data Point")
plt.show()