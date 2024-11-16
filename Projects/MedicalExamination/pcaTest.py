import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score,precision_score, recall_score, f1_score, make_scorer
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt

survey_file = "./data/survey.csv"
survey_data = pd.read_csv(survey_file, delimiter=",")
data = survey_data.drop(survey_data.columns[0], axis=1)

# Set the file into x and y 
data_extracted = data.iloc[:,-2].str.split('_').str[0].rename('Difficulty')
data.iloc[:,-2] = data_extracted
x = data.iloc[:,:-1]
y = data.iloc[:,-1].rename("Result")

X_dropped = x.dropna()
y = y[X_dropped.index]
# Step 1: Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_dropped)

# Step 2: Apply PCA
pca = PCA()
X_pca = pca.fit_transform(X_scaled)

# Step 3: Determine the number of components to keep
explained_variance_ratio = np.cumsum(pca.explained_variance_ratio_)
n_components = np.argmax(explained_variance_ratio >= 0.95) + 1  # Keep enough components to retain 95% variance

# Plot the explained variance to see the number of components needed
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, marker='o')
plt.axhline(y=0.95, color='r', linestyle='--', label='95% Variance')
plt.xlabel("Number of Components")
plt.ylabel("Cumulative Explained Variance")
plt.legend()
plt.grid(True)
plt.show()

print(f"Number of components to retain 95% variance: {n_components}")

# Step 4: Transform the data using the selected number of components
pca = PCA(n_components=n_components)
X_reduced = pca.fit_transform(X_scaled)

# Step 5: Evaluate the performance with reduced dimensions
# Using RandomForest as an example classifier for cross-validation
clf = RandomForestClassifier(
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
scoring = make_scorer(roc_auc_score, needs_proba = True, multi_class = "ovr", average ="macro")
cv_scores = cross_val_score(clf, X_reduced, y, cv=5, scoring=scoring)

print(f"Cross-validated accuracy with {n_components} components: {cv_scores.mean():.4f}")

# Get loadings (contributions of each original feature to each principal component)
loadings = pd.DataFrame(pca.components_.T, columns=[f'PC{i+1}' for i in range(n_components)], index=x.columns)

# Show top questions contributing to each principal component
top_questions = loadings.abs().nlargest(5, 'PC1')  # Top 5 questions for the first principal component
print("Top questions contributing to the first principal component:\n", top_questions)

