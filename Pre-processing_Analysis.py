import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

# Load the merged dataset
df = pd.read_csv('merged_data.csv')

# Display basic information
print(df.head())
print(df.info())

# Check for missing values
print(df.isnull().sum())

# Handle missing values if any
df.dropna(inplace=True)

# Perform an ANOVA test to compare expression levels between cancerous and normal tissues
model = ols('expression_level ~ tissue_type', data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

print(anova_table)

# Visualize the expression levels using boxplots
plt.figure(figsize=(10, 6))
sns.boxplot(x='tissue_type', y='expression_level', data=df)
plt.title('TP53 Expression Levels in Cancerous vs. Normal Tissues')
plt.xlabel('Tissue Type')
plt.ylabel('Expression Level')
plt.show()

# Encode categorical variable
df['tissue_type'] = df['tissue_type'].map({'normal': 0, 'cancerous': 1})

# Split the data
X = df[['expression_level']]
y = df['tissue_type']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Train a Random Forest Classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

# Predict and evaluate
y_pred = clf.predict(X_test)
print('Accuracy:', accuracy_score(y_test, y_pred))
print(classification_report(y_test, y_pred))

# Feature importance
importances = clf.feature_importances_
print('Feature Importances:', importances)

# Generate a summary report
report = f"""
Statistical Analysis (ANOVA):
{anova_table}

Random Forest Classifier:
Accuracy: {accuracy_score(y_test, y_pred)}
Classification Report:
{classification_report(y_test, y_pred)}

Feature Importances: {importances}
"""

print(report)

# Save the report to a file
with open('tp53_analysis_report.txt', 'w') as f:
    f.write(report)
