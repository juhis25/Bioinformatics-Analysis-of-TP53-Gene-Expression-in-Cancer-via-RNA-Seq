import pandas as pd

# Load the TSV files
clinical_df = pd.read_csv('clinical.tsv', sep='\t')
exposure_df = pd.read_csv('exposure.tsv', sep='\t')
family_history_df = pd.read_csv('family_history.tsv', sep='\t')
follow_up_df = pd.read_csv('follow_up.tsv', sep='\t')
pathology_detail_df = pd.read_csv('pathology_detail.tsv', sep='\t')

# Merge the DataFrames on 'patient_id'
merged_df = clinical_df.merge(exposure_df, on='patient_id', how='left') \
                       .merge(family_history_df, on='patient_id', how='left') \
                       .merge(follow_up_df, on='patient_id', how='left') \
                       .merge(pathology_detail_df, on='patient_id', how='left')

# Save the merged DataFrame to a CSV file
merged_df.to_csv('merged_data.csv', index=False)

print("Merged data has been saved to merged_data.csv")
