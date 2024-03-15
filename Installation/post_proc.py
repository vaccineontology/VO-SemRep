"""
Remove duplicate concepts from MRCONSO.RRF
"""

import pandas as pd

# Read MRCONSO.RRF into a DataFrame
# mrconso_df = pd.read_csv('MRCONSO.RRF', sep='|', header=None)
dupe_df = pd.read_csv('IDOBRU_VO_UMLS_meta/MRSTY_MERGED.RRF', sep='|', header=None)
# print(mrconso_df.columns)
# print(mrconso_df.head())

# # Define column names (adjust based on your file's structure)
# mrconso_df.columns = ['CUI', 'LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE', 'STR', 'SRL', 'SUPPRESS', 'CVF', 'DUM']
# # print(mrconso_df.head())

# dupe_df.columns = ['DUPE_CUI'] + list(range(1, dupe_df.shape[1]))

# If 'DUPE_CUI' is not the actual column name in your file, change it accordingly.

# Filter out rows from mrconso_df where 'CUI' matches any 'DUPE_CUI' in dupe_df
# filtered_mrconso_df = mrconso_df[~mrconso_df['CUI'].isin(dupe_df['DUPE_CUI'])]
filtered_mrconso_df = dupe_df.drop_duplicates(keep='first')

# Write the cleaned DataFrame back to an RRF file
filtered_mrconso_df.to_csv('IDOBRU_VO_UMLS_meta/FILTERED_MRSTY.RRF', sep='|', index=False, header=False)


# 13620680 MRCONSO