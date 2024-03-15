import pandas as pd

folder_name = 'IDOBRU_VO_UMLS_meta'
# umls-2023AB
input_merge = 'VO_UMLS_meta'

# MERGE MRCONSO.RRF
print('Merging MRCONSO.RRF...')
vo_conso = pd.read_csv('MRCONSO.RRF', sep='|', header=None)
umls_conso = pd.read_csv(input_merge + '/MRCONSO_MERGED.RRF', sep='|', header=None)

concat_conso = pd.concat([vo_conso, umls_conso], ignore_index=True)

concat_conso.to_csv(folder_name + '/MRCONSO_MERGED.RRF', sep='|', index=False, header=False)

# MERGE MRSTY.RRF
print('Merging MRSTY.RRF...')
vo_sty = pd.read_csv('MRSTY.RRF', sep='|', header=None)
umls_sty = pd.read_csv(input_merge + '/MRSTY_MERGED.RRF', sep='|', header=None)

concat_sty = pd.concat([vo_sty, umls_sty], ignore_index=True)

concat_sty.to_csv(folder_name + '/MRSTY_MERGED.RRF', sep='|', index=False, header=False)

# # MERGE MRSAB.RRF
# print('Merging MRSAB.RRF...')
# vo_sab = pd.read_csv('MRSAB.RRF', sep='|', header=None)
# umls_sab = pd.read_csv(input_merge + '/MRSAB_MERGED.RRF', sep='|', header=None)

# concat_sab = pd.concat([vo_sab, umls_sab], ignore_index=True)

# concat_sab.to_csv(folder_name + '/MRSAB_MERGED.RRF', sep='|', index=False, header=False)

# # MERGE MRRANK.RRF
# print('Merging MRRANK.RRF...')
# vo_rank = pd.read_csv('MRRANK.RRF', sep='|', header=None)
# umls_rank = pd.read_csv(input_merge + '/MRRANK_MERGED.RRF', sep='|', header=None)

# # need to adjust the first column
# new_pred = range(len(vo_rank) + len(umls_rank), 0, -1)
# concat_rank = pd.concat([vo_rank, umls_rank], ignore_index=True)

# concat_rank[0] = [str(i).zfill(4) for i in new_pred]    

# concat_rank.to_csv(folder_name + '/MRRANK_MERGED.RRF', sep='|', index=False, header=False)



