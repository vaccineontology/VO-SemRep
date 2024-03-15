"""
Script that finds duplicate concepts from VO in UMLS.
"""

import pip._vendor.requests as requests
import json
from pprint import pprint
from rdflib import ConjunctiveGraph
import umls_references as umls
import pandas as pd
import math


def vo_umls_refs():
    graph = ConjunctiveGraph()
    graph.parse('vo.owl')

    mondo_umls_refs = umls.test_vo(graph)
    # pprint(response)

    pprint(mondo_umls_refs)

    mrsty_df = pd.read_csv('umls-2023AB/MRSTY.RRF', sep='|', header=None)
    output = {}
    # pprint(mrsty_df.head())
    file = open('CUI_DUPE_MYSTY.RRF', 'w', encoding='utf-8')
    for item in mondo_umls_refs:
        result = mrsty_df[mrsty_df[0] == item[1]]
        # pprint(result[3].to_string(index=False, dtype=False))
        # pprint(item[0])
        # pprint(result)
        if item[0] not in output:
            output[item[0]] = (item[1], result[1].astype(str).values[0], result[2].astype(str).values[0], result[3].astype(str).values[0])
            file.write(item[0] + '|' + result[1].astype(str).values[0] + '|' 
                       + result[2].astype(str).values[0] + '|' + result[3].astype(str).values[0] 
                       + '|' + result[4].astype(str).values[0] + '|' 
                       + ('' if math.isnan(result[5].values[0]) 
                          else str(result[5].astype(int).values[0])) + '|\n')

        # pprint(output[item[0]])
    file.close()

    return output
# # API APPROACH
# file = open('output_test.txt', 'w', encoding='utf-8')
# pprint("GENERATING...")
# SEM_TYPES = {}
# for item in mondo_umls_refs:
#     res = requests.get("https://uts-ws.nlm.nih.gov/rest/content/current/CUI/" + item[1] + "?apiKey=fabafb43-2caf-49b2-94bb-a27a5eb3b105")
#     response = res.json()['result']
#     semtypes = response['semanticTypes']

#     if semtypes[0]['name'] not in SEM_TYPES:
#         sem_type = requests.get("%s" % semtypes[0]['uri'])
#         sem_type_info = sem_type.json()['result']
#         SEM_TYPES[semtypes[0]['name']] = (sem_type_info['ui'], sem_type_info['treeNumber'], sem_type_info['name'], sem_type_info['abbreviation'])

#     sem_types = SEM_TYPES[semtypes[0]['name']]
#     item.extend([semtypes[0]['name'], sem_types[0], sem_types[1], sem_types[2], sem_types[3]])
#     file.write("|".join(item))
# file.close()

# pprint(mondo_umls_refs)
    # pprint(semtypes[0]['uri'])
    # pprint(sem_type_info)
    # pprint(response)
#     file.write("[" + item[0] + ": " + item[1] + " | " + semtypes[0]['name'] + " | " + sem_type_info['abbreviation']  + " | " + sem_type_info['treeNumber']  + " | " + sem_type_info['ui']+ "]\n")
# file.close()

if __name__ == "__main__":
    vo_umls_refs()


