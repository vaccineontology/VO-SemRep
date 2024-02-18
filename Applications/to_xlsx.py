import pandas as pd
import numpy as np
from certainty_estimator.predict_certainty import CertaintyEstimator
from collections import deque

def convert_semrep_output_to_excel(input_file, output_excel):

    lines = []
    cache_size = 10  # Adjust as needed depending on your file's structure
    line_cache = deque(maxlen=cache_size)

    # Open the SemRep output file
    with open(input_file, 'r') as infile:
        for line in infile:
            # Cache the cleaned current line
            clean_line = line.strip()  # Remove newline characters
            line_cache.append(clean_line)
            
            # Check for '|' in the line
            if "|" in clean_line:
                data = clean_line.split('|')  # Split the line into a list
                data_id = data[0]  # Extract the ID
                
                # Search the cache for the corresponding source line
                source_line = None
                for cached_line in reversed(line_cache):
                    if cached_line.startswith(data_id) and '|' not in cached_line:
                        source_line = cached_line
                        break
                
                # If the source line is found, combine it with the data
                if source_line is not None:
                    # Split the source line at the first space and take everything after
                    source_line = ' '.join(source_line.split(' ')[1:])
                    
                    # Add source line into new column
                    line_combined = [source_line] + data

                    # Add to lines
                    lines.append(line_combined)

    # Define column names
    column_names = ["Text", "ID", "SemRep Format", "Subject CUI",
                    "Subject Name", "Subject STypes", "Relation SType",
                    "Norm. Gene ID", "Norm. Gene Name", "Predicate",
                    "Object CUI", "Object Name", "Object STypes",
                    "Relation SType", "Norm. Gene ID", "Norm. Gene Name"]
    
    # Convert the list of lists into a pandas DataFrame
    df = pd.DataFrame(lines, columns=column_names)

    # Certainty-Level predict
    estimator_asp = CertaintyEstimator("aspect-level")
    asp = estimator_asp.predict(df["Text"])

    estimator_sent = CertaintyEstimator("sentence-level")
    sentence = estimator_sent.predict(df["Text"])
    
    df["Aspects-Level"] = asp
    df["Sentence-Level"] = sentence

    # Write the DataFrame to an Excel file
    df.to_excel(output_excel, index=False)

def main():
    # Example usage: *Hard-coded filenames
    input_file_path = 'VO_output/sciminer_med_VO.txt'
    output_excel_path = 'VO_output/sciminer_med_VO2.xlsx'

    convert_semrep_output_to_excel(input_file_path, output_excel_path)

if __name__ == "__main__":
    main()