import argparse
import time
from Bio import Entrez, Medline
import subprocess

"""
Ideal goal: 
    Fetch for keywords:
        Brucella (1464) // Brucella vaccine (217)
        E. coli (45,821) (?) // E. coli vaccine (1165)
        Protein vaccine [2018-2024] (6536) // 

ensure that metamap servers are running
"""


def import_ids(filename):
    """In case user wants to import IDs manually.
        Though, this is implemented if the textfile in 
        formatted in complete PubMed data.
        Example: Each abstract has more than just PMID and TI,
        so it only tries to extract the PMID value
    """
    # 'filename.txt' should be replaced with your file name
    with open(filename, 'r') as f:  
        data = f.readlines()

    ids = []
    for line in data:
        if line.startswith('PMID-'):
            id = line.replace('PMID-', '').strip()
            ids.append(id)

    print(ids)
    return ids


def extract_attributes(record):
    """Only store PMID, Title, and Abstract in
        the input SemRep textfile.
    """
    pmid = record.get("PMID", "?")
    title = record.get("TI", "?")
    abstract = record.get("AB", "?")

    return f"PMID- {pmid}\nTI  - {title}\nAB  - {abstract}\n\n"


def append_file(source_file, destination_file):
    """Used to append from output of one 
        abstract to a collection of SemRep output.
    """
    # Reads result from one abstract
    with open(source_file, 'r') as source:
        data = source.read()

    # TODO: Writes to excel file

    # Appends to the whole result 
    with open(destination_file, 'a') as destination:
        destination.write(f"\n{data}")


def main(args):

    Entrez.email = args.email  # Replace with your email address

    # METHOD 1: Query for abstract lists
    # if args.method.lower() == "search":
    # Specify your search term and other parameters
    search_term = input("Search term: ").lower() # Ex: "brucella vaccine"
    try:
        input_year = int(input("Enter the relative relevant year (1,2,5,etc.): "))
        if input_year <= 0:
            print("Please enter a positive integer.")
            exit(1)
    except ValueError:
        print("Invalid input. Please enter a whole number.")
        exit(1)
    reldate_days = 365*input_year 
    batch_size = 10 # Hard-coded value
    retmax = 1000 # Hard-coded value
    # Perform the search
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed",
            term=search_term,
            retmax=retmax,
            reldate=reldate_days,
            datetype="pdat",
            usehistory="y"
        )
    )
    count = int(search_results["Count"])
    print(f"Found {count} results")
    # if args.max_res != -1:
    #     print(f"Fetching {args.max_res} results..")
    #     count = args.max_res

    # List of Pubmed IDS
    pubmed_ids = search_results["IdList"]

    # # METHOD 2: Import ID from a list
    # if args.method.lower() == "import":
    #     filename = input("Input textfile name (.txt): ")
    #     pubmed_ids = import_ids(filename)
    #     count = len(pubmed_ids)

    # Open a file for writing the abstracts
    # try:
    #     output_file = open("input_abstract.txt", "w", encoding='utf-8')
    # except IOError:
    #     print(f"Error opening file, {output_file}")
    #     return 1
    
    # Reset test_output file
    with open("test_output.txt", 'w') as i:
        print(f"Resetting {i} to empty file...")

    failed = []
    # print(pubmed_ids)
    # Fetch records in batches and write to the file
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print(f"Going to download record {start + 1} to {end}")

        # Fetch records using Entrez.efetch
        stream = Entrez.efetch(
            db="pubmed",
            rettype="medline",
            retmode="text",
            retmax=batch_size,
            id=pubmed_ids[start:end]
        )

        # Parse the XML content
        records = Medline.parse(stream)
       
        # Iterate through records and write selected attributes to the output file
        counter = start
        for record in records:
        # print(records)
        # item = [record for record in records]
        # print(item)
            print(f"Processing abstract {counter} [{pubmed_ids[counter]}]...")
            formatted_data = extract_attributes(record)
            # Write to input_abstract.txt
            with open("input_abstract.txt", "w", encoding='utf-8') as f:
                f.write(formatted_data)

            # Run SemRep with one abstract at a time
            try:
                # replace 'your_command_for_semrep' with actual command to run SemRep
                # hard-coded command
                command = 'semrep.v1.8 -M -Z 2022AB -V umls_vo1 input_abstract.txt dummy_output.txt'
                
                result = subprocess.run(
                    command,
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    text=True,
                    timeout=20)

                cmd_out = result.stdout
                # print(cmd_out)

                # check if result is not an error
                if 'ERROR' not in cmd_out:
                    # assuming output.txt is the file to which you want to append the result
                    append_file("dummy_output.txt", "test_output.txt")
                    # output_file.write(result)
                else: # Found: Never execute because SemRep error message is not included in cmd_out
                    print('Error in SemRep execution for abstract: ', pubmed_ids[start])

            except subprocess.CalledProcessError as e:
                failed.append(pubmed_ids[start])
                print('Unexpected error for abstract: ', formatted_data)
                print('Error details: ', e)
                # output_file.write(formatted_data)
            except subprocess.TimeoutExpired as e:
                print(f"Timeout expired. The command took longer than {20} seconds to complete.")
                print(f"Will write incomplete output to file.")
                append_file("dummy_output.txt", "test_output.txt")
            counter += 1
        stream.close()

        # Delay between requests
        time.sleep(1.5)

    # Close the output file
    # output_file.close()

    print(f"{len(failed)} of total {retmax} failed to parse.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='fetch abstracts from PubMed and run SemRep')
    # parser.add_argument('method', help='search/import')
    parser.add_argument('email', help='name@email.com')
    args = parser.parse_args()
    main(args)
