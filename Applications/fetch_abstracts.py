import argparse
import time
from Bio import Entrez, Medline
import config as cfg


def extract_attributes(record):
    pmid = record.get("PMID", "?")
    title = record.get("TI", "?")
    abstract = record.get("AB", "?")

    return f"PMID- {pmid}\nTI  - {title}\nAB  - {abstract}\n\n"


def main(args):

    Entrez.email = cfg.email  # Replace with your email address

    # Specify your search term and other parameters
    search_term = f"{args.search_term}" 
    reldate_days = cfg.reldate_days  # 5 years 
    batch_size = cfg.batch_size 

    # Perform the search
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed",
            term=search_term,
            reldate=reldate_days,
            datetype="pdat",
            usehistory="y"
        )
    )

    count = int(search_results["Count"])
    print(f"Found {count} results")
    if args.max_res != -1:
        print(f"Fetching {args.max_res} results..")
        count = args.max_res

    # Open a file for writing the abstracts
    try:
        output_file = open(f"{args.output_file}", "w", encoding='utf-8')
    except IOError:
        print(f"Error opening file, {output_file}")
        return 1
    
    # Fetch records in batches and write to the file
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print(f"Going to download record {start + 1} to {end}")

        # Fetch records using Entrez.efetch
        stream = Entrez.efetch(
            db="pubmed",
            rettype="medline",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )

        # Parse the XML content
        records = Medline.parse(stream)
       
        # Iterate through records and write selected attributes to the output file
        for record in records:
            formatted_data = extract_attributes(record)
            output_file.write(formatted_data)
        stream.close()

        # Delay between requests
        time.sleep(3)

    # Close the output file
    output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='fetch abstracts from PubMed')
    parser.add_argument('output_file', help='FILE.txt')
    parser.add_argument('search_term', help='brucella')
    args = parser.parse_args()
    main(args)
