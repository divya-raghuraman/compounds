import requests
import pandas as pd
import urllib.parse
import os
import time

def find_similar_compounds(smile):
    smile = urllib.parse.quote(smile)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{smile}/cids/JSON?Threshold=85&MaxRecords=100"   
    response = requests.get(url)
    print(url)
    
    if response.status_code == 200:
        data = response.json()
        if "IdentifierList" in data and "CID" in data["IdentifierList"]:
            return data["IdentifierList"]["CID"]
        else:
            return []
    else:
        print(f"Failed to retrieve data for SMILES: {smile}")
        return []

    
def find_gene_data(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/GeneID/JSON'
    print(url)
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if "InformationList" in data and "Information" in data["InformationList"]:
            genes = data["InformationList"]["Information"][0].get("GeneID", [])
            if genes:
                return genes
            else:
                return ["No data available"]
        else:
            return ["No data available"]
    else:
        print(f"Failed to retrieve data for CID: {cid}")
        return ["No data available"]
    
    
def find_gene_data_info(gene_id):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/{gene_id}/summary/JSON'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "GeneSummaries" in data and "GeneSummary" in data["GeneSummaries"]:
            return data["GeneSummaries"]["GeneSummary"]
        else:
            return []
    

def get_compound_properties(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title,MolecularFormula,MolecularWeight/JSON'
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            properties = data["PropertyTable"]["Properties"][0]
            return properties
        else:
            return {}
    else:
        print(f"Failed to retrieve data for CID: {cid}")
        return {}

def get_compound_based_on_name(name):
    search_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON'
    search_response = requests.get(search_url)

    if search_response.status_code == 200:
        search_data = search_response.json()
        if "IdentifierList" in search_data and "CID" in search_data["IdentifierList"]:
            cid = search_data["IdentifierList"]["CID"][0]

            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON'
            response = requests.get(url)

            if response.status_code == 200:
                data = response.json()
                if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                    smiles_list = [item['CanonicalSMILES'] for item in data["PropertyTable"]["Properties"]]
                    return smiles_list
                else:
                    return []
            else:
                return []
        else:
            return []
    else:
        return []
    
def get_smiles(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON'
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            smiles_list = [item['CanonicalSMILES'] for item in data["PropertyTable"]["Properties"]]
            return smiles_list
        else:
            return []
    else:
        return []

def find_cid(smile):
    smile = urllib.parse.quote(smile)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/smiles/{smile}/cids/JSON?identity_type=same_isotope"
    response = requests.get(url)
    print(url)
    
    if response.status_code == 200:
        data = response.json()
        if "IdentifierList" in data and "CID" in data["IdentifierList"]:
            return data["IdentifierList"]["CID"][0]
        else:
            return None
    else:
        print(f"Failed to retrieve data for SMILES: {smile}")
        return None


drugs = [
    "Pridopidine",
    "Branaplam",
    "Riluzole",
    "Olanzapine",
    "Quetiapine",
    "Clonazepam",
    "Amantadine",
    "Lonafarnib",
    "Triheptanoin",
    "Edaravone"

]

csv_file = "similar_compounds.csv"
if not os.path.exists(csv_file):
    pd.DataFrame(columns=["Drug", "Molecular Formula", "Molecular Weight", "SMILES ID", "Associated Genes"]).to_csv(csv_file, index=False)

all_smiles = []

for name in drugs:
    smiles = get_compound_based_on_name(name)
    all_smiles.extend(smiles)

all_cids = []

for smile in all_smiles:
    main_cid = find_cid(smile)
    cids = find_similar_compounds(smile)
    all_cids.extend(cids)

print("All CIDs fetched:", len(all_cids))

all_cids.pop(0)

done = set()
count = 0

try:
    for cid in all_cids:
        if cid in done or cid is None:
            continue

        time.sleep(2)  
        properties = get_compound_properties(cid)
        genes = find_gene_data(cid)

        if properties:
            gene_symbols = []
            gene_names = []
            gene_summaries = []

            for gene in genes:
                gene_info = find_gene_data_info(gene)
                if gene_info:
                    gene_symbols.append(gene_info[0].get('Symbol', None))
                    synonyms = gene_info[0].get('Synonym', [])
                    gene_names.append(synonyms[0] if synonyms else None)
                    gene_summaries.append(gene_info[0].get('Description', None))

            gene_symbols = [symbol for symbol in gene_symbols if symbol]
            gene_names = [name for name in gene_names if name]
            gene_summaries = [summary for summary in gene_summaries if summary]

            if "No data available" not in genes and "No data available" not in gene_symbols and "No data available" not in gene_names and "No data available" not in gene_summaries:
                row = {
                    "Drug": properties.get("Title", "No data available"),
                    "Molecular Formula": properties.get("MolecularFormula", "No data available"),
                    "Molecular Weight": properties.get("MolecularWeight", "No data available"),
                    "SMILES ID": get_smiles(cid)[0],
                    "Associated Genes": ", ".join(gene_symbols)
                }

                pd.DataFrame([row]).to_csv(csv_file, mode='a', header=False, index=False)
                print(f"#{count + 1}: Added CID {cid} with all genes to CSV.")
                count += 1
                done.add(cid)
            else:
                print(f"CID {cid} skipped: No meaningful gene data available.")



except Exception as e:
    print(f"An error occurred: {e}")
    df = pd.read_csv(csv_file)
    df.to_excel("similar_compounds.xlsx", index=False)
    df_cleaned = df[~df.isin(["No data available"]).any(axis=1)]
    df_cleaned.to_excel("similar_compounds_cleaned.xlsx", index=False)

print("Data saved")
df = pd.read_csv(csv_file)
df.to_excel("similar_compounds.xlsx", index=False)

print("Data saved")
