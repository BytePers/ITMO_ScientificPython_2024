#FIRST TASK
import re

import requests, sys, json

#Function 1 “get_<db_name>” that accepts list of IDs and outputs API request.
def get_Uniprot(accessions_list: list):
  return requests.get("https://rest.uniprot.org/uniprotkb/accessions", params={'accessions': accessions_list})

#Function 2 “parse_response_<db_name>” that accepts response and outputs parsed data.
def parse_response_Uniprot(resp = None):
    resp = resp.json()
    resp = resp["results"]
    output = {}
    for val in resp:
        acc = val['primaryAccession']
        species = val['organism']['scientificName']
        gene = val['genes']
        seq = val['sequence']
        output[acc] = {'organism':species, 'geneInfo':gene, 'sequenceInfo':seq, 'type':'protein'}

    return output

#Function 1 “get_<db_name>” that accepts list of IDs and outputs API request.
def get_ENSEMBL(accessions_list):

    json_ids = json.dumps({"ids": accessions_list})
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.post(server + ext, headers=headers, data=f'{json_ids}')

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    return decoded

#Function 2 “parse_response_<db_name>” that accepts response and outputs parsed data.
def parse_response_ENSEMBL(decoded):
    output = {}
    for val in decoded:
        acc = decoded[val]['id']
        species = decoded[val]['species']
        gene = decoded[val]['description']
        type = decoded[val]["object_type"]
        output[acc] = {'organism':species, 'geneInfo':gene, 'type':type}

    return output

#SECOND TASK
#Function 3
def parse_ID(ID_list):

    for ID in ID_list:
        if not re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', 'r' + ID):
            break
        if ID_list.index(ID) == len(ID_list) - 1:
            return parse_response_Uniprot(get_Uniprot(accessions_list))

    for ID in ID_list:
        if not re.search(r'^((ENS[FPTG]\\d{11}(\\.\\d+)?)|(FB\\w{2}\\d{7})|(Y[A-Z]{2}\\d{3}[a-zA-Z](\\-[A-Z])?)|([A-Z_a-z0-9]+(\\.)?(t)?(\\d+)?([a-z])?))$', 'r' + ID):
            break
        if ID_list.index(ID) == len(ID_list)-1:
            return parse_response_ENSEMBL(get_ENSEMBL(accessions_list))

    return 'Check the correctness of the entered identifiers'

    #Checks the IDs for the regular expression match.
    #Calls the get and parse functions of corresponding databases.
    # Outputs the parsed data.


accessions_list = list(input('Enter the accession IDs separated by spaces:').split())
print(parse_ID(accessions_list))