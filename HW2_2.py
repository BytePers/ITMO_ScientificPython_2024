import requests, re, sys, json, subprocess
from Bio import SeqIO
from tabulate import tabulate

def get_Uniprot(ID: list):
  return requests.get("https://rest.uniprot.org/uniprotkb/accessions", params={'accessions': ID})

#Function 2 “parse_response_<db_name>” that accepts response and outputs parsed data.
def parse_response_Uniprot(resp = None):
    resp = resp.json()
    resp = resp["results"]
    output = ''
    for val in resp:
        acc = val['primaryAccession']
        species = val['organism']['scientificName']
        gene = val['genes']
        seq = val['sequence']
        output = f'accession: {acc}, ' + f'organism:{species}, ' + 'type:protein'
    return output

#Function 1 “get_<db_name>” that accepts list of IDs and outputs API request.
def get_ENSEMBL(ID):
    json_ids = json.dumps({"ids": [ID]})
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
        type = decoded[val]["object_type"]
        output = f'accession: {acc}, ' + f'organism:{species}, ' + f'type:{type}'

    return output

#SECOND TASK
#Function 3
def parse_ID(ID, database):

    match database:
        case 'uniprot':
            return parse_response_Uniprot(get_Uniprot(ID))
        case 'ensembl':
            return parse_response_ENSEMBL(get_ENSEMBL(ID))

    return 'Check the correctness of the entered identifiers'

    #Checks the IDs for the regular expression match.
    #Calls the get and parse functions of corresponding databases.
    # Outputs the parsed data.


def SeqKit(filename):
    seqkit = subprocess.run(("/home/ilya/anaconda3/envs/seqkit/bin/seqkit", "stats", filename, "-a"),
                                capture_output=True, # When set to True, will capture the standard output and standard error
                                text=True) # When set to True, will return the stdout and stderr as string, otherwise as bytes

    #Parse an output from the seqkit. If collection raised an error – capture it and return as a final output.
    if seqkit.stdout:
        seqkit_out = seqkit.stdout.strip().split('\n')
        # split names and values
        prop_names = seqkit_out[0].split()[1:]
        prop_vals = seqkit_out[1].split()[1:]
        # using zip
        seq_result = dict(zip(prop_names, prop_vals))
        return seq_result

    elif seqkit.stderr:
        return seqkit.stderr


def main(fileList):
    #clear output file before running the script
    with open('./parse_output.txt', 'w') as f:
        # Write new data to the file
        f.write('')
    for filename in fileList:
        seq_result = SeqKit(filename)
        parse_response = []
        ext = filename.split('.')[-1]
        sequences = SeqIO.parse(filename, ext)
        # if seq_result == seqkit.stderr:
        #     print('1', seq_result)
        #
        big_sep = '■■■■■■■■■■■'*16
        sep = '-----------'*10
        if type(seq_result) == dict:
            seq_table = list(seq_result.items())
            table = tabulate(seq_table, headers=("Parameter", "Value"), tablefmt='fancy_grid', colalign=("left", "right"))
        try:
            if seq_result['type'] and seq_result['type'] == 'Protein':
                parse_response.append('**Database: Uniprot**' + '\n' + sep)
                for seq in sequences:
                    parse_response.append('\n' + seq.description)
                    match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', 'r' + seq.id)
                    parse_response.append(str(parse_ID(match[0], 'uniprot')) + '\n' + sep)
            elif seq_result['type']:
                parse_response.append('**Database: Ensembl**' + '\n' + sep)
                for seq in sequences:
                    parse_response.append('\n' + seq.description)
                    match = re.search(r'(ENS[a-zA-Z]{0,4}(E|FM|G|GT|P|R|T)\d{11}|MGP[a-zA-Z\_]*\d{7})', 'r' + seq.id)
                    parse_response.append(str(parse_ID(match[0], 'ensembl')) + '\n' + sep)
            final_response = big_sep + '\n'+f"****filename: {'.'.join(filename.split('.')[-2:])};" + f" data_type: {seq_result['type']}****\n" + table + '\n' + '\n'.join(parse_response) + '\n'
        except TypeError:
            parse_response.append(str())
            final_response = big_sep + '\n'+f"****filename: {'.'.join(filename.split('.')[-2:])}****" + '\n' + f"{seq_result}"
        with open('./parse_output.txt', 'a') as output:
            output.write(final_response)

fileList = ['./hw_file1.fasta', './hw_file2.fasta', './hw_file3.fasta']
main(fileList)





