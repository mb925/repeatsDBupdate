import requests
import json
from ftplib import FTP
import config as cfg

aa_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
aa_d_r = {aa_d[key]: key for key in aa_d}


def download(url):
    r = requests.get(url, headers={"Accept": "application/json"})
    if not r.ok:
        # r.raise_for_status()
        return None
    else:
        return r.text


# uniprot from PDB and PDB from uniprot
# https://www.ebi.ac.uk/pdbe/api/mappings/P29373
def map_accession(acc):
    service = 'https://www.ebi.ac.uk/pdbe/api/mappings/'
    data = download(service + acc)
    if data:
        return data
    else:
        print('Accession', acc, 'error in retrieving exons from Protein API')
        return None


# uniprot sequence from Protein API
# https://www.ebi.ac.uk/proteins/api/proteins/Q04724
def download_sequence(acc):
    service = 'https://www.ebi.ac.uk/proteins/api/proteins/'
    data = download(service + acc)
    if data:
        return json.loads(data)
    else:
        print('Accession', acc, 'error in retrieving sequence from Protein API')


def residue_list(pdb):
    service = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/'
    data = download(service + pdb)
    if data:
        entities = json.loads(data)[pdb]["molecules"]
        res = {
            e["entity_id"]:
                {c["struct_asym_id"]:
                     [r for r in c["residues"]]
                 for c in e["chains"]}
            for e in entities}
        return res
    else:
        print('Residue mapping', pdb, 'error')
        return None


def process_pdb(pdb, folder=cfg.data['collections'] + '/output_collections/mappings/'): # download_sifts
    ftp_folder = pdb[1:3]
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd('pub/databases/msd/sifts/split_xml/' + ftp_folder)
    flname = folder + '{}.xml.gz'.format(pdb)
    print('download ' + flname)
    ftp.retrbinary('RETR {}.xml.gz'.format(pdb), open(flname, 'wb').write)
