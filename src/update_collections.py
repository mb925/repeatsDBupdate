# packages needed
import gzip
import os
import json
import shutil

import pandas as pd
import xml.etree.ElementTree as ET
import config as cfg
from src.pdb_mapping import process_pdb

## add ORCID and date info in the json files, write entries collection and a list of the new pdbs
def add_orchid():

    pdbs = {}
    i = 0
    with open(cfg.data['collections'] + '/output_collections/new_entries.mjson', 'w') as fp:
        for subdir, dirs, files in os.walk(cfg.data['data']):
            for file in files:
                print(os.path.join(subdir, file))
                with open(os.path.join(subdir, file), 'rt') as f:
                  for line in f:
                    if (len(os.path.join(subdir, file).split('/')) > 7):
                      line = json.loads(line)
                      line[0]['date'] = os.path.join(subdir, file).split('/')[7]
                      line[0]['annotator'] = os.path.join(subdir, file).split('/')[6]
                      pdbs[i] = line[0]['pdb_id']
                      i = i + 1
                      print(line[0])
                      fp.write(json.dumps(line[0]) + '\n')

    df = pd.DataFrame.from_dict(pdbs, "index")
    df.to_csv(cfg.data['collections'] + '/output_collections/new_pdbs_list.tsv', sep='\t', index=False, header=False)


## create a new Entries collection copy ---but-- if the pdb is in new pdbs list, discard it (it will be overwritten by the new annotated ones)
def create_old_entries_copy():
    new_entries = pd.read_csv(cfg.data['collections'] + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[0].to_list()
    print(new_entries)
    entries_dictionary = dict.fromkeys(new_entries, "")
    with open(cfg.data['collections'] + '/output_collections/old_entries.mjson', 'w') as fp:
        with open(cfg.data['collections'] + "/input_collections/entries_20201015.mjson", 'rt') as f:
            for line in f:
                line = json.loads(line)
                line_text = json.dumps(line)
                if line['pdb_id'] not in entries_dictionary:
                  fp.write(line_text + '\n')


# concat the remaining pdbs to Entries collection
def create_entries_collection():
    shutil.copyfile(cfg.data['collections'] + '/output_collections/old_entries.mjson',
                    cfg.data['collections'] + '/output_collections/entries_20220516.mjson')

    with open(cfg.data['collections'] + '/output_collections/entries_20220516.mjson', 'a') as fp:
        with open(cfg.data['collections'] + '/output_collections/new_entries.mjson', 'rt') as f:
            for line in f:
                fp.write(line)

def mapping(): # generate all xmls
    ## transform PDB indexes of the new collection
    i = 0
    new_entries = pd.read_csv(cfg.data['collections'] + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[0].to_list()
    print(new_entries)
    for pdb in new_entries:
        if os.path.isfile(cfg.data['collections'] + '/output_collections/mappings/' + pdb + '.xml.gz'):
            print("File exists")
        else:
            print("File does not exist, downloading")
            process_pdb(pdb)
            print(pdb)

        i += 1
        print(i)

        df = pd.DataFrame(new_entries)
        # keeping track of where we are if the process breaks
        df.to_csv(cfg.data['collections'] + '/output_collections/to_map_pdbs_list.tsv', sep='\t', index=False, header=False)



## create a new Seqres collection copy ---but-- if the pdb is in new pdbs list, discard it (it will be overwritten by the new annotated ones)
def create_old_seqres_copy():
    new_entries = pd.read_csv(cfg.data['collections'] + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[0].to_list()
    entries_dictionary = dict.fromkeys(new_entries, "")

    with open(cfg.data['collections'] + '/output_collections/old_seqres.mjson', 'w') as fp:
        with open(cfg.data['collections'] + "/input_collections/seqres_20201016.mjson", 'rt') as f:
            for line in f:
                line = json.loads(line)
                line_text = json.dumps(line)
                if line['pdb_id'] not in entries_dictionary:
                    fp.write(line_text + '\n')


def create_new_seqres_collection(): # todo: not found EC ID e dssp

    with open(cfg.data['collections'] + '/output_collections/new_seqres.mjson', 'w') as fp:

        for file in os.listdir(cfg.data['collections'] + '/output_collections/mappings/'):
            print(file)
            if file.endswith(".gz"):

                inp = gzip.open(cfg.data['collections'] + '/output_collections/mappings/' + file, 'r')

                tree = ET.parse(inp)
                root = tree.getroot()

                for x in root[2]:
                    # I could take only pdb chain start end needed
                    # but anyway seqres collection was already freaking big to begin with, it won't hurt if I store even more residues
                    # the bright side is that next time these pdbs will be already downloaded and parsed


                    pdb_id = x.attrib['segId'].split('_')[0]
                    pdb_chain = x.attrib['segId'].split('_')[0]

                    for i in x[0]: # residue
                        entry = {}
                        entry['pdb_id'] = pdb_id
                        entry['pdb_chain'] = pdb_chain
                        entry['seqres_index'] = i.attrib['dbResNum']
                        entry['residue_name_3lett'] = i.attrib['dbResName']
                        for r in i:
                            if r.attrib['dbSource'] == 'PDB':
                                entry['pdb_residue_id'] = i.attrib['dbResNum']
                                entry['observed'] = False
                                if  entry['pdb_residue_id']:
                                    entry['observed'] = True
                                entry['residue_name'] = r.attrib['dbResName']

                            if r.attrib['dbSource'] == 'UniProt':
                                entry['uniprot_index'] = r.attrib['dbResNum']
                                entry['uniprot_residue_name'] = r.attrib['dbResName']
                                entry['uniprot_id'] = r.attrib['dbAccessionId']
                            if r.attrib['dbSource'] == 'CATH':
                                entry['cath'] = r.attrib['dbAccessionId']
                            if r.attrib['dbSource'] == 'InterPro':
                                if 'interpro' not in entry:
                                    entry['interpro'] = []
                                entry['interpro'].append(r.attrib['dbAccessionId'])
                            if r.attrib['dbSource'] == 'SCOP':
                                entry['scop'] = r.attrib['dbAccessionId']
                            if r.attrib['dbSource'] == 'Pfam':
                                entry['pfam'] = r.attrib['dbAccessionId']

                        print(entry)
                        fp.write(json.dumps(entry) + '\n')


# concat the old and new SeqRes collection
def create_seqres_collection():
    shutil.copyfile(cfg.data['collections'] + '/output_collections/old_seqres.mjson', cfg.data['collections'] + '/output_collections/seqres_20220516.mjson')
    with open(cfg.data['collections'] + '/output_collections/seqres_20220516.mjson', 'a') as fp:
        with open(cfg.data['collections'] + '/output_collections/new_seqres.mjson', 'rt') as f:
          for line in f:
            fp.write(line)
# todo: not really sure why this collection is smaller than the previous one...

