# packages needed
import gzip

import os
import json
import shutil
import pandas as pd
import config as cfg


## ENTRIES
## add ORCID and date info in the json files, write entries collection and a list of the new pdbs
# todo origin: predicted restano tali, reviewed restano tali, quelli nuovi se prima non c'erano diventano annotated, altrimenti reviewed
from src.uniprot_exon_data import uniprot_exons


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
    new_entries = pd.read_csv(cfg.data['collections'] + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[
        0].to_list()
    print(new_entries)
    entries_dictionary = dict.fromkeys(new_entries, "")
    with open(cfg.data['collections'] + '/output_collections/old_entries.mjson', 'w') as fp:
        with open(cfg.data['collections'] + "/input_collections/entries_20201015.mjson", 'rt') as f:
            for line in f:
                line = json.loads(line)
                line_text = json.dumps(line)
                if line['pdb_id'] in entries_dictionary and 'manually' in line['origin']:
                    pass
                else:
                    fp.write(line_text + '\n')


# concat the remaining pdbs to Entries collection
def create_entries_collection():
    shutil.copyfile(cfg.data['collections'] + '/output_collections/old_entries.mjson',
                    cfg.data['collections'] + '/output_collections/entries_20220516.mjson')

    with open(cfg.data['collections'] + '/output_collections/entries_20220516.mjson', 'a') as fp:
        with open(cfg.data['collections'] + '/output_collections/new_entries.mjson', 'rt') as f:
            for line in f:
                fp.write(line)

def update_entries_collection():
    add_orchid()
    create_old_entries_copy()
    create_entries_collection()
