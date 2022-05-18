# packages needed
import gzip
import os
import json
import pandas as pd



## add ORCID and date info in the json files, write entries collection and a list of the new pdbs
from src.pdb_mapping import process_pdb


def add_orchid(rootdir):
    entries = []
    pdbs = {}
    i = 0
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            print(subdir)
            with open(os.path.join(subdir, file), 'rt') as f:
              for line in f:
                if (len(os.path.join(subdir, file).split('/')) > 5):
                  line = json.loads(line)
                  line[0]['date'] = os.path.join(subdir, file).split('/')[4]
                  line[0]['annotator'] = os.path.join(subdir, file).split('/')[3]
                  pdbs[i] = line[0]['pdb_id']
                  entries.append(line[0])
                  i = i + 1
                  print(line[0])

    with open(rootdir + '/output_collections/entries.mjson', 'w') as f:
        json.dump(entries, f)

    df = pd.DataFrame.from_dict(pdbs, "index")
    df.to_csv(rootdir + '/output_collections/new_pdbs_list.tsv', sep='\t', index=False, header=False)


## create a new Entries collection copy ---but-- if the pdb is in new pdbs list, discard it (it will be overwritten by the new annotated ones)
def create_old_entries_copy(rootdir):
    new_entries = pd.read_csv(rootdir + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[0].to_list()
    print(new_entries)
    entries_dictionary = dict.fromkeys(new_entries, "")
    with open(rootdir + '/output_collections/old_entries.mjson', 'w') as fp:
        with open(rootdir + "/input_collections/entries_20201015.mjson", 'rt') as f:
            for line in f:
                line = json.loads(line)
                line_text = json.dumps(line)
                if line['pdb_id'] not in entries_dictionary:
                  fp.write(line_text + '\n')


# concat the remaining pdbs to Entries collection
def create_entries_collection(rootdir):
    new_entries_coll = []
    with open(rootdir + '/output_collections/entries.mjson', 'rt') as f:
      for line in f:
        new_entries_coll = json.loads(line)


    old_entries_coll = []
    with open(rootdir + '/output_collections/old_entries.mjson', 'rt') as f:
      for line in f:
        old_entries_coll = json.loads(line)

    entries_collection = [*old_entries_coll, *new_entries_coll]
    with open(rootdir + '/output_collections/entries_20220516.mjson', 'w') as f:
        json.dump(entries_collection, f)


## create a new Seqres collection copy ---but-- if the pdb is in new pdbs list, discard it (it will be overwritten by the new annotated ones)
def create_old_seqres_copy(rootdir):
    new_entries = pd.read_csv(rootdir + '/output_collections/new_pdbs_list.tsv', sep='\t', header=None)[0].to_list()
    entries_dictionary = dict.fromkeys(new_entries, "")

    with open(rootdir + '/output_collections/old_seqres.mjson', 'w') as fp:
        with open(rootdir + "/input_collections/seqres_20201016.mjson", 'rt') as f:
            for line in f:
                line = json.loads(line)
                line_text = json.dumps(line)
                if line['pdb_id'] not in entries_dictionary:
                    fp.write(line_text + '\n')
# todo check why the last id is incomplete in output file


def mapping(rootdir): # todo: finish generating all xmls and see if I can make the process faster (I don't need tsv for example)
    ## transform PDB indexes of the new collection
    new_entries = pd.read_csv(rootdir + '/output_collections/to_map_pdbs_list3.tsv', sep='\t', header=None)[0].to_list()
    for pdb in new_entries:
        process_pdb(pdb, rootdir) # todo: execute if you do not have the file already
        new_entries.remove(pdb)
        print(len(new_entries))
        df = pd.DataFrame(new_entries)
        df.to_csv(rootdir + '/output_collections/to_map_pdbs_list.tsv', sep='\t', index=False, header=False) # keep track of where we are because the process breaks (idk why, connection timeout)


def create_seqres_collection(rootdir):
    for file in os.listdir(rootdir + 'output_collections/mappings/'):
        if file.endswith(".gz"):
            print(file)
            a_file = gzip.open(rootdir + 'output_collections/mappings/' + file, "rb")
            xml_contents = a_file.read()
            print(xml_contents)
