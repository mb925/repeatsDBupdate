import os
import gzip
import json
import config as cfg

uniprot_list=[]
counter=0
with open(f'{cfg.entries_path["entries_collection"]}/UniProt_list_from_db_files.txt', 'w+') as uniprot_list_from_db:
    print(f'Parsing all the .db files... {len(os.walk("/mnt/projects/RepeatsDB/repeatsdb3/data/db_files").__next__()[2])}')
    for i in os.listdir('/mnt/projects/RepeatsDB/repeatsdb3/data/db_files'):
        counter +=1
        if counter % 1000 == 0:
            print(f'COUNT: ........... {counter:,} .......... {len(uniprot_list)}')
        folder=i.split('.')[0][1:3]
        pdb_id=i.split('.')[0][:4]
        chain_id=i.split('.')[0][4:]
        with gzip.open(f'/mnt/db/biodbs/seqres/{folder}/{pdb_id}_seqres.mjson.gz', 'rt') as open_pdb:
            for line in open_pdb:
                data = json.loads(line)
                if data['pdb_id'] == pdb_id and data['pdb_chain'] == chain_id:
                    try:
                        if data['uniprot_id'] not in uniprot_list:
                            uniprot_list.append(data['uniprot_id'])
                        else:
                            pass
                    except KeyError:
                        pass
    print(f'Parsed all .db files: {counter}. The number of UniProts associated to these PDBs is: {len(uniprot_list)}')
    print(f'Writing list on file...')
    for i in uniprot_list:
        uniprot_list_from_db.write(f'{i}\n')
