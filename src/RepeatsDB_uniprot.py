import gzip
import argparse, sys, logging
# from new_function import *
import itertools
import json
import os
import pickle
import config as cfg
# from uniprot_exon_data import *
from collections import defaultdict

from src.main_function_uniprot import add_repeatsdb_annotation, obtain_uniprot_sequence, \
    obtain_segment_interpro_pfam, obtain_missing_residues, obtain_vector_of_repeated_elements, unify_all_units, \
    make_dataframe_to_align, consensus, missing_residues_mapping, obtain_period_of_consensus_all, obtain_json, \
    obtain_unip_ind_from_PDB


def retrieve_uniprot_list_from_seqres():
    uniprots = []
    # unip_id = []
    # unip = []
    uniprot_list = {}
    with open(cfg.data_path['data'] + 'good_code/output.json', 'r') as list_json:
        for _ in range(5):
            next(list_json)
        data = json.load(list_json)
        for ln in data:
            uniprots.append(ln['_id'])

    for j in uniprots:
        unip_id = []
        for line in data:
            if line['_id'] == j:
                for i in line['pdb_id']:
                    unip_id.append(f'{i["pdb_id"]}{i["pdb_chain"]}')
        uniprot_list[f'{j}'] = f'{unip_id}'
    # print(uniprot_list)
    return uniprot_list

def subset_uniprot_list(list_all_uniprot):
    list_db_files = []
    for i in os.listdir(cfg.data_path['data'] + 'pre_submission/db_files_15_10_20'):
        list_db_files.append(i.split('.')[0])
    unip_pdb = defaultdict(list)
    for k in list_all_uniprot:
        if set(list_all_uniprot[k]) & set(list_db_files):
            unip_pdb[k].extend(list_all_uniprot[k])
    # with open(cfg.data_path['data'] + 'good_code/list_RDB3_uniprot_and_pdbs.json', 'w') as f:
    #     json.dump(unip_pdb, f)
    return unip_pdb


def get_pdb_up_dict(pdb_to_uniprot, ups_from_rdb):
    res_dict = defaultdict(list)
    with gzip.open(pdb_to_uniprot, 'rt') as f:
        for row in f.readlines()[2:]:
            fields = row.strip().split('\t')
            if fields[2] in ups_from_rdb:
                res_dict[fields[2]].append([fields[0], fields[1], fields[7], fields[8]])
    return res_dict


def get_up_interpro_dict(uniprot_to_interpro, ups_from_rdb):
    pickle_file = f'{uniprot_to_interpro}.pickle' # todo: right one has to be generated
    if os.path.exists(pickle_file):
        logging.debug(f'  loading mapping from pickle file {pickle_file} ...')
        with open(pickle_file, 'rb') as handle:
            res_dict = pickle.load(handle)
        return res_dict
    else:
        res_dict = defaultdict(list)
        count = 0
        line_count = 0
        with gzip.open(uniprot_to_interpro, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                line_count += 1
                if line_count % 10000000 == 0:
                    logging.debug(f'  ... {line_count:,} ...')
                if fields[0] in ups_from_rdb:
                    count += 1
                    if fields[3].startswith('PF'):
                        res_dict[fields[0]].append([fields[1], fields[3], fields[4], fields[5]])
                    else:
                        res_dict[fields[0]].append([fields[1], '[]', fields[4], fields[5]])
        dict_out = dict(res_dict)
        logging.debug(f'  saving mapping to pickle file {pickle_file} ...')
        with open(pickle_file, 'wb') as handle:
            pickle.dump(dict_out, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return dict(res_dict)


def main():
    list_all_uniprot = retrieve_uniprot_list_from_seqres()
    unip_pdb = {"A0A011": ["3vkcB", "3vkcA", "3vkaB", "3vkdB", "3vkbA", "3vkaA", "3vk5B", "3vkdA", "3vkbB", "3vk5A"]}
    ups_from_rdb = unip_pdb.keys()

    # print(unip_pdb)
    # exit()
    #
    # with open('./UniProt_with_entries_RDB3.json', 'w') as unip_and_pdb:
    #     json.dump(unip_pdb, unip_and_pdb)
    #
    # ups_from_rdb = unip_pdb.keys()

    pdb_uniprot_file = open(cfg.data_path['data'] + 'good_code/PDB_for_seqres_filtering.txt', 'w+')

    # initargs = '--pdb-to-uniprot ../pdb_chain_uniprot.tsv.gz --up-rdb ./UniProt_with_entries_RDB3.json --uniprot-to-interpro ../p2ipr.dat.gz'.split()
    initargs = '--pdb-to-uniprot ../pdb_chain_uniprot.tsv.gz --uniprot-to-interpro ./../data/good_code/protein2ipr_small2.dat.gz'.split()

    # with open('./UniProt_with_entries_RDB3.json') as data:
    #     unip_pdb = json.load(data)
    #     ups_from_rdb = unip_pdb.keys()

    # subset_uniprot_list(list_all_uniprot)

    args = parse_arguments(initargs)

    setup_logging(args)

    logging.info(f'Loaded {len(unip_pdb)} UniProts with associated PDBs')

    logging.info(f'Loading UniProt to Interpro mappings from {args.uniprot_to_interpro} ... ')
    up_pfam_interpro_dict = get_up_interpro_dict(args.uniprot_to_interpro, ups_from_rdb)
    logging.info(f'Mapped {len(up_pfam_interpro_dict)} UniProts to Interpro')

    # p = re.compile(' [A-Z][A-Z]=.*')
    # all_interest_pdb=[]
    count = 0
    with open(cfg.data_path['data'] + 'good_code/uniprot_collection_2020_10_15_prova.mjson', 'w') as uniprot_collection:
        # for unip in pdb_up_dict.keys():
        for unip in unip_pdb.keys():
            # unip='P62871'
            # unip = 'P31224' #none our PDBs associated to this uniprot
            # unip = 'P46675' # 3 regions
            # unip = 'P18754' #missing residue
            # unip = 'P0A749' # pdb with 2 regions
            # unip = 'P20436' # missing residues
            # unip = 'P03528' # stupid regions
            ###################
            # unip = 'Q11NN9' #not in blast database... (deprecated)
            # unip = 'Q11TI6' #not in blast database... (deprecated)
            ###################
            print(unip)
            count += 1
            if count % 100 == 0:
                print(f'COUNT: ........... {count:,} ............. {unip}')
            if unip == 'P69905':
                print('UniProt P69905 associated to Hemoglobin is still in the database...')
                print('Removing PDB chains associated to this UniProt...')
                vector, class_and_sub, infos_pdb, unit_pdb, units_uniprot_index, useful_vector, pdbs_of_the_uniprot, complete_pdb = obtain_unip_ind_from_PDB(unip_pdb[unip], unip)
                pdbs_removed = 0
                for i in complete_pdb:
                    for j in os.listdir('/mnt/projects/RepeatsDB/repeatsdb3/data/new_source_db_files'):
                        if i.rstrip() == j.split('.')[0]:
                            os.system(f'rm /mnt/projects/RepeatsDB/repeatsdb3/data/new_source_db_files/{i.rstrip()}.db')
                            # print(i)
                            pdbs_removed += 1
                print(f'Removed {pdbs_removed}')
                pass
            else:
                #### 1) see if the uniprot is curated or not
                unit_unip, ins_unip, reg_unip, curator = add_repeatsdb_annotation(unip)

                #### 2) obtain uniprot TITLE and SEQUENCE
                # unip_info = obtain_uniprot_sequence(unip) # todo: fix
                unip_info = 'MNASPQLDHHTELHAAPPLWRPGRVLARLREHQPGPVHIIDPFKVPVTEAVEKAAELTRLGFAAVLLASTDYESFESHMEPYVAAVKAATPLPVVLHFPPRPGAGFPVVRGADALLLPALLGSGDDYFVWKSFLETLAAFPGRIPREEWPELLLTVALTFGEDPRTGDLLGTVPVSTASTEEIDRYLHVARAFGFHMVYLYSRNEHVPPEVVRHFRKGLGPDQVLFVSGNVRSGRQVTEYLDSGADYVGFAGALEQPDWRSALAEIAGRRPAAPARPGSGR\nTITLE'

                if not unip_info:
                    logging.warning(f'No information finding for UniProt {unip.rstrip()} ...')
                    # pass
                else:
                    unip_seq = unip_info.split('\n')[1] # todo: fix
                    unip_title = unip_info.split('\n')[0] # todo: fix
                    uniprot_vector = list(range(1, len(unip_seq)))  #a vector of the same lenght of uniprot sequence

                    #### 3) obtain PFAM and INTERPRO informations
                    interpro_pfam_coord = obtain_segment_interpro_pfam(unip, up_pfam_interpro_dict)

                    #### 4) find the pdb units in uniprot index (units_uniprot_index)
                    vector, class_and_sub, infos_pdb, unit_pdb, units_uniprot_index, useful_vector, pdbs_of_the_uniprot, complete_pdb = obtain_unip_ind_from_PDB(unip_pdb[unip], unip)
                    # exit()
                    for pdbs in pdbs_of_the_uniprot:
                        pdb_uniprot_file.write(f'{pdbs.rstrip()}\n')

                    if all([not elem for elem in units_uniprot_index]) == True:
                        logging.warning(f'This UniProt {unip.rstrip()} is empty ...')
                        pass
                    else:
                        #### 5) find missing residues in the units of each PDB chains
                        missing_residues = obtain_missing_residues(units_uniprot_index, useful_vector)
                        list_missing_residues = list(itertools.chain.from_iterable(missing_residues)) # questa lista è una lista di tutti i missing residues

                        #### 6)
                        ########questo è il vettore delle unità complete, con i missing residues presenti; è praticamente un vettore start-end delle unità
                        new_final_vector_miss = obtain_vector_of_repeated_elements(useful_vector)

                        ########questo è il vettore senza missing residues presenti, con i buchi
                        new_final_vector = unify_all_units(units_uniprot_index)

                        #### 7) make the dataframe for the consensus
                        R_vector = make_dataframe_to_align(vector, uniprot_vector, new_final_vector_miss)

                        #### 8) obtain consensus one (O), majority (M) and all (A)
                        O_set, M_set, A_set, vector_consensus = consensus(R_vector, useful_vector)

                        #### 9) see where the missing resides map on the consensus
                        true_missing_residues = missing_residues_mapping(vector_consensus, list_missing_residues, R_vector)

                        if not O_set:
                            logging.warning(f'No consensus finding for UniProt {unip.rstrip()} ...')
                            pass
                        else:
                            period_one = obtain_period_of_consensus_all(O_set, units_uniprot_index)
                            uniprot_json = obtain_json(unip, unip_seq, unip_title, interpro_pfam_coord, O_set, M_set, A_set,
                                                       class_and_sub, infos_pdb, unit_unip, ins_unip, reg_unip, curator,
                                                       period_one, true_missing_residues, useful_vector)
                            # print(uniprot_json)
                            # exit()
                            logging.info(f'Saving the UniProt {unip.rstrip()} ...')
                            json.dump(uniprot_json, uniprot_collection)
                            ############ to save the collection directly in mongo (remember to change the localhost) ################
                            # logging.info('Importing the collection to mongodb...')
                            # to_mongo(uniprot_json)
                            # exit()
                            # uniprot_exons(unip.rstrip(), folder='./', dest=None)
            # download_exon(unip.rstrip(), folder='', dest=None, save=True)
            # exit()

    pdb_uniprot_file.close()
    # uniprot_list.close()
    return 0


def parse_arguments(initargs=None):
    if initargs is None:
        initargs = []
    parser = argparse.ArgumentParser(description="Code to create the collection 'RepeatsDB-UniProt'.")
    parser.add_argument('--pdb-to-uniprot', required=True, help="UniProt to Interpro mapping file")
    parser.add_argument('--uniprot-to-interpro', required=True, help="PDB to UniProt mapping file")
    # parser.add_argument('--up-rdb', required=True, help="File with UniProt IDs used by RepeatsDB")

    if not initargs or len(sys.argv) > 1:
        args = parser.parse_args()
    else:
        args = parser.parse_args(initargs)
    return args


def setup_logging(args):
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s | %(levelname)-5.5s | %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S', stream=sys.stdout)
    logging.info(f'{os.path.basename(__file__)} started')
    logging.debug(f'Arguments: {vars(args)}')
    return


if __name__ == '__main__':
    sys.exit(main())
