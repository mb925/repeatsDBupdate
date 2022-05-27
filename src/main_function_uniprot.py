import gzip, json, os, sys
import logging
import pandas as pd
import numpy as np
import collections
import pymongo
from operator import itemgetter
import itertools
from itertools import groupby, count
import subprocess
import warnings
from datetime import date
import operator
from Bio import BiopythonWarning, PDB
from Bio.PDB import *
from Bio import BiopythonWarning
from urllib.request import urlretrieve
from uniprot_exon_data import *
import config as cfg

warnings.simplefilter('ignore', BiopythonWarning)

def add_repeatsdb_annotation(unip):
    unit_unip, reg_unip, ins_unip, curator = [], [], [], []
    if os.path.exists(cfg.data_path['data'] + '/curated_uniprot/' + unip + '.txt'):
        with open(cfg.data_path['data'] + '/curated_uniprot/' + unip + '.txt') as rdb_ann:
            for line in rdb_ann:
                if line.startswith('CURATOR'):
                    curator.append(line.split('\t')[1].rstrip())
                if line.startswith('UNIT'):
                    unit = (int(line.split('\t')[1]), int(line.split('\t')[2].rstrip()))
                    unit_unip.append(unit)
                if line.startswith('INS'):
                    ins = (int(line.split('\t')[1]), int(line.split('\t')[2].rstrip()))
                    ins_unip.append(ins)
                if line.startswith('REG'):
                    reg = (int(line.split('\t')[1]), int(line.split('\t')[2]), '{}.{}'.format(line.split("\t")[3], line.split("\t")[4].rstrip()))
                    reg_unip.append(reg)
    return unit_unip, ins_unip, reg_unip, curator


def obtain_uniprot_sequence(unip):
    cmd = f'blastdbcmd -db /mnt/db/blastdb/uniprot.fasta -entry {unip} -outfmt'.split()
    # cmd = f'blastdbcmd -db /local/db/blastdb/uniprot.fasta -entry {unip} -outfmt'.split()
    cmd.append('%t\n%s')
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as proc:
        stdout, stderr = proc.communicate()
    if stderr:
        logging.error(f'{stderr.rstrip()}. Check again...')
        stdout = []
    # remove OS and other infos from title with 'p'
    # std = p.sub('', stdout, count=1)
    # print(stdout)
    return stdout


def obtain_segment_interpro_pfam(unip, up_interpro_dict):
    for key in up_interpro_dict.keys():
        if key == unip:
            return up_interpro_dict[key]


def obtain_unip_ind_from_PDB(dict_uniprot, unip):
    print(dict_uniprot)
    vector, infos_pdb, class_sub, units_uniprot, useful_vector = [], [], [], [], []
    all_pdb_ass_to_uniprot=[]
    complete_pdb=[]
    for i in range(len(dict_uniprot)):

        if dict_uniprot[i][0] not in all_pdb_ass_to_uniprot:
            all_pdb_ass_to_uniprot.append(dict_uniprot[i][0])
        else:
            pass
        PDB_name = dict_uniprot[i][:4]
        chain_ID = dict_uniprot[i][4:]

        complete_pdb.append('{}{}'.format(PDB_name, chain_ID))
        resseqs = obtain_resseqs(PDB_name, chain_ID, unip)
        reg_pdb, class_and_sub, unit_pdb, pdb_not_rdb, origin = obtain_reg_from_db(PDB_name, chain_ID)
        uniprot_index, residue_name, start_end, units_uniprot_index, reg_unip = from_pdb_reg_to_unip(PDB_name, chain_ID, reg_pdb, resseqs, unit_pdb, unip)

        if start_end:
            infos_pdb.append('{}{}\t{}\t{}\t{}\t{}\t{}'.format(PDB_name, chain_ID, start_end[0], start_end[len(start_end) - 1], residue_name, origin, list(set(class_and_sub))))
        else:
            if not uniprot_index:
                pass
            else:
                infos_pdb.append('{}{}\t{}\t{}\t{}\t{}\t{}'.format(PDB_name, chain_ID, uniprot_index[0], uniprot_index[len(uniprot_index)-1], residue_name, origin, list(set(class_and_sub))))
        if not reg_pdb:
            pass
        else:
            vector.append(uniprot_index)
        for item in class_and_sub:
            if item not in class_sub:
                class_sub.append(item)
        units_uniprot.append(units_uniprot_index)
        if not units_uniprot_index:
            pass
        else:
            if not class_and_sub:
                pass
            else:
                new_uni_ind=[]
                for uui in units_uniprot_index:
                    if not uui:
                        pass
                    else:
                        new_uni_ind.append(uui)
                if (not new_uni_ind) and (all([not elem for elem in units_uniprot_index ])==True):
                    pass
                elif not new_uni_ind:
                    useful_vector.append('{}\t{}\t{}\t{}\t{}\t{}'.format(PDB_name, chain_ID, units_uniprot_index[0][0], units_uniprot_index[len(units_uniprot_index) - 1][len(units_uniprot_index[len(units_uniprot_index) - 1]) - 1], list(set(class_and_sub)), reg_unip))
                else:
                    sorted_new_uni_ind=sorted(new_uni_ind, key=itemgetter(0))
                    useful_vector.append('{}\t{}\t{}\t{}\t{}\t{}'.format(PDB_name, chain_ID, sorted_new_uni_ind[0][0], sorted_new_uni_ind[len(sorted_new_uni_ind)-1][len(sorted_new_uni_ind[len(sorted_new_uni_ind)-1])-1], list(set(class_and_sub)), reg_unip))

    return vector, class_sub, infos_pdb, unit_pdb, units_uniprot, useful_vector, list(set(all_pdb_ass_to_uniprot)), complete_pdb


def obtain_resseqs(PDB_name, chain_ID, unip):
    residue=[]
    folder=PDB_name[1:3]
    try:
        with gzip.open(f'/mnt/db/biodbs/seqres/{folder}/{PDB_name}_seqres.mjson.gz', 'rt') as open_pdb:
            for line in open_pdb:
                data = json.loads(line)
                try:
                    if data['pdb_id'] == PDB_name and data['pdb_chain'] == chain_ID and data['uniprot_id'] == unip.rstrip():
                        residue.append(data['pdb_residue_id'])
                except KeyError:
                    pass
        while 'null' in residue: residue.remove('null')

    except FileNotFoundError:
        if f'{PDB_name}.pdb' in os.listdir('/home/sara/Documents/new_SRUL/data/SRUL_REPEATSDB/DATABASE_PDB'):
            pass
        else:
            url = urlretrieve(f'http://files.rcsb.org/download/{PDB_name}.pdb', f'/home/sara/Documents/new_SRUL/data/SRUL_REPEATSDB/DATABASE_PDB/{PDB_name}.pdb')
        structure = PDBParser().get_structure(PDB_name, f'/home/sara/Documents/new_SRUL/data/SRUL_REPEATSDB/DATABASE_PDB/{PDB_name}.pdb')
        model = structure[0]
        chain = model[chain_ID]
        residue = [residue.id[1] for residue in chain.get_residues()]

    resseqs_1 = list(map(str, residue))
    return resseqs_1


def obtain_reg_from_db(PDB_name, chain_ID):
    class_and_sub, reg_pdb, unit_pdb = [], [], []
    pdb_not_rdb = []
    # if PDB_name + chain_ID + '.db' in os.listdir('/mnt/projects/RepeatsDB/repeatsdb3/data/db_files'):
    #     file_db = open(f'/mnt/projects/RepeatsDB/repeatsdb3/data/db_files/{PDB_name + chain_ID}.db')
    if PDB_name + chain_ID + '.db' in os.listdir(cfg.data_path['data'] + '/pre_submission/db_files_15_10_20'):
        file_db = open(cfg.data_path['data'] + '/pre_submission/db_files_15_10_20/' +PDB_name + chain_ID + '.db')
    else:
        file_db = []
        pdb_not_rdb.append(f'{PDB_name}\t{chain_ID}')
        origin = 'Not analyzed'

    for line in file_db:
        if line.startswith('SOURCE'):
            origin = line.split('\t')[1].rstrip()
        if line.startswith('UNIT'):
            unit = (line.split('\t')[1].split(' ')[0], line.split('\t')[1].split(' ')[1].rstrip())
            unit_pdb.append(unit)
        if line.startswith('REG'):
            class_and_sub.append('{}.{}'.format(line.split('\t')[1].split(' ')[2], line.split('\t')[1].split(' ')[3].rstrip()))
            region = (line.split('\t')[1].split(' ')[0], line.split('\t')[1].split(' ')[1])
            reg_pdb.append(region)
    return reg_pdb, class_and_sub, unit_pdb, pdb_not_rdb, origin


def from_pdb_reg_to_unip(PDB_name, chain_ID, reg_pdb, resseqs, unit_pdb, unip):
    reg_in_unip_ind=[]
    uniprot_index, residue_name, pdb_residue_id, info_pdb = [], [], [], []
    if not reg_pdb and not unit_pdb:
        ress=[]
    else:
        ress = list(map(str, resseqs))
    folder = PDB_name[1:3]
    try:
        with gzip.open(f'/mnt/db/biodbs/seqres/{folder}/{PDB_name}_seqres.mjson.gz', 'rt') as open_pdb:
            for line in open_pdb:
                data = json.loads(line)
                try:
                    if data['pdb_id'] == PDB_name and data['pdb_chain'] == chain_ID and data['uniprot_id'] == unip.rstrip():
                        pdb_residue_id.append(data['uniprot_index'])
                        residue_name.append(data['residue_name'])
                except KeyError:
                    pass
                if not reg_pdb:
                    pass
                elif not resseqs:
                    pass
                else:
                    for i in range(len(reg_pdb)):
                        if (int(reg_pdb[i][0]) > int(ress[len(ress)-1])) and (int(reg_pdb[i][1]) > int(ress[len(ress)-1])):
                            result=[]
                        elif (int(reg_pdb[i][0]) < int(ress[0])) and (int(reg_pdb[i][1]) < int(ress[0])):
                            result=[]
                        elif (str(reg_pdb[i][0]) in ress) and (str(reg_pdb[i][1]) not in ress):
                            result = ress[ress.index(str(reg_pdb[i][0])):]
                        elif (str(reg_pdb[i][1]) in ress) and (str(reg_pdb[i][0]) not in ress):
                            result = ress[: ress.index(str(reg_pdb[i][1])) + 1]
                        elif (int(reg_pdb[i][0]) < int(ress[0])) and (int(reg_pdb[i][1]) > int(ress[len(ress)-1])):
                            result=ress
                        else:
                            try:
                                if data['pdb_id'] == PDB_name and data['pdb_chain'] == chain_ID:
                                    try:
                                        if data['uniprot_id'] == unip.rstrip():
                                            if (data['pdb_residue_id'] == reg_pdb[i][0]) or (data['pdb_residue_id'] == reg_pdb[i][1]):
                                                reg_in_unip_ind.append(int(data['uniprot_index']))
                                    except KeyError:
                                        pass
                                if not ress:
                                    pass
                                else:
                                    result = ress[ress.index(str(reg_pdb[i][0])): ress.index(str(reg_pdb[i][1])) + 1]
                            except ValueError:
                                int_list = [int(x) for x in ress]
                                missing_values = []
                                for miss in range(min(int_list), max(int_list)):
                                    if not miss in int_list:
                                        missing_values.append(miss)


                                missing_res = []
                                for k, g in groupby(enumerate(missing_values), lambda i_x: i_x[0] - i_x[1]):
                                    mist = (list(map(itemgetter(1), g)))
                                    missing_res.append(mist)
                                try:
                                    if (int(reg_pdb[i][0]) in missing_values) and (
                                            int(reg_pdb[i][1]) in missing_values):
                                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(
                                            str(missing_res[0][len(missing_res[0]) - 1] + 1))]
                                    elif int(reg_pdb[i][0]) < int(ress[0]):
                                        result = ress[ress.index(str(ress[0])): ress.index(str(missing_values[0] - 1))]
                                    elif int(reg_pdb[i][1]) > int(ress[len(ress) - 1]):
                                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(
                                            str(ress[len(ress) - 1]))]
                                    else:
                                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(
                                            str(reg_pdb[i][1]))]
                                except ValueError:

                                    result = ress[ress.index(str(missing_values[0]-1)): ress.index(str(reg_pdb[i][1]))]
                        if not result:
                            pass
                        else:
                            for j in result:
                                if data['pdb_residue_id'] == j:
                                    if data['pdb_id'] == PDB_name and data['pdb_chain'] == chain_ID:
                                        try:
                                            if data['uniprot_id'] == unip.rstrip():
                                                if j == result[0] or j == result[len(result)-1]:
                                                    reg_in_unip_ind.append(data['uniprot_index'])
                                                uniprot_index.append(data['uniprot_index'])
                                        except KeyError:
                                            pass
    except FileNotFoundError:
        pass

    reg_unip=list(set(reg_in_unip_ind))
    sorted_reg_unip=sorted(reg_unip, key=int)
    it = iter(sorted_reg_unip)
    region_unip=list(zip(it, it))
    if not reg_unip:
        units_uniprot_index = []
        pass
    else:
        units_uniprot_index = obtain_unit_uniprot(folder, PDB_name, ress, unit_pdb, chain_ID, unip)

    residue_name = "".join(residue_name)
    return uniprot_index, residue_name, pdb_residue_id, units_uniprot_index, region_unip


def obtain_unit_uniprot(folder, PDB_name, ress, unit_pdb, chain_ID, unip):
    uniprot_index = []
    for i in range(len(unit_pdb)):
        unit_upr = []
        try:
            if (int(unit_pdb[i][0]) > int(ress[len(ress)-1])) and (int(unit_pdb[i][1]) > int(ress[len(ress)-1])):
                result = []
            elif (int(unit_pdb[i][0]) < int(ress[0])) and (int(unit_pdb[i][1]) < int(ress[0])):
                result = []
            elif (int(unit_pdb[i][0]) < int(ress[0])) and (int(unit_pdb[i][1]) > int(ress[len(ress)-1])):
                result = ress
            elif (str(unit_pdb[i][0]) in ress) and (str(unit_pdb[i][1]) not in ress):
                result = ress[ress.index(str(unit_pdb[i][0])):]
            elif (str(unit_pdb[i][1]) in ress) and (str(unit_pdb[i][0]) not in ress):
                result = ress[: ress.index(str(unit_pdb[i][1])) + 1]
            else:
                result = ress[ress.index(str(unit_pdb[i][0])): ress.index(str(unit_pdb[i][1])) + 1]
        except ValueError:
            try:
                result = ress[ress.index(str(unit_pdb[i][0])): ress.index(str(unit_pdb[i][1]) + 'A') + 1]
            except ValueError:
                int_list = [int(x) for x in ress]
                missing_values = []
                for miss in range(min(int_list), max(int_list)):
                    if not miss in int_list:
                        missing_values.append(miss)
                missing_res=[]
                for k, g in groupby(enumerate(missing_values), lambda i_x: i_x[0] - i_x[1]):
                    mist=(list(map(itemgetter(1), g)))
                    missing_res.append(mist)
                try:
                    if (int(unit_pdb[i][0]) in missing_values) and (int(unit_pdb[i][1]) in missing_values):
                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(str(missing_res[0][len(missing_res[0])-1] + 1))]
                    elif int(unit_pdb[i][0]) < int(ress[0]):
                        result = ress[ress.index(str(ress[0])): ress.index(str(missing_values[0] - 1))]
                    elif int(unit_pdb[i][1]) > int(ress[len(ress)-1]):
                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(str(ress[len(ress)-1]))]
                    else:
                        result = ress[ress.index(str(missing_values[0] - 1)): ress.index(str(unit_pdb[i][1]))]
                except ValueError:
                    result = ress[ress.index(str(unit_pdb[i][0])): ress.index(str(missing_values[0] - 1))]
        if not result:
            pass
        else:
            with gzip.open(f'/mnt/db/biodbs/seqres/{folder}/{PDB_name}_seqres.mjson.gz', 'rt') as open_pdb:
                for line in open_pdb:
                    data = json.loads(line)
                    try:
                        for j in result:
                            if data['pdb_id'] == PDB_name and data['pdb_chain'] == chain_ID and data['uniprot_id'] == unip.rstrip():
                                if data['pdb_residue_id'] == j:
                                    unit_upr.append(data['uniprot_index'])
                    except KeyError:
                        pass
            uniprot_index.append(unit_upr)
    return uniprot_index


def obtain_missing_residues(final_vector, useful_vector):
    new=[]
    new_final_vector, missing_res=[], []
    for i in range(len(final_vector)):
        if not final_vector[i]:
            pass
        else:
            for j in final_vector[i]:
                new_v=list(set(j))
                if not j:
                    pass
                else:
                    int_list = [int(x) for x in new_v]
                    missing_values = []
                    for miss in range(min(int_list), max(int_list)):
                        if not miss in int_list:
                            missing_values.append(miss)
                    if not missing_values:
                        new_final_vector.append(j)
                    else:
                        prova =list(range(j[0], j[len(j)-1]+1))
                        new_final_vector.append(prova)
                        missing_res.append(missing_values)
    missing_res.sort()
    m_r_s = list(missing_res for missing_res, _ in itertools.groupby(missing_res))
    return m_r_s


def make_dataframe_to_align(final_v, uniprot_vector, new_final_vector):

    if len(new_final_vector) == 1:
        final_vector = new_final_vector
    else:
        final_vector = list(new_final_vector)
    new_vector = []

    for i in range(len(final_vector)):
        if not final_vector[i]:
            pass
        else:
            new_list = [x if x in final_vector[i] else '-' for x in uniprot_vector]
            for ind, t in enumerate(new_list):
                if t == '-':
                    pass
                else:
                    new_list[ind] = 'R'
            sequence = "".join(new_list)
            new_vector.append(sequence)
    return new_vector


def obtain_vector_of_repeated_elements(useful_vector):
    new=[]
    for i in useful_vector:
        prova=[]
        tuple_reg = []
        for j in i.split('\t')[5].split(','):
            if i.split('\t')[5] == '[]':
                tuple_reg.append(i.split('\t')[2])
                tuple_reg.append(i.split('\t')[3])
            else:
                tuple_reg.append(j.replace('(', '').replace(')', '').replace('[', '').replace(']', '').replace(' ', ''))
        it = iter(tuple_reg)
        reg_t=list(zip(it, it))
        for t in reg_t:
            vector = list(range(int(t[0]), int(t[1])+1))
            prova.append(vector)
        if len(prova) > 1:
            new_prova=list(itertools.chain.from_iterable(prova))
            new.append(new_prova)
        else:
            new.append(prova[0])
    return new


def unify_all_units(final_v):
    new=[]
    for i in final_v:
        if not i:
            pass
        else:
            new_list = list(itertools.chain.from_iterable(i))
            new.append(new_list)
    return new


def consensus(R_vector, infos_pdb):
    if not R_vector:
        O_set, M_set, A_set = [], [], []
        return O_set, M_set, A_set
    df = pd.DataFrame(R_vector, columns=['sequence'])
    df['sequence'].replace('\n', ' ')
    df['sequence'].replace('', np.nan, inplace=True)
    majority_array, all_array, one_array = [], [], []
    A_set, M_set, O_set = [], [], []
    vector_consensus = []
    if len(df) == 3:
        consensus_majority = 2
    else:
        consensus_majority = round(len(df) / 2) + 1
    consensus_one = 1
    consensus_all = len(df)
    for i in range(len(df.sequence[0])):
        counter = 0
        for j in range(len(df)):
            if df.sequence[j][i] == '-':
                counter = counter + 1
        vector_consensus.append(len(df) - counter)
    for ind, val in enumerate(vector_consensus):
        if val == 0:
            pass
        elif val == consensus_one:
            vector_consensus[ind] = 'O'
        elif val == consensus_all:
            vector_consensus[ind] = 'A'
        elif val == consensus_majority:
            vector_consensus[ind] = 'M'
        else:
            vector_consensus[ind] = 'M'
    for i, v in enumerate(vector_consensus):
        if v == 'A' or v == 'M' or v == 'O':
            one_array.append(i)
        if v == 'M' or v == 'A':
            majority_array.append(i)
        if v == 'A':
            all_array.append(i)
        else:
            pass
    if (not majority_array) and (not all_array) and (len(infos_pdb)==1):
        majority_array = all_array = one_array
    for k, g in groupby(enumerate(majority_array), lambda ix: ix[0] - ix[1]):
        M_set.append(list(map(itemgetter(1), g)))
    for k, g in groupby(enumerate(all_array), lambda ix: ix[0] - ix[1]):
        A_set.append(list(map(itemgetter(1), g)))
    for k, g in groupby(enumerate(one_array), lambda ix: ix[0] - ix[1]):
        O_set.append(list(map(itemgetter(1), g)))

    return O_set, M_set, A_set, vector_consensus


def missing_residues_mapping(vector_consensus, list_missing_residues, R_vector):
    new_missign_res, m_r_s=[], []
    list_missing_residues.sort()
    x = collections.Counter(list_missing_residues)
    new_list=x.most_common()
    for j in new_list:
        if j[1] == len(R_vector):
            new_missign_res.append(j[0])
    for i,v in enumerate(vector_consensus):
        if (i in list_missing_residues) and (v == 'O'):
            new_missign_res.append(i)
    new_missign_res.sort()
    for k, g in groupby(enumerate(list(set(new_missign_res))), lambda i_x: i_x[0] - i_x[1]):
        m_r_s.append(list(map(itemgetter(1), g)))
    n=[]
    for k, g in groupby(enumerate(sorted(set(list(itertools.chain.from_iterable(m_r_s))))), lambda i_x: i_x[0] - i_x[1]):
        n.append(list(map(itemgetter(1), g)))
    return n


def obtain_period_of_consensus_all(O_vector, units_uniprot_index):
    period_one=[]
    len_units=0
    num_units=0
    #### calcolare la somma di tutte le lunghezze delle unità in consensus one ######
    for i in O_vector:
        for j in units_uniprot_index:
            for k in j:
                if not k:
                    pass
                else:
                    if (k[0] in i) or (k[len(k)-1] in i):
                        len_units=len_units+len(k)
                        num_units +=1
        period_one.append(len_units/num_units)
    return period_one


def obtain_json(unip, unip_seq, unip_title, interpro_coord, O_set, M_set, A_set, class_and_sub, infos_pdb, unit_unip, ins_unip, reg_unip, curator, period_one, miss_residues, useful_vector):

    uniprot_json = {'uniprot_id': unip.rstrip(), 'uniprot_sequence': unip_seq.rstrip(), 'uniprot_name': unip_title.rstrip()}

    if not interpro_coord:
        pass
    else:
        children_interpro = []
        for b in interpro_coord:
            children_interpro.append({'id': b[0], 'start': b[2], 'end': b[3].rstrip()})
        uniprot_json["interpro"] = children_interpro

        children_pfam = []
        for b in interpro_coord:
            if b[1].startswith('PF'):
                children_pfam.append({'id': b[1], 'start': b[2], 'end': b[3].rstrip()})
        if not children_pfam:
            pass
        else:
            uniprot_json['pfam'] = children_pfam

        pfam_cons=[]
        for b in interpro_coord:
            if b[1].startswith('PF'):
                for c in O_set:
                    if (int(b[2]) in c) or (int(b[3].rstrip()) in c):
                        pfam_cons.append(b[1])
        if not pfam_cons:
            pass
        elif len(pfam_cons) == 1:
            uniprot_json['pfam_consensus'] = pfam_cons
        else:
            uniprot_json['pfam_consensus'] = list(set(pfam_cons))

    children_pdb = []
    for w in range(len(infos_pdb)):
        children_pdb.append({'id': infos_pdb[w].split('\t')[0], 'seqres_sequence': infos_pdb[w].split('\t')[3], 'start': infos_pdb[w].split('\t')[1], 'end': infos_pdb[w].split('\t')[2], 'origin': infos_pdb[w].split('\t')[4]})
    uniprot_json["pdb_chains"] = children_pdb

    children_consensus_one = []
    if (len(O_set) == 1) and (len(class_and_sub) > 1):
        try:
            for ln, p in zip(range(len(O_set)), range(len(period_one))):
                if not miss_residues:
                    children_consensus_one.append(
                        {'start': O_set[ln][0]+1, 'end': O_set[ln][len(O_set[ln]) - 1]+1, 'classification': class_and_sub, 'period': round(period_one[p], 2), 'ambiguity': 'multiple classification'})
                else:
                    children_consensus_one.append(
                        {'start': O_set[ln][0] + 1, 'end': O_set[ln][len(O_set[ln]) - 1] + 1,
                         'classification': class_and_sub, 'period': round(period_one[p], 2),
                         'ambiguity': ['multiple classification', 'missing residues'],
                         'missing_residues': miss_residues})
        except TypeError:
            for ln in (range(len(O_set))):
                if not miss_residues:
                    children_consensus_one.append(
                        {'start': O_set[ln][0]+1, 'end': O_set[ln][len(O_set[ln]) - 1]+1, 'classification': class_and_sub, 'period': round(period_one, 2), 'ambiguity': 'multiple classification'})
                else:
                    children_consensus_one.append(
                        {'start': O_set[ln][0] + 1, 'end': O_set[ln][len(O_set[ln]) - 1] + 1,
                         'classification': class_and_sub, 'period': round(period_one, 2),
                         'ambiguity': ['multiple classification', 'missing residues'],
                         'missing_residues': miss_residues})
    else:
        try:
            for ln, p in zip(range(len(O_set)), range(len(period_one))):
                class_in_cons_one, new_class_in_cons_one = [], []
                for useful in useful_vector:
                    if int(useful.split('\t')[2]) <= int(O_set[ln][0]+1) <= int(useful.split('\t')[3]):
                        class_in_cons_one.append(useful.split('\t')[4].replace('[', '').replace(']', '').replace("'", ''))
                uniq_class_in_cons_one=list(set(class_in_cons_one))
                for k in uniq_class_in_cons_one:
                    try:
                        a,b = k.split(',', 1)
                        new_class_in_cons_one.append(a)
                        new_class_in_cons_one.append(b.replace(' ', ''))
                    except:
                        new_class_in_cons_one.append(k)
                if (not miss_residues) or (miss_residues[0][0] > O_set[len(O_set) - 1][len(O_set[len(O_set) - 1]) - 1]) or (int(miss_residues[0][0]) < int(O_set[ln][0])):
                    children_consensus_one.append(
                        {'start': O_set[ln][0]+1, 'end': O_set[ln][len(O_set[ln]) - 1]+1, 'classification': list(set(new_class_in_cons_one)), 'period': round(period_one[p], 2)})
                else:
                    children_consensus_one.append(
                        {'start': O_set[ln][0] + 1, 'end': O_set[ln][len(O_set[ln]) - 1] + 1,
                         'classification': new_class_in_cons_one, 'period': round(period_one[p], 2),
                         'ambiguity': 'missing residues', 'missing_residues':miss_residues})
        except TypeError:
            for ln in(range(len(O_set))):
                if (not miss_residues) or (miss_residues[0][0] > O_set[len(O_set) - 1][len(O_set[len(O_set) - 1]) - 1]):
                    children_consensus_one.append(
                        {'start': O_set[ln][0]+1, 'end': O_set[ln][len(O_set[ln]) - 1]+1, 'classification': class_and_sub, 'period': round(period_one, 2)})
                elif (not miss_residues) or (int(miss_residues[0][0]) < int(O_set[ln][0])):
                    children_consensus_one.append(
                        {'start': O_set[ln][0] + 1, 'end': O_set[ln][len(O_set[ln]) - 1] + 1,
                         'classification': class_and_sub, 'period': round(period_one, 2)})
                else:
                    children_consensus_one.append(
                        {'start': O_set[ln][0] + 1, 'end': O_set[ln][len(O_set[ln]) - 1] + 1,
                         'classification': class_and_sub, 'period': round(period_one, 2),
                         'ambiguity': 'missing residues', 'missing_residues': miss_residues})
    uniprot_json["repeatsdb_consensus_one"] = children_consensus_one

    children_consensus_majority = []
    for fl in range(len(M_set)):
        children_consensus_majority.append(
            {'start': M_set[fl][0]+1, 'end': M_set[fl][len(M_set[fl]) - 1]+1, 'classification': class_and_sub})
    uniprot_json["repeatsdb_consensus_majority"] = children_consensus_majority

    children_consensus_all = []
    for al in range(len(A_set)):
        children_consensus_all.append(
            {'start': A_set[al][0]+1, 'end': A_set[al][len(A_set[al]) - 1]+1, 'classification': class_and_sub})
    uniprot_json["repeatsdb_consensus_all"] = children_consensus_all

    if (not unit_unip) and (not reg_unip) and (not curator):
        pass
    else:
        children_repeatsdb_annotation = {'units': unit_unip, 'insertions': ins_unip, 'regions': reg_unip, 'curator': curator}
        uniprot_json['repeatsdb_annotation'] = children_repeatsdb_annotation

    exon_information = uniprot_exons(unip.rstrip(), folder='./', dest=None)
    if not exon_information:
        pass
    else:
        exon_info=[]
        for exon in exon_information:
            exon_info.append({'id': exon, 'isoform': exon_information[exon][0]['isoform'], 'location': exon_information[exon][0]['location'], 'start': int(exon_information[exon][0]['start']), 'end': int(exon_information[exon][0]['end'])})
        uniprot_json["exon_mapping"] = exon_info
    return uniprot_json


def to_mongo(uniprot_json):
    myclient = pymongo.MongoClient("mongodb://localhost:27017/")
    mydb = myclient["RepeatsDB"]
    mycol = mydb["uniprot"]
    x = mycol.insert_one(uniprot_json)  # .inserted_id

