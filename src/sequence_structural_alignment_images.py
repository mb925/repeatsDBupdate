import gzip
import os, time
import subprocess
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
from string import whitespace
import subprocess
import pandas as pd
from urllib.request import urlretrieve
import shutil
from shutil import copyfile
from operator import itemgetter
from itertools import groupby
import sys

import config as cfg

warnings.simplefilter('ignore', BiopythonWarning)

path = cfg.data_path['data'] + 'multiple_structural_align/'

def mustang_alignment(name, regions, folder): # todo
    cmd = './../MUSTANG_v3.2.3/bin/mustang-3.2.3 ' \
          '-f ./../data/multiple_structural_align/temp/{0}_{1}_{2}/{0}_{1}_{2}_mustang_description_file ' \
          '-o ./../data/multiple_structural_align/data/{3}/{0}_{1}_{2} ' \
          '-F fasta html ' \
          '-s ON'.format(name.split(".")[0], regions[k][0], regions[k][1], folder).split()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")

    while proc.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    return_code = proc.returncode
    return return_code


def mtmalign_alternative_mustang(name, regions):
    # print('we are here')
    with open(
            f'../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[0]}_{regions[1]}/{name.split(".")[0]}_{regions[0]}_{regions[1]}_mtmalign_description_file',
            'w+') as mtmalign_conf_file:
        for j in os.listdir(f'./../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[0]}_{regions[1]}'):
            if j.endswith('.pdb'):
                mtmalign_conf_file.write(f'{j}\n')
    os.chdir(f'./../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[0]}_{regions[1]}')
    cmd = './../mTM-align/src/mTM-align -i {}_{}_{}_mtmalign_description_file'.format(name.split(".")[0], regions[0], regions[1]).split() # todo
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")

    while proc.poll() is None:
        # Process hasn't exited yet, let's wait some
        time.sleep(0.5)

    return_code = proc.returncode
    return return_code


count = 0
for name in os.listdir('/mnt/projects/RepeatsDB/repeatsdb3/data/db_files'):
    # with open('/home/sara/Documents/pre_submission/multiple_structural_align_mustang/prova.txt') as problems:
    #     for prob in problems:
    #     name = prob.split(' ')[2]
    #     name='6crzB.db'
    # name = '3r8yE.db'
    #     name='5ojfA.db'

    # for name in os.listdir('/home/sara/Documents/pre_submission/cartella_x'):
    folder = name.split('.')[0][1:3]
    pdbname = name.split('.')[0][:4]
    pdbchain = name.split('.')[0][4:]
    count += 1

    if (name == '1swgA.db') or (name == '1swfB.db') or (name == '4gyxC.db') or (name == '1swfD.db') or (
            name == '1swgB.db') or (name == '4gyxA.db') or (name == '1swfC.db') or (name == '4gyxB.db') or (
            name == '1swgD.db'):
        pass
    else:
        if count % 100 == 0:
            print(f'COUNT: ...{count:,}... {name}')
        try:
            os.makedirs('./../data/multiple_structural_align/temp/chain_pdb')
        except FileExistsError:
            pass
        try:
            os.makedirs(f'./../data/multiple_structural_align/data/{folder}')
        except FileExistsError:
            pass

        file_db = open(f'/mnt/projects/RepeatsDB/repeatsdb3/data/db_files/{name.rstrip()}')
        io = PDBIO()
        # try:
        if f'pdb{pdbname}.ent.gz' not in os.listdir(f'/mnt/db/pdb/{folder}/'):
            print(f'No .pdb file in our database for {name}')
            resseqs = []
        else:
            # url = urlretrieve(f'http://files.rcsb.org/download/{pdbname}.pdb',
            #                   f'/home/sara/Documents/new_SRUL/data/SRUL_REPEATSDB/DATABASE_PDB/{pdbname}.pdb')
            with gzip.open(f'/mnt/db/pdb/{folder}/pdb{pdbname}.ent.gz', 'rb') as pdb_file:
                with open(f'./../data/multiple_structural_align/temp/{pdbname}.pdb', 'wb') as pdb_out:
                    shutil.copyfileobj(pdb_file, pdb_out)
            pdb = PDBParser().get_structure(pdbname,
                                            f'./../data/multiple_structural_align/temp/{pdbname}.pdb')
            ppb = PPBuilder()
            io.set_structure(pdb)
            for chain in pdb.get_chains():
                if pdbchain == chain.get_id():
                    io.set_structure(chain)
                    io.save(
                        f'./../data/multiple_structural_align/temp/chain_pdb/{pdb.get_id()}_{pdbchain}.pdb')
                else:
                    pass
            chain_pdb = PDBParser().get_structure(f'{pdbname}_{pdbchain}',
                                                  f'./../data/multiple_structural_align/temp/chain_pdb/{pdbname}_{pdbchain}.pdb')
            chain_pdb_file = open(
                f'./../data/multiple_structural_align/temp/chain_pdb/{pdbname}_{pdbchain}.pdb').readlines()
            residues_chain = [residue.id[1] for residue in chain_pdb.get_residues()]
            resseqs = list(map(str, residues_chain))
        # except:
        #     print(f'HTTP Error 404: Not Found {name}')
        #     resseqs = []

        # pdb = PDBParser().get_structure(pdbname,
        #                                 f'/home/sara/Documents/new_SRUL/data/SRUL_REPEATSDB/DATABASE_PDB/{pdbname}.pdb')
        # ppb = PPBuilder()
        # io.set_structure(pdb)
        # for chain in pdb.get_chains():
        #     if pdbchain == chain.get_id():
        #         io.set_structure(chain)
        #         io.save(f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/chain_pdb/{pdb.get_id()}_{pdbchain}.pdb')
        #     else:
        #         pass
        # chain_pdb = PDBParser().get_structure(f'{pdbname}_{pdbchain}',
        #                                       f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/chain_pdb/{pdbname}_{pdbchain}.pdb')
        # chain_pdb_file = open(f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/chain_pdb/{pdbname}_{pdbchain}.pdb').readlines()
        # residues_chain = [residue.id[1] for residue in chain_pdb.get_residues()]
        # resseqs = list(map(str, residues_chain))
        residue_length, regions, units, insertion = [], [], [], []
        if not resseqs:
            pass
        else:
            for j in file_db:
                if j.startswith('REG'):
                    region = ('{}'.format(j.split('\t')[1].split(' ')[0]), '{}'.format(j.split('\t')[1].split(' ')[1]))
                    regions.append(region)
                elif j.startswith('UNIT'):
                    start = j.split('\t')[1].split(' ')[0]
                    end = j.split('\t')[1].split(' ')[1].rstrip()
                    unit_one = (f'{start}', f'{end}')
                    units.append(unit_one)

                elif j.startswith('INS'):
                    start = j.split('\t')[1].split(' ')[0]
                    end = j.split('\t')[1].split(' ')[1].rstrip()
                    ins_one = (f'{start}', f'{end}')
                    insertion.append(ins_one)

                else:
                    pass
            for k in range(len(regions)):
                try:
                    os.makedirs(
                        f'./../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}')
                except FileExistsError:
                    pass
                # print('ciao')
                # if f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}.html' and f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}.afasta' \
                #         and f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}.pdb' in os.listdir(f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/data/{folder}'):
                #     # print(f'The .html, .pdb and .afasta for the protein {name} already exist.')
                #     pass
                # else:
                units_in_region, ins_in_region = [], []
                for i in range(len(units)):
                    if int(regions[k][0]) <= int(units[i][0]) < int(regions[k][1]):
                        units_in_region.append(units[i])
                if not insertion:
                    pass
                else:
                    for r in range(len(insertion)):
                        for u in range(len(units_in_region)):
                            if int(units_in_region[u][0]) <= int(insertion[r][0]) < int(units_in_region[u][1]):
                                ins_in_region.append(insertion[r])
                for l in range(len(units_in_region)):
                    unit_db = open(f'./../data/multiple_structural_align/temp/'
                                   f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}/'
                                   f'{name.split(".")[0]}_{units_in_region[l][0]}_{units_in_region[l][1]}.pdb',
                                   'w+')
                    # if not ins_in_region:
                    # for l in range(len(units_in_region)):
                    #     unit_db = open(f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/'
                    #                    f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}/{name.split(".")[0]}_{units_in_region[l][0]}_{units_in_region[l][1]}.pdb',
                    #                    'w+')
                    if (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                        resseqs.index(str(units_in_region[l][0]))]) and \
                            (resseqs[resseqs.index(str(units_in_region[l][1]))] == resseqs[
                                resseqs.index(str(units_in_region[l][1]))]):
                        result = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(
                            str(units_in_region[l][1])) + 1]
                    elif (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                        resseqs.index(str(units_in_region[l][0]) + 'A')]) and \
                            (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                                resseqs.index(str(units_in_region[l][0]) + 'A')]):
                        result = resseqs[resseqs.index(str(units_in_region[l][0]) + 'A'): resseqs.index(
                            str(units_in_region[l][0]) + 'A') + 1]
                    elif (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                        resseqs.index(str(units_in_region[l][0]) + 'A')]) and \
                            (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                                resseqs.index(str(units_in_region[l][0]))]):
                        result = resseqs[resseqs.index(str(units_in_region[l][0]) + 'A'): resseqs.index(
                            str(units_in_region[l][0])) + 1]
                    elif (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                        resseqs.index(str(units_in_region[l][0]))]) and \
                            (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[
                                resseqs.index(str(units_in_region[l][0]) + 'A')]):
                        result = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(
                            str(units_in_region[l][0]) + 'A') + 1]
                    else:
                        result = []
                        print(resseqs)
                        print(name)
                        exit()
                    # print(result)
                    if not ins_in_region:
                        for g in range(len(result)):
                            for line in chain_pdb_file:
                                if line[22:26] != '':
                                    if line[0:6] == "ATOM  ":
                                        if line[22:26].translate(dict.fromkeys(map(ord, whitespace))) == result[g]:
                                            unit_db.write(line)
                        unit_db.close()
                        # pass
                    else:
                        for t in range(len(ins_in_region)):
                            if (ins_in_region[t][0] in result) or (ins_in_region[t][1] in result):
                                try:
                                    result_1 = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(
                                        str(ins_in_region[t][0]))]
                                    result_2 = resseqs[resseqs.index(str(ins_in_region[t][1])) + 1: resseqs.index(
                                        str(units_in_region[l][1])) + 1]
                                    result = result_1 + result_2
                                except ValueError:
                                    int_list = [int(x) for x in result]
                                    missing_values = []
                                    missing_res = []
                                    for miss in range(min(int_list), max(int_list)):
                                        if not miss in int_list:
                                            missing_values.append(miss)
                                    for w, g in groupby(enumerate(missing_values), lambda i_x: i_x[0] - i_x[1]):
                                        mist = (list(map(itemgetter(1), g)))
                                        missing_res.append(mist)
                                    for m in missing_res:
                                        if ins_in_region[t][0] in m:
                                            result_1 = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(
                                                str(m[0])) - 1]
                                            result_2 = resseqs[
                                                       resseqs.index(str(ins_in_region[t][1])) + 1: resseqs.index(
                                                           str(units_in_region[l][1])) + 1]
                                            result = result_1 + result_2
                                        elif ins_in_region[t][1] in m:
                                            result_1 = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(
                                                str(ins_in_region[t][0]))]
                                            result_2 = resseqs[resseqs.index(str(m[0])) + 1: resseqs.index(
                                                str(units_in_region[l][1])) + 1]
                                            result = result_1 + result_2
                                        else:
                                            pass
                        # print(result)

                        for g in range(len(result)):
                            for line in chain_pdb_file:
                                if line[22:26] != '':
                                    if line[0:6] == "ATOM  ":
                                        if line[22:26].translate(dict.fromkeys(map(ord, whitespace))) == result[g]:
                                            unit_db.write(line)
                        unit_db.close()

                ################################ MUSTANG ################################
                description_file = open(f'./../data/multiple_structural_align/temp/'
                                        f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}/'
                                        f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}_mustang_description_file',
                                        'w+')
                description_file.write(
                    f'>./../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}\n')
                for file in os.listdir(f'./../data/multiple_structural_align/temp/'
                                       f'{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}'):
                    if file.endswith('.pdb') and file.startswith(pdbname):
                        description_file.write(f'+ {file}\n')
                description_file.close()
                output_code = mustang_alignment(name, regions, folder) # todo: check
                # print(output_code)
                if output_code == 0:
                    shutil.rmtree(
                        f'./../data/multiple_structural_align/temp/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}')
                    pass
                else:
                    ############################# MTM-align #############################
                    print(
                        f'This protein {name} has problem with mustang ({output_code}). Aligning the protein with mTM-align......{count}')
                    mtm_output_code = mtmalign_alternative_mustang(name, regions[k])
                    # print(mtm_output_code)
                    if mtm_output_code == 0:

                        f = f'./mTM_result/cc.pdb'
                        with open(f, 'r') as contact_file:
                            contact = 0
                            for ln in contact_file:
                                if ln.startswith('ATOM'):
                                    contact += 1
                                else:
                                    pass
                            if contact > 0:
                                cur = os.getcwd()
                                shutil.copy('./mTM_result/result.fasta',
                                         f'./mTM_result/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}.afasta')
                                shutil.copy('./mTM_result/result.pdb',
                                         f'./mTM_result/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}.pdb')
                            else:
                                print(f'NO contacts for this protein {name}')
                                pass
                    else:
                        print(
                            f'This protein {name} has problem both mustang {output_code} and mTM-align {mtm_output_code}......{count}')
                        pass

                    shutil.rmtree(cur)
            exit()
            # print(regions)
            # os.system(f'mustang '
            #           f'-f /home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}_mustang_description_file '
            #           f'-o /home/sara/Documents/pre_submission/multiple_structural_align_mustang/data/{folder}/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]} '
            #           f'-F fasta html '
            #           f'-s ON > 2')
            # shutil.rmtree(f'/home/sara/Documents/pre_submission/multiple_structural_align_mustang/temp/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}')
    # print(os.getcwd())
    # exit()
