import os, sys
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning
from string import whitespace
import subprocess
import pandas as pd
from urllib.request import urlretrieve
import matplotlib.pyplot as plt
import seaborn as sns
import config as cfg
from datetime import date
import gzip
import shutil

warnings.simplefilter('ignore', BiopythonWarning)


def save_pdb_structure_chain_extract_residues(pdbname, pdbchain, folder):
    io = PDBIO()
    if os.path.exists(f'/mnt/db/pdb/{folder}/pdb{pdbname}.ent.gz'):
        with gzip.open(f'/mnt/db/pdb/{folder}/pdb{pdbname}.ent.gz', 'rb') as pdb_file:
            with open(f'{cfg.matrix_path["matrix"]}/temp/{pdbname}.pdb', 'wb') as pdb_out:
                shutil.copyfileobj(pdb_file, pdb_out)
        pdb = PDBParser().get_structure(pdbname, f'{cfg.matrix_path["matrix"]}/temp/{pdbname}.pdb')
        ppb = PPBuilder()
        io.set_structure(pdb)
        for chain in pdb.get_chains():
            if pdbchain == chain.get_id():
                io.set_structure(chain)
                io.save(f'{cfg.matrix_path["matrix"]}/temp/{pdb.get_id()}_{pdbchain}.pdb')
            else:
                pass
        chain_pdb = PDBParser().get_structure(f'{pdbchain}_{pdbname}', f'{cfg.matrix_path["matrix"]}/temp/{pdb.get_id()}_{pdbchain}.pdb')
        chain_pdb_file = open(f'{cfg.matrix_path["matrix"]}/temp/{pdb.get_id()}_{pdbchain}.pdb').readlines()
        residues_chain = [residue.id[1] for residue in chain_pdb.get_residues()]
    else:
        print(f'This file {pdbname}{pdbchain} has no .PDB FILE')
        residues_chain, chain_pdb_file = [], []
        exit()

    resseqs = list(map(str, residues_chain))
    return resseqs, chain_pdb_file


def main():
    count = 0
    #todo: use json files
    print(f'Create matrices for all the .db files... {len(os.walk("/mnt/projects/RepeatsDB/repeatsdb3/data/db_files").__next__()[2])}')
    for name in os.listdir('/mnt/projects/RepeatsDB/repeatsdb3/data/db_files'):
        # name='5b0mF.db'
        # name = '2ozqA.db'
        name = '1bpoA.db' #pdb with 2 regions
        # name='3o3qD.db'
        count += 1
        if count % 100 == 0:
            print(f'COUNT: ...{count:,}... {name}')
        folder = name.split('.')[0][1:3]
        pdbname = name.split('.')[0][:4]
        pdbchain = name.split('.')[0][4:]

        try:
            os.makedirs(f'{cfg.matrix_path["matrix"]}/matrix_{date.today()}/{folder}')
        except FileExistsError:
            pass

        file_db = open(f'/mnt/projects/RepeatsDB/repeatsdb3/data/db_files/{name.rstrip()}')

        if not file_db:
            print(f'This file {name.rstrip()} has no .DB FILE')
            exit()

        try:
            os.makedirs(f'{cfg.matrix_path["matrix"]}/temp')
        except FileExistsError:
            pass

        resseqs, chain_pdb_file = save_pdb_structure_chain_extract_residues(pdbname, pdbchain, folder)

        residue_length, regions, units, insertion = [], [], [], []
        for j in file_db:
            if j.startswith('REG'):
                region = ('{}'.format(j.split('\t')[1].split(' ')[0]), '{}'.format(j.split('\t')[1].split(' ')[1]))
                regions.append(region)
            elif j.startswith('UNIT'):
                start = j.split('\t')[1].split(' ')[0]
                end = j.split('\t')[1].split(' ')[1].rstrip()
                unit_one= ('{}'.format(start), '{}'.format(end))
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
                os.mkdir(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}')
            except FileExistsError:
                pass
            units_in_region, ins_in_region, highlighted_unit = [], [], []
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
                unit_db=open(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}/{name.split(".")[0]}_{units_in_region[l][0]}_{units_in_region[l][1]}.pdb', 'w+')
                if (resseqs[resseqs.index(str(units_in_region[l][0]))]== resseqs[resseqs.index(str(units_in_region[l][0]))]) and (resseqs[resseqs.index(str(units_in_region[l][1]))] == resseqs[resseqs.index(str(units_in_region[l][1]))]):
                    result = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(str(units_in_region[l][1]))+1]
                elif (resseqs[resseqs.index(str(units_in_region[l][0]))]== resseqs[resseqs.index(str(units_in_region[l][0]) +'A')]) and (resseqs[resseqs.index(str(units_in_region[l][0]))]== resseqs[resseqs.index(str(units_in_region[l][0])+'A')]):
                    result = resseqs[resseqs.index(str(units_in_region[l][0]) + 'A'): resseqs.index(str(units_in_region[l][0]) + 'A') + 1]
                elif (resseqs[resseqs.index(str(units_in_region[l][0]))]== resseqs[resseqs.index(str(units_in_region[l][0])+'A')]) and (resseqs[resseqs.index(str(units_in_region[l][0]))]== resseqs[resseqs.index(str(units_in_region[l][0]))]):
                    result = resseqs[resseqs.index(str(units_in_region[l][0]) + 'A'): resseqs.index(str(units_in_region[l][0])) + 1]
                elif (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[resseqs.index(str(units_in_region[l][0]))]) and (resseqs[resseqs.index(str(units_in_region[l][0]))] == resseqs[resseqs.index(str(units_in_region[l][0]) + 'A')]):
                    result = resseqs[resseqs.index(str(units_in_region[l][0])): resseqs.index(str(units_in_region[l][0]) + 'A') + 1]
                else:
                    result=[]
                    print(resseqs)
                    print(name)
                    exit()
                for g in range(len(result)):
                    for line in chain_pdb_file:
                        if line[22:26] != '':
                            if line[0:6] == "ATOM  ":
                                if line[22:26].translate(dict.fromkeys(map(ord, whitespace))) == result[g]:
                                    unit_db.write(line)
                unit_db.close()
                if not ins_in_region:
                    pass
                else:
                    for t in range(len(ins_in_region)):
                        if ins_in_region[t][0] in result:
                            if units_in_region[l] in highlighted_unit:
                                pass
                            else:
                                highlighted_unit.append(units_in_region[l])
                        else:
                            pass
            tm_score=[]
            #just for particular cases#
            # k=0
            ###########################
            for l in os.listdir(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}'):
                for m in os.listdir(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}'):
                    print(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}')
                    cmd = f'./../TMalign-master/TMalign {cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}/{l} {cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}/{m}'.split()
                    with subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding="utf-8") as proc:
                        out = proc.stdout.readlines()
                    len_l=int(l.split('_')[1]) - int(l.split('_')[1])
                    len_m = int(m.split('_')[1]) - int(m.split('_')[1])
                    if len_l >= len_m:
                        for line in out:
                            print(line)
                            if line.startswith('TM-score') and 'Structure_2' in line:
                                tm_score.append('{0}\t{1}\t{2}'.format(l, m, float(line.split('=')[1].split('(')[0])))
                    else:
                        for line in out:
                            if line.startswith('TM-score') and 'Structure_2' in line:
                                tm_score.append('{0}\t{1}\t{2}'.format(l, m, float(line.split('=')[1].split('(')[0])))
            tm_matrix = pd.DataFrame([ln.split('\t') for ln in tm_score])
            tm_matrix.columns = ['unit_1', 'unit_2', 'TMscore']
            tm_matrix.to_csv(
                f'{cfg.matrix_path["matrix"]}/data_10_09/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}_matrix.tsv', sep='\t')

            tm_matrix.iloc[:,0] = tm_matrix.iloc[:,0].map(lambda x: x.lstrip('{}'.format(name.split('.')[0])).rstrip('.pdb'))
            tm_matrix.iloc[:,1] = tm_matrix.iloc[:,1].map(lambda x: x.lstrip('{}'.format(name.split('.')[0])).rstrip('.pdb'))
            tm_matrix.iloc[:,0] = tm_matrix.iloc[:,0].str.replace('^_', ' ').str.replace('_', '-')
            tm_matrix.iloc[:,1] = tm_matrix.iloc[:,1].str.replace('^_', ' ').str.replace('_', '-')

            prova = tm_matrix.pivot(index='unit_1', columns = 'unit_2', values = 'TMscore').reset_index().fillna(0)
            prova['sort'] = '0'
            for p in range(len(prova.unit_1)):
                start_unit=prova.unit_1[p].split('-')[0]
                if start_unit == ' ':
                    start_unit = prova.unit_1[p].split('-')[1]
                prova.sort[p] = prova.sort[p].replace(prova.sort[p], start_unit)
            prova['sort'] = prova['sort'].astype(int)
            prova.sort_values('sort', inplace=True, ascending=True)
            prova.reset_index(inplace=True)
            prova=prova.drop(['sort', 'index'], axis=1)
            prova.loc[-1] = prova.columns.values
            prova.index = prova.index + 1
            for s in range(1, len(prova.unit_1)):
                if prova.loc[0][s].split('-')[0] == ' ':
                    prova.loc[0][s] = prova.loc[0][s].replace(prova.loc[0][s], prova.loc[0][s].split('-')[1])
                else:
                    prova.loc[0][s] = prova.loc[0][s].replace(prova.loc[0][s], prova.loc[0][s].split('-')[0])
            prova.loc[0]['unit_1'] = prova.loc[0]['unit_1'].replace(prova.loc[0]['unit_1'], '0')
            prova.loc[0] = prova.loc[0].astype(float)
            prova = prova.sort_values(0, axis=1)
            prova = prova.drop(prova[prova.unit_1 == 0].index)
            prova.set_index('unit_1', inplace=True)
            prova = prova.astype(float)

            f, ax = plt.subplots(figsize=(16, 7))
            for xtick in ax.xaxis.get_major_ticks():
                xtick.label.set_fontsize(14)
            for ytick in ax.yaxis.get_major_ticks():
                ytick.label.set_fontsize(14)

            if len(prova.columns) > 10:
                font_size_cell = 10
            elif len(prova.columns) > 15:
                font_size_cell = 7
            else:
                font_size_cell = 17
            sns.heatmap(prova, annot=True, linewidths=.5, cmap="BuPu", vmin=0, vmax=1, cbar=False,
                        annot_kws={'fontsize':font_size_cell}, square=True)
            for h in highlighted_unit:
                label=f'{int(h[0])}-{int(h[1])}'
                print(label)
                wanted_index=prova.columns.get_loc(f' {label}')
                ax.get_xticklabels()[wanted_index].set_color("#f9a516")
                ax.get_yticklabels()[wanted_index].set_color("#f9a516")

            plt.yticks(rotation='horizontal')
            plt.xticks(rotation='vertical')

            plt.xlabel('')
            plt.ylabel('')
            plt.savefig(f'{cfg.matrix_path["matrix"]}/matrix_{date.today()}/{folder}/{name.split(".")[0]}_{regions[k][0]}_{regions[k][1]}_matrix.svg', format='svg', dpi=1200, bbox_inches = 'tight', transparent=True) #######change the directory in wich you want to save the matrix if there is a second region (os.path.basename(os.path.normpath(mustang_align)) != chains)
            plt.close()
            print(f'CREATE MATRIX {name.split(".")[0]} for the region {regions[k][0]}_{regions[k][1]}')

            shutil.rmtree(f'{cfg.matrix_path["matrix"]}/temp/temp_{regions[k][0]}_{regions[k][1]}')
        shutil.rmtree(f'{cfg.matrix_path["matrix"]}/temp')
        # exit()
    total_len=len(os.walk(f'{cfg.matrix_path["matrix"]}/matrix_{date.today()}').__next__()[2])
    print(f'Created {total_len} matrices')

if __name__ == '__main__':
    sys.exit(main())
