# mapping script
import requests
import argparse
import json
from ftplib import FTP
import gzip
from xml.dom.minidom import parseString

aa_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
aa_d_r = {aa_d[key]: key for key in aa_d}


def parse_sifts(siftsfile, chain, acc=None):
    def get_one(r, attribute, acc=acc):
        for ref in r.getElementsByTagName("crossRefDb"):
            if ref.getAttribute('dbSource') == attribute:
                if acc and ref.getAttribute('dbAccessionId') == acc:
                    return int(ref.getAttribute('dbResNum'))
                if ref.getAttribute('dbResNum') != "null":
                    return ref.getAttribute('dbResNum')
        return None

    # open and read gzipped xml file
    infile = gzip.open(siftsfile)
    content, references = infile.read(), {}
    # parse xml file content
    dom = parseString(content)
    entities = dom.getElementsByTagName("entity")
    for e in entities:
        if e.getAttribute('entityId') == chain:
            segments = e.getElementsByTagName("segment")
            for s in segments:
                # select necessary references
                refs = {get_one(r, "UniProt", acc=acc): int(r.getAttribute("dbResNum"))
                        for r in s.getElementsByTagName("residue") if get_one(r, "PDB")}
                references.update(refs)
    return references


def download_sifts(pdb, chain, acc, folder=''):
    ftp_folder = pdb[1:3]
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd('pub/databases/msd/sifts/split_xml/' + ftp_folder)
    flname = folder + '{}.xml.gz'.format(pdb)
    ftp.retrbinary('RETR {}.xml.gz'.format(pdb), open(flname, 'wb').write)
    return parse_sifts(flname, chain, acc)


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


def uniprot_summary(uniprot, folder='', filter=None):
    dest = folder + uniprot + '_summary.tsv'
    outfl = open(dest, 'w+')
    # headers
    headers, seq_uniprot, maxlen, maxcoverage = [uniprot], [], 0, ('', 0)
    mapped_pdbs = json.loads(map_accession(uniprot))[uniprot]["PDB"]

    sequence = download_sequence(uniprot)["sequence"]["sequence"]
    entities = {}
    # filter
    if filter:
        filtered_mapped_pdbs = {p: mapped_pdbs[p] for p in mapped_pdbs if p in filter}
        mapped_pdbs = filtered_mapped_pdbs
    # compile mappings
    for key in mapped_pdbs:
        e = mapped_pdbs[key]
        residues = residue_list(key)
        for match in e:
            if key + match["chain_id"] not in entities:
                entities[key + match["chain_id"]] = []
            # select right chain and residues
            pdb_start, pdb_end, matched_sequence = match["start"]["residue_number"], match["end"]["residue_number"], {}
            # correspond from match["unp_start"] to match["unp_end"]
            if match["chain_id"] in residues[match["entity_id"]]:
                matched_residues = [
                    r for r in residues[match["entity_id"]][match["chain_id"]]
                    if pdb_start <= r["residue_number"] <= pdb_end
                ]
                uniprot_positions = range(match["unp_start"], match["unp_end"] + 1)
                if len(uniprot_positions) != len(matched_residues):
                    # mapping gap, investigate it with SIFTS
                    sifts = download_sifts(key, match["chain_id"], uniprot, folder=folder)
                    matched = {r["residue_number"]: r for r in residues[match["entity_id"]][match["chain_id"]]}
                    seq = {key: matched[sifts[key]] for key in sifts}
                else:
                    seq = {uniprot_positions[index]: matched_residues[index] for index in range(len(uniprot_positions))}
                entities[key + match["chain_id"]].append({
                    'start': match["unp_start"],
                    'end': match["unp_end"],
                    'seq': seq
                })
                headers.append(key + match["chain_id"])
    # output file uniprot summary
    outfl.write('\t'.join(headers) + '\n')
    # write output file
    for i in range(1, len(sequence) + 1):
        outfl.write(str(i) + ',' + aa_d_r[sequence[i - 1]])
        for h in headers[1:]:
            entries, seqs = entities[h], {}  # is a list
            for e in entries:
                seqs.update(e['seq'])
            # outside ranges
            if i not in seqs:
                outfl.write("\t-")
            else:
                res = seqs[i]
                outfl.write("\t" + str(res['author_residue_number']) + ',' + res['residue_name'])
        outfl.write('\n')


def process_pdb(pdb_input, rootdir):
    pdblist = pdb_input.split(',')
    pdb_filter = None
    if len(pdblist) > 1:
        pdb = pdblist[0]
        pdb_filter = pdblist
    else:
        pdb = pdblist[0]
    pdb_mapped = json.loads(map_accession(pdb))[pdb]["UniProt"]
    if pdb_mapped != {}:
        uniprots = pdb_mapped.keys()
        for u in uniprots:
            print('PDB', pdb, 'is mapped to UniProt', u)
            uniprot_summary(u, folder=rootdir + '/output_collections/mappings/', filter=pdb_filter)
    else:
        print('PDB was not mapped to any UniProt ACC')
        return None
