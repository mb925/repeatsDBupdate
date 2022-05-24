import requests
import argparse
import json


def download(url):
    r = requests.get(url, headers={"Accept": "application/json"})
    if not r.ok:
        # r.raise_for_status()
        return None
    else:
        return r.text


def parse_exons(data):
    exons, isoform = [], 1
    locations = [c["genomicLocation"] for c in data["gnCoordinate"] if "gnCoordinate" in data if "genomicLocation" in c]
    for l in locations:
        if 'exon' in l:
            for e in l['exon']:
                rev = '-1' if l['reverseStrand'] else '1'
                # proteinlocation does not have start and end if length = 1
                exons.append({
                    'isoform': str(isoform),
                    'location': str(l['chromosome'] + '_' + str(l['start']) + '_' + str(
                        l['end']) + '_' + rev) if 'chromosome' in l else str(l['start']) + '_' + str(
                        l['end']) + '_' + rev,
                    'id': e['id'],
                    'start': str(e['proteinLocation']['begin']['position']) if 'begin' in e['proteinLocation'] else
                    e['proteinLocation']['position']['position'],
                    'end': str(e['proteinLocation']['end']['position']) if 'end' in e['proteinLocation'] else
                    e['proteinLocation']['position']['position']
                })
        isoform += 1
    return exons


# https://www.ebi.ac.uk/proteins/api/coordinates/P37840
def download_exon(acc, folder='', dest=None, save=True):
    service = "https://www.ebi.ac.uk/proteins/api/coordinates/" + acc
    down = download(service)
    if dest:
        folder = dest
    if down:
        data = json.loads(down)
        parsed = parse_exons(data)
        if save:
            open(folder + acc + '_exons.json', 'w+').write(json.dumps(parsed))
        return parsed
    else:
        if save:
            open(folder + acc + '_exons.json', 'w+').write(json.dumps([]))
        return None


def uniprot_exons(acc, folder='', dest=None):
    def format_index(i):
        return i.split(',')[0] if i != '-' else '-'

    if not dest:
        dest = folder + acc + "_exons.tsv"
    data = download_exon(acc, folder=folder, save=False)
    json_object = {}
    if data is None:
        pass
    else:
        for e in data:
            children_object = [
                {'isoform': e['isoform'], 'location': e['location'], 'start': str(e['start']), 'end': str(e['end'])}]
            json_object[f'{e["id"]}'] = children_object

    return json_object
