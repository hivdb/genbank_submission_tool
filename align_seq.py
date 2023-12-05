import requests
import json
# from hashlib import sha512
from pathlib import Path
from Bio import SeqIO
import sys
import more_itertools
from tqdm import tqdm
from file_format import dump_csv
from file_format import dump_json
from file_format import load_json


WS = Path(__file__).resolve().parent


def align_seq(seq_list):
    # TODO use the original seq id

    seq_list = [
        {
            # 'header': sha512(seq.encode('utf-8')).hexdigest(),
            'header': header,
            'sequence': seq,
        }
        for header, seq in seq_list.items()
    ]

    aligned = []

    with tqdm(total=len(seq_list)) as pbar:
        for seq in more_itertools.batched(seq_list, 20):
            aligned.extend(
                query_sierra(seq)[
                    'data']['viewer']['sequenceAnalysis'])
            pbar.update(20)

    return aligned


def query_sierra(sequences):
    url = 'https://hivdb.stanford.edu/graphql'

    graphql = open(WS / 'sierra.graphql').read()

    resp = requests.post(
        url,
        data=json.dumps({
            'operationName': 'align',
            'query': graphql,
            'variables': {
                'sequences': sequences
            }
        })
    )

    return resp.json()


def dump_seq_meta_info(report_file, aligned_file):

    aligned_seq = load_json(aligned_file)

    report = []
    for seq in aligned_seq:
        for gene in seq['alignedGeneSequences']:
            mutations = [
                i['shortText']
                for i in gene['mutations']]
            insertions = [
                i['shortText']
                for i in gene['mutations']
                if i['isInsertion']
                ]
            deletions = [
                i['shortText']
                for i in gene['mutations']
                if i['isDeletion']
                ]

            report.append({
                'Isolate': seq['inputSequence']['header'],
                'virus': seq['strain']['name'],
                'subtype': seq['bestMatchingSubtype'].get(
                    'displayWithoutDistance', ''),
                'gene': gene['gene']['name'],
                'length': gene['gene']['length'],
                'firstNA': gene['firstNA'],
                'lastNA': gene['lastNA'],
                'mutations': ', '.join(mutations),
                'num_mutations': len(mutations),
                'insertions': ', '.join(insertions),
                'deletions': ', '.join(deletions),
                'aligned_NA': gene['adjustedAlignedNAs']
            })

    dump_csv(report_file, report)


if __name__ == '__main__':
    unaligned_file = Path(sys.argv[1]).resolve()
    aligned_file = Path(sys.argv[2]).resolve()
    seq_meta_file = Path(sys.argv[3]).resolve()

    if aligned_file.exists():
        print("Aligned sequence is already there.")
        dump_seq_meta_info(seq_meta_file, aligned_file)
        exit()

    seq_list = {}
    for i in SeqIO.parse(unaligned_file, 'fasta'):
        seq_list[i.id] = str(i.seq)

    aligned = align_seq(seq_list)

    dump_json(aligned_file, aligned)

    dump_seq_meta_info(seq_meta_file, aligned_file)
