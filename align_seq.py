import requests
import json
# from hashlib import sha512
from pathlib import Path
from Bio import SeqIO
import sys
import more_itertools
from tqdm import tqdm
from file_format import dump_json


WS = Path(__file__).resolve().parent


def align_seq(seq_list):

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
    url = 'https://hivdb.stanford.edu/hcv/graphql'

    graphql = open('./sierra.graphql').read()

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


if __name__ == '__main__':
    unaligned_file = Path(sys.argv[1]).resolve()
    aligned_file = Path(sys.argv[2]).resolve()

    if aligned_file.exists():
        print("Aligned sequence is already there.")
        exit()

    seq_list = {}
    for i in SeqIO.parse(unaligned_file, 'fasta'):
        seq_list[i.id] = str(i.seq)

    aligned = align_seq(seq_list)

    dump_json(aligned_file, aligned)
