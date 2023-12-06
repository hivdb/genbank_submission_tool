import requests
import json
# from hashlib import sha512
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import sys
import more_itertools
from tqdm import tqdm
from file_format import dump_csv
from file_format import dump_json
from file_format import load_json
from itertools import groupby


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


def get_gaps(deletions):

    gaps = []
    for _, pos_list in groupby(
            enumerate(deletions),
            lambda pair: pair[-1] - pair[0]):
        pos_list = [
            pos
            for idx, pos in pos_list
        ]
        start = pos_list[0]
        stop = pos_list[-1]
        gaps.append((start, stop, stop - start + 1))
    return gaps


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
                i['position']
                for i in gene['mutations']
                if i['isDeletion']
                ]
            stops = [
                i['shortText']
                for i in gene['mutations']
                if i['hasStop']
                ]

            gaps = get_gaps(deletions)

            # adjustedAlignedNAs dont have insertion
            # seq_NA = gene['adjustedAlignedNAs']
            # length = len(seq_NA)

            gene_AA_length = gene['gene']['length']
            gene_NA_length = gene_AA_length * 3

            AA_positions = gene['prettyPairwise']['positionLine']
            AA_positions = [
                int(i)
                for i in AA_positions
                if i.strip()
            ]

            trimmed_seq_NA = gene['prettyPairwise']['alignedNAsLine']
            if '-' in trimmed_seq_NA[0]:
                trimmed_seq_NA = trimmed_seq_NA[1:]
                AA_positions.pop(0)
            if '-' in trimmed_seq_NA[-1]:
                trimmed_seq_NA = trimmed_seq_NA[:-1]
                AA_positions.pop(-1)

            trimmed_seq_NA = ''.join(trimmed_seq_NA)

            start_AA_pos = min(AA_positions)
            start_NA_pos = start_AA_pos * 3 - 2

            stop_AA_pos = max(AA_positions)
            stop_NA_pos = stop_AA_pos * 3

            aligned_NA_length = stop_NA_pos - start_NA_pos + 1

            prepend = '.' * (start_NA_pos - 1)
            append = '.' * (gene_NA_length - stop_NA_pos)
            aligned_NA = prepend + trimmed_seq_NA + append

            untrans_reason = ''

            if (start_NA_pos % 3) != 1:
                untrans_reason = 'start NA pos not aligned to codon'
            elif aligned_NA_length % 3 != 0:
                untrans_reason = 'stop NA pos not aligned to codon'
            else:
                try:
                    Seq(trimmed_seq_NA).translate()
                except TranslationError:
                    untrans_reason = 'untranslatable'

            report.append({
                'Isolate': seq['inputSequence']['header'],
                'virus': seq['strain']['name'],
                'subtype': seq['bestMatchingSubtype'].get(
                    'displayWithoutDistance', ''),
                'gene': gene['gene']['name'],
                'gene_AA_length': gene_AA_length,
                'gene_NA_length': gene_NA_length,
                'NA_length': aligned_NA_length,
                'start_NA_pos': start_NA_pos,
                'stop_NA_pos': stop_NA_pos,
                # 'firstNA': gene['firstNA'],
                # 'lastNA': gene['lastNA'],
                'mutations': ', '.join(mutations),
                'num_mutations': len(mutations),
                'insertions': ', '.join(insertions),
                'deletions': ', '.join([str(i) for i in deletions]),
                'gaps': '\n'.join([
                    f'[{start}-{stop}] ({length})'
                    for start, stop, length in gaps
                ]),
                'hasStop': ', '.join(stops),
                'N in seq': 'yes' if '.' in trimmed_seq_NA else '',
                'del in seq': 'yes' if '-' in trimmed_seq_NA else '',
                'translatable': 1 if not untrans_reason else 0,
                'untrans_reason': untrans_reason,
                'trimmed_NA': trimmed_seq_NA,
                'aligned_NA': aligned_NA,
            })

    dump_csv(report_file, report)

    num_untranslate = [
        i
        for i in report
        if i['translatable'] != 1
    ]

    print(f'# untranslatable: {len(num_untranslate)}')


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
