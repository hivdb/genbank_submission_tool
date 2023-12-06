# from hashlib import sha512
from pathlib import Path
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import sys
from file_format import dump_csv
from file_format import load_json
from itertools import groupby


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


def get_NA_seq(gene):
    # adjustedAlignedNAs has errors, they dont have insertion
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

    # Remove heading and tailing untranslatable codon
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

    return {
        'gene_AA_length': gene_AA_length,
        'gene_NA_length': gene_NA_length,
        'aligned_NA_length': aligned_NA_length,
        'start_AA_pos': start_AA_pos,
        'start_NA_pos': start_NA_pos,
        'stop_AA_pos': stop_AA_pos,
        'stop_NA_pos': stop_NA_pos,
        'N in seq': 'yes' if '.' in trimmed_seq_NA else '',
        'del in seq': 'yes' if '-' in trimmed_seq_NA else '',
        'translatable': 1 if not untrans_reason else 0,
        'untrans_reason': untrans_reason,
        'trimmed_NA': trimmed_seq_NA,
        'aligned_NA': aligned_NA,
    }


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

            record = {
                'Isolate': seq['inputSequence']['header'],
                'virus': seq['strain']['name'],
                'subtype': seq['bestMatchingSubtype'].get(
                    'displayWithoutDistance', ''),
                'gene': gene['gene']['name'],
                'mutations': ', '.join(mutations),
                'num_mutations': len(mutations),
                'insertions': ', '.join(insertions),
                'deletions': ', '.join([str(i) for i in deletions]),
                'gaps': '\n'.join([
                    f'[{start}-{stop}] ({length})'
                    for start, stop, length in gaps
                ]),
                'stops': ', '.join(stops),
            }

            record.update(get_NA_seq(gene))

            report.append(record)

    dump_csv(report_file, report)

    num_untranslate = [
        i
        for i in report
        if i['translatable'] != 1
    ]

    print(f'# untranslatable: {len(num_untranslate)}')


if __name__ == '__main__':
    aligned_file = Path(sys.argv[1]).resolve()
    seq_meta_file = Path(sys.argv[2]).resolve()

    dump_seq_meta_info(seq_meta_file, aligned_file)
