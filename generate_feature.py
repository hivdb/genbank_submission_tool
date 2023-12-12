from file_format import load_csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from itertools import groupby
from datetime import datetime
import re


GENE_MAP = {
    'CA': 'Capsid',
    'PR': 'protease',
    'RT': 'reverse transcriptase',
    'IN': 'integrase',
}


def format_date(collect_date):
    date = datetime.strptime(collect_date, '%Y-%m-%d')
    return datetime.strftime(date, '%d-%b-%Y')


def parse_gap(gap):
    match = re.match(
        r'\[(?P<start>\d+)\-(?P<stop>\d+)\]'
        r' \((?P<length>\d+)\)',
        gap)

    match = match.groupdict()
    return int(match['start']), int(match['stop']), int(match['length'])


def get_NA_length(items):
    seq = ''.join([
            i['aligned_NA']
            for i in items
        ])

    seq = seq.lstrip('.').rstrip('.')

    NA_length = len(seq)

    return NA_length


def feature_gene(feature_table, items):

    NA_length = get_NA_length(items)

    feature_table.append('\t'.join([
        '<1',
        f'>{NA_length}',
        'gene'
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'gene', 'pol'
    ]))


def feature_gaps(feature_table, items):
    gap_list = []

    prev_length = 0
    for i in items:
        if not i['gaps']:
            prev_length += int(i['gene_NA_length'])
            continue

        gaps = i['gaps'].split('\n')
        for gap in gaps:
            start, stop, gap_length = parse_gap(gap)
            start = start * 3 - 2 + prev_length
            stop = stop * 3 + prev_length
            gap_length = gap_length * 3

            gap_list.append((start, stop, gap_length))

        prev_length += int(i['gene_NA_length'])

    prev_length = 0
    start = 1
    stop = 0
    for i in items:
        stop = int(i['start_NA_pos']) + prev_length - 1

        if stop >= start and start != 1:
            # print(items[0]['Isolate'], start, stop)
            gap_length = stop - start + 1
            gap_list.append((start, stop, gap_length))

        start = int(i['stop_NA_pos']) + 1 + prev_length

        prev_length += int(i['gene_NA_length'])

    for start, stop, gap_length in gap_list:

        feature_table.append('\t'.join([
            f'{start}', f'{stop}', 'gap',
        ]))
        feature_table.append('\t'.join([
            '', '', '', 'estimated_length', f'{gap_length}'
        ]))


def feature_nonfunction(feature_table, items):

    NA_length = get_NA_length(items)

    has_stop = [
        i
        for i in items
        if i['stops']
    ]

    non_trans = [
        i
        for i in items
        if not int(i['translatable'])
    ]

    genes = [
        GENE_MAP[i['gene']]
        for i in items
    ]

    if has_stop or non_trans:
        feature_table.append('\t'.join([
            '<1',
            f'>{NA_length}',
            'misc_feature'
        ]))
        feature_table.append('\t'.join([
            '', '', '', 'note',
            f'nonfunctional pol protein due to mutation; contains {" and ".join(genes)}'
        ]))
        return True
    else:
        return False


def feature_cds(feature_table, items):

    NA_length = get_NA_length(items)
    genes = [
        GENE_MAP[i['gene']]
        for i in items
    ]

    feature_table.append('\t'.join([
        '<1',
        f'>{NA_length}',
        'CDS',
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'gene', 'pol'
    ]))

    start = items[0]['start_NA_pos']
    feature_table.append('\t'.join([
        '', '', '', 'codon_start', start
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'product', ' and '.join(genes)
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'transl_table', '1'
    ]))


def generate_feature_table(seq_info, feature_file):
    seq_info = load_csv(seq_info)
    seq_info.sort(key=itemgetter('Isolate'))
    seq_info = {
        isolate: list(items)
        for isolate, items in groupby(seq_info, key=itemgetter('Isolate'))
    }

    feature_table = []
    for isolate, items in seq_info.items():
        items.sort(key=lambda x: ['CA', 'PR', 'RT', 'IN'].index(x['gene']))

        feature_table.append(f'>Feature {isolate}')

        feature_gene(feature_table, items)

        feature_gaps(feature_table, items)

        nonfunction = feature_nonfunction(feature_table, items)

        if not nonfunction:
            feature_cds(feature_table, items)

        feature_table.append('')

    with open(feature_file, 'w') as f:
        for i in feature_table:
            f.write(f"{i}\n")


if __name__ == '__main__':
    seq_info = sys.argv[1]
    feature_file = sys.argv[2]
    generate_feature_table(seq_info, feature_file)
