from file_format import load_csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from itertools import groupby
from datetime import datetime
import re
from pathlib import Path


GENE_MAP = {
    # 'CA': 'capsid',
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


def feature_gene(feature_table, items, gene_name):

    NA_length = get_NA_length(items)

    feature_table.append('\t'.join([
        '<1',
        f'>{NA_length}',
        'gene'
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'gene', gene_name.lower()
    ]))


def feature_gaps(feature_table, items):
    gap_list = []

    num_heading_del = int(items[0]['aligned_NA_heading'])
    # num_tailing_del = int(items[-1]['aligned_NA_tailing'])

    # Gap with in gene
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

    # gap between gene
    seq_length_acc = 0
    for i in range(len(items)):
        j = i + 1
        if j == len(items):
            break

        g1 = items[i]
        g2 = items[j]
        gap_length = int(
            g1['aligned_NA_tailing']) + int(
            g2['aligned_NA_heading'])

        if not gap_length:
            pass
        else:
            start = seq_length_acc + int(
                g1['aligned_NA_tailing_pos']) + 1
            stop = seq_length_acc + int(
                g1['gene_NA_length']) + int(
                g2['aligned_NA_heading_pos']) - 1
            gap_list.append((start, stop, gap_length))

        seq_length_acc += int(g1['gene_NA_length'])

    gap_list.sort(key=lambda x: x[0])

    for start, stop, gap_length in gap_list:
        start = start - num_heading_del
        stop = stop - num_heading_del

        feature_table.append('\t'.join([
            f'{start}', f'{stop}', 'gap',
        ]))
        feature_table.append('\t'.join([
            '', '', '', 'estimated_length', f'{gap_length}'
        ]))


def feature_nonfunction(feature_table, items, gene_name):

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
            f'nonfunctional {gene_name.lower()} protein due to mutation; contains {" and ".join(genes)}'
        ]))
        return True
    else:
        return False


def feature_cds(feature_table, items, gene_name):

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
        '', '', '', 'gene', gene_name.lower()
    ]))

    num_heading_del = int(items[0]['aligned_NA_heading'])
    start = int(items[0]['aligned_NA_heading_pos']) - num_heading_del
    feature_table.append('\t'.join([
        '', '', '', 'codon_start', str(start)
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'product', ' and '.join(genes)
    ]))

    feature_table.append('\t'.join([
        '', '', '', 'transl_table', '1'
    ]))


def generate_feature_table(seq_info, feature_file, batch_size):

    seq_info = load_csv(seq_info)
    seq_info.sort(key=itemgetter('Isolate'))
    seq_info = {
        isolate: list(items)
        for isolate, items in groupby(seq_info, key=itemgetter('Isolate'))
    }
    if not batch_size:
        generate_feature_table_per_batch(seq_info, feature_file)
    else:
        seq_info = [
            (k, v)
            for k, v in seq_info.items()
        ]
        batches = [
            dict(seq_info[i:i + batch_size])
            for i in range(0, len(seq_info), batch_size)
        ]
        for idx, seq_info in enumerate(batches):
            feature_file_name = (
                feature_file.parent / f"{idx+1}_{feature_file.name}")
            generate_feature_table_per_batch(seq_info, feature_file_name)


def generate_feature_table_per_batch(seq_info, feature_file):

    feature_table = []
    for isolate, items in seq_info.items():
        items.sort(key=lambda x: ['CA', 'PR', 'RT', 'IN'].index(x['gene']))

        ca_items = [
            i
            for i in items
            if i['gene'] == 'CA'
        ]

        pol_items = [
            i
            for i in items
            if i['gene'] != 'CA'
        ]

        # if ca_items:
        #     update_feature_table(feature_table, isolate, ca_items, 'CAPSID')

        if pol_items:
            update_feature_table(feature_table, isolate, pol_items, 'POL')

    with open(feature_file, 'w') as f:
        for i in feature_table:
            f.write(f"{i}\n")


def update_feature_table(feature_table, isolate, items, gene_name):
    feature_table.append(f'>Feature {isolate}_{gene_name}')

    feature_gene(feature_table, items, gene_name)

    feature_gaps(feature_table, items)

    nonfunction = feature_nonfunction(feature_table, items, gene_name)

    if not nonfunction:
        feature_cds(feature_table, items, gene_name)

    feature_table.append('')


if __name__ == '__main__':
    seq_info = sys.argv[1]
    feature_file = Path(sys.argv[2])

    if len(sys.argv) == 4:
        batch_size = int(sys.argv[3])
    else:
        batch_size = None
    generate_feature_table(seq_info, feature_file, batch_size)
