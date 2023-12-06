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


def generate_bankit(seq_info, modifier, organism, host, bankit_file):
    seq_info = load_csv(seq_info)
    seq_info.sort(key=itemgetter('Isolate'))
    seq_info = {
        isolate: list(items)
        for isolate, items in groupby(seq_info, key=itemgetter('Isolate'))
    }
    modifier = load_csv(modifier)
    modifier = {
        i['Isolate']: i
        for i in modifier
    }

    fasta_seq = []

    for isolate, items in seq_info.items():
        items.sort(key=lambda x: ['CA', 'PR', 'RT', 'IN'].index(x['gene']))

        seq = ''.join([
            i['aligned_NA']
            for i in items
        ])

        seq = seq.lstrip('.').rstrip('.')
        seq = seq.replace('.', 'N')

        assert not ('.' in seq), isolate
        # assert not ('-' in seq), isolate

        mod = modifier[isolate]

        description = {
            'Collection_date': format_date(mod['Collection_date']),
            'Country': mod['Country'],
            'Organism': organism,
            'Host': host,
            'Isolate': isolate,
            'Isolate source': mod['Isolation source']
        }

        description = ' '.join([
            f"[{k}={v}]"
            for k, v in description.items()
        ])

        fasta_seq.append(
            SeqRecord(
                id=isolate,
                seq=Seq(seq),
                description=description
            )
        )

    SeqIO.write(
        fasta_seq,
        str(bankit_file),
        'fasta')


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

        NA_length = sum([
            int(i['gene_NA_length'])
            for i in items
        ])

        feature_table.append('\t'.join([
            '<1',
            f'>{NA_length}',
            'gene'
        ]))

        genes = [
            GENE_MAP[i['gene']]
            for i in items
        ]

        feature_table.append('\t'.join([
            '', '', '', 'gene', 'pol'
        ]))

        # TODO gaps
        prev_length = 0
        for i in items:
            gaps = i['gaps'].split('\n')
            for gap in gaps:
                if not gap:
                    continue
                start, stop, gap_length = parse_gap(gap)
                start = start * 3 - 2 + prev_length
                stop = stop * 3 + prev_length
                gap_length = gap_length * 3
                feature_table.append('\t'.join([
                    f'{start}', f'{stop}', 'gap',
                ]))
                feature_table.append('\t'.join([
                    '', '', '', 'estimated_length', f'{gap_length}'
                ]))
            prev_length += int(i['gene_NA_length'])

        # stop codon
        has_stop = [
            i
            for i in items
            if i['hasStop']
        ]
        # TODO if the na insertion is at the end, remove it.
        non_trans = [
            i
            for i in items
            if not int(i['translatable'])
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

        feature_table.append('')

    with open(feature_file, 'w') as f:
        for i in feature_table:
            f.write(f"{i}\n")


def parse_gap(gap):
    match = re.match(
        r'\[(?P<start>\d+)\-(?P<stop>\d+)\]'
        r' \((?P<length>\d+)\)',
        gap)

    match = match.groupdict()
    return int(match['start']), int(match['stop']), int(match['length'])


if __name__ == '__main__':
    seq_info = sys.argv[1]
    source_modifier = sys.argv[2]
    organism = sys.argv[3]
    host = sys.argv[4]
    bankit_file = sys.argv[5]
    feature_file = sys.argv[6]
    generate_bankit(seq_info, source_modifier, organism, host, bankit_file)
    generate_feature_table(seq_info, feature_file)
