from file_format import load_csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from itertools import groupby
from datetime import datetime


GENE_ORDER = ['CA', 'PR', 'RT', 'IN']


def format_date(collect_date):
    date = datetime.strptime(collect_date, '%Y-%m-%d')
    return datetime.strftime(date, '%d-%b-%Y')


def reformat_sequence(seq, isolate):
    # Trim heading and tailing .
    seq = seq.lstrip('.').rstrip('.')

    # Use N for unsequenced part
    seq = seq.replace('.', 'N')
    seq = seq.replace('-', 'N')

    assert not ('.' in seq), isolate
    assert not ('-' in seq), isolate

    return seq


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
        items.sort(key=lambda x: GENE_ORDER.index(x['gene']))

        seq = ''.join([
            i['aligned_NA']
            for i in items
        ])

        seq = reformat_sequence(seq, isolate)

        mod = modifier[isolate]

        description = {
            'Collection_date': format_date(mod['Collection_date']),
            'Country': mod['Country'],
            'Organism': organism,
            'Host': host,
            'Isolate': isolate,
            'Isolation source': mod['Isolation source'],
            'Subtype': mod['Subtype'] if 'Subtype' in mod else items[0]['subtype']
        }

        if mod['note']:
            description['note'] = mod['note']

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


if __name__ == '__main__':
    seq_info = sys.argv[1]
    source_modifier = sys.argv[2]
    organism = sys.argv[3]
    host = sys.argv[4]
    bankit_file = sys.argv[5]
    generate_bankit(seq_info, source_modifier, organism, host, bankit_file)
