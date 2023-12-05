from pathlib import Path
from file_format import load_csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from itertools import groupby
from datetime import datetime


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

        # seq = seq.lstrip('.').rstrip('.').lstrip('-').rstrip('-')
        # seq = seq.replace('.', 'N')

        # assert not ('.' in seq), isolate
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


if __name__ == '__main__':
    seq_info = sys.argv[1]
    source_modifier = sys.argv[2]
    organism = sys.argv[3]
    host = sys.argv[4]
    bankit_file = sys.argv[5]
    feature_file = sys.argv[6]
    generate_bankit(seq_info, source_modifier, organism, host, bankit_file)
