from file_format import load_csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
from itertools import groupby
from datetime import datetime
from pathlib import Path


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


def generate_bankit(
        seq_info, modifier, organism, host, bankit_file, batch_size):
    seq_info = load_csv(seq_info)
    seq_info.sort(key=itemgetter('Isolate'))
    seq_info = {
        isolate: list(items)
        for isolate, items in groupby(seq_info, key=itemgetter('Isolate'))
    }

    print('# isolates:', len(seq_info))

    if not batch_size:
        generate_bankit_per_batch(
            seq_info, modifier, organism, host, bankit_file)
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
            bankit_file_name = (
                bankit_file.parent / f"{idx+1}_{bankit_file.name}")
            generate_bankit_per_batch(
                seq_info, modifier, organism, host, bankit_file_name)


def generate_bankit_per_batch(
        seq_info, modifier, organism, host, bankit_file):

    modifier = load_csv(modifier)
    modifier = {
        i['Isolate']: i
        for i in modifier
    }

    fasta_seq = []

    for isolate, items in seq_info.items():
        items.sort(key=lambda x: GENE_ORDER.index(x['gene']))

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

        if ca_items:
            fasta_seq.append(
                get_seq_record(
                    ca_items, isolate, modifier, organism, host, 'CAPSID')
            )

        if pol_items:
            fasta_seq.append(
                get_seq_record(
                    pol_items, isolate, modifier, organism, host, 'POL')
            )

    SeqIO.write(
        fasta_seq,
        str(bankit_file),
        'fasta')


def get_seq_record(items, isolate, modifier, organism, host, gene_name):

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

    return SeqRecord(
                id=isolate + '_' + gene_name,
                seq=Seq(seq),
                description=description
            )


if __name__ == '__main__':
    seq_info = sys.argv[1]
    source_modifier = sys.argv[2]
    organism = sys.argv[3]
    host = sys.argv[4]
    bankit_file = Path(sys.argv[5]).resolve()

    if len(sys.argv) == 7:
        batch_size = int(sys.argv[6])
    else:
        batch_size = None
    generate_bankit(
        seq_info, source_modifier, organism, host, bankit_file, batch_size)
