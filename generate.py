from pathlib import Path
import pandas as pd
import numpy as np
import re
import itertools


WS = Path(__file__).resolve().parent


def work(seq, modifier):
    seq = pd.read_csv(seq)
    modifier = pd.read_csv(modifier)

    df = pd.merge(
        seq, modifier,
        on='sequence_header',
        how='left',
        indicator=True)


if __name__ == '__main__':
    seq_file = WS / 'sample' / 'aligned_sequences.csv'
    modifier_file = WS / 'sample' / 'source_modifier.csv'
    work(seq_file, modifier_file)
