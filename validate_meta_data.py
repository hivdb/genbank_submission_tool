from file_format import load_csv
import sys
from itertools import groupby


def validate(source_file):
    meta_data = load_csv(source_file)

    isolate = [
        i['Isolate']
        for i in meta_data
    ]

    isolate.sort()

    for key, items in groupby(isolate):
        items = list(items)
        if not key:
            print(f'Empty isolate: {key}')
        if len(items) > 1:
            print(f"{key} is duplidated isolate")


if __name__ == '__main__':
    validate(sys.argv[1])
