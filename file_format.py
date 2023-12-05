import json
from pathlib import Path
import csv


def load_json(json_path):
    with json_path.open() as fd:
        return json.load(fd)


def dump_json(json_path, obj):
    json_path = Path(json_path)
    json_path.parent.mkdir(exist_ok=True, parents=True)
    with open(json_path, 'w') as fd:
        json.dump(obj, fd, indent=4, ensure_ascii=False)


def dump_csv(file_path, table, headers=[], remain=True):

    table_headers = []
    for rec in table:
        for key in rec.keys():
            if key not in table_headers:
                table_headers.append(key)

    if not headers:
        headers = table_headers
    else:
        remain_headers = [
            i
            for i in table_headers
            if i not in headers
        ]
        if remain:
            headers = headers + remain_headers
        table = [
            {
                k: v
                for k, v in i.items()
                if k in headers
            }
            for i in table
        ]

    file_path.parent.mkdir(exist_ok=True, parents=True)

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(table)


def load_csv(file_path):

    table = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            table.append(record)

    return table
