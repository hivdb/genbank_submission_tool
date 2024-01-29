from file_format import load_csv
from file_format import dump_csv
import sys
from pathlib import Path


def generate_sc(
        modifier, treatment, sc_file):
    modifier = load_csv(modifier)
    patients = {
        i['PatientID'] if i.get('PatientID') else i['Isolate']: i['Isolate']
        for i in modifier
    }

    if not treatment.exists():
        return

    treatment = load_csv(treatment)

    trigger_error = False
    for i in treatment:
        if i['PatientID'] not in patients:
            print(f"Error: Patient not found {i['PatientID']}")
            trigger_error = True

    if trigger_error:
        return

    for i in treatment:
        i.update({
            'StructuredCommentPrefix': 'HIV-DataBaseData',
            'StructuredCommentSuffix': 'HIV-DataBaseData',
            'Isolate': patients[i['PatientID']],
        })

        notes = ''
        if not i['weeks']:
            notes += 'duration unknown;'
        if not i['treatment_order']:
            notes += 'treatment order unknown;'

        i['notes'] = notes

    dump_csv(sc_file, treatment, headers=[
        'Isolate',
        'PatientID',
        'treatment',
        'weeks',
        'treatment_order',
        'notes',
        'StructuredCommentPrefix',
        'StructuredCommentSuffix',
    ], remain=False)


if __name__ == '__main__':
    modifier = Path(sys.argv[1]).resolve()
    treatment_file = Path(sys.argv[2]).resolve()
    sc_file = Path(sys.argv[3]).resolve()

    generate_sc(
        modifier, treatment_file, sc_file)
