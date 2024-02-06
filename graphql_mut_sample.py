import requests
import json
from pathlib import Path

WS = Path(__file__).resolve().parent


def query_sierra():
    url = 'https://hivdb.stanford.edu/graphql'

    graphql = open('./sierra.mutation.graphql').read()

    resp = requests.post(
        url,
        data=json.dumps({
            'operationName': 'align',
            'query': graphql,
            'variables': {
                'algorithms': [],
                'customAlgoriths': [],
                'includeGenes': ['PR', 'RT', 'IN'],
                'patternNames': ["RT:V108I+RT:P225H"],
                'patterns': [["RT:V108I", "RT:P225H"]],
            }
        })
    )

    return resp.json()
