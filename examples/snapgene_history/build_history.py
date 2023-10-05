"""
A script to clone a gene into a plasmid. The same operation will be done in Snapgene,
and we will try to extract the snapgene history to match the documentation by ShareYourCloning.
"""
# %%
import sys
sys.path.append('../../')

from pydantic_models import RepositoryIdSource, PrimerAnnealingSettings, PrimerModel, PCRSource, RestrictionEnzymeDigestionSource  # noqa
from fastapi.testclient import TestClient  # noqa
from main import app  # noqa
import json  # noqa
# %%

# lists to collect sources and sequences
sources = list()
sequences = list()

# Load the app client, where we can make mock http requests

client = TestClient(app)

# 1. Request Ase1 to genbank (a pombe gene) ==========================

input_source = RepositoryIdSource(
    repository='genbank',
    repository_id='NM_001018957.2',
)

payload = client.post('/repository_id', json=input_source.model_dump()).json()

# There should only be one in the list, so we can simply sum them
source = payload['sources'][0]
source['id'] = 1
source['output'] = 2
sequence = payload['sequences'][0]
sequence['id'] = 2
sources.append(source)
sequences.append(sequence)

# 2. Amplify Ase1 by PCR ==========================

primer_fwd = PrimerModel(
    sequence='TAAGCAGTCGACatgcaaacagtaatgatggatg',
    id=3,
    name='ase1_forward'
)

primer_rvs = PrimerModel(
    sequence='TAAGCAGGCGCGCCaaagccttcttctccccatt',
    id=4,
    name='ase1_rvs'
)

input_source = PCRSource(
    input=[2],
    primer_annealing_settings=PrimerAnnealingSettings(minimum_annealing=13)
)

primers = [primer_fwd.model_dump(), primer_rvs.model_dump()]

data = {'source': input_source.model_dump(), 'sequences': sequences, 'primers': primers}
payload = client.post("/pcr", json=data).json()

source = payload['sources'][0]
source['id'] = 5
source['output'] = 6
sequence = payload['sequences'][0]
sequence['id'] = 6
sources.append(source)
sequences.append(sequence)

# 3. Digest PCR product ==========================

input_source = RestrictionEnzymeDigestionSource(
    input=[6],
    restriction_enzymes=['AscI', 'SalI'],
)

data = {'source': input_source.model_dump(), 'sequences': [sequences[-1]]}
payload = client.post("/restriction", json=data).json()

source = payload['sources'][1]
source['id'] = 7
source['output'] = 8
sequence = payload['sequences'][1]
sequence['id'] = 8
sources.append(source)
sequences.append(sequence)

# 4. Load plasmid from file ==========================

with open('addgene-plasmid-39296-sequence-49545.dna', 'rb') as f:
    payload = client.post("/read_from_file", files={"file": f}).json()

source = payload['sources'][0]
source['id'] = 9
source['output'] = 10
sequence = payload['sequences'][0]
sequence['id'] = 10
sources.append(source)
sequences.append(sequence)


# 5. Digest the plasmid ==========================

input_source = RestrictionEnzymeDigestionSource(
    input=[10],
    restriction_enzymes=['AscI', 'SalI'],
)

data = {'source': input_source.model_dump(), 'sequences': [sequences[-1]]}
payload = client.post("/restriction", json=data).json()

source = payload['sources'][1]
source['id'] = 11
source['output'] = 12
sequence = payload['sequences'][1]
sequence['id'] = 12
sources.append(source)
sequences.append(sequence)

output_dict = {
    'sequences': sequences,
    'sources': sources,
    'primers': primers
}


with open('history.json', 'w') as fp:
    json.dump(output_dict, fp)
