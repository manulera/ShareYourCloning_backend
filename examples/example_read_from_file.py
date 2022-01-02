import requests
from example_settings import main_url

# Read a gb file with gbk extension
files = {'file': open(
    'plasmids/addgene-plasmid-39296-sequence-49545.gbk', 'rb')}

resp = requests.post(main_url + 'read_from_file', files=files)
print(resp.json())

# Read a dna file (snapgene)

files = {'file': open(
    'plasmids/addgene-plasmid-39296-sequence-49545.dna', 'rb')}

resp = requests.post(main_url + 'read_from_file', files=files)
print(resp.json())
