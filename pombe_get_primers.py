from bs4 import BeautifulSoup
import asyncio
from httpx import AsyncClient, Response
import re

default_settings_primer_design = {
    'length': 80,
    'increment': 40,
    'add_seq': 400,
    'opt_len': 22,
    'min_len': 20,
    'max_len': 28,
    'opt_tm': 60.0,
    'min_tm': 57.0,
    'max_tm': 63.0,
    'min_gc': 30,
    'max_gc': 60,
    '.submit': 'Submit',
}

genes = ['ase1']


async def async_post(url, headers, data, params=None) -> Response:
    async with AsyncClient(timeout=20.0) as client:
        return await client.post(url, headers=headers, data=data, params=params)


async def main(gene):

    # A first request to access the primers
    data = {
        'gene': gene,
        'length': 80,
        'plasmid': 'pFA6a',
        'increment': 40,
        '.submit': 'Submit',
        '.cgifields': 'plasmid',
    }

    url = 'http://bahlerweb.cs.ucl.ac.uk/cgi-bin/PPPP/pppp_deletion.pl'
    headers = None
    resp = await async_post(url, headers, data)
    # Parse with BeautifulSoup
    soup = BeautifulSoup(resp.text, 'html.parser')
    # select forward and reverse primers by default
    forward_primer = soup.find('input', {'name': 'for_sel', 'checked': True})
    reverse_primer = soup.find('input', {'name': 'rev_sel', 'checked': True})
    forward_primer_seq = re.sub(r'[^a-zA-Z]', '', forward_primer['value'])
    reverse_primer_seq = re.sub(r'[^a-zA-Z]', '', reverse_primer['value'])
    # Make a second request to get the checking primers
    data = {
        'gene': gene,
        'for_sel': forward_primer['value'],
        'rev_sel': reverse_primer['value'],
    }
    data.update(default_settings_primer_design)

    url = 'http://bahlerweb.cs.ucl.ac.uk/cgi-bin/PPPP/pppp_checking.pl'
    resp = await async_post(url, headers, data)
    # Parse with BeautifulSoup
    soup = BeautifulSoup(resp.text, 'html.parser')
    # Find a pre tag with the text "Left Primer:"
    left_check_primer = soup.find('pre', string=re.compile(r'Left Primer:\s+Sequence:')).get_text().strip()
    right_check_primer = soup.find('pre', string=re.compile(r'Right Primer:\s+Sequence:')).get_text().strip()

    pattern = r'Sequence:\s+(\S+)'
    left_check_primer_seq = re.search(pattern, left_check_primer).group(1)
    right_check_primer_seq = re.search(pattern, right_check_primer).group(1)
    print(forward_primer_seq)
    print(reverse_primer_seq)
    print(left_check_primer_seq)
    print(right_check_primer_seq)


if __name__ == '__main__':
    asyncio.run(main('ase1'))
