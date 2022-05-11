# %%
from bs4 import BeautifulSoup
with open('final_plasmid_history.svg') as input_handle:
    soup = BeautifulSoup(input_handle, 'xml')


with open('trimmed_svg.svg','w') as out_handle:
    for text in soup.find_all('text'):
        out_handle.write(str(text)+'\n')
