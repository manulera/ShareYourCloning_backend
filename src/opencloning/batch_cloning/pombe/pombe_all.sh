set -e
# Clear the output directory except for the gene_list.txt
find batch_cloning_output -type f ! -name 'gene_list.txt' -delete
find batch_cloning_output -type d -empty -delete

python pombe_get_primers.py --genes batch_cloning_output/gene_list.txt
python pombe_clone.py --genes batch_cloning_output/gene_list.txt
python pombe_summary.py
python pombe_gather.py
