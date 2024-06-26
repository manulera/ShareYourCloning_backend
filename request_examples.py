genome_region_examples = {
    'full': {
        'summary': 'All parameters provided',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'assembly_accession': 'GCF_000002945.1',
            'locus_tag': 'SPAPB1A10.09',
            'gene_id': 2543372,
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'id_omitted': {
        'summary': 'Gene ID omitted (filled in response)',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'assembly_accession': 'GCF_000002945.1',
            'locus_tag': 'SPAPB1A10.09',
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'assembly_accession_omitted': {
        'summary': 'Assembly accession omitted (filled in response)',
        'value': {
            'id': 1,
            'sequence_accession': 'NC_003424.3',
            'start': 1877009,
            'end': 1881726,
            'strand': 1,
        },
    },
    'viral_sequence': {
        'summary': 'Viral sequence not associated with assembly',
        'value': {
            'id': 1,
            'sequence_accession': 'DQ208311.2',
            'start': 20,
            'end': 2050,
            'strand': -1,
        },
    },
}

oligonucleotide_hybridization_examples = {
    'default': {
        'summary': 'Typical example',
        'description': 'blah',
        'value': {
            'source': {
                'id': 1,
                'input': [],
                'output': 0,
                'forward_oligo': 2,
                'reverse_oligo': 3,
            },
            'primers': [
                {'id': 2, 'name': 'primer1', 'sequence': 'aaGCGGCCGCgtagaactttatgtgcttccttacattggt'},
                {'id': 3, 'name': 'primer2', 'sequence': 'aaGCGGCCGCaccaatgtaaggaagcacataaagttctac'},
            ],
        },
    },
}
