#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Assembly of sequences by homologous recombination.

Should also be useful for related techniques such as Gibson assembly and fusion
PCR. Given a list of sequences (Dseqrecords), all sequences are analyzed for
shared homology longer than the set limit.

A graph is constructed where each overlapping region form a node and
sequences separating the overlapping regions form edges.


::


                -- A --
    catgatctacgtatcgtgt     -- B --
                atcgtgtactgtcatattc
                            catattcaaagttct



    --x--> A --y--> B --z-->   (Graph)

    Nodes:

    A : atcgtgt
    B : catattc

    Edges:

    x : catgatctacgt
    y : actgt
    z : aaagttct


The NetworkX package is used to trace linear and circular paths through the
graph.
"""
from Bio.SeqFeature import ExactPosition as _ExactPosition
from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from pydna.utils import rc as _rc
from pydna.utils import memorize as _memorize
from pydna._pretty import pretty_str as _pretty_str
from pydna.contig import Contig as _Contig
from pydna.common_sub_strings import common_sub_strings
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
import networkx as _nx
from copy import deepcopy as _deepcopy
import itertools as _itertools
import logging as _logging
from Bio.Seq import reverse_complement

_module_logger = _logging.getLogger("pydna." + __name__)


# TODO use quicker inits for contig
# TODO remove maxnodes for init

def match_to_location_pair(match, second_length=None) -> tuple[_SimpleLocation, _SimpleLocation]:
    x_start, y_start, length = match
    x_loc = _SimpleLocation(x_start, x_start + length, 1)
    if second_length is not None:
        y_loc = _SimpleLocation(y_start, y_start + length, -1)._flip(second_length)
    else:
        y_loc = _SimpleLocation(y_start, y_start + length, 1)

    return x_loc, y_loc


class _Memoize(type):
    @_memorize("pydna.assembly.Assembly")
    def __call__(cls, *args, **kwargs):
        return super().__call__(*args, **kwargs)


class Assembly(object, metaclass=_Memoize):
    """Assembly of a list of linear DNA fragments into linear or circular
    constructs. The Assembly is meant to replace the Assembly method as it
    is easier to use. Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object. Several methods are available for analysis
    of overlapping sequences, graph construction and assembly.

    Parameters
    ----------

    fragments : list
        a list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered
    algorithm : function, optional
        The algorithm used to determine the shared sequences.
    max_nodes : int
        The maximum number of nodes in the graph. This can be tweaked to
        manage sequences with a high number of shared sub sequences.



    Examples
    --------

    >>> from pydna.assembly import Assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    Assembly
    fragments....: 33bp 34bp 35bp
    limit(bp)....: 14
    G.nodes......: 6
    algorithm....: common_sub_strings
    >>> x.assemble_circular()
    [Contig(o59), Contig(o59)]
    >>> x.assemble_circular()[0].seq.watson
    'acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGt'


    """

    def __init__(self, frags=None, limit=25, algorithm=common_sub_strings, use_fragment_order=True):
        # TODO: allow for the same fragment to be included more than once
        G = _nx.MultiDiGraph()
        G.add_nodes_from((i + 1, {'seq': f}) for (i, f) in enumerate(frags))
        G.add_nodes_from((-(i + 1), {'seq': f.reverse_complement()}) for (i, f) in enumerate(frags))
        # all combinations of fragments are compared.
        # see https://docs.python.org/3.10/library/itertools.html
        # itertools.combinations('ABCD', 2)-->  AB AC AD BC BD CD
        fragment_pairs = list(_itertools.combinations(enumerate(frags), 2))
        if use_fragment_order:
            temp = list(enumerate(frags))
            fragment_pairs = list(zip(temp, temp[1:] + temp[:1]))

        for (index_first, first), (index_secnd, secnd) in fragment_pairs:
            index_first += 1
            index_secnd += 1
            matches_fwd = algorithm(str(first.seq).upper(), str(secnd.seq).upper(), limit)
            for pair in map(match_to_location_pair, matches_fwd):
                G.add_edge(index_first, index_secnd, f'{index_first}{pair[0]}:{index_secnd}{pair[1]}', locations=pair)
                G.add_edge(-index_first, -index_secnd, f'{index_first}_rc{pair[0]}:{index_secnd}_rc{pair[1]}', locations=[pair[1]._flip(len(secnd)), pair[0]._flip(len(first))])

            matches_rvs = algorithm(str(first.seq).upper(), reverse_complement(str(secnd.seq).upper()), limit)
            for pair in map(match_to_location_pair, matches_rvs):
                G.add_edge(index_first, -index_secnd, f'{index_first}{pair[0]}:{index_secnd}_rc{pair[1]}', locations=pair)
                G.add_edge(-index_first, index_secnd, f'{index_first}_rc{pair[0]}:{index_secnd}{pair[1]}', locations=[pair[1]._flip(len(secnd)), pair[0]._flip(len(first))])

        self.G = G
        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order

        return

    def validate_assembly(self, assembly):
        # The assembly must include all fragments (this could be made optional)
        if len(assembly) != len(self.fragments) - 1:
            return False
        return True
        # Here we check whether subsequent pairs of fragments are compatible, for instance:
        # Compatible (overlap of 1 and 2 occurs before overlap of 2 and 3):
        #    -- A --
        #  gtatcgtgt     -- B --
        #    atcgtgtactgtcatattc
        #                catattcaa
        # Incompatible (overlap of 1 and 2 occurs after overlap of 2 and 3):
        #               -- A --
        #  -- B --    gtatcgtgt
        #  catattcccccccatcgtgtactgt
        #  catattcaa
        fragment_pairs = zip(assembly, assembly[1:])
        for (u1, v1, key1), (u2, v2, key2) in fragment_pairs:
            end_of_1 = self.G[u1][v1][key1]['locations'][1].parts[-1].end
            start_of_2 = self.G[u2][v2][key2]['locations'][0].parts[0].end
            # TODO: double-check exact match
            if end_of_1 > start_of_2:
                return False

        return True

    def get_linear_assemblies(self):
        # The constrain for linear assemblies is that the first node is in the initial orientation

        # Copy the graph since we will add the begin and end nodes
        G = _nx.MultiDiGraph(self.G)
        G.add_nodes_from(['begin', 'end'])
        if self.use_fragment_order:
            G.add_edge('begin', 1, locations=(None, None))
        for node in G.nodes:
            if not self.use_fragment_order:
                if type(node) == int and node > 0:
                    G.add_edge('begin', node, locations=(None, None))
            G.add_edge(node, 'end', locations=(None, None))

        return list(filter(self.validate_assembly, map(lambda x : x[1:-1], _nx.all_simple_edge_paths(G, 'begin', 'end'))))

    def expand_circular_assembly(self, pairs):
        # Feels like this should be built-in to networkx, or that find_cycle would return all possible cycles
        combine = list()
        for u, v, _, _ in pairs:
            combine.append([key for key in self.G[u][v]])
        return list(_itertools.product(*combine))

    def get_circular_assemblies(self):
        # The constrain of circular sequence is that the first node is the first fragment in its initial orientation
        return self.expand_circular_assembly(_nx.cycles.find_cycle(self.G, source=1, orientation='original'))

    def execute_assembly(self, assembly):

        pos = 0
        first_fragment = assembly[0][0]
        out_dseq = self.G.nodes[first_fragment]['seq']
        for left, right, key in assembly:
            if right == first_fragment:
                out_dseq = out_dseq.looped()
                break
            left_location, right_location = self.G[left][right][key]['locations']
            right_seq = self.G.nodes[right]['seq']

            new_dseq = out_dseq.seq[:left_location.parts[-1].end + pos] + right_seq.seq[right_location.parts[-1].end:]

            overlap_length = len(left_location)

            # This could be a separate Dseqrecord method
            new_features = right_seq[right_location.parts[0].start:].features
            for feature in new_features:
                feature = feature._shift(left_location.parts[-1].end - overlap_length)
            out_dseq = _Dseqrecord(new_dseq, features=out_dseq.features + new_features)

            pos += left_location.parts[-1].end

        return out_dseq







    def __repr__(self):
        # https://pyformat.info
        return _pretty_str(
            "Assembly\n"
            "fragments..: {sequences}\n"
            "limit(bp)..: {limit}\n"
            "G.nodes....: {nodes}\n"
            "algorithm..: {al}".format(
                sequences=" ".join(
                    "{}bp".format(len(x["mixed"])) for x in self.fragments
                ),
                limit=self.limit,
                nodes=self.G.order(),
                al=self.algorithm.__name__,
            )
        )


example_fragments = (
    _Dseqrecord("AacgatCAtgctcc", name="a"),
    _Dseqrecord("TtgctccTAAattctgc", name="b"),
    _Dseqrecord("CattctgcGAGGacgatG", name="c"),
)


linear_results = (
    _Dseqrecord("AacgatCAtgctccTAAattctgcGAGGacgatG", name="abc"),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatG", name="ac_rc"),
    _Dseqrecord("AacgatG", name="ac"),
)


circular_results = (
    _Dseqrecord("acgatCAtgctccTAAattctgcGAGG", name="abc", circular=True),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatTTA", name="abc_rc", circular=True),
)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
