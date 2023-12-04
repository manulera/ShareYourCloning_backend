"""Slightly different assembly implementation"""

from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from pydna.utils import shift_location as _shift_location
from pydna.utils import memorize as _memorize
from pydna._pretty import pretty_str as _pretty_str
from pydna.common_sub_strings import common_sub_strings
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
import networkx as _nx
from copy import deepcopy as _deepcopy
import itertools as _itertools
import logging as _logging
from Bio.Seq import reverse_complement
from dna_functions import create_location

def circular_permutation_min_abs(lst):
    """Returns the circular permutation of lst with the smallest absolute value first."""
    min_abs_index = min(range(len(lst)), key=lambda i: abs(lst[i]))
    return lst[min_abs_index:] + lst[:min_abs_index]

def add_edges_from_match(match, index_first, index_secnd, first, secnd, graph: _nx.MultiDiGraph):
    x_start, y_start, length = match
    # We use this function in case there are origin-spanning features. Circular case should be handled properly.
    locs = [create_location(x_start, x_start + length, 1, len(first), first.circular),
            create_location(y_start, y_start + length, 1, len(secnd), secnd.circular)]
    rc_locs = [locs[0]._flip(len(first)), locs[1]._flip(len(secnd))]

    combinations = (
        (index_first, index_secnd, locs),
        (index_secnd, index_first, locs[::-1]),
        (-index_first, -index_secnd, rc_locs),
        (-index_secnd, -index_first, rc_locs[::-1]),
    )
    for u, v, l in combinations:
        graph.add_edge(u, v, f'{u}{l[0]}:{v}{l[1]}', locations=l)


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

    def __init__(self, frags=None, limit=25, algorithm=common_sub_strings, use_fragment_order=True, use_all_fragments=False):
        # TODO: allow for the same fragment to be included more than once
        G = _nx.MultiDiGraph()
        G.add_nodes_from((i + 1, {'seq': f}) for (i, f) in enumerate(frags))
        G.add_nodes_from((-(i + 1), {'seq': f.reverse_complement()}) for (i, f) in enumerate(frags))
        edge_pairs = _itertools.combinations(enumerate(frags), 2)

        for (index_first, first), (index_secnd, secnd) in edge_pairs:
            index_first += 1
            index_secnd += 1

            matches_fwd = algorithm(str(first.seq).upper(), str(secnd.seq).upper(), limit)
            for match in matches_fwd:
                add_edges_from_match(match, index_first, index_secnd, first, secnd, G)

            matches_rvs = algorithm(str(first.seq).upper(), reverse_complement(str(secnd.seq).upper()), limit)
            for match in matches_rvs:
                add_edges_from_match(match, index_first, -index_secnd, first, secnd, G)

        self.G = G
        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order
        self.use_all_fragments = use_all_fragments

        return

    def validate_assembly(self, assembly):
        """Function used to filter paths returned from the graph, see conditions tested below.
        """
        is_circular = assembly[0][0] == assembly[-1][1]

        if self.use_all_fragments and (len(assembly) - (1 if is_circular else 0) != len(self.fragments) - 1):
            return False

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
        #
        # Redundant: overlap of 1 and 2 is at the same spot as overlap of 2 and 3
        #    -- A --
        #  gtatcgtgt
        #   catcgtgtactgtcatattc
        #   catcgtgtactgtcatattc
        #   -- B ---
        if is_circular:
            edge_pairs = zip(assembly, assembly[1:] + assembly[:1])
        else:
            edge_pairs = zip(assembly, assembly[1:])

        for (u1, v1, key1), (u2, v2, key2) in edge_pairs:
            left_edge = self.G[u1][v1][key1]['locations']
            right_edge = self.G[u2][v2][key2]['locations']

            # Incompatible as described in figure above
            if left_edge[1].parts[-1].end >= right_edge[0].parts[0].end:
                return False

        return True

    def get_linear_assemblies(self):
        # The constrain for linear assemblies is that the first node is in the initial orientation

        # Copy the graph since we will add the begin and end nodes
        G = _nx.MultiDiGraph(self.G)
        G.add_nodes_from(['begin', 'end'])
        if self.use_fragment_order:
            G.add_edge('begin', 1)
            G.add_edge('begin', -1)
            G.add_edge(len(self.fragments), 'end')
            G.add_edge(-len(self.fragments), 'end')
        else:
            for node in filter(lambda x: type(x) == int, G.nodes):
                if not self.use_fragment_order and node > 0:
                    G.add_edge('begin', node)
                G.add_edge(node, 'end')

        return list(filter(self.validate_assembly, map(lambda x : x[1:-1], _nx.all_simple_edge_paths(G, 'begin', 'end'))))

    def cycle2circular_assemblies(self, cycle):
        # Feels like this should be built-in to networkx, or that find_cycle would return all possible cycles
        combine = list()
        for u, v in zip(cycle, cycle[1:] + cycle[:1]):
            combine.append([(u, v, key) for key in self.G[u][v]])
        return list(_itertools.product(*combine))

    def get_circular_assemblies(self):
        # The constrain of circular sequence is that the first node is the first fragment in its initial orientation
        sorted_cycles = map(circular_permutation_min_abs, _nx.cycles.simple_cycles(self.G))
        sorted_cycles = filter(lambda x: x[0] > 0, sorted_cycles)
        assemblies = sum(map(self.cycle2circular_assemblies, sorted_cycles),[])

        return filter(self.validate_assembly, assemblies)

    def edge_representation2fragment_representation(self, assembly):
        """
        Turn this kind of edge representation fragment 1, fragment 2, right edge on 1, left edge on 2
        a = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b', 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
        Into this: fragment 1, left edge on 1, right edge on 1
        b = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
        """

        is_circular = assembly[0][0] == assembly[-1][1]
        if is_circular:
            temp = assembly[-1:] + assembly
        else:
            temp = [(None, assembly[0][0], None)] + assembly + [(assembly[-1][1], None, None)]
        edge_pairs = zip(temp, temp[1:])
        alternative_representation = list()
        for (u1, v1, key1), (u2, v2, key2) in edge_pairs:
            start_location = None if u1 is None else self.G[u1][v1][key1]['locations'][1]
            end_location = None if v2 is None else self.G[u2][v2][key2]['locations'][0]
            alternative_representation.append((v1, start_location, end_location))

        return alternative_representation

    def execute_assembly(self, assembly_edges):

        fragment_overlaps = [len(self.G[u][v][key]['locations'][1]) for u, v, key in assembly_edges]
        assembly_fragments = self.edge_representation2fragment_representation(assembly_edges)

        fragments = list()
        for node, start_location, end_location in assembly_fragments:
            seq = self.G.nodes[node]['seq']
            start = 0 if start_location is None else start_location.parts[0].start
            end = None if end_location is None else end_location.parts[-1].end
            fragments.append(seq[start:end])

        out_dseqrecord = _Dseqrecord(fragments[0])

        for fragment, overlap in zip(fragments[1:], fragment_overlaps):
            new_features = [f._shift(len(out_dseqrecord)-overlap) for f in fragment.features]
            out_dseqrecord = _Dseqrecord(out_dseqrecord.seq + fragment.seq[overlap:], features=out_dseqrecord.features + new_features)

        # For circular assemblies, close the loop and wrap origin-spanning features
        if assembly_fragments[0][1] != None:
            overlap = fragment_overlaps[-1]
            out_dseqrecord = _Dseqrecord(out_dseqrecord.seq[:-overlap], features=out_dseqrecord.features, circular=True)
            for feature in out_dseqrecord.features:
                if feature.location.parts[0].start >= len(out_dseqrecord) or feature.location.parts[-1].end > len(out_dseqrecord):
                    # Wrap around the origin
                    feature.location = _shift_location(feature.location, 0, len(out_dseqrecord))

        return out_dseqrecord

    def assemble_linear(self):
        assemblies = self.get_linear_assemblies()
        return list(map(self.execute_assembly, assemblies))

    def assemble_circular(self):
        assemblies = self.get_circular_assemblies()
        return list(map(self.execute_assembly, assemblies))

    def __repr__(self):
        # https://pyformat.info
        return _pretty_str(
            "Assembly\n"
            "fragments..: {sequences}\n"
            "limit(bp)..: {limit}\n"
            "G.nodes....: {nodes}\n"
            "algorithm..: {al}".format(
                sequences=" ".join(
                    "{}bp".format(len(x)) for x in self.fragments
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


'CattctgcGAGGacgatCAtgctcc'

linear_results = (
    _Dseqrecord("AacgatCAtgctccTAAattctgcGAGGacgatG", name="abc"),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatG", name="ac_rc"),
    _Dseqrecord("AacgatG", name="ac"),
)


circular_results = (
    _Dseqrecord("acgatCAtgctccTAAattctgcGAGG", name="abc", circular=True),
)