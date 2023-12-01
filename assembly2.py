"""Slightly different assembly implementation"""

from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from pydna.utils import rc as _rc
from pydna.utils import memorize as _memorize
from pydna._pretty import pretty_str as _pretty_str
from pydna.common_sub_strings import common_sub_strings
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
import networkx as _nx
from copy import deepcopy as _deepcopy
import itertools as _itertools
import logging as _logging
from Bio.Seq import reverse_complement


def circular_permutation_min_abs(lst):
    """Returns the circular permutation of lst with the smallest absolute value first."""
    min_abs_index = min(range(len(lst)), key=lambda i: abs(lst[i]))
    return lst[min_abs_index:] + lst[:min_abs_index]

def add_edges_from_match(match, index_first, index_secnd, len_first, len_scnd, graph: _nx.MultiDiGraph):
    x_start, y_start, length = match
    # TODO: adapt to circular (maybe SimpleLocation would not work for origin-spanning locations)
    locs = [_SimpleLocation(x_start, x_start + length, 1), _SimpleLocation(y_start, y_start + length, 1)]
    rc_locs = [locs[0]._flip(len_first), locs[1]._flip(len_scnd)]

    combinations = (
        (index_first, index_secnd, locs),
        (index_secnd, index_first, locs[::-1]),
        (-index_first, -index_secnd, rc_locs),
        (-index_secnd, -index_first, rc_locs[::-1]),
    )
    for u, v, l in combinations:
        graph.add_edge(u, v, f'{u}{l[0]}:{v}{l[1]}', locations=l)

    return

def fragments_in_assembly(assembly):
    return [assembly[0][0]] + list(map(lambda x: x[1], assembly))

# def assembly_is_sorted(assembly):
#     fragment_order = list(map(abs, fragments_in_assembly(assembly)))
#     return fragment_order == sorted(fragment_order)

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
        fragment_pairs = _itertools.combinations(enumerate(frags), 2)

        # if use_fragment_order:
        #     temp = list(enumerate(frags))
        #     fragment_pairs = list(zip(temp, temp[1:] + temp[:1]))

        for (index_first, first), (index_secnd, secnd) in fragment_pairs:
            index_first += 1
            index_secnd += 1
            matches_fwd = algorithm(str(first.seq).upper(), str(secnd.seq).upper(), limit)
            for match in matches_fwd:
                add_edges_from_match(match, index_first, index_secnd, len(first), len(secnd), G)

            matches_rvs = algorithm(str(first.seq).upper(), reverse_complement(str(secnd.seq).upper()), limit)
            for match in matches_rvs:
                add_edges_from_match(match, index_first, -index_secnd, len(first), len(secnd), G)


        self.G = G
        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order
        self.use_all_fragments = use_all_fragments

        return

    def validate_assembly(self, assembly):

        is_circular = assembly[0][0] == assembly[-1][1]

        # The assembly must include all fragments (this could be made optional)
        # TODO: handle circular assemblies
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

    def execute_assembly(self, assembly):
        pos = 0
        first_fragment = assembly[0][0]
        out_dseq = self.G.nodes[first_fragment]['seq']
        for left, right, key in assembly:

            left_location, right_location = self.G[left][right][key]['locations']
            right_seq = self.G.nodes[right]['seq']
            print(left_location, right_location)
            print('a', pos)
            # TODO: shouldn' this be start?
            pos += left_location.parts[-1].end
            print(pos)
            if right == first_fragment:
                # Go to circular joining
                break

            new_dseq = out_dseq.seq[:pos] + right_seq.seq[right_location.parts[-1].end:]

            overlap_length = len(left_location)

            # This could be a separate Dseqrecord method
            new_features = right_seq[right_location.parts[0].start:].features
            new_features = [f._shift(pos-overlap_length) for f in right_seq[right_location.parts[0].start:].features]

            out_dseq = _Dseqrecord(new_dseq, features=out_dseq.features + new_features)
            pos += left_location.parts[-1].end

        else:
            return out_dseq

        # Close the circle in circular assemblies
        # Cut the dseqrecord at the joining point, including the overlap at the beginning and end.
        # E.g. GGGxxxxxxxxxxxxxxGGG
        print(out_dseq.seq, len(out_dseq.seq), left_location)
        print(pos)
        print(right_location.parts[-1].start,left_location.parts[-1].end + pos)
        out_dseq = out_dseq[right_location.parts[-1].start:left_location.parts[-1].end + pos]
        overlap_length = len(left_location)
        new_features = _deepcopy(out_dseq.features)
        out_dseq = _Dseqrecord(out_dseq.seq[:-overlap_length], features=new_features)
        # Wrap origin-spanning features
        print('a', len(out_dseq))
        for feature in new_features:
            print(feature.location.start)
            if feature.location.start >= len(out_dseq):
                feature = feature._shift(len(out_dseq))

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


'CattctgcGAGGacgatCAtgctcc'

linear_results = (
    _Dseqrecord("AacgatCAtgctccTAAattctgcGAGGacgatG", name="abc"),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatG", name="ac_rc"),
    _Dseqrecord("AacgatG", name="ac"),
)


circular_results = (
    _Dseqrecord("acgatCAtgctccTAAattctgcGAGG", name="abc", circular=True),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatTTA", name="abc_rc", circular=True),
)