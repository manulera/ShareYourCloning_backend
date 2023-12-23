"""Slightly different assembly implementation"""

from pydna.utils import shift_location as _shift_location, flatten
from pydna._pretty import pretty_str as _pretty_str
from pydna.common_sub_strings import common_sub_strings as common_sub_strings_str
from pydna.common_sub_strings import terminal_overlap as terminal_overlap_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.dseq import Dseq as _Dseq
import networkx as _nx
import itertools as _itertools
from Bio.SeqFeature import SimpleLocation, Location
from dna_functions import sum_is_sticky
from Bio.Seq import reverse_complement
from Bio.Restriction.Restriction import RestrictionBatch, AbstractCut

def ends_from_cutsite(cutsite: tuple[tuple[int,int],AbstractCut], seq: _Dseq):
    if cutsite is None:
        raise ValueError('None is not supported')
    cut_watson, cut_crick = cutsite[0]
    enz = cutsite[1]
    if enz.ovhg < 0:
        # TODO check the edge in circular
        return (
            ("5'", str(seq[cut_watson:cut_crick].reverse_complement()).lower()),
            ("5'", str(seq[cut_watson:cut_crick]).lower()),
        )
    elif enz.ovhg > 0:
        return (
            ("3'", str(seq[cut_crick:cut_watson]).lower()),
            ("3'", str(seq[cut_crick:cut_watson].reverse_complement()).lower()),
        )

    return ('blunt', ''), ('blunt', '')

def restriction_ligation_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, enzymes=RestrictionBatch):
    """Find overlaps. Like in stiky and gibson, the order matters"""
    cuts_x = seqx.seq.get_cutsites(*enzymes)
    cuts_y = seqy.seq.get_cutsites(*enzymes)
    matches = list()
    for cut_x, cut_y in _itertools.product(cuts_x, cuts_y):
        overlap = sum_is_sticky(
            ends_from_cutsite(cut_x, seqx.seq)[0],
            ends_from_cutsite(cut_y, seqy.seq)[1]
        )
        if overlap:
            left_x = cut_x[0][0] if cut_x[1].ovhg < 0 else cut_x[0][1]
            left_y = cut_y[0][0] if cut_y[1].ovhg < 0 else cut_y[0][1]
            matches.append((left_x, left_y, overlap))
    return matches


def common_sub_strings(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=25):
    return common_sub_strings_str(str(seqx.seq).upper(), str(seqy.seq).upper(), limit)

def terminal_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=25):
    return terminal_overlap_str(str(seqx.seq).upper(), str(seqy.seq).upper(), limit)

def gibson_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=25):
    """
    The order matters, we want alignments like:

    oooo------xxxx
              xxxx------oooo
    Product: oooo------xxxx------oooo

    Not like:

              oooo------xxxx
    xxxx------oooo
    Product (unwanted): oooo
    """
    stringx = str(seqx.seq).upper()
    stringy = str(seqy.seq).upper()
    return [
        m
        for m in common_sub_strings_str(stringx, stringy, limit)
        if (m[1] == 0 and m[0] + m[2] == len(stringx))
    ]

def sticky_end_sub_strings(seqx: _Dseqrecord, seqy: _Dseqrecord, limit=0):
    """For now, if limit 0 / False only full overlaps are considered."""
    overlap = sum_is_sticky(seqx.seq.three_prime_end(), seqy.seq.five_prime_end(), limit )
    if overlap:
        return [(len(seqx)-overlap, 0, overlap)]
    return []

def fill_left(seq: _Dseq):
    """Fill the left overhang of a sequence with the complementary sequence."""
    new_watson = seq.watson
    new_crick = seq.crick

    # Watson 5' overhang
    if seq.ovhg < 0:
        new_crick = new_crick + reverse_complement(seq.watson[:-seq.ovhg])
    # Crick 5' overhang
    elif seq.ovhg > 0:
        new_watson = reverse_complement(seq.crick[-seq.ovhg:]) + new_watson

    return _Dseq(new_watson, new_crick, 0)

def fill_right(seq: _Dseq):
    """Fill the right overhang of a sequence with the complementary sequence."""
    new_watson = seq.watson
    new_crick = seq.crick

    # Watson 3' overhang
    watson_ovhg = seq.watson_ovhg()
    if watson_ovhg < 0:
        new_watson = new_watson + reverse_complement(seq.crick[:-watson_ovhg])

    # Crick 3' overhang
    elif watson_ovhg > 0:
        new_crick = reverse_complement(seq.watson[-watson_ovhg:]) + new_crick

    return _Dseq(new_watson, new_crick, seq.ovhg)

def fill_dseq(seq: _Dseq):
    return fill_left(fill_right(seq))

def reverse_complement_assembly(assembly: list[tuple[int, int, Location, Location]], fragments: list[_Dseqrecord]) -> list[tuple[int, int, Location, Location]]:
    """Complement an assembly, i.e. reverse the order of the fragments and the orientation of the overlaps."""
    new_assembly = list()
    for u, v, locu, locv in assembly:
        f_u = fragments[abs(u)-1]
        f_v = fragments[abs(v)-1]
        new_assembly.append((-v, -u, locv._flip(len(f_v)), locu._flip(len(f_u))))
    return new_assembly[::-1]

def filter_linear_subassemblies(linear_assemblies, circular_assemblies, fragments):
    """Remove linear assemblies which are sub-assemblies of circular assemblies"""
    all_circular_assemblies = circular_assemblies + [reverse_complement_assembly(c, fragments) for c in circular_assemblies]
    filtered_assemblies = [l for l in linear_assemblies if not any(is_sublist(l, c, True) for c in all_circular_assemblies)]
    # I don't think the line below is necessary, but just in case
    # filtered_assemblies = [l for l in filtered_assemblies if not any(is_sublist(reverse_complement_assembly(l, fragments), c, True) for c in all_circular_assemblies)]
    return filtered_assemblies

def remove_subassemblies(assemblies):
    """Filter out subassemblies, i.e. assemblies that are contained within another assembly.

    For example:
        [(1, 2, '1[8:14](+):2[1:7](+)'), (2, 3, '2[10:17](+):3[1:8](+)')]
        [(1, 2, '1[8:14](+):2[1:7](+)')]
    The second one is a subassembly of the first one.
    """

    # Sort by length, longest first
    assemblies = sorted(assemblies, key=len, reverse=True)

    filtered_assemblies = list()
    for assembly in assemblies:
        # Check if this assembly is a subassembly of any of the assemblies we have already found
        if not any(is_sublist(assembly, a) for a in filtered_assemblies):
            filtered_assemblies.append(assembly)

    return filtered_assemblies

def assembly2str(assembly):
    return str(tuple(f'{u}{lu}:{v}{lv}' for u, v, lu, lv in assembly))

def assembly_is_valid(fragments, assembly, is_circular, use_all_fragments, fragments_only_once=True):
    """Function used to filter paths returned from the graph, see conditions tested below.
    """
    if is_circular is None:
        return False

    # Linear assemblies may get begin-1-end, begin-2-end, these are removed here.
    if len(assembly) == 0:
        return False

    if use_all_fragments and len(fragments) != len(set(flatten(map(abs, e[:2]) for e in assembly))):
        return False

    # Here we check whether subsequent pairs of fragments are compatible, for instance:
    # Compatible (overlap of 1 and 2 occurs before overlap of 2 and 3):
    # (1,2,[2:9],[0:7]), (2,3,[12:19],[0:7])
    #    -- A --
    # 1 gtatcgtgt     -- B --
    # 2   atcgtgtactgtcatattc
    # 3               catattcaa
    # Incompatible (overlap of 1 and 2 occurs after overlap of 2 and 3):
    # (1,2,[2:9],[13:20]), (2,3,[0:7],[0:7])
    #                 -- A --
    #  1 -- B --    gtatcgtgt
    #  2 catattcccccccatcgtgtactgt
    #  3 catattcaa
    # Redundant: overlap of 1 and 2 ends at the same spot as overlap of 2 and 3
    # (1,2,[2:9],[1:8]), (2,3,[0:8],[0:8])
    #    -- A --
    #  gtatcgtgt
    #   catcgtgtactgtcatattc
    #   catcgtgtactgtcatattc
    #   -- B ---
    if is_circular:
        # In a circular assembly, first and last fragment must be the same
        if assembly[0][0] != assembly[-1][1]:
            return False
        edge_pairs = zip(assembly, assembly[1:] + assembly[:1])
    else:
        edge_pairs = zip(assembly, assembly[1:])

    for (u1, v1, _, start_location), (u2, v2, end_location, _) in edge_pairs:
        # Incompatible as described in figure above
        fragment = fragments[abs(v1)-1]
        if not fragment.circular and start_location.parts[-1].end >= end_location.parts[0].end:
            return False

    if fragments_only_once:
        nodes_used = [f[0] for f in edge_representation2subfragment_representation(assembly, is_circular)]
        if len(nodes_used) != len(set(map(abs,nodes_used))):
            return False

    return True

def assemble(fragments, assembly, is_circular):
    """Execute an assembly, from the representation returned by get_linear_assemblies or get_circular_assemblies."""

    subfragment_representation = edge_representation2subfragment_representation(assembly, is_circular)

    # Length of the overlaps between consecutive assembly fragments
    fragment_overlaps = [len(e[-1]) for e in assembly]

    subfragments = get_assembly_subfragments(fragments, subfragment_representation)

    out_dseqrecord = _Dseqrecord(subfragments[0])

    for fragment, overlap in zip(subfragments[1:], fragment_overlaps):
        # Shift the features of the right fragment to the left by `overlap`
        new_features = [f._shift(len(out_dseqrecord)-overlap) for f in fragment.features]
        # Join the left sequence including the overlap with the right sequence without the overlap
        # we use fill_right / fill_left so that it works for ligation of sticky ends
        out_dseqrecord = _Dseqrecord(fill_right(out_dseqrecord.seq) + fill_left(fragment.seq)[overlap:], features=out_dseqrecord.features + new_features)

    # For circular assemblies, close the loop and wrap origin-spanning features
    if is_circular:
        overlap = fragment_overlaps[-1]
        # Remove trailing overlap
        out_dseqrecord = _Dseqrecord(fill_dseq(out_dseqrecord.seq)[:-overlap], features=out_dseqrecord.features, circular=True)
        for feature in out_dseqrecord.features:
            if feature.location.parts[0].start >= len(out_dseqrecord) or feature.location.parts[-1].end > len(out_dseqrecord):
                # Wrap around the origin
                feature.location = _shift_location(feature.location, 0, len(out_dseqrecord))

    return out_dseqrecord

def edge_representation2subfragment_representation(assembly, is_circular):
    """
    Turn this kind of edge representation fragment 1, fragment 2, right edge on 1, left edge on 2
    a = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b', 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    Into this: fragment 1, left edge on 1, right edge on 1
    b = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    """

    if is_circular:
        temp = list(assembly[-1:]) + list(assembly)
    else:
        temp = [(None, assembly[0][0], None, None)] + list(assembly) + [(assembly[-1][1], None, None, None)]
    edge_pairs = zip(temp, temp[1:])
    subfragment_representation = list()
    for (u1, v1, _, start_location), (u2, v2, end_location, _) in edge_pairs:
        subfragment_representation.append((v1, start_location, end_location))

    return subfragment_representation

def get_assembly_subfragments(fragments: list[_Dseqrecord], subfragment_representation):
    """From the fragment representation returned by edge_representation2subfragment_representation, get the subfragments that are joined together.

        Subfragments are the slices of the fragments that are joined together

        For example:
        ```
          --A--
        TACGTAAT
          --B--
         TCGTAACGA

        Gives: TACGTAA / CGTAACGA
        ```
        To reproduce:
        ```
        a = Dseqrecord('TACGTAAT')
        b = Dseqrecord('TCGTAACGA')
        f = Assembly([a, b], limit=5)
        a0 = f.get_linear_assemblies()[0]
        print(assembly2str(a0))
        a0_subfragment_rep =edge_representation2subfragment_representation(a0, False)
        for f in get_assembly_subfragments([a, b], a0_subfragment_rep):
            print(f.seq)

        # prints TACGTAA and CGTAACGA
        ```

        Subfragments: `cccccgtatcgtgt`, `atcgtgtactgtcatattc`
    """
    subfragments = list()
    for node, start_location, end_location in subfragment_representation:
        seq = fragments[node-1] if node > 0 else fragments[-node-1].reverse_complement()
        start = 0 if start_location is None else start_location.parts[0].start
        end = None if end_location is None else end_location.parts[-1].end
        # Special case, some of it could be handled by better Dseqrecord slicing in the future
        if seq.circular and start_location == end_location:
            # This could be definitely be done better, but for now it works:
            dummy_cut = ((start, end), type('DynamicClass', (), {'ovhg': start-end})())
            open_seq = seq.apply_cut(dummy_cut, dummy_cut)
            subfragments.append(_Dseqrecord(fill_dseq(open_seq.seq), features=open_seq.features))
            continue
        subfragments.append(seq[start:end])
    return subfragments

def is_sublist(sublist, my_list, my_list_is_cyclic=False):
    """Returns True if sublist is a sublist of my_list (can be treated as cyclic), False otherwise.

    Examples
    --------
    >>> is_sublist([1, 2], [1, 2, 3], False)
    True
    >>> is_sublist([1, 2], [1, 3, 2], False)
    False

    # See the case here for cyclic lists
    >>> is_sublist([3, 1], [1, 2, 3], False)
    False
    >>> is_sublist([3, 1], [1, 2, 3], True)
    True
    """
    n = len(sublist)
    if my_list_is_cyclic:
        my_list = my_list + my_list
    for i in range(len(my_list) - n + 1):
        # Just in case tuples were passed
        if list(my_list[i:i+n]) == list(sublist):
            return True
    return False

def circular_permutation_min_abs(lst):
    """Returns the circular permutation of lst with the smallest absolute value first.

    Examples
    --------
    >>> circular_permutation_min_abs([1, 2, 3])
    [1, 2, 3]
    >>> circular_permutation_min_abs([3, 1, 2])
    [1, 2, 3]
    """
    min_abs_index = min(range(len(lst)), key=lambda i: abs(lst[i]))
    return lst[min_abs_index:] + lst[:min_abs_index]

def add_edges_from_match(match, index_first, index_secnd, first: _Dseqrecord, secnd: _Dseqrecord, graph: _nx.MultiDiGraph):
    """Add edges to the graph from a match returned by an `algorithm` function (see pydna.common_substrings).

    TODO: this is now outdated, and combinations are handled differently.

    The edges added to the graph have the following format: (index_first, index_secnd, key, locations), where:
    - index_first and index_secnd are the indices of the fragments in the input list of fragments.
      The index is positive if the fragment is in the forward orientation, negative if it is in the reverse orientation.
    - key is a string that represents the location of the overlap. In the format: 'u[start:end](strand):v[start:end](strand)'.
    - locations is a list of two FeatureLocation objects, representing the location of the overlap in the first and second fragment.

    All possible combinations of fragments and orientations are added to the graph, see the example below where fragments 1 and 2 share
    an overlap represented by ===. There are two possible fragments recombined, first part of 1 and second part of 2, and first part of 2
    with second part of 1. Edges representing the joining of reverse complements of both fragments are also added.

    ```

    1 ---         ---
          \\     /
           =====
          /     \\
    2 ---         ---
    ```

    ```
    example_fragments = (
        Dseqrecord("AacgatCAtgctcc", name="a"),
        Dseqrecord("TtgctccTAAattctgc", name="b"),
    )

    graph = nx.MultiDiGraph()
    # Nodes represent these fragments in their current orientation
    graph.add_nodes_from([1, 2])

    matches = common_sub_strings(str(example_fragments[0].seq).upper(), str(example_fragments[1].seq).upper(), 5)

    add_edges_from_match(matches[0], 1, 2, example_fragments[0], example_fragments[1], graph)

    for edge in graph.edges:
        u, v, key = edge
        print('u:', u)
        print('v:', v)
        print('key:', key)
        locations = graph.get_edge_data(u, v, key)['locations']
        print('locations: ',locations)
        print()
    ```

    Prints this:
    ```
    u: 1
    v: 2
    key: 1[8:14](+):2[1:7](+)
    locations:  [SimpleLocation(ExactPosition(8), ExactPosition(14), strand=1), SimpleLocation(ExactPosition(1), ExactPosition(7), strand=1)]

    u: 2
    v: 1
    key: 2[1:7](+):1[8:14](+)
    locations:  [SimpleLocation(ExactPosition(1), ExactPosition(7), strand=1), SimpleLocation(ExactPosition(8), ExactPosition(14), strand=1)]

    u: -1
    v: -2
    key: -1[0:6](-):-2[10:16](-)
    locations:  [SimpleLocation(ExactPosition(0), ExactPosition(6), strand=-1), SimpleLocation(ExactPosition(10), ExactPosition(16), strand=-1)]

    u: -2
    v: -1
    key: -2[10:16](-):-1[0:6](-)
    locations:  [SimpleLocation(ExactPosition(10), ExactPosition(16), strand=-1), SimpleLocation(ExactPosition(0), ExactPosition(6), strand=-1)]
    ```


    """
    x_start, y_start, length = match
    # We use shift_location with 0 to wrap origin-spanning features
    locs = [_shift_location(SimpleLocation(x_start, x_start + length, 1), 0, len(first)),
            _shift_location(SimpleLocation(y_start, y_start + length, 1), 0, len(secnd))]
    rc_locs = [locs[0]._flip(len(first)), locs[1]._flip(len(secnd))]

    # For an homology-like assembly, we could do as below, and not do the other combinations,
    # but for a non-symmetrical assembly, such as a sticky end assembly, 1 -> 2 is not the same as 2 -> 1.
    # combinations = (
    #         (index_first, index_secnd, locs),
    #         (index_secnd, index_first, locs[::-1]),
    #         (-index_first, -index_secnd, rc_locs),
    #         (-index_secnd, -index_first, rc_locs[::-1]),
    #     )

    combinations = (
        (index_first, index_secnd, locs),
        (-index_secnd, -index_first, rc_locs[::-1]),
    )

    for u, v, l in combinations:
        graph.add_edge(u, v, f'{u}{l[0]}:{v}{l[1]}', locations=l)

class Assembly:
    """Assembly of a list of linear DNA fragments into linear or circular
    constructs. The Assembly is meant to replace the Assembly method as it
    is easier to use. Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object. Several methods are available for analysis
    of overlapping sequences, graph construction and assembly.

    The assembly contains a directed graph, where nodes represent fragments and
    edges represent overlaps between fragments. :
    - The node keys are integers, representing the index of the fragment in the
    input list of fragments. The sign of the node key represents the orientation
    of the fragment, positive for forward orientation, negative for reverse orientation.
    - The edges contain the locations of the overlaps in the fragments. For an edge (u, v, key):
        - u and v are the nodes connected by the edge.
        - key is a string that represents the location of the overlap. In the format:
        'u[start:end](strand):v[start:end](strand)'.
        - Edges have a 'locations' attribute, which is a list of two FeatureLocation objects,
        representing the location of the overlap in the first and second fragment.

    If fragment 1 and 2 share a subsequence of 6bp, [8:14](+) in fragment 1 and [1:7](+) in fragment 2,
    there will be 4 edges representing that overlap in the graph, for all possible
    orientations of the fragments (see add_edges_from_match for details):
    - `(1, 2, '1[8:14](+):2[1:7](+)')`
    - `(2, 1, '2[1:7](+):1[8:14](+)')`
    - `(-1, -2, '-1[0:6](-):-2[10:16](-)')`
    - `(-2, -1, '-2[10:16](-):-1[0:6](-)')`

    An assembly can be represented as a tuple of graph edges, like this:
    - Linear: ((1, 2, '1[8:14](+):2[1:7](+)'), (2, 3, '2[10:17](+):3[1:8](+)'))
    - Circular: ((1, 2, '1[8:14](+):2[1:7](+)'), (2, 3, '2[10:17](+):3[1:8](+)'), (3, 1, '3[12:17](+):1[1:6](+)'))
    Note that the first and last fragment are the same in a circular assembly.

    The following constrains are applied to remove duplicate assemblies:
    - Circular assemblies: the first subfragment is not reversed, and has the smallest index in the input fragment list.
      use_fragment_order is ignored.
    - Linear assemblies:
        - If use_fragment_order is False, the first fragment is always in the forward orientation.
        - If use_fragment_order is True, the first fragment is always the first fragment in the input list,
        in forward or reverse order, and the last one is the last fragment in the input list, in forward or reverse order.
        This leads

    Parameters
    ----------

    fragments : list
        a list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered
    algorithm : function, optional
        The algorithm used to determine the shared sequences.
    use_fragment_order : bool, optional
        Legacy pydna behaviour: only assemblies that start with the first fragment and end with the last are considered.
    use_all_fragments : bool, optional
        Constrain the assembly to use all fragments.

    Examples
    --------

    from assembly2 import Assembly, example_fragments
    asm = Assembly(example_fragments, limit=5, use_fragment_order=False)
    print('Linear ===============')
    for assembly in asm.get_linear_assemblies():
        print(' ', assembly)
    print('Circular =============')
    for assembly in asm.get_circular_assemblies():
        print(' ', assembly)

    # Prints
    Linear ===============
      ((1, 2, '1[8:14](+):2[1:7](+)'), (2, 3, '2[10:17](+):3[1:8](+)'))
      ((2, 3, '2[10:17](+):3[1:8](+)'), (3, 1, '3[12:17](+):1[1:6](+)'))
      ((3, 1, '3[12:17](+):1[1:6](+)'), (1, 2, '1[8:14](+):2[1:7](+)'))
      ((1, 3, '1[1:6](+):3[12:17](+)'),)
      ((2, 1, '2[1:7](+):1[8:14](+)'),)
      ((3, 2, '3[1:8](+):2[10:17](+)'),)
    Circular =============
      ((1, 2, '1[8:14](+):2[1:7](+)'), (2, 3, '2[10:17](+):3[1:8](+)'), (3, 1, '3[12:17](+):1[1:6](+)'))

    """

    def __init__(self, frags: list[_Dseqrecord], limit=25, algorithm=common_sub_strings, use_fragment_order=True, use_all_fragments=False):
        # TODO: allow for the same fragment to be included more than once?
        G = _nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        G.add_nodes_from((i + 1, {'seq': f}) for (i, f) in enumerate(frags))
        G.add_nodes_from((-(i + 1), {'seq': f.reverse_complement()}) for (i, f) in enumerate(frags))

        # Iterate over all possible combinations of fragments
        edge_pairs = _itertools.combinations(filter(lambda x : x>0, G.nodes), 2)
        for index_first, index_secnd in edge_pairs:
            combinations = (
                (index_first, index_secnd),
                (index_secnd, index_first),
                (index_first, -index_secnd),
                (-index_secnd, index_first)
            )
            for u, v in combinations:
                u_seq = G.nodes[u]['seq']
                v_seq = G.nodes[v]['seq']
                matches = algorithm(u_seq, v_seq, limit)
                for match in matches:
                    add_edges_from_match(match, u, v, u_seq, v_seq, G)

        self.G = G
        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm
        self.use_fragment_order = use_fragment_order
        self.use_all_fragments = use_all_fragments

        return


    def format_assembly_edge(self, assembly_edge):
        """Go from the (u, v, key) to the (u, v, locu, locv) format."""
        u, v, key = assembly_edge
        locu, locv = self.G.get_edge_data(u, v, key)['locations']
        return u, v, locu, locv


    def get_linear_assemblies(self):
        """Get linear assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid). Subassemblies are removed (see remove_subassemblies)."""

        # Copy the graph since we will add the begin and end mock nodes
        G = _nx.MultiDiGraph(self.G)
        G.add_nodes_from(['begin', 'end'])

        if self.use_fragment_order:
            # Path must start with the first fragment and end with the last
            G.add_edge('begin', 1)
            G.add_edge('begin', -1)
            G.add_edge(len(self.fragments), 'end')
            G.add_edge(-len(self.fragments), 'end')
        else:
            # Path must start with forward fragment
            for node in filter(lambda x: type(x) == int, G.nodes):
                if not self.use_fragment_order and node > 0:
                    G.add_edge('begin', node)
                G.add_edge(node, 'end')

        assemblies = [tuple(map(self.format_assembly_edge, x[1:-1])) for x in _nx.all_simple_edge_paths(G, 'begin', 'end')]
        return remove_subassemblies([a for a in assemblies if assembly_is_valid(self.fragments, a, False, self.use_all_fragments)])

    def cycle2circular_assemblies(self, cycle):
        """Convert a cycle in the format [1, 2, 3] (as returned by _nx.cycles.simple_cycles) to a list of all possible circular assemblies.

        There may be multiple assemblies for a given cycle,
        if there are several edges connecting two nodes, for example two overlaps between 1 and 2, and single overlap between 2 and 3 should
        return 3 assemblies. If there was a built-in function in networkx that returned cycles like in all_simple_edge_paths, this would not
        be necessary.
        """
        combine = list()
        for u, v in zip(cycle, cycle[1:] + cycle[:1]):
            combine.append([(u, v, key) for key in self.G[u][v]])
        return [tuple(map(self.format_assembly_edge, x)) for x in _itertools.product(*combine)]

    def get_circular_assemblies(self):
        """Get circular assemblies, applying the constrains described in __init__, ensuring that paths represent
        real assemblies (see assembly_is_valid)."""
        # The constrain of circular sequence is that the first node is the fragment with the smallest index in its initial orientation,
        # this is ensured by the circular_permutation_min_abs function + the filter below
        sorted_cycles = map(circular_permutation_min_abs, _nx.cycles.simple_cycles(self.G))
        sorted_cycles = filter(lambda x: x[0] > 0, sorted_cycles)
        # cycles.simple_cycles returns lists [1,2,3] not assemblies, see self.cycle2circular_assemblies
        assemblies = sum(map(self.cycle2circular_assemblies, sorted_cycles),[])
        return [a for a in assemblies if assembly_is_valid(self.fragments, a, True, self.use_all_fragments)]

    def format_insertion_assembly(self, assembly):
        edge_pair_index = list()

        for i, ((u1, v1, _, start_location), (u2, v2, end_location, _)) in enumerate(zip(assembly, assembly[1:] + assembly[:1])):
            fragment = self.fragments[abs(v1)-1]
            if not fragment.circular and (start_location.parts[-1].start > end_location.parts[0].end or start_location == end_location):
                edge_pair_index.append(i)

        if len(edge_pair_index) != 1:
            return None

        shift_by = len(assembly) - edge_pair_index[0] - 1
        return assembly[shift_by:] + assembly[:shift_by]

    def get_insertion_assemblies(self):
        # We find cycles first
        assemblies = sum(map(self.cycle2circular_assemblies, _nx.cycles.simple_cycles(self.G)),[])
        # We select those that contain exactly only one suitable edge
        assemblies = [b for a in assemblies if (b := self.format_insertion_assembly(a)) is not None]
        # First fragment should be in the + orientation
        assemblies = list(filter(lambda x: x[0][0] > 0, assemblies))
        return [a for a in assemblies if assembly_is_valid(self.fragments, a, False, self.use_all_fragments, False)]

    def assemble_linear(self):
        """Assemble linear constructs, from assemblies returned by self.get_linear_assemblies."""
        assemblies = self.get_linear_assemblies()
        return [assemble(self.fragments, a, False) for a in assemblies]

    def assemble_circular(self):
        """Assemble circular constructs, from assemblies returned by self.get_circular_assemblies."""
        assemblies = self.get_circular_assemblies()
        return [assemble(self.fragments, a, True) for a in assemblies]

    def assemble_insertion(self):
        """Assemble insertion constructs, from assemblies returned by self.get_insertion_assemblies."""
        assemblies = self.get_insertion_assemblies()
        return [assemble(self.fragments, a, False) for a in assemblies]

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

class PCRAssembly(Assembly):

    def __init__(self, frags: tuple[_Dseqrecord, _Dseqrecord, _Dseqrecord], limit=25, algorithm=common_sub_strings):

        # TODO: allow for the same fragment to be included more than once?
        G = _nx.MultiDiGraph()
        # Add positive and negative nodes for forward and reverse fragments
        forward_primer, template, reverse_primer = frags
        G.add_node(1, seq=forward_primer)
        G.add_node(2, seq=template)
        G.add_node(-2, seq=template.reverse_complement())
        G.add_node(-3, seq=reverse_primer.reverse_complement())

        combinations = (
                (1, 2),
                (1, -2),
                (2, -3),
                (-2, -3)
        )
        for u, v in combinations:
            u_seq = G.nodes[u]['seq']
            v_seq = G.nodes[v]['seq']
            matches = algorithm(u_seq, v_seq, limit)
            # For now we use the same algorithm function, and filter after those in which
            # the primer anneals until its 3'-most base, but a separate algorithm could be used.
            if u == 1:
                matches = filter(lambda x: x[0] + x[2] == len(forward_primer), matches)
            elif v == -3:
                matches = filter(lambda x: x[1] == 0, matches)
            for match in matches:
                add_edges_from_match(match, u, v, u_seq, v_seq, G)

        # These two are constrained
        self.use_fragment_order=True
        self.use_all_fragments=True

        self.G = G
        self.fragments = frags
        self.limit = limit
        self.algorithm = algorithm

        return

    def get_linear_assemblies(self):
        """Adds extra constrains to prevent clashing primers."""
        assemblies = super().get_linear_assemblies()
        # Error if clashing primers
        for a in assemblies:
            edge_pairs = zip(a, a[1:])
            for (u1, v1, _, start_location), (u2, v2, end_location, _) in edge_pairs:
                # Incompatible as described in figure above
                if start_location.parts[-1].end > end_location.parts[0].start:
                    raise ValueError('Clashing primers in assembly ' + assembly2str(a))

        return assemblies


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