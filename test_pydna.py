import unittest


class PydnaTest(unittest.TestCase):

    def test_cut_product_order(self):

        from Bio.Restriction.Restriction import CommOnly, RestrictionBatch
        from pydna.dseqrecord import Dseqrecord
        from pydna.dseq import Dseq

        enzyme_names = ['BamHI', 'EcoRV']
        enzymes = [CommOnly.format(e) for e in enzyme_names]
        first_product_sequence = list()
        # With an unpacked list of enzymes ===================

        # seq_object = Dseq('AAAGGATCCAAAAGATATCAAAAA', linear=False)
        # first_product_sequence.append(str(seq_object.cut(*enzymes)[0]))
        # first_product_sequence.append(str(seq_object.cut(*enzymes[::-1])[0]))

        # seq_object = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', linear=False)
        # first_product_sequence.append(str(seq_object.cut(*enzymes)[0].seq))
        # first_product_sequence.append(str(seq_object.cut(*enzymes[::-1])[0].seq))

        # With a RestrictionBatch ==============================

        seq_object = Dseq('AAAGGATCCAAAAGATATCAAAAA', linear=False)
        first_product_sequence.append(str(seq_object.cut(RestrictionBatch(enzymes))[0]))
        first_product_sequence.append(str(seq_object.cut(RestrictionBatch(enzymes[::-1]))[0]))

        seq_object = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', linear=False)
        first_product_sequence.append(str(seq_object.cut(RestrictionBatch(enzymes))[0].seq))
        first_product_sequence.append(str(seq_object.cut(RestrictionBatch(enzymes[::-1]))[0].seq))

        # Check if all products are identical
        assert first_product_sequence.count(first_product_sequence[0]) == len(first_product_sequence)