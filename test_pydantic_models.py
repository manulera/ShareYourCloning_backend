from unittest import TestCase
from pydantic_models import AssemblyJoin, SimpleSequenceLocation, AssemblyJoinComponent, AssemblySource
from assembly2 import assembly2str
from Bio.SeqFeature import SimpleLocation


class DummyFragment:
    def __init__(self, id):
        self.id = id


class AssemblyJoinTest(TestCase):
    def test_str(self):
        join = AssemblyJoin(
            left=AssemblyJoinComponent(
                sequence=1,
                reverse_complemented=False,
                location=SimpleSequenceLocation(start=1, end=10),
            ),
            right=AssemblyJoinComponent(
                sequence=2,
                reverse_complemented=False,
                location=SimpleSequenceLocation(start=20, end=30),
            ),
        )
        self.assertEqual(str(join), '1[1:10]:2[20:30]')
        self.assertEqual(join.__repr__(), '1[1:10]:2[20:30]')

        join.left.reverse_complemented = True

        self.assertEqual(str(join), '-1[1:10]:2[20:30]')
        self.assertEqual(join.__repr__(), '-1[1:10]:2[20:30]')

    def test_join_tuple(self):
        join_tuple_1 = (1, 2, SimpleLocation(0, 10), SimpleLocation(10, 20))
        join1 = AssemblyJoin.from_join_tuple(join_tuple_1)

        self.assertEqual(join1.left.sequence, 1)
        self.assertEqual(join1.left.reverse_complemented, False)
        self.assertEqual(join1.left.location.start, 0)
        self.assertEqual(join1.left.location.end, 10)

        self.assertEqual(join1.right.sequence, 2)
        self.assertEqual(join1.right.reverse_complemented, False)
        self.assertEqual(join1.right.location.start, 10)
        self.assertEqual(join1.right.location.end, 20)

        self.assertEqual(join1.to_join_tuple([DummyFragment(1), DummyFragment(2)]), join_tuple_1)

        join_tuple_2 = (-1, -2, SimpleLocation(0, 10), SimpleLocation(10, 20))
        join2 = AssemblyJoin.from_join_tuple(join_tuple_2)

        self.assertEqual(join2.left.sequence, 1)
        self.assertEqual(join2.left.reverse_complemented, True)
        self.assertEqual(join2.left.location.start, 0)
        self.assertEqual(join2.left.location.end, 10)

        self.assertEqual(join2.right.sequence, 2)
        self.assertEqual(join2.right.reverse_complemented, True)
        self.assertEqual(join2.right.location.start, 10)
        self.assertEqual(join2.right.location.end, 20)

        self.assertEqual(join2.to_join_tuple([DummyFragment(1), DummyFragment(2)]), join_tuple_2)


class AssemblySourceTest(TestCase):
    def test_get_assembly_plan(self):

        join = AssemblyJoin(
            left=AssemblyJoinComponent(
                sequence=4,
                reverse_complemented=False,
                location=SimpleSequenceLocation(start=1, end=10),
            ),
            right=AssemblyJoinComponent(
                sequence=5,
                reverse_complemented=False,
                location=SimpleSequenceLocation(start=20, end=30),
            ),
        )

        asm = AssemblySource(
            id=1,
            input=[2, 3],
            assembly=[join],
        )

        assert assembly2str(asm.get_assembly_plan([DummyFragment(4), DummyFragment(5)])) == "('1[1:10]:2[20:30]',)"
        assert assembly2str(asm.get_assembly_plan([DummyFragment(5), DummyFragment(4)])) == "('2[1:10]:1[20:30]',)"

    def test_from_assembly(self):
        assembly = [
            (1, 2, SimpleLocation(0, 10), SimpleLocation(10, 20)),
            (2, 3, SimpleLocation(0, 10), SimpleLocation(10, 20)),
        ]

        fragments = [DummyFragment(4), DummyFragment(5), DummyFragment(6)]
        assembly_source = AssemblySource.from_assembly(assembly=assembly, fragments=fragments, id=0, circular=False)

        assert [str(join) for join in assembly_source.assembly] == ['4[0:10]:5[10:20]', '5[0:10]:6[10:20]']
