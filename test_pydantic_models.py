from unittest import TestCase
from pydantic_models import AssemblyJoin, SimpleSequenceLocation, AssemblyJoinComponent, AssemblySource
from assembly2 import assembly2str
from Bio.SeqFeature import SimpleLocation


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

        self.assertEqual(join1.to_join_tuple(), join_tuple_1)

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

        self.assertEqual(join2.to_join_tuple(), join_tuple_2)


class AssemblySourceTest(TestCase):
    def test_get_assembly_plan(self):

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

        asm = AssemblySource(
            id=1,
            input=[2, 3],
            assembly=[join],
        )

        print(assembly2str(asm.get_assembly_plan()))
