from unittest import TestCase
from pydantic_models import AssemblyJoin, SimpleSequenceLocation


class AssemblyJoinTest(TestCase):
    def test_str(self):
        join = AssemblyJoin(
            left_fragment=1,
            right_fragment=2,
            left_location=SimpleSequenceLocation(start=1, end=10),
            right_location=SimpleSequenceLocation(start=20, end=30),
        )
        self.assertEqual(str(join), '1[1:10]:2[20:30]')
        self.assertEqual(join.__repr__(), '1[1:10]:2[20:30]')
