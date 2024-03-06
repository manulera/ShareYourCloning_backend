from fastapi.testclient import TestClient
import unittest
import shutil
from pydantic_models import ManuallyTypedSource
import os

# activate the stubs
os.environ['RECORD_STUBS'] = '1'
from main import app

client = TestClient(app)

class StubRouteTest(unittest.TestCase):

    # DO this before each test
    def setUp(self):
        # remove the stubs folder
        shutil.rmtree('stubs', ignore_errors=True)

    # Remove the stubs folder after each test
    def tearDown(self):
        shutil.rmtree('stubs', ignore_errors=True)

    def test_stub_route(self):
        source = ManuallyTypedSource(
            user_input='ATGC',
        )

        response = client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)

        self.assertTrue(os.path.exists('stubs/manually_typed/'))
        self.assertTrue(len(os.listdir('stubs/manually_typed/')), 1)
        # Also works for 422 response
        source.user_input = 'io'
        response = client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 422)
        self.assertTrue(len(os.listdir('stubs/manually_typed/')), 2)

if __name__ == '__main__':
    unittest.main()
