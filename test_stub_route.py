from fastapi.testclient import TestClient
import unittest
import shutil
import os
from pydantic_models import ManuallyTypedSource
from pytest import MonkeyPatch


class StubRouteTest(unittest.TestCase):
    # DO this before each test
    def setUp(self):
        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('RECORD_STUBS', '1')
        from main import app

        client = TestClient(app)
        self.client = client
        # remove the stubs folder
        shutil.rmtree('stubs', ignore_errors=True)

    # Remove the stubs folder after each test
    def tearDown(self):
        shutil.rmtree('stubs', ignore_errors=True)

    def test_stub_route(self):
        source = ManuallyTypedSource(
            user_input='ATGC',
        )
        print('reached the test')
        response = self.client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)

        self.assertTrue(os.path.exists('stubs/manually_typed/'))
        self.assertTrue(len(os.listdir('stubs/manually_typed/')), 1)
        # Also works for 422 response
        source.user_input = 'io'
        response = self.client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 422)
        self.assertEqual(len(os.listdir('stubs/manually_typed/')), 2)


if __name__ == '__main__':
    unittest.main()
