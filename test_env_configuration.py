import unittest
from pytest import MonkeyPatch
from importlib import reload
from fastapi.testclient import TestClient
from utils import TemporaryFolderOverride
import shutil


class TestServeFrontend(unittest.TestCase):
    # DO this before each test
    def setUp(self):

        self.folder_override = TemporaryFolderOverride('frontend')
        self.folder_override.__enter__()

        shutil.move('test_files/dummy_frontend', './frontend')

        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('SERVE_FRONTEND', '1')
        import main

        reload(main)
        client = TestClient(main.app)
        self.client = client

    # DO this after each test
    def tearDown(self):
        self.folder_override.__exit__(None, None, None)

    def test_serve_frontend(self):
        response = self.client.get('/')
        assert response.status_code == 200
        assert '<div id="root">Dummy example</div>' in response.text
