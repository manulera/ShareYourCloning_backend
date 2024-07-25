import unittest
from pytest import MonkeyPatch
from importlib import reload
from fastapi.testclient import TestClient
from utils import TemporaryFolderOverride
import shutil
import glob
import os
import tempfile


class TestTemporaryFolderOverride(unittest.TestCase):

    def test_temporary_folder_override(self):
        # File that did not exist before
        with tempfile.TemporaryDirectory() as tempdir:
            override_folder = os.path.join(tempdir, 'dummy_test_file_folder_override')
            with TemporaryFolderOverride(override_folder):
                assert os.path.exists(override_folder)
            assert not os.path.exists(override_folder)

        # File that existed before
        with tempfile.TemporaryDirectory() as tempdir:
            override_folder = os.path.join(tempdir, 'dummy_test_file_folder_override')
            os.mkdir(override_folder)
            # Add a file to the folder
            with open(os.path.join(override_folder, 'dummy_file.txt'), 'w') as f:
                f.write('Dummy file')
            with TemporaryFolderOverride(override_folder):
                assert os.path.exists(override_folder)
                assert not os.path.exists(os.path.join(override_folder, 'dummy_file.txt'))

            assert os.path.exists(override_folder)
            assert os.path.exists(os.path.join(override_folder, 'dummy_file.txt'))


class TestServeFrontend(unittest.TestCase):
    # DO this before each test
    def setUp(self):

        self.folder_override = TemporaryFolderOverride('frontend')
        self.folder_override.__enter__()

        for f in glob.glob('test_files/dummy_frontend/*'):
            print(f)
            if os.path.isdir(f):
                shutil.copytree(f, f'./frontend/{f.split("/")[-1]}', dirs_exist_ok=True)
            else:
                shutil.copy2(f, './frontend')

        # Has to be imported here to get the right environment variable
        MonkeyPatch().setenv('SERVE_FRONTEND', '1')
        import main

        reload(main)
        client = TestClient(main.app)
        self.client = client

    # DO this after each test
    def tearDown(self):
        self.folder_override.__exit__(None, None, None)
        MonkeyPatch().setenv('SERVE_FRONTEND', '0')
        import main

        reload(main)

    def test_serve_frontend(self):
        # The index is served at the root
        response = self.client.get('/')
        assert response.status_code == 200
        assert '<noscript>You need to enable JavaScript to run this app.</noscript>' in response.text

        # The rest of files can be accessed:
        for file in ['config.json', 'favicon.ico', 'robots.txt', 'logo192.png', 'logo512.png', 'manifest.json']:
            response = self.client.get(f'/{file}')
            assert response.status_code == 200

        # Can still make a normal request (e.g. enzymes)
        response = self.client.get('/restriction_enzyme_list')
        assert response.status_code == 200
        assert 'EcoRI' in response.text
