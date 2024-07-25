import tempfile
import shutil
import os


class TemporaryFolderOverride:
    def __init__(self, target_folder):
        self.target_folder_exists = os.path.exists(target_folder)
        self.target_folder = target_folder
        self.backup_folder = None

    def __enter__(self):
        # Create a temporary directory
        if self.target_folder_exists:
            self.backup_folder = tempfile.mkdtemp()
            shutil.move(self.target_folder, self.backup_folder)

        # Create a new empty folder in place of the original
        os.mkdir(self.target_folder)

    def __exit__(self, exc_type, exc_value, traceback):
        # Remove the temporary folder
        shutil.rmtree(self.target_folder)

        if self.target_folder_exists:
            # Restore the original folder
            shutil.copytree(self.backup_folder, self.target_folder)
            shutil.rmtree(self.backup_folder)
