[![Python tests](https://github.com/manulera/ShareYourCloning_backend_fastAPI/actions/workflows/ci.yml/badge.svg)](https://github.com/manulera/ShareYourCloning_backend_fastAPI/actions/workflows/ci.yml)
# ShareYourCloning Backend API

This API is part of a bigger application, before going further, please go to the [main project readme](https://github.com/manulera/ShareYourCloning), where you can find an introduction.

This python API is built with [FastAPI](https://fastapi.tiangolo.com/) and is for *in silico* cloning.

## Summary

Read [main project readme](https://github.com/manulera/ShareYourCloning) first.

This API provides a series of entry points. The API documentation can be accessed [here](https://shareyourcloning.api.genestorian.org/docs)

## Getting started

### Prerequisites

You should have python 3.9 installed in your machine. For the management of the dependencies I used `pipenv`, if you don't have it:

```
pip install pipenv
```

### Local installation

In the project directory:

```bash
# This should install the dependencies and create a virtual environment
pipenv install

# Activate the virtual environment
pipenv shell

```

Now you should be able to run the api in the debug mode by doing:

```bash
# These commands are also in the file run_app.sh
export FLASK_APP=app.py
export FLASK_ENV=development
flask run
```

If you go to the url of the flask application (by default [http://127.0.0.1:5000/](http://127.0.0.1:5000/)), you should see a greeting json response. That means the application is working.

You can see some example requests to the API in `examples/request_to_restriction_sites.py`. However, at this point, the goal is to use it as a backend for [ShareYourCloning](https://github.com/manulera/ShareYourCloning) UI.

## My settings for vscode

If you are interested in the settings that I use for vscode, you can create a folder in the directory of the project called `.vscode`, and create a `settings.json` as below.

You will have to change `path/to/python/environment/bin/` by the location of the bin folder of the virtual environment created by pipenv. For that, run `pipenv shell` in the project directory to activate the virtual environment (after you have installed the dependencies), and then run `which python`.

```json
{
    "files.exclude": {
        "**/.git": true,
        "**/.svn": true,
        "**/.hg": true,
        "**/CVS": true,
        "**/.DS_Store": true,
        "**/*.pyc": true,
        "**/__pycache__": true
    },
    "python.linting.enabled": true,
    "python.linting.flake8Enabled": true,
    "python.linting.flake8Path": "path/to/python/environment/bin/flake8",
    "python.defaultInterpreterPath": "path/to/python/environment/bin/python",
    "python.pythonPath": "path/to/python/environment/bin/python"

}
```


### Tests
* The type field is set correctly
* The circular field is set correctly