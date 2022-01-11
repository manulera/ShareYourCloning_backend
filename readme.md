[![Python tests](https://github.com/manulera/ShareYourCloning_backend_fastAPI/actions/workflows/ci.yml/badge.svg)](https://github.com/manulera/ShareYourCloning_backend_fastAPI/actions/workflows/ci.yml)
# ShareYourCloning Backend API

This API is part of a bigger application, before going further, please go to the [main project readme](https://github.com/manulera/ShareYourCloning), where you can find an introduction.

This python API is built with [FastAPI](https://fastapi.tiangolo.com/) and is for *in silico* cloning.

## Summary

Read [main project readme](https://github.com/manulera/ShareYourCloning) first.

This API provides a series of entry points. The API documentation can be accessed [here](https://shareyourcloning.api.genestorian.org/docs). You can use the documentation page to try some request directly on the browser. Otherwise, the API is open for you to make requests from a python script or command line at: [https://shareyourcloning.api.genestorian.org/](https://shareyourcloning.api.genestorian.org/).

## Getting started

If you want to quickly set up a local instance of the frontend and backend of the application, check [getting started in 5 minutes](https://github.com/manulera/ShareYourCloning#timer_clock-getting-started-in-5-minutes) in the main repository.

### Local installation

You should have python 3.9 installed in your machine. For the management of the dependencies I used `pipenv`, if you don't have it:

```
pip install pipenv
```

In the project directory:

```bash
# This should install the dependencies and create a virtual environment
pipenv install

# Activate the virtual environment
pipenv shell

```
Now you should be able to run the api by running:

```bash
# The --reload argument will reload the API if you make changes to the code
uvicorn main:app --reload
```

Then you should be able to open the API docs at [http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs) to know that your API is working.

### Running locally with docker

You can build the docker image and run it:

```bash
docker build -t shareyourcloningapi .
docker run -d --name apicontainer -p 8000:80 shareyourcloningapi
```

If you don't want to download the repository and build the image, you can fetch the latest image from dockerhub (same image that is used in [https://shareyourcloning.api.genestorian.org/](https://shareyourcloning.api.genestorian.org/))

```
docker pull manulera/shareyourcloningapi
docker run -d --name apicontainer -p 8000:80 manulera/shareyourcloningapi
```

The api will be running at `http://localhost:8000`, so you should be able to access the docs at [http://localhost:8000/docs](http://localhost:8000/docs0).

### Connecting to the frontend

If you want to receive requests from the [frontend](https://github.com/manulera/ShareYourCloning_frontend), or from another web application that is not served at `http://localhost:3000`, you must include the url of the frontend application it in the CORS exceptions, by adding it to the list `origins` in `main.py`:

```python
# at the beginning of main.py file
origins = ["http://localhost:3000", "https://shareyourcloning.netlify.app"]
```

Finally, if you are running your api at an address other than `http://127.0.0.1:8000/docs`, you have to configure your frontend to send request to your api address ([see here](https://github.com/manulera/ShareYourCloning_backend#connecting-to-the-frontend)).

## Contributing

Check [contribution guidelines in the main repository](https://github.com/manulera/ShareYourCloning/blob/master/CONTRIBUTING.md).

## Acknowledgements

Thanks to [@maratumba](https://github.com/maratumba) for recommending the usage of FastAPI and for giving some general guidelines for the development.

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
}
```
