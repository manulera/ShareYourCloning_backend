# Backend for ShareYourCloning
# https://github.com/manulera/ShareYourCloning_backend

# BUILDER IMAGE
FROM python:3.11-slim-bookworm as builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y gcc git g++

RUN useradd -ms /bin/bash backend
USER backend
WORKDIR /home/backend

ENV PIP_DISABLE_PIP_VERSION_CHECK=on

# Poetry
# https://python-poetry.org/docs/configuration/#using-environment-variables
# make poetry install to this location
ENV POETRY_HOME="/home/backend/.bin"
# do not ask any interactive question
ENV POETRY_NO_INTERACTION=1
# never create virtual environment automatically, only use env prepared by us
ENV POETRY_VIRTUALENVS_CREATE=false
# this is where our dependencies and virtual environment will live
ENV VIRTUAL_ENV="/home/backend/venv"
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$POETRY_HOME/bin:$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir poetry

COPY ./src ./src
COPY ./poetry.lock .
COPY ./pyproject.toml .
COPY ./README.md .

RUN poetry build
RUN pip install dist/shareyourcloning-*.whl

# FINAL IMAGE
FROM python:3.11-slim-bookworm

# directly output things to stdout/stderr, without buffering
ENV PYTHONUNBUFFERED=1

# create a user to run the app
RUN useradd -ms /bin/bash backend
USER backend
WORKDIR /home/backend

# Accept the COMMIT_SHA as a build argument
ARG COMMIT_SHA
ARG VERSION

# Store the COMMIT_SHA and VERSION in files
RUN echo "${COMMIT_SHA}" > commit_sha.txt
RUN echo "${VERSION}" > version.txt

ENV VIRTUAL_ENV="/home/backend/venv"
COPY --from=builder $VIRTUAL_ENV $VIRTUAL_ENV

ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# For example, ROOT_PATH="/syc"
ENV ROOT_PATH=""
# Only add --root-path if ROOT_PATH is not empty, otherwise uvicorn will throw an error
CMD uvicorn shareyourcloning.main:app --host 0.0.0.0 --port 8000 ${ROOT_PATH:+--root-path ${ROOT_PATH}}
