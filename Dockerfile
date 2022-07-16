# TODO: stage the build to not have pipenv in production, which is not needed
FROM python:3.9

WORKDIR /api

RUN pip install poetry

COPY ./poetry.lock /api/poetry.lock
COPY ./pyproject.toml /api/pyproject.toml
COPY ./dna_functions /api/dna_functions
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev

RUN rm poetry.lock pyproject.toml
RUN pip uninstall --yes poetry

COPY . /api

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]
