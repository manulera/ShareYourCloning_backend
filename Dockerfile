# TODO: use alpine in prod, but you need gcc for some of the deps
FROM python:3.11

WORKDIR /api

RUN pip install poetry

COPY ./poetry.lock /api/poetry.lock
COPY ./pyproject.toml /api/pyproject.toml
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev

RUN rm poetry.lock pyproject.toml
RUN pip uninstall --yes poetry

COPY . /api

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]
