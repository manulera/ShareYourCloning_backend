# TODO: stage the build to not have pipenv in production, which is not needed
FROM python:3.9

WORKDIR /api

# Use the pipenv files to create a requirements.txt file with the dependencies
# since we do not want to use virtual environments in production
# The pipenv lock already exclude dev dependencies

RUN pip install pipenv

COPY ./Pipfile /api/Pipfile
COPY ./Pipfile.lock /api/Pipfile.lock

RUN pipenv lock --keep-outdated --requirements > requirements.txt
# Keeping the pipfile files would give an error
RUN rm Pipfile Pipfile.lock
RUN pip install --upgrade --no-cache-dir -r requirements.txt

COPY . /api

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]
