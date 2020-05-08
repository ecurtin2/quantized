FROM python:3.8 as base
ADD requirements-dev.txt .
RUN pip install -r requirements-dev.txt
ADD . /code



