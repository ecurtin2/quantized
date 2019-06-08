FROM python:3.7
RUN apt-get update -y \
  && apt-get install -y pandoc

COPY . /code
WORKDIR /code
RUN pip install .[all]
