FROM python:3.7

RUN apt-get update &&\
    apt-get install --no-install-recommends -y \
    build-essential \
    graphviz \
    libopenmpi-dev
RUN pip install --no-cache-dir ford

COPY entrypoint.sh /
RUN ["chmod", "+x", "entrypoint.sh"]
ENTRYPOINT ["/entrypoint.sh"]