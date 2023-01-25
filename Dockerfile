FROM python:3.8 as build

ARG BUILD_VERSION=3.0.3
ARG BUILD_DATE
ARG VCS_REF

LABEL com.oncoimmunity.build-date=$BUILD_DATE
LABEL com.oncoimmunity.vcs-url="https://github.com/OncoImmunity/spladder"
LABEL com.oncoimmunity.vcs-ref=$VCS_REF
LABEL com.oncoimmunity.version=$BUILD_VERSION

WORKDIR /usr/src/app

RUN pip install --upgrade pip wheel

COPY . .

RUN pip install -r requirements.txt

RUN make install

ENV PATH=/root/.local/bin:$PATH PYTHONUNBUFFERED=1

# procps is required by Nextflow
RUN apt-get update -y \
    && apt-get install -y -o Acquire::Retries=10 --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

CMD ["spladder"]
