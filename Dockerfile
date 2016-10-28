FROM ubuntu:xenial

RUN apt-get update && apt-get install -yq --no-install-recommends \
    libncurses5 libncurses5-dev \
    zlib1g-dev \
    python python-pip curl \
    && rm -rf /var/lib/apt/lists/*

RUN curl -sSL https://github.com/vcftools/vcftools/releases/download/v0.1.14/vcftools-0.1.14.tar.gz | tar -xzf \
      && cd vcftools \
      && ./autogen.sh \
      && ./configure \
      && make \
      && cp src/cpp/vcftools /bin \
      && rm -rf -- vcftools

RUN pip install pybedtools==0.7.8 \
    && pip install plink==1.07

CMD ["python", "scripts/gwas_to_genes.py", "--disease", "autism"]
