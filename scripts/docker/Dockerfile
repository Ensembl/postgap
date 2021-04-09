FROM ubuntu:latest AS builder

RUN apt update && apt install -y --no-install-recommends\
   libvcflib1\
   libbigwig-dev\
   libgsl-dev\
   libhts-dev\
   wget\
   python2-dev\
   libsqlite3-dev\
   gcc\
   g++\
   ca-certificates\
   git\
   make\
   libstdc++-10-dev\
   libcurl4-openssl-dev\
   zlib1g-dev

RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py --no-check-certificate && python2 get-pip.py
RUN pip install pysqlite requests h5py scipy pybedtools

RUN git clone https://github.com/Ensembl/WiggleTools
WORKDIR WiggleTools
RUN   make LIBS='-lwiggletools -lBigWig -lcurl -lz -lhts -lm -lgsl -lpthread' 
RUN   make test

WORKDIR ../postgap
COPY C C
RUN mkdir bin
RUN gcc -Wall -O3 C/ld_vcf.c -I htslib -o bin/ld_vcf -Lhtslib -Wl,-rpath,htslib -lhts -lm

FROM ubuntu:latest
RUN apt update && apt install -y --no-install-recommends \
    libbigwig0 \
    libhts3\
    tabix\
    bcftools\
    vcftools\
    bedtools\
    python2\
    sqlite3\
    libgsl23

COPY   --from=builder /WiggleTools/bin/wiggletools /usr/local/bin/
COPY   --from=builder /postgap/bin/ld_vcf /usr/local/bin/
COPY   --from=builder /usr/local/lib/python2.7 /usr/local/lib/python2.7
COPY   --from=builder /usr/local/lib/python2.7/ /usr/local/lib/python2.7/

COPY POSTGAP.py POSTGAP.py
COPY lib postgaplib
COPY configuration configuration
ENV PYTHONPATH postgaplib

ENTRYPOINT ["python2","POSTGAP.py"]
CMD ["--help"]

