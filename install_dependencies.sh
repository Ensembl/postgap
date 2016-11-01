#!/bin/bash

mkdir -p bin

## Tabix, bgzip
git clone git@github.com:samtools/htslib.git
git checkout 1.3.2
cd htslib
make
cp tabix bgzip ../bin
cd ..

## vcfkeepsamples
git clone git@github.com:vcflib/vcflib.git
git checkout v1.0.0-rc0
cd vcflib
make
cp vcfkeepsamples ../bin
cd ..

## bcftools
git clone git@github.com:samtools/bcftools.git
gti checkout 1.3.1
cd bcftools
make
cp bcftools ../bin
cd ..

# vcftools
git clone git@github.com:vcftools/vcftools.git
git checkout v0.1.14
cd vcftools
./autogen.sh
./configure
make
cp src/cpp/vcftools ../bin
cd ..

# bedtools
git clone https://github.com/arq5x/bedtools2.git
git checkout v2.26.0
cd bedtools2
make
cp bin/bedtools ../bin
cd ..

# pybedtools v0.7.8
pip install pybedtools

# Plink v1.07
git clone https://github.com/chrchang/plink-ng.git
cd plink-ng
mv Makefile.std Makefile
./plink_first_compile
cp plink ../bin
cd ..
