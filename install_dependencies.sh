#!/bin/bash

mkdir -p bin

## Tabix, bgzip
git clone https://github.com/samtools/htslib.git
cd htslib
git checkout 1.3.2
make
cp tabix bgzip ../bin
cd ..

## vcfkeepsamples
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
git checkout v1.0.0-rc1
make
cp bin/vcfkeepsamples ../bin
cd ..

## bcftools
git clone https://github.com/samtools/bcftools.git
cd bcftools
git checkout 1.3.1
make
cp bcftools ../bin
cd ..

# vcftools
git clone https://github.com/vcftools/vcftools.git
cd vcftools
git checkout v0.1.14
./autogen.sh
./configure
make
cp src/cpp/vcftools ../bin
cd ..

# bedtools
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git checkout v2.26.0
make
cp bin/bedtools bin/intersectBed ../bin
cd ..

# pybedtools v0.7.8, requests, pandas
pip install pybedtools==0.7.4 requests pandas

# Plink v1.07
git clone https://github.com/chrchang/plink-ng.git
cd plink-ng/1.9
mv Makefile.std Makefile
sed -e 's/^NO_LAPACK =/NO_LAPACK = 1/' -i Makefile
sed -e 's@^// #define NOLAPACK@#define NOLAPACK@' -i plink_common.h
./plink_first_compile
cp plink ../../bin
cd ../..
