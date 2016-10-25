#!/bin/bash

mkdir -p bin

## Tabix, bgzip
git clone git@github.com:samtools/htslib.git
cd htslib
make
cp tabix bgzip ../bin
cd ..

## vcfkeepsamples
git clone git@github.com:vcflib/vcflib.git
cd vcflib
make
cp vcfkeepsamples ../bin
cd ..

## bcftools
git clone git@github.com:samtools/bcftools.git
cd bcftools
make
cp bcftools ../bin
cd ..

# vcftools
git clone git@github.com:vcftools/vcftools.git
cd vcftools
./autogen.sh
./configure
make
cp src/cpp/vcftools ../bin
cd ..

# pybedtools
pip install pybedtools
