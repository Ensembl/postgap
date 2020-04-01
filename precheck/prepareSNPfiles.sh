#prework about SNP files

#cd yalan/variants/

#1. download all the SNP files from ensembl in one goal
wget ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz

#2. unzip and delete .gz files
gzip -d homo_sapiens-chr*.vcf.gz

#3. extract SNPs coordinates
chrlist=({1..22} 'X' 'Y' 'MT')
for n in "${chrlist[@]}"; do
	# use only the first 50 lines of SNP files for test:
	#head -n50 homo_sapiens-chr$n.vcf | grep -v ^## | cut -f 1-3 > test-chr$n.vcf

	# use the whole SNP files:
	#1. grep - to remove annotations at the top that starts with ##,
	#		   but keep headers, which lines start with #
	#2. cut - to select the first 3 columns
	grep -v ^## homo_sapiens-chr$n.vcf | cut -f 1-3 > variants_coords-chr$n.tsv
	#3. head & tail - remain the header, then remove SNPs without rsID from the rest lines
	#4. sort - to sort by rsID, by numbers after "rs"
	(head -n1 variants_coords-chr$n.tsv && tail -n+2 variants_coords-chr$n.tsv | grep "rs" | sort -t "s" -k2n) > variants_rsID_coords-chr$n.tsv
	#5. delete intermediate files
	rm homo_sapiens-chr$n.vcf
	rm variants_coords-chr$n.tsv

	echo "chr$n is done!"
	done

#4. merge chrX, chrY and chrMT together into one called chrO - other chromosomes, as they are all small
cat variants_rsID_coords-chrX.tsv <(tail -n+2 variants_rsID_coords-chrY.tsv) <(tail -n+2 variants_rsID_coords-chrMT.tsv) > variants_rsID_coords-chrO.tsv

# to submit the job:
#$gsub -m 10 -n snpfiles -d /homes/yalan/yalan/variants/ -q production-rh74 'sh prepareSNPfiles.sh' -r
