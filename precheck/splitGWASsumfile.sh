#split GWAS summary file by chromosomes

#cd yalan/variants/

chrlist=({1..22} 'O')
for n in "${chrlist[@]}"; do
	# test using file "dubois_test.tsv" *which column in the GWAS summary file is rsID
	(head -n1 dubois_test.tsv && awk 'NR==FNR {a[$3]; next} $3 in a {print $0}' variants_rsID_coords-chr$n.tsv dubois_test.tsv) > dubois-chr$n.tsv
	echo "chr$n is done!"
	done



# to submit the job:
#$gsub -m 10 -n splitgwas -d /homes/yalan/yalan/variants/ -q production-rh74 'sh splitGWASsumfile.sh' -r
