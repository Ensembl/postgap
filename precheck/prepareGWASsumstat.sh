#preworks before running postgap.py

#cd yalan/variants/

filename="dubois_test.tsv"

# save the first line into to a variable for the following checking step
checkline=$(head -n1 $filename)

# create an array to split filename for renaming new generated GWAS files
IFS='.' read -r -a fnarr <<< "$filename"
filetype="${fnarr[-1]}" #eg:tsv
fnstart="${filename//."${fnarr[-1]}"}" #eg:dubois_test

#1. check the very important three colnames: variant_id, beta, p-value
if grep -q -P '^(?=.*variant_id)(?=.*beta)(?=.*p-value)' <<< "$checkline"; then
	echo "colnames checking passed"
else
	echo -e "Please make sure all the following information is included in your data, and please specify them with assigned column names:\nrsid -> variant_id\nbeta -> beta\np -> p-value"
	#exit
fi

#2. check if the delimiter is tab
if grep -q $'\t' <<< "$checkline"; then
	echo "delimiter checking passed"
else
	echo "$filename is not separated by tab"
	seplist=(' ' ',' ';' ':' '|')
	t=0
	for n in "${seplist[@]}"; do
		if grep -q "$n" <<< "$checkline"; then
			echo "$filename is separated by $n"
			# replace delimiter with tab, and save into a new file
			mv $filename backup_$filename
			cat backup_$filename | tr "$n" '\t' > $filename
			echo "delimiter replacement is done"
			rm backup_$filename
			break
		else
			echo "$filename is not separated by $n"
			((t++))
		fi
	done
	if [[ "$t" == "${#seplist[@]}" ]]; then
		echo "unknown delimiter"
		#exit
	fi
fi

#3. split GWAS summary file by chromosomes
chrlist=({1..22} 'O')
for n in "${chrlist[@]}"; do
	# change the column number!!! *which column in the GWAS summary file is rsID
	(head -n1 $filename && awk 'NR==FNR {a[$3]; next} $3 in a {print $0}' variants_rsID_coords-chr$n.tsv $filename) > $fnstart-chr$n.$filetype
	echo "chr$n is done!"
done
