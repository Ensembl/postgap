#check the very important three colnames: variant_id, beta, p-value

#cd yalan/variants/

filename="dubois_test.tsv"

if grep -q "variant_id" <(head -n1 $filename) && grep -q "beta" <(head -n1 $filename) && grep -q "p-value" <(head -n1 $filename); then
	echo colnames checking passed
else
	echo -e "Please make sure all the following information is included in your data, and please specify them with assigned column names:\nrsid -> variant_id\nbeta -> beta\np -> p-value"
	#exit
fi
