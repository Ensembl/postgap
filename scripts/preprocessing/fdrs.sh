cat ~/lustre2/CTTV24/databases/Fantom5.txt | sed -e 's/\t[^\t;]*;/\t/' -e 's/;.*$//' | python STOPGAP_FDR.py > ~/lustre2/CTTV24/databases/Fantom5.fdrs
cat ~/lustre2/CTTV24/databases/DHS.txt | python STOPGAP_FDR.py > ~/lustre2/CTTV24/databases/DHS.fdrs
