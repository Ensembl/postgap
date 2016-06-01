curl 'http://www.ebi.ac.uk/gwas/api/search?q=text%3A%22diabetes%22' | jq ' .response.docs | map(select(.resourcename == "association")) | .[].beta '
