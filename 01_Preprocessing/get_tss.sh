#!/usr/bin/bash

# Rf. https://gist.github.com/brentp/3131928

module load racs-eb MySQL
module load bedtools/2.25.0

UPSTREAM=400
INSTREAM=100
ORG=hg19

mysql -A -h ensembldb.ensembl.org -u anonymous -e \
    "select chrom, txStart, txEnd, X.geneSymbol, strand from knownGene as K, kgXref as X WHERE txStart != txEnd AND X.kgID = K.name" \
    | awk -v ups=$UPSTREAM -v ins=$INSTREAM 'BEGIN{OFS=FS="\t"}
           $5 == "-" { print $1,$3-ins,$3+ups,$4 } 
           $5 == "+" { print $1,$2-ins,$2+ups,$4 }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -nms > tss.bed