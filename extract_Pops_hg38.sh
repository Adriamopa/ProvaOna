#!/bin/bash

wd="$(pwd)"

condition="$1"

## Test
#condition="Standard" "Unphased" "Shapeit_Phase"

pops="EUR AFR EAS"

for inv in $(ls $wd/../VCFs/Inversion_Regions/$condition)
do

	for pop in $pops
	do
		commSamps=$(cat <(zgrep "#CHROM" $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz | sed 's/\t/\n/g') \
		<(grep $pop $wd/../VCFs/30X/Panel30x | cut -f1) \
		| sort | uniq -d)

		bcftools view -s $(echo $commSamps | tr " " ",") -Oz $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz -o $wd/../VCFs/Inversion_Regions/$condition/$inv/${pop}.vcf.gz

	done

done









