#!/bin/bash

wd="$(pwd)"

pops="EUR AFR EAS"

for inv in $(ls $wd/../VCFs/Phased_Phase3)
do
	for pop in $pops
	do
		commSamps=$(cat <(zgrep "#CHROM" $wd/../VCFs/Phased_Phase3/$inv/Output.vcf.gz | sed 's/\t/\n/g') \
		<(grep $pop $wd/../VCFs/Phase3/Phase3.panel | cut -f1) \
		| sort | uniq -d)

		bcftools view -s $(echo $commSamps | tr " " ",") -Oz $wd/../VCFs/Phased_Phase3/$inv/Output.vcf.gz -o $wd/../VCFs/Phased_Phase3/$inv/${pop}.vcf.gz
	done
done









