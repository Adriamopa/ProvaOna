#!/bin/bash

wd="$(pwd)"

inv="$1"
 

## Test
# condition="Standard" or condition="Unphased"
# inv="HsInv0015" 

mkdir -p $wd/../VCFs/Phased_Phase3
mkdir -p $wd/../VCFs/Phased_Phase3/$inv

chr="$(grep "$inv[[:space:]]" $wd/../VCFs/Phase3/Inversions_coordinates_hg19.csv | cut -f2 | sed 's/chr//g')"
BP1=$(($(grep "$inv[[:space:]]" $wd/../VCFs/Phase3/Inversions_coordinates_hg19.csv | cut -f3)-200000))
BP2=$(($(grep "$inv[[:space:]]" $wd/../VCFs/Phase3/Inversions_coordinates_hg19.csv | cut -f6)+200000))
Midpoint=$(($(($(grep "$inv[[:space:]]" $wd/../VCFs/Phase3/Inversions_coordinates_hg19.csv | cut -f4)+$(grep "$inv[[:space:]]" $wd/../VCFs/Phase3/Inversions_coordinates_hg19.csv | cut -f5)))/2))

if [ $chr == "X" ]
then
	vcf="$wd/../VCFs/Phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
else
	vcf="$wd/../VCFs/Phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
fi

nline=$(head -1 $wd/../VCFs/Phase3/GTypesINVs.csv | sed 's/\t/\n/g' | cat -n | grep "$inv$" | awk '{print $1}')

## For the inclusion of genotypes
cat <(cut -f1,$nline $wd/../VCFs/Phase3/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) <(cut -f1 $wd/../VCFs/Phase3/Phase3.panel) | sort | uniq -d > $wd/../VCFs/Phased_Phase3/$inv/Samples
cut -f1,$nline $wd/../VCFs/Phase3/GTypesINVs.csv | grep -f $wd/../VCFs/Phased_Phase3/$inv/Samples | sort -k1 > $wd/../VCFs/Phased_Phase3/$inv/InvGenotypes.txt
samps=$(cat $wd/../VCFs/Phased_Phase3/$inv/Samples)

## The VCF
bcftools view -r $chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/../VCFs/Phased_Phase3/$inv/Region.vcf

#bcftools annotate --rename-chrs $wd/../VCFs/Phase3/chr_names.txt  $wd/../VCFs/Phased_Phase3/$inv/Region.vcf -Ov -o $wd/../VCFs/Phased_Phase3/$inv/RegionMod.vcf	

std_line=$(echo -e "$chr\t$Midpoint\t$inv\tA\tG\t100\tPASS\t.\tGT")
genos=$(cut -f2 $wd/../VCFs/Phased_Phase3/$inv/InvGenotypes.txt)

paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/../VCFs/Phased_Phase3/$inv/Region.vcf
vcf-sort $wd/../VCFs/Phased_Phase3/$inv/Region.vcf > $wd/../VCFs/Phased_Phase3/$inv/RegionSorted.vcf

#cat $wd/../VCFs/Phased_Phase3/$inv/RegionSorted.vcf | sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' > $wd/../VCFs/Phased_Phase3/$inv/RegionFilled.vcf

# rm $wd/../VCFs/Phased_Phase3/$inv/Samples $wd/../VCFs/Phased_Phase3/$inv/InvGenotypes.txt $wd/../VCFs/Phased_Phase3/$inv/Region.vcf

java -jar $wd/../Programes/imputation_methods/beagle.22Jul22.46e.jar \
	gt=$wd/../VCFs/Phased_Phase3/$inv/RegionSorted.vcf \
	out=$wd/../VCFs/Phased_Phase3/$inv/Output \
	map=$wd/../VCFs/Phase3/Maps/plink.chr${chr}.GRCh37.map
	
gunzip $wd/../VCFs/Phased_Phase3/$inv/Output.vcf.gz

sed -i -E 's/END=[0-9]+/./g' $wd/../VCFs/Phased_Phase3/$inv/Output.vcf
sed -i -E "/##source/ a ##contig=<ID=${chr},length=10000000000000000>" $wd/../VCFs/Phased_Phase3/$inv/Output.vcf

#bcftools annotate --rename-chrs $wd/../VCFs/Phase3/chr_names_inv.txt  $wd/../VCFs/Phased_Phase3/$inv/OutputMod.vcf -Ov -o $wd/../VCFs/Phased_Phase3/$inv/Output.vcf	


bgzip $wd/../VCFs/Phased_Phase3/$inv/Output.vcf

tabix -f -p vcf $wd/../VCFs/Phased_Phase3/$inv/Output.vcf.gz

# rm $wd/../VCFs/Phased_Phase3/$inv/RegionSorted.vcf		
