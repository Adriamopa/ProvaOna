#!/bin/bash

wd="$(pwd)"

condition="$1"
inv="$2"


## Test
# condition="Standard" or condition="Unphased"
# inv="HsInv0015" 

mkdir -p $wd/../VCFs/Inversion_Regions/$condition
mkdir -p $wd/../VCFs/Inversion_Regions/$condition/$inv/

chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
BP1=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f5)-200000))
BP2=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f8)+200000))
Midpoint=$(($(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f6)+$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f7)))/2))

if [ $chr == "X" ]
then
	vcf="$wd/../VCFs/30X/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
else
	vcf="$wd/../VCFs/30X/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
fi

nline=$(head -1 $wd/../VCFs/30X/GTypesINVs.csv | sed 's/\t/\n/g' | cat -n | grep "$inv$" | awk '{print $1}')

## For the inclusion of genotypes
cat <(cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) <(cut -f1 $wd/../VCFs/30X/Panel30x) | sort | uniq -d > $wd/../VCFs/Inversion_Regions/$condition/$inv/Samples
cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -f $wd/../VCFs/Inversion_Regions/$condition/$inv/Samples | sort -k1 > $wd/../VCFs/Inversion_Regions/$condition/$inv/InvGenotypes.txt
samps=$(cat $wd/../VCFs/Inversion_Regions/$condition/$inv/Samples)

## The VCF
bcftools view -r chr$chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/../VCFs/Inversion_Regions/$condition/$inv/Region.vcf

bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names.txt  $wd/../VCFs/Inversion_Regions/$condition/$inv/Region.vcf -Ov -o $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionMod.vcf	

std_line=$(echo -e "$chr\t$Midpoint\t$inv\tA\tG\t100\tPASS\t.\tGT")
genos=$(cut -f2 $wd/../VCFs/Inversion_Regions/$condition/$inv/InvGenotypes.txt)

paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionMod.vcf
vcf-sort $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionMod.vcf > $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionSorted.vcf

cat $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionSorted.vcf | sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' > $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionFilled.vcf

rm $wd/../VCFs/Inversion_Regions/$condition/$inv/Samples $wd/../VCFs/Inversion_Regions/$condition/$inv/InvGenotypes.txt $wd/../VCFs/Inversion_Regions/$condition/$inv/Region.vcf $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionMod.vcf $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionSorted.vcf


if ! [[ $condition =~ [Uu]nphased ]]
then

	java -jar $wd/../Programes/imputation_methods/beagle.22Jul22.46e.jar \
		gt=$wd/../VCFs/Inversion_Regions/$condition/$inv/RegionSorted.vcf \
		out=$wd/../VCFs/Inversion_Regions/$condition/$inv/Output \
		map=$wd/../VCFs/30X/Maps/plink.chr${chr}.GRCh38.map
		
	rm  $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionFilled.vcf
	
	gunzip $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz
	
	sed -i -E 's/END=[0-9]+/./g' $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf
	sed -i -E "/##source/ a ##contig=<ID=${chr},length=10000000000000000>" $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf
	
	#bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names_inv.txt  $wd/../VCFs/Inversion_Regions/$condition/$inv/OutputMod.vcf -Ov -o $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf	

	
	bgzip $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf
	
	tabix -f -p vcf $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz
	
	#rm $wd/../VCFs/Inversion_Regions/$condition/$inv/OutputMod.vcf

else
	cat $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionFilled.vcf | awk 'BEGIN{FS=OFS="\t"} {if ($1 ~ /^#/) print; else {$8="."; print}}' | bgzip > $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz
	
	tabix -f -p vcf $wd/../VCFs/Inversion_Regions/$condition/$inv/Output.vcf.gz
	
	rm $wd/../VCFs/Inversion_Regions/$condition/$inv/RegionFilled.vcf

fi			
