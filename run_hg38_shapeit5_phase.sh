#!/bin/bash

wd="$(pwd)"

inv="$1"

# For testing
#inv="HsInv0015"

mkdir -p $wd/../VCFs/Inversion_Regions/Shapeit_Phase
mkdir -p $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/

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

cat <(cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) <(cut -f1 $wd/../VCFs/30X/Panel30x) | sort | uniq -d > $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Samples
cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -f $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Samples | sort -k1 > $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/InvGenotypes.txt
samps=$(cat $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Samples)

bcftools view -r chr$chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Region.vcf

bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names.txt  $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Region.vcf -Ov -o $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionMod.vcf

std_line=$(echo -e "$chr\t$Midpoint\t$inv\tA\tG\t100\tPASS\t.\tGT")
genos=$(cut -f2 $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/InvGenotypes.txt)

paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionMod.vcf

## Change annotation from non SNP variants (<INS>, <INV>, <DEL> ...) to trivial SNP (G T) for SHAPEIT to work
bcftools view $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionMod.vcf | sed -E ' /^[^#]/ s/[[:alpha:]]+\t<[^ ]+>/G\tT/g ' > $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionCurated.vcf

vcf-sort $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionCurated.vcf > $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionSorted.vcf

bcftools +fill-AN-AC $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionSorted.vcf -Ou -o $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotated.bcf
#bcftools norm -d all $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionSorted.vcf -Oz -o $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotatedNormalized.vcf.gz

tabix -f -p bcf $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotated.bcf

rm $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Samples $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/InvGenotypes.txt $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Region.vcf $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionMod.vcf $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionCurated.vcf $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionSorted.vcf

# Phase common variants (caution with --filter-maf parameter)
SHAPEIT5_phase_common	--input $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotated.bcf \
			--region "$chr:$BP1-$BP2"\
			--filter-maf 0.05 \
			--map $wd/../VCFs/30X/Maps/chr${chr}.b38.gmap.gz \
			--output $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionScaffold.bcf \
			--thread 5 \
			--log $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/OutputCommon.log
			
# Phase rare variants (pbwt-modulo parameter care with it)
SHAPEIT5_phase_rare	--input $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotated.bcf \
			--scaffold $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionScaffold.bcf \
			--input-region "$chr:$BP1-$BP2" \
			--scaffold-region "$chr:$BP1-$BP2" \
			--map $wd/../VCFs/30X/Maps/chr${chr}.b38.gmap.gz \
			--pbwt-modulo 0.005 \
			--output $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionPhased.bcf \
			--thread 5 \
			--log $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/OutputRare.log


bcftools annotate --remove FORMAT/PP $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionPhased.bcf -Oz -o $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Output.vcf.gz

tabix -f -p vcf $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/Output.vcf.gz
			
rm $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionAnnotated.bcf* $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionScaffold.bcf* $wd/../VCFs/Inversion_Regions/Shapeit_Phase/$inv/RegionPhased.bcf*

