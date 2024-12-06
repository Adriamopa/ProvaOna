#!/bin/bash

wd="$(pwd)"
condition="$1"

mkdir -p $wd/../VCFs/Inversion_Regions/
mkdir -p $wd/../VCFs/Inversion_Regions/$condition

parallel -j16 -u "bash $wd/Extract_regions_pops/run_hg38_beagle_phase.sh $condition" < $wd/../VCFs/30X/Invs2Imp
