#!/bin/bash

wd="$(pwd)"

mkdir -p $wd/../VCFs/Inversion_Regions/
mkdir -p $wd/../VCFs/Inversion_Regions/Shapeit_Phase

parallel -j4 -u "bash $wd/Extract_regions_pops/run_hg38_shapeit5_phase.sh" < $wd/../VCFs/30X/Invs2Imp
