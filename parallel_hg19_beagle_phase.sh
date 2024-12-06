#!/bin/bash

wd="$(pwd)"

mkdir -p $wd/../VCFs/Phased_Phase3/

parallel -j16 -u "bash $wd/Extract_regions_pops/run_hg19_beagle_phase.sh" < $wd/../VCFs/Phase3/InvsWOTags.txt
parallel -j16 -u "bash $wd/Extract_regions_pops/run_hg19_beagle_phase.sh" < $wd/../VCFs/Phase3/InvsWTags.txt
