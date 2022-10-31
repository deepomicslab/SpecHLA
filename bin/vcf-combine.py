#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: vcf-combine.py
#
#        USAGE: python vcf-combine.py input1.vcf input2.vcf ... > output.vcf  
#
#  DESCRIPTION: Simply combine multiple VCF files into one file 
#
#      OPTIONS: ---
# REQUIREMENTS: PyVCF module (https://github.com/jamescasbon/PyVCF)
#				Bedtools (https://github.com/arq5x/bedtools2)
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (cauyrd@gmail.com), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: Oct 15, 2015
#     REVISION: ---
#===============================================================================
import sys
import vcf
vcf_header = vcf.Reader(filename=sys.argv[1])
vcf_output = vcf.Writer(sys.stdout,vcf_header)
for vcffile in sys.argv[1:]:
	vcf_reader = vcf.Reader(open(vcffile, 'r'))
	for each in vcf_reader:
		vcf_output.write_record(each)
