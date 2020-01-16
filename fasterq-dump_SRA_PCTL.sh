#!/bin/bash

#author: Michael Skaro
#The pupose of this script is to download the publically availbale prostate cancer transcriptomics data. This data
# was made available on the SRA in December of 2019. 

# I am not sure the faster dumpwill work, so we will piece meal the operations to make sure each file is 
# downloaded successfully before moving onto the next project. Hence we have done 28 individual faster1-dumps instead
# of one large array of dumps. 

#date: 1/13/2020


fasterq-dump	ERS025246;
fasterq-dump	ERS025234;
fasterq-dump	ERS025232;
fasterq-dump	ERS025236;
fasterq-dump	ERS025243;
fasterq-dump	ERS025248;
fasterq-dump	ERS025245;
fasterq-dump	ERS025223;
fasterq-dump	ERS025241;
fasterq-dump	ERS025240;
fasterq-dump	ERS025244;
fasterq-dump	ERS025222;
fasterq-dump	ERS025225;
fasterq-dump	ERS025238;
fasterq-dump	ERS025237;
fasterq-dump	ERS025247;
fasterq-dump	ERS025227;
fasterq-dump	ERS025230;
fasterq-dump	ERS025229;
fasterq-dump	ERS025224;
fasterq-dump	ERS025242;
fasterq-dump	ERS025233;
fasterq-dump	ERS025226;
fasterq-dump	ERS025228;
fasterq-dump	ERS025239;
fasterq-dump	ERS025235;
fasterq-dump	ERS025221;
fasterq-dump	ERS025231;




#END

