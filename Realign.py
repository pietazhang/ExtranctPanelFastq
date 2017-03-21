#! /usr/bin/python
import os, re, sys
import argparse


def Realign(args):
    fq1 = os.path.abspath(args.oPrefix+".R1.fastq.gz")
    fq2 = os.path.abspath(args.oPrefix+".R2.fastq.gz")
    
    seqmulePath = args.Seqmule
    cmd = "source /home/zhuying/perl5/perlbrew/etc/bashrc;"
    #/home/zhuying/grandbox/app/seqmule/bin/seqmule_modified_gc pipeline -a $fq1 -b $fq2 -prefix M16A0281 -capture M16A0281.target.bed -e -threads 12 -rg M16A0281 -advanced /home/zhuying/grandbox/app/seqmule/misc/predefined_config/bwa_gatk.config
    configFile = os.path.abspath("%s/misc/predefined_config/bwa_gatk.config"%(args.Seqmule))
    prefix = os.path.basename(args.oPrefix)
    bedFile = os.path.abspath(args.oPrefix+".target.bed")
    #print(args.oPrefix+".target.bed")
    #print(bedFile)
    if os.path.dirname(args.oPrefix):

        os.chdir(os.path.dirname(args.oPrefix))
        #print("pwd")
        #os.system("pwd")
    cmd += "%s/bin/seqmule_modified_gc pipeline -a %s -b %s -prefix %s -capture %s -e -threads 12 -rg %s -advanced %s"%(seqmulePath, fq1, fq2, prefix, bedFile, prefix, configFile)
    #print(cmd)
    os.system(cmd)