#! /usr/bin/env python
#coding=utf-8
#Author:    zhangpeng
#email:     zhangpeng@grandomics.com
#Version:   1.0
import re, os, sys
import CheckGene
import argparse
import GetTargetFastqFromBam
import Realign

def GetOpt():

    parser = argparse.ArgumentParser()
    parser.add_argument("-G", "--genelist", type=str, action="store",
        default="", help="gene list file.")
    parser.add_argument("-RE", "--refgene", type=str, action="store",
        default = "", 
        help="bed file for refgene")
    parser.add_argument("-CA", "--chipName", type=str, action="store",
        default = "Medexome", choices=["Medexome", "TSO", "Agilent"],
        help="name of chip.")
    parser.add_argument("-O", "--oPrefix", type=str, action="store",
        default = "./out", help="prefix for output.")
    parser.add_argument("-S", "--Seqmule", type=str, action="store",
        default = "", help="dicrectory for sequmule software."
        )
    parser.add_argument("-Q1", "--fq1", type=str, action="store",
        default="", help="path to fq1 file.")
    parser.add_argument("-Q2", "--fq2", type=str, action="store",
        default="", help="path to fq2 file.")
    parser.add_argument("-b", "--bam", type=str, action="store",
        default="", help="bamfile; if not provided, bam will be recreated.")
    #parser.add_argument("-GI", "--GeneInfo", type=str, action="store",
    #    default = "", help="gene info file."
    #    )
    args = parser.parse_args()
    if not args.refgene or not args.genelist:
        parser.print_help()
        sys.exit(0)

    if not args.genelist or not args.refgene or not args.Seqmule or not args.fq1 or not args.fq2:
        parser.print_help()
        sys.exit(0)
    return(args)

if __name__ == "__main__":
    args = GetOpt()
    dictTargetGene = CheckGene.CheckGene(args)
    #print dictTargetGene
    GetTargetFastqFromBam.Pipe(args, dictTargetGene)
    Realign.Realign(args)
    
