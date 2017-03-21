#! /usr/bin/env python

import os, sys, re
import argparse
import copy
import logging



def GetTargetGene(geneListFile):
    f = open(geneListFile, "r")
    qGene = re.split("\s+", f.read().strip("\n"))
    resGene = {i.strip("\n"):[] for i in qGene}
    f.close()
    return resGene

def CheckIfMainSymbol(dictTargetGene, geneInfoFile):
    dictTarGeneButNotSym = copy.deepcopy(dictTargetGene)
    print(geneInfoFile)
    f = open(geneInfoFile, "r")
    f.readline()
    for line in f:
        qTmp = re.split("\t", line.strip("\n"))
        if qTmp[0] != "9606":
            continue;
        ######################
        geneMainSym = qTmp[2]
        qAlias = re.split("|", qTmp[4])
        if geneMainSym in dictTarGeneButNotSym:
            dictTarGeneButNotSym.pop(geneMainSym)
        for aliasGene in qAlias:
            if aliasGene in dictTarGeneButNotSym:
                dictTarGeneButNotSym[aliasGene].append(geneMainSym)


    for gene in dictTarGeneButNotSym.keys():
        aliasGene = ",".join(dictTarGeneButNotSym[gene])
        logging.warning("%s not a main gene symbol, substituted by %s"%(gene, aliasGene))
        for alias in dictTarGeneButNotSym[gene]:
            dictTargetGene[alias] = []
    #return dictTargetGene[alias]


    f.close()
def CheckCoverage(dictTargetGene, chipName):
    coverageFile = os.path.dirname(os.path.abspath(sys.argv[0]))+"/config/coverage.%s.txt"%(chipName)
    dictRawTargetGene = copy.deepcopy(dictTargetGene)
    if not os.path.exists(coverageFile):
        sys.stderr.write("cannot find %s"%coverageFile)
        sys.exit(0)
    f = open(coverageFile, "r")
    for line in f:
        geneName, cov = re.split("\t", line.strip("\n"))
        if float(cov) < 0.8 and geneName in dictRawTargetGene:
            logging.warning("%s coverage is %s < 0.8"%(geneName, cov))
            dictRawTargetGene.pop(geneName)
    f.close()
    for geneName in dictRawTargetGene:
        logging.warning("%s is covered by %s"%(geneName, chipName))

def CheckGene(args):
    logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(levelname)s %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename='%s.geneCheck.log'%args.oPrefix,
                filemode='w')
    #ReadGeneList
    dictTargetGene = GetTargetGene(args.genelist)
    geneInfoFile = os.path.dirname(os.path.abspath(sys.argv[0]))+"/config/Homo_sapiens.gene_info"
    CheckIfMainSymbol(dictTargetGene, geneInfoFile)
    CheckCoverage(dictTargetGene, args.chipName)
    return dictTargetGene

    #Check if the target gene all in main symbol

    #if not in main symbol check the get alias symbol
        #if not find logging
    #if the genelist of chip exists, check if all symbols in chip genelist
        #if not logging
    #get the bed msg for target gene
        # if gene not in bedfile of refgene
            #logging
        #else
            #record in final result


