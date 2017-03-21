#! /usr/bin/env python
#coding=utf-8
import os, sys, re
import gzip



def GetTargetReg(dictTargetGene, bedFileForRefGene, oPrefix):
    dictTarReg = {}
    f = open(bedFileForRefGene, "r")
    for line in f:
        line = re.sub("chr", "", line)
        qTmp = re.split(r"\t", line.strip("\n"))
        if qTmp[-1] in dictTargetGene:
            if qTmp[0] not in dictTarReg:
                dictTarReg[qTmp[0]] = []
            dictTarReg[qTmp[0]].append([int(qTmp[1]), int(qTmp[2])])
    f.close();
    oF = open("%s.target.bed"%oPrefix, "w")
    for chr in dictTarReg:
        dictTarReg[chr] = sorted(dictTarReg[chr], key=lambda x:(x[0], x[1]))
        for i in dictTarReg[chr]:
            oF.write("%s\t%d\t%d\n"%(chr, i[0], i[1]))
    oF.close()
    return(dictTarReg)

def Align(fq1, fq2, seqmulePath, oPrefix):
    #/home/zhuying/grandbox/app/seqmule/exe/bwa/bwa mem -M  -T 0  -A 1  -B 4  -O 6  -E 1  -L 5  -U 17  -t 12 -R '@RG\tID:M16A0281_M16A0281\tSM:M16A0281\tPL:ILLUMINA\tLB:LIBRARY' /home/zhuying/grandbox/app/seqmule/bin/secondary/../../database/bwa/human_g1k_v37.fasta M16A0281_result/M16A0281.1.fastq.gz M16A0281_result/M16A0281.2.fastq.gz 
    #|  /home/zhuying/grandbox/app/seqmule/exe/samtools/samtools view -b -S  -F 256  -@ 12 - 
    #| /home/zhuying/grandbox/app/seqmule/exe/samtools/samtools sort  -@ 12 - M16A0281_result/M16A0281_bwamem.sort
    bwaPath = "%s/exe/bwa/bwa"%(seqmulePath)
    ref = "%s/database/bwa/human_g1k_v37.fasta"%(seqmulePath)
    samtools = "%s/exe/samtools/samtools"%(seqmulePath)
    rg = os.path.basename(oPrefix)
    command = "%s mem -M  -T 0  -A 1  -B 4  -O 6  -E 1  -L 5  -U 17  -t 12 %s %s %s"%(bwaPath, ref, fq1, fq2)
    command += "| %s view -b -S  -F 256  -@ 12 - "%(samtools)
    command += "| %s sort  -@ 12 - %s_bwamem.sort"%(samtools, oPrefix)
    #print(command)
    os.system(command)
    return("%s_bwamem.sort.bam"%(oPrefix))

def GetFastqFromBam(args, bamFile, dictTarReg, oPrefix, fq1, fq2):
    dictReadName = {};
    samtoolsPath = "%s/exe/samtools/samtools"%(args.Seqmule)
    #print("%s view %s"%(samtoolsPath, bamFile))
    #sys.exit(0)
    f = os.popen("%s view %s"%(samtoolsPath, bamFile))
    lastChr = "";
    lastIndexOfTag = 0;
    flank = 200
    for line in f:
        qTmp = re.split("\t", line.strip("\n"))
        #if last
        staPos = int(qTmp[3]);
        length = len(qTmp[9])
        if qTmp[2] != lastChr:
            lastChr = qTmp[2]; lastIndexOfTag = 0;
        if qTmp[2] not in dictTarReg:
            continue
        else:
            while(dictTarReg[qTmp[2]][lastIndexOfTag][1]+flank < staPos and lastIndexOfTag < len(dictTarReg[qTmp[2]]) - 1):
                lastIndexOfTag += 1;
            if staPos <= dictTarReg[qTmp[2]][lastIndexOfTag][1]+flank and staPos+length >= dictTarReg[qTmp[2]][lastIndexOfTag][0] -flank:
                dictReadName[qTmp[0]] = 1;

    #f.close();
    #print(dictReadName.keys())
    #exit(0)

    oFile1 = gzip.open("%s.R1.fastq.gz"%(args.oPrefix), "w");
    oFile2 = gzip.open("%s.R2.fastq.gz"%(args.oPrefix), "w");
    f1 = gzip.open(fq1, "r");
    f2 = gzip.open(fq2, "r")
    while (1):
        readName = f1.readline().strip("\n"); readName2 = f2.readline();
        fqseq1 = f1.readline(); fqseq2 = f2.readline();
        f1.readline(); f2.readline();
        fqQua1 = f1.readline(); fqQua2 = f2.readline();
        if not readName:
            break;
        sReadName = re.sub("\s+.*$", "", readName);
        sReadName = sReadName[1:]
        if sReadName in dictReadName:
            oFile1.write(readName+"\n");
            oFile1.write(fqseq1)
            oFile1.write("+\n")
            oFile1.write(fqQua1);
            oFile2.write(readName2);
            oFile2.write(fqseq2);
            oFile2.write("+\n");
            oFile2.write(fqQua2)
    f1.close();
    f2.close();
    oFile1.close();
    oFile2.close();

def Pipe(args, dictTargetGene):
    if not args.bam:
        bamFile = Align(args.fq1, args.fq2, args.Seqmule, args.oPrefix)
    else:
        bamFile = args.bam
    dictTarReg = GetTargetReg(dictTargetGene, args.refgene, args.oPrefix)
    #print(bamFile)
    GetFastqFromBam(args, bamFile, dictTarReg, args.oPrefix, args.fq1, args.fq2)

if __name__ == '__main__':
    #geneList, bamFile, fq1, fq2, oPrefix = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    #geneList, samDir = sys.argv[1], sys.argv[2]
    #samName = os.path.basename(samDir)
    #bamFile = samDir + "/"+samName+"_result/%s_bwamem.sort.bam"%(samName)
    #fq1 = "%s/%s_result/%s.1.fastq.gz"%(samDir, samName, samName)
    #fq2 = "%s/%s_result/%s.2.fastq.gz"%(samDir, samName, samName)
    geneList = "DM17A0067.gene.list"
    bamFile = "DM17A0067/LDM17A0029_result/LDM17A0029_bwamem.sort.bam"
    fq1 = "DM17A0067/LDM17A0029_result/LDM17A0029.1.fastq.gz"
    fq2 = "DM17A0067/LDM17A0029_result/LDM17A0029.2.fastq.gz"
    bedFileFromRefGene = "/home/zhangpeng/database/UniqCDS/refGene.hg19.uni.CDS";
    bedFileFromMedExome = "/home/chungong/bed/medexome/MedExome_hg19_capture_targets.bed"
    dictGeneRegion = GetGeneReg(bedFileFromRefGene, bedFileFromMedExome);
    #print dictGeneRegion;
    dictTarReg = GetTargetReg(geneList, dictGeneRegion);
    GetFastqFromBam(bamFile, dictTarReg, oPrefix, fq1, fq2);
