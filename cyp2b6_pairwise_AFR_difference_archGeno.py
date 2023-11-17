"""
This script was written by Kelsey Witt in August 2021. The goal of this script is to calculate pairwise distances between all combinations of African
1000 genomes haplotypes and archaic haplotypes. Archaic genomes are pseudophased - if an individual is heterozygous the alleles are randomly assigned
to the two haplotypes except the Vindija Neanderthal, which is sorted so that the superarchaic variants are on chromosome "2" and the more Neanderthal-like
SNPs are on chromosome "1". The output is a table that specifies the IDs of the two individuals being compared, the type of comparison (within modern
humans, within archaic humans, and between modern and archaic humans), and the pairwise distance.
"""

import gzip
import random

geneSt = 41495203
archEnd = 41524278
geneEnd = 41524301
geneLen = geneEnd-geneSt

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

afrCols = []
afrNames = []
col_counter = 9
posList = []

with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] == "AFR": 
            afrCols.append(col_counter)
            afrNames.append(line_col[0])
        col_counter += 1
print(afrCols)
print(afrNames)
afrHaps = [[] for i in range(len(afrCols)*2)]

archNames = ["Denisovan", "Altai", "Chagyrskaya", "Vindija"]
archCols = [9,10,11,12]
archHaps = [[] for i in range(len(archNames))]

mInfile = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
aInfile = "./archaic_cyp2b6_merged.vcf.gz"

with gzip.open(mInfile, 'rt') as f:
    mLine = f.readline()
    g = gzip.open(aInfile, 'rt')
    aLine = g.readline()
    while '#' in mLine:
        mLine = f.readline()
    while '#' in aLine:
        aLine = g.readline()
    spline = mLine[:-1].split()
    mpos, _, ref, alt = spline[1:5]
    aSpline = aLine[:-1].split()
    apos, _, aRef, aAlt = aSpline[1:5]
    while int(mpos) < geneSt:
        mLine = f.readline()
        spline = mLine[:-1].split()
        mpos, _, ref, alt = spline[1:5]
    while int(mpos) <= archEnd:
        while int(apos) < int(mpos) and apos != archEnd: #Assign archaic genotypes as shown, assign refs for modern
            print(str(apos),"apos only")
            if len(aRef) == len(aAlt) == 1:
                for i in range(len(afrCols)):
                    afrHaps[i*2].append("0")
                    afrHaps[i*2+1].append("0")
                for i in range(len(archCols)):
                    aData = aSpline[archCols[i]].split(sep=":")
                    aGeno = aData[0]
                    archHaps[i].append(aGeno)
            aLine = g.readline()
            aSpline = aLine[:-1].split()
            apos, _, aRef, aAlt = aSpline[1:5]
        while int(mpos) < int(apos): #Assign refs for ancients, assign modern genotypes as shown
            if len(ref) == len(alt) == 1:
                for i in range(len(afrCols)):
                    afrCol = afrCols[i]
                    genotypes = spline[afrCol].split(sep="|")
                    afrHaps[i*2].append(genotypes[0])
                    afrHaps[i*2+1].append(genotypes[1])
                for i in range(len(archNames)):
                    archHaps[i].append("0/0")
            mLine = f.readline()
            spline = mLine[:-1].split()
            mpos, _, ref, alt = spline[1:5]
        if apos == mpos:
            if len(ref) == len(alt) == len(aAlt) == 1 and alt == aAlt:
                print(str(mpos))
                for i in range(len(afrCols)):
                    afrCol = afrCols[i]
                    genotypes = spline[afrCol].split(sep="|")
                    afrHaps[i*2].append(genotypes[0])
                    afrHaps[i*2+1].append(genotypes[1])
                for i in range(len(archCols)):
                    aData = aSpline[archCols[i]].split(sep=":")
                    print(aData)
                    aGeno = aData[0]
                    print(aGeno)
                    archHaps[i].append(aGeno)
            mLine = f.readline()
            aLine = g.readline()
            spline = mLine[:-1].split()
            mpos, _, ref, alt = spline[1:5]
            aSpline = aLine[:-1].split()
            apos, _, aRef, aAlt = aSpline[1:5]
    while int(mpos) <= geneEnd:
        if len(ref) == len(alt) == 1:
            for i in range(len(afrCols)):
                afrCol = afrCols[i]
                genotypes = spline[afrCol].split(sep="|")
                afrHaps[i*2].append(genotypes[0])
                afrHaps[i*2+1].append(genotypes[1])
            for i in range(len(archNames)):
                archHaps[i].append("0/0")
        mLine = f.readline()
        spline = mLine[:-1].split()
        mpos, _, ref, alt = spline[1:5]

outfile = "pairwise_AFR_haplotype_distances_geno.txt"
with open(outfile, 'w') as j:
    outCols = ["Ind1","Ind2","Comparison","PairwiseDiff"]
    outLine = "\t".join(outCols)+"\n"
    j.write(outLine)
    for k in range(0,len(afrHaps)-1):
        for l in range(k+1,len(afrHaps)):
            ind1Hap = afrHaps[k]
            ind1Label = afrNames[k//2] + "_" + str(k%2)
            ind2Hap = afrHaps[l]
            ind2Label = afrNames[l//2] + "_" + str(l%2)
            numMismatch = 0
            for x, y in zip(ind1Hap, ind2Hap):
                if y != x:
                    numMismatch += 1
            pairwiseDiff = numMismatch/geneLen
            outCols = [ind1Label, ind2Label, "Within Modern Humans", str(pairwiseDiff)]
            outLine = "\t".join(outCols)+"\n"
            j.write(outLine)
    for m in range(0,len(archHaps)):
        for n in range(0,len(afrHaps)):
            ind1Hap = archHaps[m]
            ind1Label = archNames[m]
            ind2Hap = afrHaps[n]
            ind2Label = afrNames[n//2] + "_" + str(n%2)
            numMismatch = 0
            for x, y in zip(ind1Hap, ind2Hap):
                numMismatch += 1-(x.count(str(y))/2)
            pairwiseDiff = numMismatch/geneLen
            outCols = [ind1Label, ind2Label, "Comparing Modern/Archaic Humans", str(pairwiseDiff)]
            outLine = "\t".join(outCols)+"\n"
            j.write(outLine)
    for p in range(0,len(archHaps)-1):
        for q in range(p+1,len(archHaps)):
            ind1Hap = archHaps[p]
            ind1Label = archNames[p]
            ind2Hap = archHaps[q]
            ind2Label = archNames[q]
            numMismatch = 0
            for x, y in zip(ind1Hap, ind2Hap):
                if y != x:
                    if y=="0/1" or x=="0/1":
                        numMismatch += 0.5
                    else:
                        numMismatch += 1
            pairwiseDiff = numMismatch/geneLen
            outCols = [ind1Label, ind2Label, "Within Archaic Humans", str(pairwiseDiff)]
            outLine = "\t".join(outCols)+"\n"
            j.write(outLine)