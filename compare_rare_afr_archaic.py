"""
This script was written by Kelsey Witt in April of 2021. It takes outputs from the script "archaic_allele_freqs_pharmacogenees_mp.py"
and combined them to identify SNVs that are rare in Africa and shared with archaic humans. It requires two input files, "X_AFRrare_snp_freqs.csv"
and "archaic_SNPs_at_rare_pos_X.csv" where X is a pharmacogene, and makes an output file called "archaic_mp_nonafr_snp_freqs.csv". The output
is a table that contains the pharmacogene name, chromosome and position, the reference and alternate allele, archaic genotypes and archaic allele,
as well as the population frequency of the archaic allele for each of the non-African populations in the 1000 Genomes dataset
"""


pharmGenes = ["cyp1a2", "cyp2a6","cyp2b6","cyp2c19","cyp2c8","cyp2c9","cyp2d6","cyp2e1","cyp2j2","cyp3a4","cyp3a5"]
populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]

outfile = "archaic_mp_nonafr_snp_freqs.csv"
f = open(outfile, 'w')
headercols = ["pharmacogene","chromosome","position","rs value","ref","alt","rare allele","Den","Altai","Chagyrskaya","Vindija","arch_alt"] + populations
headerline = ",".join(headercols)+"\n"
f.write(headerline)
for pharmGene in pharmGenes:
    afrInfile = pharmGene + "_AFRrare_snp_freqs.csv"
    print(afrInfile)
    archInfile = "archaic_SNPs_at_rare_pos_" + pharmGene + ".csv"
    h = open(archInfile)
    with open(afrInfile) as g:
        h.readline()
        next(g)
        for afrLine in g:
            afrSharedArchaic = False
            archLine = h.readline()
            archSpline = archLine.split(sep=",")
            denGeno,altGeno,chagGeno,vinGeno = archSpline[3:7]
            if denGeno != "*" or altGeno != "*" or chagGeno != "*" or vinGeno != "*": #ensures there's data
                archPos = archSpline[2]
                afrSpline = afrLine.split(sep=",")
                afrPos = afrSpline[2]
                if archPos == afrPos:
                    nonAfrGeno = afrSpline[6]
                    for archaic in (denGeno,altGeno,chagGeno,vinGeno):
                        if nonAfrGeno in archaic:
                            afrSharedArchaic = True
                    if afrSharedArchaic:
                        afrStart=afrSpline[0:7]
                        afrEnd=afrSpline[7:]
                        archAlt=archSpline[8][:-1]
                        outLine = afrStart+[denGeno,altGeno,chagGeno,vinGeno,archAlt]+afrEnd
                        dataline = ",".join(outLine) +"\n"
                        f.write(dataline)   
f.close()