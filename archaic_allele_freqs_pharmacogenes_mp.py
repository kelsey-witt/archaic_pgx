"""
This script was written by Kelsey Witt in March of 2021. This script uses the specified pharmacogenes and their coordinates
and identifies single nucleotide variants for each pharmacogene that are at low frequency in Africa. The first output
file, pharmGene + "_AFRrare_snp_freqs.csv", outputs a table of positions that includes reference, alternate, and rare
alleles as well as frequencies in each non-African population.

The second half of the script looks at archaic genome files that have merged archaic human genotypes for each of the 
pharmacogenes. These files were created using bcftools. The second half of the script outputs all genotypes where
at least one archaic individual has an alternate allele, and the output file "archaic_SNPs_at_rare_pos_" + pharmGene + ".csv"
has a list of archaic genotypes and alternate alleles.
"""

import gzip
import re

pharmGenes = {"cyp1a2": ["15",75041183,75048941], "cyp2a6": ["19",41349442,41356352], "cyp2b6": ["19",41497203,41524301], \
    "cyp2c8": ["10",96796528,96829254], "cyp2c9": ["10",96698414,96749148], "cyp2c19": ["10",96522437,96612962], "cyp2d6": ["22",42522500,42526883], \
        "cyp2e1": ["10",135340866,135352620], "cyp2j2": ["1",60358979,60392470], "cyp3a4": ["7",99354582,99381811], \
            "cyp3a5": ["7",99245811,99277649]}

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

populations = ["AFR", "CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
afrPops = ["YRI", "LWK", "GWD", "MSL", "ESN"]
archaics = ["Denisovan","NAltai","NChagyrskaya","NVindija"]
Pop_Tracker = {}
rareAfrAlleles = []
archaicGenotypes = {}

for pop in populations:
    Pop_Tracker[pop] = []

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[2] == "AFR" and line_col[1] in afrPops: 
            Pop_Tracker["AFR"].append(col_counter)
        elif line_col[2] != "AFR":
            pop = line_col[1]
            Pop_Tracker[pop].append(col_counter)
        col_counter += 1

def decode_split_line(line):
    #line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

def snp_is_biallelic(line):
    ref_allele = line[3]
    alt_allele = line[4]
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return True
    else:
        return False

def parse_arch_line(line):
    position, _, ref_allele, arch_allele = line[1:5]
    if arch_allele == '.':
        arch_allele = ref_allele
    arch_info = line[9]
    return [position, ref_allele, arch_allele, arch_info]

def parse_mod_line(line):
    position, rsVal, ref_allele, mod_alt_allele = line[1:5]
    mod_info = line[7]
    return [position, ref_allele, mod_alt_allele, mod_info, rsVal]

def calc_afr_allele_freq(line):
    afr_allele_total = len(Pop_Tracker["AFR"]) * 2
    afr_allele_1_sum = 0
    for afr_ind in Pop_Tracker["AFR"]:
        afr_allele_1_sum += line[afr_ind].count("1")
    afr_allele_1_fq = float(afr_allele_1_sum/afr_allele_total)
    afr_allele_1_round = round(afr_allele_1_fq)
    return afr_allele_1_round

def identify_non_afr_alleles(line):
    non_afr_allele = "null"
    afr_alt_af = calc_afr_allele_freq(line)
    if afr_alt_af < 0.01 or afr_alt_af > 0.99:
        if afr_alt_af < 0.01: #allele of interest is alternate
            non_afr_allele = line[4]
        else: #allele of interest is reference
            non_afr_allele = line[3]
    return non_afr_allele

def arch_passes_quality_check(arch_info):
    arch_data_parse = arch_info.split(":")
    if len(arch_data_parse) == 11: #Altai Nea or Den info line
        genotype_quality = arch_data_parse[2]
        mapping_quality = arch_data_parse[9]
        if genotype_quality == ".":
            genotype_quality = 0.0
        if mapping_quality == ".":
            mapping_quality = 0.0
        if float(genotype_quality) >= 40.0 and float(mapping_quality) >= 30.0:
            return True
        else:
            return False
    elif len(arch_data_parse) == 8: #Vindija Nea info line
        genotype_quality = arch_data_parse[7]
        if float(genotype_quality) >= 40.0:
            return True
        else:
            return False

def calc_pop_freq(population, nonafr_allele):
    pop_allele_count = 0
    allele_total = len(Pop_Tracker[population]) * 2
    for ind in Pop_Tracker[population]:
        pop_allele_count += mod_line[ind].count(nonafr_allele)
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return rounded_freq

for pharmGene in pharmGenes:
    rareAfrAlleles = []
    archaicGenotypes = {}
    outfile = pharmGene + "_AFRrare_snp_freqs.csv"
    f = open(outfile, 'w')
    headercols = ["pharmacogene","chromosome","position","rs value","ref","alt","rare allele"] + populations
    headerline = ",".join(headercols)+"\n"
    f.write(headerline)
    chromosome, pharmSt, pharmEnd = pharmGenes[pharmGene][0:3]
    modern_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    with gzip.open(modern_file, 'rt') as mf:
        for mod_line in mf:
            mod_line = decode_split_line(mod_line)
            #print(mod_line)
            if '#' not in mod_line:
                mod_vars = parse_mod_line(mod_line)
                mod_position, mod_ref, mod_alt = mod_vars[0:3]
                mod_rs = mod_vars[4]
                if int(mod_position) >= pharmSt and int(mod_position) <= pharmEnd:
                    freq_by_pop = []
                    if snp_is_biallelic(mod_line):
                        mod_alt_allele = mod_line[4]
                        for pop in populations[1:]:
                                alleleCounter = 0
                                allele_total = len(Pop_Tracker[pop]) * 2
                                for ind in Pop_Tracker[pop]:
                                    alleleCounter += mod_line[ind].count("1")
                        non_afr_allele = identify_non_afr_alleles(mod_line)
                        if non_afr_allele != "null":
                            if non_afr_allele == mod_ref:
                                non_afr_genotype = "0"
                            elif non_afr_allele == mod_alt:
                                non_afr_genotype = "1"
                            rareAfrAlleles.append(int(mod_position))
                            archaicGenotypes[int(mod_position)] = ["*","*","*","*",mod_ref,mod_alt]
                            for pop in populations:
                                arch_freq = calc_pop_freq(pop, non_afr_genotype)
                                freq_by_pop.append(str(arch_freq))
                            datacols = [pharmGene,chromosome,mod_position,mod_rs,mod_ref,mod_alt,non_afr_genotype]+freq_by_pop
                            dataline = ",".join(datacols) +"\n"
                            f.write(dataline)
    f.close()
    for i in range(0,len(archaics)):
        archaic = archaics[i]
        archInfile = archaic + "_" + pharmGene + ".vcf.gz"
        archOutfile = "archaic_SNPs_at_rare_pos_" + pharmGene + ".csv"
        f = open(archOutfile, 'w')
        headercols = ["pharmacogene","chromosome","position","Den","Altai","Chag","Vind","ref","alts"]
        headerline = ",".join(headercols)+"\n"
        f.write(headerline)
        with gzip.open(archInfile, 'rt') as af:
            for archline in af:
                if '#' not in archline:
                    archline = decode_split_line(archline)
                    arch_vars = parse_arch_line(archline) 
                    arch_position, ref_allele, arch_alt_allele, arch_info = arch_vars[0:4]
                    if int(arch_position) in rareAfrAlleles:
                        if arch_alt_allele not in archaicGenotypes[int(arch_position)][5]:
                            archaicGenotypes[int(arch_position)][5] += "_" + arch_alt_allele
                        arch_genotype = arch_info[0:3]
                        archaicGenotypes[int(arch_position)][i] = arch_genotype
    for position in rareAfrAlleles:
        pos = int(position)
        datacols = [pharmGene,chromosome,str(pos)]+archaicGenotypes[pos]
        dataline = ",".join(datacols) +"\n"
        f.write(dataline)
    f.close()