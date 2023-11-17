# archaic_pgx
This repository contains the code needed to repeat the analyses used in "The impact of modern admixture on archaic human ancestry in human populations", published in GBE in 2023 by Wroblewski, Witt, and colleaguess. A link to the publication will be added once the manuscript is officially published.

The aim of this study was to examine pharmacogenetic variation in archaic humans and compare archaic and modern human pharmacogene variation to attempt to identify signals of archaic introgression. This project uses the 1000 Genomes Project Phase 3 data and the published high-coverage archaic genomes.

Some of the data in this project was analyzed by other co-authors and can be found at the [Claw Lab Github](https://github.com/the-claw-lab/aDNA_PGx_2021)

## Data Used
* [1000 Genomes Project Phase 3 VCF files](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/)
* [Denisovan and Altai, Chagyrskaya, and Vindija Neanderthal Genomes](https://www.eva.mpg.de/genetics/genome-projects/)


## Analyses and scripts
### Identifying archaic SNVs in pharmacogenes
archaic_allele_freqs_pharmacogenes_mp.py identifies alleles that are rare in Africans in 1000 genomes populations and alternate alleles present in archaic genotypes for pharmacogenes.
* Usage: python3 archaic_allele_freqs_pharmacogenes_mp.py
* Inputs: 1000 Genomes vcf and panel file, archaic VCF subsets representing each pharmacogene. Archaic files were created using bcftools to extract a region of interest from the Max Planck archaic genome VCFs and then combining them across all four archaic individuals. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 21), modern_file (vcf file, line 128), the output file for modern SNVs that are rare in 1000 Genomes African populations (outFile, line 122), the archaic merged vcf files (archInfile, line 163), and the output for pharmacogene haplotypes in archaic humans (archOutfile, line 164). A dictionary specifies which genomic regions are targeted (pharmGenes, line 16-19).
* Output: For each pharmacogene, two files are generated. The first output file, "pharmGene_AFRrare_snp_freqs.csv", outputs a table of positions that includes reference, alternate, and rare alleles as well as frequencies in each non-African population. The second output file, "archaic_SNPs_at_rare_pos_pharmGene.csv"
has a list of archaic genotypes and alternate alleles.

compare_rare_afr_archaic.py takes the two above outputs and combines them to produce a list of pharmacogene SNVs where the archaic allele is found in modern humans and rare in 1000 Genomes African populations.
* Usage: python3 compare_rare_afr_archaic.py
* Inputs: The two previous output files, "pharmGene_AFRrare_snp_freqs.csv" and "archaic_SNPs_at_rare_pos_pharmGene.csv", are specified as afrInfile (line 19) and archInfile (line 21) respectively. The output file is definied as "outfile" (line 13)
* Outputs: A csv file with the following columns: pharmacogene, chromosome, position, rsID, reference allele, alt allele, rare in Africa allele, Denisovan genotype, Chagyrskaya genotype, VIndija genotype, archaic allele, and the allele frequencies for each non-African 1000 Genomes population.

### Calculating pairwise distance between archaic and modern humans
cyp2b6_pairwise_AFR_difference_arch_geno.py calculates the pairwise distance between all combinations of African 1000 genomes individual haplotypes and archaic human haplotypes for CYP2B6, though it can be modified for other genes.
* Usage: python3 cyp2b6_pairwise_AFR_difference_arch_geno.py
* Inputs: The 1000 genomes VCF (mInfile, line 40) and an archaic merged genotype file for the region of interest (aInfile, line 41)
* Output: A tab-separated file with a line for each comparison that includes the individuals being compared, the type of comparison (within modern humans, within archaic humans, or between modern and archaic humans), and the pairwise distance.

For any questions, please contact Kelsey Witt at kwittdi@clemson.edu.
