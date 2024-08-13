# Real_Allele
-Uses biopython version 1.81

- The Real_Allele Script is designed to determine allele frequencies and assign allele numbers from a consensus FASTA file containing DNA sequences.  
  
  
- Sequencing artifacts will be filtered out based on a criterion of zygosity proportion, and minimum read count.  

- Alleles are assigned based on their frequency, allele 1 will be assigned to the most common allele.  
- Creates a TSV dictionary file containing sequences and their counts.

  
- The FASTA file is needed to be named in this format: Gene_Exon.fasta (A minimum of 2 fields with underscores is required)
- The sequences within the file also need to be named to this convention:  
RHIMI-ARA13-xx-DN-020-JB_S86_17441  
"RHIMI-," "-xx-DN-020-JB", "-xx-VE-020-JB" will be removed.  
ARA13 is the sample name.  
17441 is the read count. 
This assumes three fields, all separated by "_," with the last field being the read count.  
The output name will have five elements: ARA13_AQP1_E0101_17441_A1.  
  
  
- Test files are included with the repository. First, run the script on an initial file:  
python allele_driver.py -f TEST_V1.fasta  
--a command will use the TSV dictionary file and apply the previous allele naming scheme to the new fasta file, keeping the alleles defined in the dictionary uniform:  
python allele_driver.py -f TEST_V2.fasta -a TEST_V1_alleles.tsv  
Help command python allele_driverV5.py -h   


- allele_driverV5.py is the version of the script used in this publication: https://www.sciencedirect.com/science/article/pii/S1877959X24000372
 
