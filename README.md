# Real_Allele
-Uses biopython version 1.81

- The Real_Allele Script is designed for data sets as long as it formatted in the specifc ways listed below, Allele_Driver is the script used for the manuscript listed below.

-This script takes a consensus FASTA file for any gene and reports the alleles

  
-The script assumes that the input file will be named in this format


Gene_Exon.fasta (A minimum of 2 fields with underscores is required)
  
-The script will sort and organize by the frequency of an allele in the fasta.
  
--a command will use the TSV dictionary file and apply the previous allele naming scheme to the new fasta file, keeping the alleles defined in the dictionary uniform.  
  
-The sequences within the file also need to be named to this convention:  
RHIMI-ARA13-xx-DN-020-JB_S86_17441  
"RHIMI-," "-xx-DN-020-JB", "-xx-VE-020-JB" will be removed.  
ARA13 is the sample name  
17441 is the read count 
  
The parameters find the sample with the highest read count, and any duplicate samples will be removed if the read count is below the set parameter (e.g., multiple versions of the same sample with a read count below 50% will be deleted.) 
  
This assumes three fields, all separated by "_," with the last field being the read count.  
  
The output name will have five elements: ARA13_AQP1_E0101_17441_A1  
  
After the filtration step, the remaining samples are renamed to fit the new naming scheme of "Sample_Gene_Exon_Readcount_Allele."  
  
The dictionary tsv file will list all defined alleles. The allele with the most frequency will be at the top of the dictionary and named "Allele 1."  
  
-Test files are included with the repository. First, you would run the script on an initial file:  
python allele_driver.py -f TEST_V1.fasta  
-Then, you would run the new file from the same locus using the allele file produced in the first script:  
python allele_driver.py -f TEST_V2.fasta -a TEST_V1_alleles.tsv  
Help command python allele_driverV5.py -h  
