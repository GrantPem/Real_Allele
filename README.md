# allelinatorV5
-Uses biopython version 1.81  
-This script takes a FASTA file and reports the alleles   
-The script assumes that the input file will be named in this format  
AQ01_01.fasta (2 fields with underscores is important)  
-The script will sort and organizers via the frequency of an allele present in the fasta.  
-a command will use the TSV dictionary file an apply the previous allele naming scheme to the new fasta file keeping the allele's defined in the dictionary uniform.   
-The sequences within the file also need to be named to this convention:
 RHIMI-ARA13-xx-DN-020-JB_S86_17441 or 
"RHIMI-", "-xx-DN-020-JB", "-xx-VE-020-JB" will be removed
 ARA13 is the sample name
 17441 is the read count
The parameters choose the sample with the highest read count and any of the duplicate samples will be removed if the read count is below the parameter (e.g., multiple version of the same sample below 50% will be deleted.)
This assumes 3 fields, all separated by "_" with the last field being the read count

-Test files are included with the repository. First you would run the script on an initial file:
python allele_driver.py -f TEST_V1.fasta

-Then you would run the new file from the same locus, using the allele file produced in the first script
python allele_driver.py -f TEST_V2.fasta -a TEST_V1_alleles.tsv


Help command `python allele_driverV5.py -h`
