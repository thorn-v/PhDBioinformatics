#===========================================================================
                                  ABOUT
  @author         :  Veronica
  @email          :  thornv2@mcmaster.ca
  @createdOn      :  Dec 4, 2024
  @updated        :  May 17, 2025
  description    :  

#===========================================================================

If not using samples from SRA, skip to 

# SRA

1. First, obtain whole genome sequencing data from the NCBI SRA by using the specific search query: `"Aspergillus fumigatus"[orgn] AND ("biomol dna"[Properties] AND "strategy wgs"[Properties] AND "filetype fastq"[Properties])`. This query filters the search results to include only whole genome sequencing (WGS), DNA datasets in FASTQ format (ex. add `"AND ("biomol dna"[Properties] AND "strategy wgs"[Properties] AND "filetype fastq"[Properties])` to the search).

2. Click **Send to** in the top right corner and select **Run Selector** > **Go**. On the newly opened page, download both the **Total Metadata** and **Accession List** from the middle, "select" box. Transfer that list to the server (using scp)

3. Using the accession list, we are going to download the raw sequencing files remotely from NCBI using the **sra-toolkit**. As there are many accessions we will use **GNU Parallels** in conjunction with a bash function download in parallel, taking advantage of a whole node.
run:
```
sbatch sraSetUp.sh
```
In the ecript, the `-j` option should correspond to the number in `--cpus-per-task=48`. `accs.txt` should be the downloaded list of accessions with no header. If the file has a header, delete it.

It does not finish in time, add add `--resume` before `--joblog` and re-run.
4. 

