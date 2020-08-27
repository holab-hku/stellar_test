# Stellar - Test

Stellar test was written specifically to test the stellar software. It doesnot uses stellar independently but has the stellar logic inbuild to produce the required results. It takes in raw fastq files(R1 and R2) both human and mouse i.e. in total 4 fastq files and uses STAR aligner to:
  1. Run stellar on the individual mouse and human dataset input. These experiments are called experiment 3(mouse) and experiment 4(human).
  2. Use the two results to formulate chimeric data where mouse and human barcodes are mixed in proportion to mimic chimera dataset. These chimeric dataset are called experiment5.

On each dataset of exp3-5, stellar is run -> divide the raw fastq files into two individual human and mouse sets while finally aligning these inidividual sets using STAR. 


## Usage?
<ol>
  <li> <strong>Bring preliminary files from the storage directory to your personal account:</strong> <code>cp -r /storage/holab/datamart {your personalised path}</code></li>
  <li> <strong>Put <i>STARsolo</i> and <i>samtools</i> to path:</strong> It is recommended to fetch the latest version of these two softwares to get the enhanced experience but a copy of executables is present in the storage directory:
    <code>/storage/holab/datamart/STAR</code> & <code>/storage/holab/datamart/samtools-1.10</code><br/><br/>
    <i>To put these into path:</i><br/>
    - Go to the executable directory and then do <code>export PATH="$PATH:`pwd`"</code><br/>
    - Test if it was succesfully added to the PATH variable: <code>STAR --version</code>
   
  </li>
  
  <li> <strong>Create a directory and clone this repository:</strong> <code>git clone </code> </li>
  <li> <strong>Running the main script:</strong></li><br>
</ol>

  ```ruby
  perl pipeline.pl {R2 human fastq file} {R1 human fastq file} {R2 mouse fastq file} {R1 mouse fastq file} {path_to_chimera_index_files} {path_to_whitelist_barcode} {path_to_human_gtf} {path_to_mouse_gtf} {path_to_mouse_index_files} {path_to_human_index_files} {threshold} {output_directory} 

  Required:
    -{R2 human fastq file}: Path to the human R2 fastq file for 10x data
    -{R1 human fastq file}: Path to the human R1 fastq file for 10x data
    -{R2 mouse fastq file}: Path to the mouse R2 fastq file for 10x data
    -{R1 mouse fastq file}: Path to the mouse R1 fastq file for 10x data
    -{path_to_chimera_index_files}: Path to the chimera genome index - human + mouse mixed
    -{path_to_whitelist_barcode}: Path to the barcode whitelist file which is needed for 10X data. Depends on the chemistry(v2/v3) used to produce the 10X data
    -{path_to_human_gtf}: Path to human annotation gtf file
    -{path_to_mouse_gtf}: Path to mouse annotation gtf file
    -{path_to_mouse_index_files}: Path to mouse genome index file - only mouse
    -{path_to_human_index_files}: Path to human genome index file - only human
    -{threshold}: The threshold on the basis of which a particular barcode is classified as either mouse/human/unspecified. For eg if it is 80% then a barcode must have gene pct >=80% for either human/mouse for it to be classified otherwise it'll be regarded as 'unspecified' 
    -{output_directory}: The path of the directory where all the results of the current instance are put.
      
  ```
  
  ### Result Structure:
  This should be the tree view of the resulting directories/files after succesffuly running this script.
  
  ```bash
  .
  ├── output_dir: #
  │   ├── 5.1: #
  │   │   ├── analysis: #
  │   │   │   ├── 5.1_barcodes_classification.csv
  │   │   │   ├── 5.1.human_R1.fastq.gz
  │   │   │   ├── 5.1.human_R2.fastq.gz
  │   │   │   ├── 5.1.mouse_R1.fastq.gz
  │   │   │   ├── 5.1.mouse_R2.fastq.gz
  │   │   │   ├── Aligned.sortedByCoord.out.bam
  │   │   │   ├── barcodes.tsv.gz
  │   │   │   ├── features.tsv.gz
  │   │   │   ├── HumanFilteredAligned.out.sam.bam
  │   │   │   ├── matrix.csv
  │   │   │   ├── matrix.mtx.gz
  │   │   │   └── MouseFilteredAligned.out.sam.bam
  │   │   ├── Log.final.out
  │   │   ├── Log.out
  │   │   ├── Log.progress.out
  │   │   ├── SJ.out.tab
  │   │   └── Solo.out
  │   │       ├── Barcodes.stats
  │   │       └── Gene
  │   │           ├── Features.stats
  │   │           ├── filtered
  │   │           ├── raw
  │   │           │   ├── barcodes.tsv
  │   │           │   ├── features.tsv
  │   │           │   └── matrix.mtx
  │   │           ├── Summary.csv
  │   │           └── UMIperCellSorted.txt
  │   ├── 5.2
  │   │   ├── analysis
  │   │   │   ├── 5.2_barcodes_classification.csv
  │   │   │   ├── 5.2.human_R1.fastq.gz
  │   │   │   ├── 5.2.human_R2.fastq.gz
  │   │   │   ├── 5.2.mouse_R1.fastq.gz
  │   │   │   ├── 5.2.mouse_R2.fastq.gz
  │   │   │   ├── Aligned.sortedByCoord.out.bam
  │   │   │   ├── barcodes.tsv.gz
  │   │   │   ├── features.tsv.gz
  │   │   │   ├── HumanFilteredAligned.out.sam.bam
  │   │   │   ├── matrix.csv
  │   │   │   ├── matrix.mtx.gz
  │   │   │   └── MouseFilteredAligned.out.sam.bam
  │   │   ├── Log.final.out
  │   │   ├── Log.out
  │   │   ├── Log.progress.out
  │   │   ├── SJ.out.tab
  │   │   └── Solo.out
  │   │       ├── Barcodes.stats
  │   │       └── Gene
  │   │           ├── Features.stats
  │   │           ├── filtered
  │   │           ├── raw
  │   │           │   ├── barcodes.tsv
  │   │           │   ├── features.tsv
  │   │           │   └── matrix.mtx
  │   │           ├── Summary.csv
  │   │           └── UMIperCellSorted.txt
  │   ├── 5.3
  │   │   ├── analysis
  │   │   │   ├── 5.3_barcodes_classification.csv
  │   │   │   ├── 5.3.human_R1.fastq.gz
  │   │   │   ├── 5.3.human_R2.fastq.gz
  │   │   │   ├── 5.3.mouse_R1.fastq.gz
  │   │   │   ├── 5.3.mouse_R2.fastq.gz
  │   │   │   ├── Aligned.sortedByCoord.out.bam
  │   │   │   ├── barcodes.tsv.gz
  │   │   │   ├── features.tsv.gz
  │   │   │   ├── HumanFilteredAligned.out.sam.bam
  │   │   │   ├── matrix.csv
  │   │   │   ├── matrix.mtx.gz
  │   │   │   └── MouseFilteredAligned.out.sam.bam
  │   │   ├── Log.final.out
  │   │   ├── Log.out
  │   │   ├── Log.progress.out
  │   │   ├── SJ.out.tab
  │   │   └── Solo.out
  │   │       ├── Barcodes.stats
  │   │       └── Gene
  │   │           ├── Features.stats
  │   │           ├── filtered
  │   │           ├── raw
  │   │           │   ├── barcodes.tsv
  │   │           │   ├── features.tsv
  │   │           │   └── matrix.mtx
  │   │           ├── Summary.csv
  │   │           └── UMIperCellSorted.txt
  │   ├── exp5files: The experiment 5 data is formulated in this directory
  │   │   ├── 5.1.human._R1.fastq.gz: #Human R1 fastq produced from 5.1.human.sam.bam
  │   │   ├── 5.1.human._R2.fastq.gz: #Human R1 fastq produced from 5.1.human.sam.bam
  │   │   ├── 5.1.human.sam.bam: #The human BAM file containing 50% of human barcodes from original human data.
  │   │   ├── 5.1.mouse._R1.fastq.gz: #Mouse R1 fastq produced from 5.1.mouse.sam.bam
  │   │   ├── 5.1.mouse._R2.fastq.gz #Mouse R2 fastq produced from 5.1.mouse.sam.bam
  │   │   ├── 5.1.mouse.sam.bam: #The mouse BAM file containing 50% of mouse barcodes from original mouse data.
  │   │   ├── 5.1._R1.fastq.gz: #Concat 5.1.human._R1.fastq.gz & 5.1.mouse._R1.fastq.gz
  │   │   ├── 5.1._R2.fastq.gz #Concat 5.1.human._R2.fastq.gz & 5.1.mouse._R2.fastq.gz
  │   │   ├── 5.2.human._R1.fastq.gz
  │   │   ├── 5.2.human._R2.fastq.gz
  │   │   ├── 5.2.human.sam.bam
  │   │   ├── 5.2.mouse._R1.fastq.gz
  │   │   ├── 5.2.mouse._R2.fastq.gz
  │   │   ├── 5.2.mouse.sam.bam
  │   │   ├── 5.2._R1.fastq.gz
  │   │   ├── 5.2._R2.fastq.gz
  │   │   ├── 5.3.human._R1.fastq.gz
  │   │   ├── 5.3.human._R2.fastq.gz
  │   │   ├── 5.3.human.sam.bam
  │   │   ├── 5.3.mouse._R1.fastq.gz
  │   │   ├── 5.3.mouse._R2.fastq.gz
  │   │   ├── 5.3.mouse.sam.bam
  │   │   ├── 5.3._R1.fastq.gz
  │   │   └── 5.3._R2.fastq.gz
  │   ├── human
  │   │   ├── analysis
  │   │   │   ├── Aligned.sortedByCoord.out.bam
  │   │   │   ├── barcodes.tsv.gz
  │   │   │   ├── features.tsv.gz
  │   │   │   ├── human_barcodes_classification.csv
  │   │   │   ├── HumanFilteredAligned.out.sam.bam
  │   │   │   ├── human_R1.fastq.gz
  │   │   │   ├── human_R2.fastq.gz
  │   │   │   ├── matrix.csv
  │   │   │   └── matrix.mtx.gz
  │   │   ├── Log.final.out
  │   │   ├── Log.out
  │   │   ├── Log.progress.out
  │   │   ├── SJ.out.tab
  │   │   └── Solo.out
  │   │       ├── Barcodes.stats
  │   │       └── Gene
  │   │           ├── Features.stats
  │   │           ├── filtered
  │   │           ├── raw
  │   │           │   ├── barcodes.tsv
  │   │           │   ├── features.tsv
  │   │           │   └── matrix.mtx
  │   │           ├── Summary.csv
  │   │           └── UMIperCellSorted.txt
  │   ├── Log_STARAnalysis_5.1.txt
  │   ├── Log_STARAnalysis_5.2.txt
  │   ├── Log_STARAnalysis_5.3.txt
  │   ├── Log_STARAnalysis_Human.txt
  │   ├── Log_STARAnalysis_Mouse.txt
  │   └── mouse
  │       ├── analysis
  │       │   ├── Aligned.sortedByCoord.out.bam
  │       │   ├── barcodes.tsv.gz
  │       │   ├── features.tsv.gz
  │       │   ├── matrix.csv
  │       │   ├── matrix.mtx.gz
  │       │   ├── mouse_barcodes_classification.csv
  │       │   ├── MouseFilteredAligned.out.sam.bam
  │       │   ├── mouse_R1.fastq.gz
  │       │   └── mouse_R2.fastq.gz
  │       ├── Log.final.out
  │       ├── Log.out
  │       ├── Log.progress.out
  │       ├── SJ.out.tab
  │       └── Solo.out
  │           ├── Barcodes.stats
  │           └── Gene
  │               ├── Features.stats
  │               ├── filtered
  │               ├── raw
  │               │   ├── barcodes.tsv
  │               │   ├── features.tsv
  │               │   └── matrix.mtx
  │               ├── Summary.csv
  │               └── UMIperCellSorted.txt
  └── results
      ├── 3h
      ├── 3m
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 3mLog_STARAnalysis.txt
      ├── 4h
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 4hLog_STARAnalysis.txt
      ├── 4m
      ├── 5.1h
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 5.1.Log_STARAnalysis.txt
      ├── 5.1m
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 5.2h
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 5.2.Log_STARAnalysis.txt
      ├── 5.2m
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 5.3h
      │   ├── Aligned.sortedByCoord.out.bam
      │   ├── Log.final.out
      │   ├── Log.out
      │   ├── Log.progress.out
      │   ├── SJ.out.tab
      │   └── Solo.out
      │       ├── Barcodes.stats
      │       └── Gene
      │           ├── Features.stats
      │           ├── filtered
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── raw
      │           │   ├── barcodes.tsv
      │           │   ├── features.tsv
      │           │   └── matrix.mtx
      │           ├── Summary.csv
      │           └── UMIperCellSorted.txt
      ├── 5.3.Log_STARAnalysis.txt
      └── 5.3m
          ├── Aligned.sortedByCoord.out.bam
          ├── Log.final.out
          ├── Log.out
          ├── Log.progress.out
          ├── SJ.out.tab
          └── Solo.out
              ├── Barcodes.stats
              └── Gene
                  ├── Features.stats
                  ├── filtered
                  │   ├── barcodes.tsv
                  │   ├── features.tsv
                  │   └── matrix.mtx
                  ├── raw
                  │   ├── barcodes.tsv
                  │   ├── features.tsv
                  │   └── matrix.mtx
                  ├── Summary.csv
                  └── UMIperCellSorted.txt
  ```
  Summary: There are two principal directories being produced inside the root directory that the user specify. These are 'output_dir' and 'results'. The output_dir is where the initial STARSolo results on the raw input fastq, and formulation of exp5 data are put. The results dir hold the final results we want i.e. splitting an entire fastq dataset into human + mouse set and runninf STARSolo on it individually
  
  ### Example Run
  ```ruby
  nohup perl pipeline.pl ../fastq/exp_5.1_all_human_R2.fastq.gz ../fastq/exp_5.1_all_human_R1.fastq.gz ../fastq/exp_5.1_all_mouse_R2.fastq.gz ../fastq/exp_5.1_all_mouse_R1.fastq.gz /home/msnaveed/sra_local_repo/chimera_index/v3 /home/msnaveed/sra_local_repo/10x_genomics/3M-february-2018.txt /home/msnaveed/sra_local_repo/chimera_genome/human_genome/v3/*.gtf /home/msnaveed/sra_local_repo/chimera_genome/altered_mouse_genome/v3/*.gtf /home/msnaveed/sra_local_repo/chimera_genome/mouse_index/v3 /home/msnaveed/sra_local_repo/chimera_genome/human_index/v3 80 LEXP5 &
  ```
  #### Note
  - Make sure Rscript and Perl is also installed and is already on the PATH variable.
  - Use <code>nohup ... &</code> so that the command doesn't break upon logging off and continue to run in the backgroun even after signing off.
  - Within each script, code is commented for convenience
  - To learn on how to create genome index files visit the holab page[https://holab-hku.github.io/10X-workshop/processing-10x-rna-seq-data.html]

## File Hierarchy

### Stellar Files
  - **pipeline.pl**: The main script that contains the logic of running the whole stellar testing
  - **convert_to_genes_cells_matrices.r**: Converts dense matrix to sparse using R package Seurat. The output file is matrix.csv
  - **runSTARAnalysis.sh**: Runs STAR command and then runs convert_to_genes_cells_matrices.r on it to produce a matrix.csv
  - **classify_barcodes_pipeline.pl**: This script uses matrix.cvs to create a labelled file classifying each barcode in the raw input data as either human/mouse or unspecified
  - **convertBamtToFastq.pl**: Converts a bam file to fastq file. Need to run this individually for each R1 and R2 for 10X data
  - **runFinalSTAR.sh**: Runs a STAR command individually on each human and mouse dataset produced from original dataset
  
  <i>Note: To know more about the arguments of these files take, reads the header comments within each file.</i>
  
### DataMart Files
Present in the holab storage directory: <code>/storage/holab/datamart</code>
  - **gencode.v33.unique_gene_names.gtf**: Human gtf fie with unique gene names
  - **gencode.vM24.unique_gene_names.gtf**: Mouse gtf with unique gene names
  - **gencode.v33.vM24.concat.unique_gene_names.gtf**: Concatenated human + mouse gtf file with unique gene names
  - **human_index**: Indexed files generated using human genome
  - **mouse_index**: Indexed files generated using mouse genome
  - **chimera_index**: Indexed files generated using human + mouse genome. Concatenated gtf file and two fasta files were used
  - **3M-february-2018.txt**: Most updated barcode whitelist for v3 chemistry
  - **737K-august-2016.txt**: Barcode whitelist for v2 chemistry
  - **processed_gtf_files**: This directory contains the files for producing concatenated(human+mouse) gtf and raw human and mouse gtfs<br/>
      |<br/>
      |__ **create_gtf.pl**: Contains the logic to obtain individual and concatenated versions of human and mouse gtfs because of common gene names between mouse and human gtfs<br/>
      |<br/>
      |__  **gencode.v33.chr_patch_hapl_scaff.annotation.gtf**:  Original human gtf otained from the web<br/>
      |<br/>
      |__  **gencode.vM24.chr_patch_hapl_scaff.annotation.gtf**: Original mouse gtf otained from the web but m was appended to each line to mkae the chromosome name to mchr to differentiate from human<br/>
      
  - **STAR**: Executable directory of STAR version 2.7.3a
  - **samtools-1.10**: Executable directory of samtools version samtools 1.10
  
  <i>Note: It is highly recommended to get the latest version of samtools and starsolo and better to obtain the latest human and mouse annotation gtf files as well.</i>

### Sample Data for testing
A small sample data is put into the holab storage directory: <code>/storage/holab/stellar_test_sample</code>
  - **exp_5.1_all_human_R1.fastq.gz**: Human R1 fastq
  - **exp_5.1_all_human_R2.fastq.gz**: Human R2 fastq
  - **exp_5.1_all_mouse_R1.fastq.gz**: Mouse R1 fastq
  - **exp_5.1_all_mouse_R2.fastq.gz**: Mouse R2 fastq

## Time Complexity
The time it takes to run the script is directly related to the size of the fastq files fed. On running the script on human and mouse each R1(around 7GB) and R2(around 6GB) fastq files, it took around 2hr 5min to complete the whole process.


## Caution
This stellar test was performed on v2 chemistry 10X dataset. If there's any error thrown by STAR command, it is likely that some parameters might need to be changed within the command. These changes should be made in the STAR commands present in runSTARAnalysis.sh & runFinalSTAR.sh.
