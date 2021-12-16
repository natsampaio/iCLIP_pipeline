# iCLIP_pipeline
MDA5 iCLIP pipeline - full sample processing, from fastq files to mapped and deduplicated reads.

Input requirements: 
- FATSQ file of Illumina-sequenced iCLIP libraries
- barcodes.txt: tab-delimited list of barcodes and respective sample codes; sample provided
- pipeline_iclip.yml: yaml file containing environment details (channels and dependencies), locations of required files, and editable variables for pipeline; sample provided
- Genome index: STAR index for desired genome for mapping; should be created with STAR 2.7 or above
