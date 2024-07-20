# FXS-carrier-identifier
Overview

FXS carrier identifier was developed by Qiwei Guo's laboratory at Xiamen University and can use for analyzing the data of fragile X syndrome carrier screening using nanopore sequence. FXS carrier identifier performs on the Ubuntu system, please install related softwares mentioned in Summary.R file in advance.


Getting Started
1. Download the FXS carrier identifier (STR_pip.py, SummaryplotSTR.R and SummaryplotSTR_motif.R) and its attached files. The attached files can be modified according to your needs.
2. Use Guppy to translate the raw signal into nucleotide sequence. The command is as follows: guppy_basecaller -i fast5_files -s output --config dna_r9.4.1_450bps_sup.cfg -r -x auto
3. Pass in the parameters by shell script: 
python STR_pip.py \
  --fastq fastq file  \
	-n sample name \
	-o output file \
  --qvalue 4 \
  --lowper 800 \
  --lowper1 0.5 \
  --gc_content 90 \
  --barseq tar60.fa \
  --tarfile tar.bed \
  --ifbar bar \
  --barcodemis 0 \
  --barcodefa barcode.fa \
  --STRmotif AGG \
	--SE SE.tsv \
5. Run the shell script: sh shell script
6. The FXS carrier identifier will generate a result file and display a part of results to the terminal.
