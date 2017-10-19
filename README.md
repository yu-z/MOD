# MOD

## What is MOD?

* MOD stands for Mearsurement of Overlapping of DMRs. It is a pipeline to compare DNA methylation profiles among mutiple Bisulfite-Seq libraries. For details, see Ref. [1].

* The pipline includes several procedures:
   * BS-seq read mapping
   * DMR (Differentially methylated regions) calling by mutilple control libraries
   * S-MOD (Statistical Measurement of Overlapping of DMRs)
   * Q-MOD (Quantitative Mearsurement of Overlapping of DMRs)
   
* It is a very early version, we are working on the improvement of computational efficiency and making it more user freindly. 

## Run MOD

* It is a very early version, we are working on the improvement of computational efficiency and making it more user freindly. So you have to be a perl programmer to use this pipeline for current version.
The following paragraphs describe how to use the pipeline:

   * 1. BS-seq reads mapping (/1_BS-seq_mapping)
	This step processes raw fastq file from BS-seq. You have to install BSMAP firstly before run this step. 
	We have two scripts for mapping BS-seq reads using BSMAP:
	multiple_BSMAP.pl for single end libraries; multiple_BSMAP_PE.pl for pair end libraries. 
	The two scripts depend on the other three scripts:  methratio_alt.py, BSmap_to_cytosine.pl, Cytosine_to_100bp.pl; the two files also depend on two large file: TAIR10 genome, Cytosine files. So you have to check the route is right in above two scripts. 
	The step will generate 100bp files, seprated into CHH, CHG and CG files. You can use Parse_BSMAP_log_By100bp.pl to get the mapping statistcs result. 
		
   * 2. DMR calling (/2_DMR_calling)
	multiple_DMR_WTvsMUTANT.pl  processes a folder containing 100bp files from above step. It depends on the R script. You can use DMR_Allchr_commnad_noFtest.R if you have a lot of libraries to run. 
	This step will generate files containing DMRs from comaprision with test libraries and multiple libraries. 
	
   * 3. Mearsurement of Overlapping of DMRs (/DOM)
  This step proccess the the folder from second step. You can get S-MOD matrix using script S-MOD.pl and QMOD using Q-MOD.pl

## Feedbacks:

There is a comprehensive manual in the "doc" directory. You can e-mail (zhangy9@sustc.edu.cn) the author if you find errors in the manual or bugs in the source code, or have any suggestions/questions about the manual and code. Thank you!

## Citations

If you use MOD in your published work, we kindly ask you to cite the following paper which describes the central algorithms used in MOD:
* [1] xxx


