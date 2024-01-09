# Tawfik_Project_Course
Contains script for project course

The data consists of treated and mock data from patients with SarsCov-2 and controll 

The preproccesing of the data that was preformed in Galaxy using the following tools and options: 

- FastQC Trimming (using Fastp) <br/>
- Alignment (using Hisat2) <br/>
- MutliQC Quantifications (using feature counts) <br/>
- differential expression analysis (edgeR) <br/>

The code in r contains the following:
- Reading in data <br/>
- Formatting data <br/>
- Creating a meta data file <br/>
-  Unsupervised clustering of the samples using dimensionality reduction <br/>
- Statistical test for differentially expressed genes (t-test) <br/>
- Adjusting p-values to correct for multiple testing. <br/>
- Calculating logFC <br/>
- Visualization and downstream analysis <br/>
- You need to use at least three different plots to visualize your data. <br/>
- Incorporate at least two approaches/analyses/plots you have learnt in <br/>
