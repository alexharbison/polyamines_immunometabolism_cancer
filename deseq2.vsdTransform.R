# This program can be run as
# cat counts.txt | Rscript deseq2.vsdTransform.R arg
# Produces a table with variance stabilized counts on the standard output.

# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

# Load DESeq2 library.
library(DESeq2)

# Read the data from the standard input.
countData = read.table("stdin", header=TRUE, sep="\t", row.names=1 )

# Build the dataframe from the conditions
samples = names(countData)
colData = data.frame(samples=samples)

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~1)

# Generate variance stabilized normalized matrix and write to file
vsd <- vst(dds, blind=FALSE)
vsd.df <- assay(vsd)
vsd.df.deseq2 = data.frame("id"=rownames(vsd.df),vsd.df)
vsd.df.ssGSEA = data.frame("NAME"=rownames(vsd.df),"DESCRIPTION"="na",vsd.df)

write.table(vsd.df.deseq2, file=paste0(args,"_","vst_matrix_deseq2.txt"), sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

write.table(vsd.df.ssGSEA, file=paste0(args,"_","vst_matrix_ssGSEA_ready_deseq2.gct"), sep="\t", row.name=FALSE, col.names=TRUE,quote=FALSE)
