# This script is run from the command line as follows:
# cat counts.txt | Rscript deseq2.r NxM > results.txt
# N is the number of samples in the first condition and M is the number of samples in the second condition.
# Condition one is the reference condition ("N").
# Produces a table with differentially expressed genes on the standard output.

# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Experimental design must be specified as: NxM at the command line", call.=FALSE)
}

first = args[1]

# Load DESeq2 library.
library(DESeq2)

# Extract the experimental design from the command line.
design = unlist(strsplit(first, 'x'))

# Find the desing counts.
cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])

# Set up the conditions based on the experimental setup.
cond_1 = rep("cond1", cond1_num)
cond_2 = rep("cond2", cond2_num)

# Read the data from the standard input.
countData = read.table("stdin", header=TRUE, sep="\t", row.names=1 )

# Build the dataframe from the conditions
samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

#Set the reference to be compared
dds$condition = relevel(dds$condition,"cond1")

# Run differential expression program.
dds = DESeq(dds)

# Format results.
res = results(dds)

# Sort the results by the padj and fold change.
sorted = res[with(res, order(padj, -log2FoldChange)), ]

# Create dataframe.
sorted.df = data.frame("id"=rownames(sorted),sorted)

# Write to txt file.
write.table(sorted.df, file="", sep="\t", col.names=NA, quote=FALSE)
