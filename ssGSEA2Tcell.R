# Define T cell infiltration strata using single sample GSEA (ssGSEA) results.
# Input requires Bindea T cells CD8TCELLS_BINDEA and CYTOTOXIC_TCELLS_BINDEA ssGSEA results in .gct format.
# Run this code from the command line as:
# cat ssGSEA.results.gct | Rscript ssGSEA2Tcell.R TumorType date
# TumorType is a string describing the samples you are working with (e.g., HNSC)
# Date is a string in any format without special characters or spaces (e.g., 210101)

# Load libraries ----------------------------------------------------------
library(data.table)
library(readr)
library(CePa)
library(tidyr)
library(tidyverse)

# Read the command line arguments.
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Arguments must be specified as: 'tumorType' 'date' at the command line minus the quotes", call.=FALSE)
}

# Load ssGSEA data from the standard input and create T cell variables ----------------
ssGSEAdata = read.gct("stdin")
ssGSEAdata.mat.scale <- scale(t(ssGSEAdata))

# Use this codeblock to make immune cell variables based on entire dataset --------------------------
ssGSEAdata.mat.scale.tbl <-  as.data.frame(ssGSEAdata.mat.scale) %>%
  rownames_to_column(var = "SampleID")  %>%
  as_tibble() %>%
  dplyr::mutate(Tstatus = as.factor((if_else((CD8TCELLS_BINDEA>quantile(CD8TCELLS_BINDEA,0.75) |
                                                CYTOTOXIC_TCELLS_BINDEA>quantile(CYTOTOXIC_TCELLS_BINDEA,0.75)),"Thi","Tlo")))
  ) %>%
  dplyr::select(SampleID,Tstatus) 

write_delim(ssGSEAdata.mat.scale.tbl,
            file=paste0(args[1],"_","Tcell_status_table","_",args[2],".txt"),
              delim = "\t",
              col_names = T)


