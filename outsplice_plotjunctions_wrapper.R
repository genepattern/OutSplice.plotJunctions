## The Regents of the University of California and The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2018) by the
## Regents of the University of California abd the 
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# Load any packages used to in our code to interface with GenePattern.
# Note the use of suppressMessages and suppressWarnings here.  The package
# loading process is often noisy on stderr, which will (by default) cause
# GenePattern to flag the job as failing even when nothing went wrong. 
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(TxDb.Hsapiens.UCSC.hg38.knownGene)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(Homo.sapiens)))
suppressMessages(suppressWarnings(library(OutSplice)))


# Print the sessionInfo so that there is a listing of loaded packages, 
# the current version of R, and other environmental information in our
# stdout file.  This can be useful for reproducibility, troubleshooting
# and comparing between runs.
args = commandArgs(trailingOnly=TRUE)
print("=============")
print(args)
sessionInfo()

is.emptyString=function(a){return (trimws(a)=="")}

# Get the command line arguments.  We'll process these with optparse.
# https://cran.r-project.org/web/packages/optparse/index.html
arguments <- commandArgs(trailingOnly=TRUE)

# Declare an option list for optparse to use in parsing the command line.
option_list <- list(
  # Note: it's not necessary for the names to match here, it's just a convention
  # to keep things consistent.
  make_option("--data.file", dest="data.file"),
  make_option("--number", dest="number"),
  make_option("--junctions", dest="junctions"),
  make_option("--tail", dest="tail", default=NULL),
  make_option("--p.value", dest="p.value", type="double"),
  make_option("--out.file.prefix", dest="out.file.prefix", default="OutputFile"),
  make_option("--gene", dest="gene", default='FALSE')  ,
  make_option("--symbol", dest="symbol", default=NULL),
  make_option("--tumor.color", dest="tumor.color", default="red"),
  make_option("--normal.color", dest="normal.color", default="blue")
)

# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=FALSE, args=arguments)
opts <- opt
data_file = opts$data.file

outdir = getwd()

if (!(endsWith(outdir, '/'))){
   outdir = paste(outdir, '/', sep="")
}

pdf <- paste0( opts$out.file.prefix, "", '.pdf')
pdf_output <- paste0('./', pdf)
print(data_file)
if (opts$tail == "NULL"){
    tail = NULL
} else {
    tail=opts$tail
}
print(paste("tail is ", tail))
PlotJunctionData(data_file, NUMBER=opts$number, junctions=opts$junctions, tail=tail, p_value = opts$p.value, GENE=opts$gene, SYMBOL=opts$symbol, makepdf=TRUE, pdffile = pdf_output, tumcol='red', normcol='blue')





