library("airway")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("stats")
library("DESeq2")

dir <- system.file("extdata", package="airway", mustWork = TRUE)
list.files(dir)

dir.exists("C:/Users/Ricardo/R/library/airway/extdata")

csvfile <- file.path(dir, "sample_table.csv")
(sampleTable <- read.csv(csvfile, row.names = 1))

#pasting sampleTable$Run_subset.bam name into filenames and check it
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)

#Provide R interface to BAM files and details about how the BAM files should be treated
bamfiles <- BamFileList(filenames, yieldSize = 2000000)

#See the sequence name in bamfiles
seqinfo(bamfiles[1])

gtffile <- file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- txdbmaker::makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))

#Granges object of the exons for a gene
(ebg <- exonsBy(txdb, by="gene"))

#register(SerialParam())
registered()
BiocParallel::register(BiocParallel::SnowParam(workers = 4), default = TRUE)
bp_result <- bplapply(1:4, function(x) {
  library(stats)
  return(sessionInfo())
}, BPPARAM = bpparam())
print(bp_result)

#Summarized Experiments
se <- summarizeOverlaps(features=ebg, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand = TRUE, fragment=TRUE)
dim(se)
#assay contains all the numerical data 
assayNames(se)
head(assay(se))
colSums(assay(se),3)
#data obtain from txdb and exonsBy showing where the genes is located in, seqnames and range + metadata of columns
rowRanges((se))
str(metadata(rowRanges(se)))

colData(se)
#because colData(se) is still empty, we can input the data from sampleTable
(colData(se) <- DataFrame(sampleTable))

se$cell
se$dex

se$dex <- factor(se$dex)

round(colSums(assay(se)),1)

dds <- DESeqDataSet(se, design = ~ cell + dex)
View(dds)
