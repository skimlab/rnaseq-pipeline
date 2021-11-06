
library(tximport) 
library(GenomicFeatures)
    
# tx2gene links transcript IDs to gene IDs for summarization
TxDb <- makeTxDbFromGFF(file = "/data1/CCSB/data/human/gtf/Homo_sapiens.GRCh38.91.chr.gtf")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
head(tx2gene)

#obtain the quantification file list
samples <- read.table(file.path("/data1/CCSB/data/human/quant/cdna","samples.txt"), header=TRUE)
files <- file.path("/data1/CCSB/data/human/quant/cdna", samples$run, "quant.sf")
  
#generate the gene-level summurization 
txi <- tximport(files, type="salmon", tx2gene=tx2gene)



