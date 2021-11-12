#!/home/liy/anaconda3/envs/R/bin/Rscript
suppressMessages(library(argparser,quietly=TRUE))
p <- arg_parser("___afsamlon___by:liy(liyin59375@gmail.com)")
p <- add_argument(p, "input1", help="转录本对应文件", type="character")
p <- add_argument(p, "input2", help="样品信息表", type="character")
p <- add_argument(p, "-s", help="输入1或-1选择前/后-对照组/实验组,默认1,即前为实验组",type="character",default="1")
p <- add_argument(p, "-o", help="设置输出路径", type="character",default="./")

argv <- parse_args(p)


#### 1.

suppressMessages(library(DESeq2,quietly=TRUE))
suppressMessages(library(tidyverse,quietly=TRUE))
suppressMessages(library("tximport",quietly=TRUE))
suppressMessages(library("readr",quietly=TRUE))
options(warn=-1)


#### 2.

a <- argv$input1 ##读取数据
b <- argv$input2
tx2gene <- read.csv(a, sep = "\t")
samples <- read.csv(b, sep = "\t")
colnames(samples) <- c('sample','treatment')
files <- file.path('.', paste(substring(samples$sample, 1,4),".sf", sep = ""))
names(files) <- samples$sample ##得到'.sf文件路径

txi <- tximport(files, type="salmon", tx2gene=tx2gene)##tximport导入结果
dds <- DESeqDataSetFromTximport(txi, colData=samples, design= ~ treatment)##tximport导入的结果传递给DESeq2
dds <- DESeq(dds)##计算并输出


#### 3.
nn <- unique(samples$treatment)%>%as.data.frame()

if(argv$s==1)
{res <- results(dds, contrast = c("treatment",nn[1,1],nn[2,1]))
suppressMessages(sig<-print(paste(nn[1,1],"vs",nn[2,1])))
res$signal<-sig
}else
{res <- results(dds, contrast = c("treatment",nn[2,1],nn[1,1]))
suppressMessages(sig<-print(paste(nn[2,1],"vs",nn[1,1])))
res$signal<-sig
}
if(argv$o!='./'){
  suppressMessages(dir.create(argv$o))
  setwd(argv$o)}
write.table(res,"result.csv", sep = ",", row.names = TRUE,quote = T)
