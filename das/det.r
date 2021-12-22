#!/usr/bin/env Rscript

############################自动化差异分析脚本######################
### 制    作:liy                                                 ###
### 联系方式:liyin59375@gmail.com                                ###
####################################################################

library(argparser, quietly=TRUE)

p <- arg_parser("--差异分析并绘制火山图--")
p <- add_argument(p, "input", help="输入文件", type="character")
p <- add_argument(p, "-W", help="输出pdf,宽度设置", type="character",default="7.5")
p <- add_argument(p, "-H", help="输出pdf,高度设置", type="character",default="6")
p <- add_argument(p, "-F", help="log2FoldChange", type="character",default="1")
p <- add_argument(p, "-P", help="padj", type="character",default="0.05")
p <- add_argument(p, "-o", help="设置输出路径", type="character",default="./")

argv <- parse_args(p)


########
#### 1.加载包
library(ggplot2)
library(ggrepel)
suppressMessages(library(tidyverse))
options(warn=-1)
##############################################################


Dat<-read.csv(argv$input,header = T,row.names=1,stringsAsFactors = F,check.names = F,sep=',')
Dat$id <- rownames(Dat)
#确定是上调还是下调，用于给图中点上色）

FF <- argv$F%>%as.numeric()
PP <- argv$P%>%as.numeric()
FF1 = ifelse( argv$F>1,-2,-1)

Dat$threshold = factor(ifelse(Dat$padj < PP & abs(Dat$log2FoldChange) >= FF, ifelse(Dat$log2FoldChange>= FF ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Dat<-na.omit(Dat) 


picture <- ggplot(Dat,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p-adj)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名x称
  geom_vline(xintercept=c(FF1,FF),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(PP),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05P



#### 2.保存数据
wi<-argv$W%>%as.numeric()
he<-argv$H%>%as.numeric()
if(argv$o=='./'){
 dir.create("DET_result")
 setwd("DET_result")
}else{
 suppressMessages(dir.create(argv$o))
 setwd(argv$o)}


ggsave(picture,filename ="piciure.pdf",width =wi,height =he)

det.stat <- list()
det.stat$'---------cutoff-------' <- '' 
det.stat$log2FC<-FF
det.stat$padj<-PP
B <- Dat$threshold%>% table()%>%as.data.frame()
det.stat$'---------------------' <- '' 
det.stat$Up<- B[1,2]
det.stat$Down<- B[2,2]
det.stat$NoSignifi <- B[3,2]
det.stat$NoSignifi <- B[3,2]
det.stat <- t(t(det.stat))

gene <- filter(Dat,
               abs(log2FoldChange)>FF&padj<PP)%>%
  pull(id)%>%as.data.frame()
colnames(gene)<-'id'
Dat2 <- dplyr::select(Dat,id,log2FoldChange,pvalue,padj,threshold)

write.csv(gene,file='DET.csv',row.names = F)
write.table(det.stat,file='2.基因表达结果统计.txt',row.names = T,col.names = F)
write.csv(Dat2,file='1.基因表达结果.csv',row.names = F)
