#!/usr/bin/env Rscript


cat('-------欢迎使用富集分析脚本--------\n')

############################自动化富集分析脚本######################
### 制    作:liy                                                 ###
### 联系方式:liyin59375@gmail.com                                ###
####################################################################

suppressMessages(library(argparser,quietly=TRUE))
p <- arg_parser("自动分析脚本__liy")
p <- add_argument(p, "input", help="输入文件:差异分析结果", type="character")
p <- add_argument(p, "input2", help="输入文件:数据库文件(go/kegg)", type="character")
p <- add_argument(p, "-W", help="输出pdf,宽度设置", type="character",default="7.5")
p <- add_argument(p, "-H", help="输出pdf,高度设置", type="character",default="6")
p <- add_argument(p, "-o", help="设置输出路径", type="character",default="./")

argv <- parse_args(p)

### 加载R包
###############################################

suppressMessages(library(tidyverse,quietly=TRUE))
suppressMessages(library(clusterProfiler,quietly=TRUE))
options(stringsAsFactors = F)
suppressMessages(library(do,quietly=TRUE))
suppressMessages(library(lubridate))
options(warn=-1)
data1 <- lubridate::now()
data1 <-minute(as.POSIXct(data1))



####
cat('\n-----分析进度:10%\n')
####①此处可修改#################################

gene_list<-read.csv(argv$input)
colnames(gene_list)='gene'
gene2des<-read.csv(argv$input2,sep=',')
argv$s=ifelse('GO'%in%colnames(gene2des),'go','kegg')

if(argv$s=='go'){
  class <- dplyr::select(gene2des,GO,class)%>%unique()
  colnames(class) <- c('ID','class')
  argv$gene2_GID <- 'GID' ##描述文件基因列表头
  argv$gene2_term <- 'GO' ##描述文件term号表头/或KO
  argv$gene2_Name <- 'description' ##描述文件基因描述表头
  argv$result_wd <- 'GO富集分析' ##结果目录
  argv$result_file <- 'go_rich_result.csv' ##结果文件名
  argv$p_headline <- 'Top20 of GO enrichment'##图片标题
  argv$class<-'class'##go分类
}  else
{
  argv$gene2_GID <- 'GID' ##描述文件基因列表头
  argv$gene2_term <- 'KO' ##描述文件term号表头/或KO
  argv$gene2_Name <- 'description' ##描述文件基因描述表头
  argv$result_wd <- 'KEGG富集分析' ##结果目录
  argv$result_file <- 'kegg_rich_result.csv' ##结果文件名
  argv$p_headline <- 'Top20 of KEGG enrichment'##图片标题
}

####
cat('\n-----分析进度:30%\n')
#### 1.富集部分
rich <- enricher(gene = gene_list$gene,
                 TERM2GENE = gene2des[c(argv$gene2_term, argv$gene2_GID)],
                 TERM2NAME = gene2des[c(argv$gene2_term, argv$gene2_Name)],
                 pvalueCutoff = 1,###原值为0.05
                 pAdjustMethod = 'BH',
                 qvalueCutoff = 1,###原值为0.2
                 maxGSSize = 500)
###
cat('\n-----分析进度:50%\n')
#### 2.文件生成:输出显著富集结果 和 前20显著的结果
rich_results <- as.data.frame(rich)

if(argv$s=='go'){
  suppressMessages(rich_results<- left_join(rich_results,class))
}

if(argv$o!='./'){
 suppressMessages(dir.create(argv$o))
 setwd(argv$o)}
dir.create(argv$result_wd)
setwd(argv$result_wd)

rich_results1 <- rich_results[rich_results$pvalue<0.05 &rich_results$qvalue<0.2,]
rich_20 <- rich_results[1:20,]
richa <- rich_20[c("BgRatio","GeneRatio")]
richa <- Replace(data=richa,from='/.*$',to="")
a1 <- as.numeric(richa$BgRatio)
b1 <- as.numeric(richa$GeneRatio)
rich_20$Rich.factor <- a1/b1

write.csv(rich_results ,"all_result.csv", row.names = FALSE, quote = TRUE)
write.csv(rich_results1 , argv$result_file, row.names = FALSE,quote = TRUE)
write.csv(rich_20,"rich_20.csv",row.names = FALSE,quote = TRUE)


####
cat('\n-----分析进度:70%\n')
#### 3.画图部分

rich2.csv <- read.csv("rich_20.csv",encoding = "UTF-8")
labels=(levels(factor(rich2.csv$Description))[as.factor(rich2.csv$Description)])

rich2.csv$number <- factor(rev(1:nrow(rich2.csv)))
names(labels) = rev(1:20)
rich2.csv$shortname<-labels

if(argv$s=='go'){
  p <- ggplot(data=rich2.csv, aes(x=number, y=Rich.factor)) +
    geom_point(mapping = aes(size=Count,colour=-log10(qvalue),shape=class))+
    coord_flip() + theme_test() +
    scale_color_gradient(low = "darkgreen",high = "red")+
    scale_x_discrete(labels=labels) +
    labs(title = argv$p_headline,x=" ",y="Rich factor",
  colour="-log10(qvalue)",size="Gene number")+theme_bw()
}  else
{
  p <- ggplot(data=rich2.csv, aes(x=number, y=Rich.factor)) +
    geom_point(mapping = aes(size=Count,colour=-log10(qvalue)))+
    coord_flip() + theme_test() +
    scale_color_gradient(low = "darkgreen",high = "red")+
    scale_x_discrete(labels=labels) +
    labs(title = argv$p_headline,x=" ",y="Rich factor",
  colour="-log10(qvalue)",size="Gene number")+theme_bw()
  }



####
cat('\n-----分析进度:90%\n')
#### 4.统计部分
a <-data.frame(1)
a$name<-argv$s
a$pdf<-'.pdf'
a$stat<-'stat.txt'
a<-tidyr::unite(a,"pdf_result",name,pdf,remove=F)
a<-tidyr::unite(a,"stat_result",name,stat,remove=F)
picture <- p + guides(colour = guide_colorbar(order = 1),size = guide_legend(order = 2))
wi<-argv$W%>%as.numeric()
he<-argv$H%>%as.numeric()
ggsave(picture,filename = a$pdf_result,width =wi,height =he)


rich.stat <- list()
rich.stat$'----------------------' <- ''
rich.stat$det_gene_num<-dim(gene_list)[1]
rich.stat$total_num<- dim(rich_results)[1]
rich.stat$rich_num<- dim(rich_results1)[1]
rich.stat$pvalue<- '0.05'
rich.stat$qvalue<- '0.2'
if(argv$s=='go'){
 rich.stat$'---class_stat(rich)---' <- ''
 B <- rich_results1$class %>% table()%>%as.data.frame()
 rich.stat$BP <- B[1,2]
 rich.stat$CC <- B[2,2]
 rich.stat$MF <- B[3,2]
}
rich.stat <- t(t(rich.stat))
write.table(rich.stat,a$stat_result,row.names = T,col.names = F,sep = "\t",quote = F)


data2 <- lubridate::now()
data2 <-minute(as.POSIXct(data2))
data3<-data2-data1 +1
cat('\n-----分析进度:100%\n')
cat(paste("\n本次分析共计用时不到",data3,"min\n"))
cat(paste("\n-----恭喜您完成",argv$s,"富集分析,祝您生活愉快,٩(^ᴗ^)۶\n"))
