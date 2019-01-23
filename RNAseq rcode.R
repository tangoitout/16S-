msource("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

##################################in all###########################################################
b1<-read.table("b1.txt", header=F)
b2<-read.table("b2.txt", header=F)
b3<-read.table("b3.txt", header=F)
c1<-read.table("c1.txt", header=F)
c2<-read.table("c2.txt", header=F)
c3<-read.table("c3.txt", header=F)
d1<-read.table("d1.txt", header=F)
d2<-read.table("d2.txt", header=F)
d3<-read.table("d3.txt", header=F)
e1<-read.table("e1.txt", header=F)
e2<-read.table("e2.txt", header=F)
e3<-read.table("e3.txt", header=F)
f1<-read.table("f1.txt", header=F)
f2<-read.table("f2.txt", header=F)
f3<-read.table("f3.txt", header=F)

mergex=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2],d1[,2],d2[,2],e1[,2],e2[,2],e3[,2],f1[,2],f2[,2],f3[,2])
rownames(mergex) <- a1[,1]
countData=mergex
condition<-as.factor(c("DK","DK","DK","WT","WT","WT","WT+RbdG","WT+RbdG","NT","NT","NT","NTDrbdG","NTDrbdG","NTDrbdG"))
sampleTable1<-data.frame(condition=as.factor(condition))
colnames(sampleTable1)<-colnames(countData)
deseq <- DESeqDataSetFromMatrix(countData,DataFrame(condition), ~ condition)
d.deseq<-DESeq(deseq)
##PCA Plot###############################################################
vsdB<-varianceStabilizingTransformation(d.deseq)
jpeg(filename=paste('pca.jpg'),width=1024,height=728)
plotPCA(vsdB,intgroup=c("condition"))

##################heatmap#####################################################
vsdB_table<-as.data.frame(assay(vsdB))
vsdB_table_rowsum<-transform(vsdB_table,sum=rowSums(vsdB_table))
colnames(vsdB_table_rowsum)
selected<-order(vsdB_table_rowsum$sum,decreasing=TRUE)[1:100]
k<-as.matrix(vsdB_table[selected,])
colnames(k)<-c("DK","DK","DK","WT","WT","WT","WT+RbdG","WT+RbdG","NT","NT","NT","NTDrbdG","NTDrbdG","NTDrbdG")
####different clustering method)#######################3
Cluster_Method<-c( "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
for (i in 1:length(Cluster_Method)){
  #make a function to extract the cluster method
  myclust<-function(x){
    hclust(x,method=Cluster_Method[i])
  }
  #make heatmap by jpeg
  jpeg(filename=paste(Cluster_Method[i],'.jpg'),width=1024,height=728)
  heatmap.2(k,
            trace='none',
            dendrogram="column",
            scale="row",
            key=T,
            keysize=1.5,
            cexCol=0.9,
            hclustfun=myclust,
            labRow=NA,
            xlab='Strains',
            ylab='Genes',
            main=Cluster_Method[i],
            col=redgreen(75))
  dev.off()
}

heatmap.2(k,dendrogram="column",col=redgreen(75), scale="row",key=T,keysize=1.5,trace="none",cexCol=0.9)

desity.info="none",,labRow=NA
#############seperate senario################################################
###############################wt/dk############################################################################################
mergex1=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2],d1[,2],d2[,2],e1[,2],e2[,2],e3[,2],f1[,2],f2[,2],f3[,2])
rownames(mergex1) <- a1[,1]
countData1=mergex1
condition1<-as.factor(c("DK","DK","DK","WT","WT","WT"))
sampleTable1<-data.frame(condition1=as.factor(condition1))
colnames(sampleTable1)<-colnames(countData1)
deseq <- DESeqDataSetFromMatrix(countData,DataFrame(condition), ~ condition)
d.deseq<-DESeq(deseq)
condition1<-factor(c(rep("WT",3),rep("dk",3)))
#results#
res<-results(d.deseq,alpha=0.05)
res

#ordered by small p

resOrdered<-res[order(res$padj),]
head(resOrdered)

#adjusted q-value<0.05
sum(res$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig<-subset(resOrdered,padj<0.05)

head(regSig)

#export to excel
write.csv(as.data.frame(regSig),file="all.rapid_control_padj0.05.csv")



condition2<-factor(c(rep("WTP",3),rep("WTG",3)))
condition3<-factor(c(rep("NT",3),rep("NTDG",3)))
condition4<-factor(c(rep("WT",3),rep("NT",3)))
sampleTable<-data.frame(condition=as.factor(condition))
sampleTable
colnames(sampleTable)<-colnames(countData)
deseq <- DESeqDataSetFromMatrix(countData,DataFrame(condition), ~ condition)
deseq 
head(deseq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<3*3>
  countdata=read.table("xxxx.txt",sep=" ", header=TURE, row.names= null, check.rows = false, check.names = true)
condition<-factor(c(rep("wt",3),rep("dk",3),rep("slow",3),rep("control",3),rep("rapid",3),rep("slow",3),rep("control",3),rep("rapid",3),rep("slow",3)))
timepoints<-factor(c)rep("t1",8,rep("t2",8),rep("t3",8)))
sampleTable<-data.frame(condition=as.factor(condition),timepoints=as.factor(timepoints))
rownames(sampleTable)<-colnames(countdata)
sampleTable
deseq<-DESeqDataDetFromMatrix(countData=countdata,colData=sampleTable, design=~condition+timepoints)
--------------------------------------------------------------------------------------------------------------------

d.deseq<-DESeq(deseq)
##PCA Plot
vsdB<-varianceStabilizingTransformation(d.deseq)
plotPCA(vsdB,intgroup=c("condition","timepoints"))
plotPCA(vsdB,intgroup=c("timepoints"))
plotPCA(vsdB,intgroup=c("condition"))

#heatmap
library(gplots)
vsdB_table<-as.data.frame(assay(vsdB))
vsdB_table_rowsum<-transform(vsdB_table,sum=rowSums(vsdB_table))
colnames(vsdB_table_rowsum)
selected<-order(vsdB_table_rowsum$sum,decreasing=TRUE)[1:100]
vsdB_table[selected,]
heatmap.2(as.matrix(vsdB_table[selected,]),Rowv=F,dendrogram="column",col=redgreen(75), scale="row",key=T,keysize=1.5,desity.info="none",trace="none",cexCol=0.9,labRow=NA)
#results#
res<-results(d.deseq,alpha=0.05)
res

#ordered by small p

resOrdered<-res[order(res$padj),]
head(resOrdered)

#adjusted q-value<0.05
sum(res$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig<-subset(resOrdered,padj<0.05)

head(regSig)

#export to excel
write.csv(as.data.frame(regSig),file="all.rapid_control_padj0.05.csv")

get#ordered reslts by smallest P

t3.slow_controls=read.table("t3.slow_controls_padj0.05_overrepresentation.txt")
t3.slow=c(as.vector(t3.slow_control[,1]))
############################################################################################################################
# Download the data from github (click the "raw" button, save as a text file called "results.txt").
# https://gist.github.com/stephenturner/806e31fce55a8b7175af
res <- read.table("results.txt", header=TRUE)
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.5 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.5 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))



plot(abs(res$log2FoldChange),abs(res$log2FoldChange), col=res$Colour)

res$Colour="black"
# Set new column values to appropriate colours
res$Colour[abs(res$log2FoldChange)>=0.2]="red"
res$Colour[abs(res$log2FoldChange)>=0.3]="blue"
# Plot all points at once, using newly generated colours
plot(data$col_name1,data$col_name2, ylim=c(0,5), col=data$Colour, ylim=c(0,10))


# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))

###########################################################################################
all.timecourse=read.table("xxxx.txt")
all.timecourse=c(as.vector(all.timecourse[,1]))

grid.newpage()
venn.t2_t3_rapid_slow<-draw.qua.venn(
  area1=1374,
  area2=216,
  area3=4148,
  area4=3182,
  n12=198,
  n13=1155,
  n14=1056,
  n23=184,
  n24=183,
  n34=2910,
  n123=178,
  n134=1022,
  n234=176,
  n1234=172,
  category=("aaa","bbb","ccc","ddd"),
  fill=c("cornflower blue","pink","light green","mediumorchid"),
  fontface=rep("plain",15),
  fontfamily=rep("sans",15),
  cat.fontface=rep("plain",4),
  cat.fontfamily=rep("sans",4),
  scale=F,
  lty="blank",
  cex=2,
  cat.cex=2
)






)
