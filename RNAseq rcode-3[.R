library(DESeq2)
library(gplot)
library(calibrate)
annot=read.table("ann.txt",fill=TRUE)

ann<-read.csv("ann.csv",header=T)
rawdata<-read.csv("heatmap.csv",header=T)
dat <- merge(x=rawdata, y=ann, by.x = "geneID",by.y="geneID", all.x = TRUE)
dat2<-na.omit(dat)
data1<-dat[order(dat$PMID),]
cdata<-head(data1, 55)
head(cdata$V1)
##combine the two data set
vars <- c(cdata$V1:cdata$V1214])
cdata[vars] <- lapply((cdata[vars]+1), log2)

##############################
head(regSig)
finaldata<-as.matrix(regSig)
regSig[,1]
rawdata<-read.table("e1.txt", header=T)
ann<-read.csv("ann.csv",header=T)
x <- merge(x=rata, y=ann, by.x = "geneID",by.y="geneID", all.x = TRUE)
x1<-na.omit(x)
write.csv(as.data.frame(dat),file="dat.csv")
#############PCA##########Analysis######################


log transformation
head
##################################in all###########################################################
annotation
A1<-read.table("c1.txt", header=F)
A2<-read.table("c2.txt", header=F)
A3<-read.table("c3.txt", header=F)
A4<-read.table("e1.txt", header=F)
A5<-read.table("e2.txt", header=F)
A6<-read.table("e3.txt", header=F)
mergex=cbind(A1[,2],A2[,2],A3[,2],A4[,2],A5[,2],A6[,2])  
rownames(mergex) <- A1[,1]
write.csv(as.data.frame(mergex),file="mergex.csv")

countData<-read.csv("dat1.csv",header=T)
x2 <- merge(x=rawdatax, y=ann, by.x = "geneID",by.y="geneID", all.x = TRUE)
x1<-na.omit(x2)
write.csv(as.data.frame(x1),file="PP.csv")

x1<-read.csv("PP.csv",header=T)
countData=x1
condition<-as.factor(c("wt","wt","wt","nt","nt","nt"))
sampleTable1<-data.frame(condition=as.factor(condition))
rownames(sampleTable1)<-colnames(countData)
deseq <- DESeqDataSetFromMatrix(countData,DataFrame(condition), ~ condition)
d.deseq<-DESeq(deseq)
##PCA Plot###############################################################
vsdB<-varianceStabilizingTransformation(d.deseq)
#jpeg(filename=paste('pca.jpg'),width=1024,height=728)
plotPCA(vsdB,intgroup=c("condition"))

A7<-read.table("A7.txt", header=F)
A8<-read.table("A8.txt", header=F)
A9<-read.table("A9.txt", header=F)
A10<-read.table("A10.txt", header=F)
A11<-read.table("A11.txt", header=F)
A12<-read.table("A12.txt", header=F)
A13<-read.table("A13.txt", header=F)
A14<-read.table("A14.txt", header=F)
A15<-read.table("A15.txt", header=F)
c1<-read.table("A1.txt", header=F)
c2<-read.table("A2.txt", header=F)
c3<-read.table("A3.txt", header=F)

x2<-read.csv("dat1.csv",header=F)
mergex=cbind(x2[,2],x2[,3],x2[,4],x2[,5],x2[,6])  
rownames(mergex) <- x2[,1]
countData2=mergex
condition<-as.factor(c("cim_72", "cim_165", "cim_208","cim_220","cim_231"))
sampleTable1<-data.frame(condition=as.factor(condition))
colnames(sampleTable1)<-colnames(countData2)
deseq <- DESeqDataSetFromMatrix(countData2,DataFrame(condition), ~ condition)
d.deseq<-DESeq(deseq)

mergex=cbind(rawdata[,2],rawdata[,3],rawdata[,4],rawdata[,5])  
rownames(x1) <- x1[,1]

countData=mergex
condition<-as.factor(colnames(c("x1","x1","x2","x3","x3")))
sampleTable1<-data.frame(condition=as.factor(condition))
colnames(sampleTable1)<-colnames()
deseq <- DESeqDataSetFromMatrix(countData,DataFrame(condition), ~ condition)
condition<-colnames(x1)
deseq <- DESeqDataSetFromMatrix(x1,DataFrame(x1), ~ condition)

d.deseq<-DESeq(deseq)
##PCA Plot###############################################################
vsdB<-varianceStabilizingTransformation(d.deseq)
vsdB<-varianceStabilizingTransformation(x1)

#jpeg(filename=paste('pca.jpg'),width=1024,height=728)
plotPCA(x1,intgroup=c("condition"))

##################heatmap#####################################################




vsdB_table<-as.data.frame(assay(vsdB))
vsdB_table_rowsum<-transform(vsdB_table,sum=rowSums(vsdB_table))
colnames(vsdB_table_rowsum)
selected<-order(vsdB_table_rowsum$sum,decreasing=TRUE)[1:100]
k<-as.matrix(vsdB_table[selected,])
colnames(k)<-c("cim_72","cim_165","cim_208","cim_220","cim_231")


k<-as.matrix(countData2)
colnames(k)<-c("cim_72","cim_165","cim_208","cim_220","cim_231")
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

test<-heatmap.2(k,dendrogram="column",col=redgreen(75), scale="row",key=T,keysize=1.5,trace="none",cexCol=0.9)

y=k[rev(test$rowInd), test$colInd](get output)
#############seperate senario################################################


y<-c("Flagella","c-di-GMP","c-di-GMP","pili-TM","c-di-GMP","QS-hfq","pili-MSHA","c-di-GMP","Flagella","QS-hfq","
Flagella","Flagella","Flagella","Cys Pathway Rbd","c-di-GMP","Cys Pathway","Flagella","Flagella","c-di-GMP","
pili-pilA","Flagella","pili-MSHA","Rbd ","Flagella","c-di-GMP","Brp","Flagella","Flagella","Brp","
QS-hfq ","c-di-GMP","Flagella","mfp","c-di-GMP","Brp","c-di-GMP","Rbd","Flagella","Flagella","
c-di-GMP","pili-MSHA","Flagella","pili-pilA","Flagella","c-di-GMP","Rbd ","QS-hfq","Brp","Flagella","x")      


###############################wt/dk############################################################################################
mergex1=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2])
rownames(mergex1) <- b1[,1]
annot<-select(org.Hs.eg.db, symbols, "PMID", "ALIAS")
countData1=mergex1
condition1<-as.factor(c("DK","DK","DK","WT","WT","WT"))
sampleTable1<-data.frame(condition1=as.factor(condition1))
colnames(sampleTable1)<-colnames(countData1)
deseq1 <- DESeqDataSetFromMatrix(countData1,DataFrame(condition1), ~ condition1)
d.deseq1<-DESeq(deseq1)


#results#
res<-results(d.deseq,alpha=0.5)
res
#ordered by small p
resOrdered<-res[order(res$padj),]
head(resOrdered)
#adjusted q-value<0.05
sum(res$padj<0.5,na.rm=TRUE)
#Check result with a pass p value thred hold
regSig<-subset(resOrdered,padj<0.5)
head(regSig)
#export to excel
write.csv(as.data.frame(regSig),file="AB.csv")

#jpeg(filename=paste('wtdk.jpg'),width=1024,height=728)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log2(pvalue), pch=20, main="austin", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="green"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log2(pvalue), labs=b1[,1], cex=.8))


###############################wt/wt+g############################################################################################
mergex2=cbind(d1[,2],d2[,2],d3[,2],c1[,2],c2[,2],c3[,2])
rownames(mergex2) <- d1[,1]
countData2=mergex2
condition2<-as.factor(c("WT+G","WT+G","WT+G","WT","WT","WT"))
sampleTable2<-data.frame(condition2=as.factor(condition2))
colnames(sampleTable2)<-colnames(countData2)
deseq2<- DESeqDataSetFromMatrix(countData2,DataFrame(condition2), ~ condition2)
d.deseq2<-DESeq(deseq2)
#results#
res2<-results(d.deseq2,alpha=0.05)
res2

#ordered by small p

resOrdered2<-res2[order(res2$padj),]
head(resOrdered2)

#adjusted q-value<0.05
sum(res2$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig2<-subset(resOrdered2,padj<0.05)

head(regSig2)

#export to excel
write.csv(as.data.frame(regSig2),file="wt_WT+G.csv")

#jpeg(filename=paste('wt_WT+G.jpg'),width=1024,height=728)
# Make a basic volcano plot
with(res2, plot(log2FoldChange, -log2(pvalue), pch=20, main="WT_WT+G", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res2, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))
with(subset(res2, abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="orange"))
with(subset(res2, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="green"))
with(subset(res2, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log2(pvalue), labs=b1[,1], cex=.8))

###############################wt/dk############################################################################################
mergex1=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2])
rownames(mergex1) <- b1[,1]
countData1=mergex1
condition1<-as.factor(c("DK","DK","DK","WT","WT","WT"))
sampleTable1<-data.frame(condition1=as.factor(condition1))
colnames(sampleTable1)<-colnames(countData1)
deseq1 <- DESeqDataSetFromMatrix(countData1,DataFrame(condition1), ~ condition1)
d.deseq1<-DESeq(deseq1)
#results#
res1<-results(d.deseq1,alpha=0.05)
res1

#ordered by small p

resOrdered1<-res1[order(res1$padj),]
head(resOrdered1)

#adjusted q-value<0.05
sum(res1$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig1<-subset(resOrdered1,padj<0.05)

head(regSig1)

#export to excel
write.csv(as.data.frame(regSig1),file="wtdk.csv")

jpeg(filename=paste('wtdk.jpg'),width=1024,height=728)
# Make a basic volcano plot
with(res1, plot(log2FoldChange, -log2(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res1, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))
with(subset(res1, abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="orange"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="green"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log2(pvalue), labs=b1[,1], cex=.8))

###############################wt/dk############################################################################################
mergex1=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2])
rownames(mergex1) <- b1[,1]
countData1=mergex1
condition1<-as.factor(c("DK","DK","DK","WT","WT","WT"))
sampleTable1<-data.frame(condition1=as.factor(condition1))
colnames(sampleTable1)<-colnames(countData1)
deseq1 <- DESeqDataSetFromMatrix(countData1,DataFrame(condition1), ~ condition1)
d.deseq1<-DESeq(deseq1)
#results#
res1<-results(d.deseq1,alpha=0.05)
res1

#ordered by small p

resOrdered1<-res1[order(res1$padj),]
head(resOrdered1)

#adjusted q-value<0.05
sum(res1$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig1<-subset(resOrdered1,padj<0.05)

head(regSig1)

#export to excel
write.csv(as.data.frame(regSig1),file="wtdk.csv")

jpeg(filename=paste('wtdk.jpg'),width=1024,height=728)
# Make a basic volcano plot
with(res1, plot(log2FoldChange, -log2(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res1, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))
with(subset(res1, abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="orange"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="green"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log2(pvalue), labs=b1[,1], cex=.8))

###############################wt/dk############################################################################################
mergex1=cbind(b1[,2],b2[,2],b3[,2],c1[,2],c2[,2],c3[,2])
rownames(mergex1) <- b1[,1]
countData1=mergex1
condition1<-as.factor(c("DK","DK","DK","WT","WT","WT"))
sampleTable1<-data.frame(condition1=as.factor(condition1))
colnames(sampleTable1)<-colnames(countData1)
deseq1 <- DESeqDataSetFromMatrix(countData1,DataFrame(condition1), ~ condition1)
d.deseq1<-DESeq(deseq1)
#results#
res1<-results(d.deseq1,alpha=0.05)
res1

#ordered by small p

resOrdered1<-res1[order(res1$padj),]
head(resOrdered1)

#adjusted q-value<0.05
sum(res1$padj<0.05,na.rm=TRUE)

#Check result with a pass p value thred hold
regSig1<-subset(resOrdered1,padj<0.05)

head(regSig1)

#export to excel
write.csv(as.data.frame(regSig1),file="wtdk.csv")

jpeg(filename=paste('wtdk.jpg'),width=1024,height=728)
# Make a basic volcano plot
with(res1, plot(log2FoldChange, -log2(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res1, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))
with(subset(res1, abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="orange"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="green"))
with(subset(res1, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log2(pvalue), labs=b1[,1], cex=.8))

plot(abs(res1$log2FoldChange),abs(res1$log2FoldChange), col=res1$Colour)

res$Colour="black"
# Set new column values to appropriate colours
res$Colour[abs(res$log2FoldChange)>=0.2]="red"
res$Colour[abs(res$log2FoldChange)>=0.3]="blue"
# Plot all points at once, using newly generated colours
plot(data$col_name1,data$col_name2, ylim=c(0,5), col=data$Colour, ylim=c(0,10))


# Label points with the textxy function from the calibrate plot



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
",