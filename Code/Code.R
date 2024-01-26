library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggisoband)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(ggnewscale)
library(WGCNA)
library(future)

##########################################################################################
##########################################################################################
#############################################Load expression matrix and cell-level meta data
#############################################"ExpressionMatrix.csv" and "CellMetaData.csv" are in folder namely "materials"
matrix <- fread("ExpressionMatrix.csv",sep=",",stringsAsFactors=F,header=T)
gene_name=as.vector(as.matrix(matrix[,1]))
matrix <- matrix[,-1]
matrix <- as.matrix(matrix)
row.names(matrix)=gene_name

metadata <- fread("CellMetaData.csv",sep=",",stringsAsFactors=F,header=T)
metadata = data.frame(metadata,row.names = as.vector(as.matrix(metadata[,1])))
metadata=metadata[,-1]
row.names(metadata)=colnames(matrix)

Merge <- CreateSeuratObject(counts = matrix, min.cells = -1, min.features = -1, meta.data=metadata, project = "10X")



##########################################################################################
##########################################################################################
#############################################Plotting UMAP of control donor
#############################################"UMAP_Control.csv" is in folder namely "materials"
mixc=subset(Merge, subset = Disease == 'Control')

mixc <- NormalizeData(mixc, normalization.method = "LogNormalize", scale.factor = 10000)
mixc <- FindVariableFeatures(mixc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mixc)
mixc <- ScaleData(mixc, features = all.genes)
mixc <- RunPCA(mixc, features = VariableFeatures(object = mixc))
mixc <- FindNeighbors(mixc, dims = 1:20)
mixc <- RunUMAP(mixc, dims = 1:20)
UMAP_control=fread("UMAP_Control.csv",sep=",",stringsAsFactors=F,header=T)
UMAP_control=UMAP_control[,-1]
UMAP_control=as.matrix(UMAP_control)
row.names(UMAP_control)=colnames(mixc)
mixc@reductions$umap@cell.embeddings=UMAP_control

pdf('1b.pdf')
DimPlot(mixc, reduction = "umap",label = T,group.by = 'CellType')
dev.off()
pdf('1e.pdf')
DimPlot(mixc, reduction = "umap",label = F,group.by = 'Region')
dev.off()



##########################################################################################
##########################################################################################
#############################################Heatmap with scaling by column showed the relative number of predefined cell type markers overlapped with top 100 genes ranked by PEM scores
#############################################"bglist.csv" is in folder namely "materials"
cluster_v=c()
for(i in 1:dim(mixc)[2]){
if(mixc@meta.data[i,'SubCellType'] == 'Oligo-1'){cluster_v=c(cluster_v,0)}
if(mixc@meta.data[i,'SubCellType'] == 'D2-MSN-1'){cluster_v=c(cluster_v,1)}
if(mixc@meta.data[i,'SubCellType'] == 'D1-MSN'){cluster_v=c(cluster_v,2)}
if(mixc@meta.data[i,'SubCellType'] == 'Oligo-2'){cluster_v=c(cluster_v,3)}
if(mixc@meta.data[i,'SubCellType'] == 'Microglia'){cluster_v=c(cluster_v,4)}
if(mixc@meta.data[i,'SubCellType'] == 'Astro-1'){cluster_v=c(cluster_v,5)}
if(mixc@meta.data[i,'SubCellType'] == 'OPC'){cluster_v=c(cluster_v,6)}
if(mixc@meta.data[i,'SubCellType'] == 'Endothelial-1'){cluster_v=c(cluster_v,7)}
if(mixc@meta.data[i,'SubCellType'] == 'Astro-2'){cluster_v=c(cluster_v,8)}
if(mixc@meta.data[i,'SubCellType'] == 'Endothelial-2'){cluster_v=c(cluster_v,9)}
if(mixc@meta.data[i,'SubCellType'] == 'InterN-1'){cluster_v=c(cluster_v,10)}
if(mixc@meta.data[i,'SubCellType'] == 'Da'){cluster_v=c(cluster_v,11)}
if(mixc@meta.data[i,'SubCellType'] == 'D1/D2 Mixed-1'){cluster_v=c(cluster_v,12)}
if(mixc@meta.data[i,'SubCellType'] == 'D1/D2 Mixed-2'){cluster_v=c(cluster_v,13)}
if(mixc@meta.data[i,'SubCellType'] == 'Endothelial-3'){cluster_v=c(cluster_v,14)}
if(mixc@meta.data[i,'SubCellType'] == 'InterN-2'){cluster_v=c(cluster_v,15)}
if(mixc@meta.data[i,'SubCellType'] == 'InterN-3'){cluster_v=c(cluster_v,16)}
if(mixc@meta.data[i,'SubCellType'] == 'D2-MSN-2'){cluster_v=c(cluster_v,17)}
if(mixc@meta.data[i,'SubCellType'] == 'InterN-4'){cluster_v=c(cluster_v,18)}
if(mixc@meta.data[i,'SubCellType'] == 'Ependy'){cluster_v=c(cluster_v,19)}
}
mixc <- AddMetaData(object = mixc,metadata = factor(as.vector(cluster_v),levels=c(0:19)),col.name = "cluster")

calculatePEM <- function(SeuratObject,GroupName,subGroupName){
CP_ClusterNum = length(subGroupName)
CP_cluster_mean <- matrix(nr=dim(SeuratObject)[1],nc=CP_ClusterNum)
for(i in 1:CP_ClusterNum){
CP_tmp <- FetchData(object = SeuratObject, vars = row.names(SeuratObject) , cells = colnames(SeuratObject)[which(SeuratObject@meta.data[,GroupName] == subGroupName[i])])
CP_tmp = as.matrix(CP_tmp)
CP_tmp <- colMeans(CP_tmp)
CP_tmp <- t(CP_tmp)
CP_tmp <- matrix(CP_tmp)
row.names(CP_tmp) <- row.names(SeuratObject)
CP_cluster_mean[,i] <- CP_tmp
}
row.names(CP_cluster_mean) <- row.names(SeuratObject)
colnames(CP_cluster_mean) <- subGroupName
CP_cluster_sum <- colSums(CP_cluster_mean)
CP_psedo <- mean(CP_cluster_sum/dim(SeuratObject)[1])
CP_cluster_mean <- CP_cluster_mean + 0.1*CP_psedo
CP_cluster_sum <- colSums(CP_cluster_mean)
CP_sum <- sum(CP_cluster_sum)
CP_PEM <- matrix(nr=dim(SeuratObject)[1],nc=CP_ClusterNum)
row.names(CP_PEM) <- row.names(SeuratObject)
colnames(CP_PEM) <- subGroupName
for(i in 1:dim(SeuratObject)[1]){
for(j in 1:CP_ClusterNum){
CP_PEM[i,j] = log10((CP_sum*CP_cluster_mean[i,j])/(CP_cluster_sum[j]*sum(CP_cluster_mean[i,])))
}}
CP_PEM[is.na(CP_PEM)] <- -Inf
return(CP_PEM)
}
PEM = calculatePEM(mixc,'cluster',c(0:19))

bglist <- fread("bglist.csv",sep=",",stringsAsFactors=F,header=T)
bglist <- as.matrix(bglist)
num10=100
num_cluster=20
m6 <- as.character(matrix(c(0:(num_cluster-1)),nrow=1))
m7 <- paste('cluster',m6,sep='')
PEMsort <- matrix(nr=num10,nc=3*num_cluster)
colnames(PEMsort) <- as.character(1:(3*num_cluster))
for(i in 1:num_cluster){
a <- data.frame(PEM[,i])
b <- row.names(a)[order(a[,1],decreasing=T)]
a<- a[order(a[,1],decreasing=T),1]
a <- matrix(a)
row.names(a) <- b
PEMsort[,3*i-2] <- row.names(a)[1:num10]
PEMsort[,3*i-1] <- a[1:num10,1]
colnames(PEMsort)[3*i-2] <- paste('cluster',m6[i],seq='')
colnames(PEMsort)[3*i-1] <- 'PEM score'
colnames(PEMsort)[3*i] <- 'cell type'
for(j in 1:num10){
if(length(which(bglist[,1]==PEMsort[j,3*i-2])) != 0) { 
  PEMsort[j,3*i] = bglist[which(bglist[,1]==PEMsort[j,3*i-2])[1],2]
   }
  else { PEMsort[j,3*i] = '' }
}
}
m33 <- list(c('^Neuron','^Astro','^Oligo','^Vascular','^Microglia','^OPC','^Macrophage','^Ependy','^Endothelial','^Mural'),c('Neuron','Astro','Oligo','Vascular','Microglia','OPC','Macrophage','Ependy','Endothelial','Mural'))
m44 <- matrix(nr=10,nc=num_cluster)
row.names(m44) <- m33[[2]]
colnames(m44) <- m7
for(i in 1:10){
for(j in 1:num_cluster){
gene0 = PEMsort[,3*j]
gene0 = gene0[grepl(m33[[1]][i],gene0)]
if(dim(matrix(gene0))[1] > 0){
m44[i,j] = dim(matrix(gene0))[1]
}
else {m44[i,j] = 0}
}
}
pdf('1c.pdf')
pheatmap(m44,scale="column",cluster_cols=F,cluster_rows=F,angle_col = c("315"),color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()



##########################################################################################
##########################################################################################
#############################################Markers and expression statistics used for defining cell types in control basal ganglia
m7=c("D1-MSN","D2-MSN-1","D2-MSN-2","D1/D2 Mixed-1","D1/D2 Mixed-2","InterN-1","InterN-2","InterN-3","InterN-4","Da","Oligo-1","Oligo-2","Astro-1","Astro-2","Microglia","Endothelial-1","Endothelial-2","Endothelial-3","OPC","Ependy")
gene_list = c('SYT1','MYT1L','GRIA1','CADPS','GABRB3','GAD1','GAD2','TH','SLC6A3','PPP1R1B','DRD1','DRD2','ELAVL2','CLSTN2','PLP1','PLEKHH1','AQP4','ACSBG1','APBB1IP','IKZF1','LEF1','ITIH5','PCDH15','FAM19A1','SPAG17','RGS22')
num_cluster=length(m7)
nums=length(gene_list)
ave_exp=matrix(nr=nums,nc=num_cluster)
ave_per=matrix(nr=nums,nc=num_cluster)
row.names(ave_exp) = gene_list
colnames(ave_exp) = m7
row.names(ave_per) = gene_list
colnames(ave_per) = m7
ge = gene_list
for(j in 1:length(ge)){
for(i in 1:num_cluster){
cells = colnames(mixc)[which(mixc@meta.data[,'SubCellType']==m7[i])]
ave_exp[j,i] = mean(mixc[['RNA']]@data[ge[j],cells])
ave_per[j,i]=sum(mixc[['RNA']]@data[ge[j],cells]>0)/length(cells)
}
ave_exp[j,] = t(scale(ave_exp[j,],center=T,scale=T))
}
write.csv(ave_exp,file='ave_exp.csv')
write.csv(ave_per,file='ave_per.csv')
dat = cbind(melt(ave_exp),melt(ave_per)[,3])
colnames(dat) = c('Y','X','exp','per')
dat[,1] = factor(dat[,1],level = rev(gene_list))
dat[,2] = factor(dat[,2],level = m7)
p <- ggplot(dat,aes(x = X, y = Y)) +
geom_point(aes(size = per, fill = exp), shape = 21 ) +
scale_fill_gradient(low = "#FFFFFF", high = "#FF0000")+
xlab("") + ylab("") +  labs(size = "Expression Ratio",fill = 'Z-Score') +
theme_bw() +
theme(panel.grid=element_line(colour=NA),axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=9),axis.text.y = element_text(color = "black",size=9))
pdf('1d.pdf')
print(p)
dev.off()



##########################################################################################
##########################################################################################
#############################################Cell ratios of cell type composition per region or regional composition per cell type
regions <- c('caudate','globus pallidus','putamen','substantia nigra')
types <- c("D1-MSN","D2-MSN-1","D2-MSN-2","D1/D2 Mixed-1","D1/D2 Mixed-2","InterN-1","InterN-2","InterN-3","InterN-4","Da","Oligo-1","Oligo-2","Astro-1","Astro-2","Microglia","Endothelial-1","Endothelial-2","Endothelial-3","OPC","Ependy")
cl_di  <- matrix(nr=length(types),nc=4)
row.names(cl_di) <- types
colnames(cl_di) <- regions
for(i in 1:length(types)){
for(j in 1:4){
tee <- mixc@meta.data[,c('Region','SubCellType')]
tee <- tee[which(tee[,1]==regions[j]),]
tee <- tee[which(tee[,2]==types[i]),]
cl_di[i,j] <- dim(tee)[1]
}}
cs <- colSums(cl_di)
tee2 <- cl_di
for(i in 1:4){
tee2[,i] <- tee2[,i]/cs[i]
}

cll <- tee2[,1]
cll2 <- names(cll)
cll3 <- rep(colnames(tee2)[1],dim(tee2)[1])
for(i in 2:dim(tee2)[2]){
clll <- tee2[,i]
clll2 <- names(clll)
clll3 <- rep(colnames(tee2)[i],dim(tee2)[1])
cll <- c(cll,clll)
cll2 <- c(cll2,clll2)
cll3 <- c(cll3,clll3)
}

pdf('2a.pdf')
ggplot(data=NULL, aes(x=cll3, y=cll)) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(fill = factor(cll2,level = row.names(cl_di)))) +
theme(axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=9)) +
labs(x ='',y='Cell ratio') +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(limits=colnames(cl_di)) +
labs(fill = "CellType")
dev.off()

tee2 = cl_di/rowSums(cl_di)
cll <- tee2[,1]
cll2 <- names(cll)
cll3 <- rep(colnames(tee2)[1],dim(tee2)[1])
for(i in 2:dim(tee2)[2]){
clll <- tee2[,i]
clll2 <- names(clll)
clll3 <- rep(colnames(tee2)[i],dim(tee2)[1])
cll <- c(cll,clll)
cll2 <- c(cll2,clll2)
cll3 <- c(cll3,clll3)
}

pdf('2b.pdf')
ggplot(data=NULL, aes(x=cll2, y=cll)) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(fill = cll3)) +
theme(axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=9)) +
labs(x ='',y='Cell ratio') +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(limits=row.names(cl_di)) +
labs(fill = "Region")
dev.off()



##########################################################################################
##########################################################################################
#############################################Heatmap showed cell type- and region-specific gene expression signatures
regions <- c('caudate','globus pallidus','putamen','substantia nigra')
types <- c("D1-MSN","D2-MSN","D1/D2 Mixed","Interneuron","Da","Oligo","Astro","Microglia","Endothelial","OPC","Ependy")
iiiii=1
store = matrix(nc=5,nr=30000)
colnames(store) = as.character(matrix(c(1:5)))
plan("multiprocess", workers = 32)
options(future.globals.maxSize=999999999999)
for(j in 1:length(types)){
for(i in 1:length(regions)){
meta=mixc@meta.data[,c('Region','CellType')]
if(sum(meta[,1]==regions[i] & meta[,2]==types[j]) < 20){
next
}
label=rep(0,dim(mixc)[2])
label[meta[,1]==regions[i] & meta[,2]==types[j]]=1
mixc <- AddMetaData(object = mixc,metadata = as.vector(label),col.name = "label")
markerstmp <- FindMarkers(mixc, ident.1 = 1,ident.2 =0,min.pct = 0.25,min.diff.pct=0.1,logfc.threshold=0.5,group.by = 'label')
markerstmp <- markerstmp[markerstmp[,'p_val_adj']<0.01,]
markerstmp <- markerstmp[markerstmp[,'avg_log2FC']>0,]
ds = 50
store[(1+(iiiii-1)*ds):(iiiii*ds),1] = row.names(markerstmp)[1:ds]
store[(1+(iiiii-1)*ds):(iiiii*ds),2] = rep(types[j],ds)
store[(1+(iiiii-1)*ds):(iiiii*ds),3] = rep(regions[i],ds)
store[(1+(iiiii-1)*ds):(iiiii*ds),4] = markerstmp[1:ds,'p_val_adj']
store[(1+(iiiii-1)*ds):(iiiii*ds),5] = markerstmp[1:ds,'avg_log2FC']
colnames(store)=c('gene','Cluster','Region','pvalue','log2foldchange')
iiiii = iiiii+1
print(paste(types[j],regions[i],sep=' '))
}}
write.csv(store,file = 'Signature_50.csv',row.names=F)

dise <- c('C','GP','P','SN')
Signature <- fread("Signature_50.csv",sep=",",stringsAsFactors=F,header=T)
Signature <- as.matrix(Signature)
Signature <- Signature[!is.na(Signature[,1]),]
y_lab=c()
iiiii=1
exp_ma=matrix(0,nr=length(Signature[,1]),nc=35)
for(j in 1:length(types)){
for(i in 1:length(regions)){
meta=mixc@meta.data[,c('Region','CellType')]
if(sum(meta[,1]==regions[i] & meta[,2]==types[j]) < 20){
next
}
cells=row.names(mixc@meta.data)[meta[,1]==regions[i] & meta[,2]==types[j]]
y_lab=c(y_lab,paste(types[j],dise[i],sep=' '))
exp_ma[,iiiii]=rowMeans(mixc[['RNA']]@data[Signature[,1],cells])
iiiii=iiiii+1
}}
row.names(exp_ma)=paste(Signature[,1],1:1750)
colnames(exp_ma)=y_lab

pdf('2c.pdf')
annotation_row <- data.frame(CellType=factor(Signature[,2],levels=c("D1-MSN","D2-MSN","D1/D2 Mixed","Interneuron","Da","Oligo","Astro","Microglia","Endothelial","OPC","Ependy")),Region=factor(Signature[,3],levels=c('caudate','putamen','globus pallidus','substantia nigra')))
row.names(annotation_row) <- paste(Signature[,1],1:1750)
pheatmap(exp_ma,scale="row",border_color=NA,cluster_rows=F,cluster_cols = F,show_rownames = F,show_colnames = T,annotation_row =annotation_row,annotation_names_row=T,annotation_names_col=F,color = c(colorRampPalette(colors = c("blue","white"))(50/2),colorRampPalette(colors = c("white","red"))(50/2)))
dev.off()



##########################################################################################
##########################################################################################
#############################################Number of DEGs detected by inter-regional differential expression analysis within the same cell type
iiiii=1
store = matrix(nc=10000,nr=30000)
colnames(store) = as.character(matrix(c(1:10000)))
plan("multiprocess", workers = 32)
options(future.globals.maxSize=999999999999)
reg <- c('caudate','globus pallidus','putamen','substantia nigra')
CellT <- c("D1-MSN","D2-MSN","D1/D2 Mixed","Interneuron","Oligo","Astro","Microglia","OPC","Endothelial")
for(i in 1:length(CellT)){
mix8 = subset(mixc, subset = CellType == CellT[i])
cott <- table(mix8@meta.data[,'Region'])
cott <- cott[cott > 20]
diseaseee2 <- names(cott)
all_dis <- reg[sort(match(diseaseee2,reg))]
bingn <- dim(matrix(all_dis))[1]
if(bingn < 2){
next
}
com <- combn(all_dis,2)
pairs_num = length(com)/2
for(j in 1:pairs_num){
mix8 = subset(mixc, subset = CellType == CellT[i])
mix8 = subset(mix8, subset = Region == com[1,j] | Region == com[2,j] )
markerstmp <- FindMarkers(mix8, ident.1 = com[1,j],ident.2 =com[2,j],min.pct = 0.25,min.diff.pct=0.1,logfc.threshold=0.5,group.by = 'Region')
markerstmp <- markerstmp[markerstmp[,'p_val_adj']<0.01,]
if(dim(markerstmp)[1] > 0){
ds = dim(markerstmp)[1]
store[1:ds,(1+(iiiii-1)*3)] = row.names(markerstmp)[1:ds]
store[1:ds,(2+(iiiii-1)*3)] = markerstmp[1:ds,'p_val_adj']
store[1:ds,(3+(iiiii-1)*3)] = markerstmp[1:ds,'avg_log2FC']
colnames(store)[((1+(iiiii-1)*3):(3+(iiiii-1)*3))] = c(paste(com[1,j],' vs ',com[2,j],' in ',CellT[i],sep=''),'pvalue','log2foldchange')
iiiii = iiiii+1
print(paste(com[1,j],' vs ',com[2,j],' in ',CellT[i],sep=''))
}
}
}
store[is.na(store)] = ''
store=store[,1:(3*iiiii+10)]
write.csv(store,file = 'DE_control_CellType_Region.csv',row.names=F)

DEseq2 <- fread("DE_control_CellType_Region.csv",sep=",",stringsAsFactors=F,header=T)
DEseq2 <- as.matrix(DEseq2)
x_label = c('C vs GP','C vs P','C vs SN','GP vs P','GP vs SN','P vs SN')
pdf('2e.pdf')
reg <- c('caudate','globus pallidus','putamen','substantia nigra')
for(z in 1:length(CellT)){
tmp_inx = grepl(CellT[z],colnames(DEseq2))
tmp_inx = which(tmp_inx == TRUE)
if(length(tmp_inx) == 0){
next}
tmp_inx=sort(c(tmp_inx,tmp_inx+1,tmp_inx+2))
allDE3 <- DEseq2[,tmp_inx]
len1 <- dim(allDE3)[2]/3
com <- combn(reg,2)
com <- paste(com[1,],'vs',com[2,],sep=' ')
countss <- rep(0,length(com))
countss2 <- rep(0,length(com))
for(i in 1:length(com)){
tmp_inx = grepl(com[i],colnames(allDE3))
tmp_inx = which(tmp_inx == TRUE)
if(length(tmp_inx) == 0){
next}
geneb1 = as.numeric(as.vector(as.matrix(allDE3[,tmp_inx+2])))
geneb1 = as.numeric(geneb1[!is.na(geneb1)])
numm = sum(geneb1>0)
numm3 = sum(geneb1<0)
countss[i] = numm
countss2[i] = numm3
}
p <- ggplot(data=NULL) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(x=x_label, y=countss)) +
theme(plot.title = element_text(hjust = 0.5,size=18),axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=15),axis.text.y = element_text(color = "black",size=15)) +
labs(x ='',y='',title =paste('DE Gene Counts of',CellT[z],sep=' ')) +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(limits=x_label) +
new_scale_color() +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(x=x_label, y=(-countss2))) +
geom_hline(aes(yintercept=0))#+coord_cartesian(ylim = c(-50, 50))
print(p)
}
dev.off()

##########################################################################################
##########################################################################################
#############################################Plotting UMAP of BD donor
#############################################"UMAP_BD.csv" is in folder namely "materials"
mixd=subset(Merge, subset = Disease == 'BD')

mixd <- NormalizeData(mixd, normalization.method = "LogNormalize", scale.factor = 10000)
mixd <- FindVariableFeatures(mixd, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mixd)
mixd <- ScaleData(mixd, features = all.genes)
mixd <- RunPCA(mixd, features = VariableFeatures(object = mixd))
mixd <- FindNeighbors(mixd, dims = 1:20)
mixd <- RunUMAP(mixd, dims = 1:20)
UMAP_BD=fread("UMAP_BD.csv",sep=",",stringsAsFactors=F,header=T)
UMAP_BD=UMAP_BD[,-1]
UMAP_BD=as.matrix(UMAP_BD)
row.names(UMAP_BD)=colnames(mixd)
mixd@reductions$umap@cell.embeddings=UMAP_BD

pdf('3a.pdf')
DimPlot(mixd, reduction = "umap",label = T,group.by = 'CellType')
dev.off()



##########################################################################################
##########################################################################################
#############################################Cell ratio of regional composition per cell type (for BD donor)
regions <- c('caudate','globus pallidus','putamen','substantia nigra')
types <- c("D1-MSN-C","D1-MSN-P","D2-MSN-C","D2-MSN-P","D1/D2 Mixed-1","D1/D2 Mixed-2","InterN-1","InterN-2","InterN-3","Oligo-1","Oligo-2","Astro-1","Astro-2","Microglia","Endothelial-1","Endothelial-2","Endothelial-3","OPC","Ependy")
cl_di  <- matrix(nr=length(types),nc=4)
row.names(cl_di) <- types
colnames(cl_di) <- regions
for(i in 1:length(types)){
for(j in 1:4){
tee <- mixd@meta.data[,c('Region','SubCellType')]
tee <- tee[which(tee[,1]==regions[j]),]
tee <- tee[which(tee[,2]==types[i]),]
cl_di[i,j] <- dim(tee)[1]
}}
cs <- colSums(cl_di)
tee2 <- cl_di
for(i in 1:4){
tee2[,i] <- tee2[,i]/cs[i]
}

tee2 = cl_di/rowSums(cl_di)
cll <- tee2[,1]
cll2 <- names(cll)
cll3 <- rep(colnames(tee2)[1],dim(tee2)[1])
for(i in 2:dim(tee2)[2]){
clll <- tee2[,i]
clll2 <- names(clll)
clll3 <- rep(colnames(tee2)[i],dim(tee2)[1])
cll <- c(cll,clll)
cll2 <- c(cll2,clll2)
cll3 <- c(cll3,clll3)
}

pdf('3b.pdf')
ggplot(data=NULL, aes(x=cll2, y=cll)) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(fill = cll3)) +
theme(axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=9)) +
labs(x ='',y='Cell ratio') +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(limits=row.names(cl_di)) +
labs(fill = "Region")
dev.off()



##########################################################################################
##########################################################################################
#############################################Number of DEGs detected by differential expression analysis of each cell type between BD and controls in the
#############################################caudate (C), putemen (P), globus pallidus (GP), substantia nigra (SN), respectively.
Control=mixc
Disease=mixd
Disease <- AddMetaData(object = Disease,metadata = as.vector(rep('BD',dim(Disease)[2])),col.name = "Disease")
Control <- AddMetaData(object = Control,metadata = as.vector(rep('Control',dim(Control)[2])),col.name = "Disease")

#for caudate
pairs=c(c('D1-MSN-C','D1-MSN'),c('D2-MSN-C','D2-MSN-1'),c('D1/D2 Mixed-1','D1/D2 Mixed-1'),c('D1/D2 Mixed-2','D1/D2 Mixed-2'),c('InterN-1','InterN-1'),c('InterN-2','InterN-2'),c('InterN-3','InterN-3'),c('Oligo-1','Oligo-1'),c('Astro-1','Astro-1'),c('Astro-2','Astro-2'),c('Microglia','Microglia'),c('Endothelial-1','Endothelial-1'),c('Endothelial-2','Endothelial-2'),c('Endothelial-3','Endothelial-3'),c('OPC','OPC'),c('Ependy','Ependy'))
cellt=pairs[seq(1,length(pairs),2)]
cellt[1]='D1-MSN'
cellt[2]='D2-MSN'
ana_region='caudate'
Disease = subset(Disease, subset = Region == ana_region)
Control = subset(Control, subset = Region == ana_region)

#for globus pallidus
pairs=c(c('InterN-1','InterN-1'),c('Oligo-1','Oligo-1'),c('Oligo-2','Oligo-2'),c('Astro-1','Astro-1'),c('Microglia','Microglia'),c('Endothelial-1','Endothelial-1'),c('Endothelial-2','Endothelial-2'),c('Endothelial-3','Endothelial-3'),c('OPC','OPC'))
cellt=pairs[seq(1,length(pairs),2)]
ana_region='globus pallidus'
Disease = subset(Disease, subset = Region == ana_region)
Control = subset(Control, subset = Region == ana_region)

#for putamen
pairs=c(c('D1-MSN-P','D1-MSN'),c('D2-MSN-P','D2-MSN-1'),c('D1/D2 Mixed-1','D1/D2 Mixed-1'),c('D1/D2 Mixed-2','D1/D2 Mixed-2'),c('InterN-1','InterN-1'),c('InterN-2','InterN-2'),c('InterN-3','InterN-3'),c('Oligo-1','Oligo-1'),c('Oligo-2','Oligo-2'),c('Astro-1','Astro-1'),c('Astro-2','Astro-2'),c('Microglia','Microglia'),c('Endothelial-1','Endothelial-1'),c('Endothelial-2','Endothelial-2'),c('Endothelial-3','Endothelial-3'),c('OPC','OPC'))
cellt=pairs[seq(1,length(pairs),2)]
cellt[1]='D1-MSN'
cellt[2]='D2-MSN'
ana_region='putamen'
Disease = subset(Disease, subset = Region == ana_region)
Control = subset(Control, subset = Region == ana_region)

#for substantia nigra
pairs=c(c('InterN-1','InterN-1'),c('Oligo-1','Oligo-1'),c('Oligo-2','Oligo-2'),c('Astro-1','Astro-1'),c('Microglia','Microglia'),c('Endothelial-1','Endothelial-1'),c('Endothelial-2','Endothelial-2'),c('Endothelial-3','Endothelial-3'),c('OPC','OPC'))
cellt=pairs[seq(1,length(pairs),2)]
ana_region='substantia nigra'
Disease = subset(Disease, subset = Region == ana_region)
Control = subset(Control, subset = Region == ana_region)


iiiii=1
store = matrix(nc=1000,nr=30000)
colnames(store) = as.character(matrix(c(1:1000)))
plan("multiprocess", workers = 32)
options(future.globals.maxSize=999999999999)
for(i in 1:length(cellt)){
if(table(Disease@meta.data[,'SubCellType'])[pairs[2*i-1]] < 20 | table(Control@meta.data[,'SubCellType'])[pairs[2*i]] < 20){
next
}
DiseaseSub = subset(Disease, subset = SubCellType == pairs[2*i-1])
ControlSub = subset(Control, subset = SubCellType == pairs[2*i])
mixcd = merge(DiseaseSub,ControlSub)
mixcd <- NormalizeData(mixcd, normalization.method = "LogNormalize", scale.factor = 10000)
markerstmp <- FindMarkers(mixcd, ident.1 = 'BD',ident.2 ='Control',min.pct = 0.25,min.diff.pct=0.1,logfc.threshold=0.5,group.by = 'Disease')
markerstmp <- markerstmp[markerstmp[,'p_val_adj']<0.01,]
if(dim(markerstmp)[1] > 0){
ds = dim(markerstmp)[1]
store[1:ds,(1+(iiiii-1)*3)] = row.names(markerstmp)[1:ds]
store[1:ds,(2+(iiiii-1)*3)] = markerstmp[1:ds,'p_val_adj']
store[1:ds,(3+(iiiii-1)*3)] = markerstmp[1:ds,'avg_log2FC']
colnames(store)[((1+(iiiii-1)*3):(3+(iiiii-1)*3))] = c(paste('BD',' vs ','control',' in ',ana_region,' in ',cellt[i],sep=''),'pvalue','log2foldchange')
iiiii = iiiii+1
print(paste('BD',' vs ','control',' in ',ana_region,' in ',cellt[i],sep=''))
}
}
store[is.na(store)] = ''
write.csv(store,file = paste('DE_BDVScontrol_CellType_in_',ana_region,'.csv',sep=''),row.names=F)

DEseq2 <- fread(paste('DE_BDVScontrol_CellType_in_',ana_region,'.csv',sep='') ,sep=",",stringsAsFactors=F,header=T)
DEseq2 <- as.matrix(DEseq2)
x_label =cellt
pdf(paste('3c-g_',ana_region,'.pdf',sep=''))
reg <- cellt
for(z in 1:1){
allDE3 <- DEseq2
len1 <- dim(allDE3)[2]/3
com <- reg
countss <- rep(0,length(com))
countss2 <- rep(0,length(com))
for(i in 1:length(com)){
tmp_inx = grepl(com[i],colnames(allDE3))
tmp_inx = which(tmp_inx == TRUE)
if(length(tmp_inx) == 0){
next}
geneb1 = as.numeric(as.vector(as.matrix(allDE3[,tmp_inx+2])))
geneb1 = as.numeric(geneb1[!is.na(geneb1)])
numm = sum(geneb1>0)
numm3 = sum(geneb1<0)
countss[i] = numm
countss2[i] = numm3
}
p <- ggplot(data=NULL) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(x=x_label, y=countss)) +
theme(plot.title = element_text(hjust = 0.5,size=18),axis.text.x = element_text(angle = 45, hjust =1, vjust = 1,color = "black",size=15),axis.text.y = element_text(color = "black",size=15)) +
labs(x ='',y='',title =paste('DE Gene Counts of BD VS Control in ',ana_region,sep='')) +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(limits=x_label) +
new_scale_color() +
geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.7),aes(x=x_label, y=(-countss2))) +
geom_hline(aes(yintercept=0)) 
#+
#coord_cartesian(ylim = c(-50, 50))
print(p)
}
dev.off()



##########################################################################################
##########################################################################################
#############################################WGCNA analysis
#############################################Gene modules identified in DEGs between human BD and control basal ganglia nuclei using WGCNA
#############################################"WGCNA_genes.csv" and "geneTree.rds" are in folder namely "materials"
DEGs <- fread("WGCNA_genes.csv",sep=",",stringsAsFactors=F,header=T)
DEGs <- names(table(as.vector(as.matrix(DEGs))))
DEGs=DEGs[-1]
DEGs=DEGs[-1]

enableWGCNAThreads(32)
allowWGCNAThreads(32)
options(stringsAsFactors = FALSE)
counts <- Merge[['RNA']]@data[DEGs,]
cells <- colnames(counts)

nSets = 1
multiExpr = vector(mode = 'list',length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(as.matrix(counts))))
exprSize = checkSets(multiExpr)
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
multiExpr[[1]] = list(data = as.data.frame(t(as.matrix(counts[unlist(gsg$goodGenes),unlist(gsg$goodSamples)]))))
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
exprSize = checkSets(multiExpr)

cor = WGCNA::cor
type = "signed"
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(multiExpr[[1]]$data, powerVector = powers, verbose = 5,networkType=type)

softPower <- sft$powerEstimate
softPower <- 3
cor <- WGCNA::cor
adjacency = adjacency(multiExpr[[1]]$data, power = softPower,type="signed")#corFnc = "bicor"
TOM = TOMsimilarity(adjacency,TOMType = type)
dissTOM = 1-TOM

#geneTree = hclust(as.dist(dissTOM))
geneTree = readRDS(file='geneTree.rds')
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,minClusterSize = minModuleSize, cutHeight = NULL,deepSplit = 3)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(multiExpr[[1]]$data, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25
merge = mergeCloseModules(multiExpr[[1]]$data, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors
mergedMEs = merge$newMEs
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
pdf('5a.pdf')
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()



##########################################################################################
##########################################################################################
#############################################WGCNA analysis
#############################################Heatmap indicating the correlation between module eigengenes and cell type identities
CellT=c('D1-MSN','D2-MSN','D1/D2 Mixed','Interneuron','Da','Oligo','Astro','Microglia','Endothelial','OPC','Ependy')
j=1
datTraits <- Merge@meta.data[cells,'CellType']
datTraits <- match(datTraits,CellT[j])
datTraits[is.na(datTraits)] <- '0'
datTraits <- as.numeric(datTraits)
datTraits <- as.matrix(as.vector(datTraits))
row.names(datTraits) <- cells
colnames(datTraits) <- CellT[j]
for(i in 2:length(CellT)){
j=i
datTraits2 <- Merge@meta.data[cells,'CellType']
datTraits2 <- match(datTraits2,CellT[j])
datTraits2[is.na(datTraits2)] <- '0'
datTraits2 <- as.numeric(datTraits2)
datTraits2 <- as.matrix(as.vector(datTraits2))
row.names(datTraits2) <- cells
colnames(datTraits2) <- CellT[j]
datTraits <- cbind(datTraits,datTraits2)
}

oMEs = orderMEs(mergedMEs)
nGenes = ncol(multiExpr[[1]]$data)
nSamples = nrow(multiExpr[[1]]$data)
moduleTraitCor = cor(oMEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf('5b.pdf')
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = as.matrix(moduleTraitCor),
xLabels = colnames(datTraits),
yLabels = names(oMEs),
ySymbols = names(oMEs),
colorLabels = FALSE,
colors = blueWhiteRed(80),
#textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.35,
zlim = c(-1,1),
main = paste("Module-trait relationships"),cex.lab.x=1,cex.lab.y=1)
dev.off()



##########################################################################################
##########################################################################################
#############################################WGCNA analysis
#############################################Heatmap indicating the correlation between module eigengenes and other traits
datTraits <- Merge@meta.data[cells,2]
datTraits <- as.numeric(datTraits)
datTraits <- as.matrix(as.vector(datTraits))
row.names(datTraits) <- cells
colnames(datTraits) <- 'nCount_RNA'
datTraits2 <- Merge@meta.data[cells,3]
datTraits2 <- as.numeric(datTraits2)
datTraits2 <- as.matrix(as.vector(datTraits2))
colnames(datTraits2) <- 'nFeature_RNA'
datTraits <- cbind(datTraits,datTraits2)

disease <- Merge@meta.data[cells,'Disease']
disease[which(disease == 'Control')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'Control'
datTraits <- cbind(datTraits,disease)

disease <- Merge@meta.data[cells,'Disease']
disease[which(disease == 'BD')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'BD'
datTraits <- cbind(datTraits,disease)

disease <- Merge@meta.data[cells,'Region']
disease[which(disease == 'caudate')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'caudate'
datTraits <- cbind(datTraits,disease)
disease <- Merge@meta.data[cells,'Region']
disease[which(disease == 'globus pallidus')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'globus pallidus'
datTraits <- cbind(datTraits,disease)
disease <- Merge@meta.data[cells,'Region']
disease[which(disease == 'putamen')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'putamen'
datTraits <- cbind(datTraits,disease)
disease <- Merge@meta.data[cells,'Region']
disease[which(disease == 'substantia nigra')] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
colnames(disease) <- 'substantia nigra'
datTraits <- cbind(datTraits,disease)

oMEs = orderMEs(mergedMEs)
nGenes = ncol(multiExpr[[1]]$data)
nSamples = nrow(multiExpr[[1]]$data)
moduleTraitCor = cor(oMEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 1), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf('5c.pdf')
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = as.matrix(moduleTraitCor),
xLabels = colnames(datTraits),
yLabels = names(oMEs),
ySymbols = names(oMEs),
colorLabels = FALSE,
colors = blueWhiteRed(80),
#textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.35,
zlim = c(-1,1),
main = paste("Module-trait relationships"),cex.lab.x=1,cex.lab.y=1)
dev.off()




##########################################################################################
##########################################################################################
#############################################WGCNA analysis
#############################################Bar plot showing the eigengene expression of the gene modules among cell types and regions
pdf('5d.pdf')
for(z in 1:dim(mergedMEs)[2]){
print(z)
CellT=c('D1-MSN','D2-MSN','D1/D2 Mixed','Interneuron','Da','Oligo','Astro','Microglia','Endothelial','OPC','Ependy')
diseaseee <- c('caudate','globus pallidus','putamen','substantia nigra')
i=1
j=1
infor = Merge@meta.data[,c("CellType","Region")]
infor = infor[infor[,1]==CellT[i],]
infor = infor[infor[,2]==diseaseee[j],]
cellss <- row.names(infor)
m1 <- mergedMEs[cellss,z]
m2 <- substring(colnames(mergedMEs)[z],3)
m3 <- rep(CellT[i],dim(infor)[1])
m4 <- rep(diseaseee[j],dim(infor)[1])
for(i in 1:length(CellT)){
for(j in 1:4){
if(i != 1 | j != 1){
infor = Merge@meta.data[,c("CellType","Region")]
infor = infor[infor[,1]==CellT[i],]
infor = infor[infor[,2]==diseaseee[j],]
cellss <- row.names(infor)
m1 <- c(m1,mergedMEs[cellss,z])
m3 <- c(m3,rep(CellT[i],dim(infor)[1]))
m4 <- c(m4,rep(diseaseee[j],dim(infor)[1]))
}}}
m4 <- factor(m4, levels= c('caudate','globus pallidus','putamen','substantia nigra'))
p <- ggplot(data=NULL, aes(x=m3, y=m1)) +
theme_bw() +
theme(panel.grid=element_line(colour=NA)) +
geom_boxplot(width=0.5,position=position_dodge(0.7),aes(fill = m4),outlier.size = 0.001,outlier.shape=15,outlier.colour = 'black') +
theme(axis.text.x = element_text(angle = 70, hjust =1, vjust = 1,color = "black",size=9)) +
labs(x ='',y=paste('Eigengene','(',m2,') ','Expression',sep='')) +
scale_x_discrete(limits=CellT) +
labs(fill = "Region")
print(p)
}
dev.off()



##########################################################################################
##########################################################################################
#############################################WGCNA analysis
#############################################Network of hub genes identified in module green
diseaseee <- c('BD','Control')
cor <- WGCNA::cor
datKME = signedKME(multiExpr[[1]]$data,mergedMEs,outputColumnName="")
j=1
disease <- Merge@meta.data[row.names(multiExpr[[1]]$data),'Disease']
disease[which(disease == diseaseee[j])] <- 1
disease[which(disease != 1)] <- 0
disease <- as.matrix(as.vector(as.numeric(disease)))
GS1= as.numeric(cor(disease,multiExpr[[1]]$data, use="p"))
wggenes = colnames(multiExpr[[1]]$data)
names(GS1) <- wggenes
names(moduleColors) <- wggenes
MM=c()
md='green'
mogene = wggenes[which(moduleColors == md)]
write.csv(mogene,paste('Module_',md,'_Genes.csv',sep=''),row.names=F)
for(i in 1:length(mogene)){
MM <- c(MM,datKME[mogene[i],md])
}
names(MM) <- mogene
topMM <- MM[order(MM,decreasing=TRUE)[1:60]]
topMMgene <- names(topMM)
topMM_GS <- GS1[topMMgene]
topGS <- topMM_GS[order(abs(topMM_GS),decreasing=TRUE)[1:30]]
topGS_MM <- topMM[names(topGS)]
topGS_gene <- names(topGS)
topGS_color <- moduleColors[topGS_gene]
tmplen <- length(topGS_MM)
outputcsv = matrix(nr=tmplen,nc=7)
colnames(outputcsv) = c('gene','module','group','GS','pvalue of GS','MM','pvalue of MM')
outputcsv[1:tmplen,1] = topGS_gene
outputcsv[1:tmplen,2] = topGS_color
outputcsv[1:tmplen,3] = rep(diseaseee[j],tmplen)
outputcsv[1:tmplen,4] = topGS
outputcsv[1:tmplen,5] = corPvalueStudent(topGS, nSamples)
outputcsv[1:tmplen,6] = topGS_MM
outputcsv[1:tmplen,7] = corPvalueStudent(topGS_MM, nSamples)
write.csv(outputcsv,file='Hubgene.csv',row.names=F)

hub <- topGS_gene
colnames(TOM)<-colnames(adjacency)
row.names(TOM) <- row.names(adjacency)
exportNetworkToCytoscape(adjacency[hub,hub],
             edgeFile = 'edges.txt',
             nodeFile = 'nodes.txt',
             nodeAttr = topGS_MM,
             weighted = TRUE, threshold = quantile(adjacency[hub,hub],0.9),
             nodeNames = row.names(adjacency[hub,hub]))




##########################################################################################
##########################################################################################
#############################################SCENIC analysis
library(SCENIC)
library(SCopeLoomR)

dbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc",dbDir='/home/user/SCENIC',dbs = dbs,nCores=32) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#For example, caudate
regions <- c('caudate','globus pallidus','putamen','substantia nigra')
mix8 <- subset(Merge, subset = Region==regions[1])

exprMat=mix8[['RNA']]@counts
exprMat=as.matrix(exprMat)

cellInfo=mix8@meta.data[,c('CellType','SubCellType','nCount_RNA','nFeature_RNA','Disease')]
colnames(cellInfo)[3:4]=c('nUMI','nGene')
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

write.csv(t(exprMat_filtered), file = "subMatrix.csv")

#run following code in Python
import os, sys 
import loompy as lp
import numpy as np
import scanpy as sc

x=sc.read_csv("subMatrix.csv")
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("subMatrix.loom",x.X.transpose(),row_attrs,col_attrs)
#run Python end

#run following code in bash
dir=/home/user/SCENIC
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

input_loom=subMatrix.loom
pyscenic grn --num_workers 32 --output subMatrix.tsv --method grnboost2 $input_loom $tfs
pyscenic ctx subMatrix.tsv $feather --annotations_fname $tbl --expression_mtx_fname $input_loom --mode "dask_multiprocessing" --output reg.csv --num_workers 32 --mask_dropouts
pyscenic aucell $input_loom reg.csv --output pySCENIC_results.loom --num_workers 32
#run bash end

loom <- open_loom('pySCENIC_results.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
close_loom(loom)

sub_regulonAUC <- regulonAUC[,match(colnames(exprMat_filtered),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(exprMat_filtered))

TF=c()
Target=c()
for(i in 1:length(regulons)){
TF_name=names(regulons)[i]
TF_target=regulons[[TF_name]]
TF=c(TF,rep(TF_name,length(TF_target)))
Target=c(Target,TF_target)
}
Out=data.frame(TF=TF,Target=Target)
#corresponding to Figure 6c
write.csv(Out,file='TF_Target.csv',row.names=F

Keep_TF=row.names(table(TF))[table(TF) >=20]
Keep_sub_regulonAUC=sub_regulonAUC[Keep_TF,]

Keep_TF_newName=gsub('[(+)]','',row.names(table(TF))[table(TF) >=20])
Keep_TF_newName=paste(Keep_TF_newName,'(',table(TF)[table(TF) >=20],'g)',sep='')
row.names(Keep_sub_regulonAUC)=Keep_TF_newName

Disease_name=c('BD','Control')
Disease_bar=c()
out=list()
for(i in 1:length(Disease_name)){
cellInfo_sub=cellInfo[which(cellInfo[,'Disease']==Disease_name[i]),]
dim1=dim(Keep_sub_regulonAUC)[1]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo_sub), cellInfo_sub$CellType),function(cells) rowMeans(matrix(Keep_sub_regulonAUC@assays@data@listData$AUC[,cells],nr=dim1)))
row.names(regulonActivity_byCellType)=row.names(Keep_sub_regulonAUC)
keep_SubCellType=names(table(cellInfo_sub[,'CellType']))[table(cellInfo_sub[,'CellType']) > 20]
regulonActivity_byCellType=regulonActivity_byCellType[,keep_SubCellType]
if(Disease_name[i]=='BD'){
relative=c("D1-MSN","D2-MSN","D1/D2 Mixed","Interneuron","Oligo","Astro","Microglia","Endothelial","OPC","Ependy")
keep_SubCellType=relative[sort(match(keep_SubCellType,relative))]
regulonActivity_byCellType=regulonActivity_byCellType[,keep_SubCellType]
Disease_bar=c(Disease_bar,rep(Disease_name[i],length(keep_SubCellType)))
}
if(Disease_name[i]=='Control'){
relative=c("D1-MSN","D2-MSN","D1/D2 Mixed","Interneuron","Da","Oligo","Astro","Microglia","Endothelial","OPC","Ependy")
keep_SubCellType=relative[sort(match(keep_SubCellType,relative))]
regulonActivity_byCellType=regulonActivity_byCellType[,keep_SubCellType]
Disease_bar=c(Disease_bar,rep(Disease_name[i],length(keep_SubCellType)))
}
out[[Disease_name[i]]]=regulonActivity_byCellType
}
inter_celltype=intersect(colnames(out[[Disease_name[1]]]),colnames(out[[Disease_name[2]]]))
name1_out=out[[Disease_name[1]]][,inter_celltype]
name2_out=out[[Disease_name[2]]][,inter_celltype]
out=name1_out-name2_out
out_row_max=apply(out,1,function(x){max(abs(x))})
out=out[out_row_max>=0.01,]

pdf('6a.pdf')
bk1 <- c(seq(range(out)[1],0,by=0.0001))
bk2 <- c(seq(0,range(out)[2],by=0.0001))
pheatmap(out,main=paste(Disease_name[1],'-',Disease_name[2]),scale='none',border_color=NA,cluster_rows=T,treeheight_row=25,fontsize_row=5,cluster_cols = F,show_rownames = T,show_colnames = T,annotation_names_row=T,annotation_names_col=T,color = c(colorRampPalette(colors = c("blue","white"))(length(bk1)),colorRampPalette(colors = c("white","red"))(length(bk2))))
dev.off()
write.csv(out,file='TF_CellType_minus_heatmap.csv')



##########################################################################################
##########################################################################################
#############################################CellChat analysis
Sys.setenv(RETICULATE_PYTHON="Python_Work_Dictionary")
library(CellChat)

Disease=mixd
Control=mixc

#For example, caudate
diseaseee <- c('caudate','globus pallidus','putamen','substantia nigra')
Controlsub <- subset(Control, subset = Region==diseaseee[1])
Diseasesub <- subset(Disease, subset = Region==diseaseee[1])

keep_SubCellType1=names(table(Controlsub@meta.data[,'CellType']))[table(Controlsub@meta.data[,'CellType']) > 20]
keep_SubCellType2=names(table(Diseasesub@meta.data[,'CellType']))[table(Diseasesub@meta.data[,'CellType']) > 20]
inter_celltype=intersect(keep_SubCellType1,keep_SubCellType2)
inter_gene=intersect(row.names(Controlsub),row.names(Diseasesub))
Controlcell=!is.na(match(Controlsub@meta.data[,'CellType'],inter_celltype))
Diseasecell=!is.na(match(Diseasesub@meta.data[,'CellType'],inter_celltype))

cellchat_control <- createCellChat(as.matrix(Controlsub@assays$RNA@data[inter_gene,Controlcell]), meta = Controlsub@meta.data[Controlcell,], group.by = "CellType")
cellchat_control <- setIdent(cellchat_control, ident.use = "CellType", levels = inter_celltype)
table(cellchat_control@idents)

cellchat_disease <- createCellChat(as.matrix(Diseasesub@assays$RNA@data[inter_gene,Diseasecell]), meta = Diseasesub@meta.data[Diseasecell,], group.by = "CellType")
cellchat_disease <- setIdent(cellchat_disease, ident.use = "CellType", levels = inter_celltype)
table(cellchat_disease@idents)

cellchat_control@DB <- CellChatDB.human
cellchat_control <- subsetData(cellchat_control)
future::plan("multiprocess", workers = 16)
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_control <- projectData(cellchat_control, PPI.human)
cellchat_control <- computeCommunProb(cellchat_control, raw.use = FALSE, population.size = TRUE)
cellchat_control <- filterCommunication(cellchat_control)
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)

cellchat_disease@DB <- CellChatDB.human
cellchat_disease <- subsetData(cellchat_disease)
future::plan("multiprocess", workers = 16)
cellchat_disease <- identifyOverExpressedGenes(cellchat_disease)
cellchat_disease <- identifyOverExpressedInteractions(cellchat_disease)
cellchat_disease <- projectData(cellchat_disease, PPI.human)
cellchat_disease <- computeCommunProb(cellchat_disease, raw.use = FALSE, population.size = TRUE)
cellchat_disease <- filterCommunication(cellchat_disease)
cellchat_disease <- computeCommunProbPathway(cellchat_disease)
cellchat_disease <- aggregateNet(cellchat_disease)

cc.list <- list(BD=cellchat_disease,Control=cellchat_control)
cellchat <- mergeCellChat(cc.list, add.names = names(cc.list), cell.prefix = TRUE)

pdf('7a.pdf')
gg1 <- compareInteractions(cellchat, show.legend = F,group = c('1','2'),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F,group = c('1','2'),measure = "weight")
p <- gg1 + gg2
print(p)
dev.off()

pdf('7b-c.pdf')
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count",color.edge=c("#2166ac","#b2182b"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.edge=c("#2166ac","#b2182b"))
dev.off()

pdf('7d.pdf')
rankNet(cellchat, mode = "comparison",measure = c("weight"),stacked = T,do.stat =T,font.size = 5)
dev.off()
