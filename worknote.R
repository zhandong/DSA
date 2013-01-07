
### load the benchmark data and combine the replicates by taking their mean value

lbl.data = read.table('GSE19830_series_matrix.txt',row.names=1,header=T,sep='\t')
lbl.data = as.matrix(lbl.data)

lbl.sdata=apply(2^lbl.data,1,function(r){
	tapply(r,rep(1:14,each=3),mean)
})
lbl.sdata= t(lbl.sdata)


### load the true weight
weight= read.table("liver_brain_lung_weight.txt",header=T)/100

### load the tissue specific markers and conver to list data type

tissue_markers = read.table('tissue specific genes.txt',header=F,sep='\t')

tm = list()
tm[[1]]=rownames(lbl.data) %in% tissue_markers[tissue_markers[,2]== "liver",1]
tm[[2]]=rownames(lbl.data) %in% tissue_markers[tissue_markers[,2]== "brain",1]
tm[[3]]=rownames(lbl.data) %in% tissue_markers[tissue_markers[,2]== "lung",1]



