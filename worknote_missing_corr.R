rm(list=ls())


### load the benchmark data and combine the replicates by taking their mean value


lbl.data = read.table('GSE19830_series_matrix.txt',row.names=1,header=T,sep='\t')
lbl.data = as.matrix(lbl.data)

lbl.sdata=apply(2^lbl.data,1,function(r){
	tapply(r,rep(1:14,each=3),mean)
})
lbl.sdata= t(lbl.sdata)


### load the tissue specific genes in ts
tsgenes=read.table('../Fig 2/tissue specific genes.txt')

truew= read.table("../../../liver_brain_lung_weight.txt",header=T)
truew = truew[-(1:3),]
truew = as.matrix(truew)/100



splist = list()
splist[[1]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='liver',1]
splist[[2]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='brain',1]
splist[[3]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='lung',1]


g12=rbind(colMeans(lbl.sdata[splist[[2]],4:14]),colMeans(lbl.sdata[splist[[3]],4:14]))

g12=(g12-rowMeans(g12))/sd(t(g12))

g12=colMeans(g12)

g3=lbl.sdata[,4:14]
g3=(g3-rowMeans(g3))/sd(t(g3))
tt=c()
for(i in 1:nrow(g3)){
	tt=c(tt,cor(g12,g3[i,]))
}

tti=sort(tt,index.return=TRUE)


index3=tti$ix[tti$x< -0.99]
index3= apply(lbl.sdata[index3,4:14],1,max) >1000
index3 = names(which(index3==TRUE))

splist[[1]]=rownames(lbl.sdata) %in% index3

splist[[1]]=rownames(lbl.sdata) %in% rownames(lbl.sdata[tti$ix[1:10],])

splist[[1]]=sample(rownames(lbl.sdata),10)


w=t(weight_estimate_1(lbl.sdata[,4:14],splist))

plot(truew,w); abline(0,1)

points(truew,w,col="red")
points(truew,w,col="green")

####### 




