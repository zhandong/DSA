### DSA with g_n missing

### load the benchmark data and combine the replicates by taking their mean value

lbl.data = read.table('GSE19830_series_matrix.txt',row.names=1,header=T,sep='\t')
lbl.data = as.matrix(lbl.data)

lbl.sdata=apply(2^lbl.data,1,function(r){
	tapply(r,rep(1:14,each=3),mean)
})
lbl.sdata= t(lbl.sdata)

### load the tissue specific genes in ts
tsgenes=read.table('../Fig 2/tissue specific genes.txt')


splist = list()
splist[[1]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='liver',1]
splist[[2]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='brain',1]
splist[[3]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='lung',1]


#### estimate weight


truew= read.table("../../../liver_brain_lung_weight.txt",header=T)
truew = truew[-(1:3),]
truew = as.matrix(truew)/100

w=t(weight_estimate_1(lbl.sdata[,4:14],splist))

#### if we don't know which one is good. 
splist = list()
splist[[1]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='liver',1]
splist[[2]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='brain',1]
#splist[[3]]=rownames(lbl.sdata) %in% tsgenes[tsgenes[,2]=='lung',1]

tmp=weight_estimate_missing(lbl.sdata[,4:14],splist,"lm")

tmp2=matrix(unlist(tmp),ncol=3,byrow=TRUE)

a=seq(12,nrow(tmp2),by=12)

mse=(tmp2[a,1])
fc=apply(log2(lbl.sdata),1,function(r){r[3]-mean(r[1],r[2])})
plot(fc,mse)
for(i in seq(12,nrow(tmp2),by=12)){
	 tmpw= tmp2[(i-11):(i-1),]
	 err=c(err,sum((tmpw-truew)^2)/(3*11))
	}
plot(fc,log(err))
plot(fc,log(mse))
plot(log(mse),log(err))

mseSorted=sort(mse,decreasing=T,index.return=TRUE)
aw=matrix(0,nrow=11,ncol=3)

for( i in mseSorted$ix[1:2])
{
		 tmpw= tmp2[((i-1)*12+1):(i*12-1),]
		 
		 aw= aw+tmpw
		 plot(tmpw,truew)
	
}
plot(aw/100,truew)








