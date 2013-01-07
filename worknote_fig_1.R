source('utils.R')

### load the benchmark data and combine the replicates by taking their mean value

lbl.data = read.table('GSE19830_series_matrix.txt',row.names=1,header=T,sep='\t')
lbl.data = as.matrix(lbl.data)

lbl.sdata=apply(2^lbl.data,1,function(r){
	tapply(r,rep(1:14,each=3),mean)
})
lbl.sdata= t(lbl.sdata)

### load the true weight
weight= read.table("liver_brain_lung_weight.txt",header=T)/100


### tissue reconstruction from liver brain and lung


lbl.hat.linear =lbl.sdata[,1:3]%*%t(weight) 
lbl.hat.log = 2^(log2(lbl.sdata[,1:3])%*%t(weight))

pdf("supplementary reconstruction 1 .pdf", height=8.5,width=11)
par(mfrow=c(3,4))
for(i in 4:14){
	plot.dens(log2(lbl.sdata[,i]),log2(lbl.hat.linear[,i]),cex=0.5,cex.axis=1)
	}
dev.off()	

pdf("supplementary reconstruction 2 .pdf", height=8.5,width=11)
par(mfrow=c(3,4))
for(i in 4:14){
plot.dens(log2(lbl.sdata[,i]),log2(lbl.hat.log[,i]),cex=0.5,cex.axis=1)
	}
dev.off()	

### deconvolution 
ob = lbl.sdata[,4:14]
ob= as.matrix(ob)
w = weight[-(1:3),]
w = as.matrix(w)

w2=w[-1,]
ob2 = ob[,-1]


pure.original = t(2^coef(lm.fit(w2,log2(t(ob2)))))
pure.log.hat = t(2^coef(lm.fit(w2,log2(t(ob2)))))
pure.linear.hat =QP(ob,w,l=4.3,u=22029)

par(mfrow=c(3,2))
for(i in 1:3)
{
	plot.dens((lbl.sdata[,i]),(pure.original[,i]),cex=0.5,cex.axis=1.2)
	print(cor((lbl.sdata[,i]),(pure.original[,i])))
	plot.dens((lbl.sdata[,i]),(pure.linear.hat[,i]),cex=0.5,cex.axis=1.2)
	print(cor((lbl.sdata[,i]),(pure.linear.hat[,i])))
}


#######################


out=c()
for (i in 1:300){
index=sample(1:11,6)
ob_p = ob[,index]
w_p = w[index,]
#tmp = t(coef(lm.fit(w_p,(t(ob_p)))))
tmp = QP(ob_p,w_p,l=4.3,u =22029)
print(i)
out=cbind(out,tmp)

}

pure.sample = pure.linear.hat

out[out<4.3]=NA
for (i in 1:3){
	tmp=(rowMeans(out[,seq(i,90,3)],na.rm=T))
	tmp[is.na(tmp)]=4.3
	pure.sample[,i]=tmp
}


par(mfrow=c(3,2))
for(i in 1:3)
{
	plot.dens(log2(lbl.sdata[,i]),log2(pure.original[,i]),cex=0.7,cex.axis=1.2)
	plot.dens(log2(lbl.sdata[,i]),log2(pure.linear.hat[,i]),cex=0.7,cex.axis=1.2)
#	plot.dens(log2(lbl.sdata[,i]),log2(pure.sample[,i]),cex=0.5,cex.axis=1.2)

}

par(mfrow=c(1,3))
for(i in 1:3){
plot.diff(log2(lbl.sdata[,i]),log2(pure.original[,i]), log2(pure.linear.hat[,i]), F)
}



source('mypca.R')

pdf("PCA.pdf",width=4.5, height =4.5)

src=lbl.sdata[,1:3]
sol=pure.linear.hat


pcadata =cbind(src,sol,ob[,c(2,3,1)],ob[,c(11,5,6)])
myPCA((pcadata[,1:3]),pcadata,label=rep(c("black","red","royalblue","seagreen"),c(3,3,3,3)), shape = rep(c(15:17),4))

pcadata = cbind(lbl.sdata[,1:3],pure.linear.hat,pure.original,ob[,c(2,3,1)],ob[,c(11,5,6)])
#pcadata = log2(pcadata)
myPCA((pcadata[,1:3]),pcadata,label=rep(c("black","red","purple","royalblue","seagreen"),c(3,3,3,3,3)), shape = rep(c(15:18),4))

dev.off()



tmp = cbind(lbl.sdata[,1:3],pure.original,pure.linear.hat)
myPCA(tmp[,1:3],tmp,rep(1:3,3))


##### supplementary 

TIGER = read.table('tissue specific genes.txt',header=F,sep='\t')

tissue_specific_list = list()
tissue_specific_list[[1]] = rownames(lbl.sdata)%in% as.character(TIGER[TIGER[,2]=="liver",1])
tissue_specific_list[[2]] =rownames(lbl.sdata)%in% as.character(TIGER[TIGER[,2]=="brain",1])
tissue_specific_list[[3]] = rownames(lbl.sdata)%in%as.character(TIGER[TIGER[,2]=="lung",1])

w=t(weight_estimate_1(lbl.sdata[,4:14], tissue_specific_list))

pure.linear.hat =QP(lbl.sdata[,4:14],w,l=4.3,u=22029)


par(mfrow=c(1,3))
for(i in 1:3)
{
	#plot.dens(log2(lbl.sdata[,i]),log2(pure.original[,i]),cex=0.7,cex.axis=1.2)
	plot.dens((lbl.sdata[,i]),(pure.linear.hat[,i]),cex=0.7,cex.axis=1.2)
	print(cor((lbl.sdata[,i]),(pure.linear.hat[,i])))
#	plot.dens(log2(lbl.sdata[,i]),log2(pure.sample[,i]),cex=0.5,cex.axis=1.2)

}
