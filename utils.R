plot.dens<-function(x,y,xlab="",ylab=".Rapp.history",xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),transparency=70,fc=log2(2),cex=1,cex.axis=1){

require(geneplotter)

mdata=cbind(x,y)


colMap=heat.colors(100)
colMap=densCols(mdata, colramp=colorRampPalette(colMap))
colMap=apply(col2rgb(colMap),2,function(r){rgb(r[1],r[2],r[3], transparency,maxColorValue=255)})

plot(mdata,col=colMap,pch=16,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,cex=cex,cex.axis=cex.axis)

abline(a=fc,b=1)
abline(a=-fc,b=1)

    
} 


plot.diff<-function(ref,est1,est2,log=T){
	if(log){
	x1=sort(abs(ref-est1))
	x2=sort(abs(ref-est2))
	y=rep(1,length(x1))
	y=1-cumsum(y)/length(y)
	plot(2^(x1),y,ylim=c(0,0.5))
	points(2^(x2),y,col="red")
	}
	else{
	
	x1=sort(abs(ref-est1))
	x2=sort(abs(ref-est2))
	y=rep(1,length(x1))
	y=1-(cumsum(y)/length(y))
	plot(log2(x1),y,xlim=c(0,10),ylim=c(0,1))
	points(log2(x2),y,col="red")
	}
	
}





plot.roc <- function(ref, diff1, diff2="null",add=F){
	require(ROCR)
	
	pred1 = prediction(diff1,ref)
	perf1 = performance(pred1,"tpr", "fpr")
	perf1.auc = performance(pred1,"auc")
	plot(perf1,main=perf1.auc@y.values[[1]][1], add=add, colorize=F)
	print(perf1.auc@y.values[[1]][1])
	if(diff2!="null"){
	pred2 = prediction(diff2,ref)
	perf2 = performance(pred2,"tpr", "fpr")
	plot(perf2,add=T,col='red')
	}
}
	
QP <- function(ob, weight, l, u, meq =0){
require(quadprog)	
	
sol= c()
for (id in 1:nrow(ob)){
A= weight
b= ob[id,]

Dmat = t(A)%*%A
dvec = b%*%(A)

#Amat = diag(rep(1,3))
#bvec=c(rep(3.15,3))

numC = ncol(weight)
Amat = diag(rep(1, numC))
Amat = rbind(Amat, diag(rep(-1, numC)))
Amat = t(Amat)
bvec=c(rep(l, numC),rep(-u, numC))

if(meq>0){
	Amat=cbind(rep(1,numC),Amat)
	bvec=c(1,bvec)
	}

#
#print("start")
#
#print(dim(Amat))
#print((bvec))
#print(dim(Dmat))
#print(dim(dvec))



sol=rbind(sol,solve.QP(Dmat,dvec,Amat,bvec=bvec,meq =meq)$solution)

}

	
return(sol)	
	}
	
	
blm <- function(x, H, theta0, sigma2){
	
	
	thetaU= solve((t(H)%*%H))%*%t(H)%*%x
	
	C = sigma2*solve(t(H)%*%H)
	
	m = theta0/(theta0+sum(diag(solve(t(H)%*%(H)))))
	
	return(m*thetaU)

	}	
	
	
myPCA = function(train,input=train,label,shape=10)
{
    PCAs=prcomp(t(train));
   mapped=t(input)%*%PCAs$rotation;
    plot(mapped,col=label,pch=shape,cex=2)
   #scatterplot3d(mapped[,1:3],color=label,pch=shape)
  # text(mapped[,1],mapped[,2],labels=colnames(input),cex=1,offset=0,pos=1)
    print(summary(PCAs))
    return(mapped)
}
	
	## this function I will select the number of specific gene for each cell type, by the proportion
## parameter Qp_linear is a logical variable with false linear model, true as QP model
## this estimate is better
weight_estimate_1 <- function(mix_ob, gene_list, method = "Lm")
{
   select_mix_ob <- matrix()
   for ( i in 1 : length(gene_list))
   {
       print("the sum of gene_list is ")
       print(sum(gene_list[[i]]))
       print("the length of third one")
       print(sum(gene_list[[3]]))
       #print(colMeans( mix_ob[gene_list[[i]], ]))
       if(i == 1)
       {
          if(sum(gene_list[[i]]) == 1)
            select_mix_ob <- as.matrix(mix_ob[gene_list[[i]], ])
          else
            select_mix_ob <-  as.matrix(colMeans( mix_ob[gene_list[[i]], ]))
       }
       else
       {
          if(sum(gene_list[[i]]) == 1)
            select_mix_ob <- cbind(select_mix_ob, as.matrix( mix_ob[gene_list[[i]], ]))
          else
            select_mix_ob <- cbind(select_mix_ob, as.matrix(colMeans( mix_ob[gene_list[[i]], ])))
       }
   }
   #print("dim of select_mix_ob is ")
   #print(dim(select_mix_ob))
   #print(head(select_mix_ob))
   y <- rep(1, times = nrow(select_mix_ob))
   #print("the length of y is ")
   #print(length(y))
   #stop()
   b_par <- numeric()
   ## we will use linear fit or QP method to get the estimator of paramter
   if(method != "Lm" && method != "Qp")
   {
      print("the 'method' paramter is wrong, should be 'Lm' or 'Qp'")
      stop()
   }
   if(method == "Lm")
   {
      lmob <- lm( y ~  -1 + (select_mix_ob))
      
      b_par <- coef(lmob)
   }
   else
   {
      print("before call QP function ")
     
      Qp <- QP(y/10000, (select_mix_ob)/10000,3,2^14)
      b_par <- Qp * 10000
      print(b_par)
      stop()
   }
   len <- as.integer(length(gene_list))
   par_matrix <- diag(b_par, len, len) 
   print("par_matrix is ")
   print(par_matrix) 
   estimat_weight <- par_matrix %*% t((select_mix_ob ))
   return (estimat_weight)
}




gene_select <- function(dataM, low, high)
{
   result <- list()
   dataM <- as.matrix(dataM)
   #cellName <- colnames(dataM)
   #i <-3 
   #print(sum(dataM[,i] > high))
   #print("rowSum is ")
   #print(sum((rowSums(dataM[,-i] < low) == (ncol(dataM) - 1))))
   #print(sum( (dataM[,i] > high) & (rowSums(dataM[,-i] < low) == (ncol(dataM) - 1))))
   for( i in 1 : ncol(dataM))
   {
       index <-  (dataM[,i] > high) & (rowSums(dataM[,-i] < low) == (ncol(dataM) - 1))
       result[[i]] <- index
       print("sum of index is ")
       print(sum(index))
       if( sum(index) == 0)
       {
          print("no specific gene for this cell type")
          print("cell type Id is ")
          print(i)
          stop()
       }
   }
   return (result)
}
	