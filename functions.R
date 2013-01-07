## this is the main function
## mix_file_name is the mixture array signal
## cell_gene_table is a table key mapping between cell type and cell specific genes
## weight_file_name is the file contain the weight matrix, if the weight matrix is offered by user, otherwise 
## this weight matrix will be automatically estimated by program.
## method: we have lm method and ml and QP (quadratic programming with constrain on the esitmated parameter
## out_cell_fillName which contain the 
Rf_main <- function(mix_file_name, cell_gene_table_name, weight_file_name=NULL, method = "Lm", out_cell_fileName, out_weight_fileName, log2 = T)
{
   
   #source("deconvolution.R")
   #source("three_cell_weight.R")
   
   mix_file <- as.matrix(read.delim(file = mix_file_name, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE))
   if(log2 == FALSE)
      mix_file <- log2(mix_file) 
   #return (mix_file)
   cell_gene_table <- as.matrix(read.table(file = cell_gene_table_name, header = F, sep = "\t", quote="\"", dec=".", fill = TRUE))
   if(ncol(cell_gene_table) != 2)
   {
      print("the cell_gene_table should be a table with only two columns")
      stop()
   }
   #nrnum <- nrow(cell_gene_table)
   #ncnum <- ncol(cell_gene_table)
   #cell_gene_table <- as.matrix(as.character(cell_gene_table), nrnum, ncnum)
   unique_cell_type <- unique(cell_gene_table[,2]) 
   gene_list <- list()
   for( i in 1 : length(unique_cell_type))
   {
      gene_list[[i]] <- cell_gene_table[cell_gene_table[,2] == unique_cell_type[i],1]
   } 
   names(gene_list) <- unique_cell_type
   if(is.null(weight_file_name))
   {
      estimate_weight <- weight_estimate_1(2^mix_file, gene_list, method)
      rownames(estimate_weight) <- names(gene_list)
   }
   else
      ## this function I read the weight matrix from file
   {
      estimate_weight <- (as.matrix(read.delim(file = weight_file_name, header = TRUE, sep = "\t", fill = TRUE)))
      estimate_weight <- estimate_weight/colSums(estimate_weight)
   }
   ## this step I need the deconvoluted signals
   ten_cell_decon_QP_LM <- deconvoltion_general(2^mix_file,t(estimate_weight), method )
   rownames(ten_cell_decon_QP_LM) <- rownames(mix_file)
   colnames(ten_cell_decon_QP_LM) <- unique_cell_type
   write.table(log2(ten_cell_decon_QP_LM), file = out_cell_fileName, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
               col.names = TRUE)
   write.table(estimate_weight, file = out_weight_fileName, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
               col.names = TRUE)  
   return (list(estimate_weight,log2(ten_cell_decon_QP_LM)))
}




# Lm is the linear space regression
# the input is anti-log intensity data
# lg is the log transformed space regression and take anti-log return value
deconvoltion_general <- function(data, weight, method = "LM")
{
   i <- 1
   data <- as.matrix(data)
   weight <- as.matrix(weight)
   paraM <- matrix(0, nrow(data), ncol(weight))
   
   if( method == "LM")
   {
      for(i in 1 : nrow(data))
      {
         y <- data[i,]
         lmOb <- lm(y ~ -1 + weight)
         paraM[i,] <- lmOb$coefficients
      }
   }
   else if (method == "QP_LM")
   {
      paraM <- GSM_QP(data, weight, l = min(data), u = max(data) , meq =0) 
   }
   else if (method == "LG")
   {
      for(i in 1 : nrow(data))
      {
         y <- log2(data[i,])
         lmOb <- lm(y ~ -1 + weight)
         paraM[i,] <- 2^lmOb$coefficients
      }
   }
   else if ( method == "QP_LG")
   {
      data <- log2(data)
      paraM <- GSM_QP(data, weight, l = min(data), u = max(data) , meq =0) 
      paraM <- 2^paraM
   }
   else
   {
       print("method notation is wrong")
       stop()
   }
   return (paraM)
}


weight_estimate_missing <- function(mix_ob,gene_list,method="Lm",topN=10){
	
	
	tmp=mix_ob[unlist(gene_list),]
	tmp=(tmp-rowMeans(tmp))/sd(t(tmp))
	tmp=colMeans(tmp)
	
	tt=c()
   for(i in 1:nrow(mix_ob)){
	   tt=c(tt,cor(tmp,mix_ob[i,]))
	}

	tti=sort(tt,index.return=TRUE)
   
	
	
	gene_list[[length(gene_list)+1]]=rownames(mix_ob)[tti$ix[1:10]]		
	res=weight_estimate_1(mix_ob,gene_list,method="Lm")
	
	return(res);
	
	}





#
#weight_estimate_missing <- function(mix_ob,gene_list,method="Lm"){
#	
#	gene_list[[length(gene_list)+1]]=rep(FALSE,nrow(mix_ob))
#		
#	res <- c()
#	
#	for(i in 1:nrow(mix_ob)){
#		
#		print(i)
#		
#		gene_list[[length(gene_list)]][i]=TRUE
#		
#		#return(list(mix_ob,gene_list))
#		
#		res <- c(res,weight_estimate_2(mix_ob,gene_list,method="Lm"))
#		
#		gene_list[[length(gene_list)]][i]=FALSE
#
#		
#	}
#	
#	
#	return(res);
#	
#	
#	}




## this function I will select the number of specific gene for each cell type, by the proportion
## parameter Qp_linear is a logical variable with false linear model, true as QP model
## this estimate is better
weight_estimate_1 <- function(mix_ob, gene_list, method = "Lm")
{
   select_mix_ob <- matrix()
   for ( i in 1 : length(gene_list))
   {
       print("the length of gene_list is ")
       print(length(gene_list[[i]]))
       #print("the length of third one")
       #print(length(gene_list[[3]]))
       #print(colMeans( mix_ob[gene_list[[i]], ]))
       if(i == 1)
       {
          if(length(gene_list[[i]]) == 1)
            select_mix_ob <- as.matrix(mix_ob[gene_list[[i]], ])
          else
            select_mix_ob <-  as.matrix(colMeans( mix_ob[gene_list[[i]], ]))
       }
       else
       {
          if(length(gene_list[[i]]) == 1)
            select_mix_ob <- cbind(select_mix_ob, as.matrix( mix_ob[gene_list[[i]], ]))
          else
            select_mix_ob <- cbind(select_mix_ob, as.matrix(colMeans( mix_ob[gene_list[[i]], ])))
       }
   }
   #print("dim of select_mix_ob is ")
   #print(dim(select_mix_ob))
   #print(head(select_mix_ob))
   y <- rep(1, times = nrow(select_mix_ob))
   #print("dim of select_mix_ob is ")
   #print(dim(select_mix_ob))
   #y <- matrix(1, nrow(select_mix_ob),ncol(select_mix_ob))
   #print("y is ")
   #print(y)
   #print("the length of y is ")
   #print(length(y))
   #stop()
   b_par <- numeric()
   ## we will use linear fit or QP method to get the estimator of paramter
   if(method != "Lm" && method != "QP_LM")
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
      #y <- t(matrix(1, nrow(select_mix_ob),ncol(select_mix_ob)))
      y <- rep(1, times = nrow(select_mix_ob))
      #print("before call QP function ")
      #print("the length of y is ")
      #print(length(y))
      #print(y)
      #print("dimentionof select_mix_ob is ")
      #print(dim(select_mix_ob))
      y <- t(as.matrix(y))
      #result <- list(y, select_mix_ob)
      #return (result)
      Qp <- GSM_QP(y, (select_mix_ob), 0.0, 2^34)
      print("Qp is ")
      print(Qp) 
      #b_par <- Qp * 10000
      b_par <- Qp 
      print(b_par)
      #stop()
   }
   len <- as.integer(length(gene_list))
   
  # return (list(b_par, len))
   par_matrix <- diag(c(b_par), len, len) 
   #estimat_weight <- par_matrix %*% t((select_mix_ob/10000 ))
   estimat_weight <- par_matrix %*% t((select_mix_ob ))
   return (estimat_weight)
}



## this function I will select the number of specific gene for each cell type, by the proportion
## parameter Qp_linear is a logical variable with false linear model, true as QP model
## this estimate is better
weight_estimate_2 <- function(mix_ob, gene_list, method = "Lm")
{
   select_mix_ob <- matrix()
   for ( i in 1 : length(gene_list))
   {
       #print("the length of gene_list is ")
       #print(sum(gene_list[[i]]))
       #print("the length of third one")
       #print(length(gene_list[[3]]))
       #print(colMeans( mix_ob[gene_list[[i]], ]))
       if(i == 1)
       {
          if(length(gene_list[[i]]) == 1)
            select_mix_ob <- as.matrix(mix_ob[gene_list[[i]], ])
          else
            select_mix_ob <-  as.matrix(colMeans( mix_ob[gene_list[[i]], ]))
       }
       else
       {
          if(length(gene_list[[i]]) == 1)
            select_mix_ob <- cbind(select_mix_ob, as.matrix( mix_ob[gene_list[[i]], ]))
          else
            select_mix_ob <- cbind(select_mix_ob, as.matrix(colMeans( mix_ob[gene_list[[i]], ])))
       }
   }
   #print("dim of select_mix_ob is ")
   #print(dim(select_mix_ob))
   #print(head(select_mix_ob))
   y <- rep(1, times = nrow(select_mix_ob))
   #print("dim of select_mix_ob is ")
   #print(dim(select_mix_ob))
   #y <- matrix(1, nrow(select_mix_ob),ncol(select_mix_ob))
   #print("y is ")
   #print(y)
   #print("the length of y is ")
   #print(length(y))
   #stop()
   b_par <- numeric()
   mse <- numeric()
   ## we will use linear fit or QP method to get the estimator of paramter
   if(method != "Lm" && method != "QP_LM")
   {
      print("the 'method' paramter is wrong, should be 'Lm' or 'Qp'")
      stop()
   }
   if(method == "Lm")
   {
      lmob <- lm( y ~  -1 + (select_mix_ob))
      
      #print(str(lmob))
      mse <- sum(lmob$residuals^2)
      b_par <- coef(lmob)
   }
   else
   {
      #y <- t(matrix(1, nrow(select_mix_ob),ncol(select_mix_ob)))
      y <- rep(1, times = nrow(select_mix_ob))
      #print("before call QP function ")
      #print("the length of y is ")
      #print(length(y))
      #print(y)
      #print("dimentionof select_mix_ob is ")
      #print(dim(select_mix_ob))
      y <- t(as.matrix(y))
      #result <- list(y, select_mix_ob)
      #return (result)
      Qp <- GSM_QP(y, (select_mix_ob), 0.0, 2^34)
      print("Qp is ")
      print(Qp) 
      #b_par <- Qp * 10000
      b_par <- Qp 
      print(b_par)
      #stop()
   }
   len <- as.integer(length(gene_list))
   
  # return (list(b_par, len))
   par_matrix <- diag(c(b_par), len, len) 
   #estimat_weight <- par_matrix %*% t((select_mix_ob/10000 ))
   estimat_weight <- par_matrix %*% t((select_mix_ob ))
   return (data.frame(estimat_weight,mse))
}

# this function is about the quadratic programming
GSM_QP <- function(ob, weight, l, u, meq =0)
{
  require("quadprog")
  sol= c()
  ob <- as.matrix(ob)
  weight <-as.matrix(weight) 
  for (id in 1:nrow(ob))
  {
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
      if(meq>0)
      {
          Amat=cbind(rep(1,numC),Amat)
          bvec=c(1,bvec)
      }
      sol=rbind(sol,solve.QP(Dmat,dvec,Amat,bvec=bvec,meq =meq)$solution)
   }
   return (sol)
}

weight_plot <- function(estimat_weight, weight, prefix)
{
   print("dim of estimat_weight is ")
   print((estimat_weight))
   print("the weight is ")
   print(dim(weight))
   print(weight)
   colorVec <- c("red", "blue", "orange", "green","purple")
   pdf(paste(prefix,"_","weight_estimate.pdf", sep = ""), height = 4.5, width = 4.5)
   #plot(estimat_weight[3,], weight[3,], xlab = "Estimated Weight", col = "yellow")
   #dev.off()
   for(i in 1 : nrow(estimat_weight))
   {
      if ( i == 1)
      {
         plot(weight[i,], estimat_weight[i,],  xlim = c(0, 1), ylim = c(0, 1), ylab = "Estimated Weight",xlab = "True Weight", col = colorVec[i], pch = 20, cex = 1.5) 

      }
      else
      {
         points(weight[i,], estimat_weight[i,], col = colorVec[i], pch = 20, cex = 1.5)
      } 
   }
   abline(0,1)
   dev.off()
   
} 
   
## this function is used for the plot 
GSM_decon_plot <- function(decon, sources, names, prefix)
{
  
    for( i in 1 : length(names))
    {
       print(paste(prefix, names[i], "-", names[i], ".pdf", sep = ""))
       pdf(paste(prefix, names[i], "-", names[i], ".pdf"))
       cov_eff <- cor(sources[,i], decon[,i])
       print("the cov_eff is ")
       print(cov_eff)
       my_title <- paste(c("R-Squared = "), format(cov_eff, width=4), sep = "")
       plot(sources[,i], decon[,i], xlab = paste("Pure-",names[i], sep = ""), ylab =  paste("Decon-",names[i], sep = ""), main = my_title)
       abline(0,1)
       dev.off()
    }
}





