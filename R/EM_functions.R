############################################################
#Internal functions and methods for SCINA
#Maintainer: Ze Zhang, <ze.zhang@utsouthwestern.edu>
############################################################

#Dependencies:
#' @import MASS
#' @importFrom gplots heatmap.2

#Description: 
#' An internal function called by SCINA.
#' Checking input parameters' formats, integrity, and whether their values are within the designed value range.

#Parameters:
#' See more details in EM_model.R

#Outputs:
#' @return qual A binary value, 0 (FALSE) represents the SCINA is not able to use this input sets. 1 (TRUE) represents that the inputs' quality is satisfied.
#' @return sig Signatures treated up to user requirements.
#' @return para Other parameters, marked in log_file if they are out of range and are changed to defaults.

check.inputs=function(exp, signatures, max_iter, convergence_n, convergence_rate, sensitivity_cutoff, rm_overlap, log_file)
{
  # Initialize parameters.
  quality=1
  def_max_iter=1000
  def_conv_n=10
  def_conv_rate=0.99
  def_dummycut=0.33
  allgenes=row.names(exp)
  # Check sequence matrices.
  if (any(is.na(exp))){
    cat('NA exists in expression matrix.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    quality=0
  }
  # Check signatures.
  if (any(is.na(signatures))){
    cat('Null cell type signature genes.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    quality=0
  }else{
    signatures=sapply(signatures,function(x) unique(x[(!is.na(x)) & (x %in% allgenes)]),simplify=F)
    # Remove duplicate genes.
    if(rm_overlap==1){
      tmp=table(unlist(signatures))
      signatures=sapply(signatures,function(x) x[x %in% names(tmp[tmp==1])],simplify=F)
    }
    # Check if any genes have all 0 counts
    signatures=sapply(signatures,function(x) x[apply(exp[x,,drop=F],1,sd)>0],simplify=F)
  }
  # Clean other parameters.
  if (is.na(convergence_n)){
    cat('Using convergence_n=default',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    convergence_n=def_conv_n
  }
  if (is.na(max_iter)){
    cat('Using max_iter=default',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    max_iter=def_max_iter
  }else{
    if (max_iter<convergence_n){
      cat('Using max_iter=default due to smaller than convergence_n.',file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      max_iter=convergence_n 
    }
  }
  if (is.na(convergence_rate)){
    cat('Using convergence_rate=default.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    convergence_rate=def_conv_rate
  }
  if (is.na(sensitivity_cutoff)){
    cat('Using sensitivity_cutoff=default.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    sensitivity_cutoff=def_dummycut
  }
  # Return cleaned parameters.
  return(list(qual=quality,sig=signatures,
              para=c(max_iter,convergence_n,convergence_rate,sensitivity_cutoff)))
}
 
                      
#Description:
#' An internal function called by SCINA.
#' Designed to calculate the probability matrices for each cell type from their distribution parameters.
                      
density_ratio=function(e,mu1,mu2,inverse_sigma1,inverse_sigma2){
  tmp1=colSums((e-mu1)*(inverse_sigma1%*%(e-mu1)))
  tmp2=colSums((e-mu2)*(inverse_sigma2%*%(e-mu2)))
  tmp=exp(-1/2*(tmp1+log(1/det(inverse_sigma1))-tmp2-log(1/det(inverse_sigma2))))
  tmp[tmp>1e200]=1e200
  tmp[tmp<1e-200]=1e-200
  return(tmp)
}

                      
#Description:
#' A function to plot SCINA results in a heatmap.
                      
#Parameters:
#' @param exp See more details in \code{\link{SCINA}}
#' @param results An output object returned from SCINA.
#' @param signatures See more details in \code{\link{SCINA}}
                      
#Output:
#' @return A heatmap showing signature genes' expression level and SCINA predicted cell types.

#Usage:
#' plotheat.SCINA(exp,results,signatures)

plotheat.SCINA=function(exp,results,signatures){
  # Remove nonexist signature genes.
  allgenes=row.names(exp)
  signatures=sapply(signatures,function(x) unique(x[(!is.na(x)) & (x %in% allgenes)]),simplify=F)
  # Build side color bars.
  col_row=topo.colors(length(signatures))
  col_col=cm.colors(length(unique(results$cell_labels)))
  # Build matrices to heatmap.2 function.
  exp2plot=exp[unlist(signatures),order(factor(results$cell_labels,levels=c(names(signatures),'unknown')),decreasing=F)]
  exp2plot=as.matrix(exp2plot)
  col_colside=col_col[as.factor(results$cell_labels[order(factor(results$cell_labels,levels=c(names(signatures),'unknown')),decreasing=F)])]
  col_rowside=col_row[unlist(sapply(1:length(signatures),function(i) rep(i, length(signatures[[i]]))))]
  # Plot heatmaps.
  par(mar=rep(1,4))
  plot=heatmap.2(exp2plot,trace='none',col=c('white','mistyrose1','lightpink1','lightpink2','lightpink3','lightpink4'),Rowv=FALSE,Colv=FALSE,dendrogram='none',xlab='Cells',ylab='Genes',labCol=NA,ColSideColors=col_colside,RowSideColors=col_rowside,key = F,margins = c(6.5,6.5))
  legend_text=c(paste('Gene identifiers',names(signatures),sep='_'),names(signatures),'unknown')
  legend_cor=c(col_row,unique(col_colside))
  legend('topright',legend=legend_text,fill=legend_cor,cex=0.5)
  return(invisible(plot))
}


#Description:                                    
#' A function to convert signatures uploaded via .csv files to lists used by SCINA.

#Parameter:
#' @param file_path The path where the .csv file stores. The first row of the file should be cell type names. Each column is occupied by the signature genes/protein markers for the cell type in the first row. Please find more details in \code{\link{SCINA}}.

#Output: 
#' @return signatures A list of signature gene lists as an input for SCINA.

#Usage:                                    
#' signatures=preprocess.signatures('./example/example_signatures.csv')
                                    
preprocess.signatures=function(file_path){
  csv_signatures=read.csv(file_path,stringsAsFactors = F,header = T,fill = F,quote ='')
  signatures=as.list(csv_signatures)
  signatures=sapply(signatures,function(x) x[(!is.na(x)) & x!=""],simplify = F)
  return(signatures)
}

