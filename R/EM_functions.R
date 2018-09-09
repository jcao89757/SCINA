EM_check=function(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap)
{
  # initialize parameters
  quality=1
  def_max_iter=1000
  def_conv_n=10
  def_conv_rate=0.99
  def_dummycut=0.33
  allgenes=row.names(exp)

  # clean exp
  if (any(is.na(exp)))
  {
    cat('NA exists in expression matrix.',file='SCINA.log',append=T,sep='\n')
    quality=0
  }

  # clean signatures
  if (any(is.na(signatures)))
  {
    cat('Null cell type signature genes.',file='SCINA.log',append=T,sep='\n');quality=0
  }else
  {
    #na.omit(NAs) in each signatures,remove unvalid genes (not in exp)
    signatures=sapply(signatures,function(x) unique(x[(!is.na(x)) & (x %in% allgenes)]),simplify=F)

    #remove duplicate genes
    if(rm_overlap==1)
    {
      tmp=table(unlist(signatures))
      signatures=sapply(signatures,function(x) x[x %in% names(tmp[tmp==1])],simplify=F)
    }
    #check if any genes have all 0 counts
    signatures=sapply(signatures,function(x) x[apply(exp[x,,drop=F],1,sd)>0],simplify=F)
  }
  # clean other parameters
  if (is.na(convergence_n))
  {
    cat('Using convergence_n=default',file='SCINA.log',append=T,sep='\n')
    convergence_n=def_conv_n
  }
  if (is.na(max_iter))
  {
    cat('Using max_iter=default',file='SCINA.log',append=T,sep='\n')
    max_iter=def_max_iter
  }else
  {
    if (max_iter<convergence_n)
    {
      cat('Using max_iter=default due to smaller than convergence_n.',file='SCINA.log',append=T,sep='\n')
      max_iter=def_max_iter # this theorectically, still doesn't ensure max_iter>convergence_n (Tao)
    }
  }

  if (is.na(convergence_rate))
  {
    cat('Using convergence_rate=default.',file='SCINA.log',append=T,sep='\n')
    convergence_rate=def_conv_rate
  }

  if (is.na(sensitivity_cutoff))
  {
    cat('Using sensitivity_cutoff=default.',file='SCINA.log',append=T,sep='\n')
    sensitivity_cutoff=def_dummycut
  }
  # return cleaned data
  return(list(qual=quality,sig=signatures,
    para=c(max_iter,convergence_n,convergence_rate,sensitivity_cutoff)))
}

density_ratio=function(e,mu1,mu2,inverse_sigma1,inverse_sigma2){
  tmp1=colSums((e-mu1)*(inverse_sigma1%*%(e-mu1)))
  tmp2=colSums((e-mu2)*(inverse_sigma2%*%(e-mu2)))
  tmp=exp(-1/2*(tmp1+log(1/det(inverse_sigma1))-tmp2-log(1/det(inverse_sigma2))))
  tmp[tmp>1e200]=1e200
  tmp[tmp<1e-200]=1e-200
  return(tmp)
}

#' Provide proper gene signatures depending on labeled cell population sequencing results.
#'
#'@param batched_mat #' A list of lists, contains all labeled cells, each cell type has an expression matrix with rows represent genes and columns represent cells. Each matrix is assigned in a daughter list with its cell type as the name of the list.
#'@param cells #' A list of cell types need to generate gene signatures from.
#'@param n #' An integer larger than 0, indicating the number of cells in each cell type for the purpose of signature generation.
#'@param cutoff #' A float larger than 0 to control signature gene amounts for each of the cell types. A cutoff close to 0 will results in more genes in each signature. Default 1.5.
#'@return signatures #' See more details in \code{\link{SCINA}}

findgenesigs=function(batched_mat,cells,n,cutoff=1.5){
  mean_cell=c()
  for (cell in cells){
    tmp=log(batched_mat[[cell]][,1:n]+1)
    mean_cell=cbind(mean_cell,rowMeans(tmp,na.rm = T))
  }
  signatures=list()
  for(i in 1:dim(mean_cell)[1]){
    cell_max=which.max(mean_cell[i,])
    if (mean_cell[i,cell_max,drop=F]>max(mean_cell[i,-cell_max,drop=F])+cutoff) {
      signatures[[cells[cell_max]]]=c(signatures[[cells[cell_max]]],rownames(mean_cell)[i])
    }
  }
  return(signatures)
}

#' A function to plot SCINA results in a heatmap.
#'
#' @param exp See more details in \code{\link{SCINA}}
#' @param results An output object returned from SCINA.
#' @param signatures See more details in \code{\link{SCINA}}
#' @return Plot A heatmap showing signature genes' expression level and SCINA predicted cell types.
#' @examples
#' plot.SCINA(exp,results,signatures)
#'
plot.SCINA=function(exp,results,signatures){
  col_row=topo.colors(length(signatures))
  col_col=cm.colors(length(unique(results$cell_labels)))
  exp2plot=exp[unlist(signatures),order(factor(results$cell_labels,levels=c(names(signatures),'unknown')),decreasing=F)]
  exp2plot=as.matrix(exp2plot)
  col_colside=col_col[as.factor(results$cell_labels[order(factor(results$cell_labels,levels=c(names(signatures),'unknown')),decreasing=F)])]
  col_rowside=col_row[unlist(sapply(1:length(signatures),function(i) rep(i, length(signatures[[i]]))))]
  par(mar=rep(1,4))
  plot=heatmap.2(exp2plot,trace='none',col=c('white','mistyrose1','lightpink1','lightpink2','lightpink3','lightpink4'),Rowv=FALSE,Colv=FALSE,dendrogram='none',xlab='Cells',ylab='Genes',labCol=NA,ColSideColors=col_colside,RowSideColors=col_rowside,key = F,margins = c(6.5,6.5))
  legend_text=c(paste('Gene identifiers',names(signatures),sep='_'),names(signatures),'unknown')
  legend_cor=c(col_row,unique(col_colside))
  legend('topright',legend=legend_text,fill=legend_cor,cex=0.5)
  return(invisible(plot))
}

#' A function to convert signatures uploaded via .csv files to lists used by SCINA.
#'
#' @param file_path The path where the .csv file stores. The first row of the file should be cell type names. Each column is occupied by the signature genes/protein markers for the cell type in the first row. Please find more details in \code{\link{SCINA}}.
#' @return A list of signature gene lists as an input for SCINA.
#' @examples
#' signatures=preprocess.signatures('./example/example_signatures.csv')
preprocess.signatures=function(file_path){
  csv_signatures=read.csv(file_path,stringsAsFactors = F,header = T,fill = F,quote ='')
  signatures=as.list(csv_signatures)
  return(signatures)
}
