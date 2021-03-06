#######################################################################
#SCINA: A Semi-Supervised Category Identification and Assignment Tool.
#Maintainer: Ze Zhang, <ze.zhang@utsouthwestern.edu>
#######################################################################

#Description: 
#' An automatic cell type detection and assignment algorithm for single cell RNA-Seq (scRNA-seq) and Cytof/FACS data. 
#' SCINA is capable of assigning cell type identities to a pool of cells profiled by scRNA-Seq or Cytof/FACS data with prior knowledge of identifiers, 
#' such as genes and protein symbols,that are highly or lowly expressed in each category. 
#' See Ze Z, Danni L, et al (2018) for more details.

#Dependencies:
#' @import MASS
#' @importFrom gplots heatmap.2

#Parameters:
#' @param exp A normalized matrix representing the target dataset. Columns correpond to objects (cell barcodes for example), rows correspond to attributes or variables (gene symbols for example).
#' @param signatures A list contains multiple signature identifier lists. Each signature identifier list (genes for example) represents prior knowledge for one category (cell type for example), containing genes or protein symbols with high degree of detection.
#' @param max_iter An integer > 0. Default is 100. Max iterations allowed for the EM algorithm.
#' @param covergence_n An integer > 0. Default is 10. Stop the EM algorithm if during the last n rounds of iterations, cell type assignment keeps steady above the convergence_rate.
#' @param convergence_rate A float between 0 and 1. Default is 0.99. Percentage of cells for which the type assignment remains stable for the last n rounds. 
#' @param sensitivity_cutoff A float between 0 and 1. Default is 1. The cutoff to remove signatures whose cells types are deemed as non-existent at all in the data by the SCINA algorithm. 
#' @param rm_overlap A binary value, default 1 (TRUE), denotes that shared symbols between signature lists will be removed. If 0 (FALSE) then allows different cell types to share the same identifiers.
#' @param allow_unknown A binary value, default 1 (TRUE). If 0 (FALSE) then no cell will be assigned to the 'unknown' category.
#' @param log_file A string names the record of the running status of the SCINA algorithem, default 'SCINA.log'.

#Outputs:
#' @return cell_labels A vector contains cell type mapping results for each cell.
#' @return probabilities A probability matrix indicating the predicted probability for each cell belongs to each cell type respectively.

#Usage:
#' SCINA(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,allow_unknown,log_file)

SCINA=function(exp,signatures,max_iter=100,convergence_n=10,convergence_rate=0.99,
               sensitivity_cutoff=1,rm_overlap=1,allow_unknown=1,log_file='SCINA.log'){
  cat('Start running SCINA.',file=log_file,append=F)
  cat('\n',file=log_file,append=T)
  #Create a status file for the webserver
  status_file=paste(log_file,'status',sep='.')
  all_sig=unique(unlist(signatures))
  # Create low-expression signatures.
  invert_sigs=grep('^low_',all_sig,value=T)
  if(!identical(invert_sigs, character(0))){
    cat('Converting expression matrix for low_genes.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    invert_sigs_2add=unlist(lapply(invert_sigs,function(x) strsplit(x,'_')[[1]][2]))
    invert_sigs=invert_sigs[invert_sigs_2add%in%row.names(exp)]
    invert_sigs_2add=invert_sigs_2add[invert_sigs_2add%in%row.names(exp)]
    sub_exp=-exp[invert_sigs_2add,,drop=F]
    row.names(sub_exp)=invert_sigs
    exp=rbind(exp,sub_exp)
    rm(sub_exp,all_sig,invert_sigs,invert_sigs_2add)
  }
  # Check input parameters.
  quality=check.inputs(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,log_file)
  if(quality$qual==0){
    cat('EXITING due to invalid parameters.',file=log_file,append=T)
    cat('\n',file=log_file,append=T)
    cat('0',file=status_file,append=F)
    stop('SCINA stopped.')
  }
  signatures=quality$sig
  max_iter=quality$para[1]
  convergence_n=quality$para[2]
  convergence_rate=quality$para[3]
  sensitivity_cutoff=quality$para[4]
  # Initialize variables.
  exp=as.matrix(exp)
  exp=exp[unlist(signatures),,drop=F]
  labels=matrix(0,ncol=convergence_n, nrow=dim(exp)[2])
  unsatisfied=1
  if(allow_unknown==1){
    tao=rep(1/(length(signatures)+1),length(signatures))
  }else{tao=rep(1/(length(signatures)),length(signatures))}
  theta=list()
  for(i in 1:length(signatures)){
    theta[[i]]=list()
    theta[[i]]$mean=t(apply(exp[signatures[[i]],,drop=F],1,function(x) quantile(x,c(0.7,0.3))))
    tmp=apply(exp[signatures[[i]],,drop=F],1,var)
    theta[[i]]$sigma1=diag(tmp,ncol = length(tmp))
    theta[[i]]$sigma2=theta[[i]]$sigma1
  }
  for(marker_set in 1:length(theta)){
    if(is_empty(theta[[marker_set]]$sigma1) == TRUE){
      theta <- theta[-marker_set]
    }
  }
  sigma_min=min(sapply(theta,function(x) min(c(diag(x$sigma1),diag(x$sigma2)))))/100
  remove_times=0
  # Run SCINA algorithm.
  while(unsatisfied==1){
    prob_mat=matrix(tao,ncol=dim(exp)[2],nrow=length(tao))
    row.names(prob_mat)=names(signatures)
    iter=0
    labels_i=1
    remove_times=remove_times+1
    while(iter<max_iter){
      iter=iter+1
      # E step: estimate variables.
      for(i in 1:length(signatures)){
        theta[[i]]$inverse_sigma1=theta[[i]]$inverse_sigma2=chol2inv(chol(theta[[i]]$sigma1))
      }
      for (r in 1:dim(prob_mat)[1]){
        prob_mat[r,]=tao[r]*density_ratio(exp[signatures[[r]],,drop=F],theta[[r]]$mean[,1],
                                          theta[[r]]$mean[,2],theta[[r]]$inverse_sigma1,theta[[r]]$inverse_sigma2)
      }
      prob_mat=t(t(prob_mat)/(1-sum(tao)+colSums(prob_mat)))
      # M step: update sample distributions.
      tao=rowMeans(prob_mat)
      for(i in 1:length(signatures)){
        theta[[i]]$mean[,1]=(exp[signatures[[i]],]%*%prob_mat[i,])/sum(prob_mat[i,])
        theta[[i]]$mean[,2]=(exp[signatures[[i]],]%*%(1-prob_mat[i,]))/sum(1-prob_mat[i,])
        keep=theta[[i]]$mean[,1]<=theta[[i]]$mean[,2]
        if(any(keep)){
          theta[[i]]$mean[keep,1]=rowMeans(exp[signatures[[i]][keep],,drop=F])
          theta[[i]]$mean[keep,2]=theta[[i]]$mean[keep,1]
        }
        tmp1=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,1])^2)
        tmp2=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,2])^2)
        diag(theta[[i]]$sigma1)=diag(theta[[i]]$sigma2)=
          colSums(tmp1*prob_mat[i,]+tmp2*(1-prob_mat[i,]))/dim(prob_mat)[2]
        diag(theta[[i]]$sigma1)[diag(theta[[i]]$sigma1)<sigma_min]=sigma_min
        diag(theta[[i]]$sigma2)[diag(theta[[i]]$sigma2)<sigma_min]=sigma_min
      }
      labels[,labels_i]=apply(rbind(1-colSums(prob_mat),prob_mat),2,which.max)-1
      # Compare estimations with stop rules.
      if(mean(apply(labels,1,function(x) length(unique(x))==1))>=convergence_rate){
        cat('Job finished successfully.',file=log_file,append=T)
        cat('\n',file=log_file,append=T)
        cat('1',file=status_file,append=F)
        break
      }
      labels_i=labels_i+1
      if(labels_i==convergence_n+1){
        labels_i=1
      }
      if(iter==max_iter){
        cat('Maximum iterations, breaking out.',file=log_file,append=T)
        cat('\n',file=log_file,append=T)
      }
    }
    #Build result matrices.
    colnames(prob_mat)=colnames(exp)
    row.names(prob_mat)=names(signatures)
    row.names(labels)=colnames(exp)
    # Attempt to remove unused signatures. 
    dummytest=sapply(1:length(signatures),function(i) mean(theta[[i]]$mean[,1]-theta[[i]]$mean[,2]==0))
    if(all(dummytest<=sensitivity_cutoff)){
      unsatisfied=0
    }else{
      rev=which(dummytest>sensitivity_cutoff);
      cat(paste('Remove dummy signatures:',rev,sep=' '),file=log_file,append=T)
      cat('\n',file=log_file,append=T)
      signatures=signatures[-rev]
      tmp=1-sum(tao)
      tao=tao[-rev]
      tao=tao/(tmp+sum(tao))
      theta=theta[-rev]
    }
  }
  return(list(cell_labels=c("unknown",names(signatures))[1+labels[,labels_i]],probabilities=prob_mat))
}
