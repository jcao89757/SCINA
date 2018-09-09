#' @title A semi-supervised category identification and assignment tool.
#'
#' @description An automatic cell type detection and assignment algorithm for single cell RNA-Seq (scRNA-seq) and Cytof/FACS data. SCINA is capable of assigning cell type identities to a pool of cells profiled by scRNA-Seq or Cytof/FACS data with prior knowledge of identifiers, such as genes and protein symbols,
#' that are highly or lowly expressed in each category. See Ze Z, Danni L, et al (2018) <doi:> for more details.
#'
#' @param exp A normalized matrix representing the target dataset. Columns correpond to objects (cell barcodes for example), rows correspond to attributes or variables (gene symbols for example).
#' @details 
#' For efficiency of data transfer and computation, the user is encouraged to upload the subset of the gene expression matrix that contains only the genes that appeared in the signature list.
#' @param signatures A list contains multiple signature identifier lists. Each signature identifier list (genes for example) represents prior knowledge for one category (cell type for example), containing genes or protein symbols with high degree of detection.
#' @details
#' For any symbols in signature lists, if the category is identified with symbol X's low detection level, please specify the symbol as 'low_X'. The name for the list is the category.
#' Details for 'low_X' (take scRNA-Seqs as an example):
#' (a) There are 4 cell types, the first one highly express one gene A, and the other three lowly express the same gene. Then it is better to specify A as the high marker for cell type 1, but it is not a good idea to specify A as the low expression marker for cell type 2,3,4.\cr
#' \if{html}{\figure{sub1.png}{options: width=0.5in}}
#' \if{latex}{\figure{sub1.png}{options: width=0.5in}}
#' (b) There are 4 cell types, the first one lowly express one gene A, and the other three highly express the same gene. Then is it better to specify A as the low marker for cell type 1, but it is not a good idea to specify A as the high expression marker for cell type 2,3,4.\cr
#' \if{html}{\figure{sub2.png}{options: width=0.5in}}
#' \if{latex}{\figure{sub2.png}{options: width=0.5in}}
#' (c) There are 4 cell types, the first one lowly express one gene A, the second and third one moderately express gene A, and the last one highly express gene A. Then is it better to specify A as the low marker for cell type 1, and as the high expression marker for cell type 4.\cr
#' \if{html}{\figure{sub3.png}{options: width=0.5in}}
#' \if{latex}{\figure{sub3.png}{options: width=0.5in}}
#' (d) The same specification can be applied to protein markers in CyTOF anlysis.\cr
#' \if{html}{\figure{sub4.png}{options: width=0.5in}}
#' \if{latex}{\figure{sub4.png}{options: width=0.5in}}
#' For immune cell types, an example signature lists can be found in \code{\link[DisHet]{eTME_signatures}.
#' @param max_iter An integer > 0. Default is 100. Max iterations allowed for the EM algorithm.
#' @param covergence_n An integer > 0. Default is 10. Stop the EM algorithm if during the last n rounds of iterations, cell type assignment keeps steady above the convergence_rate.
#' @param convergence_rate A float between 0 and 1. Default is 0.99. Percentage of cells for which the type assignment remains stable for the last n rounds. 
#' @param sensitivity_cutoff A float between 0 and 1. Default is 1. The cutoff to remove signatures whose cells types are deemed as non-existent at all in the data by the SCINA algorithm. 
#' @details
#' Small sensitivity_cutoff leads to more signatures to be removed, and 1 denotes that no signature is removed.
#' @param rm_overlap A binary value, default 1 (TRUE), denotes that shared symbols between signature lists will be removed. If 0 (FALSE) then allows different cell types to share the same identifiers.
#' @param allow_unknown A binary value, default 1 (TRUE). If 0 (FALSE) then no cell will be assigned to the 'unknown' category.
#'
#' @return cell_labels A vector contains cell type mapping results for each cell.
#' @return probabilities A probability matrix indicating the predicted probability for each cell belongs to each cell type respectively.
#' @import MASS
#' @importFrom gplots heatmap.2
#' @importFrom Rdpack reprompt
#' @export SCINA
#' @export plot.SCINA
#' @export preprocess.signatures
#' @examples

#' load(system.file("example/example_expmat.RData",package="SCINA"))
#' load(system.file("example/example_signatures.RData",package="SCINA"))
#' exp=exp_test$exp_data
#' results=SCINA(exp,signatures,max_iter=120,convergence_n=12,convergence_rate=0.999,sensitivity_cutoff=0.9)
#' table(exp_test$true_label,results$cell_labels)

#' @references
#' \insertRef{@Rpack:bibtex}{Rdpack}
#'
SCINA=function(exp,signatures,max_iter=100,convergence_n=10,convergence_rate=0.99,
               sensitivity_cutoff=1,rm_overlap=1,allow_unknown=1)
{
  cat('Start running SCINA.',file='SCINA.log',append=F)
  #check for low-expression signatures
  all_sig=unique(unlist(signatures))
  invert_sigs=grep('^low_',all_sig,value=T)
  if(!identical(invert_sigs, character(0))){
    cat('Converting expression matrix for low_genes.',file='SCINA.log',append=T)
    invert_sigs_2add=unlist(lapply(invert_sigs,function(x) strsplit(x,'_')[[1]][2]))
    invert_sigs=invert_sigs[invert_sigs_2add%in%row.names(exp)]# To avoid users adding 'low_' to non-existing genes.
    invert_sigs_2add=invert_sigs_2add[invert_sigs_2add%in%row.names(exp)]
    sub_exp=-exp[invert_sigs_2add,,drop=F]
    row.names(sub_exp)=invert_sigs
    exp=rbind(exp,sub_exp)
    rm(sub_exp,all_sig,invert_sigs,invert_sigs_2add)
  }
  #invalid parameters check
  quality=EM_check(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap)
  if(quality$qual==0){cat('EXITING due to invalid parameters.',file='SCINA.log',append=T,sep='\n');stop('SCINA stopped.')}
  signatures=quality$sig
  max_iter=quality$para[1]
  convergence_n=quality$para[2]
  convergence_rate=quality$para[3]
  sensitivity_cutoff=quality$para[4]
  # initialize variables
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
  sigma_min=min(sapply(theta,function(x) min(c(diag(x$sigma1),diag(x$sigma2)))))/100
  remove_times=0
  # SCINA
  while(unsatisfied==1)
  {
    prob_mat=matrix(tao,ncol=dim(exp)[2],nrow=length(tao))
    row.names(prob_mat)=names(signatures)
    iter=0
    labels_i=1
    remove_times=remove_times+1

    while(iter<max_iter)
    {
      iter=iter+1
      #2. E_step
      for(i in 1:length(signatures))
      {
        theta[[i]]$inverse_sigma1=theta[[i]]$inverse_sigma2=chol2inv(chol(theta[[i]]$sigma1))
      }
      for (r in 1:dim(prob_mat)[1])
      {
        prob_mat[r,]=tao[r]*density_ratio(exp[signatures[[r]],,drop=F],theta[[r]]$mean[,1],
                                          theta[[r]]$mean[,2],theta[[r]]$inverse_sigma1,theta[[r]]$inverse_sigma2)
      }
      prob_mat=t(t(prob_mat)/(1-sum(tao)+colSums(prob_mat)))
      #3. M-step
      tao=rowMeans(prob_mat)
      for(i in 1:length(signatures)){
        #update means
        theta[[i]]$mean[,1]=(exp[signatures[[i]],]%*%prob_mat[i,])/sum(prob_mat[i,])
        theta[[i]]$mean[,2]=(exp[signatures[[i]],]%*%(1-prob_mat[i,]))/sum(1-prob_mat[i,])
        #remain u1>u2
        keep=theta[[i]]$mean[,1]<=theta[[i]]$mean[,2]
        if(any(keep)){
          theta[[i]]$mean[keep,1]=rowMeans(exp[signatures[[i]][keep],,drop=F])
          theta[[i]]$mean[keep,2]=theta[[i]]$mean[keep,1]
        }
        #update covars
        tmp1=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,1])^2)
        tmp2=t((exp[signatures[[i]],,drop=F]-theta[[i]]$mean[,2])^2)
        diag(theta[[i]]$sigma1)=diag(theta[[i]]$sigma2)=
          colSums(tmp1*prob_mat[i,]+tmp2*(1-prob_mat[i,]))/dim(prob_mat)[2]
        if (sum(prob_mat[i,])==0) {theta[[i]]$sigma1=theta[[i]]$sigma2}
        diag(theta[[i]]$sigma1)[diag(theta[[i]]$sigma1)<sigma_min]=sigma_min
        diag(theta[[i]]$sigma2)[diag(theta[[i]]$sigma2)<sigma_min]=sigma_min
      }

      labels[,labels_i]=apply(rbind(1-colSums(prob_mat),prob_mat),2,which.max)-1
      if(mean(apply(labels,1,function(x) length(unique(x))==1))>=convergence_rate)
      {cat('Job finished successfully.',file='SCINA.log',append=T,sep='\n');break}
      labels_i=labels_i+1
      if(labels_i==convergence_n+1){labels_i=1}
      if(iter==max_iter){cat('Maximum iterations, breaking out.',file='SCINA.log',append=T,sep='\n')}
    }
    colnames(prob_mat)=colnames(exp)
    row.names(prob_mat)=names(signatures)
    row.names(labels)=colnames(exp)
    #4. test dummy signature lists
    dummytest=sapply(1:length(signatures),function(i) mean(theta[[i]]$mean[,1]-theta[[i]]$mean[,2]==0))
    if(all(dummytest<=sensitivity_cutoff))
    {
      unsatisfied=0
    }else
    {
      rev=which(dummytest>sensitivity_cutoff);
      cat(paste('Remove dummy signatures:',rev,sep=' '),file='SCINA.log',append=T,sep='\n');
      # note the changes below
      signatures=signatures[-rev]
      tmp=1-sum(tao)
      tao=tao[-rev]
      tao=tao/(tmp+sum(tao))
      theta=theta[-rev]
    }
  }
  return(list(cell_labels=c("unknown",names(signatures))[1+labels[,labels_i]],probabilities=prob_mat))
}
