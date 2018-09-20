SCINA=function(exp,signatures,max_iter=100,convergence_n=10,convergence_rate=0.99,
               sensitivity_cutoff=1,rm_overlap=1,allow_unknown=1,log_file='SCINA.log')
  {
  cat('Start running SCINA.',file=log_file,append=F)
  all_sig=unique(unlist(signatures))
  invert_sigs=grep('^low_',all_sig,value=T)
  if(!identical(invert_sigs, character(0))){
    cat('Converting expression matrix for low_genes.',file=log_file,append=T)
    invert_sigs_2add=unlist(lapply(invert_sigs,function(x) strsplit(x,'_')[[1]][2]))
    invert_sigs=invert_sigs[invert_sigs_2add%in%row.names(exp)]# To avoid users adding 'low_' to non-existing genes.
    invert_sigs_2add=invert_sigs_2add[invert_sigs_2add%in%row.names(exp)]
    sub_exp=-exp[invert_sigs_2add,,drop=F]
    row.names(sub_exp)=invert_sigs
    exp=rbind(exp,sub_exp)
    rm(sub_exp,all_sig,invert_sigs,invert_sigs_2add)
  }
  quality=check.inputs(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,log_file)
  if(quality$qual==0){cat('EXITING due to invalid parameters.',file=log_file,append=T,sep='\n');stop('SCINA stopped.')}
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
        diag(theta[[i]]$sigma1)[diag(theta[[i]]$sigma1)<sigma_min]=sigma_min
        diag(theta[[i]]$sigma2)[diag(theta[[i]]$sigma2)<sigma_min]=sigma_min
      }

      labels[,labels_i]=apply(rbind(1-colSums(prob_mat),prob_mat),2,which.max)-1
      if(mean(apply(labels,1,function(x) length(unique(x))==1))>=convergence_rate)
      {cat('Job finished successfully.',file=log_file,append=T,sep='\n');break}
      labels_i=labels_i+1
      if(labels_i==convergence_n+1){labels_i=1}
      if(iter==max_iter){cat('Maximum iterations, breaking out.',file=log_file,append=T,sep='\n')}
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
      cat(paste('Remove dummy signatures:',rev,sep=' '),file=log_file,append=T,sep='\n');
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
