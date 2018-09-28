#Sep04,2018
#Interface for EM model, input parameters at here,receive results out
#exp supposed to be csv. signatures list of signature genes,max_iter=1000,convergence_n=10,convergence_rate=0.9
setwd('/Users/zhangze/Desktop/SCINA_test')
source('EM_model.R')
source('EM_functions.R')
suppressMessages(library('gplots'))
args=commandArgs(trailingOnly = T)
exp_csv=args[1];signatures_csv=args[2];max_iter=args[3];convergence_n=args[4];convergence_rate=args[5];sensitivity_cutoff=args[6];rm_overlap=args[7];allow_unknown=args[8];output=args[9];job_id=args[10]#output:RData with path
max_iter=as.numeric(max_iter);convergence_n=as.numeric(convergence_n);convergence_rate=as.numeric(convergence_rate);sensitivity_cutoff=as.numeric(sensitivity_cutoff);rm_overlap=as.numeric(rm_overlap);allow_unknown=as.numeric(allow_unknown)
exp=read.csv(exp_csv,row.names=1,stringsAsFactors = F)
signatures=preprocess.signatures(signatures_csv)
results=SCINA(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,allow_unknown,log_file=paste(job_id,'SCINA.log',sep='_'))
jpeg(paste(job_id,'output_plot.jpg',sep='_'),width =2444,height =2444,res=600 )
par(mar=rep(3,4))
plotheat.SCINA(exp,results,signatures)
dev.off()
base_name=basename(output)
base_path=gsub(base_name,'',output)
save(results,file=paste(base_path,job_id,'_',base_name,sep=''))
