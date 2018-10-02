###########################
#Sep04,2018
#Maintainer: Ze Zhang 
###########################

#@' Description: Shell command line tool for EM model, wrap up input parameters to the script, receive results and saved in the system.

#@' Parameters: See more in EM_model.R. The parameter ranges are the same, exp and signatures are supposed to be .csv files.
#@' Parameters: output The RData file to store the output results and plots. Paths may be included.
#@' Parameters: job_id Any string could work. Used to avoid overwritting.

#@' Usage: Rscript EM_interface.R example_expmat_test.csv example_signatures.csv 100 10 0.99 1 1 1 test.RData 1


setwd('$HOME/SCINA_test')
# A test folder on Ze Zhang's system. Contains all input .csv files and scripts.
source('EM_model.R')
source('EM_functions.R')
suppressMessages(library('gplots'))
#Initiate input parameters
args=commandArgs(trailingOnly = T)
exp_csv=args[1];signatures_csv=args[2];max_iter=args[3];convergence_n=args[4];convergence_rate=args[5]
sensitivity_cutoff=args[6];rm_overlap=args[7];allow_unknown=args[8];output=args[9];job_id=args[10]
max_iter=as.numeric(max_iter);convergence_n=as.numeric(convergence_n);convergence_rate=as.numeric(convergence_rate)
sensitivity_cutoff=as.numeric(sensitivity_cutoff);rm_overlap=as.numeric(rm_overlap);allow_unknown=as.numeric(allow_unknown)
#Read data
exp=read.csv(exp_csv,row.names=1,stringsAsFactors = F)
signatures=preprocess.signatures(signatures_csv)
#Run SCINA
results=SCINA(exp,signatures,max_iter,convergence_n,convergence_rate,sensitivity_cutoff,rm_overlap,allow_unknown,log_file=paste(job_id,'SCINA.log',sep='_'))
#Draw output figures
jpeg(paste(job_id,'output_plot.jpg',sep='_'),width =2444,height =2444,res=600 )
par(mar=rep(3,4))
plotheat.SCINA(exp,results,signatures)
dev.off()
#Save results
base_name=basename(output)
base_path=gsub(base_name,'',output)
save(results,file=paste(base_path,job_id,'_',base_name,sep=''))
