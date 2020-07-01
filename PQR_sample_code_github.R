library(xts)
library(quantreg)
library(rqpd)
library(fBasics)
library(compiler)
library(parallel)
enableJIT(3)

load("datarv_xts_sample.RData")
load("dataret_xts_sample.RData")
#  return ~ RV --------
length_insample<-1001
qq <- seq(0.05, 0.95, 0.05)

RV_coef=lapply(1:ncol(dataret_xts),
               function (i) {asset=i;
               y=dataret_xts[2:length_insample,asset];
               RV=sqrt(datarv_xts[1:length_insample-1,asset]);
               coef(rq((y)~(RV),tau=qq))})


# save coefficients #
ret_RV_tau=sapply(1:length(qq), 
                  function(j) sapply(1:ncol(dataret_xts),
                                     function(i) RV_coef[[i]][2,j]))
colnames(ret_RV_tau) <- qq

#  Panel Quantile Regressions RV  
m <- length_insample-1
n <- ncol(dataret_xts)
s <- as.factor(rep(1:n,rep(m,n)))

ret=vec(vec(dataret_xts[2:length_insample,1:n]))
RV=vec(vec(sqrt(datarv_xts[1:(length_insample-1),1:n])))

# panel QR ret ~ RV 
q_panel<- seq(0.05, 0.95, 0.05)
panel_fit_ret_RV=lapply(1:length(q_panel), 
                        function(i) rqpd(ret ~ RV | s, 
                                         panel(taus=q_panel[i], tauw=rep(1/length(q_panel[i]), length(q_panel[1])), lambda=0)))
names(panel_fit_ret_RV)<-q_panel

# save estimated coef
ret_RV_panel_tau=sapply(1:length(q_panel),
                        function(i) panel_fit_ret_RV[[i]]$coefficients[1])



# bootstrap 
cl <- makeCluster(detectCores())
clusterExport(cl, list("boot.rqpd", "panel_fit_ret_RV"))
system.time(
  bootstrap_QR_ret_RV<-parLapply(cl, 1:length(q_panel),
                                 function(i) boot.rqpd(panel_fit_ret_RV[[i]]$ids,panel_fit_ret_RV[[i]]$X,
                                                       panel_fit_ret_RV[[i]]$Z, panel_fit_ret_RV[[i]]$y, panel_fit_ret_RV[[i]]$panel,
                                                       panel_fit_ret_RV[[i]]$control, R=500, bsmethod = "wxy"))
)
stopCluster(cl)



n_row<-length(panel_fit_ret_RV$`0.05`$coefficients)
n_col<-length(q_panel)

coef_panel<-matrix(nrow = n_row,ncol = n_col)
coef_panel<-sapply(1:length(q_panel),
                   function(i) panel_fit_ret_RV[[i]]$coef)
colnames(coef_panel)<-as.character(c(q_panel))


boot_SE<- matrix(nrow = n_row,ncol = n_col)
boot_SE<-sapply(1:length(q_panel),
                function(i)sqrt(diag(cov(bootstrap_QR_ret_RV[[i]]))))
row.names(boot_SE)<-row.names(coef_panel)
colnames(boot_SE)<-as.character(c(q_panel))

t_stat<- matrix(nrow = n_row,ncol = n_col)
t_stat<- matrix(paste("(",round(coef_panel/boot_SE,2),")",sep=""),ncol = n_col)
row.names(t_stat)<-row.names(coef_panel)
colnames(t_stat)<-as.character(c(q_panel))


summary_table <- matrix(nrow=2*n_row,ncol=n_col)
for (i in 1:n_row) {
  summary_table[2*i-1,]<-round(coef_panel[i,],2)
  summary_table[2*i,]<-t_stat[i,]
}
row.names(summary_table)<-rep(c(row.names(t_stat)), times = c(rep(2,n_row)))
colnames(summary_table)<-as.character(c(q_panel))
summary_table= data.frame(row = rownames(summary_table),data.frame(summary_table))
colnames(summary_table)=as.character(c(" ",q_panel))

keeps <- c(" ","0.05", "0.1","0.25","0.5","0.75","0.9","0.95") # quantiles to keep
rows_to_keep<-12 #number of rows to keep for latex table
summary_table_RV<-summary_table[1:rows_to_keep,keeps]

##### 2,5-97,5% conficence intervals
CI_QR_ret_RV<-list(RV=NULL)
CI_QR_ret_RV$RV=sapply(1:length(q_panel), function(i) quantile(bootstrap_QR_ret_RV[[i]][,1],c(0.975,0.025)))
# ### PLOT RV 
boxplot(ret_RV_tau,border= "dimgray", ylim=c(-2.5,2), xlab=expression(Quantiles),ylab=expression(Coefficients),
        cex.lab=1.65, cex.axis=1.65, cex.main=1.5, cex.sub=1.5)
lines(seq(1,19,1),ret_RV_panel_tau,col="black",lwd=4)
lines(CI_QR_ret_RV$RV[1,],type = "l",col= "black",lty = 2, lwd=2)
lines(CI_QR_ret_RV$RV[2,],type = "l",col= "black",lty = 2, lwd=2)
