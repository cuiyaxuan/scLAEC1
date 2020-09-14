# scLAEC1  
## A local affinity enhancement clustering algorithm for cell typing using single-cell RNA-seq data.  
#The program need to install parallel,SNFtool,apcluster package.  
install.packages(parallel)  
install.packages(SNFtool)  
install.packages(apcluster)  
library(scLAEC1)  
#Import the dataset  
Express=read.table("E:/Biase3celltypes.txt",header = T)  
Express=Express[,-1]  
#Preprocess the data(You can omit these steps)  
len=ncol(Express)  
Express=Express[apply(Express,1,function(x) sum(x>1)>len*0.05),]  
Express=Express[apply(Express,1,function(x) sum(x>1)<len*0.95),]  
Express=apply(Express,2,as.numeric)  
#Distribute the number of cores in a computer  
library(parallel)  
detectCores()  
cl <- makeCluster(5)  
#Load the file  
library(SNFtool)  
library(apcluster)  
alpha=0.5  
clusterExport(cl,"KNN_SMI",envir = environment())  
clusterExport(cl,"Express",envir = environment())  
clusterExport(cl,"alpha",envir = environment())  
#Parallel computation.You can change the parameters  
parLapply(cl, c(1,2,3,4,5),K=10,T=50,X1=100,X2=200,X3=300,X4=400,X5=500, select_features1)  
stopCluster(cl)  
#Read the resulting group file  
group1=read.table("./data1.txt",header = T,quote = "",sep=' ')  
group2=read.table("./data2.txt",header = T,quote = "",sep=' ')  
group3=read.table("./data3.txt",header = T,quote = "",sep=' ')  
group4=read.table("./data4.txt",header = T,quote = "",sep=' ')  
group5=read.table("./data5.txt",header = T,quote = "",sep=' ')  
group1=t(group1)  
group2=t(group2)  
group3=t(group3)  
group4=t(group4)  
group5=t(group5)  
max(group1)  
max(group2)  
max(group3)  
max(group4)  
max(group5)  
#Merge the results  
tt=result(group1,group2,group3,group4,group5)  
#Spectral clustering  
#Selection of digital  
#number is the median in max(group1),max(group2),max(group3),max(group4),max(group5)  
res=spectralClustering(tt, number)  
library(mclust)  
#Compute ARI  
a<-c(x,x,.........,y,y,y,y)  
a<-c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)  
b=as.vector(res)  
adjustedRandIndex(a,b)  






