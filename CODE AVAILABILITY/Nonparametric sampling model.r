################################
#rarefaction on discovery rate
################################



################################
#XT is the vector of species sequence abundance at time T
#a size of total number of species being sequenced up to time T
#Xt1 is species sequence abundance vector at time t1
#Xt2 is species sequence abundance vector at time t2
#in all these vectors, species identity and order are not important
Dt2<-function(Xt1,Xt2,XT)
{
NT=sum(XT)
ST=length(XT[XT>0])
one=NT-XT
Nt1=sum(Xt1)
Nt2=sum(Xt2)
###
two=0:(Nt1-1)
three=0:(Nt2-1)
v=0
for(i in 1:ST)
{
p1=one[i]-two
p2=NT-two
p=prod(p1/p2)
q1=one[i]-three
q2=NT-three
q=prod(q1/q2)
v=v+p-q
}#i
return(v)
}#
################################


################################
#Xt1 is the vector of species sequence abundance at time t1
#Xt2 is the vector of species sequence abundance at time t2
#XT is the vector of species sequence abundance at time T
#in all these vectors, species identity and order are not important
Dt3<-function(Xt1,Xt2,XT)
{
ST=length(XT[XT>0])
NT=sum(XT)
St1=length(Xt1[Xt1>0])
Nt1=sum(Xt1)
Nt2=sum(Xt2)
########
one=ST-St1
p=1/one
dN=Nt2-Nt1
f1=length(which(Xt1==1))
C=1-f1/Nt1+1e-300
v=0
for(i in 0:dN)
{
b=lgamma(dN+1)-lgamma(i+1)-lgamma(dN-i+1)
part=log(one)+log(1-(1-p)^i)+b+i*log(1-C)+(dN-i)*log(C)
v=v+exp(part)
}#i
return(v)
}
################################




#################################
#################################
#################################
#empirical test
dat=read.table("D:\\other study\\濒危物种研究\\N\\allCRanimals-N.csv",sep=",",header=TRUE)
years=as.numeric(as.character(dat[,1]))
mat=dat[,-1]
mat[is.na(mat)]=0
###############
Xt=list()
for(i in 1:length(years))
{
if(i==1)
{
this=mat[i,]
}else
{
this=colSums(mat[1:i,])
}
if(length(which(this==0))>0)
{
this=this[-which(this==0)]
}
Xt[[i]]=this
}#i
#################################

#################################
obs=v2=v3=vector()
for(i in 2:length(years))
{
v2=c(v2,Dt2(Xt[[i-1]],Xt[[i]],Xt[[30]]))
v3=c(v3,Dt3(Xt[[i-1]],Xt[[i]],Xt[[30]]))
obs=c(obs,length(Xt[[i]])-length(Xt[[i-1]]))
}#
#v2 is adjusted discovery rate based on Dt2 model
#v3 is adjusted discovery rate based on Dt3 model
#obs is the observed discovery rate, not sample size-adjusted
##############
length(which(obs>v2))
length(which(obs>v3))
###################