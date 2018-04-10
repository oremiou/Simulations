set.seed(42)
rm(list = ls())
library(netmeta)


############################################################################################################
##### SCENARIO 1- STAR NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO########################################
##### STAR NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO #
N.sim=1000
N.treat=10
N.studies.comp=1 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c(rep(1,N.studies.comp*(N.treat-1)))
t2=c()
for (j in 2:N.treat)  {t2=c(t2,rep(j,N.studies.comp))}

###generate data
for (i in 1:N.sim)
{
 seTE=c(rep(0,N.studies.comp*(N.treat-1))) 
 TE=c(rep(0,N.studies.comp*(N.treat-1))) 
 data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=0.2*rchisq(N.studies.comp*(N.treat-1),df=1)+0.5
  data1[[i]]$TE=rnorm(N.studies.comp*(N.treat-1),0,data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.studies.comp*(N.treat-1)))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}


##### netmeta
for (i in 1:N.sim)
{
   net1[[i]]=netmeta(TE=TE,comb.fixed=TRUE, comb.random=F, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
}  
#netgraph(net1[[1]], plastic=F, multiarm = F, points=T, cex.points=5, col=1)  


#### find st. sing. findings, all treatment effects
sum1=c()
sum2=c()
sum3=c()
sum4=c()
sum5=c()

for (i in 1:N.sim)
{
sum1[i]=sum(net1[[i]]$lower.fixed[upper.tri(net1[[i]]$lower.fixed)]>0) ### number of comparisons with CIlower>0
sum2[i]=sum(net1[[i]]$upper.fixed[upper.tri(net1[[i]]$upper.fixed)]<0) ### number of comparisons with CIupper<0
sum3[i]=sum1[i]+sum2[i] ### number of comparisons statistically significant
}
sum(sum3>0)/N.sim # percent of NMAs with at least one stat. sign. finding


for (i in 1:N.sim)
{
sum4[i]=sum(net1[[i]]$lower.fixed[1,2:N.treat]>0) ### number of NMAs with at least one statistically basic parameter (lower)
sum5[i]=sum(net1[[i]]$upper.fixed[1,2:N.treat]<0) ### number of NMAs with at least one statistically basic parameter (lower)
}
sum((sum4+sum5)>0)/N.sim ## percent of NMAs with at least one st. sign. basic parameter. 


SUCRA_max=c()
SUCRA_min=c()
treat.max=c()
for (i in 1:N.sim)
{
  SUCRA_max=  c(SUCRA_max, max(netrank(net1[[i]])[[1]]))
  SUCRA_min=  c(SUCRA_min, min(netrank(net1[[i]])[[1]]))
  treat.max=c(treat.max, c(which(netrank(net1[[i]])[[1]]==max(netrank(net1[[i]])[[1]]))))
}

mean(SUCRA_min)
min(SUCRA_min)
max(SUCRA_min)

mean(SUCRA_max)
min(SUCRA_max)
max(SUCRA_max)

###############


############################################################################################################
##### SCENARIO 2 - STAR NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO ,HETEROGENEITY, 3 ARMS PER COMPARISON #########################################################################################
##### STAR NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO ,HETEROGENEITY, 3 ARMS PER COMPARISON#
N.sim=1000
N.treat=10
N.studies.comp=3 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c(rep(1,N.studies.comp*(N.treat-1)))
t2=c()
for (j in 2:N.treat)  {t2=c(t2,rep(j,N.studies.comp))}
###generate data
tau_sq=rlnorm(N.sim,-3,1.5)
tau=sqrt(tau_sq)
for (i in 1:N.sim)
{
  seTE=c(rep(0,N.studies.comp*(N.treat-1))) 
  TE_true=c(rep(0,N.studies.comp*(N.treat-1))) 
  TE=c(rep(0,N.studies.comp*(N.treat-1))) 
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=0.2*rchisq(N.studies.comp*(N.treat-1),df=1)+0.5
  data1[[i]]$TE_true=rnorm(N.studies.comp*(N.treat-1),0,tau[i])
  data1[[i]]$TE=rnorm(N.studies.comp*(N.treat-1),data1[[i]]$TE_true,data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.studies.comp*(N.treat-1)))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}


##### netmeta
for (i in 1:N.sim)
{
  net1[[i]]=netmeta(TE=TE,comb.fixed=TRUE, comb.random=F, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
}  
#netgraph(net1[[1]], plastic=F, multiarm = F, points=T, cex.points=5, col=1)  


#### find st. sing. findings, all treatment effects
sum1=c()
sum2=c()
sum3=c()
sum4=c()
sum5=c()

for (i in 1:N.sim)
{
  sum1[i]=sum(net1[[i]]$lower.fixed[upper.tri(net1[[i]]$lower.fixed)]>0) ### number of comparisons with CIlower>0
  sum2[i]=sum(net1[[i]]$upper.fixed[upper.tri(net1[[i]]$upper.fixed)]<0) ### number of comparisons with CIupper<0
  sum3[i]=sum1[i]+sum2[i] ### number of comparisons statistically significant
}
sum(sum3>0)/N.sim # percent of NMAs with at least one stat. sign. finding


for (i in 1:N.sim)
{
  sum4[i]=sum(net1[[i]]$lower.fixed[1,2:N.treat]>0) ### number of NMAs with at least one statistically basic parameter (lower)
  sum5[i]=sum(net1[[i]]$upper.fixed[1,2:N.treat]<0) ### number of NMAs with at least one statistically basic parameter (lower)
}
sum((sum4+sum5)>0)/N.sim ## percent of NMAs with at least one st. sign. basic parameter. 


SUCRA_max=c()
SUCRA_min=c()
treat.max=c()
for (i in 1:N.sim)
{
  SUCRA_max=  c(SUCRA_max, max(netrank(net1[[i]])[[1]]))
  SUCRA_min=  c(SUCRA_min, min(netrank(net1[[i]])[[1]]))
  treat.max=c(treat.max, c(which(netrank(net1[[i]])[[1]]==max(netrank(net1[[i]])[[1]]))))
}

mean(SUCRA_min)
min(SUCRA_min)
max(SUCRA_min)

mean(SUCRA_max)
min(SUCRA_max)
max(SUCRA_max)

###############

##### SCENARIO 3 -FULLY CONNECTED NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO - HOMOGENEITY#######################
##### FULLY CONNECTED NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO #

N.sim=1000
N.treat=10
N.studies.comp=1 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c()
t2=c()
for (i in 1:(N.treat-1)){
  for (k in (i+1):N.treat){
    for(j in 1:N.studies.comp){
    t1=c(t1,i)
    t2=c(t2,k)  
    }
   }
}

N.stud=length(t1)
for (i in 1:N.sim)
{
  seTE=c(rep(0,N.stud)) 
  TE=c(rep(0,N.stud)) 
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=0.2*rchisq(N.studies.comp*(N.treat-1),df=1)+0.5
  data1[[i]]$TE=rnorm(N.stud,0,data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}

for (i in 1:N.sim)
{
  net1[[i]]=netmeta(TE=TE, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
}  
#netgraph(net1[[1]], plastic=F, multiarm = F, points=T, cex.points=5, col=1)  




#### find st. sing. findings, all treatment effects
sum1=c()
sum2=c()
sum3=c()
sum4=c()
sum5=c()

for (i in 1:N.sim)
{
  sum1[i]=sum(net1[[i]]$lower.fixed[upper.tri(net1[[i]]$lower.fixed)]>0) ### number of comparisons with CIlower>0
  sum2[i]=sum(net1[[i]]$upper.fixed[upper.tri(net1[[i]]$upper.fixed)]<0) ### number of comparisons with CIupper<0
  sum3[i]=sum1[i]+sum2[i] ### number of comparisons statistically significant
}
sum(sum3>0)/N.sim # percent of NMAs with at least one stat. sign. finding


for (i in 1:N.sim)
{
  sum4[i]=sum(net1[[i]]$lower.fixed[1,2:N.treat]>0) ### number of NMAs with at least one statistically basic parameter (lower)
  sum5[i]=sum(net1[[i]]$upper.fixed[1,2:N.treat]<0) ### number of NMAs with at least one statistically basic parameter (lower)
}
sum((sum4+sum5)>0)/N.sim ## percent of NMAs with at least one st. sign. basic parameter. 


SUCRA_max=c()
SUCRA_min=c()
treat.max=c()
for (i in 1:N.sim)
{
  SUCRA_max=  c(SUCRA_max, max(netrank(net1[[i]])[[1]]))
  SUCRA_min=  c(SUCRA_min, min(netrank(net1[[i]])[[1]]))
  treat.max=c(treat.max, c(which(netrank(net1[[i]])[[1]]==max(netrank(net1[[i]])[[1]]))))
}

mean(SUCRA_min)
min(SUCRA_min)
max(SUCRA_min)

mean(SUCRA_max)
min(SUCRA_max)
max(SUCRA_max)


############################################################################################################

##### SCENARIO 4 - FULLY CONNECTED NETWORK - ZERO EFFECTS FOR ALL TREATMENTS VS PLACEBO - HETEROGENEITY  #########################################################################################
N.sim=1000
N.treat=10
N.studies.comp=3 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c()
t2=c()
for (i in 1:(N.treat-1)){
  for (k in (i+1):N.treat){
    for(j in 1:N.studies.comp){
      t1=c(t1,i)
      t2=c(t2,k)  
    }
  }
}
###generate data
tau_sq=rlnorm(N.sim,-3,1.5)
tau=sqrt(tau_sq)

N.stud=length(t1)
for (i in 1:N.sim)
{
  seTE=c(rep(0,N.stud)) 
  TE_true=c(rep(0,N.stud)) 
  TE=c(rep(0,N.stud))
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=0.2*rchisq(N.stud,df=1)+0.5
  data1[[i]]$TE_true=rnorm(N.stud,0,tau[i])
  data1[[i]]$TE=rnorm(N.stud,data1[[i]]$TE_true,data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}






#####################


##### SCENARIO 5 - STAR NETWORK - EQUAL EFFECTS FOR ALL TREATMENTS VS PLACEBO ######
N.sim=1000
N.treat=10
N.studies.comp=1 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c(rep(1,N.studies.comp*(N.treat-1)))
t2=c()
for (j in 2:N.treat)  {t2=c(t2,rep(j,N.studies.comp))}
N.stud=length(t1)
# all treatment effects vs. placebo are equal to 1. No relative effects between treatments

for (i in 1:N.sim)
{
  seTE=c(rep(0,N.studies.comp*(N.treat-1))) 
  TE=c(rep(0,N.studies.comp*(N.treat-1))) 
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=0.2*rchisq(N.stud,df=1)+0.5
  data1[[i]]$TE=rnorm(N.studies.comp*(N.treat-1),2,data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.studies.comp*(N.treat-1)))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}



for (i in 1:N.sim)
{
  net1[[i]]=netmeta(TE=TE,comb.fixed=T, comb.random=F, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
} 

netleague(net1[[1]])
#netgraph(net1[[1]], plastic=F, multiarm = F, points=T, cex.points=5, col=1)  





#############
####### proxeiro
sum.1v9=c()
sum.1v8=c()
sum.1v7=c()
sum.1v6=c()
sum.1v5=c()
sum.1v4=c()

best.treat=c()
fourth.best.treat=c()
fifth.best.treat=c()
sixth.best.treat=c()
seventh.best.treat=c()
eighth.best.treat=c()


worst.treat=c()
for (i in 1:N.sim){
  best.treat[i]=order(netrank(net1[[i]])[[1]])[10]
  worst.treat[i]=order(netrank(net1[[i]])[[1]])[2]
  eighth.best.treat[i]=order(netrank(net1[[i]])[[1]])[3]
  seventh.best.treat[i]=order(netrank(net1[[i]])[[1]])[4]
  sixth.best.treat[i]=order(netrank(net1[[i]])[[1]])[5] 
    fifth.best.treat[i]=order(netrank(net1[[i]])[[1]])[6]
    fourth.best.treat[i]=order(netrank(net1[[i]])[[1]])[7]
  
  sum.1v9=c(sum.1v9,net1[[i]]$upper.fixed[best.treat[i],worst.treat[i]]<0)
  sum.1v8=c(sum.1v8,net1[[i]]$upper.fixed[best.treat[i],eighth.best.treat[i]]<0)
  sum.1v7=c(sum.1v7,net1[[i]]$upper.fixed[best.treat[i],seventh.best.treat[i]]<0)
  sum.1v6=c(sum.1v6,net1[[i]]$upper.fixed[best.treat[i],sixth.best.treat[i]]<0)
  sum.1v5=c(sum.1v5,net1[[i]]$upper.fixed[best.treat[i],fifth.best.treat[i]]<0) 
  sum.1v4=c(sum.1v4,net1[[i]]$upper.fixed[best.treat[i],fourth.best.treat[i]]<0)
}

  sum(sum.1v9)/N.sim
  sum(sum.1v8)/N.sim
  sum(sum.1v7)/N.sim
  sum(sum.1v6)/N.sim
  sum(sum.1v5)/N.sim
  sum(sum.1v4)/N.sim
  
  # power
  count.power=c()
  for (i in 1:N.sim){
    count.power[i]= sum(net1[[i]]$lower.fixed[1,2:10]>0)
   # sum(net1[[i]]$upper.fixed[1,2:10]<0)
  }
  sum(count.power)/(N.sim*N.treat)

  
  
##################################################################################################################################

  
  
  
  
  
  
##### SCENARIO 4  ##########################################################################################
##### FULLY CONNECTED NETWORK - EQUAL EFFECTS FOR ALL TREATMENTS VS PLACEBO #

N.sim=100
N.treat=10
N.studies.comp=1 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c()
t2=c()
for (i in 1:(N.treat-1)){
  for (k in (i+1):N.treat){
    for(j in 1:N.studies.comp){
      t1=c(t1,i)
      t2=c(t2,k)  
    }
  }
}

diff=as.numeric(t1!=1)+1
mu=c(2,0)


###generate data

N.stud=length(t1)
for (i in 1:N.sim)
{
  seTE=c(rep(0,N.stud)) 
  TE=c(rep(0,N.stud)) 
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=runif(N.stud,0.9,1.1)
  data1[[i]]$TE=rnorm(N.stud,mu[diff],data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}


for (i in 1:N.sim)
{
  net1[[i]]=netmeta(TE=TE,comb.fixed=TRUE, comb.random=F, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
} 
#netgraph(net1[[1]], plastic=F, multiarm = F, points=T, cex.points=5, col=1)  


sum.1v9=c()
sum.1v8=c()
sum.1v7=c()
sum.1v6=c()
sum.1v5=c()
sum.1v4=c()

best.treat=c()
fourth.best.treat=c()
fifth.best.treat=c()
sixth.best.treat=c()
seventh.best.treat=c()
eighth.best.treat=c()


worst.treat=c()
for (i in 1:N.sim){
  best.treat[i]=order(netrank(net1[[i]])[[1]])[10]
  worst.treat[i]=order(netrank(net1[[i]])[[1]])[2]
  eighth.best.treat[i]=order(netrank(net1[[i]])[[1]])[3]
  seventh.best.treat[i]=order(netrank(net1[[i]])[[1]])[4]
  sixth.best.treat[i]=order(netrank(net1[[i]])[[1]])[5] 
  fifth.best.treat[i]=order(netrank(net1[[i]])[[1]])[6]
  fourth.best.treat[i]=order(netrank(net1[[i]])[[1]])[7]
  
  sum.1v9=c(sum.1v9,net1[[i]]$upper.fixed[best.treat[i],worst.treat[i]]<0)
  sum.1v8=c(sum.1v8,net1[[i]]$upper.fixed[best.treat[i],eighth.best.treat[i]]<0)
  sum.1v7=c(sum.1v7,net1[[i]]$upper.fixed[best.treat[i],seventh.best.treat[i]]<0)
  sum.1v6=c(sum.1v6,net1[[i]]$upper.fixed[best.treat[i],sixth.best.treat[i]]<0)
  sum.1v5=c(sum.1v5,net1[[i]]$upper.fixed[best.treat[i],fifth.best.treat[i]]<0) 
  sum.1v4=c(sum.1v4,net1[[i]]$upper.fixed[best.treat[i],fourth.best.treat[i]]<0)
}

sum(sum.1v9)/N.sim
sum(sum.1v8)/N.sim
sum(sum.1v7)/N.sim
sum(sum.1v6)/N.sim
sum(sum.1v5)/N.sim
sum(sum.1v4)/N.sim

# power
count.power=c()
for (i in 1:N.sim){
  count.power[i]= sum(net1[[i]]$lower.fixed[1,2:10]>0)
  # sum(net1[[i]]$upper.fixed[1,2:10]<0)
}
sum(count.power)/(N.sim*N.treat)
############################################################################################################



##### SCENARIO 5  ##########################################################################################
##### FULLY CONNECTED NETWORK - DIFFERENT EFFECTS FOR ALL TREATMENTS VS PLACEBO #

N.sim=100
N.treat=10
N.studies.comp=1 # studies per comparison
data1=list()
net1=list()

#define the treatments in the studies
t1=c()
t2=c()
for (i in 1:(N.treat-1)){
  for (k in (i+1):N.treat){
    for(j in 1:N.studies.comp){
      t1=c(t1,i)
      t2=c(t2,k)  
    }
  }
}

mu=0:(N.treat-1)/N.treat*2


###generate data

N.stud=length(t1)
for (i in 1:N.sim)
{
  seTE=c(rep(0,N.stud)) 
  TE=c(rep(0,N.stud)) 
  data1[[i]]=data.frame(TE,seTE)
  data1[[i]]$seTE=runif(N.stud,0.9,1.1)
  data1[[i]]$TE=rnorm(N.stud,(mu[t2]-mu[t1]),data1[[i]]$seTE)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$t1=t1
  data1[[i]]$t2=t2
}


for (i in 1:N.sim)
{
  net1[[i]]=netmeta(TE=TE,comb.fixed=TRUE, comb.random=F, seTE=seTE,treat1=t1, treat2=t2, studlab=studlab,data=data1[[i]]) 
} 

# power
count.power=c()
for (i in 1:N.sim){
  count.power[i]= sum(
    net1[[i]]$lower.fixed[upper.tri(net1[[i]]$lower.fixed)]
    >0)
  # sum(net1[[i]]$upper.fixed[1,2:10]<0)
}
sum(count.power)/(N.sim*N.treat*(N.treat-1)/2)
############################################################################################################








####### proxeiro ######
N.sim=100
se=runif(N.sim,0.5,1.5)
y=rnorm(N.sim,0,se)
CIU=y+1.96*se
CID=y-1.96*se
(sum(CIU<0)+sum(CID>0))/N.sim

### sanity check
CIU=data1[[1]]$TE+1.96*data1[[1]]$seTE
CID=data1[[1]]$TE-1.96*data1[[1]]$seTE
(sum(CIU<0)+sum(CID>0))/length(data1[[1]]$TE)

############################
