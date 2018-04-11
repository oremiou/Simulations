library(rjags)
##### fit the model ####
# Start the clock!
ptm <- proc.time()

model1.string <-  "
model{
for(i in 1:NS){
dm[i]<-d[t2[i]]-d[t1[i]]
prec[i]<-1/(SE[i]*SE[i])
y[i]~dnorm(dm[i],prec[i])}
d[1]<-0
for(i in 2:NT){
d[i]~dnorm(md,sd)}
md~dnorm(0,0.01)
sd<-1/(td*td)
td~dunif(0,5)

for (i in 1:NT){
for (j in i:NT){
D[j,i]<-d[j]-d[i]}}

#TreatmeNT hierarchy
  order[1:NT]<- NT+1- rank(d[1:NT])
for(k in 1:NT) {
# this is when the outcome is positive - omit  'NT+1-' when the outcome is negative
most.effective[k]<-equals(order[k],1)
for(j in 1:NT) {
effectiveness[k,j]<- equals(order[k],j)}}
for(k in 1:NT) {
for(j in 1:NT) {
cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])}}

#SUCRAS#
for(k in 1:NT) {
SUCRA[k]<- sum(cumeffectiveness[k,1:(NT-1)]) /(NT-1)}}
"
####### definitions 
jags.m=list()
samps=list()
A1=list()
count1=N.treat*(N.treat-1)/2+N.treat
lower=list()
upper=list()
lower_D=list()
upper_D=list()
count2=N.treat*(N.treat-1)/2
SUCRA=list()
SUCRA_max=c()
SUCRA_min=c()
best.worst=c()
best.2worst=c()
sum7=c()
sum8=c()
count3=N.treat*(N.treat-1)/2
lower3=list()
sortedSUCRA=list()
best=c()
worst=c()
sec.worst=c()
third.worst=c()
fourth.worst=c()
count4=N.treat*(N.treat-1)/2+N.treat+N.treat-1

for (i in 1:N.sim){
model1.spec<-textConnection(model1.string) 
data <- list(y = data1[[i]]$TE,SE=data1[[i]]$seTE, NS=length(data1[[i]]$studlab), t1=data1[[i]]$t1,t2=data1[[i]]$t2, NT=N.treat)
jags.m[[i]] <- jags.model(model1.spec, data = data, n.chains = 2, n.adapt = 5000)
print(i)

params <- c("d[2]","d[3]","d[4]","d[5]","d[6]","d[7]","d[8]","d[9]","d[10]","md","td")

for (i in 1:9){
  for (j in (i+1):10){
    params=c(params, paste("D[",j,",",i,"]",sep=""))
    closeAllConnections()
  }
}

for (i in 1:10){
  params=c(params, paste("SUCRA[",i,"]" ,sep=""))
}
}

for (i in 1:N.sim){
samps[[i]] <- coda.samples(jags.m[[i]], params, n.iter = 30000)

#plot(samps)
burn.in <- 15000
A1[[i]]=summary(window(samps[[i]],start = burn.in))

lower[[i]]= A1[[i]]$quantiles[(count1+1):(length(A1[[i]]$quantiles[,1])-2),1]
upper[[i]]= A1[[i]]$quantiles[(count1+1):(length(A1[[i]]$quantiles[,1])-2),5]
SUCRA[[i]]=A1[[i]]$statistics[(count2+1):(count2+N.treat),1]
lower_D[[i]]=A1[[i]]$quantiles[1:count2,1]
upper_D[[i]]=A1[[i]]$quantiles[1:count2,5]
best[i]=(which(SUCRA[[i]]==sort(SUCRA[[i]])[10]))-1+10*((which(SUCRA[[i]]==sort(SUCRA[[i]])[10]))==1)
worst[i]=(which(SUCRA[[i]]==sort(SUCRA[[i]])[2]))-1+10*((which(SUCRA[[i]]==sort(SUCRA[[i]])[2]))==1)
sec.worst[i]=(which(SUCRA[[i]]==sort(SUCRA[[i]])[3]))-1+10*((which(SUCRA[[i]]==sort(SUCRA[[i]])[3]))==1)
third.worst[i]=(which(SUCRA[[i]]==sort(SUCRA[[i]])[4]))-1+10*((which(SUCRA[[i]]==sort(SUCRA[[i]])[4]))==1)
fourth.worst[i]=(which(SUCRA[[i]]==sort(SUCRA[[i]])[5]))-1+10*((which(SUCRA[[i]]==sort(SUCRA[[i]])[5]))==1)
best.worst[i]=which(rownames(A1[[i]]$quantiles)==paste("D[",max(best[i],worst[i]),",",min(best[i],worst[i]),"]",sep=""))
sum7[i]=sum((A1[[i]]$quantiles[best.worst[i],1])>0)+sum(A1[[i]]$quantiles[best.worst[i],5]<0)
best.2worst[i]=which(rownames(A1[[i]]$statistics)==paste("D[",max(best[i],sec.worst[i]),",",min(best[i],sec.worst[i]),"]",sep=""))
sum8[i]=sum((A1[[i]]$quantiles[best.2worst[i],1])>0)+sum(A1[[i]]$quantiles[best.2worst[i],5]<0)
#lower3[[i]]= A1[[i]]$quantiles[1:count3,1]


samps[[i]]=c(0)
A1[[i]]=c(0)
jags.m[[i]]=c(0)
print(i)
closeAllConnections()
}


#plot(samps[[10]][[1]][,66])

# Stop the clock
proc.time() - ptm


#### when no treatment effects ######


 sum1=c(rep(100,N.sim))
 for (i in 1:N.sim)
 {
 sum1[i]=sum(lower[[i]]>0)+sum(upper[[i]]<0)
 }
 
 sum(sum1>0)/N.sim ## percent of NMAs with at least one st. sign. basic parameter. 


 for (i in 1:N.sim){
   SUCRA_max[i]=max(SUCRA[[i]])
   SUCRA_min[i]=min(SUCRA[[i]])
 }
 
 mean(SUCRA_max)
 min(SUCRA_max)
 max(SUCRA_max)
 
 mean(SUCRA_min)
 min(SUCRA_min)
 max(SUCRA_min)
 
 
 sum6=c(rep(100,N.sim))
 for (i in 1:N.sim)
 {
   sum6[i]=sum(lower_D[[i]]>0)+sum(upper_D[[i]]<0)
 }
 
 sum(sum6>0)/N.sim  # percent of NMAs with at least one stat. sign. finding

 
 
 
 
 ####### when non-zero treatment effects, SCENARIO 5 ###########
 
 ### false positive rates
 sum(sum7>0)/N.sim ### % of networks where the best vs worst treatment (not placebo) are SS
 sum(sum8>0)/N.sim ### % of networks where the best vs 2nd worst treatment (not placebo) are SS
 
 ### false positive rates in all treatment-treatment comparisons
 sum11=c(rep(100,N.sim))
 sum12=c(rep(100,N.sim))
 low=list()
 up=list()
 low1=list()
 up1=list()
 sum11=c(rep(100,N.sim))
 sum12=c(rep(100,N.sim))
 for (i in 1:N.sim)
 { low[[i]]=as.matrix(lower_D[[i]])
 up[[i]]=as.matrix(upper_D[[i]])
 low1[[i]]=low[[i]][substr(rownames(low[[i]]),nchar(rownames(low[[i]]))-1,nchar(rownames(low[[i]]))-1)!="1"]
 up1[[i]]=up[[i]][substr(rownames(up[[i]]),nchar(rownames(up[[i]]))-1,nchar(rownames(up[[i]]))-1)!="1"]
 }
 
 for (i in 1:N.sim)
 {  sum11[i]=sum(low1[[i]]>0)
 sum12[i]=sum(up1[[i]]<0)
 }
 
 sum(sum11+sum12)/(N.sim*(N.treat-1)*(N.treat-2)/2) # false positive rates in all comparisons
 
 
 ### power to detect effects vs placebo
 sum9=c(rep(100,N.sim))
 for (i in 1:N.sim)
 {  sum9[i]=sum(lower[[i]]>0)}
 sum(sum9)/(N.sim*(N.treat-1))  # percent of treat vs placebo comparisons that were SS in all datasets 
 sum(sum9!=0)/N.sim# percent of networks with at least 1 SS in comparisons vs placebo 
 
 

 
#### proxeiro #####

net1[[11]]$TE.fixed[11,]
net1[[11]]$lower.fixed[11,]
net1[[11]]$upper.fixed[11,]
 
 gelman.diag(samps)
 gelman.plot(samps)
 
 