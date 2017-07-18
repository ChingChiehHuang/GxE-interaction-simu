library(iSKAT)
library(rareGE)
library(survey)

setwd('/home/cchuang/Simu/Cosi1/simulation')

ncase = 1000
ncontrol = 1000

significance.level = 1/100000

beta_GE1 = log(1.6) #1:3
beta_GE2 = log(1.3) #4:6
beta_G = 0
beta_E = 0
prevalence = 0.05
Pe = 0.5
n =  15000     # the number of replications
MAF <- 0.9

Haplo <- read.table('/home/cchuang/Simu/Cosi1/out.hap-1')[,3:740]
SNP <- read.table('/home/cchuang/Simu/Cosi1/out.pos-1',header = T)

SNP_pruned <- SNP[which(SNP[,7]< MAF),]

E <- sample(c(0,1),size=n,replace=TRUE,p=c(0.5,0.5))

for(k in 1:100){
  person <- matrix(0,n,nrow(SNP_pruned))
  control <- matrix(NA,n,nrow(SNP_pruned)+1)
  case <- matrix(NA,n,nrow(SNP_pruned)+1)
  Y  <- c()
  
  j <- 1
  repeat{
    for(i in 1:nrow(SNP_pruned)){
      if(sum(sample(Haplo[,which(SNP[,7]< MAF)][,i],size=2))==4){
        person[j,i] <- 0
      }else if(sum(sample(Haplo[,which(SNP[,7]< MAF)][,i],size=2))==3){
        person[j,i] <- 1
      }else{
        person[j,i] <- 2
      }
    }
    p <- exp(log(prevalence/(1-prevalence))+sum(beta_GE1*(person[j,1:3]*E[j]))+sum(beta_GE2*(person[j,4:6]*E[j])))/
      (1+exp(log(prevalence/(1-prevalence))+sum(beta_GE1*(person[j,1:3]*E[j]))+sum(beta_GE2*(person[j,4:6]*E[j]))))
    Y[j] <- sample(c(0,1),size=1,replace=TRUE,prob=c(1-p,p))
    
    if(Y[j]==0){
      control[j,] <- c(E[j],person[j,])
    }else{
      case[j,] <- c(E[j],person[j,])
    }
    
    if(nrow(matrix(case[-which(apply(case,1,function(x)all(is.na(x)))),],ncol=nrow(SNP_pruned)+1))>=ncase){
      break
    }
    j <- j+1
    
  }
  
  control <- control[-which(apply(control,1,function(x)all(is.na(x)))),]
  case <- case[-which(apply(case,1,function(x)all(is.na(x)))),]
  control <- control[1:ncase,]
  dim(control)
  dim(case)
  colnames(control)<- colnames(case)<- c('E',paste0('SNP',SNP_pruned[,1]))
  
  write.csv(case,file = paste0('case',k,'.csv'),col.names = TRUE ,
            row.names = FALSE,append=TRUE)
  write.csv(control,file = paste0('control',k,'.csv'),col.names = TRUE ,
            row.names = FALSE,append=TRUE)
  
}

setwd('/home/cchuang/Simu/Cosi1')

Age <- sample(18:65,size=2000,replace=TRUE)

SEBRIA_P <- c()
rareGE_RAN_P <- c()
rareGE_FIX_P <- c()
GESAT_P <- c()
iSKAT_P <- c()
MixGE_RAN_P <- c() 
MixGE_FIX_P <- c()

for(l in 1:100){
  
  case <- read.csv(paste('/home/cchuang/Simu/Cosi1/simulation/case',l,'.csv',sep=''),header=T)
  control <- read.csv(paste('/home/cchuang/Simu/Cosi1/simulation/control',l,'.csv',sep=''),header=T)
  y <- as.numeric(c(rep(1,1000),rep(0,1000)))
  E <- as.numeric(c(case[,1],control[,1]))
  gene_raw <- as.matrix(rbind(case[,-1],control[,-1]))
  
  #filering step
  out <- c()
  for(i in 1: ncol(gene_raw)){
    Z <- gene_raw[,i]*E
    if(nrow(summary(glm(y~Age+E+gene_raw[,i]+E*gene_raw[,i]))$coefficient)< 5){
      out[i] <- 0
    }else{ 
      out[i] <- summary(glm(y~Age+E+gene_raw[,i]+E*gene_raw[,i]))$coefficient[5,4]
    }
  }
  gene <- gene_raw[,out < 0.1]
  
  #rareGE_RAN_P
  rareGE_RAN_P[l] <- rareGE(phenotype= y, genotypes=gene, covariates=cbind(E,Age),family = "binomial",B = 1000)$pINT_RAN
  
  #rareGE_FIX_P
  rareGE_FIX_P[l] <- rareGE(phenotype= y, genotypes=gene, covariates=cbind(E,Age),family = "binomial",B = 1000)$pINT_FIX
  
  #SEBRIA_P
  Z_score <- apply(gene,2,function(x)1-summary(glm(E~x))$coefficients[2,4])
  indicator <- function(x) ifelse(abs(x)>=min(abs(Z_score)[rank(abs(Z_score))> 0.9*length(Z_score)]),1,0)
  W <- indicator(Z_score)*sign(Z_score)+0.0001
  coef <- summary(lm(y~E+Age+gene+E*as.matrix(gene)%*%W ))$coefficients
  SEBRIA_P[l]  <- coef[nrow(coef),4]
  
  
  #MixGE_RAN_P
  W <- cbind(rep(1/ncol(gene),ncol(gene)))
  model1 <- glm(y ~ Age + E + gene%*%W + diag(E)%*%gene%*%W )
  mu_hat <- predict(model1)
  res <- model1$residuals
  S_tau <- t(res)%*%diag(E)%*%gene%*%t(gene)%*%diag(E)%*%(res)
  
  S1 <- c()
  for(a in 1:length(y)){
    S1[a] <-  mu_hat[a]*(1-mu_hat[a])
  }
  D1 <- diag(S1)                
  sloD1 <- solve(D1)
  M <- cbind(Age,E,gene%*%W,diag(E)%*%gene%*%W)
  P <-  sloD1-sloD1%*%M%*%solve((t(M)%*%sloD1%*%M))%*%t(M)%*%sloD1
  eig2 <- eigen(P)
  evals2 <- eig2$values[eig2$values != 0]
  MixGE_RAN_P[l] <- pchisqsum(S_tau, df=rep(1, length(evals2)),a= evals2,lower.tail = F,method="sad")
  
  #MixGE_FIX_P
  
  model2 <- glm(y ~ Age + E + gene%*%W  )
  mu_hat2 <- predict(model2)
  U_pi <- t(diag(E)%*%gene%*%W)%*%(y-mu_hat2)
  
  S <- c()
  for(b in 1:length(y)){
    S[b] <-  mu_hat2[b]*(1-mu_hat2[b])
  }
  D <- diag(S)
  V <- cbind(Age,E,gene%*%W)
  sigma <- t(diag(E)%*%gene%*%W)%*%(D-D%*%V%*%solve((t(V)%*%D%*%V))%*%t(V)%*%D)%*%(diag(E)%*%gene%*%W)
  
  fix  <- (U_pi*U_pi)/sigma
  MixGE_FIX_P[l] <- pchisq(fix,df=1,lower.tail=FALSE) 
  
  #GESAT_P
  GESAT_P[l] <- GESAT(Z=gene, Y=as.matrix(y), E=as.matrix(E), X=as.matrix(Age), out_type="D")$pvalue 
  #iSKAT_P
  iSKAT_P[l] <- iSKAT(Z=gene, Y=as.matrix(y), E=as.matrix(E), X=as.matrix(Age), out_type="D"
                      , scale.Z=TRUE, MAF_cutoff=1, weights.beta=NULL)$pvalue
}


Output <- cbind(SEBRIA_P, rareGE_RAN_P, rareGE_FIX_P, GESAT_P, iSKAT_P, MixGE_RAN_P, MixGE_FIX_P)
colnames(Output) <- c('SEBRIA_P', 'rareGE_RAN_P', 'rareGE_FIX_P', 'GESAT_P', 'iSKAT_P', 'MixGE_RAN_P', 'MixGE_FIX_P')

write.csv(Output,file = 'Output.csv',col.names = TRUE ,row.names = FALSE,append=TRUE)
Output
