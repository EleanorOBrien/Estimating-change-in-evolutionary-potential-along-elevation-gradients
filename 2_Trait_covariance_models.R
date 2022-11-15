#Code to estimate additive genetic variance, heritability and maternal effect variance in each trait, and genetic and maternal covariances among pairs of traits in flies produced in laboraotory crosses


#load packages

library(MCMCglmm)


#load data

#Pedigree files
ped_Pal <- read.csv("Paluma_pedigree.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE) 
ped_Dan <- read.csv("Danbulla_pedigree.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE) 

#Phenotype files
dat_Pal <- read.csv("Paluma_Traits.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE) 
dat_Dan <- read.csv("Danbulla_Traits.csv", header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE) 


#Ensure variable types are specified correctly

dat_Pal$ID <- as.factor(dat_Pal$ID)
dat_Pal$Mother <- as.factor(dat_Pal$Mother)
dat_Pal$Batch_CT <- as.factor(dat_Pal$Batch_CT)
dat_Pal$Batch_HT <- as.factor(dat_Pal$Batch_HT)
dat_Pal$Site_number <- as.factor(dat_Pal$Site_number)
dat_Pal$Cross_WS <- as.factor(dat_Pal$Cross_WS)
dat_Pal$Replicate_CT_HT <- as.factor(dat_Pal$Replicate)


dat_Dan$ID <- as.factor(dat_Dan$ID)
dat_Dan$Mother <- as.factor(dat_Dan$Mother)
dat_Dan$Batch_CT <- as.factor(dat_Dan$Batch_CT)
dat_Dan$Batch_HT <- as.factor(dat_Dan$Batch_HT)
dat_Dan$Site_number <- as.factor(dat_Dan$Site_number)
dat_Dan$Cross_WS <- as.factor(dat_Dan$Cross_WS)
dat_Dan$Replicate_CT_HT <- as.factor(dat_Dan$Replicate)

ped_Pal$ID <- as.factor(ped_Pal$ID)
ped_Pal$Dam <- as.factor(ped_Pal$Dam)
ped_Pal$Sire <- as.factor(ped_Pal$Sire)

ped_Dan$ID <- as.factor(ped_Dan$ID) 
ped_Dan$Dam <- as.factor(ped_Dan$Dam)
ped_Dan$Sire <- as.factor(ped_Dan$Sire)


#Convert wing size to mm
dat_Pal$WS_mm <- (dat_Pal$Centroid_Size)*0.000861
dat_Dan$WS_mm <- (dat_Dan$Centroid_Size)*0.000861


#Calculate raw means, variances and standard deviation of each trait

Raw_ave_PalCT <- mean(dat_Pal$Offspring_CT, na.rm = TRUE)
Raw_Vp_PalCT <- var(dat_Pal$Offspring_CT, na.rm = TRUE)
Raw_SD_PalCT <- sqrt(Raw_Vp_PalCT)

Raw_ave_PalHT <- mean(dat_Pal$TimeHT_MINS, na.rm = TRUE)
Raw_Vp_PalHT <- var(dat_Pal$TimeHT_MINS, na.rm = TRUE)
Raw_SD_PalHT <- sqrt(Raw_Vp_PalHT)

Raw_ave_PalWS <- mean(dat_Pal$WS_mm, na.rm = TRUE)
Raw_Vp_PalWS <- var(dat_Pal$WS_mm, na.rm = TRUE)
Raw_SD_PalWS <- sqrt(Raw_Vp_PalWS)


Raw_ave_DanCT <- mean(dat_Dan$Offspring_CT, na.rm = TRUE)
Raw_Vp_DanCT <- var(dat_Dan$Offspring_CT, na.rm = TRUE)
Raw_SD_DanCT <- sqrt(Raw_Vp_DanCT)

Raw_ave_DanHT <- mean(dat_Dan$TimeHT_MINS, na.rm = TRUE)
Raw_Vp_DanHT <- var(dat_Dan$TimeHT_MINS, na.rm = TRUE)
Raw_SD_DanHT <- sqrt(Raw_Vp_DanHT)

Raw_ave_DanWS <- mean(dat_Dan$WS_mm, na.rm = TRUE)
Raw_Vp_DanWS <- var(dat_Dan$WS_mm, na.rm = TRUE)
Raw_SD_DanWS <- sqrt(Raw_Vp_DanWS)

#Standardise traits to mean of 0 and variance of 1

dat_Pal$Z_CT <- ((dat_Pal$Offspring_CT) - Raw_ave_PalCT)/Raw_SD_PalCT
dat_Pal$Z_HT <- ((dat_Pal$TimeHT_MINS) - Raw_ave_PalHT)/Raw_SD_PalHT
dat_Pal$Z_WS <- ((dat_Pal$WS_mm) - Raw_ave_PalWS)/Raw_SD_PalWS

dat_Dan$Z_CT <- ((dat_Dan$Offspring_CT) - Raw_ave_DanCT)/Raw_SD_DanCT
dat_Dan$Z_HT <- ((dat_Dan$TimeHT_MINS) - Raw_ave_DanHT)/Raw_SD_DanHT
dat_Dan$Z_WS <- ((dat_Dan$WS_mm) - Raw_ave_DanWS)/Raw_SD_DanWS




#Function to generate relationship matrices


ISO.A <- function(pedigree, f, r.d, r.d.d){
  
  Ped.num <- matrix(0, dim(pedigree)[1], 3)
  colnames(Ped.num) <- colnames(pedigree)
  
  Ped.num[, 1] <- match(pedigree[, 1], pedigree[,1], nomatch = 0)
  Ped.num[, 2] <- match(pedigree[, 2], pedigree[,1], nomatch = 0)
  Ped.num[, 3] <- match(pedigree[, 3], pedigree[,1], nomatch = 0)
  
  Ped.num <- as.data.frame(Ped.num)
  
  nid <- dim(Ped.num)[1]
  ISO.f <- which(Ped.num[, 2] == 0 & Ped.num[, 3] == 0)
  d.parents <- which(Ped.num[, 2] != 0 & Ped.num[, 3] == 0)
  d.off <- which(Ped.num[, 2] != 0 & Ped.num[, 3] != 0)
  
  A <- diag(dim(Ped.num)[1])
  rownames(A) <- pedigree[, 1]
  colnames(A) <- pedigree[, 1]
  
  for (i in ISO.f){
    sibs <- d.parents[which(Ped.num[d.parents, 2] == i)]
    A[sibs, i] <- r.d 
    A[i, sibs] <- A[sibs, i]
    A[sibs,sibs]  <- r.d.d
  }
  
  diag(A)[d.parents] <- 1 + f
  
  for (i in 1:dim(Ped.num)[1]){
    for (j in d.off){
      parents <- as.numeric(Ped.num[j,2:3])
      A[i,j] <- (A[parents[1], i] + A[parents[2],i])/2
      A[j,i] <-  A[i,j]
      A[j,j] <- 1 + A[parents[1],parents[2]]/2
    }
  }
  list(A = A)
}


#Multivariate models to estimate genetic and maternal covariance among traits within each gradient

#Limit analysis of trait covariances to just those with an ID (i.e. exclude field fitness)

dat_Pal_Gcov <- subset(dat_Pal, Trait!="Field_fitness")
dat_Dan_Gcov <- subset(dat_Dan, Trait!="Field_fitness")


#Paluma

#phenotypic variance of standardised traits
Vp_Pal_Z_CT <- var(dat_Pal_Gcov$Z_CT, na.rm = TRUE)
Vp_Pal_Z_HT <- var(dat_Pal_Gcov$Z_HT, na.rm = TRUE)
Vp_Pal_Z_WS <- var(dat_Pal_Gcov$Z_WS, na.rm = TRUE)

Vp_Pal <- diag(c(Vp_Pal_Z_CT, Vp_Pal_Z_HT, Vp_Pal_Z_WS))

A_Pal <- ISO.A(pedigree = ped_Pal, f = 0.256, r.d = 0, r.d.d = 0.511)$A #if this doesn't work, update packages


Ainv_Pal <- as(round(solve(A_Pal), 8),"CsparseMatrix")

#Prior specification

n.t_Pal <- dim(Vp_Pal)[1] #number of traits



prior_Pal <- list(R = list(V = Vp_Pal/3, n = n.t_Pal),
                  G = list(G1 = list(V = Vp_Pal/3, nu = n.t_Pal),
                           G2 = list(V = Vp_Pal/3, nu = n.t_Pal)))



#MCMC multivariate model

date()
mv.Pal <- MCMCglmm(cbind(Z_CT, Z_HT, Z_WS) ~ trait -1,
                   random = ~ us(trait):ID + us(trait):Mother,
                   rcov = ~ idh(trait):units,
                   #setting rcov to idh(trait):units instead of us(trait):units sets residual covariances to 0, which is appropriate given traits were measured on different individuals.  
                   family = c("gaussian", "gaussian", "gaussian"),
                   ginverse = list(ID = Ainv_Pal),
                   nitt=2500000, thin = 2000, burnin = 500000,
                   pr = T, #saves the BLUPs in m.1$Sol
                   prior = prior_Pal,
                   data = dat_Pal_Gcov)
date()


# save model object
saveRDS(mv.Pal, "mv.Pal")

#model diagnostics

plot(mv.Pal$VCV)
summary(mv.Pal)
autocorr(mv.Pal$VCV)

#G Matrix
G_Pal <- matrix(colMeans(mv.Pal$VCV[,1:9]), n.t_Pal, n.t_Pal)

# M (maternal effects) matrix
M_Pal <- matrix(colMeans(mv.Pal$VCV[,10:18]), n.t_Pal, n.t_Pal)


#mean heritability

#first set up an index to extract the variances 

index <- diag(matrix(1:9, n.t_Pal, n.t_Pal))
index_MPal <- diag(matrix(10:18, n.t_Pal, n.t_Pal))
index_RPal <- diag(matrix(19:21, n.t_Pal, n.t_Pal))

#calculate the heritabilities
h_sq_Pal <- mv.Pal$VCV[,index]/(mv.Pal$VCV[,index] + mv.Pal$VCV[,index_MPal] + mv.Pal$VCV[,index_RPal])

plot(h_sq_Pal)

colMeans(h_sq_Pal)
HPDinterval(h_sq_Pal)

m_varPal <- mv.Pal$VCV[,index_MPal]/(mv.Pal$VCV[,index] + mv.Pal$VCV[,index_MPal] + mv.Pal$VCV[,index_RPal])

colMeans(m_varPal)
HPDinterval(m_varPal)

#genetic correlations 

post.dist.g.cor_Pal <- posterior.cor(mv.Pal$VCV[,1:9])


#maternal effect correlations

post.dist.m.cor_Pal <- posterior.cor(mv.Pal$VCV[,10:18])

#for estimates of genetic and maternal effect correlations and HPDintervals

G.cor_Pal <- matrix(colMeans(post.dist.g.cor_Pal),n.t_Pal, n.t_Pal)

mean_plus_HPD_Pal <- list(lower = matrix(HPDinterval(post.dist.g.cor_Pal)[,1], n.t_Pal, n.t_Pal), mean = G.cor_Pal, upper = matrix(HPDinterval(post.dist.g.cor_Pal)[,2], n.t_Pal, n.t_Pal))

M.cor_Pal <- matrix(colMeans(post.dist.m.cor_Pal),n.t_Pal, n.t_Pal)

mean_plus_HPD_M_Pal <- list(lower = matrix(HPDinterval(post.dist.m.cor_Pal)[,1], n.t_Pal, n.t_Pal), mean = M.cor_Pal, upper = matrix(HPDinterval(post.dist.m.cor_Pal)[,2], n.t_Pal, n.t_Pal))




#Danbulla
Vp_Dan_Z_CT <- var(dat_Dan_Gcov$Z_CT, na.rm = TRUE)
Vp_Dan_Z_HT <- var(dat_Dan_Gcov$Z_HT, na.rm = TRUE)
Vp_Dan_Z_WS <- var(dat_Dan_Gcov$Z_WS, na.rm = TRUE)

Vp_Dan <- diag(c(Vp_Dan_Z_CT, Vp_Dan_Z_HT, Vp_Dan_Z_WS))

A_Dan <- ISO.A(pedigree = ped_Dan, f = 0.256, r.d = 0, r.d.d = 0.511)$A

any(eigen(A_Dan)$values <= 0) 

Ainv_Dan <- as(round(solve(A_Dan), 8),"CsparseMatrix")


#Prior specification
#Danbulla
n.t_Dan <- dim(Vp_Dan)[1] #number of traits

#The only difference then for a multivariate model is the degrees of freedom for the IW distribtion specified by the nu below. Typically the starting point is nu = number of traits but depending on the model you can go as far as nu = 0.002 + (number of traits - 1) 

prior_Dan <- list(R = list(V = Vp_Dan/3, nu = n.t_Dan),
                  G = list(G1 = list(V = Vp_Dan/3, nu = n.t_Dan),
                           G2 = list(V = Vp_Dan/3, nu = n.t_Dan)))

#Danbulla
START <- date()
mv.Dan <- MCMCglmm(cbind(Z_CT, Z_HT, Z_WS) ~ trait -1,
                   random = ~ us(trait):ID + us(trait):Mother,
                   rcov = ~ idh(trait):units,
                   #setting rcov to idh(trait):units instead of us(trait):units sets residual covariances to 0, which is appropriate given traits were measured on different individuals.  
                   family = c("gaussian", "gaussian", "gaussian"),
                   ginverse = list(ID = Ainv_Dan),
                   nitt=2500000, thin = 2000, burnin = 500000,
                   pr = T, #saves the BLUPs in m.1$Sol
                   prior = prior_Dan,
                   data = dat_Dan_Gcov) 

END <- date()

saveRDS(mv.Dan, "mv.Dan")

#model diagnostics

plot(mv.Dan$VCV)
summary(mv.Dan)
autocorr(mv.Dan$VCV)

#G Matrix
G_Dan <- matrix(colMeans(mv.Dan$VCV[,1:9]), n.t_Dan, n.t_Dan)

# M (maternal effects) matrix
M_Dan <- matrix(colMeans(mv.Dan$VCV[,10:18]), n.t_Dan, n.t_Dan)


#mean heritability

#first set up an index to extract the variances 

index <- diag(matrix(1:9, n.t_Dan, n.t_Dan))
index_MDan <- diag(matrix(10:18, n.t_Dan, n.t_Dan))
index_RDan <- diag(matrix(19:21, n.t_Dan, n.t_Dan))

#calculate the heritabilities
h_sq_Dan <- mv.Dan$VCV[,index]/(mv.Dan$VCV[,index] + mv.Dan$VCV[,index_MDan] + mv.Dan$VCV[,index_RDan])

plot(h_sq_Dan)

colMeans(h_sq_Dan)
HPDinterval(h_sq_Dan)

m_varDan <- mv.Dan$VCV[,index_MDan]/(mv.Dan$VCV[,index] + mv.Dan$VCV[,index_MDan] + mv.Dan$VCV[,index_RDan])

colMeans(m_varDan)
HPDinterval(m_varDan)

#genetic correlations 

post.dist.g.cor_Dan <- posterior.cor(mv.Dan$VCV[,1:9])


#maternal effect correlations

post.dist.m.cor_Dan <- posterior.cor(mv.Dan$VCV[,10:18])

#for estimates of genetic and maternal effect correlations and HPDintervals

G.cor_Dan <- matrix(colMeans(post.dist.g.cor_Dan),n.t_Dan, n.t_Dan)

mean_plus_HPD_Dan <- list(lower = matrix(HPDinterval(post.dist.g.cor_Dan)[,1], n.t_Dan, n.t_Dan), mean = G.cor_Dan, upper = matrix(HPDinterval(post.dist.g.cor_Dan)[,2], n.t_Dan, n.t_Dan))

M.cor_Dan <- matrix(colMeans(post.dist.m.cor_Dan),n.t_Dan, n.t_Dan)

mean_plus_HPD_M_Dan <- list(lower = matrix(HPDinterval(post.dist.m.cor_Dan)[,1], n.t_Dan, n.t_Dan), mean = M.cor_Dan, upper = matrix(HPDinterval(post.dist.m.cor_Dan)[,2], n.t_Dan, n.t_Dan))

