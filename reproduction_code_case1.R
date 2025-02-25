###### package ######"

library(foreach)
library(doParallel)
library(splines)
library(lcmm)
library(survival)
library(MASS)
library(lme4)
library(dplyr)

###### Data load ############

Donnees_expo_cas1 <- read.delim("C:/Users/mlf/Desktop/Doc_stage/Code_m_Wagner/reproduction_code/Donnees_expo_cas1.txt")
donnees_outcome_cas1 <- base_fin <- read.delim("C:/Users/mlf/Desktop/Doc_stage/Code_m_Wagner/reproduction_code/base_fin.txt")

#########################################
############# Parameters ################
#########################################

# individual number
nb_ind = 1000
nb_ind = length(unique(as.factor(Donnees_expo_cas1$id)))





########################################
######## exposition estimation #########
########################################

# Model to estimate exposition
## Linear mixed effects model for the exposure
## Time is taken into account in splines (random)
Donnees_expo_cas1$edu = as.factor(Donnees_expo_cas1$edu)
m_bmi        = hlme(BMI_obs ~ edu+age0+B1+B2+B3+B4+B5, random = ~B1+B2+B3+B4+B5, subject = "id", data = Donnees_expo_cas1)

# Don't necessite this
m_bmi$conv #converge ok
#if((z==1) & (m_bmi$conv==2)) 
#{ z = 1   
#} else if((z>1) & (m_bmi$conv==2))
#{ z = z + 1  
#} else if (m_bmi$conv == 1)
#{
  
# New tab to put the estimation
  tab_bmi      = data.frame(id = rep(seq(1,nb_ind), each = 25),int = rep(1,25*nb_ind))
  tab_bmi$time = rep(seq(-24,0), nb_ind)
  covar        = subset(Donnees_expo_cas1, select=c(id, edu, age0), numvis == 1)
  tab_bmi2     = merge(tab_bmi, covar, by = "id", all = T)
  
  ## splines (without intercept)
  K1  = quantile(tab_bmi2$time,probs = c(0.20),na.rm=T)
  K2  = quantile(tab_bmi2$time,probs = c(0.40),na.rm=T)
  K3  = quantile(tab_bmi2$time,probs = c(0.60),na.rm=T)
  K4  = quantile(tab_bmi2$time,probs = c(0.80),na.rm=T)
  b5  = quantile(tab_bmi2$time,probs = c(0.05),na.rm=T)
  b95 = quantile(tab_bmi2$time,probs = c(0.95),na.rm=T)    
  splines           = as.matrix(ns(tab_bmi2$time, knots = c(K1,K2,K3,K4), Boundary.knots = c(b5 , b95)))
  colnames(splines) = c("B1","B2","B3","B4","B5")
  
  # merge the splines
  tab_bmi2               = cbind(tab_bmi2, splines)
  
  ## merge the coefficients of the lmm by individual (random effect)
  colnames(m_bmi$predRE) = c("id", "intBMI", "RE_B1", "RE_B2", "RE_B3", "RE_B4", "RE_B5")
  tab_bmi3               = merge(tab_bmi2, m_bmi$predRE, by="id", all=T)
  
  ## Computation of individual predictions of exposure
  tab_bmi3$edu      = as.numeric(tab_bmi3$edu)
  tab_bmi3$BMI_pred = m_bmi$best[1] + tab_bmi3$edu*m_bmi$best[2] + 
    tab_bmi3$age0*m_bmi$best[3] + 
    tab_bmi3$B1*m_bmi$best[4] +
    tab_bmi3$B2*m_bmi$best[5] + 
    tab_bmi3$B3*m_bmi$best[6] + 
    tab_bmi3$B4*m_bmi$best[7] + 
    tab_bmi3$B5*m_bmi$best[8] + 
    tab_bmi3$intBMI + 
    tab_bmi3$B1*tab_bmi3$RE_B1 +  
    tab_bmi3$B2*tab_bmi3$RE_B2 + 
    tab_bmi3$B3*tab_bmi3$RE_B3 + 
    tab_bmi3$B4*tab_bmi3$RE_B4 +  
    tab_bmi3$B5*tab_bmi3$RE_B5
  
  ###### Okay jusque là !!
  
  
  
  
  ############ Computation of the history exposure #######################
  K11 = quantile(tab_bmi3$time,probs = c(0.33),na.rm=T)	
  K22 = quantile(tab_bmi3$time,probs = c(0.66),na.rm=T)
  b5  = quantile(tab_bmi3$time,probs = c(0.05),na.rm=T)
  b95 = quantile(tab_bmi3$time,probs = c(0.95),na.rm=T)
  tab_bmi3 = tab_bmi3[order(tab_bmi3[,1], tab_bmi3[,3]), ]
  
  ## splines recompile
  B2K      = as.matrix(ns(tab_bmi3$time,knots = c(K11, K22), Boundary.knots = c(b5, b95), intercept = T)) 
  nz       = ncol(B2K)
  tab      = tab_bmi3
  
  ## add the new splines column 
  tab[,((ncol(tab_bmi3)+1):(ncol(tab_bmi3)+nz))] = B2K
  colnames(tab) = c(colnames(tab)[1:(ncol(tab)-nz)],paste("H",1:nz,sep = "")) #rename the splines column
  
  
  ## faire la prédiction * les nouveaux splines (principe du WCIE pour 
  # récupérer le poids à chaque temps)
  tab$coco1 = tab$BMI_pred*tab$H1
  tab$coco2 = tab$BMI_pred*tab$H2
  tab$coco3 = tab$BMI_pred*tab$H3
  tab$coco4 = tab$BMI_pred*tab$H4
  
  tab$Hist_K1 = 0 
  tab$Hist_K2 = 0 
  tab$Hist_K3 = 0 
  tab$Hist_K4 = 0
  
  
######## Ici calcul la somme des prédiction * splines pour chaques noeuds ##########
  
  ## H1
  i = 1   ## lines
  j = -24 ## times
  k = 1 ## subject
  x = 0 ## sum 
  for (i in 1:dim(tab)[1]){
    if (tab$id[i]==k){  #vérifie qu'il s'agit du bon individu
      x = x + tab$coco1[i] #somme des (prédiction * splines)
      tab$Hist_K1[i] = x   #affect la valeur à la ligne correspondante
      i = i+1  
    }else{  # permet de changer d'individu
      x = tab$coco1[i]
      tab$Hist_K1[i] = x #remet à 0 le compteur quand chgmt d'ind
      k=k+1}
  }  
  
  ## H2
  i = 1   
  j = -24 
  k = 1 
  x = 0  
  for (i in 1:dim(tab)[1]){
    if (tab$id[i]==k){
      x = x + tab$coco2[i]
      tab$Hist_K2[i] = x
      i = i+1
    }else{
      x = tab$coco2[i]
      tab$Hist_K2[i] = x
      k=k+1
    }
  }
  
  ## H3
  i = 1   
  j = -24 
  k = 1 
  x = 0 
  for (i in 1:dim(tab)[1]){
    if (tab$id[i]==k){
      x = x + tab$coco3[i]
      tab$Hist_K3[i] = x
      i = i+1
    }else{
      x = tab$coco3[i]
      tab$Hist_K3[i] = x
      k=k+1}
  }  
  
  ## H4
  i = 1   
  j = -24 
  k = 1 
  x = 0 
  for (i in 1:dim(tab)[1]){
    if (tab$id[i]==k){
      x = x + tab$coco4[i]
      tab$Hist_K4[i] = x
      i = i+1
    }else{
      x = tab$coco4[i]
      tab$Hist_K4[i] = x
      k=k+1}
  }  
  
  # récupère les splines à t0
  try         = subset(tab, select=c(id,Hist_K1:Hist_K4), time==0)
  
  # d'où sort le w<sq ?
  base_final2 = merge(donnees_outcome_cas1, try, by="id", all=T)#w<sq
  

  ## Linear mixed effects model for the outcome 
  base_final2$primo = as.factor(base_final2$primo)
  base_final2$edu   = as.factor(base_final2$edu)    
  m_2K_acute = hlme(TICS_obs ~ time + primo + edu + age0 +
                      Hist_K1 + Hist_K2 + Hist_K3+ Hist_K4 +
                      time:(edu + age0 + Hist_K1 + Hist_K2 + Hist_K3+ Hist_K4), random = ~time, subject = "id", data = base_final2)
  
  ## dataset with estimated parameters of the simulation replicate z
  Bestmod           = as.data.frame(t(estimates(m_2K_acute)))
  Bestmod$conv_tics = m_2K_acute$conv
  Bestmod$iter_tics = m_2K_acute$niter
  Bestmod$conv_bmi  = m_bmi$conv
  Bestmod$iter_bmi  = m_bmi$niter
  # Bestmod$simu      = z   numéro de la simu pas besoin pour nous
  
  # matrice de variance covariance
  Vmod = VarCov(m_2K_acute)
  
  # ecart type 
  sd = data.frame(t(sqrt(diag(Vmod))))
  
  # matrice de variance covariance de l'historique des splines
  varCov_Hist      = as.matrix(Vmod[c("Hist_K1","Hist_K2","Hist_K3","Hist_K4"),c("Hist_K1","Hist_K2","Hist_K3","Hist_K4")])
  # matrice de varcovar de l'intéraction temps avec splines
  varCov_Hist_time = as.matrix(Vmod[c("time:Hist_K1","time:Hist_K2","time:Hist_K3","time:Hist_K4"),c("time:Hist_K1","time:Hist_K2","time:Hist_K3","time:Hist_K4")])
  
  
  ### Computation of effects (init + slope) every year 
  ## intercept
  effect$eff_H1 = m_2K_acute$best["Hist_K1"]
  effect$eff_H2 = m_2K_acute$best["Hist_K2"]
  effect$eff_H3 = m_2K_acute$best["Hist_K3"]
  effect$eff_H4 = m_2K_acute$best["Hist_K4"]
  ## slope
  effect$eff_H1_int = m_2K_acute$best["time:Hist_K1"]
  effect$eff_H2_int = m_2K_acute$best["time:Hist_K2"]
  effect$eff_H3_int = m_2K_acute$best["time:Hist_K3"]
  effect$eff_H4_int = m_2K_acute$best["time:Hist_K4"]
  ## vecteur B1(t), B2(t), B3(t)
  vecH = as.matrix(subset(effect, select=c("H1","H2","H3","H4")))
  
  
  #### Calcul des effets + IC95% avec la méthode Delta classique (sans bootstrap)    
  ## intercept
  effect$eff_int    = effect$H1*effect$eff_H1 + effect$H2*effect$eff_H2 + effect$H3*effect$eff_H3 + effect$H4*effect$eff_H4
  effect$sd_eff_int = 0
  for (b in 1:dim(effect)[1])
  {
    effect$sd_eff_int[b] = sqrt(t(vecH[b,]) %*% varCov_Hist %*% vecH[b,])
  }
  
  ## slope
  effect$eff_slope    = effect$H1*effect$eff_H1_int + effect$H2*effect$eff_H2_int + effect$H3*effect$eff_H3_int + effect$H4*effect$eff_H4_int
  effect$sd_eff_slope = 0
  for (c in 1:dim(effect)[1])
  {
    effect$sd_eff_slope[c] = sqrt(t(vecH[c,]) %*% varCov_Hist_time %*% vecH[c,])
  }
  
  ## theoretical effects
  effect$eff_int_THEO   = effect$eff_int_norm
  effect$eff_slope_THEO = effect$eff_slope_norm
  
  box      = subset(effect, select=c(time_retro,eff_int,sd_eff_int, eff_slope,sd_eff_slope,eff_int_THEO,eff_slope_THEO))
  # box$simu = z pas nécessaire comme pas de simu
  
  box$binf_int  = box$eff_int   - 1.96*box$sd_eff_int
  box$bsup_int  = box$eff_int   + 1.96*box$sd_eff_int
  box$binf_slop = box$eff_slope - 1.96*box$sd_eff_slope
  box$bsup_slop = box$eff_slope + 1.96*box$sd_eff_slope
  
  # Indicateur - coverage (0/1)
  box$ind_int   = 1
  box$ind_int   = ifelse(box$bsup_int < box$eff_int_THEO, 0, 1)
  box$ind_int   = ifelse(box$binf_int > box$eff_int_THEO, 0, box$ind_int)
  box$ind_slope = 1
  box$ind_slope = ifelse(box$bsup_slop < box$eff_slope_THEO, 0, 1)
  box$ind_slope = ifelse(box$binf_slop > box$eff_slope_THEO, 0, box$ind_slope)
  
  # View(box)
  
  
  if (m_2K_acute$conv==1)
  {
    #############################
    ####       BOOTSTRAPS   #####
    #############################
    m       = m_bmi             
    tab_MED = Donnees_expo_cas1  
    tab_MED$int = 1     
    tab_MED$edu = as.numeric(tab_MED$edu) 
    mu    = as.matrix(estimates(m)) 
    # View(mu) View(Sigma)
    Sigma = as.matrix(VarCov(m)) 
    boot  = as.matrix(mvrnorm(nb_boot, mu, Sigma))
    #View(boot) View(tab_pred) 
    ## Calcul des effets aléatoires
    parb = ncol(boot)        
    Xrandom = as.matrix(m$predRE)
    Xrandom[,c(2:ncol(Xrandom))] = NA
    
    tab_pred=data.frame(time = rep(-24,0, length = 25))
    ## View()
    #ii=1
    #jj=1
    for (ii in 1:nb_boot)  ## boucle bootstrap
    {
      pos  = grep("^cholesky", colnames(boot))
      matC = matrix(0,nrow=6,ncol=6)
      matC[,1]       = boot[ii,pos[1]]
      matC[,2]       = boot[ii,pos[c(2,3)]]
      matC[,3]       = boot[ii,pos[c(4,5,6)]]
      matC[c(1:4),4] = boot[ii,pos[c(7,8,9,10)]]
      matC[c(1:5),5] = boot[ii,pos[c(11,12,13,14,15)]]
      matC[c(1:6),6] = boot[ii,pos[c(16,17,18,19,20,21)]]
      matC[lower.tri(matC)] = 0
      matB = t(matC)%*% matC
      beta_b = boot[ii,1:dim(as.matrix(fixef(m)[[2]]))[1]]
      matH = matrix(nrow=nb_ind,ncol=5)
      colnames(matH) = c("id","chiB1_2K_0","chiB2_2K_0","chiB3_2K_0","chiB4_2K_0")
      
      # View(Xrandom) View(tab_MED) View(matZ)
      for (jj in 1:nb_ind)  ## boucle sujet
      {
        matZ = as.matrix(subset(tab_MED, id==Xrandom[jj,1] & is.na(BMI_obs)==F, select=c("int","B1","B2","B3","B4","B5")))
        ni   = nrow(matZ)
        pos  = grep("^stder", colnames(boot))
        
        matE = as.matrix(diag(ni)*boot[ii,pos])
        matY = as.matrix(subset(tab_MED, id==Xrandom[jj,1] & is.na(BMI_obs)==F, select=c("BMI_obs")))
        matX = subset(tab_MED, id==Xrandom[jj,1] & is.na(BMI_obs)==F, select=c("int","edu","age0","B1","B2","B3","B4","B5"))
        matXbeta = as.matrix(apply(matX,1,function(x) sum(x*beta_b)))
        
        V  = matZ%*%matB%*%t(matZ)+matE
        RE = t(matB%*%t(matZ)%*%(solve(V))%*%(matY-matXbeta))
        Xrandom[jj,c(2:ncol(Xrandom))] = RE
        # head(Xrandom) head(matX)
        
        ### Calcul des prédictions individuelles 
        tpred = data.frame(numvis=seq(1,25),time = seq(-24,0, length = 25),int=1, 
                           edu=matX$edu[1], age0=matX$age0[1], 
                           RE_int=Xrandom[jj,2], 
                           RE_Z1=Xrandom[jj,3],
                           RE_Z2=Xrandom[jj,4],
                           RE_Z3=Xrandom[jj,5],
                           RE_Z4=Xrandom[jj,6],
                           RE_Z5=Xrandom[jj,7])
        
        Z  = as.matrix(ns(tpred$time,knots = c(-19,-14,-10,-5), Boundary.knots = c(-23,-1)))
        nz = ncol(Z)
        tpred[,((ncol(tpred)+1):(ncol(tpred)+nz))] = Z
        colnames(tpred) = c(colnames(tpred)[1:(ncol(tpred)-nz)],paste("Z",1:nz,sep = ""))
        
        ### COMPUTATION OF INDIVIDUAL PREDICTIONS OF BMI
        tpred$predBMI = beta_b[1] + 
          beta_b[2]*tpred$edu +
          beta_b[3]*tpred$age0 + 
          beta_b[4]*tpred$Z1 + 
          beta_b[5]*tpred$Z2 + 
          beta_b[6]*tpred$Z3 + 
          beta_b[7]*tpred$Z4 + 
          beta_b[8]*tpred$Z5 + 
          tpred$RE_int + 
          tpred$RE_Z1*tpred$Z1+ 
          tpred$RE_Z2*tpred$Z2+ 
          tpred$RE_Z3*tpred$Z3+ 
          tpred$RE_Z4*tpred$Z4+ 
          tpred$RE_Z5*tpred$Z5
        
        ## CALCUL DE L'HISTOIRE DE L'EXPOSITION
        ##################################################################
        ####            Base de splines  -   2  KNOTS                 ####
        ##################################################################
        ####################
        ####  -24 to 0  ####
        ####################
        ### BASE DE SPLINES 
        quantile(tpred$time,probs = c(0.05,0.33,0.66,0.95),na.rm=T)	# -23.2 ; 0
        B2K = as.matrix(ns(tpred$time,knots = c(-16.08,-8.16),Boundary.knots = c(-22.8,-1.2), intercept = T)) 
        tpred[,((ncol(tpred)+1):(ncol(tpred)+ncol(B2K)))] = B2K
        colnames(tpred) = c(colnames(tpred)[1:(ncol(tpred)-ncol(B2K))],paste("BH",1:ncol(B2K),sep = ""))
        
        ### Somme des histoires
        tpred$chiB1_2K_0 = sum(tpred$predBMI*tpred$BH1)
        tpred$chiB2_2K_0 = sum(tpred$predBMI*tpred$BH2)
        tpred$chiB3_2K_0 = sum(tpred$predBMI*tpred$BH3)
        tpred$chiB4_2K_0 = sum(tpred$predBMI*tpred$BH4)
        
        matH[jj,1] = Xrandom[jj,1]
        matH[jj,2] = as.matrix(tpred$chiB1_2K_0[1])
        matH[jj,3] = as.matrix(tpred$chiB2_2K_0[1])
        matH[jj,4] = as.matrix(tpred$chiB3_2K_0[1])
        matH[jj,5] = as.matrix(tpred$chiB4_2K_0[1])
        ## sujet suivant
        jj = jj + 1
      }
      
      
      ## FUSION - merge
      endBoot     = merge(base_final2,matH, by="id",all.x=T)
      endBoot$edu = as.factor(endBoot$edu)
      
      ## -24 ans à 0
      m_2KBoot = hlme(TICS_obs ~ time+primo+edu+age0+
                        edu*time+age0*time +
                        chiB1_2K_0 + chiB2_2K_0 + chiB3_2K_0 +chiB4_2K_0 +
                        time:(chiB1_2K_0 + chiB2_2K_0 + chiB3_2K_0 +chiB4_2K_0),
                      random = ~time, subject = "id", data = endBoot)
      
      ### Récupérer Paramètres estimés
      bestBoot0       = as.data.frame(t(estimates(m_2KBoot)))
      bestBoot0$boots = ii
      bestBoot0$conv  = m_2KBoot$conv
      
      #### Récupérer Variance empirique (issues du summary)
      VmodBoot = VarCov(m_2KBoot)
      sdBoot0  = data.frame(boots=ii,t(sqrt(diag(VmodBoot))),conv=m_2KBoot$conv)
      
      
      ################
      ##  2 KNOTS   ##
      ################
      final = matrix(nrow=25,ncol=14)
      rownames(final) = c(-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0)
      colnames(final) = c("B1","B2","B3","B4","chiB1","chiB2","chiB3","chiB4","chiB1slop","chiB2slop","chiB3slop","chiB4slop","effectU","effectUslop")
      # splines
      datnew <- data.frame(time = seq(-24,0, length = 25)) 
      #quantile(datnew$time, probs=c(0.05,0.33,0.66,0.95))
      B <- as.matrix(ns(datnew$time,knots = c(-16.08,-8.16),Boundary.knots = c(-22.8,-1.20), intercept = T))  # boundaries à 95%
      final[,c(1:4)]= B
      ## Effect alpha(u), gamma(u)
      ### ChiBk
      final[,5] = m_2KBoot$best[6]
      final[,6] = m_2KBoot$best[7]
      final[,7] = m_2KBoot$best[8]
      final[,8] = m_2KBoot$best[9]
      ### ChiBkSlop
      final[,9]  = m_2KBoot$best[12]
      final[,10] = m_2KBoot$best[13]
      final[,11] = m_2KBoot$best[14]
      final[,12] = m_2KBoot$best[15]
      ### Calcul des effects
      for (kk in 1:dim(final)[1])
      {
        final[kk,13] = final[kk,1]*final[kk,5]+final[kk,2]*final[kk,6]+final[kk,3]*final[kk,7]+final[kk,4]*final[kk,8]
        final[kk,14] = final[kk,1]*final[kk,9]+final[kk,2]*final[kk,10]+final[kk,3]*final[kk,11]+final[kk,4]*final[kk,12]
      }
      
      boxBoot0  = data.frame(simu=z,boots=ii,conv=m_2KBoot$conv, time=seq(-24,0,length=25), eff_UBoot = final[,13], eff_SBoot = final[,14])
      
      
      if (ii==1)
      { 
        bestBoot   = bestBoot0     # paramètres individuels estimés du modèle boostrap
        boxBoot    = boxBoot0      # effets de l'int et slope estimés (et calculés) à partir du modèle boostrap
        sdBoot     = sdBoot0       
        compt_conv = m_2KBoot$conv
        var_intra  = VarCov(m_2KBoot)    ### conservation de la matrix VarCov pour b=1
      } else if(ii>1) 
      {
        boxBoot    = rbind(boxBoot,boxBoot0)
        bestBoot   = rbind(bestBoot,bestBoot0)
        sdBoot     = rbind(sdBoot,sdBoot0)
        compt_conv = compt_conv + m_2KBoot$conv
        var_intra  = var_intra + VarCov(m_2KBoot)  ### sum des 500 matrix VarCov
      }
      ### COMPTEUR
      ii=ii+1
    }
    
    #View(var_intra) View(boxBoot) View(bestBoot) View(boxBoot)
    ##########################
    ###  END OF BOOTSTRAP  ###
    ##########################
    
    ### MEAN EFFECT FOR B BOOTSTRAPS View(brice)
    box$eff_intB   = 0
    box$eff_slopeB = 0
    for (i in -24:0)
    {
      brice          = subset(boxBoot,time==i)
      box$eff_intB   = ifelse(box$time==i,mean(brice$eff_UBoot),box$eff_intB)
      box$eff_slopeB = ifelse(box$time==i,mean(brice$eff_SBoot),box$eff_slopeB)
      i = i + 1
    }
    
    ### VARIANCE DECOMPOSITION FORMULA (used in the construction of 95%CI) ###
    ## FOR PARAMETERS
    ## Intra-individual variability  
    var_intra = var_intra/nb_boot
    
    
    ## Inter-individual variability 
    vec_par   = t(as.matrix(bestBoot))
    mean_par  = as.matrix(apply(vec_par, 1, function(x) mean(x)))
    var_tempo = matrix(0,nrow=dim(vec_par)[1], ncol=dim(vec_par)[1])
    var_inter = matrix(0,nrow=dim(vec_par)[1], ncol=dim(vec_par)[1])
    for (i in 1:dim(vec_par)[2])
    {
      var_tempo = (vec_par[,i]-mean_par)%*%t(vec_par[,i]-mean_par)
      var_inter = var_inter + var_tempo
      i = i + 1 
    }
    
    var_inter = var_inter[-c(20,21,22),-c(20,21,22)]
    var_inter = var_inter/nb_boot
    
    ### VARIANCE TOTALE PARAMETERS
    var_tot = var_inter + var_intra
    
    
    
    ### FOR TIME-VARYING EFFECTS
    # Intra-individual variability (sd)
    eff_varCov_tot      = var_intra[c("chiB1_2K_0","chiB2_2K_0","chiB3_2K_0","chiB4_2K_0"),c("chiB1_2K_0","chiB2_2K_0","chiB3_2K_0","chiB4_2K_0")]
    eff_varCov_tot_time = var_intra[c("time:chiB1_2K_0","time:chiB2_2K_0","time:chiB3_2K_0","time:chiB4_2K_0"),c("time:chiB1_2K_0","time:chiB2_2K_0","time:chiB3_2K_0","time:chiB4_2K_0")]
    # intercept
    box$var_intra_effI   = 0
    for (bb in 1:dim(box)[1])
    {
      box$var_intra_effI[bb] = t(vecH[bb,]) %*% eff_varCov_tot %*% vecH[bb,]
    }
    
    #slope
    box$var_intra_effS = 0
    for (cc in 1:dim(box)[1])
    {
      box$var_intra_effS[cc] = t(vecH[cc,]) %*% eff_varCov_tot_time %*% vecH[cc,]
    }
    
    ## Inter-individual variability
    ## INTERCEPT
    pos1 = which(colnames(bestBoot) == "chiB1_2K_0")
    pos2 = which(colnames(bestBoot) == "chiB2_2K_0")
    pos3 = which(colnames(bestBoot) == "chiB3_2K_0")
    pos4 = which(colnames(bestBoot) == "chiB4_2K_0")
    
    par_Int       = t(as.matrix(bestBoot[c(pos1,pos2,pos3,pos4)]))
    mean_par      = as.matrix(apply(par_Int,1,function(x) mean(x)))
    var_tempo     = matrix(0,nrow=dim(par_Int)[1], ncol=dim(par_Int)[1])
    var_inter_Int = matrix(0,nrow=dim(par_Int)[1], ncol=dim(par_Int)[1])
    
    for (i in 1:dim(par_Int)[2])
    {
      var_tempo = (par_Int[,i]-mean_par)%*%t(par_Int[,i]-mean_par)
      var_inter_Int = var_inter_Int + var_tempo
      i = i + 1 
    }
    
    var_inter_Int = var_inter_Int/nb_boot
    
    box$var_inter_effI = 0
    for (bb in 1:dim(box)[1])
    {
      box$var_inter_effI[bb] = t(vecH[bb,]) %*% var_inter_Int %*% vecH[bb,]
    }
    
    ## SLOPE
    pos1 = which(colnames(bestBoot) == "time:chiB1_2K_0")
    pos2 = which(colnames(bestBoot) == "time:chiB2_2K_0")
    pos3 = which(colnames(bestBoot) == "time:chiB3_2K_0")
    pos4 = which(colnames(bestBoot) == "time:chiB4_2K_0")
    
    par_Slope       = t(as.matrix(bestBoot[c(pos1,pos2,pos3,pos4)]))
    mean_par        = as.matrix(apply(par_Slope,1,function(x) mean(x)))
    var_tempo       = matrix(0,nrow=dim(par_Slope)[1], ncol=dim(par_Slope)[1])
    var_inter_Slope = matrix(0,nrow=dim(par_Slope)[1], ncol=dim(par_Slope)[1])
    
    for (i in 1:dim(par_Slope)[2])
    {
      var_tempo = (par_Slope[,i]-mean_par)%*%t(par_Slope[,i]-mean_par)
      var_inter_Slope = var_inter_Slope + var_tempo
      i = i + 1 
    }
    
    var_inter_Slope = var_inter_Slope/nb_boot
    
    box$var_inter_effS = 0
    for (bb in 1:dim(box)[1])
    {
      box$var_inter_effS[bb] = t(vecH[bb,]) %*% var_inter_Slope %*% vecH[bb,]
    }
    
    
    ### TOTALE VARIANCE - TIME-VARYING EEFECTS 
    box$var_tot_effI = box$var_inter_effI + box$var_intra_effI
    box$var_tot_effS = box$var_inter_effS + box$var_intra_effS
    
    box$sd_eff_intB   = sqrt(box$var_tot_effI)
    box$sd_eff_slopeB = sqrt(box$var_tot_effS)
    
    ## Confidence Intervals 
    box$binf_intB  = box$eff_intB   - 1.96*box$sd_eff_intB
    box$bsup_intB  = box$eff_intB   + 1.96*box$sd_eff_intB
    box$binf_slopB = box$eff_slopeB - 1.96*box$sd_eff_slopeB
    box$bsup_slopB = box$eff_slopeB + 1.96*box$sd_eff_slopeB
    box$conv_ticsB = compt_conv
    
    ## Coverage (0/1) ?
    ## Initial level 
    box$ind_intB = 1
    box$ind_intB = ifelse(box$bsup_intB < box$eff_int_THEO, 0, 1)
    box$ind_intB = ifelse(box$binf_intB > box$eff_int_THEO, 0, box$ind_intB)
    ## Slope
    box$ind_slopeB = 1
    box$ind_slopeB = ifelse(box$bsup_slopB < box$eff_slope_THEO, 0, 1)
    box$ind_slopeB = ifelse(box$binf_slopB > box$eff_slope_THEO, 0, box$ind_slopeB)
    
  } 
  
  if (z==1){
    best    = Bestmod
    sd_end  = sd
    box_end = box
    bestB   = bestBoot
    sdB     = sdBoot
    
  } else if(z>1) 
  {
    best    = rbind(best,Bestmod)
    sd_end  = rbind(sd_end,sd)
    box_end = rbind(box_end,box)
    bestB   = rbind(bestB,bestBoot)
    sdB     = rbind(sdB,sdBoot)
  }
  
  print(z)

#}