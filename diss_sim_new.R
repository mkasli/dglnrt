

sim.conditions <- expand.grid(n.person= c(500),
                              n.item  = c(50,100),
                              p.item  = c(.0,.2,.4,.6),
                              p.person= c(.05,.1))

# item preknowledge effect: 0.4, 0.7, 1.0 and 1.5
# these values will be arranged based on the extracted avg parameters

#condition is gonna change 1 to 90
    condition = 4
    
    n.person = sim.conditions[condition,]$n.person
    n.item   = sim.conditions[condition,]$n.item
    p.item   = sim.conditions[condition,]$p.item
    p.person = sim.conditions[condition,]$p.person


#### bu kisimdan sonrasi update olacak

# A generic function to simulate DG-LNRT data with no item preknowledge

    sim_dglnrt <- function() {
      
      # MODEL PARAMETERS
      
      #averaging parameters from 2 datasets: 
      # alpha (1.42, 0.37), alpha (2.06, 0.28) = (1,74,0.33)
      # beta (4.19, 0.38), beta (3.98,0.32) = (4.08,0.35)
      
      
      beta  <- rnorm(n.item,4.00,0.35)
      alpha <- rnorm(n.item,1.74,0.33)
      
      # Tau for unflagged examinees
      
      cor0 <- matrix(c(1,.90,.90,1),2,2)
      tau0 <- mvrnorm(n.person*(1-p.person),c(0,0),cor2cov(cor0,c(0.17,0.16)))
      
      # Tau for flagged examinees 
      
      cor1 <- matrix(c(1,.74,.74,1),2,2)
      
      # When Item preknowledge effect = 0.4
      
      tau1 <- mvrnorm(n.person*p.person,c(0.00,0.122),cor2cov(cor1,c(0.45,0.38)))
      
      # When Item preknowledge effect = 0.7
      #tau1 <- mvrnorm(n.person*p.person,c(0.00,0.214),cor2cov(cor1,c(0.45,0.38)))
      
      # When Item preknowledge effect = 1.0
      # tau1 <- mvrnorm(n.person*p.person,c(0.00,0.306),cor2cov(cor1,c(0.45,0.38)))
      
      # When Item preknowledge effect = 1.5
      # tau1 <- mvrnorm(n.person*p.person,c(0.00,0.459),cor2cov(cor1,c(0.45,0.38)))
      
      
      # A vector for item status (0: not disclosed, 1:disclosed)
      
      C    <- c(rep(0,n.item*(1-p.item)),rep(1,n.item*p.item))
      
      # Shuffle the vector, randomly selected 91 items are compromised
      
      C <- C[base::sample(1:n.item,n.item)]
      
      # RESPONSE TIME GENERATION
      
      # Note that the response time data generated is
      # already on the log scale
      
      # Unflagged Examinees
      
      rt0 <- matrix(nrow = n.person*(1-p.person), ncol = n.item)
      
      for (i in 1:n.person*(1-p.person)) {
        for (j in 1:n.item) {
          p_t = beta[j] - tau0[i, 1]
          p_c = beta[j] - tau0[i, 2]
          p   = p_t * (1 - C[j]) + p_c * C[j]
          rt0[i, j] = rnorm(1, p, 1 / alpha[j])
        }
      }
      
      # Flagged Examinees
      
      rt1 <- matrix(nrow = n.person*p.person, ncol = n.item)
      
      for (i in 1:n.person*p.person) {
        for (j in 1:n.item) {
          p_t = beta[j] - tau1[i, 1]
          p_c = beta[j] - tau1[i, 2]
          p   = p_t * (1 - C[j]) + p_c * C[j]
          rt1[i, j] = rnorm(1, p, 1 / alpha[j])
        }
      }
    
    # Combine the groups
    
    rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
                cbind(data.frame(exp(rt1)), gr = 2))
    
    
    # Random shuffle of examinees
    
    rt <- rt[sample(1:nrow(rt),nrow(rt),replace = FALSE),]
      
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      tau_t = c(tau0[, 1], tau1[, 1]),
      tau_c = c(tau0[, 2], tau1[, 2]),
      C = C
    ))
  }

##############################################################################

mod <- cmdstan_model(here::here('diss_simulation/dglnrt2.stan'))    
    
    
    data <- sim_dglnrt()
    
    d.sub <- data$rt[,1:n.item]
    
    d.sub$ID <- 1:n.person
    
    d.long <- reshape(
      data        = d.sub,
      idvar       = "ID",
      varying     = list(colnames(d.sub)[1:n.item]),
      timevar     = "Item",
      times       = 1:n.item,
      v.names      = "RT",
      direction   = "long"
    )
    
    d.long <- na.omit(d.long)
    d.long$logRT <- log(d.long$RT)
    
    d.long$i.status <- NA
    
    flagged.item <- which(data$C==1)
    
    d.long[which(d.long$Item%in%flagged.item==TRUE),]$i.status = 1
    d.long[which(d.long$Item%in%flagged.item==FALSE),]$i.status = 0
    
    # Data_rt: Object for Stan
    
    data_rt <- list(
      J              = n.item,
      I              = n.person,
      n_obs          = length(d.long$RT),
      ind_person_obs = d.long$ID,
      ind_item_obs   = d.long$Item,
      i_status_obs   = d.long$i.status,
      Y              = d.long$logRT
    )

# Model Fitting (cmdstan)
    
    fit <- mod$sample(
      data = data_rt,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      iter_warmup   = 500,
      iter_sampling = 1500,
      refresh = 10,
      adapt_delta = 0.99)

# Save the output file   
    
 
    fit$cmdstan_summary()
    
    #stanfit <- rstan::read_stan_csv(fit$output_files())
    
    filename1 <- paste ("dglnrt","sim", sep ="")
    
    save.image(filename1)

save.image(here(''))

################################################################  

# inverse-gamma
#  alpha    ~ inv_gamma(800,1550);
#  alpha    ~ inv_gamma(n.person/2=v0, v0*a); a=mean(alpha0)

# beta estimates (Eq. 23)

#time<-as.data.frame(data_rt$Y)

#beta0 = colMeans(time,na.rm=TRUE)

# alpha estimates (Eq. 24)

# dglnrt_cmdstanr_v1 file indan alindi.
#alpha0 = sqrt(1/colMeans((log(d[,5:29]) - matrix(beta0,93,25,byrow=T))^2,na.rm=TRUE))

#alpha0 = sqrt(1/colMeans((time - matrix(beta0,n.person,n.item,byrow=T))^2,na.rm=TRUE))

# hyper parameters for inverse gamma 
# Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution


#N= n.person
#a = mean(alpha0)
#v0 = N/2
#v0
#v0*a

###################################################################################
