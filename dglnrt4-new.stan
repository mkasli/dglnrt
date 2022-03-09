data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n_obs;                  
  int <lower=1> ind_person_obs[n_obs];  
  int <lower=1> ind_item_obs[n_obs];    
  int <lower=0,upper=1> i_status_obs[n_obs];
  real Y[n_obs];                       
}

parameters {
  vector[I] tau_t;            
  vector[I] tau_c;  
  vector<lower = 0>[J] beta;         
  vector <lower=0.5, upper = 1> [J]  inv_alpha;              
  real<offset = 4, multiplier = 1> mu1;
  real<lower=0, upper = 3> sigma1;
  real<lower=0> sigma_t;
  real<lower=0> sigma_c;
}

transformed parameters {
  vector[I] T;

  for (i in 1:n_obs) {
    if(tau_t[ind_person_obs[i]]>tau_c[ind_person_obs[i]])
      
      T[ind_person_obs[i]] = 0;
    
    else 
      
      T[ind_person_obs[i]] = 1;
  }
}


model{
  
  sigma_t ~ exponential(1);
  sigma_c ~ exponential(1);
  
  tau_t    ~ normal(0,sigma_t);
  tau_c    ~ normal(0,sigma_c);
  
  mu1      ~ normal(4,1);
  sigma1   ~ exponential(1);
    beta     ~ normal(mu1,sigma1);

    inv_alpha    ~ gamma(250,336);
  
  for (i in 1:n_obs) {
   real p = 0;
   if (i_status_obs[i] == 0) {
     real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
     p = p_t;
   } else {
    real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
    p = p_c;
   }
    Y[i] ~ normal(p, inv_alpha[ind_item_obs[i]]);
  }
}

generated quantities {      

  vector <lower=0> [J]  alpha = inv(inv_alpha); 
}

