data {
  int<lower=1> nSubjects;
  int<lower=1> nTrials;
  real<lower=0,upper=1> Feedback[nSubjects, nTrials]; //Feedback on the current trial
  real<lower=0, upper=9999> ProbeResp[nSubjects, nTrials]; //Probe response
}

transformed data {

}

parameters {

  real<lower=0,upper=1> v_init_mu;
  
  real<lower=0> v_init_sd;

  real<lower=0,upper=1> v_init[nSubjects];
  
  real<lower=0,upper=0.25> noise_mu;
  real<lower=0,upper=0.2> noise_sd;
  real<lower=0,upper=0.5> noise[nSubjects];
}

model {
  v_init_sd  ~ cauchy(0,1);
  noise_sd  ~ cauchy(0,1);

  
  // give the prior here: how individual-level parameters are connected to the group-level parameters
  v_init ~ normal(v_init_mu, v_init_sd);
  noise ~ normal(noise_mu, noise_sd);

  
  
  for (s in 1:nSubjects) {
    real v; 
    real pe;    
    v = v_init[s];

    for (t in 1:nTrials) {
      v =  v;


      if (ProbeResp[s,t]!=9999){
        ProbeResp[s,t] ~ normal(v,noise[s]);
      }
      else{
        ProbeResp[s,t] ~ normal(9999,0.0000000000000001); 
      }
      
    }
  }
}

