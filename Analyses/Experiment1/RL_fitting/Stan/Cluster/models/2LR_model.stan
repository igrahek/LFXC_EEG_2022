data {
  int<lower=1> nSubjects;
  int<lower=1> nTrials;
  real<lower=0,upper=1> Feedback[nSubjects, nTrials]; //Feedback on the current trial
  real<lower=0, upper=9999> ProbeResp[nSubjects, nTrials]; //Probe response
}

transformed data {

}

parameters {
  real<lower=0,upper=1> lr_pos_mu;
  real<lower=0,upper=1> lr_neg_mu;

  real<lower=0> lr_pos_sd;
  real<lower=0> lr_neg_sd;

  real<lower=0,upper=1> lr_pos[nSubjects];
  real<lower=0,upper=1> lr_neg[nSubjects];
  
  real<lower=0,upper=1> v_init_mu;
  
  real<lower=0> v_init_sd;

  real<lower=0,upper=1> v_init[nSubjects];
  
  //real<lower=0,upper=0.25> noise_mu;
  //real<lower=0,upper=0.2> noise_sd;
  real<lower=0,upper=0.5> noise;
}

model {
  lr_pos_sd  ~ cauchy(0,1);
  lr_neg_sd  ~ cauchy(0,1);
  v_init_sd  ~ cauchy(0,1);
  //noise_sd  ~ cauchy(0,1);

  
  // give the prior here: how individual-level parameters are connected to the group-level parameters
  lr_pos ~ normal(lr_pos_mu, lr_pos_sd);
  lr_neg ~ normal(lr_neg_mu, lr_neg_sd);
  v_init ~ normal(v_init_mu, v_init_sd);
  //noise ~ normal(noise_mu, noise_sd);

  
  
  for (s in 1:nSubjects) {
    real v; 
    real pe;    
    v = v_init[s];

    for (t in 1:nTrials) {
      pe = Feedback[s,t] - v;
      
      //print("LR_pos=",lr_pos[s])
      //print("LR_neg=",lr_neg[s])
      //print("Value=",v)
      //print("Feedback=",Feedback[s,t])
      //print("PE=",pe)

      //print("Estimated response=",v)
      //print("True response=",ProbeResp[s,t])

      if (pe>=0){
        v =  v + lr_pos[s] * pe;

      }
      else if (pe<0){
        v = v + lr_neg[s] * pe;

      }
      
      if (ProbeResp[s,t]!=9999){
        ProbeResp[s,t] ~ normal(v,noise);
        
      }
      
      //print("likelihood=",target())

      //else{
        //ProbeResp[s,t] ~ normal(9999,0.0000000000000001); 
      //}
      
    }
  }
}
