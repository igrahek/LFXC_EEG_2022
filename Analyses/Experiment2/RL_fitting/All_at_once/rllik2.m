function lik = rllik2(x,data,type)
    
    % Likelihood function for Q-learning on two-armed bandit with single
    % learning rate.
    %
    % USAGE: lik = rllik(x,data)
    %
    % INPUTS:
    %   x - parameters:
    %       x(1) - inverse temperature
    %       x(2) - learning rate
    %   data - structure with the following fields
    %          .c - [N x 1] choices
    %          .r - [N x 1] rewards
    %
    % OUTPUTS:
    %   lik - log-likelihood
    %
    % Sam Gershman, June 2015
    
v = x(1);  % initial values
var = x(2); % variance for the likelihood
lr = x(3); % positive learning rate

lik = 0;

for n = 1:data.N     % for every trial
    
    % Determine what we are fitting and take the appropriate subjective
    % estimate
    if strcmp(type,'reward')==1
        v_subj = data.RewRateProbeResp(n); % For reward subjective estimate
    elseif strcmp(type,'efficacy')==1
        v_subj = data.EfficacyProbeResp(n);  % For efficacy subjective estimate
    end
    
    % Take the actual reward and efficacy on the current trial
    r = data.IsRewarded(n); % For reward
    e = data.EffLvl(n); % For efficacy
    
    %------------1 learning rate RL-----------%


        if strcmp(type,'reward')==1 % if fitting reward
            rpe = r-v; % prediction error
            v = v + lr*rpe;      % update the value
        elseif strcmp(type,'efficacy')==1 % if fitting efficacy
            rpe = e-v; % prediction error
            v = v + lr*rpe;      % update the value
        end

    % calculate the likelihood, but only if there is a probe response
    % on this trial
    if isnan(v_subj) == 0             
        P = normpdf(v-v_subj,0,var);
        if isnan(P) || P==0; P = realmin; end % avoid NaNs and zeros in the logarithm
        lik = lik + log(P);

         
    else
        lik = lik;
    end
    
end

