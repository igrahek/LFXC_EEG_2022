function sse = cost(x,data,type,model)

% Fitting code by Sam Gershman, June 2015 (https://github.com/sjgershm/mfit);
% Adapted by Ivan Grahek (grahek.i@gmail.com) in order to fit a LR to the
% LFXC_EEG dataset

v = 0.5;  % initial values
SSE = 0;

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
    
    %------------Fixed intercept model-----------%
    if strcmp(model,'rl')==1      
        
        if strcmp(type,'reward')==1
            %rpe = r-v; % prediction error
            v = v;      % update the value
        elseif strcmp(type,'efficacy')==1
            %rpe = e-v; % prediction error
            v = v;      % update the value
        end
     
    end
    
    
    
    % calculate the SSE, but only if there is a probe response
    % on this trial
    if isnan(v_subj) == 0        
        SSE = SSE + (v-v_subj)^2;   
        SSE = SSE;
    else
        SSE = SSE;
    end
end
    sse = SSE;
end





