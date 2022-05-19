% This script generates data for the parameter recovery study. 

%% Import the efficacy drift
clear all
load('lfxc_eeg_40subs.mat')
EfficacyDrift1 = allSubData.EffLvl(allSubData.SubID == 1042);
EfficacyDrift2 = allSubData.EffLvl(allSubData.SubID == 1043);

RewardDrift = allSubData.IsRewarded(allSubData.SubID == 1042);

Noisevector = [0.15];

for Nnoise = 1:length(Noisevector)
    noise = Noisevector(Nnoise);

%% Set the simulation parameters
% LRsPositive = [0.01:0.02:0.4 0.01:0.02:0.4];
% LRsNegative = [0.01:0.02:0.4 0.01:0.02:0.4];

min = 0.001;
max = 0.4;
N = 200;
LRPositive = (max-min).*rand(N,1) + min;
LRNegative = (max-min).*rand(N,1) + min;

intercept = 0.5;

%% Simulate the data

% Drift 1
for s = 1:length(LRPositive)/2 % for every subject
    for t = 1:length(EfficacyDrift1) % for every trial
        if t==1 % if it's the first trial use the estimate of the initial value
            v_R = intercept;
            v_E = intercept;
        else
            v_R = v_R;
            v_E = v_E;
        end
        
        e = EfficacyDrift1(t);
        r = RewardDrift(t);
        
        data(s).EffLvl(t) = e;
        data(s).IsRewarded(t) = r;
        pe = e-v_E; %calculate efficacy pe
        rpe = r-v_R; %calculate rpe

        % Efficacy
        if pe >= 0
            v_E = v_E + LRPositive(s)*pe; % update with the positive learning rate if the pe is positive
        else
            v_E = v_E + LRNegative(s)*pe; % update with the negative learning rate if the pe is negative
        end
        
        % Save on every 6th trial
        if ismember(t,1:6:288)==1
            data(s).EfficacyProbeResp(t) = v_E + normrnd(0,noise); % save the subjective estimate
        else
            data(s).EfficacyProbeResp(t) = NaN;
        end
        
        % Reward
        if rpe >= 0
            v_R = v_R + LRPositive(s)*rpe; % update with the positive learning rate if the pe is positive
        else
            v_R = v_R + LRNegative(s)*rpe; % update with the negative learning rate if the pe is negative
        end
        
        % Save on every 6th trial
        if ismember(t,1:6:288)==1
            data(s).RewRateProbeResp(t) = v_R + normrnd(0,noise); % save the subjective estimate
        else
            data(s).RewRateProbeResp(t) = NaN;
        end
    end
   
% save
data(s).EfficacyProbeResp = data(s).EfficacyProbeResp';
data(s).EffLvl = data(s).EffLvl';

data(s).RewRateProbeResp = data(s).RewRateProbeResp';
data(s).IsRewarded = data(s).IsRewarded';

data(s).PosELR = LRPositive(s);
data(s).NegELR = LRNegative(s);
data(s).PosRLR = LRPositive(s);
data(s).NegRLR = LRNegative(s);
data(s).N = length(EfficacyDrift1);
data(s).Subject = s;
end


% Drift 2
for s = (length(LRPositive)/2+1):length(LRPositive) % for every subject
    for t = 1:length(EfficacyDrift2) % for every trial
        if t==1 % if it's the first trial use the estimate of the initial value
            v_R = intercept;
            v_E = intercept;
        else
            v_R = v_R;
            v_E = v_E;
        end
        
        e = EfficacyDrift2(t);
        r = RewardDrift(t);
        
        data(s).EffLvl(t) = e;
        data(s).IsRewarded(t) = r;
        pe = e-v_E; %calculate efficacy pe
        rpe = r-v_R; %calculate rpe

        % Efficacy
        if pe >= 0
            v_E = v_E + LRPositive(s)*pe; % update with the positive learning rate if the pe is positive
        else
            v_E = v_E + LRNegative(s)*pe; % update with the negative learning rate if the pe is negative
        end
        
        % Save on every 6th trial
        if ismember(t,1:6:288)==1
            data(s).EfficacyProbeResp(t) = v_E + normrnd(0,noise); % save the subjective estimate
        else
            data(s).EfficacyProbeResp(t) = NaN;
        end
        
        % Reward
        if rpe >= 0
            v_R = v_R + LRPositive(s)*rpe; % update with the positive learning rate if the pe is positive
        else
            v_R = v_R + LRNegative(s)*rpe; % update with the negative learning rate if the pe is negative
        end
        
        % Save on every 6th trial
        if ismember(t,1:6:288)==1
            data(s).RewRateProbeResp(t) = v_R + normrnd(0,noise); % save the subjective estimate
        else
            data(s).RewRateProbeResp(t) = NaN;
        end
    end
   
% save
data(s).EfficacyProbeResp = data(s).EfficacyProbeResp';
data(s).EffLvl = data(s).EffLvl';

data(s).RewRateProbeResp = data(s).RewRateProbeResp';
data(s).IsRewarded = data(s).IsRewarded';

data(s).PosELR = LRPositive(s);
data(s).NegELR = LRNegative(s);
data(s).PosRLR = LRPositive(s);
data(s).NegRLR = LRNegative(s);
data(s).N = length(EfficacyDrift1);
data(s).Subject = s;
end

%% Fit the model
% probe response to analyze (reward vs. efficacy)
type_array =  {'reward', 'efficacy'};

for t = 1:length(type_array)
    
        % Clear the previous results
    clear results
    clear predicted
    clear all.predicted
    
    % Which data to fit (reward rate or efficacy) ('reward' or 'efficacy')
    type = char(type_array(t));
    
    % Which model to fit ('rl')
    model = 'rl';
    
    % Which of the models is fit (fixed, intercept, one, or two LRs)
    model_fitted = 'intercept_2LR';
    
    % number of random parameter initializations
    nstarts=10;
    
% path for saving results
results_path = 'results\';
results_path_AIC = '..\Model_comparison\data\';
    
%% Fitting

a = 1.2; b = 1.2;   % parameters of beta prior
param(1).name = 'positive learning rate';
param(1).logpdf = @(x) sum(log(betapdf(x,a,b)));
param(1).lb = 0;
param(1).ub = 1;

a = 1.2; b = 1.2;   % parameters of beta prior
param(2).name = 'negative learning rate';
param(2).logpdf = @(x) sum(log(betapdf(x,a,b)));
param(2).lb = 0;
param(2).ub = 1;

mu = 0.5; sigma = 0.3;   % parameters of gaussian prior
param(3).name = 'intercept';
param(3).logpdf = @(x) sum(log(normpdf(x,mu,sigma)));
param(3).lb = 0;
param(3).ub = 1;

param(4).name = 'variance';
param(4).logpdf = @(x) 0; % uniform prior
param(4).lb = 0;
param(4).ub = 0.5;


% run optimization
disp('... Fitting model 1');

results = mfit_optimize_hierarchical(@rllik,param(1:length(param)),data,type,nstarts);

% add the subject number to the results
subs = [data.Subject];
subs = subs';
results.x(:,5) = subs; 

% add the ground truth learning rates
results.x(:,6) = LRPositive';
results.x(:,7) = LRNegative';
results.x(:,8) = repmat(noise,1,length(LRPositive))';

% % save the results
% save([sprintf('%s',type),'_',sprintf('%s',model_fitted),'.mat'],'results')

%% Save

if strcmp(type,'reward')==1
    csvwrite([results_path,['learning_rates_reward_noise_',num2str(noise),'.csv']],results.x)
else
    csvwrite([results_path,['learning_rates_efficacy_noise_',num2str(noise),'.csv']],results.x)
end

%% Plotting

% Plot the parameter estimates
figure('visible','off')
sgtitle({['Fit for ', type]})
%     ['Number of starts=' num2str(nstarts) '; method of fitting: SSE']
%     ['R2 = ' num2str(round(mean(all.predicted.R2),2)) '; AIC = ' num2str(round(mean(all.predicted.AIC),2)) '; BIC = ' num2str(round(mean(all.predicted.BIC),2))]

for i=1:length(param)
    subplot(1,length(param),i)
    histogram(results.x(:,i),100)
    xlim([-0.1 0.6])
    title([param(i).name])
    set(gcf,'units','centimeters','position',[0 0 40 20])
end
saveas(gcf,[results_path,sprintf('%s',type),'/estimates_',sprintf('%s',type),sprintf('%s',model_fitted),num2str(noise),'.png'])

% Plot the correlations

% Positive LR 
figure('visible','off')
sgtitle({['Positive LR for ', type,'; Noise = ',num2str(noise)]})
scatter(results.x(:,1),results.x(:,6))
xlabel('Estimated') 
ylabel('Ground truth') 
refline
saveas(gcf,[results_path,sprintf('%s',type),'/Pos_LR_',sprintf('%s',type),sprintf('%s',model_fitted),num2str(noise),'.png'])

% Negative LR 
figure('visible','off')
sgtitle({['Negative LR for ', type,'; Noise = ',num2str(noise)]})
scatter(results.x(:,2),results.x(:,7))
xlabel('Estimated') 
ylabel('Ground truth') 
refline
saveas(gcf,[results_path,sprintf('%s',type),'/Neg_LR_',sprintf('%s',type),sprintf('%s',model_fitted),num2str(noise),'.png'])


end

end