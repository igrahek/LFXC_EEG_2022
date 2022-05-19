% Code for fitting a single learning rate for every subject for the
% LFXC_EEG dataset.

% Code written by: Ivan Grahek (grahek.i@gmail.com) initially based on the
% mfit code by Sam Gershman (https://github.com/sjgershm/mfit)

%% Initial stuff
clear
% Which version to load
version = 'TSS_LFXC_MEGA';

% probe response to analyze (reward vs. efficacy)
type_array =  {'reward', 'efficacy'};

% Specify the root directory
root = '..\..\..\..\Analyses\Experiment2\Stats\brms\';

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
    model_fitted = 'intercept';
    
    % number of random parameter initializations
    nstarts = 20;
    
% path for saving results
results_path = '.\results\';
results_path_AIC = '..\Model_comparison\data\';

%% Importing data
% Get the list of task versions
files = dir('..\..\..\..\Analyses\Experiment2\Stats\brms\');

% Load the data
if ~any(strcmp(char({files.name}),version))

        opts = detectImportOptions(fullfile(root, ['RL_fit_data_',version,'.csv']));
        opts.VariableTypes = {'double','double','double','double','double','double','double','double','double','double','double','double','double'};
        allSubData = readtable(fullfile(root, ['RL_fit_data_',version,'.csv']),opts);
 
else
    fprintf('The dataset for this version does not exist')
end


% Make a list of all subjects
SubIDs = unique(allSubData.SubID);

for s = 1: length(SubIDs)
    %data(s).Subject = s;
    data(s).Subject = allSubData.SubID(allSubData.SubID == SubIDs(s));
    data(s).Subject = data(s).Subject(1);
    data(s).N = length(allSubData.RT(allSubData.SubID == SubIDs(s))); %length
    
    data(s).RT = allSubData.RT(allSubData.SubID == SubIDs(s));
    data(s).Efficacy_T1 = allSubData.Efficacy_T1(allSubData.SubID == SubIDs(s));
    data(s).IsRewarded_T1 = allSubData.IsRewarded_T1(allSubData.SubID == SubIDs(s));
    data(s).EffLvl = allSubData.EffLvl(allSubData.SubID == SubIDs(s));
    data(s).IsRewarded = allSubData.IsRewarded(allSubData.SubID == SubIDs(s));
    data(s).ReportedEfficacyLin = allSubData.ReportedEfficacyLin(allSubData.SubID == SubIDs(s));
    data(s).ReportedRewRateLin = allSubData.ReportedRewRateLin(allSubData.SubID == SubIDs(s));
    data(s).EfficacyProbeResp = allSubData.EfficacyProbeResp(allSubData.SubID == SubIDs(s));
    data(s).RewRateProbeResp = allSubData.RewRateProbeResp(allSubData.SubID == SubIDs(s));
    data(s).runAvgEfficacy = allSubData.runAvgEfficacy(allSubData.SubID == SubIDs(s));
    data(s).runAvgRewRate = allSubData.runAvgRewRate(allSubData.SubID == SubIDs(s));
end

% % % Min-max normalization of the efficacy response
% for s = 1:length(SubIDs)
%     data(s).EfficacyProbeResp = (data(s).EfficacyProbeResp - min(data(s).EfficacyProbeResp)) / (max(data(s).EfficacyProbeResp) - min(data(s).EfficacyProbeResp));
% end
% 
% % Min-max normalization of the reward response
% for s = 1:length(SubIDs)
%     data(s).RewRateProbeResp = (data(s).RewRateProbeResp - min(data(s).RewRateProbeResp)) / (max(data(s).RewRateProbeResp) - min(data(s).RewRateProbeResp));
% end

%% Fitting

mu = 0.5; sigma = 0.3;   % parameters of gaussian prior
param(1).name = 'intercept';
param(1).logpdf = @(x) sum(log(normpdf(x,mu,sigma)));
param(1).lb = 0;
param(1).ub = 1;

param(2).name = 'variance';
param(2).logpdf = @(x) 0; % uniform prior
param(2).lb = 0;
param(2).ub = 0.5;


% run optimization
disp('... Fitting model 1');

results = mfit_optimize_hierarchical(@rllik,param,data,type,nstarts);

% 
% % fill in missing options
% %     if nargin < 4 || isempty(nstarts); nstarts = 5; end
%     K = length(param);
%     results.K = K;
%     
%     % save info to results structure
%     results.param = param;
% %     results.likfun = cost;
%     
%     % extract lower and upper bounds
%     lb = [param.lb];
%     ub = [param.ub];
%     
%     options = optimset('Display','off','MaxFunEvals',1e6,'MaxIter',1e6);
%     warning off all
%     
%     for s = 1:length(data)
%         disp(['Subject ',num2str(s)]);
%         
%         % construct posterior function
%         f = @(x) cost(x,data(s),type,model);
% 
%         for i = 1:nstarts
%             x0 = zeros(1,K);
%             for k = 1:K
%                 x0(k) = unifrnd(param(k).lb,param(k).ub);
%             end
%             [x,sse] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
%             if i == 1 || results.sse(s) > sse
%                 results.sse(s) = sse;
%                 results.x(s,:) = x;
%             end
%         end
%         
%      
%     end

% add the subject number to the results
subs = [data.Subject];
subs = subs';
results.x(:,3) = subs; 

save([sprintf('%s',type),'_',sprintf('%s',model_fitted),'.mat'],'results')

% save the results
if strcmp(type,'reward')==1
    csvwrite([results_path,'learning_rates_reward.csv'],results.x)
else
    csvwrite([results_path,'learning_rates_efficacy.csv'],results.x)
end
%% Model checks

% Calculate the subjective estimates predicted from the model

for s = 1:length(data) % for every subject
    for n = 1:data(s).N % for every trial

            v = results.x(s,1);

        
        if strcmp(type,'reward')==1
            r = data(s).IsRewarded(n);
            rpe = r-v; %calculate rpe
            predicted(s).rpe(n) = rpe; % save rpe
            predicted(s).value(n) = v ; % save
            predicted(s).actual_value(n) = data(s).RewRateProbeResp(n); % save the actual response
            predicted(s).actual_value_interpolated(n) = data(s).ReportedRewRateLin(n); % save the actual response, but interpolated
            predicted(s).runavg(n) = data(s).runAvgRewRate(n); % save the running average (5)
            predicted(s).feedback(n) = data(s).IsRewarded(n); % save the feedback
            predicted(s).residual(n) =  predicted(s).actual_value(n) - predicted(s).value(n); % residual
            predicted(s).residual_from_mean(n) =  predicted(s).actual_value(n) - nanmean(predicted(s).actual_value); % residual
                        


        else
            r = data(s).EffLvl(n);
            rpe = r-v; %calculate rpe
            predicted(s).rpe(n) = rpe; % save rpe
            predicted(s).value(n) = v ; % save
            predicted(s).actual_value(n) = data(s).EfficacyProbeResp(n); % save the actual response
            predicted(s).actual_value_interpolated(n) = data(s).ReportedEfficacyLin(n); % save the actual response, but interpolated
            predicted(s).runavg(n) = data(s).runAvgEfficacy(n); % save the running average (5)  
            predicted(s).feedback(n) = data(s).EffLvl(n); % save the feedback
            predicted(s).residual(n) =  predicted(s).actual_value(n) - predicted(s).value(n); % residual
            predicted(s).residual_from_mean(n) =  predicted(s).actual_value(n) - nanmean(predicted(s).actual_value); % residual
        end
        
    end
    
% transpose
predicted(s).value = predicted(s).value';
predicted(s).rpe = predicted(s).rpe';
predicted(s).actual_value = predicted(s).actual_value';
predicted(s).actual_value_interpolated = predicted(s).actual_value_interpolated';
predicted(s).feedback(n) = predicted(s).feedback(n)'; % save the feedback
predicted(s).residual(n) = predicted(s).residual(n)';
predicted(s).residual_from_mean(n) = predicted(s).residual_from_mean(n)';
end


% R2 for the whole sample

% Merge the whole structure
all.predicted.actual_value = predicted(1).actual_value;
all.predicted.value = predicted(1).value;

for i = 2:length(predicted)
   all.predicted.actual_value = [all.predicted.actual_value,predicted(i).actual_value]; 
   all.predicted.value = [all.predicted.value,predicted(i).value]; 
end

% Calculate the R2
   all.predicted.SSEres = nansum((all.predicted.actual_value - all.predicted.value).^2); 
   all.predicted.SSEtot = nansum((all.predicted.actual_value - nanmean(all.predicted.actual_value)).^2);  
   all.predicted.R2 = 1 - all.predicted.SSEres/all.predicted.SSEtot;

   
% R2 per subject
for s = 1:length(data) % for every subject
 r2.SSEtot(s) = nansum(predicted(s).residual_from_mean.^2);   
 r2.SSEres(s) = nansum(predicted(s).residual.^2);   
 r2.R2(s) = 1 - r2.SSEres(s)/r2.SSEtot(s);
end

% Akaike's information criterion (Akaike, 1969)
all.predicted.AIC = data(s).N * length(data) * log(sum(all.predicted.SSEres) / length(data)) + 2 * length(param);

% Akaike's information criterion (Akaike, 1969) per subject
for s = 1:length(data) % for every subject
    AIC(s) = data(s).N * log(r2.SSEres(s) / data(s).N) + 2 * length(param);
end

% Schwarz's Bayesian criterion (or BIC) (Schwarz, 1978)
all.predicted.BIC = data(s).N * length(data) * log(sum(all.predicted.SSEres) / length(data)) + length(param) * log(data(s).N * length(data));

% Schwarz's Bayesian criterion (or BIC) (Schwarz, 1978) per subject
for s = 1:length(data) % for every subject
    BIC(s) = data(s).N * log(r2.SSEres(s) / data(s).N) + length(param) * log(data(s).N);
end

% Save the AIC and BIC
AIC_BIC.AIC = AIC';
AIC_BIC.AIC(:,2) = subs;
AIC_BIC.BIC = BIC';
AIC_BIC.BIC(:,2) = subs;

save([results_path_AIC,sprintf('%s',type),'_',sprintf('%s',model_fitted),'_AIC_BIC'],'AIC_BIC')
save([results_path_AIC,sprintf('%s',type),'_',sprintf('%s',model_fitted),'_fitting_results'],'all')

%% Plotting

% Plot the parameter estimates
figure('visible','off')
sgtitle({['Fit for ', type]
    ['Number of starts=' num2str(nstarts) '; method of fitting: SSE']
    ['R2 = ' num2str(round(mean(all.predicted.R2),2)) '; AIC = ' num2str(round(mean(all.predicted.AIC),2)) '; BIC = ' num2str(round(mean(all.predicted.BIC),2))]})

for i=1:length(param)
    subplot(1,length(param),i)
    histogram(results.x(:,i),100)
    xlim([-0 1])
    title([param(i).name])
    set(gcf,'units','centimeters','position',[0 0 40 20])
end
saveas(gcf,[results_path,sprintf('%s',type),'/estimates_',sprintf('%s',type),sprintf('%s',model_fitted),'.png'])

% Plot the fitted vs. actual values

% Dots for the actual responses and a line for the model-based estimates
for f = 1:length(data)/10 %number of plots
    figure('visible','off');
    sgtitle(['Fit for ', type, ['; nstarts=' num2str(nstarts) '; method of fitting: SSE']])
    for s = f*10-9:f*10 % for every subject
        subplot(5,2,s-((f-1)*10))
        plot(1:max(allSubData.Interval),predicted(s).actual_value,'o','DisplayName','Subjective estimate', 'MarkerFaceColor', 'b','MarkerSize',4)
        legend ('Location','NorthEastOutside')
        hold on
        plot(1:max(allSubData.Interval),predicted(s).value,'-','LineWidth',1.5,'DisplayName','Predicted values')
        xlabel('Trial','FontSize',10);
        ylabel('Estimate','FontSize',10);
        hold on
        plot(1:max(allSubData.Interval),predicted(s).runavg,'-','LineWidth',1.5,'DisplayName','Running Avg (5)')
        %title(['sub=',sprintf('%1.0f',s),'; lr=', sprintf('%1.3f',results.x(s,1)),'; vinit= ', sprintf('%1.3f',results.x(s,2))])
        title(['sub=',sprintf('%1.0f',results.x(s,2)),'; intercept=', sprintf('%1.3f',results.x(s,1)), '; R2= ', sprintf('%1.3f',r2.R2(s)),'; AIC= ', sprintf('%1.3f',AIC(s)),'; BIC= ', sprintf('%1.3f',BIC(s))])
        ylim([-0.1 1.1])
    end
        set(gcf,'units','centimeters','position',[0 0 50 20]) 
        saveas(gcf,[results_path,sprintf('%s',type),'/predicted_observed_',sprintf('%s',type),sprintf('%d',f),sprintf('%s',model_fitted),'.png'])
end




end