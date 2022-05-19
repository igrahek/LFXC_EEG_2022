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
    model_fitted = 'intercept_2LR';
    
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

save([sprintf('%s',type),'_',sprintf('%s',model_fitted),'.mat'],'results')

% % save the results
% if strcmp(type,'reward')==1
%     csvwrite([results_path,'learning_rates_reward.csv'],results.x)
% else
%     csvwrite([results_path,'learning_rates_efficacy.csv'],results.x)
% end
%% Model checks

% Calculate the subjective estimates predicted from the model

for s = 1:length(data) % for every subject
    for n = 1:data(s).N % for every trial
        if n==1 % if it's the first trial use the estimate of the initial value
            v = results.x(s,3);
        else 
            v = v;
        end
        
        if strcmp(type,'reward')==1
            r = data(s).IsRewarded(n);
            rpe = r-v; %calculate rpe
            predicted(s).rpe(n) = rpe; % save rpe
            if rpe > 0 % if the rpe is positive
                predicted(s).value(n) = v + results.x(s,1)* rpe; % save
                v = v + results.x(s,1)* rpe; % % update with the positive learning rate if the rpe is positive
            else % if the rpe is negative
                predicted(s).value(n) = v + results.x(s,2)* rpe; % save
                v = v + results.x(s,2)* rpe; % % update with the negative learning rate if the rpe is negative            
            end
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
            if rpe > 0 % if the rpe is positive
                predicted(s).value(n) = v + results.x(s,1)* rpe; % save
                v = v + results.x(s,1)* rpe; % % update with the positive learning rate if the rpe is positive
            else % if the rpe is negative
                predicted(s).value(n) = v + results.x(s,2)* rpe; % save
                v = v + results.x(s,2)* rpe; % % update with the negative learning rate if the rpe is negative            
            end
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

% save the results

% add the BIC to the learning rates results and save
bics = [BIC];
bics = bics';
results.x(:,6) = bics; 

if strcmp(type,'reward')==1
    csvwrite([results_path,'learning_rates_reward.csv'],results.x)
else
    csvwrite([results_path,'learning_rates_efficacy.csv'],results.x)
end
%% Plotting

% Plot the parameter estimates
figure('visible','off')
sgtitle({['Fit for ', type]
    ['Number of starts=' num2str(nstarts) '; method of fitting: SSE']
    ['R2 = ' num2str(round(mean(all.predicted.R2),2)) '; AIC = ' num2str(round(mean(all.predicted.AIC),2)) '; BIC = ' num2str(round(mean(all.predicted.BIC),2))]})

for i=1:length(param)
    subplot(1,length(param),i)
    histogram(results.x(:,i),100)
    xlim([-0.1 0.6])
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
        legend ('Location','NorthEastOutside','FontName', 'Helvetica','FontSize',12)
        hold on
        plot(1:max(allSubData.Interval),predicted(s).value,'-','LineWidth',1.5,'DisplayName','Predicted values')
        xlabel('Trial','FontSize',12,'FontName', 'Helvetica');
        ylabel('Estimate','FontSize',12,'FontName', 'Helvetica');
        hold on
        plot(1:max(allSubData.Interval),predicted(s).runavg,'-','LineWidth',1.5,'DisplayName','Running Average')
        %title(['sub=',sprintf('%1.0f',s),'; lr=', sprintf('%1.3f',results.x(s,1)),'; vinit= ', sprintf('%1.3f',results.x(s,2))])
        title(['sub=',sprintf('%1.0f',results.x(s,4)),'vinit=',sprintf('%1.0f',results.x(s,3)),'; lr_p_o_s=', sprintf('%1.3f',results.x(s,1)),'; lr_n_e_g= ', sprintf('%1.3f',results.x(s,2)), '; R2= ', sprintf('%1.3f',r2.R2(s)),'; AIC= ', sprintf('%1.3f',AIC(s)),'; BIC= ', sprintf('%1.3f',BIC(s))])
        ylim([-0.1 1.1])
    end
%         set(gcf,'units','centimeters','position',[0 0 50 20]) 
%         set(gcf,'units','centimeters','position',[0 0 22 16]) 
%         set(gca, 'FontName', 'Helvetica','fontsize',18)

%         saveas(gcf,[results_path,sprintf('%s',type),'\predicted_observed_',sprintf('%s',type),sprintf('%d',f),sprintf('%s',model_fitted),'.png'])
        h=gcf;
%         set(h,'PaperOrientation','landscape');
%         set(h,'PaperUnits','normalized');
%         set(h,'units','centimeters','position',[0 0 22 16]) 
        set( h,'PaperSize',[20 10], 'PaperPosition',[0 0 20 10]) 
        saveas(h,[results_path,sprintf('%s',type),'\predicted_observed_',sprintf('%s',type),sprintf('%d',f),sprintf('%s',model_fitted),'.png'])

end

% Plot the scatter of the positive and negative learning rates
figure('visible','off');
scatterhist(results.x(:,1),results.x(:,2), 'Direction','out','NBins',[10,10])
% xlim([-0.03 0.5])
% ylim([-0.03 0.5])
title({'Positive and negative learning rates for each subject' ''});
xlabel('Positive learning rate','FontSize',14);
ylabel('Negative learning rate','FontSize',14);
refline(1,0)
saveas(gcf,[results_path,sprintf('%s',type),'/Scatter_learning_rates_',sprintf('%s',type),'.png'])

% Plot the scatter of the log-positive and lognegative learning rates
figure('visible','off');
scatterhist(log(results.x(:,1)),log(results.x(:,2)), 'Direction','out','NBins',[10,10])
% xlim([-0.03 0.5])
% ylim([-0.03 0.5])
title({'Positive and negative learning rates for each subject' ''});
xlabel('log of the positive learning rate','FontSize',14);
ylabel('log of the negative learning rate','FontSize',14);
refline(1,0)
saveas(gcf,[results_path,sprintf('%s',type),'/Scatter_learning_rates_log_',sprintf('%s',type),'.png'])

% % Plot the histogram of the difference between positive and negative learning rates
% figure('visible','off');
% LR_diff = results.x(:,1) - results.x(:,2);
% histogram(LR_diff,20)
% % xlim([-0.03 0.5])
% % ylim([-0.03 0.5])
% title({'Positive minus negative learning rate' ''});
% xlabel('Positive minus negative learning rate','FontSize',14);
% ylabel('Frequency','FontSize',14);
% saveas(gcf,[results_path,sprintf('%s',type),'/Histogram_learning_rates_',sprintf('%s',type),'.png'])
% Plot the actual feedbacks
% for f = 1:length(data)/10 %number of plots
%     figure;
%     sgtitle(['Fit for ', type, ['; nstarts=' num2str(nstarts) '; method of fitting: SSE']])
%     for s = f*10-9:f*10 % for every subject
%         subplot(5,2,s-((f-1)*10))
%         plot(1:max(allSubData.Interval),predicted(s).actual_value,'o','DisplayName','Subjective estimate', 'MarkerFaceColor', 'b','MarkerSize',4)
%         legend ('Location','NorthEastOutside')
%         hold on
%         plot(1:max(allSubData.Interval),predicted(s).value,'-','LineWidth',1.5,'DisplayName','Predicted values')
%         xlabel('Trial','FontSize',10);
%         ylabel('Estimate','FontSize',10);
%         hold on
%         plot(1:max(allSubData.Interval),predicted(s).feedback,'o','DisplayName','Feedback', 'MarkerFaceColor', 'k','MarkerEdgeColor','k','MarkerSize',2)
%         %title(['sub=',sprintf('%1.0f',s),'; lr=', sprintf('%1.3f',results.x(s,1)),'; vinit= ', sprintf('%1.3f',results.x(s,2))])
%         title(['sub=',sprintf('%1.0f',results.x(s,4)),'vinit=',sprintf('%1.0f',results.x(s,3)),'; lr_p_o_s=', sprintf('%1.3f',results.x(s,1)),'; lr_n_e_g= ', sprintf('%1.3f',results.x(s,2)), '; R2= ', sprintf('%1.3f',r2.R2(s)),'; AIC= ', sprintf('%1.3f',AIC(s)),'; BIC= ', sprintf('%1.3f',BIC(s))])
%         ylim([-0.1 1.1])
%     end
%         set(gcf,'units','centimeters','position',[0 0 50 20]) 
%         saveas(gcf,[results_path,sprintf('%s',type),'/predicted_observed_',sprintf('%s',type),sprintf('%d',f),'_plotted_trials.png'])
% end

end