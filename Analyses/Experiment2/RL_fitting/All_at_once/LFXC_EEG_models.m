% Demo of hierarchical MFIT applied to a simple reinforcement learning model.
% Two models are compared:
%   Model 1: single learning rate (and inverse temperature)
%   Model 2: separate learning rates for positive and negative prediction errors
% Ground-truth data are generated from Model 1.

% ---------- Import the data ----------%

% Load the data
load('lfxc_eeg_40subs.mat')

% Prepare the data
SubIDs = unique(allSubData.SubID);

for s = 1: length(SubIDs)
    data(s).Subject = allSubData.SubID(allSubData.SubID == SubIDs(s));
    data(s).Subject = data(s).Subject(1);
    data(s).N = length(allSubData.RT(allSubData.SubID == SubIDs(s))); 
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

% ------------ fit models --------------------%

mu = 0.5; sigma = 0.3;   % parameters of gaussian prior
param(1).name = 'intercept';
param(1).logpdf = @(x) sum(log(normpdf(x,mu,sigma)));
param(1).lb = 0;
param(1).ub = 1;

param(2).name = 'variance';
param(2).logpdf = @(x) 0; % uniform prior
param(2).lb = 0;
param(2).ub = 0.5;

a = 1.2; b = 1.2;   % parameters of beta prior
param(3).name = 'positive learning rate';
param(3).logpdf = @(x) sum(log(betapdf(x,a,b)));
param(3).lb = 0;
param(3).ub = 1;

a = 1.2; b = 1.2;   % parameters of beta prior
param(4).name = 'negative learning rate';
param(4).logpdf = @(x) sum(log(betapdf(x,a,b)));
param(4).lb = 0;
param(4).ub = 1;

% select which data to fit
type = 'efficacy';

% run optimization
nstarts = 5;    % number of random parameter initializations
disp('... Fitting model 1');
results(1) = mfit_optimize_hierarchical(@rllik1,param(1:2),data,type,nstarts);
disp('... Fitting model 2');
results(2) = mfit_optimize_hierarchical(@rllik2,param(1:3),data,type,nstarts);
disp('... Fitting model 3');
results(3) = mfit_optimize_hierarchical(@rllik3,param(1:4),data,type,nstarts);

% compute predictive probability for the two models on test data
logp(:,1) = mfit_predict(testdata,results(1));
logp(:,2) = mfit_predict(testdata,results(2));
logp(:,3) = mfit_predict(testdata,results(3));
%-------- plot results -----------%

r = corr(results(1).x(:),x(:));
disp(['Correlation between true and estimated parameters: ',num2str(r)]);
figure;
plot(results(1).x(:),x(:),'+k','MarkerSize',12,'LineWidth',4);
h = lsline; set(h,'LineWidth',4);
set(gca,'FontSize',25);
xlabel('Estimate','FontSize',25);
ylabel('Ground truth','FontSize',25);

bms_results = mfit_bms(results);
figure;
bar(bms_results.xp); colormap bone;
set(gca,'XTickLabel',{'Model 1' 'Model 2' 'Model 3'},'FontSize',25,'YLim',[0 1]);
ylabel('Exceedance probability','FontSize',25);
title('Bayesian model comparison','FontSize',25);

figure;
d = logp(:,3)-logp(:,2);
m = mean(d);
se = std(d)/sqrt(S);
errorbar(m,se,'ok','MarkerFaceColor','k','MarkerSize',12,'LineWidth',4);
set(gca,'YLim',[-1 max(d)+1],'XLim',[0.5 1.5],'XTick',1,'XTickLabel',{'Model 1 vs. Model 2'},'FontSize',25);
ylabel('Relative log predictive prob.','FontSize',25);
hold on; plot([0.5 1.5],[0 0],'--r','LineWidth',3); % red line shows chance performance
title('Cross-validation','FontSize',25);
