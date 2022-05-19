% Code for plotting the AICs and BICs for all the fitted learning models 
% for the LFXC_EEG project 

% Code written by: Ivan Grahek (grahek.i@gmail.com) initially based on the
% mfit code by Sam Gershman (https://github.com/sjgershm/mfit)

%% Load all the data
clear 

% Set the directory
cd '..\data'

% Per subject

% Efficacy
AIC_BIC.int_efficacy = load('efficacy_intercept_AIC_BIC.mat');
AIC_BIC.int_1LR_efficacy = load('efficacy_intercept_1LR_AIC_BIC.mat');
AIC_BIC.int_2LR_efficacy = load('efficacy_intercept_2LR_AIC_BIC.mat');

% Reward
AIC_BIC.int_reward = load('reward_intercept_AIC_BIC.mat');
AIC_BIC.int_1LR_reward = load('reward_intercept_1LR_AIC_BIC.mat');
AIC_BIC.int_2LR_reward = load('reward_intercept_2LR_AIC_BIC.mat');

% For the whole sample

% Efficacy
int_efficacy = load('efficacy_intercept_fitting_results.mat');
int_1LR_efficacy = load('efficacy_intercept_1LR_fitting_results.mat');
int_2LR_efficacy = load('efficacy_intercept_2LR_fitting_results.mat');

% Reward
int_reward = load('reward_intercept_fitting_results.mat');
int_1LR_reward = load('reward_intercept_1LR_fitting_results.mat');
int_2LR_reward = load('reward_intercept_2LR_fitting_results.mat');

% The path for saving the plots
results_path_model_comparison = '..\results\';

results_path_efficacy = '..\results\efficacy\';

results_path_reward = '..\results\reward\';

% Calculate the number of participants
Nsub = length(int_2LR_efficacy.all.predicted.actual_value(1,:));

%% Efficacy AIC plot per subject

figure
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 15 3];
sgtitle('AIC per subject (efficacy)')


% Subjects 1-20
y = [AIC_BIC.int_efficacy.AIC_BIC.AIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_1LR_efficacy.AIC_BIC.AIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_2LR_efficacy.AIC_BIC.AIC(1:(ceil(Nsub/2)),1)];
x = (AIC_BIC.int_efficacy.AIC_BIC.AIC(1:(ceil(Nsub/2)),2));
x = categorical(x);
subplot(2,1,1)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

% Subjects 21-40
y = [AIC_BIC.int_efficacy.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_1LR_efficacy.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_2LR_efficacy.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1)];
x = (AIC_BIC.int_efficacy.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,2));
x = categorical(x);
subplot(2,1,2)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

set(gcf,'units','centimeters','position',[0 0 50 20]) 
% Save
saveas(gcf,[results_path_efficacy,'AIC_per_sub_efficacy','.png'])

%% Reward AIC plot per subject

figure
sgtitle('AIC per subject (reward)')

% Subjects 1-20
y = [AIC_BIC.int_reward.AIC_BIC.AIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_1LR_reward.AIC_BIC.AIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_2LR_reward.AIC_BIC.AIC(1:(ceil(Nsub/2)),1)];
x = (AIC_BIC.int_reward.AIC_BIC.AIC(1:(ceil(Nsub/2)),2));
x = categorical(x);
subplot(2,1,1)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

% Subjects 21-40
y = [AIC_BIC.int_reward.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_1LR_reward.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_2LR_reward.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,1)];
x = (AIC_BIC.int_reward.AIC_BIC.AIC((ceil(Nsub/2)+1):Nsub,2));
x = categorical(x);
subplot(2,1,2)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

set(gcf,'units','centimeters','position',[0 0 50 20]) 
%Save
saveas(gcf,[results_path_reward,'AIC_per_sub_reward','.png'])

%% Efficacy BIC plot per subject

figure
sgtitle('BIC per subject (efficacy)')

% Subjects 1-20
y = [AIC_BIC.int_efficacy.AIC_BIC.BIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_1LR_efficacy.AIC_BIC.BIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC(1:(ceil(Nsub/2)),1)];
x = (AIC_BIC.int_efficacy.AIC_BIC.BIC(1:(ceil(Nsub/2)),2));
x = categorical(x);
subplot(2,1,1)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

% Subjects 21-40
y = [AIC_BIC.int_efficacy.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_1LR_efficacy.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1)];
x = (AIC_BIC.int_efficacy.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,2));
x = categorical(x);
subplot(2,1,2)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])
set(gcf,'units','centimeters','position',[0 0 50 20]) 

% Save
saveas(gcf,[results_path_efficacy,'BIC_per_sub_efficacy','.png'])

%% Reward BIC plot per subject

figure
sgtitle('BIC per subject (reward)')

% Subjects 1-20
y = [AIC_BIC.int_reward.AIC_BIC.BIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_1LR_reward.AIC_BIC.BIC(1:(ceil(Nsub/2)),1) AIC_BIC.int_2LR_reward.AIC_BIC.BIC(1:(ceil(Nsub/2)),1)];
x = (AIC_BIC.int_reward.AIC_BIC.BIC(1:(ceil(Nsub/2)),2));
x = categorical(x);
subplot(2,1,1)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])

% Subjects 21-40
y = [AIC_BIC.int_reward.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_1LR_reward.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1) AIC_BIC.int_2LR_reward.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,1)];
x = (AIC_BIC.int_reward.AIC_BIC.BIC((ceil(Nsub/2)+1):Nsub,2));
x = categorical(x);
subplot(2,1,2)
bar(x,y)
set(gca, 'YDir','reverse')
legend ('Location','northeastoutside')
legend('Intercept', '1 LR', '2 LRs')
ylim([-2000 -500])
set(gcf,'units','centimeters','position',[0 0 50 20]) 

%Save
saveas(gcf,[results_path_reward,'BIC_per_sub_reward','.png'])

%% Efficacy plot for the whole sample
% BIC for efficacy 
figure
names = categorical({'Intercept','1 Learning rate','2 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate','2 Learning rate'});
BICs = [AIC_BIC.int_efficacy.AIC_BIC.BIC(:,1) AIC_BIC.int_1LR_efficacy.AIC_BIC.BIC(:,1) AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC(:,1)];
BICs_means = mean(BICs);
bar(names,BICs_means,0.2,'FaceColor','white','LineWidth', 1.4)
hold
plot(names,BICs,'o','color','black','LineWidth', 1.4)
yline(0,'LineWidth', 1.4)
% xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('BIC')
ylim([-2600 -1100])
yticks([-2400 -1900 -1400])
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca,'linewidth',2)
title('Efficacy')
% saveas(gcf,[results_path_model_comparison,'BIC_efficacy','.png'])
% figure
% sgtitle('AIC and BIC for the whole sample (efficacy)')
% % AIC
% AICs = [int_efficacy.all.predicted.AIC int_1LR_efficacy.all.predicted.AIC int_2LR_efficacy.all.predicted.AIC];
% names = categorical({'Intercept', '1 LR', '2 LRs'});
% names = reordercats(names,{'Intercept', '1 LR', '2 LRs'});
% % Plot
% subplot(1,2,1)
% bar(names,AICs)
% title('AIC')
% 
% % BIC
% BICs = [int_efficacy.all.predicted.BIC int_1LR_efficacy.all.predicted.BIC int_2LR_efficacy.all.predicted.BIC];
% names = categorical({'Intercept', '1 LR', '2 LRs'});
% names = reordercats(names,{'Intercept', '1 LR', '2 LRs'});
% % Plot
% subplot(1,2,2)
% bar(names,BICs)
% title('BIC')
% set(gcf,'units','centimeters','position',[0 0 30 10]) 
% 
% %Save
% saveas(gcf,[results_path_efficacy,'AIC_BIC_whole_sample_efficacy','.png'])


% BIC for efficacy 
figure
names = categorical({'Intercept','1 Learning rate','2 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate','2 Learning rate'});
BICs = [AIC_BIC.int_efficacy.AIC_BIC.BIC(:,1) AIC_BIC.int_1LR_efficacy.AIC_BIC.BIC(:,1) AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC(:,1)];
BICs_means = mean(BICs);
bar(names,BICs_means,0.2,'FaceColor','white','LineWidth', 1.4)
hold
plot(names,BICs,'o','color','black','LineWidth', 1.4)
yline(0,'LineWidth', 1.4)
% xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('BIC')
ylim([-2600 -1100])
yticks([-2400 -1900 -1400])
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca,'linewidth',2)
title('Efficacy')
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca,'linewidth',2)
% saveas(gcf,[results_path_model_comparison,'Delta_BIC_reward','.png'])
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure5_BNew.pdf')
% 
%% Reward plot for the whole sample

% BIC for reward 
figure
names = categorical({'Intercept','1 Learning rate','2 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate','2 Learning rate'});
BICs = [AIC_BIC.int_reward.AIC_BIC.BIC(:,1) AIC_BIC.int_1LR_reward.AIC_BIC.BIC(:,1) AIC_BIC.int_2LR_reward.AIC_BIC.BIC(:,1)];
BICs_means = mean(BICs);
bar(names,BICs_means,0.2,'FaceColor','white','LineWidth', 1.4)
hold
plot(names,BICs,'o','color','black','LineWidth', 1.4)
yline(0,'LineWidth', 1.4)
% xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('BIC')
ylim([-3000 -1000])
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca,'linewidth',2)
title('Reward')
% saveas(gcf,[results_path_model_comparison,'BIC_reward','.png'])

% 
% figure
% sgtitle('AIC and BIC for the whole sample (reward)')
% % AIC
% AICs = [int_reward.all.predicted.AIC int_1LR_reward.all.predicted.AIC int_2LR_reward.all.predicted.AIC];
% names = categorical({'Intercept', '1 LR', '2 LRs'});
% names = reordercats(names,{'Intercept', '1 LR', '2 LRs'});
% % Plot
% subplot(1,2,1)
% bar(names,AICs)
% title('AIC')
% 
% % BIC
% BICs = [int_reward.all.predicted.BIC int_1LR_reward.all.predicted.BIC int_2LR_reward.all.predicted.BIC];
% names = categorical({'Intercept', '1 LR', '2 LRs'});
% names = reordercats(names,{'Intercept', '1 LR', '2 LRs'});
% % Plot
% subplot(1,2,2)
% bar(names,BICs)
% title('BIC')
% set(gcf,'units','centimeters','position',[0 0 30 10]) 
% 
% %Save
% saveas(gcf,[results_path_reward,'AIC_BIC_whole_sample_reward','.png'])

%% Delta AIC/BIC plots

% Calculate delta AIC

% Efficacy
delta_IC.AIC.efficacy.int_min_int_2LR = AIC_BIC.int_efficacy.AIC_BIC.AIC(:,1) - AIC_BIC.int_2LR_efficacy.AIC_BIC.AIC(:,1);
delta_IC.AIC.efficacy.int_1LR_min_int_2LR = AIC_BIC.int_1LR_efficacy.AIC_BIC.AIC(:,1) - AIC_BIC.int_2LR_efficacy.AIC_BIC.AIC(:,1);
% Reward
delta_IC.AIC.reward.int_min_int_2LR = AIC_BIC.int_reward.AIC_BIC.AIC(:,1) - AIC_BIC.int_2LR_reward.AIC_BIC.AIC(:,1);
delta_IC.AIC.reward.int_1LR_min_int_2LR = AIC_BIC.int_1LR_reward.AIC_BIC.AIC(:,1) - AIC_BIC.int_2LR_reward.AIC_BIC.AIC(:,1);

% Calculate delta BIC

% Efficacy
delta_IC.BIC.efficacy.int_min_int_2LR = AIC_BIC.int_efficacy.AIC_BIC.BIC(:,1) - AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC(:,1);
delta_IC.BIC.efficacy.int_1LR_min_int_2LR = AIC_BIC.int_1LR_efficacy.AIC_BIC.BIC(:,1) - AIC_BIC.int_2LR_efficacy.AIC_BIC.BIC(:,1);
% Reward
delta_IC.BIC.reward.int_min_int_2LR = AIC_BIC.int_reward.AIC_BIC.BIC(:,1) - AIC_BIC.int_2LR_reward.AIC_BIC.BIC(:,1);
delta_IC.BIC.reward.int_1LR_min_int_2LR = AIC_BIC.int_1LR_reward.AIC_BIC.BIC(:,1) - AIC_BIC.int_2LR_reward.AIC_BIC.BIC(:,1);

% Plot
figure
sgtitle('Model fit comparison')
names = categorical({'Intercept','1 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate'});

% Plot efficacy

% AIC
Delta_AICs = [ delta_IC.AIC.efficacy.int_min_int_2LR delta_IC.AIC.efficacy.int_1LR_min_int_2LR ];
Delta_AICs_means = mean(Delta_AICs);
subplot(2,2,1)
bar(names,Delta_AICs_means,0.2,'FaceColor','white')
hold
plot(names,Delta_AICs,'o','color','black')
yline(0)
xlabel('Delta AIC compared to the 2 Learning rates model')
ylabel('Model')
title('Efficacy Delta AIC')

% BIC
Delta_BICs = [delta_IC.BIC.efficacy.int_min_int_2LR delta_IC.BIC.efficacy.int_1LR_min_int_2LR];
Delta_BICs_means = mean(Delta_BICs);
subplot(2,2,2)
bar(names,Delta_BICs_means,0.2,'FaceColor','white')
hold
plot(names,Delta_BICs,'o','color','black')
yline(0)
xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('Model')
title('Efficacy Delta BIC')

% Plot reward

% AIC
Delta_AICs = [ delta_IC.AIC.reward.int_min_int_2LR delta_IC.AIC.reward.int_1LR_min_int_2LR];
Delta_AICs_means = mean(Delta_AICs);
subplot(2,2,3)
bar(names,Delta_AICs_means,0.2,'FaceColor','white')
hold
plot(names,Delta_AICs,'o','color','black')
yline(0)
xlabel('Delta AIC compared to the 2 Learning rates model')
ylabel('Model')
title('Reward Delta AIC')

% BIC
Delta_BICs = [ delta_IC.BIC.reward.int_min_int_2LR delta_IC.BIC.reward.int_1LR_min_int_2LR];
Delta_BICs_means = mean(Delta_BICs);
subplot(2,2,4)
bar(names,Delta_BICs_means,0.2,'FaceColor','white')
hold
plot(names,Delta_BICs,'o','color','black')
yline(0)
xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('Model')
title('Reward Delta BIC')
%Save
set(gcf,'units','centimeters','position',[0 0 20 20]) 
saveas(gcf,[results_path_model_comparison,'Delta_AIC_BIC','.png'])


% BIC for efficacy only
figure
names = categorical({'Intercept','1 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate'});
Delta_BICs = [delta_IC.BIC.efficacy.int_min_int_2LR delta_IC.BIC.efficacy.int_1LR_min_int_2LR];
Delta_BICs_means = mean(Delta_BICs);
bar(names,Delta_BICs_means,0.2,'FaceColor','white','LineWidth', 1.4)
hold
plot(names,Delta_BICs,'o','color','black','LineWidth', 1.4)
yline(0,'LineWidth', 1.4)
% xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('\Delta BIC')
%Save
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gca,'linewidth',2)
saveas(gcf,[results_path_model_comparison,'Delta_BIC_efficacy','.png'])
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure5_B.pdf')

% BIC for reward only
figure
names = categorical({'Intercept','1 Learning rate'}); % Define the names of the categorical axes
names = reordercats(names,{'Intercept','1 Learning rate'});
Delta_BICs = [delta_IC.BIC.reward.int_min_int_2LR delta_IC.BIC.reward.int_1LR_min_int_2LR];
Delta_BICs_means = mean(Delta_BICs);
bar(names,Delta_BICs_means,0.2,'FaceColor','white','LineWidth', 1.4)
hold
plot(names,Delta_BICs,'o','color','black','LineWidth', 1.4)
yline(0,'LineWidth', 1.4)
% xlabel('Delta BIC compared to the 2 Learning rates model')
ylabel('\Delta BIC')
set(gca, 'FontName', 'Helvetica','fontsize',20)
set(gcf,'units','centimeters','position',[0 0 22 16])
set(gca,'linewidth',2)
saveas(gcf,[results_path_model_comparison,'Delta_BIC_reward','.png'])
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS6_B.pdf')


