% Code for plotting ERPs and topographies for the LFXC_EEG project
% Written by: Romy Froemer and Ivan Grahek


%% Importing data

% PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/';
PATH = 'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/';

% load behavioral Data
load(sprintf('%sData/Experiment1/EEG/Export/FXCallSubDataTable.mat', PATH))
% take out the last participant's behavioral data because we don't have the
% clean eeg data yet
% delete = FXCallSubDataTable.SubID==1077;
% new = FXCallSubDataTable;
% new(delete,:) = [];
% FXCallSubDataTable = new;

% load EEG Data
load(sprintf('%sData/Experiment1/EEG/Export/CDAT.mat', PATH))
load(sprintf('%sData/Experiment1/EEG/Export/SDAT.mat', PATH))
load(sprintf('%sData/Experiment1/EEG/Export/RDAT.mat', PATH))
load(sprintf('%sData/Experiment1/EEG/Export/FRDAT.mat', PATH))
load(sprintf('%sData/Experiment1/EEG/Export/FEDAT.mat', PATH))

% load the channel info
nchans = size(CDAT,1);
load(sprintf('%sData/Experiment1/EEG/Export/chanlocs.mat', PATH));
for n = 1:nchans
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end

%get unique subject names for later processing and preselection of partially existing data
[~,idx] = unique(FXCallSubDataTable.SubID,'first');
vps = FXCallSubDataTable.SubID(sort(idx));

% load & define stuff for plotting
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;

%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % for topoplot function


%% Preparing the data for plotting the ERPs


%% Defining the timining for the cue-locked activity
cTIME    = linspace(-200, 1500, size(CDAT,2));
XMIN = -200;
XMAX =  1500;
CERP=nanmean(CDAT,3);
%% Defining the timining for the everything else-locked activity
TIME    =linspace(-200, 800, size(SDAT,2));
XMIN = -200;
XMAX =  800;
SERP=nanmean(SDAT,3);
RERP=nanmean(RDAT,3);
FRERP=nanmean(FRDAT,3);
FEERP=nanmean(FEDAT,3);
%% Median split or tertile split for the running average efficacy level

% Running average split for efficacy
FXCallSubDataTable.QSRunEffLev = zeros(size(FXCallSubDataTable.runAvgEfficacy));
for nsubs = 1: length(vps)

cut= quantile(FXCallSubDataTable.runAvgEfficacy(FXCallSubDataTable.SubID==vps(nsubs)), [.25 .50 .75]);

% FXCallSubDataTable.QSRunEffLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgEfficacy<cut(1,1))=1;
% FXCallSubDataTable.QSRunEffLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgEfficacy>cut(1,3))=2;

FXCallSubDataTable.QSRunEffLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgEfficacy < 0.5)=1;
FXCallSubDataTable.QSRunEffLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgEfficacy > 0.5)=2;
end

% Running average split for reward
FXCallSubDataTable.QSRunRewLev = zeros(size(FXCallSubDataTable.runAvgRewRate));
for nsubs = 1: length(vps)

% cut= quantile(FXCallSubDataTable.runAvgRewRate(FXCallSubDataTable.SubID==vps(nsubs)), [.25 .50 .75]);
% 
% FXCallSubDataTable.QSRunRewLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgRewRate<cut(1,1))=1;
% FXCallSubDataTable.QSRunRewLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgRewRate>cut(1,3))=2;

FXCallSubDataTable.QSRunRewLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgRewRate < 0.5)=1;
FXCallSubDataTable.QSRunRewLev(FXCallSubDataTable.SubID==vps(nsubs) & FXCallSubDataTable.runAvgRewRate > 0.5)=2;
end

% For the CNV
RunEffLow= find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.isMiss~=1);
RunEffHigh= find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.isMiss~=1);

% For the efficacy feedback
RunEffLowE1 = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.isMiss~=1 );
RunEffHighE1 = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.isMiss~=1);
RunEffLowE0 = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.isMiss~=1);
RunEffHighE0 = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.isMiss~=1);

% For the efficacy feedback when running average reward is high
RunEffLowRunRewHighR1 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 );
RunEffHighRunRewHighR1 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1);
RunEffLowRunRewHighR0 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1);
RunEffHighRunRewHighR0 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1);

% For the efficacy feedback when running average reward is high only reward
% feedback first
% RunEffLowRunRewHighR1 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
% RunEffHighRunRewHighR1 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1& FXCallSubDataTable.feedbackOrder==2);
% RunEffLowRunRewHighR0 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1& FXCallSubDataTable.feedbackOrder==2);
% RunEffHighRunRewHighR0 = find(FXCallSubDataTable.QSRunRewLev==2 & FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1& FXCallSubDataTable.feedbackOrder==2);
%% Reward feedback
RunEffLowR1 = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 );
RunEffHighR1 = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1);
RunEffLowR0 = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1);
RunEffHighR0 = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1);
% E1R01 = find(FXCallSubDataTable.MSSubEffLev==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);
% E0R01 = find(FXCallSubDataTable.MSSubEffLev==0 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);
% E1R11 = find(FXCallSubDataTable.MSSubEffLev==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);
% E0R11 = find(FXCallSubDataTable.MSSubEffLev==0 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);


% Reward vs. no reward when reward feedback first
IR_rew1st = find(FXCallSubDataTable.feedbackOrder==2 & FXCallSubDataTable.IsRewarded==1 & ~isnan(FXCallSubDataTable.Resp));
NR_rew1st = find(FXCallSubDataTable.feedbackOrder==2 & FXCallSubDataTable.IsRewarded==0 & ~isnan(FXCallSubDataTable.Resp));

% Reward vs. no reward when efficacy feedback first
IR_eff1st = find(FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==1 & ~isnan(FXCallSubDataTable.Resp));
NR_eff1st = find(FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==0 & ~isnan(FXCallSubDataTable.Resp));

% Reward vs. no reward and Eff vs. No Eff when efficacy feedback first
IR_eff1st_eff = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==1 & ~isnan(FXCallSubDataTable.Resp));
NR_eff1st_eff = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==0 & ~isnan(FXCallSubDataTable.Resp));
IR_eff1st_noeff = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==1 & ~isnan(FXCallSubDataTable.Resp));
NR_eff1st_noeff = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.feedbackOrder==1 & FXCallSubDataTable.IsRewarded==0 & ~isnan(FXCallSubDataTable.Resp));

% just Reward
IR = find(FXCallSubDataTable.IsRewarded==1 & ~isnan(FXCallSubDataTable.Resp));
NR = find(FXCallSubDataTable.IsRewarded==0 & ~isnan(FXCallSubDataTable.Resp));
%% Efficacy feedback
E1 = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.isMiss~=1);
E0 = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.isMiss~=1);

E_first_Eff = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);
E_first_NoEff = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==1);

R_first_Rew_Eff = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
R_first_Rew_NoEff = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);

R_first_NoRew_Eff = find(FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
R_first_NoRew_NoEff = find(FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);


% E1R02 = find(FXCallSubDataTable.MSSubEffLev==1 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
% E0R02 = find(FXCallSubDataTable.MSSubEffLev==0 & FXCallSubDataTable.IsRewarded==0 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
% E1R12 = find(FXCallSubDataTable.MSSubEffLev==1 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
% E0R12 = find(FXCallSubDataTable.MSSubEffLev==0 & FXCallSubDataTable.IsRewarded==1 & FXCallSubDataTable.isMiss~=1 & FXCallSubDataTable.feedbackOrder==2);
%% Congruency
C = find(FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
I = find(FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
N = find(FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);
%% Response - correct vs. incorrect

% Only correct vs. incorrect
Incorrect = find(FXCallSubDataTable.Acc==0 & FXCallSubDataTable.isMiss~=1);
Correct = find(FXCallSubDataTable.Acc==1 & FXCallSubDataTable.isMiss~=1);

% Correct vs. incorrect by congruency
Incorrect_Inc = find(FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
Incorrect_Con = find(FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
Incorrect_Neu = find(FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);
Correct_Inc = find(FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
Correct_Con = find(FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
Correct_Neu = find(FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);

% Correct vs. incorrect by congruency by running efficacy 
HighEfficacy_Incorrect_Inc = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
HighEfficacy_Incorrect_Con = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
HighEfficacy_Incorrect_Neu = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);
HighEfficacy_Correct_Inc = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
HighEfficacy_Correct_Con = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
HighEfficacy_Correct_Neu = find(FXCallSubDataTable.QSRunEffLev==2 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);

LowEfficacy_Incorrect_Inc = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
LowEfficacy_Incorrect_Con = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
LowEfficacy_Incorrect_Neu = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==0 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);
LowEfficacy_Correct_Inc = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==0 & FXCallSubDataTable.isMiss~=1);
LowEfficacy_Correct_Con = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==1 & FXCallSubDataTable.isMiss~=1);
LowEfficacy_Correct_Neu = find(FXCallSubDataTable.QSRunEffLev==1 & FXCallSubDataTable.Acc==1 & FXCallSubDataTable.Congruency==2 & FXCallSubDataTable.isMiss~=1);
%% Plotting the ERPs


%% Midline 
YMIN=-10;
YMAX=10;
XMAX=800;
XMIN = -200;
figure
hold on;
title('Midline', 'fontsize', 14)
%  f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
%  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
%  f3 = fill([700 700 3000 3000],[YMIN YMAX YMAX YMIN], [.98 .98 .98]);
%  f3.EdgeColor = 'none';
%  f3.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(cTIME,CERP(Pz,:),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,SERP(Pz,:),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
Leg3=plot(TIME,RERP(Pz,:),'-','linewidth',2, 'color', [0.80 0.475, 0.655]);
Leg4=plot(TIME,FRERP(Pz,:),'-','linewidth',2, 'color', [0.80 0.8, 0.8]); %rgb(80%,47.5%,65.5%)
Leg5=plot(TIME,FEERP(Pz,:),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('Time [ms]', 'fontsize', 14)
ylabel('Amplitude [µV]', 'fontsize', 14)
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [-600 -0.8 0]);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2, Leg3, Leg4, Leg5],'Cue', 'Target','Response', 'Reward Feedback', 'Efficacy Feedback');%, Leg4
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

saveas(figure(1),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPMidline.png')
%% Cue

% CNV by Running efficacy
YMIN=-2;
YMAX=1;
XMIN = -200;
XMAX =  1500;
figure;
hold on;
% title('FCz - running average efficacy (median split)', 'fontsize', 12)
 f = fill([1000 1000 1500 1500],[YMIN YMAX YMAX YMIN], [.85 .85 .85]);
 f.EdgeColor = 'none';
 f.EdgeColor = 'none';
%  f2 = fill([700 700 1000 1000],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(cTIME,nanmean(CDAT(FCz,:,RunEffHigh),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(cTIME,nanmean(CDAT(FCz,:,RunEffLow),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(cTIME,nanmean(CDAT(Cz,:,E1R1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(cTIME,nanmean(CDAT(Cz,:,E0R1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]')
ylabel('amplitude [µV]')
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
vline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
set(vline,'Color','black','Linestyle',':','LineWidth',1);

% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'High efficacy estimate', 'Low efficacy estimate');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.2,.81,lp(3:4)])
set(gcf,'units','centimeters','position',[0 0 15 20])


% h=gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'CNV_FCz_RunAvgEfficacy','-dpdf','-r0')


saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERP_CNVlate_FCz.png')


% CNV
YMIN=-2;
YMAX=1;
XMIN = -200;
XMAX =  1500;
figure;
hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy', 'fontsize', 14)
f = fill([350 350 500 500],[YMIN (YMAX-0.5) (YMAX-0.5) YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(cTIME,nanmean(CDAT(FCz,:,RunEffHigh),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(cTIME,nanmean(CDAT(FCz,:,RunEffLow),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]','FontName', 'Helvetica', 'fontsize', 18)
ylabel({'\fontsize{20}P3b (350-500ms)','\fontsize{15}amplitude [µV]'},'FontName', 'Helvetica')
set(gca, 'fontsize', 14,'FontName', 'Helvetica');
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'High Efficacy', 'Low Efficacy');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14,'FontName', 'Helvetica')%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.29,.81,lp(3:4)])

% set(gcf,'units','centimeters','position',[0 0 15 20])
set( gcf,'PaperSize',[5 8], 'PaperPosition',[0 0 5 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPEfficacyFB_P3b_Pz.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure4_B.pdf')





% P3b
% YMIN=-6;
% YMAX=8;
% 
% XMIN = -200;
% XMAX =  1500;
% figure();
% hold on;
% title('Pz Cue Efficacy', 'fontsize', 14)
%  f = fill([250 250 550 550],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
% %  f2 = fill([700 700 1000 1000],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(cTIME,nanmean(CDAT(Pz,:,E1MS),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(cTIME,nanmean(CDAT(Pz,:,E0MS),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% % Leg3=plot(cTIME,nanmean(CDAT(Pz,:,E1R1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% % Leg4=plot(cTIME,nanmean(CDAT(Pz,:,E0R1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'high efficacy', 'low efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
%% Target 

YMIN=-6;
YMAX=6;
XMIN = -200;
XMAX =  800;
figure();
hold on;
title('FCz Target', 'fontsize', 14)
 f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
 f.EdgeColor = 'none';
 f.EdgeColor = 'none';
%  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(SDAT(FCz,:,C),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(SDAT(FCz,:,I),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
Leg3=plot(TIME,nanmean(SDAT(FCz,:,N),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(SDAT(FCz,:,E0R1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]')
ylabel('amplitude [µV]')
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2, Leg3],'Congruent', 'Incongruent', 'Neutral');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

saveas(figure(1),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPTargetN2.png')
%% Response

% incorrect vs. correct
YMIN=-6;
YMAX=6;
XMIN = -200;
XMAX =  600;
figure();
hold on;
title('FCz Accuracy', 'fontsize', 12)
 f = fill([-30 -30 70 70],[YMIN YMAX YMAX YMIN], [.85 .85 .85]);
 f.EdgeColor = 'none';
 f.EdgeColor = 'none';
%  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(RDAT(FCz,:,Incorrect),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(RDAT(FCz,:,Correct),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]')
ylabel('amplitude [µV]')
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Incorrect', 'Correct');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')

h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'CNV_FCz_RunAvgEfficacy','-dpdf','-r0')

saveas(figure(1),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPResponseERN.png')


% ERN by congruency
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure();
% hold on;
% title('FCz Correct vs. Incorrect times congruency', 'fontsize', 14)
%  f = fill([0 0 100 100],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% 
% Leg1=plot(TIME,nanmean(RDAT(FCz,:,Incorrect_Inc),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(RDAT(FCz,:,Incorrect_Con),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(RDAT(FCz,:,Incorrect_Neu),3),'-','linewidth',2, 'color', [0. 0.6, 0.8]);
% Leg4=plot(TIME,nanmean(RDAT(FCz,:,Correct_Inc),3),'-','linewidth',2, 'color', [0.5 0.6, 0.2]); % rgb(83.5%,36.9%,0%)
% Leg5=plot(TIME,nanmean(RDAT(FCz,:,Correct_Con),3),'-','linewidth',2, 'color', [0.3 0.7, 0.18]);
% Leg6=plot(TIME,nanmean(RDAT(FCz,:,Correct_Neu),3),'-','linewidth',2, 'color', [0.5 0.39, 0.3]); % rgb(83.5%,36.9%,0%)
% 
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4, Leg5, Leg6],'Incorrect Inc', 'Incorrect Con','Incorrect Neu', 'Correct Inc', 'Correct Con','Correct Neu');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')

% ERN by running efficacy by congruency
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure();
% hold on;
% title('FCz Correct vs. Incorrect', 'fontsize', 14)
%  f = fill([0 0 100 100],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% 
% Leg1=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Incorrect_Inc),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Incorrect_Con),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Incorrect_Neu),3),'-','linewidth',2, 'color', [0. 0.6, 0.8]);
% Leg4=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Correct_Inc),3),'-','linewidth',2, 'color', [0.5 0.6, 0.2]); % rgb(83.5%,36.9%,0%)
% Leg5=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Correct_Con),3),'-','linewidth',2, 'color', [0.3 0.7, 0.18]);
% Leg6=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Correct_Neu),3),'-','linewidth',2, 'color', [0.5 0.39, 0.3]); % rgb(83.5%,36.9%,0%)
% Leg7=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Incorrect_Inc),3),'-','linewidth',2, 'color', [0.2 0.17, 0.98]);
% Leg8=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Incorrect_Con),3),'-','linewidth',2, 'color', [0.35 0.1169, 0.5]); % rgb(83.5%,36.9%,0%)
% Leg9=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Incorrect_Neu),3),'-','linewidth',2, 'color', [0.9 0.7, 0.18]);
% Leg10=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Correct_Inc),3),'-','linewidth',2, 'color', [0.1 0.39, 0.3]); % rgb(83.5%,36.9%,0%)
% Leg11=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Correct_Con),3),'-','linewidth',2, 'color', [0.27 0.37, 0.98]);
% Leg12=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Correct_Neu),3),'-','linewidth',2, 'color', [0.5 0.169, 0.65]); % rgb(83.5%,36.9%,0%)
% 
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4, Leg5, Leg6, Leg7, Leg8,Leg9, Leg10, Leg11, Leg12],'HighEfficacy_Incorrect_Inc', ...
% 'HighEfficacy_Incorrect_Con',...
% 'HighEfficacy_Incorrect_Neu',...
% 'HighEfficacy_Correct_Inc',...
% 'HighEfficacy_Correct_Con',...
% 'HighEfficacy_Correct_Neu',...
% 'LowEfficacy_Incorrect_Inc',...
% 'LowEfficacy_Incorrect_Con',...
% 'LowEfficacy_Incorrect_Neu',...
% 'LowEfficacy_Correct_Inc',...
% 'LowEfficacy_Correct_Con', ...
% 'LowEfficacy_Correct_Neu');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')


% ERN by running efficacy by congruency
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure();
% hold on;
% title('FCz Correct vs. Incorrect', 'fontsize', 14)
%  f = fill([0 0 100 100],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% 
% Leg1=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Incorrect_Inc),3),'-','linewidth',2, 'color', [0.9 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(RDAT(FCz,:,HighEfficacy_Correct_Inc),3),'-','linewidth',2, 'color', [0.5 0.6, 0.6]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Incorrect_Inc),3),'-','linewidth',2, 'color', [0.3 0.17, 0.98]);
% Leg4=plot(TIME,nanmean(RDAT(FCz,:,LowEfficacy_Correct_Inc),3),'-','linewidth',2, 'color', [0.1 0.5, 0.3]); % rgb(83.5%,36.9%,0%)
% 
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'HighEffIncorrectInc', ...
% 'HighEffCorrectInc',...
% 'LowEffIncorrectInc',...
% 'LowEffCorrectInc');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
%% Reward feedback

% Reward vs no reward (P3a)
YMIN=-1.2;
YMAX=5;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('FCz (P3a) Reward feedback: Reward vs. No Reward', 'fontsize', 12)
f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FRDAT(FCz,:,NR),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(FRDAT(FCz,:,IR),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]')
ylabel({'\fontsize{20}P3a (357-457ms)','\fontsize{15}amplitude [µV]'})
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'No reward', 'Reward');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.22,.81,lp(3:4)])
set(gcf,'units','centimeters','position',[0 0 15 20])
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPRewardFB_P3a_FCz.png')
% % Running Efficacy and Reward - P3a
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('FCz (P3a) Reward Feedback: Run Efficacy X Reward', 'fontsize', 14)
% f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,RunEffHighR0),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,RunEffLowR0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FRDAT(FCz,:,RunEffHighR1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FRDAT(FCz,:,RunEffLowR1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'high efficacy no reward', 'low efficacy no reward', 'high efficacy reward', 'low efficacy  reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPRewardFB_EffSplit_P3a_FCz.png')

% Reward vs no reward (P3b)
YMIN=-1.2;
YMAX=5;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('Pz (P3b) Reward feedback: Reward vs. No Reward', 'fontsize', 12)
f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FRDAT(Pz,:,NR),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(FRDAT(Pz,:,IR),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]')
ylabel({'\fontsize{20}P3b (350-500ms)','\fontsize{15}amplitude [µV]'})
set(gca, 'fontsize', 14);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'No reward', 'Reward');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.22,.81,lp(3:4)])
set(gcf,'units','centimeters','position',[0 0 15 20])
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPRewardFB_P3b_Pz.png')


% P3b Reward vs. no reward
YMIN=-1.2;
YMAX=3;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy', 'fontsize', 14)
f = fill([350 350 500 500],[YMIN (YMAX-0.5) (YMAX-0.5) YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FRDAT(Pz,:,NR),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(FRDAT(Pz,:,IR),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]','FontName', 'Helvetica', 'fontsize', 18)
ylabel({'\fontsize{20}P3b (350-500ms)','\fontsize{15}amplitude [µV]'},'FontName', 'Helvetica')
set(gca, 'fontsize', 14,'FontName', 'Helvetica');
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'No Reward', 'Reward');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14,'FontName', 'Helvetica')%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.29,.81,lp(3:4)])

% set(gcf,'units','centimeters','position',[0 0 15 20])
set( gcf,'PaperSize',[5 8], 'PaperPosition',[0 0 5 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPEfficacyFB_P3b_Pz.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS3_A.pdf')

% 
% % Running Efficacy and Reward - P3b
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Reward Feedback: Run Efficacy X Reward', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(Pz,:,RunEffHighR0),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(Pz,:,RunEffLowR0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FRDAT(Pz,:,RunEffHighR1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FRDAT(Pz,:,RunEffLowR1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'high efficacy no reward', 'low efficacy no reward', 'high efficacy reward', 'low efficacy  reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPRewardFB_EffSplit_P3b_Pz.png')

% % reward vs no reward when reward feedback first
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  600;
% figure(2);
% hold on;
% title('FCz Reward vs. No Reward (reward feedback 1st)', 'fontsize', 12)
% %  f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% %  f.EdgeColor = 'none';
% %  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,NR_rew1st),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,IR_rew1st),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'no reward', 'is reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% 
% h=gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'CNV_FCz_RunAvgEfficacy','-dpdf','-r0')
% 
% 
% saveas(figure(2),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPRewardFB_rew1st.png')

% % reward vs no reward when efficacy feedback first
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  600;
% figure(3);
% hold on;
% title('FCz Reward vs. No Reward (efficacy feedback 1st)', 'fontsize', 12)
% %  f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% %  f.EdgeColor = 'none';
% %  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,NR_eff1st),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,IR_eff1st),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'no reward', 'is reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% 
% h=gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'CNV_FCz_RunAvgEfficacy','-dpdf','-r0')
% 
% saveas(figure(3),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPRewardFB_eff1st.png')

% % reward vs no reward and eff vs. no efficacy when efficacy feedback first
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure(4);
% hold on;
% title('FCz Reward vs. No Reward and Eff vs. No Eff when efficacy feedback 1st', 'fontsize', 14)
% %  f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% %  f.EdgeColor = 'none';
% %  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,IR_eff1st_eff),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,NR_eff1st_eff),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FRDAT(FCz,:,IR_eff1st_noeff),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FRDAT(FCz,:,NR_eff1st_noeff),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'reward efficacy', 'no reward efficacy', 'reward no efficacy', 'no reward no efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% 
% 
% saveas(figure(4),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPRewardFB_eff1st_eff_vs_noeff.png')



% % reward vs no reward for P3a
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  550;
% figure(6);
% hold on;
% title('Cz Reward vs. No Reward', 'fontsize', 12)
%  f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.85 .85 .85]);
%  f.EdgeColor = 'none';
%  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,NR),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,IR),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'no reward', 'is reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% 
% h=gcf;
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'CNV_FCz_RunAvgEfficacy','-dpdf','-r0')
% 
% 
% saveas(figure(6),'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPP3aRewardFB.png')

% Efficacy and reward
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure();
% hold on;
% title('FCz Reward Feedback Running Average Efficacy by reward', 'fontsize', 14)
% %  f = fill([300 300 400 400],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% %  f.EdgeColor = 'none';
% %  f.EdgeColor = 'none';
% %  f2 = fill([500 500 600 600],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
% %  f2.EdgeColor = 'none';
% %  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(FCz,:,RunEffHighR0),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(FCz,:,RunEffLowR0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FRDAT(FCz,:,RunEffHighR1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FRDAT(FCz,:,RunEffLowR1),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'high efficacy no reward', 'low efficacy no reward', 'high efficacy reward', 'low efficacy  reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')



% Efficacy and reward when running average reward is high
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure();
% hold on;
% title('Cz Reward Feedback X Running Average Efficacy when Running Average Reward is High (median split used for high vs. low)', 'fontsize', 14)
% %  f = fill([283 283 383 383],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% %  f.EdgeColor = 'none';
% %  f.EdgeColor = 'none';
%  f2 = fill([347 347 447 447],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FRDAT(Cz,:,RunEffLowRunRewHighR1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FRDAT(Cz,:,RunEffHighRunRewHighR1),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FRDAT(Cz,:,RunEffLowRunRewHighR0),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FRDAT(Cz,:,RunEffHighRunRewHighR0),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% % set(gca, 'fontsize', 14, 'YDir','reverse');
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'low efficacy reward', 'high efficacy reward', 'low efficacy no reward', 'high efficacy no reward');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')

% h=gcf; set(h,'Units','Inches'); pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),
% pos(4)])
% print(h,'P3a_Cz_RewFeedback_RunAvgRewHigh_Efficacy_by_Reward','-dpdf','-r0')
%% Efficacy Feedback

% P3a Efficacy vs. no efficacy
YMIN=-1.2;
YMAX=3;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('FCz (P3a) Efficacy Feedback: Efficacy vs. No efficacy', 'fontsize', 14)
f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FEDAT(FCz,:,E1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(FEDAT(FCz,:,E0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]','FontName', 'Helvetica', 'fontsize', 18)
% ylabel('amplitude [µV]')
ylabel({'\fontsize{20}P3a (357-457ms)','\fontsize{15}amplitude [µV]'},'FontName', 'Helvetica')
set(gca, 'fontsize', 18,'FontName', 'Helvetica');
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
% vline = xline(0);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% set(vline,'Color','black','Linestyle',':','LineWidth',1);

% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Performance-Based', 'Random');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 18,'FontName', 'Helvetica')%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.26,.81,lp(3:4)])
% set(gca, 'FontName', 'Helvetica')

set(gcf,'units','centimeters','position',[0 0 15 20])




% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPEfficacyFB_P3a_FCz.png')

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_A.pdf')


% % P3a Efficacy vs. no efficacy by high and low running efficacy
% YMIN=-2;
% YMAX=4;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('FCz (P3a) Efficacy Feedback: Efficacy vs. No efficacy X Running efficacy', 'fontsize', 14)
% f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(FCz,:,RunEffLowE1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(FCz,:,RunEffLowE0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FEDAT(FCz,:,RunEffHighE1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FEDAT(FCz,:,RunEffHighE0),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'efficacy low running efficacy', 'no efficacy low running efficacy', 'efficacy high running efficacy', 'no efficacy high running efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3a_FCz.png')


% P3b Efficacy vs. no efficacy
YMIN=-1.2;
YMAX=2.5;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy', 'fontsize', 14)
f = fill([350 350 500 500],[YMIN (YMAX-0.5) (YMAX-0.5) YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FEDAT(Pz,:,E1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698],'LineWidth',2);
Leg2=plot(TIME,nanmean(FEDAT(Pz,:,E0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]','FontName', 'Helvetica', 'fontsize', 18)
ylabel({'\fontsize{20}P3b (350-500ms)','\fontsize{15}amplitude [µV]'},'FontName', 'Helvetica')
set(gca, 'fontsize', 14,'FontName', 'Helvetica');
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Performance-Based', 'Random');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 14,'FontName', 'Helvetica')%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.29,.81,lp(3:4)])

% set(gcf,'units','centimeters','position',[0 0 15 20])
set( gcf,'PaperSize',[5 8], 'PaperPosition',[0 0 5 8]) 

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/ERPEfficacyFB_P3b_Pz.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_A.pdf')


% P3a Efficacy vs. no efficacy
YMIN=-1.2;
YMAX=3;
XMIN = -200;
XMAX =  800;
figure;
hold on;
% title('FCz (P3a) Efficacy Feedback: Efficacy vs. No efficacy', 'fontsize', 14)
f = fill([357 357 457 457],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
f.EdgeColor = 'none';
f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(TIME,nanmean(FEDAT(FCz,:,E1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
Leg2=plot(TIME,nanmean(FEDAT(FCz,:,E0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('time [ms]','FontName', 'Helvetica', 'fontsize', 18)
% ylabel('amplitude [µV]')
ylabel({'\fontsize{20}P3a (357-457ms)','\fontsize{15}amplitude [µV]'},'FontName', 'Helvetica')
set(gca, 'fontsize', 18,'FontName', 'Helvetica');
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
% vline = xline(0);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
% set(vline,'Color','black','Linestyle',':','LineWidth',1);

% t=gridxy(0, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2],'Performance-Based', 'Random');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 18,'FontName', 'Helvetica')%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(lh,'position',[.26,.81,lp(3:4)])
% set(gca, 'FontName', 'Helvetica')

set(gcf,'units','centimeters','position',[0 0 15 20])






% 
% 
% % P3a Efficacy vs. no efficacy by high and low running efficacy
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy X Running efficacy', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(Pz,:,RunEffLowE1),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(Pz,:,RunEffLowE0),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FEDAT(Pz,:,RunEffHighE1),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FEDAT(Pz,:,RunEffHighE0),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4],'efficacy low running efficacy', 'no efficacy low running efficacy', 'efficacy high running efficacy', 'no efficacy high running efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3b_Pz_byRunAvgEff.png')
% 
% % P3a Efficacy vs. no efficacy by feedback order
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. Feedabck order', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(Pz,:,E_first_Eff),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(Pz,:,E_first_NoEff),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% Leg3=plot(TIME,nanmean(FEDAT(Pz,:,R_first_Rew_Eff),3),'-','linewidth',2, 'color', [0.80 0.475, 0.655]); %rgb(80%,47.5%,65.5%)
% Leg4=plot(TIME,nanmean(FEDAT(Pz,:,R_first_Rew_NoEff),3),'-','linewidth',2, 'color', [0 0.62, 0.451]); %rgb(80%,47.5%,65.5%)
% Leg5=plot(TIME,nanmean(FEDAT(Pz,:,R_first_NoRew_Eff),3),'-','linewidth',2, 'color', [0.50 0.475, 0.355]); %rgb(80%,47.5%,65.5%)
% Leg6=plot(TIME,nanmean(FEDAT(Pz,:,R_first_NoRew_NoEff),3),'-','linewidth',2, 'color', [0.2 0.32, 0.1]); %rgb(80%,47.5%,65.5%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2, Leg3, Leg4,Leg5,Leg6],'Eff first Efficacy', 'Eff first No Efficacy', 'Rew Efficacy', 'Rew No Efficacy', 'No Rew Efficacy', 'No Rew No Efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3b_Pz_by_feedback.png')
% 
% % P3a Efficacy vs. no efficacy when efficacy first
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy in Eff first trials', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(Pz,:,E_first_Eff),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(Pz,:,E_first_NoEff),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'Eff first Efficacy', 'Eff first No Efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3b_Pz_when_eff_first.png')
% 
% 
% % P3a Efficacy vs. no efficacy by feedback order
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy vs. No efficacy in Reward first trials', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(Pz,:,R_first_Rew_Eff),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(Pz,:,R_first_Rew_NoEff),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'Rew Efficacy', 'Rew No Efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3b_Pz_when_rew_first.png')
% 
% % P3a Efficacy vs. no efficacy by feedback order
% YMIN=-6;
% YMAX=6;
% XMIN = -200;
% XMAX =  800;
% figure;
% hold on;
% title('Pz (P3b) Efficacy Feedback: Efficacy in No Reward first trials', 'fontsize', 14)
% f = fill([350 350 500 500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
% f.EdgeColor = 'none';
% f.EdgeColor = 'none';
% set(gcf,'Color', [1 1 1])
% Leg1=plot(TIME,nanmean(FEDAT(Pz,:,R_first_NoRew_Eff),3),'-','linewidth',2, 'color', [0. 0.447, 0.698]);
% Leg2=plot(TIME,nanmean(FEDAT(Pz,:,R_first_NoRew_NoEff),3),'-','linewidth',2, 'color', [0.835 0.369, 0]); % rgb(83.5%,36.9%,0%)
% xlim([XMIN XMAX])
% ylim([YMIN YMAX])
% xlabel('time [ms]')
% ylabel('amplitude [µV]')
% set(gca, 'fontsize', 14);
% set(gca,'TickDir','out')
% xlabh = get(gca,'XLabel');
% hline = refline([0 0]);
% set(hline,'Color','black','Linestyle',':','LineWidth',1);
% t=gridxy(0, 'Linestyle',':','LineWidth',1);
% h = legend([Leg1, Leg2],'No Rew Efficacy', 'No Rew No Efficacy');
% lh=findall(gcf,'tag','legend');
% lp=get(lh,'position');
% set(lh,  'fontsize', 14)%, 'position',[.65,.85,lp(3:4)]
% set(lh,'box','off')
% set(gcf,'units','centimeters','position',[0 0 35 30])
% 
% saveas(gcf,'C:/Users/igrahek/Documents/Studies/Efficacy/LFXC_EEG/EEG/ERPEfficacyFB_EffSplit_P3b_Pz_when_noreward_first.png')

%% Plotting the LMM TOPOGRAPHIES


%% CUE LOCKED ACTIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CNV late (1000-1500ms)
load(sprintf('%sData/Experiment1/EEG/Export/ploteslateCNV.mat', PATH), 'EmatlateCNV') % estimates


% % Running average reward rate effect
% PLOTRAD = .70;
% HEADRAD = .70;
% INTRAD  =  0.9;
% MAPLIM = [-1 1];
% figure;
% title('Late CNV - Model-based reward (estimates)', 'fontsize', 12)
% topoplot(EmatlateCNV.mbased_reward_prev,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Fz,F1,F2,FCz,FC1,FC2,Cz,C1,C2],'o','k', 3,1});set(gca,'Color','none');
% xlim([-0.6 0.6])
% ylim([-0.6 0.6])
% colorbar('EastOutside');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 
% 
% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/CNVlateRunAvgRewRateEffect.png')


% Running average efficacy effect
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-0.8 0.8];
figure;
hold on;
% title('Late CNV - Model-based efficacy (estimates)', 'fontsize', 12)
topoplot(EmatlateCNV.mbased_efficacy_prev,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Fz,F1,F2,FCz,FC1,FC2,Cz,C1,C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');

% set( gcf,'PaperSize',[8 8], 'PaperPosition',[0 0 7 7]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/CNVlateRunAvgEfficacyEffect.png')
saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure4_C1.pdf')



%% EFFICACY FEEDBACK LOCKED ACTIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% P3b - estimate X efficacy level
load(sprintf('%sData/Experiment1/EEG/Export/plotesP3beff.mat', PATH), 'EmatP3beff') % estimates


% % Efficacy vs. no efficacy
% PLOTRAD = .70;
% HEADRAD = .70;
% INTRAD  =  0.9;
% MAPLIM = [-3 3];
% figure;
% % title('P3b Efficacy effect: Efficacy minus No Efficacy', 'fontsize', 12)
% topoplot(EmatP3beff.EffLvlEff_min_NEff,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
% xlim([-0.6 0.6])
% ylim([-0.6 0.6])
% colorbar('EastOutside');
% 
% set(gcf,'units','centimeters','position',[0 0 20 8]) 
% 
% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffEfficacyEffect.png')




% Running Avg Efficacy X Efficacy level
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-4 4];
figure;
% title('P3b - Efficacy estimate X Efficacy level', 'fontsize', 12)
topoplot(EmatP3beff.Eff_min_NEffbyEff_estimate,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffEffEstbyEffLvl.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_B1.pdf')


%% P3b - efficacy - PE X LR
load(sprintf('%sData/Experiment1/EEG/Export/plotesP3beff_PE_LR.mat', PATH), 'EmatP3beff') % estimates


% % Efficacy vs. no efficacy
% PLOTRAD = .70;
% HEADRAD = .70;
% INTRAD  =  0.9;
% MAPLIM = [-5 5];
% figure;
% title('P3b - Prediction error', 'fontsize', 12)
% topoplot(EmatP3beff.PE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
% xlim([-0.6 0.6])
% ylim([-0.6 0.6])
% colorbar('EastOutside');
% 
% set(gcf,'units','centimeters','position',[0 0 20 8]) 
% 
% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffPEEffect.png')

% Learning rate
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-4 4];
figure;
% title('P3b - Learning rate', 'fontsize', 12)
topoplot(EmatP3beff.LR,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffLR.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_C2.1.pdf')


% Prediction error
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-4 4];
figure;
% title('P3b - Prediction error X Learning rate', 'fontsize', 12)
topoplot(EmatP3beff.PE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffPEbyLR.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_C1.1.pdf')


% % Running Avg Efficacy X Efficacy level
% PLOTRAD = .70;
% HEADRAD = .70;
% INTRAD  =  0.9;
% MAPLIM = [-4 4];
% figure;
% % title('P3b - Prediction error X Learning rate', 'fontsize', 12)
% topoplot(EmatP3beff.PEbyLR,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
% xlim([-0.6 0.6])
% ylim([-0.6 0.6])
% colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% % set(gcf,'units','centimeters','position',[0 0 20 8]) 
% 
% % saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3beffPEbyLR.png')
% 
% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/Figure3_C1.1.pdf')



%% P3b - reward estimate X reward level
load(sprintf('%sData/Experiment1/EEG/Export/plotesP3brew.mat', PATH), 'EmatP3brew') % estimates


% Efficacy vs. no efficacy
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-4 4];
figure;
title('P3b Reward effect: Reward minus No Reward', 'fontsize', 12)
topoplot(EmatP3brew.IsRewardedRew_min_NRew,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');

set(gcf,'units','centimeters','position',[0 0 20 8]) 

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3brewRewardEffect.png')


% Running Avg Efficacy X Efficacy level
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure;
% title('P3b - Reward estimate X Reward level', 'fontsize', 12)
topoplot(EmatP3brew.Rew_min_NRewbyRew_estimate,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3brewRewEstbyRewardEffect.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS3_B1.pdf')




%% P3b - reward - PE X LR
load(sprintf('%sData/Experiment1/EEG/Export/plotesP3brew_PE_LR.mat', PATH), 'EmatP3brew') % estimates


% Efficacy vs. no efficacy
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure;
% title('P3b - Prediction error', 'fontsize', 12)
topoplot(EmatP3brew.PE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');

% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3brewPEEffect.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS3_C1.1.pdf')


% Learning rate
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure;
% title('P3b - Learning rate', 'fontsize', 12)
topoplot(EmatP3brew.LR,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3brewLR.png')

saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/Figures/FigureS3_C2.1.pdf')


% Running Avg Efficacy X Efficacy level
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-2 2];
figure;
% title('P3b - Prediction error X Learning rate', 'fontsize', 12)
topoplot(EmatP3brew.PEbyLR,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz,P1,P2,POz,PO3,PO4,CPz,CP1,CP2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside','fontsize', 14,'FontName', 'Helvetica');
% set(gcf,'units','centimeters','position',[0 0 20 8]) 

% saveas(gcf,'C:/Users/igrahek/Dropbox (Brown)/CLPS-ShenhavLab/EEG_Studies/Experiments/LFXC_EEG/Manuscript/plots/EEG/P3brewPEbyLR.png')

