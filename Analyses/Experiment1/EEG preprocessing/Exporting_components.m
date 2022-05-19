

%% Export time window average or time series to access in R for analyses  
% Author: Romy Froemer

% here we export Cue, Stimulus, response, reward and efficacy feedback
% related components. 

% A) Cue
% A1) late CNV
% A2) early CNV
% A3) Cue P3b

% B) Stimulus (N2)

% C) Reward Feedback
% C1) mean FRN
% C2) P3a
% C3) P3b
% C4) Peak to peak FRN

% D) Efficacy Feedback
% D1) mean FRN
% D2) P3a
% D3) P3b
% D4) Control for early effect onset on P3b: 200-300 ms covariate
% D5) Peak to peak FRN

% E) Response
% E1) ERN
% E2) Pe

%% general settings
PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/'; %edit to your needs
EXPORTPATH ='Data/Export/'; % Export folder (segmented data will be saved there)
srate=500; % sampling rate
bl=200; %length of prestimulus interval in ms

%% 1) load EEG data (if necessary)

load(sprintf('%sData/Export/FXCallSubDataTable.mat', PATH))
% load EEG Data
load(sprintf('%sData/Export/CDAT.mat', PATH)) % Cue
load(sprintf('%sData/Export/SDAT.mat', PATH)) % Stimulus
load(sprintf('%sData/Export/RDAT.mat', PATH)) % Response
load(sprintf('%sData/Export/FRDAT.mat', PATH))% Reward feedback 
load(sprintf('%sData/Export/FEDAT.mat', PATH))% Efficacy feedback

%% get channel locations for peak detection
nchans = size(CDAT,1);
load(sprintf('%sData/Export/chanlocs.mat', PATH));
for n = 1:nchans
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end

%% A) Cue
% for cue- locked potentials we go with the literature
BFS = CDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

%% A1) export mean CNV --> late 1000 - 1500
mint=1000; %export start in ms
maxt=1500; %export end in ms
ERPNAME='CNV'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

%% A2) export mean CNV --> early 700-1000
mint=700; %export start in ms
maxt=1000; %export end in ms
ERPNAME='CNV'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% A3) export mean Cue P3b
mint=300; %export start in ms
maxt=400; %export end in ms
ERPNAME='CueP3b'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');



%% B) stimulus locked: N2 
TIME    =linspace(-200, 800, size(SDAT,2));

GM = nanmean(SDAT,3);

% Peak detection
start= find(TIME>200,1);

[Y,I] = min(GM(FCz,start:end)); % find negative peak at FCz

Peak = TIME(start+I); % N2 peak at 331 something ms

% Get the data
BFS = SDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

% export mean N2--> 140-240 ms 
mint=Peak-50; %export start in ms
maxt=Peak+50; %export end in ms
ERPNAME='N2'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',round(mint),round(maxt)) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,round(mint),round(maxt)), sprintf('%s%d%d',ERPNAME,round(mint),round(maxt)), '-v6');


%% C) Reward feedback

% grand mean and time vector for peak detection
TIME    =linspace(-200, 800, size(FRDAT,2));
GM = nanmean(FRDAT,3);
% get relevant data

BFS = FRDAT; % set BFS to whatever your data is
nchans = size(BFS,1);


%% C1) Outdated peak FRN: find negative peak between 230 and 350; 
start= find(TIME>230,1);
maxtime = find(TIME<350,1,'last');
[Y,I] = min(GM(FCz,start:maxtime)); % find negative peak at FCz

Peak = round(TIME(start+I)); 

% export mean FRN 
mint=Peak - 50; %export start in ms
maxt=Peak + 50; %export end in ms
ERPNAME='FRNR'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% C2) now P3a, basically the same thing, only that we try to find a positive peak and later
% find positive peak between 330 and 500
start= find(TIME>330,1);
maxtime = find(TIME<500,1,'last');
[Y,I] = max(GM(FCz,start:maxtime)); % find negative peak at FCz

Peak = round(TIME(start+I)); 

% export mean P3a 
mint=Peak - 50; %export start in ms
maxt=Peak + 50; %export end in ms
ERPNAME='P3aR'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% C3) export mean P3b for reward Feedback
mint=350; %export start in ms
maxt=500; %export end in ms
ERPNAME='P3bR'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

%% C4) Peak to peak FRN Reward
P2PFCZR = p2p(FRDAT, 200, 300, FCz, 75);
%PATH = '/Volumes/daten/romy/Confidence/';
save(sprintf('%s%s/P2PFCZR.mat', PATH,EXPORTPATH), 'P2PFCZR', '-v6');



%% D) efficacy feedback
TIME    =linspace(-200, 800, size(FEDAT,2));
GM = nanmean(FEDAT,3);

% get relevant data

BFS = FEDAT; % set BFS to whatever your data is
nchans = size(BFS,1);


%% D1) Outdated mean FRN (peak detection)find negative peak between 230 and 350
start= find(TIME>230,1);
maxtime = find(TIME<350,1,'last');
[Y,I] = min(GM(FCz,start:maxtime)); % find negative peak at FCz

Peak = round(TIME(start+I)); % N2

% export mean FRN 
mint=Peak - 50; %export start in ms
maxt=Peak + 50; %export end in ms
ERPNAME='FRNE'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% D2) P3a
% find positive peak between 330 and 500 
start= find(TIME>330,1);
maxtime = find(TIME<500,1,'last');
[Y,I] = max(GM(FCz,start:maxtime)); % find negative peak at FCz

Peak = round(TIME(start+I)); 

% export mean P3a
mint=Peak - 50; %export start in ms
maxt=Peak + 50; %export end in ms
ERPNAME='P3aE'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% D3) mean Cue P3b for efficacy Feedback
mint=350; %export start in ms
maxt=500; %export end in ms
ERPNAME='P3bE'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

%% D4) get 200-300 ms Feedback locked activity for early effect onset control:

mint=200; %export start in ms
maxt=300; %export end in ms
ERPNAME='P3bEcontrol'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

%% D5) Peak to peak FRN Efficacy
P2PFCZE = p2p(FEDAT, 200, 300, FCz, 75);
%PATH = '/Volumes/daten/romy/Confidence/';
save(sprintf('%s%s/P2PFCZE.mat', PATH,EXPORTPATH), 'P2PFCZE', '-v6');


%% E) response locked 

BFS = RDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

TIME    =-200:2: 798;

GM = nanmean(RDAT,3);
%% E1) ERN peak
start= find(TIME>-50,1);
maxtime = find(TIME<150,1,'last');
[Y,I] = min(GM(FCz,start:maxtime)); % find negative peak at FCz

Peak = round(TIME(start+I)); 

% export mean ERN 
mint=Peak - 50; %export start in ms
maxt=Peak + 50; %export end in ms
ERPNAME='ERN'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',abs(mint),maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,abs(mint),maxt), '-v6');

%% E2) export Pe

BFS = RDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

% export mean Pe
mint=396-50; %export start in ms
maxt=396+50; %export end in ms
ERPNAME='Pe'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

