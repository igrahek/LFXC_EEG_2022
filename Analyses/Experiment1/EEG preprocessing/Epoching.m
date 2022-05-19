%% Making BFS for FXC EEG

% 1) Cue-locked -200: 1500 ms
% 2) Stim-locked: -200:800 ms
% 3) R-locked: -200: 800 ms
% 4) rew FB locked: -200: 800 ms
% 5) eff FB locked: -200: 800 ms

%% Info for segmentation analyses:

PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/';
addpath('~/Dropbox (Brown)/ShenhavLab/EEG_ressources/EEGfunctions/')

%% about the data:


%             p.EEGpars.trigs.cueCode = 10;
%             p.EEGpars.trigs.stimCode = 20;
%             p.EEGpars.trigs.respCode = 30;
%             p.EEGpars.trigs.rewFBCode = 40;
%             p.EEGpars.trigs.effFBCode = 50;

Trigs =[10,20, 30, 40, 50];
Mks = cellstr(num2str(Trigs(:)));
Mks=strcat({'S'}, {' '}, Mks); % adding S to Markers for segmentation


alltrigs= Mks;%
allIL= [-0.2, 1.5; -0.2, 0.8; -0.2, 0.8; -0.2, 0.8; -0.2, 0.8];%


nelecs= 65;


%% approach: Get all condition triggers
% just NAN the artifact trials and paste all subject files together.
%% run only when needed, adds EEGLab

% addpath('N:/Software/eeglab13_5_4b')
%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
SOURCEFILES = dir(strcat(PATH, 'Data/clean/*clean.set')); %all cleaned files
%SOURCEFILES = SOURCEFILES([1:14 16:end])
SUBJECTS = 1:numel(SOURCEFILES);

%% for response-locked segementation:
load(sprintf('%sData/Export/FXCallSubDataTable.mat', PATH))
if length(unique(FXCallSubDataTable.SubID)) < numel(SOURCEFILES)
    cd('~/Dropbox (Brown)/ShenhavLab/experiments/LFXC_EEG_2Cues/analysis')
    analyzeLFXC_EEG_new
    load(sprintf('%sData/Export/FXCallSubDataTable.mat', PATH))
    cd([PATH, 'Analyses/Matlab/'])
end
%%
% adjust table to match cleaned data:
if length(unique(FXCallSubDataTable.SubID))> numel(SOURCEFILES)
    FXCallSubDataTable= FXCallSubDataTable(FXCallSubDataTable.SubID<(str2double(SOURCEFILES(end).name(1:4))+1),:);
end
%%
cf = zeros(numel(SOURCEFILES),2);
cf(:,1)= unique(FXCallSubDataTable.SubID);
for iii = 1: numel(SOURCEFILES)
    cf(iii,2)= str2double(SOURCEFILES(iii).name(1:4));
    
end
%%
if mean(cf(:,1)==cf(:,2))==1
    
    
SubAllLatencies= [];
    allsegs= {'c'; 's'; 'r'; 'fr'; 'fe'}; %
    contmat=zeros(40,3, 4); % check number of trials found per participant (before & after AR)
    %% loop through different segmentations to save memory
    for nsegs = 1:length(allsegs)
        
        if strcmp(allsegs{nsegs}, 'c')
            DAT=nan(nelecs,850,height(FXCallSubDataTable)); % preallocate memory and preset datamatrix (larger for cue than others)
            %     elseif strcmp(allsegs{nsegs}, 'r')
            %
            %         DAT=nan(nelecs,500,height(FXCallSubDataTable)); % preallocate memory and preset datamatrix (larger for cue than others)
            %
            %     elseif strcmp(allsegs{nsegs}, 'fr')
            %         DAT=nan(nelecs,500,height(FXCallSubDataTable));
            
        else
            DAT=nan(nelecs,500,height(FXCallSubDataTable));
        end
        
        Mks ={alltrigs{nsegs}};
        
        %Mks=strcat({'S'}, {' '}, Mks); % adding S to Markers for segmentation;
        %seems to be a bva thing...
        
        [subIDs,uniqueRows] = unique(FXCallSubDataTable.SubID); % get first line for each sub to start putting data into
        
        
        for s =1:numel(SOURCEFILES)
            %%
            vpn=SOURCEFILES(s).name(1:4); % read out of dataset name
            s_id = str2double(SOURCEFILES(s).name(1:4)); % reads in subject ID and converts it to double (does not work for '01')
            fprintf('processing participant number %d\n', s_id)
            contmat(s,1,nsegs)=s_id;
            pcnt=uniqueRows(s); % get first line for each sub to start putting data into
            lns = sum(FXCallSubDataTable.SubID==s_id);
            %%
            if lns>0
                % load fixed trigger file
                
                % load dataset
                EEG = pop_loadset('filename',sprintf('%s_clean.set',vpn),'filepath',sprintf('%sData/clean/',PATH));
                %%
                if strcmp(allsegs{nsegs}, 'r') % correct trigger timing for response-locked segmentation
                    
                    
                    
                            %%
        % get rig of weird empty marker:
        if strcmp(EEG.event(1).type, 'empty')
           EEG.event = EEG.event(2:end); 
        end
        
        vpdouble = str2double(vpn);
        % if response is nan, add 1000/2
        subRTs = FXCallSubDataTable.RT(FXCallSubDataTable.SubID==vpdouble);
        subLatencies= nan(length(subRTs),1);
        subLatenciesCHECK= nan(length(subRTs),1);
        subRTs(isnan(subRTs)) = 1000;
        lcnt=1;
        fprintf('%d trials\n', (size(EEG.event,2)-1)/5)
        % now we basically go through all trials. Each trial has 5
        % triggers, so we go in steps of 5 to the end of our events.
        for eventnumber=2:5:size(EEG.event,2) % the first trigger is new segment, so we skip that 
            %This is just to make sure, we're not messing anything up
         if size(EEG.event,2)>(length(subRTs)*5+1)
             fprintf('Trigger number does not match number of trials!!!\n')
             break
         else
             % we need to get the latencies for the stimulus trigger, which
             % is the second trigger to replace the response trigger, which
             % is the 3rd one
           subLatencies(lcnt) = EEG.event(eventnumber+2).latency-EEG.event(eventnumber+1).latency;  % get original trigger latency relative to stim trigger
           EEG.event(eventnumber+2).latency =  EEG.event(eventnumber+1).latency + round(subRTs(lcnt)/2); % replace trigger latency by latency of previous trigger (stim) + half of RT (bc data points, sampling every 2ms)
           subLatenciesCHECK(lcnt)= EEG.event(eventnumber+2).latency-EEG.event(eventnumber+1).latency; % get new trigger latency relative to stim trigger 
         end
          lcnt=lcnt+1;  
        end 
        CHECK= [subLatencies, subLatenciesCHECK]; % this is to compare original and new trigger latencies
        CHECK2 = CHECK*2; % in RT format
        CHECK2(:,3) = CHECK2(:,1)- CHECK2(:,2); % and their difference
        SubAllLatencies=[SubAllLatencies;CHECK2]; % across all subjects.
                    
                    
                end
                %             if strcmp(allsegs{nsegs}, 'r') % do this only for r-locked data bc something seems weird...
                %             triggersFixed =readtable(sprintf('%s/TriggersCleaned/%s_triggersFixed.csv',PATH, vpn));
                %             % change triggercodes to fixed trigger codes
                %             for eventnum = 1:length(EEG.event)
                %
                %                 EEG.event(eventnum).type = triggersFixed.eventCodeNew(eventnum);
                %
                %             end
                %             clear triggersFixed;
                %             %EEGC=EEG; % save copy in temp memory
                %
                %             end
                %% get stim-locked data
                fprintf('%s-locked segmentation of participant number %d\n', allsegs{nsegs}, s_id)
                IL =allIL(nsegs,:); % 200 ms pre- stimulus + entire trial length
                % number of data points in segment
                
                [EEG, indices]=pop_epoch(EEG, Mks,IL);
                
                %%% want common baseline for everything? %%% think about it!!!
                EEG = pop_rmbase( EEG, [-200 0]); % Baseline correction relative to stimulus onset
                DP = size(EEG.data,2);
                
                % select all trials based on object onset markers
                % no further segmentation required --> matching with logfiles
                contmat(s,2,nsegs)=EEG.trials;
                %%
                % AR: Ansatz: Trials werden nicht gel?scht, sondern nur die
                % unzul?ssigen identifiziert und dann sp?ter nicht mit ausgelesen.
                % (Irej: index of rejected trials)
                %% I don't do the artifact rejection on the EOG channels
                [EEGn, Irej] = pop_eegthresh(EEG,1,[1:nelecs] ,-150,150, IL(1),IL(2)-0.002,1,0); % artifact rejection based on all but the ocular electrodes
                
                %% %   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);
                
                [rej, rejE] = rejtrend(EEG.data(1:nelecs, :,:),DP,50,0.3,1); % I don't do this for EOG either
                
                frej= find(rej==1);
                Out=unique(sort([Irej,frej]));
                if isempty(Out)
                    remout=0;
                else
                    remout = length(Out);
                end
                fprintf('%d trials removed alltogether\n', remout);
                
                %             if strcmp(allsegs{nsegs}, 'r')
                %                 ndif = length(FXCallSubDataTable.RT( FXCallSubDataTable.SubID==s_id & ~isnan(FXCallSubDataTable.Resp))) - size(EEG.data,3);
                %             else
                %                 ndif = 600- size(EEG.data,3);
                %             end
                %%
                EEG.data(:,:,Out) = nan; % invalidate data with artifacts
                DAT(:,:,pcnt :pcnt+size(EEG.data,3)-1)=EEG.data(:,:,1:size(EEG.data,3));% excludes practise trials
                contmat(s,3,nsegs)=(EEG.trials-size(Out,2));
                
            else
                fprintf('Subject not in table! Run analyzeFXC first!\n')
            end
        end
        
        
        fprintf('%d participants in dataset... saving %s locked data\n', s, allsegs{nsegs})
        if strcmp(allsegs{nsegs}, 'c')
            CDAT = DAT;
            save(sprintf('%sData/Export/CDAT.mat', PATH), 'CDAT', '-v7.3');
            clear CDAT;
        elseif strcmp(allsegs{nsegs}, 's')
            SDAT = DAT;
            save(sprintf('%sData/Export/SDAT.mat', PATH), 'SDAT', '-v7.3');
            clear SDAT;
        elseif strcmp(allsegs{nsegs}, 'r')
            RDAT=DAT;
            save(sprintf('%sData/Export/RDAT.mat', PATH), 'RDAT', '-v7.3');
            clear RDAT;
        elseif strcmp(allsegs{nsegs}, 'fr')
            FRDAT=DAT;
            save(sprintf('%sData/Export/FRDAT.mat', PATH), 'FRDAT', '-v7.3');
            clear FRDAT;
        else
            FEDAT=DAT;
            save(sprintf('%sData/Export/FEDAT.mat', PATH), 'FEDAT', '-v7.3');
            clear FEDAT;
        end
        
        
    end
    %
    % chanlocs = EEG.chanlocs;
    % save(sprintf('%sExport/chanlocs.mat', PATH), 'chanlocs', '-v7.3');
    
    save(sprintf('%sData/Export/contmat.mat', PATH), 'contmat', '-v7.3');
    
else
    fprintf('Subject numbers in behavioral and EEG data do not match!')
    
end