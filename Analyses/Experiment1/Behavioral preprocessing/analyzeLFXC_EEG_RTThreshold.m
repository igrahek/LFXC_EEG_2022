%% Preprocessing script for the behavioral data of the LFXC_EEG experiment

clear;

runGLMs = 0;
runCalibTest = 0;
EVCanalysis = 0;
plotHistograms = 0;
plotScatter = 0;
fitLR_grid  = 0;
fitLR_fmin  = 0;
generateFigs = 0;

% pick a version to analyze
version_exp = 104;

% pick a running average window for efficacy and reward rate 
RunAvgWindow_Eff = 5;
RunAvgWindow_Rew = 5;
RunAvgWindow_RT = 5;
RunAvgWindow_Acc = 5;

%% Load in data files
files = dir(sprintf('../../../Data/Experiment1/Behavior/raw/%d/ST_*mat',version_exp));
fileNames = {files(:).name};

for curSub = 1:length(fileNames)
    file = [sprintf('../../../Data/Experiment1/Behavior/raw/%d/',version_exp),fileNames{curSub}];
    data.raw(curSub) = load(file);
end

% Initialize
data.all = [];
for curSub = 1:length(fileNames)
    subID = [];
    version = [];
    phase = [];
    block = [];
    trialNum = [];
    effLvl = [];
    rewLvl = [];
    accuracy = [];
    isRewarded = [];
    RT = [];
    AccRT = [];
    resp = [];
    cong = [];
    sliderEnd = [];
    sliderChanges = [];
    probeRt = [];
    probeType = [];
    feedbackOrder = [];
    efficacyDrift = [];
    RTthreshold = [];
    
    subID = str2num(data.raw(curSub).p.subID);  
    version = data.raw(curSub).p.version;
   
    phase = [ones(data.raw(curSub).p.numTrials.perform,1)];
    
    block = [sort(repmat([1:data.raw(curSub).p.numBlocks.perform],1,data.raw(curSub).p.numTrialsPerBlock.perform))]';
    
    trialNum = [repmat([1:data.raw(curSub).p.numTrialsPerBlock.perform*data.raw(curSub).p.numBlocks.perform],1,1)]';
      
    effLvl = [data.raw(curSub).results.perform.effLevel]';
    
    rewLvl = [data.raw(curSub).results.perform.rewLevel]';
    
    accuracy =  [double(data.raw(curSub).results.perform.xaccuracy)]'; %IG added double because there was an error on line 71 when putting in NaNs because accuracy was logical
    
    RT = [cell2mat(data.raw(curSub).results.perform.rt)]';
    
    AccRT = [cell2mat(data.raw(curSub).results.perform.rt)]';
    
    resp = [cell2mat(data.raw(curSub).results.perform.resp)]';
    
    isRewarded = [data.raw(curSub).results.perform.xisRewarded]';
    
    feedbackOrder = [data.raw(curSub).results.perform.feedbackOrder]';
    
    RTthreshold = [repmat(data.raw(curSub).p.rtThreshold,1,288)]';
       
    accuracy(resp == -1) = nan;
    RT(resp == -1) = nan;
    resp(resp == -1) = nan;
    

    congPerform = [];
    curSubIsCong = [];  % this needs to be cleared out on each loop
    curSubCong = zeros(size(data.raw(curSub).p.stimuli, 2),1);
    curSubWordResp = zeros(size(data.raw(curSub).p.stimuli, 2),1);
    
    for y = 1:size(data.raw(curSub).p.stimuli.perform, 2)
        
        congPerform(y) =  data.raw(curSub).p.stimuli.perform(1,y).IsCongruent;
        % get the target and the distractor information
        curSubTarget(y) = data.raw(curSub).p.stimuli.perform(1, y).CorrectAns;
        try
            if strcmp(data.raw(curSub).p.stimuli.perform(1,y).Text,'RED')
                wordCode = 1;
            elseif strcmp(data.raw(curSub).p.stimuli.perform(1,y).Text,'YELLOW')
                wordCode = 2;
            elseif strcmp(data.raw(curSub).p.stimuli.perform(1,y).Text,'GREEN')
                wordCode = 3;
            elseif strcmp(data.raw(curSub).p.stimuli.perform(1,y).Text,'BLUE')
                wordCode = 4;
            elseif strcmp(data.raw(curSub).p.stimuli.perform(1,y).Text,'XXXXX')
                wordCode = 5;
            end
            
            if wordCode == 5
                curSubWordResp(y) = nan;
            else
                curSubWordResp(y) = find(data.raw(curSub).p.colorOrder==wordCode);
            end
        catch
        end
        curSubTarget = curSubTarget';
        
    end    
    
    cong = [congPerform]';
    
    
    %% Add the learning probe data    
    try
        sliderEnd = data.raw(curSub).results.performanceProbe.sliderEndPosition';
        
        %sliderStart = [nan(1,length(data.raw(curSub).results.sliderStartPosition))]'; % the starting point of the slider
        
        probeRt= [nan(1,length(data.raw(curSub).results.performanceProbe.rt))]';
        
%         probeType = [nan(1,length(data.raw(curSub).results.performanceProbe.type))]';
        probeType = data.raw(curSub).results.performanceProbe.type';

        sliderEnd(probeType > 0) = ((sliderEnd(probeType > 0) - (data.raw(curSub).p.xCenter - 300)) / 600);
        
        
        
    catch
       
        
    end
    
    

    
    curSubData = [subID*ones(length(phase),1),...
        version*ones(length(phase),1),...
        phase,...
        block,...
        trialNum,...
        effLvl,...
        rewLvl,...
        accuracy,...
        curSubTarget,...
        curSubWordResp',... 
        RT,...
        AccRT,...
        resp,...
        isRewarded,...
        cong,...
        sliderEnd,...
        probeRt,...
        probeType,...
        feedbackOrder,...
        RTthreshold];
    

    
    data.all = [data.all; curSubData];
    
    
end


%% Table: Long-form data
% % Assuming congruency is 3-level (incong =0, cong=1, neutral=2):
allSubHdr = {'SubID','Vers','Phase','Block','Trial','EffLvl','RewLvl','Acc','TargetResp','DistractorResp','RT','AccRT','Resp','IsRewarded','Congruency','SliderEnd',...
    'ProbeRt','ProbeType','feedbackOrder','RTthreshold'}; %,'ChosenEff','ModeChosenEff','AltEff','ChosenTrialAmt','AltTrialAmt','ChoiceRT','ChoseRight','LeftEffOption',...
    %'RightEffOption','LeftTrialsOption','RightTrialsOption'
% 'TargetResp','RollingEffLvl','Fast','AccRT','IsRewarded','IsMiss','runavgRRgen','runavgRRcond','ZRT','AccZRT','RtThreshold'};
allSubData = array2table(data.all);
allSubData.Properties.VariableNames = allSubHdr;

allSubData.AccRT = allSubData.Acc .* allSubData.RT;
allSubData.AccRT(allSubData.AccRT == 0) = nan;
allSubData.zAccRT = AS_nanzscore(allSubData.AccRT);

allSubData.isMiss = isnan(allSubData.RT);

%%%% make efficacy categorical
allSubData.EffLvlCat = categorical(allSubData.EffLvl);

%%% nan non-probe trials
allSubData.SliderEnd(allSubData.SliderEnd == 0) = nan;
allSubData.ProbeType(allSubData.ProbeType == 0) = nan;

%%% make separate columns for probes
allSubData.EfficacyProbeResp = nan(height(allSubData),1);
allSubData.RewRateProbeResp = nan(height(allSubData),1);


for row = 1:height(allSubData)
    if allSubData.ProbeType(row) == 1
        
        allSubData.EfficacyProbeResp(row) = allSubData.SliderEnd(row);
        
    elseif allSubData.ProbeType(row) == 2
        
        allSubData.RewRateProbeResp(row) = allSubData.SliderEnd(row);
        
    end
end


subs = unique(allSubData.SubID);




%%% make feedback efficacy (t-1...t-10) columns

allSubData.Efficacy_T1 = nan(height(allSubData),1);
allSubData.Efficacy_T2 = nan(height(allSubData),1);
allSubData.Efficacy_T3 = nan(height(allSubData),1);
allSubData.Efficacy_T4 = nan(height(allSubData),1);
allSubData.Efficacy_T5 = nan(height(allSubData),1);
allSubData.Efficacy_T6 = nan(height(allSubData),1);
allSubData.Efficacy_T7 = nan(height(allSubData),1);
allSubData.Efficacy_T8 = nan(height(allSubData),1);
allSubData.Efficacy_T9 = nan(height(allSubData),1);
allSubData.Efficacy_T10 = nan(height(allSubData),1);

for currSubIndex = 1:length(subs)
    
    currSub = subs(currSubIndex);
    
    currEfficacy = allSubData.EffLvl(allSubData.SubID == currSub);
    
    allSubData.Efficacy_T1(allSubData.SubID == currSub) = [NaN; currEfficacy(1:end-1)];
    allSubData.Efficacy_T2(allSubData.SubID == currSub) = [NaN;NaN; currEfficacy(1:end-2)];
    allSubData.Efficacy_T3(allSubData.SubID == currSub) = [NaN;NaN;NaN; currEfficacy(1:end-3)];
    allSubData.Efficacy_T4(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN; currEfficacy(1:end-4)];
    allSubData.Efficacy_T5(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-5)];
    allSubData.Efficacy_T6(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-6)];
    allSubData.Efficacy_T7(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-7)];
    allSubData.Efficacy_T8(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-8)];
    allSubData.Efficacy_T9(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-9)];
    allSubData.Efficacy_T10(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currEfficacy(1:end-10)];
    
end





%%% make isRewarded lag columns for lagged regression
allSubData.IsRewarded_T1 = nan(height(allSubData),1);
allSubData.IsRewarded_T2 = nan(height(allSubData),1);
allSubData.IsRewarded_T3 = nan(height(allSubData),1);
allSubData.IsRewarded_T4 = nan(height(allSubData),1);
allSubData.IsRewarded_T5 = nan(height(allSubData),1);
allSubData.IsRewarded_T6 = nan(height(allSubData),1);
allSubData.IsRewarded_T7 = nan(height(allSubData),1);
allSubData.IsRewarded_T8 = nan(height(allSubData),1);
allSubData.IsRewarded_T9 = nan(height(allSubData),1);
allSubData.IsRewarded_T10 = nan(height(allSubData),1);

for currSubIndex = 1:length(subs)
    
    currSub = subs(currSubIndex);
    
    
    
    currIsRew = allSubData.IsRewarded(allSubData.SubID == currSub);
    
    allSubData.IsRewarded_T1(allSubData.SubID == currSub) = [NaN; currIsRew(1:end-1)];
    allSubData.IsRewarded_T2(allSubData.SubID == currSub) = [NaN;NaN; currIsRew(1:end-2)];
    allSubData.IsRewarded_T3(allSubData.SubID == currSub) = [NaN;NaN;NaN; currIsRew(1:end-3)];
    allSubData.IsRewarded_T4(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN; currIsRew(1:end-4)];
    allSubData.IsRewarded_T5(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-5)];
    allSubData.IsRewarded_T6(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-6)];
    allSubData.IsRewarded_T7(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-7)];
    allSubData.IsRewarded_T8(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-8)];
    allSubData.IsRewarded_T9(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-9)];
    allSubData.IsRewarded_T10(allSubData.SubID == currSub) = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN; currIsRew(1:end-10)];
    
end


%% make nearest neighbor efficacy and reward probe ratings column

% Efficacy probe
allSubData.ReportedEfficacyNN = nan(height(allSubData),1);
allSubData.ReportedEfficacyLin = nan(height(allSubData),1);


for currSubIndex = 1:length(subs)
    
    currSub = subs(currSubIndex);
    
  
        allSubData.ReportedEfficacyNN(allSubData.SubID == currSub)...
            = fillmissing(allSubData.EfficacyProbeResp(allSubData.SubID == currSub),'nearest');
        
        preProbeTrials = min(find(~isnan(allSubData.EfficacyProbeResp(allSubData.SubID == currSub))));
        allSubData.ReportedEfficacyNN(allSubData.SubID == currSub & allSubData.Trial < preProbeTrials) = NaN; 
              
        allSubData.ReportedEfficacyLin(allSubData.SubID == currSub)...
            = fillmissing(allSubData.EfficacyProbeResp(allSubData.SubID == currSub),'linear');
        
        preProbeTrials = min(find(~isnan(allSubData.EfficacyProbeResp(allSubData.SubID == currSub))));
        allSubData.ReportedEfficacyLin(allSubData.SubID == currSub & allSubData.Trial < preProbeTrials) = NaN; 
        

end

    

% Reward rate probe

allSubData.ReportedRewRateNN = nan(height(allSubData),1);
allSubData.ReportedRewRateLin = nan(height(allSubData),1);


for currSubIndex = 1:length(subs)
    
    currSub = subs(currSubIndex);
    
  
        allSubData.ReportedRewRateNN(allSubData.SubID == currSub)...
            = fillmissing(allSubData.RewRateProbeResp(allSubData.SubID == currSub),'nearest');
        
        preProbeTrials = min(find(~isnan(allSubData.RewRateProbeResp(allSubData.SubID == currSub))));
        allSubData.ReportedRewRateNN(allSubData.SubID == currSub & allSubData.Trial < preProbeTrials) = NaN; 
              
        allSubData.ReportedRewRateLin(allSubData.SubID == currSub)...
            = fillmissing(allSubData.RewRateProbeResp(allSubData.SubID == currSub),'linear');
        
        preProbeTrials = min(find(~isnan(allSubData.RewRateProbeResp(allSubData.SubID == currSub))));
        allSubData.ReportedRewRateLin(allSubData.SubID == currSub & allSubData.Trial < preProbeTrials) = NaN; 
        

end






%% adding in mean eff ratings and running averages of Acc & RT to table


subs = unique(allSubData.SubID);

allSubData.AvgSubjEffLvl = nan(height(allSubData),1);
allSubData.runAvgRewRate = nan(height(allSubData),1);
allSubData.runAvgRewPerformCong = nan(height(allSubData),1);
allSubData.runAvgEfficacy = nan(height(allSubData),1);

allSubData.runAvgAccRT = nan(height(allSubData),1);
allSubData.runAvgAcc = nan(height(allSubData),1);

for currSub = subs(1):subs(end)

        
        %% reward rate checks
        filterKernel = RunAvgWindow_Rew;
        filterEnvelope = [zeros(1,1),ones(1,filterKernel)]/(filterKernel);
        
        curSubRunavgRRgen = filter(filterEnvelope,1,allSubData(allSubData.SubID == currSub,:).IsRewarded);       
        allSubData(allSubData.SubID == currSub,:).runAvgRewRate = curSubRunavgRRgen;
                  
        %% running avg efficacy
        filterKernel2 = RunAvgWindow_Eff;
        filterEnvelope2 = [zeros(1,1),ones(1,filterKernel2)]/(filterKernel2);
        
        curSubRunAvgRewPerformCongGen = filter(filterEnvelope2,1,allSubData(allSubData.SubID == currSub,:).EffLvl);
        allSubData(allSubData.SubID == currSub,:).runAvgEfficacy = curSubRunAvgRewPerformCongGen;   
        
        %% running avg RT
        filterKernel3 = RunAvgWindow_RT;
        filterEnvelope3 = [zeros(1,1),ones(1,filterKernel3)]/(filterKernel3);
        
        curSubRunAvgRT = filter(filterEnvelope3,1,allSubData(allSubData.SubID == currSub,:).AccRT);
        allSubData(allSubData.SubID == currSub,:).runAvgAccRT = curSubRunAvgRT; 
        
        %% running avg RT
        filterKernel4 = RunAvgWindow_RT;
        filterEnvelope4 = [zeros(1,1),ones(1,filterKernel4)]/(filterKernel4);
        
        curSubRunAvgAcc = filter(filterEnvelope4,1,allSubData(allSubData.SubID == currSub,:).Acc);
        allSubData(allSubData.SubID == currSub,:).runAvgAcc = curSubRunAvgAcc; 
   
        
end
    

%% adding in current subj eff ratings to table -- tracks ratings over time in learning

allSubData.CurrSubjEffLvl = nan(height(allSubData),1);
CurrSubjEffLvl = [];

CurrSubjEffLvl = nan;
curBlock = 1;

subEffRatingMatrix = [0.1; 0.4; 0.6; 0.9];

for row = 1:height(allSubData)
    if allSubData.Phase(row) == 1 || (allSubData.Phase(row) == 2 && allSubData.Block(row) == 1 && allSubData.Trial(row) == 1)
        if curBlock ~= allSubData.Block(row)
            subEffRatingMatrix(subEffRatingMatrix(:,1) == allSubData.EffLvl(row-1),2) = CurrSubjEffLvl;
            
            curBlock = allSubData.Block(row);
            CurrSubjEffLvl = nan;
        end
        
        if allSubData.ProbeType(row) == 1
            if ~isnan(allSubData.CurrSubjEffLvl(row-1))
                CurrSubjEffLvl = (allSubData.SliderEnd(row) + CurrSubjEffLvl)/2;  %%look at this for repeating sliderEnd column
                
            else
                CurrSubjEffLvl = allSubData.SliderEnd(row);
            end
        end
        
        allSubData.CurrSubjEffLvl(row) = CurrSubjEffLvl;
    end
    
    if allSubData.Phase(row) == 2 && allSubData.ProbeType(row) == 1
        subEffRatingMatrix(subEffRatingMatrix(:,1) == allSubData.EffLvl(row),3) = allSubData.SliderEnd(row);
    end
    
    if allSubData.Phase(row) == 3
        allSubData.CurrSubjEffLvl(row) = subEffRatingMatrix(subEffRatingMatrix(:,1) == allSubData.EffLvl(row),2);
    end
    
end

%% Wirte down the size of the windows used for the calculation of the running averages

allSubData.RunAvgWindow_Eff = repmat(RunAvgWindow_Eff,height(allSubData),1);
allSubData.RunAvgWindow_Rew = repmat(RunAvgWindow_Rew,height(allSubData),1);


%% Saving data
writetable(allSubData, sprintf('../../Data/Behavior/preprocessed/data%d.csv',version_exp))
%d/ST_*mat
