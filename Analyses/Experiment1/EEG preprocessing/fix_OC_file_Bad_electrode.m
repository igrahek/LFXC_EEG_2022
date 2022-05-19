
%% fix ocular correction file and export as brain vision


PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/';
SOURCEFILESOC = dir(strcat(PATH, 'Data/OC/*.matrix'));

s=27;
%%
load(sprintf('%sData/Export/chanlocs.mat', PATH));
for n = 1:65
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end


%% load the suspect
EEG = pop_loadbv(sprintf('%sData/OC/',PATH),sprintf('%s',SOURCEFILESOC(s).name));

%% nb1: fix C6
           x=EEG.data(O2,:);
           y = mean(EEG.data([PO4 POz Oz],:));   
           %H =[x',y'];
           % we only do that though when the difference between the
           % interpolation and the actual measurement exceeds 40 microvolt
           %z= x-y;
           %fdif= find(abs(z>40));
           EEG.data(O2,:)=y; 
% %% now fix T8
%            x=EEG.data(T8,:);
%            y = mean(EEG.data([FT8 C6 TP8],:));   
%            H =[x',y'];
%            % we only do that though when the difference between the
%            % interpolation and the actual measurement exceeds 40 microvolt
%            %z= x-y;
%            %fdif= find(abs(z>40));
%            EEG.data(T8,:)=y;
%            
 %% export as vhdr
 % this is hpw it's done!
 pop_writebva(EEG,sprintf('%sData/OC/%s', PATH,SOURCEFILESOC(s).name));
 
 
 %% fix 1081 bad electrodes
 
 PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/';
SOURCEFILESOC = dir(strcat(PATH, 'Data/raw/*.vhdr'));

s=40;
%%
load(sprintf('%sData/Export/chanlocs.mat', PATH));
for n = 1:65
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end


%% load the suspect
EEG = pop_loadbv(sprintf('%sData/raw/',PATH),sprintf('%s',SOURCEFILESOC(s).name));
%% find bad channel(s)
%     EEG.data(end+1,:) = 0;
%     EEG.nbchan = size(EEG.data,1);
%     EEG.chanlocs(end+1).labels = 'Cz';
% %%
%   
%     %load chanloc information
%     %EEG = pop_chanedit(EEG, 'lookup',sprintf('%sOccular_Correction/werfen_BESA.elp',PATH));    
%     EEG=pop_chanedit(EEG, 'lookup',sprintf('/Applications/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
%     %re-reference to average
%     EEG = pop_reref( EEG, []);
% 
% chStd = std(EEG.data')'; % --> TP9 has super huge std, so we need to fix this.
%% let's reload the raw data
EEG = pop_loadbv(sprintf('%sData/raw/',PATH),sprintf('%s',SOURCEFILESOC(s).name));
%% nb1: fix C6
%             x1=EEG.data(TP9,:);
%             x2=EEG.data(TP7,:);
%            y = mean(EEG.data([T7  P7],:));   
%             H =[x1',x2',y', EEG.data(TP8,:)', EEG.data(TP10,:)'];
%            % we only do that though when the difference between the
%            % interpolation and the actual measurement exceeds 40 microvolt
%            %z= x-y;
%            %fdif= find(abs(z>40));
%            %%
%            EEG.data(TP9 ,:)=y; 

%% due to the excentricity of the electrodes that are bad, I will just replace the left ones with the right ones.
% we will use some asymmetry information, but that's better than fucking up
% our average reference completely.
EEG.data(TP9 ,:)=EEG.data(TP10 ,:);
EEG.data(TP7 ,:)=EEG.data(TP8 ,:);

%            x=EEG.data(T8,:);
%            y = mean(EEG.data([FT8 C6 TP8],:));   
%            H =[x',y'];
%            % we only do that though when the difference between the
%            % interpolation and the actual measurement exceeds 40 microvolt
%            %z= x-y;
%            %fdif= find(abs(z>40));
%            EEG.data(T8,:)=y;
%            
 %% export as vhdr
 % this is hpw it's done!
 pop_writebva(EEG,sprintf('%sData/raw/%s', PATH,SOURCEFILESOC(s).name));

