
%% fix ocular correction file and export as brain vision

PATH = '~/Dropbox (Brown)/ShenhavLab/EEG_ressources/Experiments/LFXC_EEG/';
SOURCEFILESOC = dir(strcat(PATH, 'Data/FixFiles/*.vhdr'));

eeglab;
%%
load(sprintf('%sData/Export/chanlocs.mat', PATH));
for n = 1:65
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end




%%

for s =1:numel(SOURCEFILESOC)%
%% load the suspect
EEG = pop_loadbv(sprintf('%sData/FixFiles/',PATH),sprintf('%s',SOURCEFILESOC(s).name));
fn   = SOURCEFILESOC(s).name(1:4);

%%
if str2double(fn)== 1042

    chan1=LO2;
    chan2=IO2;

else
    
    chan1=LO1;
    chan2=IO1;
    
end
%% insert fix here
           
TEMP =EEG.data(chan1,:);

EEG.data(chan1,:)=EEG.data(chan2,:);
EEG.data(chan2,:)=TEMP;
 %% export as vhdr
 % this is hpw it's done!
 if length(SOURCEFILESOC(s).name)>9
     
     DATTYPE='Data/OC/';
 else
     DATTYPE='Data/raw/';
 end
 
 pop_writebva(EEG,sprintf('%s%s%s', PATH, DATTYPE, SOURCEFILESOC(s).name));

end