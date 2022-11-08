%% Takes a set of OEP files, cleans artefacts, saves data as MATLAB files.
% Run from directory containing the files to be edited (i.e. the 32
% .continuous file)

TTLs_all = returnTTLs_not_from_zero([pwd '\']); %returns the TTLs recorded by OEP, with timestamps matching .continuous files
TTLs=TTLs_all.digital{2};  % The number in brackets is the channel the TTL on interest is associated with
TTLs(2:2:end)=[];  %theres a TTL for 'on' and one for 'off' - we only need the 'on'
file_prefix='100_';  %amend this to match the format of your TTL files
file_postfix='_2';

%% OPTIONAL: use this code snippet to look at an example file and work out 
%  what time window around TTL needs to be edited

% fn_=[file_prefix '32' , file_postfix '.continuous'];
% [data_ex, ts_ex,  ~] = load_open_ephys_data([fn_]);
% figure
% plot(ts_ex, data_ex)
% hold on
% plot([TTLs, TTLs], [-200, 200], 'r')
%% Now set the window:
%ms after TTL to edit in format [before after]. 
%Values in this window around the TTL will be set back to the mean:
removal_window=[1,3];  


%output files will always be put into a folder called 'artefact_free'
if ~exist('artefact_free', 'dir')
   mkdir('artefact_free')
end

%% Go through each file and edit around each TTL, then save the new version of the file
for ch=1:32
    fn_raw=[file_prefix 'CH' num2str(ch), file_postfix '.continuous'];
    %    fn_raw=[file_prefix '' num2str(ch), file_postfix '.continuous']; %for new format files
    edited_file=removeArtefact(fn_raw, TTLs, removal_window);
    
    fn_edited=[file_prefix 'CH' num2str(ch) file_postfix '.mat'];
    save([pwd '\artefact_free\' fn_edited],'edited_file' ,'-v7.3');
end

% Save metadata in same folder
[~, ts_raw,  info] = load_open_ephys_data([fn_raw]);
OEPinfo.info=info;
OEPinfo.t_first=ts_raw(1);
OEPinfo.aretfact_removal_win=removal_window;
save([pwd '\artefact_free\OEPInfo.mat'],'OEPinfo' );