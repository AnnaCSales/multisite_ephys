function [data_edited] = removeArtefact(data_fn, TTLs, window)
%Removes artefact associated with TTLs. Takes filename for an OEP data file
%(i.e. one continuous file) and a list of TTL times. 
% Also expects a variable called 'window'. 'window' is
% the number of ms around the TTL to be edited, and is specifed as a two
% element vector :
% window=[time_before_TTL, time_after_TTL] (both positive numbers)
% For every TTL, the function simply sets values in the window to the mean
% value of the neural data in the entire file. This will prevent TTL
% artefacts from being confused with threshold crossing events/spikes.
% Anna Sales, UoB, November 2020


[data_raw, ts_raw,  info] = load_open_ephys_data([data_fn]);

if size(TTLs,2)>1
    TTLs=TTLs';  %want a column vector
end

mean_val=mean(data_raw);
window_s=window*0.001;
fprintf('\n Editing file %s \n', data_fn)
data_edited=data_raw;
for h=1:length(TTLs)
    if ~mod(h, 25)
        fprintf('\nEditing TTL # %d of %d', h, length(TTLs))
    end
    this_TTL=TTLs(h);
    in_window=find(ts_raw>this_TTL-window_s(1) & ts_raw<this_TTL+window_s(2)); 
    data_edited(in_window)=mean_val;
       
end


end