% Script to recreate the output from openEphys and plot example of raw data over
% multiple channels. Run this from the directory containing the data to be
% plotted (i.e. a collection of 32 .continuous files)

% Load in the channel map
chanmap_fn='C:\MATLABcode\code_for_spikes\DesktopClustering\chanMapPoly3_25s_correct.mat';
chanMap_struct=load(chanmap_fn);
chanMap=chanMap_struct.chanMap;
nchans=length(chanMap);  %' number of channels

%pick a group of eight channels - we will plot out data on these channels.
%Low numbers are further down the probe.
chans_chosen=4:11;
chans_to_plot=chanMap(chans_chosen); %the channels to plot - use a heatmap to work out which ones are important here

%% get the TTLs and select a time period to use

TTLs = returnTTLs_not_from_zero([pwd '\']);
% process footshock TTLs - these are in pairs (start, end) but we only want
% the one marking the onset:
fsts=TTLs.digital{2};
fsts(2:2:end)=[];

% Pull out sets of shocks delivered at difference frequencies - we can
% include these on the plot.
inds2=find( diff(fsts)>0.48 & diff(fsts) <0.51)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);  % 2Hz stims

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);  %0.25Hz stims


% Obreja protocol times 
obj_times=[800,1875];
fsts_Ob=fsts(fsts>obj_times(1)& fsts<obj_times(2));

% Start to read in data - first work out the format of the continuous file names
% in the folder.
[file_prefix, file_postfix]=returnOEPformat([pwd '\']);

% read in the first few samples of an example file, to check first
% timestamp and the sampling rate
ex_chan_name=[pwd '\' file_prefix 'CH' num2str(chans_to_plot(1)) file_postfix '.continuous'];
[data_ex, ts_ex,  info_ex] = load_open_ephys_data(ex_chan_name, 'Indices',1:10);
fs=info_ex.header.sampleRate;
first_ts=ts_ex(1);

%% Pick a time window containing a few footshocks, using the variables extracted above
time_window=[1700, 1703.5]; % in seconds
time_window_zeroed=time_window-first_ts;  %plot relative to a zero start point.
ind_window=round(time_window_zeroed*fs); %indices of samples in selected window
data_=[];  % for storing extracted data
fsts_inc=find(fsts>time_window(1) & fsts <time_window(2));  %TTLs included in selected window

%set up a figure for plotting
raw_data_fig=figure('Color', 'w', 'Units', 'Centimeters', 'Position', [1,1.5,25,19])

% now read in the data
for chan=1:length(chans_to_plot)
   this_chan=chans_to_plot(chan); %these vectors start with the lowermost channel, so plot from bottom up
   this_chan_x=chanMap_struct.xcoords(chans_chosen(chan)); % x-coord on probe for this chan
   fprintf('Extracting data for channel %d \n', chan);
%  chan_name=[pwd '\' file_prefix num2str(chan) file_postfix '.continuous']; % NEW FORMAT
   chan_name=[pwd '\' file_prefix 'CH' num2str(this_chan) file_postfix '.continuous'];
   
   if chan==1  % load in data. For first channel, read in timestamps, then don't bother
       [data_(chan,:), ts_,  ~] = load_open_ephys_data(chan_name, 'Indices',ind_window(1):ind_window(2));      
   else
       [data_(chan,:), ~,  ~] = load_open_ephys_data(chan_name, 'Indices',ind_window(1):ind_window(2)); 
   end
   
   %plot data on this channel
   subplot('Position', [0.065 0.05+(chan-1)*0.118 0.85 0.085]);
   plot(ts_, data_(chan,:))
   ylabel('\muV');
   
   if chan==1
    xlabel('Time (s)');
   end
   
   ylim([-55,55])
   hold on
   plot(repmat(fsts(fsts_inc), 1, 2)', [-55, 55], 'r');
   
   if this_chan_x<22 & this_chan_x >5;
       textcol=rgb('MediumSeaGreen');
   elseif this_chan_x>30
       textcol=rgb('Red');
   elseif this_chan_x==0
       textcol=rgb('MediumBlue');
   end
   
   txt=text(time_window(1)+0.1, 45, ['Chan. ' num2str(chans_chosen(chan))], 'FontSize', 9, 'Color', textcol);
   box off
   xlim(time_window)
end