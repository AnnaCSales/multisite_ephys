% Script to estimate conduction velocities for specified units, as in
% Figure 5. The details of which channels to use for which unit, and
% whether to track trough or peak are in the powerpoint file
% 'Unit_feature_tracking' in the Figure 5 folder.

%% Load in the spikeStruct, get TTLs etc

% NB for 8th Dec, good fibres are 1-4, for 9th Dec, 3 and 4

 datapath=pwd; %dont forget backslash

%for 8th Dec
% Path to ADC file with pedal on
ADC_fn=[pwd '\100_ADC2_3.continuous'];   %for 8th Dec
% Path to continuous file with ECG recording.
ECG_fn=[pwd '\100_CH34_3.continuous'];

% ADC_fn=[pwd '\113_ADC2.continuous'];   %for 9th Dec
% % Path to continuous file with ECG recording.
%  ECG_fn=[pwd '\113_CH34.continuous'];

% distance from RF to recording point in cm
 rf_dist=4.7; %for 8th Dec
% rf_dist=5.0  % for 9th Dec

%% Times for digitimer TTL
TTLchan1=2;  %this is the channel to take TTLs from
fsts=spikeStruct.TTLs.digital{TTLchan1};
fsts(2:2:end)=[];

% fsts=spikeStruct.rescuedTTLs;
 %fsts=fsts-2400; %for edited data times
    
inds2=find( diff(fsts)>0.4 & diff(fsts) <.6)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);

figure;plot(fsts, 1.2*ones(numel(fsts)), 'rX');

% Obreja times (exclude stims for thresholding)
obj_times=[800,1875];
fsts_Ob=fsts(fsts>obj_times(1)& fsts<obj_times(2));

% Low frequency stims, preceeding 2Hz (for baseline CV calculation)

% inds_slow=find( diff(fsts)>3.9 )+1;
inds_slow=find( diff(fsts)>7.9 )+1;
inds_slow=[inds_slow(1)-1; inds_slow];
fsts_slow=fsts(inds_slow);
% find 2Hz period, exclude anything that happened after that.
first_2Hz=fsts_2Hz(1); 
fsts_slow(fsts_slow>=first_2Hz)=[];
fsts_slow(fsts_slow<850)=[];

nlabels=length(spikeStruct.TTLs.manual.TTL_labels);
%% Times of pedal down
pedal_ts=pedalOnOffs(ADC_fn);
%% Pull out useful things for plotting from the spikeStruct:

nclusts=spikeStruct.nclusts;  %number of clusters
fs=spikeStruct.sample_rate;   % sampling rate
newWFs=spikeStruct.allchanWFs;  %waveforms across channels, for each cluster
c_channel=spikeStruct.c_channel;  %centre channel for the cluster
av_waveform=spikeStruct.av_waveform;  %av waveform on the centre channel
plot_pos=spikeStruct.plot_pos;     % depth of each cluster
bl_start=spikeStruct.baseline_st;  %baseline info
bl_end=spikeStruct.baseline_end;
min_t=spikeStruct.timeRange(1);    %time range of recording
max_t=spikeStruct.timeRange(2);
sq_=ceil(nclusts^0.5) ; %for calculating the number of subplots required

%set up some labels for plots below.
for pos=1:1:length(plot_pos)
    unit_test(pos)=find(plot_pos==pos); %The unit that is in the pos-th position on the plot
    chan_=spikeStruct.c_channel(unit_test(pos));   
    tt=['Clu ', int2str(unit_test(pos)), '  c chan= ', int2str(spikeStruct.c_channel(unit_test(pos)))];
    ticklabs{pos}=tt;
end




%% Now plot waveforms as in Fig 5 and estimate CVs, using selected channels only

%set up parameters:
unit=6;  %which cluster are we analysing
plotparams.max_min=-1; %1 to track peaks, -1 to track troughs
plotparams.channels=[5, 7:11,12, 13,15,18]; %channels to include in plot.

%set this up so that we don't consider anything after the 2Hz events have
%started:
plotparams.events=fsts_2Hz(1);%fsts_2Hz; %events to exclude - send empty vector if none
plotparams.exclusionWindow=2280; %window after TTL in which to ignore spikes, in seconds.

%other parameters (don't change)
plotparams.fn='dataALL.bin'; %filename of raw data, with path if not in pwd
plotparams.unit=unit; %the unit to plot
plotparams.skipCloseChans=0;
plotparams.doPlot=1;
plotparams.twin=[];


[norm_xv_mean, norm_yv_mean, v_mag_mean] = mean_spk_CV_and_dir_select_chans(spikeStruct, plotparams)
fprintf('\n\nAction potential vector from mean spike = %2.2f i + %2.2f j m/s', norm_xv_mean, norm_yv_mean)