% This is a which uses raw data to produce heatmaps 
% showing spiking responses at constant latency. The idea is to explore 
% processing options that could be implemented in real time on OEP. This would 
% help to identify channels of interest during an experiment (i.e. those
% that are most likely to show C fibre activity. This version uses
% pre-recorded data in .continuous format to test different approaches
% NB there is no filtering in this script - if using raw data, it will need
% to be bandpass filtered first.

% Run this script from the folder containing the 32 .continuous files you
% want to visualise.

% Dependencies
% - a channel map file, in kilosort format, 
% - returnTTLs_not_from_zero (helper function to return the times of TTLs on each
% digital input.)
% -returnOEPformat - helper function to get pre and postfixes of OEP
% filenames
% - OEP utilities for reading in OEP data
% Anna Sales, April 2021

%% Specify path to prerecorded data
datapath=[pwd '/Data/021220_rat2/021220/rec2_cFibres/'];

% specify a path to a channel map file (script expects this in kilosort
% format, i.e a struct with a field called 'chanMap', which lists OEP channels
% from bottom of probe upwards
chanmap_fn=[pwd '/Code/Clustering/chanMapPoly3_25s_correct.mat'];
chanMap_struct=load(chanmap_fn);
chanMap=chanMap_struct.chanMap;
nchans=length(chanMap);  %' number of channels

%% Retrieve the TTLs marking electrical stimulation (footshock timestamps - 'fsts')

% either extract from the raw data:
TTLs = returnTTLs_not_from_zero(datapath); 
% process footshock TTLs - these are in pairs (start, end) but we only want
% the one marking the onset: NB for 091220 recording need to load TTLs manually from
% 'allTTLs' in folder, and comment out the two lines below.
 fsts=TTLs.digital{2};
 fsts(2:2:end)=[];

% Pull out sets of shocks delivered at different frequencues
inds2=find( diff(fsts)>0.48 & diff(fsts) <0.51)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);  % 2Hz stims

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);  %0.25Hz stims
%% Do a quick plot of TTLs to figure out where the start of
%  the Obreja protocol is - need to find out the ind of the first 0.25Hz stim
%  that's actually during Obreja and not one of the ones used for
%  thresholding. The first 0.25 TTL in Obreja should be in the middle of a
%  long run of other TTLs, at around 800-1200s.
figure;
slow_TTLS=plot(fsts_025Hz, ones(numel(fsts_025Hz)), 'X');
hold on;
all_TTLS=plot(fsts, 1.2*ones(numel(fsts)), 'rX');
legend({'0.25Hz TTLs' 'all TTLs'})
ylim([0.9, 1.3])
%% change the numbers in next line 
Ob_start=1200; %rough time point after which Obreja protocol starts

inds_of_interest=find(fsts_025Hz>Ob_start);  
first_real_025_ind=inds_of_interest(1);  %which one of the 0.25Hz events is the start of the protocol

%% Read in data (NB don't run this if windowed data has been previously saved!

% work out the format of the continuous file names
% in the folder.
[file_prefix, file_postfix]=returnOEPformat(datapath);

% Now go through the TTLs and store a window of data for each channel.
fs=30000;  %sampling rate
window_length=300; %duration of window after each shock to keep, in ms
window_sample_length=(fs*0.001*window_length);  % number of samples in the window
tbase_window=1/fs*(1:window_sample_length);  %a time vector for the window

%specify which TTLs to use in the analysis 
ttls_of_interest=fsts_2Hz; %change depending on whether 2Hz or 0.25Hz are being plotted.

nTTLs=length(ttls_of_interest); 
ttl_list=1:nTTLs;

% idiot proof in case nTTLs is set too high
if nTTLs<=length(ttls_of_interest)  
    ttl_list(nTTLs+1:end)=[];
    ttls_to_use=ttls_of_interest(ttl_list); %randomly picked list.
else
    fprintf('\n Number of TTLs higher than number in set of interest - using all TTLs \n')
    ttls_to_use=ttls_of_interest;
end

%this variable will store the extracted data:
windowed_data=zeros(nchans,window_sample_length,nTTLs);

%now extract the data. NB this is a slow process - memory mapping makes it
%even slower so just using brute force approach here.
for chan=1:nchans

   fprintf('Extracting TTL windows on channel %d \n', chan);
%  chan_name=[datapath file_prefix num2str(chan) file_postfix '.continuous']; % NEW NAME FORMAT
    chan_name=[datapath file_prefix 'CH' num2str(chan) file_postfix '.continuous'];
   [data_full, ts_full,  ~] = load_open_ephys_data(chan_name); %read in entire chan
 
   for ttl=1:nTTLs   %pull out windows around each TTL on this channel
      this_TTL=ttls_to_use(ttl);
      ttl_win_times=[this_TTL, this_TTL+(0.001*window_length)];
      inds_in_window=find(ts_full>=ttl_win_times(1) & ts_full<=ttl_win_times(2));
      data_win = data_full(inds_in_window);
      data_win(window_sample_length+1:end)=[]; %sometimes ends up with one sample too many
      try
          windowed_data(chan,:, ttl)=data_win;
      catch
          fprintf('\n Error windowing data \n')  %unhelpful error message
      end
   end
end
%% (1) Averaging over trials, look for 'hotspots' on the probe where constant latency responses are occuring.
% Extracting spikes by simple threshold crossing, and counting them in bins.
% NB This section is for plotting the 0.25Hz events

fs=30000;  %sampling rate
window_plot_length=175 * 0.001; %window to plot in heatmap, in ms
skip_time=5 * 0.001;  %time window in ms to skip at the start of each event - for excluding artefacts.
nstds=3; % sets the threshold for spike extraction in terms of number of stds from mean.

binwin=2;  %bn width in ms, for binning spikes
binwin=0.001*binwin;
binedges=0:binwin:window_plot_length; %bin edges used for each window.
bin_cents=binedges;
bin_cents(end)=[];

%hist_data is a variable which will store binned spike counts for each TTL on
%each channel:
nTTLs=25;
hist_data=zeros(nchans, numel(bin_cents), nTTLs);

for ttl=first_real_025_ind:first_real_025_ind+nTTLs-1   
   for chan=1:32
     
      data_=windowed_data(chan,:, ttl);
      spike_times=spikes_by_threshold(data_, fs, nstds);
      [bindata, ~]=histcounts(spike_times, binedges);
      hist_data(chan,:,ttl)=bindata;
      
   end
end

exclude_=find(bin_cents<skip_time)  %exclude these bins - to kill off artefact
hist_data(:,exclude_,:)=[];

%set up some labels for plotting
labels_x=[bin_cents];
gap_ind=0.05/binwin;
xtick_ind=[1:gap_ind:length(bin_cents)];

% Average over trials. For each channel, the average should begin to
% reflect level of background noise over multiple trials.
meandata=mean(hist_data,3);
meandata=[zeros(nchans,exclude_(end)), meandata]

figure('Color', 'w', 'Units', 'Pixels', 'Position',  [367.4000 109.8000 300 400])
imagesc(meandata(fliplr(chanMap),:))
xticks(xtick_ind)
xticklabels(string(labels_x(xtick_ind)*1000))
xlabel('Time (ms)')
ylabel('Channel')
yticks(1:4:32)
% yticklabels(fliplr(chanMap)); %Channel numbers as in OEP
yticklabels([32:-4:1]); %Channel numbers as in OEP
colormap hot
% title(['Average reponse to ' num2str(nTTLs) ' trials'], 'FontWeight', 'normal')
% caxis([0 0.7])
aa=gca;
aa.FontSize=14;
colorbar

%% Flatten over channels, but with the raw voltage (no spike extraction)
%  This is for the 2Hz events.
fs=30000;

%select subset of channels to plot
chans_inc=[1,32]; %max and min channels to plot (1-32, 1-18, 1-16, 6-32, 1-32)
window_plot_length=250*0.001;  %lims 50 for rec 1, then 

if ~exist('nTTLs')
    nTTLs=size(windowed_data,3);
end

skip_time2=5 * 0.0001  ;%time window in ms to skip at the start of each event - for excluding artefacts.
binwin2=0.0002;  %binning window for voltage.(.2ms for rec 2 and 4, otherwise .1ms)
binedges2=0:binwin2:window_plot_length; %bin edges used for each window.
bin_cents2=binedges2;
bin_cents2(end)=[];
nEvents=360;

n_samples_per_bin=binwin2*fs;
window_samples=window_plot_length*fs;
% 'volt_bin_data' will hold the binned voltage data.
volt_bin_data=zeros(nEvents,(window_samples/n_samples_per_bin), nchans);
for ttl=1:nEvents
   
   for chan=1:32
     
      data_=abs(windowed_data(chan,1:window_samples, ttl));
      bin_data=reshape(data_, n_samples_per_bin, [])';
      mean_bin_data=mean(bin_data,2);    
      volt_bin_data(ttl,:,chan)=mean_bin_data;
      
   end
end

exclude_=find(bin_cents2<skip_time2);

volts_chans_sel=[];
chansSel=chanMap(chans_inc(1):chans_inc(2));
volts_chans_sel(:,:)=sum( volt_bin_data(:,:,chansSel), 3);

quiet_period=find(bin_cents2>0.07); %take mean of everything after 70ms - usually after all the unit activity
mean_val=mean(volts_chans_sel(:,quiet_period),'all');
volts_chans_sel(:,exclude_)=mean_val * ones(size(volts_chans_sel, 1), numel(exclude_));

%set up some labels for plotting
%labels for plotting
labels_x=[bin_cents2];
gap_ind=0.01/binwin2;

xtick_ind=[1:gap_ind:length(bin_cents2)];

% plot - will need to adjust some of the options below depending on which
% recording is being plotted.
figure('Color', 'w', 'Units', 'pixels', 'Position', [211 356 616 400])
volts_chans_sel(volts_chans_sel>1000)=mean_val;  %kill noise on 2511 rec.
volts_chans_sel=rescale(volts_chans_sel);
imagesc(volts_chans_sel);
xticks(xtick_ind)
xticklabels(string(1000*labels_x(xtick_ind)))
xlabel('Time (ms)')
yticks(1:40:nEvents)
yticklabels(string(0:40:nEvents))
ylabel('Trial')
title(['Channels ' num2str(chans_inc(1)) ' to ' num2str(chans_inc(2))], 'FontWeight', 'normal')

caxis([0,0.13])  % Recs 1 to 5 lims, 0.3
cb=colorbar
aa=gca;
aa.FontSize=14;
cb.Ticks=[0:0.1:cb.Ticks(end)]

xlim([301, 701]);