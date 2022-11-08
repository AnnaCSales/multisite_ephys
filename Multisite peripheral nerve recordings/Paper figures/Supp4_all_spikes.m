%% Plots spikes for a list of selected clusters around TTL - for comparing clustered data 
 % with raw heatmap data

%% specify TTLs to plot
TTLchan1=2;  %this is the channel to take TTLs from
fsts=spikeStruct.TTLs.digital{TTLchan1};
fsts(2:2:end)=[];
  
inds2=find( diff(fsts)>0.47 & diff(fsts) <0.59)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);

obj_times=[850,1875];
fsts_Ob=fsts(fsts>obj_times(1)& fsts<obj_times(2));


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

ticker=11:10:(10*(nclusts+1));
%% Plotting activity around a TTL for specified clusters - change fields as required
win=[0.25, 0.5] ; %specify a window, time before and after event to consider for pinches
binwin=0.002;

% the clusters to plot (all spikes for these clusters will be included):
sel_clu=[ 3,1,2,4];

%the TTLs to plot around:
TTL_to_plot=fsts_2Hz;  %update this as needed. 
nTTL=length(TTL_to_plot);

% set up plotting 
foot_sp_fig      = figure('color','w','NumberTitle','off', 'name','Spiking around footshock / pinch TTLs', 'units', 'centimeters', 'pos',[5 2 24 17]);
figure(foot_sp_fig)

%pull out the relevant spike times
spikes_all=[];
clusterIDs_all=[]
sel_clu_ids=spikeStruct.cids(sel_clu);
spike_inds=find(ismember(spikeStruct.clu, sel_clu_ids));

ts_=spikeStruct.st(spike_inds); % all the spike times in selected cluster group
cluIDs=spikeStruct.clu(spike_inds);
event_ts=[];
spk_count_all=[];
d=subplot(2,1,1);
nTTL=numel(TTL_to_plot);

for iTTL=1:nTTL 
    
    event_ts=TTL_to_plot(iTTL);
    
    win_st=event_ts-win(1);
    win_end=event_ts+win(2);
    
    t_ind1=find(ts_>=win_st & ts_<=win_end);
    ts_window=ts_(t_ind1);  %store all the data that's been cut.
    tbin_edges = win_st:binwin:win_end;
    
    if iTTL==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binwin/2;
        t_plot=tbin_centers-event_ts;
    end
    
 
    reps=5;
    if reps==length(ts_window);
        reps=6;  %had to put this in because the plot will mess up if ts_plot is a square!
    end
    
    ts_plot=repmat(ts_window-event_ts, 1, reps);  %NB if this is a square matrix the plot will mess up as it'll go along wrong dim
    y_marks_=linspace(-0.3,0.3, reps) + (iTTL) ; %centres at  1,2,3 etc for each trial  #
    
    if ts_window
        plot(1000*ts_plot, y_marks_', 'k', 'LineWidth', 1.5);
    end
    hold on
    yticks(0:20:nTTL);
    xlim(1000*[0.02,0.15]);
    ylim([0, nTTL+1]);
    plot( zeros(1,2), [0, nTTL+0.5], 'r', 'LineWidth', 1)  %use an odd number in the plotting as laser stamps come in pairs - that way it'll never be square

end
d.Position=[0.12, 0.1, 0.7, 0.84];
xlabel('Time (ms)')
ylabel('Trial #')
set(gca, 'FontSize', 15);
box off
  set(gca, 'YDir', 'reverse');
xlim([-30,500])
yticks([0:40:359])
