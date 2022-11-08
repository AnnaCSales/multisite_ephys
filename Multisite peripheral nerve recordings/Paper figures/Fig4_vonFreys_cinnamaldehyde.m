% A plotting script that starts with a completed spikestruct for a
% multisite peripheral nerve recording - plots responses to VF/Paintbrush
% and cinnamaldehyde.

% tailored to recording for rat #3 but can be edited for #4 if needed.
%% Load in the spikeStruct.
 datapath=pwd; %dont forget backslash
load([datapath 'spikeStruct.mat']);
%for 8th Dec
% Path to ADC file with pedal on
ADC_fn=[pwd '\100_ADC2_3.continuous']; 
% Path to continuous file with ECG recording.
ECG_fn=[pwd '\100_CH34_3.continuous'];

% distance from RF to recording point in cm
rf_dist=4.7; %for 8th Dec

%for 9 Dec:
% ADC_fn=[pwd '\113_ADC2.continuous']; 
% % Path to continuous file with ECG recording.
% ECG_fn=[pwd '\113_CH34.continuous'];
% 
% rf_dist=6.0; % for 9th Dec
%% Times for digitimer TTL
TTLchan1=2;  %this is the channel to take TTLs from
fsts=spikeStruct.TTLs.digital{TTLchan1};
fsts(2:2:end)=[];

    
inds2=find( diff(fsts)>0.4 & diff(fsts) <.6)+1;
inds2=[inds2(1)-1; inds2];
fsts_2Hz=fsts(inds2);

inds025=find( diff(fsts)>3.9 & diff(fsts) <4.1)+1;
inds025=[inds025(1)-1; inds025];
fsts_025Hz=fsts(inds025);

% Obreja times (exclude stims for thresholding)
obj_times=[850,1875];
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
labels_to_keep=1: nlabels;
labels=cell(1, length(labels_to_keep));
labels(1:length(labels_to_keep))=spikeStruct.TTLs.manual.TTL_labels(labels_to_keep)
labeltimes=spikeStruct.TTLs.manual.TTL_times(labels_to_keep)
%% Times of pedal down / up
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

ticker=11:10:(10*(nclusts+1));


%% responses to paintbrush and vonFrey - bar chart

% pull out the labels that refer to events that are worth plotting (have to
% check these indices manually)
plot_labels=spikeStruct.TTLs.manual.TTL_labels([5:10, 12:end]);  %8th Dec
plot_times=spikeStruct.TTLs.manual.TTL_times([5:10, 12:end]);

% plot_labels=spikeStruct.TTLs.manual.TTL_labels([1:5, 7:end]);  %9th Dec
% plot_times=spikeStruct.TTLs.manual.TTL_times([1:5, 7:end]);


%find groups of pedal events separated by at least 12s and match to a
%label.
pedal_starts=pedal_ts(:,1);
event_starts_inds=find(diff(pedal_starts)>12)+1;
event_starts=pedal_starts([1;event_starts_inds]);  %don't forget the first group! Need to add in the '1'
%  event_starts(7)=[];  %needed in 9th Dec
% find the label closest to each event starts (should really match one to
% one but might not
pedalLabel={};
for p=1:length(event_starts)
    [~,min_ind]=min(abs(plot_times - event_starts(p)))
    pedalLabel{p}=plot_labels{min_ind};
end
% 
figure; plot(pedal_ts(:,1), 1, 'rs'); hold on; plot(pedal_ts(:,2), 1, 'ks'); plot(event_starts, 1, '*g')
for g=1:length(pedalLabel)
   text(event_starts(g),1.2, pedalLabel{g})
end
%% Now we can actually plot for unit(s) of interest
mean_frs=[];
sem_frs=[];
units_to_plot=[3,2,1,4]; %specify the units to look at.
% units_to_plot=[3,4];
pedalLabel
events_to_plot=[1:5]; %the subset of events to analyse
figure('Color', 'w', 'Position', [427.4000 51.4000 770 702])

pltNum=ceil(sqrt(length(units_to_plot)));
for iUnit=1:length(units_to_plot)    
%     subplot(1, length(units_to_plot), iUnit)
    subplot(2,2,iUnit);
    spk_times=spikeStruct.timesSorted{units_to_plot(iUnit)};
    all_fr_info=[];
    grp_info=[];
    pp=1;
    for w=1:length(events_to_plot) %go through event groups (multiple pedal presses for each event type)
        this_label=pedalLabel(events_to_plot(w));
        event_time=event_starts(events_to_plot(w));  
        % pull out all the pedal times for this group
        if w<length(pedalLabel)
            next_event_time=event_starts(w+1);
            event_group_inds=find( pedal_ts(:,1)>=event_time & pedal_ts(:,1)<event_starts(w+1));
        else
            event_group_inds=find( pedal_ts(:,1)>=event_time);
        end

        fr=[]; %for storing firing rate info for each event in this group
        for p=1:length(event_group_inds)
            event_win=pedal_ts(event_group_inds(p),:);
            event_spks=spk_times(spk_times>event_win(1) & spk_times<= event_win(2));
            fr(p)=numel(event_spks)./range(event_win);
        end

        all_fr_info=[all_fr_info,fr];  %store fr info for each pedal press in each group
        grp_info=[grp_info, pp*ones(1,length(event_group_inds))];
        pp=pp+1;
        mean_frs(w)=nanmean(fr);
        sem_frs(w)=nanstd(fr)./sqrt(sum(~isnan(fr)));
    end

    bb=bar([1:5], mean_frs)
    bb.FaceColor=rgb('RoyalBlue');
    bb.FaceAlpha=0.8;
    bb.EdgeColor='none'
    hold on
    errorbar([1:5], mean_frs, sem_frs, 'ks')
    xticks([1:5])
    xticklabels(pedalLabel(events_to_plot))
    xtickangle(45)
    ylabel('FR (Hz)')
    title(['Fibre ' num2str(iUnit)], 'FontWeight', 'normal')
    aa=gca;
    aa.YLim(1)=0;
    aa.YLim(2)=13.5;
    aa.FontSize=16;
   box off
    %do basic ANOVA
    [p,tbl,stats] = anova1(all_fr_info,grp_info, 'off');
    [c,~,~,~]=multcompare(stats, 'Display', 'off');
    for k=1:size(c,1)
        p=c(k,6);
        if p<0.05
            sigstar(c(k,1:2), p);
        end
    end
end
%% An example of a unit firing during VFs etc

cellList=[3];  %unit to plot
timePlot=[2100, 2390]; % time to plot, in seconds.
binsize=0.5; %in seconds

labels_to_keep=5:9;  %which event labels to plot (have to enter manually)
labels=cell(1, length(labels_to_keep));
labels(1:length(labels_to_keep))=spikeStruct.TTLs.manual.TTL_labels(labels_to_keep)
labeltimes=spikeStruct.TTLs.manual.TTL_times(labels_to_keep)


fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

%work out plotting positions so that the cells are plotted as they appear
%on the probe, i.e. cells detected ventrally on the probe are at the bottom
%of the window.

plot_pos2=[];
for iUnit = 1:1:length(cellList)
    
    [~,plot_pos2(iUnit)]=find(indys==iUnit) ; % a plot position for each cluster, with lower (deeper) channels getting lower numbers
    
    %plot_pos is a vector n_clusts long with a position for each unit,
    %provided in the same order as usual.
end

nclusts=length(cellList);
ticker=11:10:(10*(nclusts+1));
ticklabs={nclusts};

%set up some labels for plots below.
for pos=1:1:length(plot_pos2)
    unit_test(pos)=find(plot_pos2==pos); %The unit that is in the pos-th position on the plot
    chan_=c_channel2(unit_test(pos));   
    tt=['Clu ', int2str(cellList(unit_test(pos)))];
    ticklabs{pos}=tt;
end 

%check if there are laser TTLs included in this period - and mark on the
%plot if there are.

fsts_inc=fsts(fsts>timePlot(1)& fsts<timePlot(2));
fs2Hz_inc=fsts_2Hz(fsts_2Hz>timePlot(1)& fsts_2Hz<timePlot(2));
fs025Hz_inc=fsts_025Hz(fsts_025Hz>timePlot(1)& fsts_025Hz<timePlot(2));

fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 3 30.2 13.2], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 
p6.BackgroundColor='w';
figure(fr_fig)

FRsubplots={};
sp_height=0.6/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.07 0.3+((plot_pos2(u)-1)*1.1*sp_height) 0.88 sp_height],  'Parent', p6);
    
    iUnit=cellList(u);
    ts_= spikeStruct.timesSorted{iUnit}; 
    ts_ind=find(ts_>=timePlot(1)& ts_<=timePlot(2));
    ts_window=ts_(ts_ind);  %store all the data that's been cut.
    tbin_edges = timePlot(1):binsize:timePlot(2);

    if u==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binsize/2;
        t_plot=tbin_centers;
    end
    
    %bin and histogram
    spk_count = histc(ts_window,tbin_edges);
    spk_count = spk_count(1:end-1) ./binsize;
    %plot
    mybar=bar(t_plot, spk_count);
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.7;
    mybar.LineWidth=1;
   
    xlabel('Time (s)')   
    ylabel('FR (Hz)', 'FontSize', 10 )

    hold on    
    aa=gca;
    aa.FontSize=14;
    ymax=aa.YLim(2)
    ymin=aa.YLim(1)
    ymax=0;
    ymin=35;
    if fsts_inc
       plot(repmat(fsts_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.05], 'LineWidth', 1);
    end
    
    aa=gca;
    lower=aa.YLim(1);
    upper=aa.YLim(2);
    lower=0;
    upper=35;
    for q=1:size(pedal_ts, 1)
        mypat=patch([pedal_ts(q,1), pedal_ts(q,2), pedal_ts(q,2), pedal_ts(q,1)],[lower,lower,upper, upper], 'g')
        mypat.FaceColor=[0 1 0.1];
        mypat.FaceAlpha=0.2;
        mypat.EdgeColor='none';
    end
    
    xlim([timePlot(1), timePlot(2)])
    ylim([0,max(spk_count)+2])
    box off
end       
    

  MARKERsubplots=subplot('Position',[0.065 0.04 0.88 0.05],  'Parent', p6);
  for r=1:length(labeltimes)
      text(labeltimes(r), 0.1, labels{r}, 'Rotation', 90 , 'Fontsize', 12)
  end
  set(gca,'xtick',[]) 
  set(gca,'ytick',[]) 
%   XAxis. TickLength = [0 0];
  xlim([timePlot(1), timePlot(2)])
 
%% Cinnamaldehyde - spontaneous firing

cellList=[4]
cin_time=spikeStruct.TTLs.manual.TTL_times(11); %8th Dec
% cin_time=spikeStruct.TTLs.manual.TTL_times(6); %9th Dec
timePlot=[60, 300]; % time to plot, in seconds.
binsize=1; %in seconds

fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

nclusts=length(cellList);
ticker=11:10:(10*(nclusts+1));
ticklabs={nclusts};

%set up some labels for plots below.
for pos=1:1:length(cellList)
    tt=['Clu ', num2str(cellList(pos))];
    ticklabs{pos}=tt;
end 



fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 3 30.5 10], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 
p6.BackgroundColor='w';
figure(fr_fig)

FRsubplots={};
sp_height=0.6/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.065 0.1+((u-1)*1.45*sp_height) 0.88 sp_height],  'Parent', p6);
    
    iUnit=cellList(u);
    ts_= spikeStruct.timesSorted{iUnit}-cin_time; 
    ts_ind=find(ts_>=-timePlot(1)& ts_<=timePlot(2));
    ts_window=ts_(ts_ind);  %store all the data that's been cut.
    tbin_edges = -timePlot(1):binsize:timePlot(2);

    if u==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binsize/2;
        t_plot=tbin_centers;
    end
    
    %bin and histogram
    spk_count = histc(ts_window,tbin_edges);
    spk_count = spk_count(1:end-1) ./binsize;
           
    %plot
    mybar=bar(t_plot, spk_count);
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.7;
    mybar.LineWidth=1;
    mybar.EdgeColor='none'; hold on
%     plot(t_plot, spk_count)
   
    xlabel('Time (s)')   
    ylabel({'FR (Hz)'}, 'FontSize', 10 )
    ylim([0,15])
    aa=gca;   
    plot([0,0], [0, aa.YLim(2)], 'r', 'LineWidth', 2)
  
    pedal_inc=find(pedal_ts(:,1)>cin_time-timePlot(1) & pedal_ts(:,2)<cin_time+timePlot(2))
    if ~isempty(pedal_inc)
        plot(pedal_ts(pedal_inc,1)-cin_time, 5, '*g')
    end
    
    xlim([-timePlot(1), timePlot(2)])
    
    if u==length(cellList)
        text(-25, aa.YLim(2)+1, 'Cinnamaldehyde', 'FontSize', 13)
    end
    aa.FontSize=16
    box off
end 
%% Cinnamaldehyde, before and after

cellList=[3,1,4]
cin_time=spikeStruct.TTLs.manual.TTL_times(11); %8th Dec
% cin_time=spikeStruct.TTLs.manual.TTL_times(6); %9th Dec
timePlot=[60, 60]; % time to plot, in seconds.
binsize=2; %in seconds

fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

nclusts=length(cellList);
ticker=11:10:(10*(nclusts+1));
ticklabs={nclusts};

%set up some labels for plots below.
for pos=1:1:length(cellList)
    tt=['Clu ', num2str(cellList(pos))];
    ticklabs{pos}=tt;
end 

figure('Color', 'w', 'Position', [3.4000 449 980 280])
pltNum=ceil(sqrt(length(units_to_plot)));
all_bef=[];
all_aft=[];
xlabs={};
pvals=[];
for u=1:length(cellList)
    subplot(1, length(cellList), u)
   
    iUnit=cellList(u);
    xlabs{u}=['Fibre ' num2str(iUnit)];
    ts_= spikeStruct.timesSorted{iUnit}-cin_time; 
    ts_ind=find(ts_>=-timePlot(1)& ts_<=timePlot(2));
    ts_window=ts_(ts_ind);  %store all the data that's been cut.
    tbin_edges = -timePlot(1):binsize:timePlot(2);

    if u==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binsize/2;
        t_plot=tbin_centers;
        t_bef=find(tbin_centers<0);
        t_aft=find(tbin_centers>=0);
    end
    
    %bin and histogram
    spk_count = histc(ts_window,tbin_edges);
    spk_count = spk_count(1:end-1) ;
          
    spk_count_bef=spk_count(t_bef);
    spk_count_aft=spk_count(t_aft);
    
    [p,h]=ranksum(spk_count_bef, spk_count_aft)
    mean_bef=mean(spk_count_bef);
    sem_bef=std(spk_count_bef)./sqrt(numel(spk_count_bef));   
    mean_aft=mean(spk_count_aft);
    sem_aft=std(spk_count_aft)./sqrt(numel(spk_count_aft))
    
    %store
    all_bef(:,u)= spk_count_bef
    all_aft(:,u)= spk_count_aft;
    
     
    mybar=bar([1], [mean_bef]); 
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.6;
    mybar.LineWidth=1;
    mybar.EdgeColor='none';
    hold on
    mybar=bar([2], [mean_aft]);
    mybar.FaceColor=[1 0 0];
    mybar.FaceAlpha=0.6;   
    mybar.LineWidth=1;
    mybar.EdgeColor='none';
    errorbar([1,2],[mean_bef, mean_aft],[sem_bef, sem_aft], 'ks')
    
   if p<0.05
       sigstar([1,2], p)
   end
   pvals(u)=p;
   xticks([1,2])
   xticklabels({'Pre-cinn' 'Post-cinn'})
   xtickangle(45)
   ylabel('FR (Hz)')
   xlim([0.2, 2.8])
   aa=gca;  
   aa.FontSize=15;
   aa.YLim(1)=0;
   box off
   if aa.YLim(2)<5
       aa.YLim(2)=0.85;
   else
       aa.YLim(2)=18;
   end
   
end 
%% Before and after cinnamaldehyde, but on the same plot
figure('Color', 'w','Position', [489 409 708 288.8000])
b1=bar([1,4,7],mean(all_bef, 1) )
b1.FaceColor=[0 0 1];
b1.FaceAlpha=0.7;
b1.LineWidth=1;
b1.EdgeColor='none'; hold on
b1.BarWidth=0.3;
errorbar([1,4,7],mean(all_bef, 1) ,std(all_bef, 1)/sqrt(size(all_bef,1)) , 'ks')

b2=bar([2,5,8],mean(all_aft, 1) )
b2.FaceColor=[1 0 0];
b2.FaceAlpha=0.7;
b2.LineWidth=1;
b2.EdgeColor='none'; hold on
b2.BarWidth=0.3;
errorbar([2,5,8],mean(all_aft, 1) ,std(all_aft, 1)/sqrt(size(all_aft,1)) , 'ks')
xticks([1.5, 4.5, 7.5])
xticklabels(xlabs)

for pv=1:length(pvals)
    if pvals(pv)<0.05
        xv=((pv-1) * 3)+1;
        sigstar([xv, xv+1], pvals(pv))
    end
end

aa=gca;
aa.FontSize=16;
box off
ylabel(['Spikes per ' num2str(binsize) 's bin'])
xlim([0.2 8.7])
ll=legend([b1, b2], {'60s baseline' '60s post-cinnamaldehyde'}, 'Fontsize', 12)
ll.EdgeColor='none';
ll.Location='northwest';

%%
%% Before and after cinnamaldehyde, but on the same plot, JPD trying to bring the three plots in to a single before and after bar chart... 28th Sept 2021

M_bef=mean(all_bef);
M_aft=mean(all_aft);

dat_jd=[(transpose(M_bef)), (transpose(M_aft))];
x = 1:2;
y = dat_jd;   
ymean = mean(y); 

%%
figure();
hold on; % plot multiple things without clearing the axes
plot(x, y, 'ok' );   % scatter of the data. 'o' for marker, 'k' for black
line(x,y,'LineStyle', '-', 'Color', 'k');

hold on;

jd_b1=bar([1],mean_bef);
jd_b1.FaceColor=[0 0 1];
jd_b1.FaceAlpha=0.7;
jd_b1.LineWidth=1;
jd_b1.EdgeColor='none'; hold on
jd_b1.BarWidth=0.3;

hold on;

jd_b2=bar([2],mean_aft)
jd_b2.FaceColor=[1 0 0];
jd_b2.FaceAlpha=0.7;
jd_b2.LineWidth=1;
jd_b2.EdgeColor='none'; hold on
jd_b2.BarWidth=0.3;

xticks([1 2])
xticklabels({'Baseline','Post Cinnamaldehyde'})


xx = dat_jd(:,1);
yy = dat_jd(:,2);
[h,p] = ttest(xx,yy)

sigstar({[1,2]})

aa=gca;
aa.FontSize=16;
box off
ylabel(['Mean Spikes per ' num2str(binsize) 's bin'])

%% Optional plot for looking at spontaneous activity over any specified time period

cellList=[1,2,3,4]
timePlot=[2000, 4000]; % time to plot, in seconds.
binsize=0.5; %in seconds

labels_to_keep=5:9;
labels=cell(1, length(labels_to_keep));
labels(1:length(labels_to_keep))=spikeStruct.TTLs.manual.TTL_labels(labels_to_keep)
labeltimes=spikeStruct.TTLs.manual.TTL_times(labels_to_keep)


fs=spikeStruct.sample_rate;  %sampling rate
c_channel=spikeStruct.c_channel;  %centre channel for the cluster

c_channel2=c_channel(cellList);
[vals, indys]=sort(c_channel2);

%work out plotting positions so that the cells are plotted as they appear
%on the probe, i.e. cells detected ventrally on the probe are at the bottom
%of the window.

plot_pos2=[];
for iUnit = 1:1:length(cellList)
    
    [~,plot_pos2(iUnit)]=find(indys==iUnit) ; % a plot position for each cluster, with lower (deeper) channels getting lower numbers
    
    %plot_pos is a vector n_clusts long with a position for each unit,
    %provided in the same order as usual.
end

nclusts=length(cellList);
ticker=11:10:(10*(nclusts+1));
ticklabs={nclusts};

%set up some labels for plots below.
for pos=1:1:length(plot_pos2)
    unit_test(pos)=find(plot_pos2==pos); %The unit that is in the pos-th position on the plot
    chan_=c_channel2(unit_test(pos));   
    tt=['Clu ', int2str(cellList(unit_test(pos)))];
    ticklabs{pos}=tt;
end 

%check if there are laser TTLs included in this period - and mark on the
%plot if there are.

fsts_inc=fsts(fsts>timePlot(1)& fsts<timePlot(2));
fs2Hz_inc=fsts_2Hz(fsts_2Hz>timePlot(1)& fsts_2Hz<timePlot(2));
fs025Hz_inc=fsts_025Hz(fsts_025Hz>timePlot(1)& fsts_025Hz<timePlot(2));

fr_fig= figure('color','w','NumberTitle','off', 'name',' Unit firing', 'units', 'centimeters', 'pos',[5 3 24 15], 'Color', 'white');
p6 = uipanel('Parent',fr_fig,'BorderType','none'); 
p6.BackgroundColor='w';
figure(fr_fig)

FRsubplots={};
sp_height=0.55/length(cellList);
if sp_height > 0.3
    sp_height=0.24
end


for u=1:length(cellList)
    
    FRsubplots{u}=subplot('Position',[0.09 0.1+((plot_pos2(u)-1)*1.7*sp_height) 0.84 sp_height],  'Parent', p6);
    
    iUnit=cellList(u);
    ts_= spikeStruct.timesSorted{iUnit}; 
    ts_ind=find(ts_>=timePlot(1)& ts_<=timePlot(2));
    ts_window=ts_(ts_ind);  %store all the data that's been cut.
    tbin_edges = timePlot(1):binsize:timePlot(2);

    if u==1 %store a time vector for plotting
        tbin_centers = tbin_edges(1:end-1)+binsize/2;
        t_plot=tbin_centers;
    end
    
    %bin and histogram
    spk_count = histc(ts_window,tbin_edges);
    spk_count = spk_count(1:end-1) ./binsize;
    %plot
    mybar=bar(t_plot, spk_count);
    mybar.FaceColor=[0 0 1];
    mybar.FaceAlpha=0.7;
    mybar.LineWidth=1;
   
    xlabel('Time (s)')   
    ylabel({ticklabs{u} 'FR (Hz)'}, 'FontSize', 10 )

    hold on    
    aa=gca;
    aa.FontSize=12;
    ymax=aa.YLim(2)
    ymin=aa.YLim(1)
    ymax=0;
    ymin=35;
    if fsts_inc
       plot(repmat(fsts_inc,1,3), [ymin, ymax, ymax], '-r', 'Color', [1 0 0 0.05], 'LineWidth', 1);
    end
    
    aa=gca;
    lower=aa.YLim(1);
    upper=aa.YLim(2);
    lower=0;
    upper=35;
    for q=1:size(pedal_ts, 1)
        mypat=patch([pedal_ts(q,1), pedal_ts(q,2), pedal_ts(q,2), pedal_ts(q,1)],[lower,lower,upper, upper], 'g')
        mypat.FaceColor=[0 1 0.1];
        mypat.FaceAlpha=0.2;
        mypat.EdgeColor='none';
    end
    
    xlim([timePlot(1), timePlot(2)])
    ylim([0,max(spk_count)+2])
end       
   