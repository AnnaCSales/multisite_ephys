function [norm_xv, norm_yv, v_mag] = mean_spk_CV_and_dir_select_chans(spikeStruct, plotparams)
% Estimates conduction velocity using the mean waveform across multiple
% channels.
% Excludes spikes in specified period ( times specified by 'events',
% duration to exclude specified by 'exclusionWindow')
% Parameters set feature to be tracked (trough +1 or peak -1), channels to
% include ('channels')


events=plotparams.events; %events to exclude - send empty vector if none
duration_s=plotparams.exclusionWindow; %window after TTL in which to ignore spikes, in seconds.
fn=plotparams.fn; %filename of raw data, with path if not in pwd
unit=plotparams.unit; %the unit to plot
max_min=plotparams.max_min; %1 to track peaks, -1 to track troughs
chans_to_use=plotparams.channels;

% First sort out the times of spikes during 'events', and those not during
% events.
wfparams.dataDir=[pwd '\'];; %where is raw data
wfparams.fileName=fn;
wfparams.nCh = 32;
wfparams.dataType='int16';
wfparams.wfWin = [-40 41];  %number of samples before and after spike peak
wfparams.nWf = 2000;   %number of waveforms to return (if they are there)
wfparams.nBad=0; %nobad channels.
fs=spikeStruct.sample_rate;

%make a table of times to consider,e.g during TTLs
spikes_duringTTL=[];
clusterIDs_duringTTL=[];  
all_spk_ts=spikeStruct.st;
all_spk_IDs=spikeStruct.clu;
% duration_s=duration/1000;
all_spk_tsNOT=all_spk_ts; %we will cut out any spikes during laser from these
all_spk_IDsNOT=all_spk_IDs;
all_inds_duringTTLs=[];

% Now extract cluster IDs for the unit we care about

clu_ids_for_unit=spikeStruct.cids(unit);


for p=1:length(events)   %get a list of spikes / cluster IDs which happened during laser.
    event_time=events(p);
    inds_duringTTL=find(all_spk_ts>event_time & all_spk_ts< (event_time+duration_s) ); 
    
    %spike times in the window of interest...
    all_inds_duringTTLs=[all_inds_duringTTLs; inds_duringTTL];     
end

fprintf('\n Extracting waveforms at times other than after specified TTLs...')
all_spk_tsNOT(all_inds_duringTTLs)=[];
all_spk_IDsNOT(all_inds_duringTTLs)=[];
wfparams.spikeTimes=round(all_spk_tsNOT * fs);

this_unit_ind=find(all_spk_IDsNOT==clu_ids_for_unit);%spike times NOT during TTLs, all clusters
wfparams.spikeClusters=all_spk_IDsNOT(this_unit_ind);  %clu IDs just for this unit
all_unit_times=round(all_spk_tsNOT * fs);
wfparams.spikeTimes=all_unit_times(this_unit_ind); %spk times for this unit only.


wf_NOTevent=getWaveForms_with_std_and_bad(wfparams);   %extract the waveforms
wfs=squeeze(wf_NOTevent.waveFormsMean);

wfs_sem=squeeze(wf_NOTevent.waveFormsSTD);
wfs_n_extracted=wfparams.nWf;
wave_time=1000/spikeStruct.sample_rate *(1:size(wfs,2)); %convert to milliseconds

xs=spikeStruct.xcoords;
ys=spikeStruct.ycoords;

%now work out where to plot, to accruately represent channel layout
cent_chan=spikeStruct.c_channel(unit);  %the centre channel for this unit
cent_y=ys(cent_chan);


%determine if we are looking for trough or peak
cent_wf=wfs(cent_chan,:)
wf_scale=round([floor(min(cent_wf))-50, ceil(max(cent_wf))+50  ]); %set scale for plotting wfs

if isfield(plotparams, 'max_min')
    max_min=plotparams.max_min; %1 to track peaks, -1 to track troughs
else
    if cent_wf(41)<0
        max_min=-1;  %trough
    else
        max_min=1;
    end
end



% %now scale channel coords from 0.1 to 0.75, to provide coordinates for the axes
% %representing each channel
plotPos=[rescale(xs(chans_to_use),0.15, 0.75), rescale(ys(chans_to_use), 0.06, 0.85)];

% work out which channel is in the bottom left of the plot - this is where
% we will show axis info

[low_val,low_ind]=min(plotPos(:,2));
lowest=find(plotPos(:,2)==low_val);
[~,leftest]=min(plotPos(lowest, 1));
bottom_left_chan=chans_to_use(lowest(leftest));

wavefig_diff=figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.39 0.75]);
for sb=1:length(chans_to_use)
    chan_=chans_to_use(sb);
    
    if chan_==cent_chan  %plot the centre channel in red to highlight
        pltcol='b';
    else
        pltcol='b';
    end
   
    %put the plot in the correct place.
    chanPlots(sb)=subplot('Position', [plotPos(sb,1), plotPos(sb,2), 0.18, 0.1]);
    seb=shadedErrorBarLight(wave_time, wfs(chan_, :), wfs_sem(chan_,:)./sqrt(wfs_n_extracted), pltcol, 1);
    seb.mainLine.LineWidth=1
    ylim(wf_scale);
    xlim([0,wave_time(end) ]);
    try
    title(['Ch. ' num2str(chan_)], 'Fontweight', 'normal') ;
    catch
        fprintf('Boo')
    end
    ax = gca;
    ax.FontSize = 12;
    
    if chan_==bottom_left_chan
        xlabel('ms');
        ylabel('\mu V');
        box off
    else
       set(gca,'xtick',[])
       set(gca,'ytick',[])
       box off
       axis off
    end
 
    
end
subplot(chanPlots(1))
aa=gca;
text(3.3, 120, ['Unit: ' num2str(unit)]);

%%
wavefig_same=figure('Color', 'w', 'Units', 'normalized', 'Position',  [0.1635 0.4898 0.46 0.3380]);

wf_scale=round([floor(min(cent_wf))-50, ceil(max(cent_wf))+50  ])
leglabels={};
trough_times=zeros(1, length(chans_to_use));
cols=[rgb('Red');rgb('Blue');rgb('MediumSeaGreen');rgb('DarkSlateGrey');rgb('Gold');rgb('Maroon');rgb('Orange');rgb('Red');rgb('Red');rgb('Red');rgb('Red');rgb('Red');rgb('Red');rgb('Red')]

for sb=1:length(chans_to_use)
    chan_=chans_to_use(sb);
    pltcol='b'

    pltcol=cols(sb, :)
  
   wf_=wfs(chan_, :);
   plot(wave_time, wf_, 'LineWidth', 1.5, 'Color', pltcol);
   hold on
   
   %*****************************%
   feat_win=35:47; %samples where feature to be tracked will be found
   
   if max_min==1 %find max
       [~, feat_ind]= max(wf_(feat_win));
   else
       [~, feat_ind]= min(wf_(feat_win));
   end
   
   feat_ind_abs=feat_win(feat_ind)
   
   trough_times(sb)=wave_time(feat_ind_abs); %the trough
   leg_labels{sb}=num2str(chan_);
   hold on
    ylim(wf_scale);
    xlim([0,wave_time(end) ]);
    title('Waveforms across channels', 'Fontweight', 'normal') ;
 
    
    
end
figure(wavefig_same)
ll=legend(leg_labels)
ll.Location='southeast'
xlabel('Time (ms)')
ylabel('\muV')
ax = gca;
ax.FontSize = 16;
%sort ascending

[~, sort_ind]=sort(trough_times);
sorted_chans=chans_to_use(sort_ind)';
trough_times_sorted=trough_times(sort_ind);
% find the distances between all channel pairs

for g=1:length(sorted_chans)
    chan1=sorted_chans(g);
    x1=spikeStruct.xcoords(chan1);
    y1=spikeStruct.ycoords(chan1);
    for h=g:length(chans_to_use)
         chan2=sorted_chans(h);
         x2=spikeStruct.xcoords(chan2);
         y2=spikeStruct.ycoords(chan2);
         
         delta_t=(trough_times_sorted(h)-trough_times_sorted(g))*0.001; %convert from ms to s
         delta_x=(x2-x1)*0.000001; %convert from um to m
         delta_y=(y2-y1)*0.000001;
         time_diff(g,h)=delta_t; 
         x_dist(g,h)=delta_x;
         y_dist(g,h)=delta_y;
         x_vel(g,h)=delta_x/delta_t; 
         y_vel(g,h)=delta_y/delta_t;
                     
    end
end


x_dist_table=array2table(x_dist);
x_dist_table.Properties.VariableNames=string(sorted_chans);
x_dist_table.Properties.RowNames=string(sorted_chans);
x_vel_table=array2table(x_vel);
x_vel_table.Properties.VariableNames=string(sorted_chans);
x_vel_table.Properties.RowNames=string(sorted_chans);
y_dist_table=array2table(y_dist);
y_dist_table.Properties.VariableNames=string(sorted_chans);
y_dist_table.Properties.RowNames=string(sorted_chans);
y_vel_table=array2table(y_vel);
y_vel_table.Properties.VariableNames=string(sorted_chans);
y_vel_table.Properties.RowNames=string(sorted_chans);

idx = ones(numel(sorted_chans),numel(sorted_chans));
idx = tril(idx) ; % Make 1's at lower triangular part 

x_vel(idx==1) = NaN ;     % replace lower triangular part with 1
x_vel(isinf(x_vel))=NaN;
mean_x_vel=nanmean(reshape(x_vel, numel(x_vel), []));
y_vel(idx==1) = NaN ;     % replace lower triangular part with 1
y_vel(isinf(y_vel))=NaN;
mean_y_vel=nanmean(reshape(y_vel, numel(y_vel), []));

norm_xv=mean_x_vel/sqrt(mean_x_vel^2 + mean_y_vel^2);
norm_yv=mean_y_vel/sqrt(mean_x_vel^2 + mean_y_vel^2);

% figure(wavefig_diff)
% for h=1:length(chans_to_use)
%     subplot(chanPlots(h));
%     place_in_order=find(sort_ind==h);
%     text(0,0, num2str(place_in_order));
% end

fprintf('Action potential vector = %2.2f i + %2.2f j m/s', norm_xv, norm_yv)
v_mag=sqrt(mean_x_vel^2  + mean_y_vel ^2);

if norm_xv<0
    startx=max(xs);
else
    startx=min(xs);
end

if norm_yv<0
    starty=max(ys);
    text_offset=-10;
else
    starty=min(ys);
    text_offset=+10;
end

endx=startx+100*v_mag*norm_xv;
endy=starty+100*v_mag*norm_yv;

figure('Color', 'w','Units', 'Normalized', 'Position',[0.3177 0.2769 0.124 0.6051] );
% subplot('Position', [0.15 0.15 0.7, 0.7])
pp=plot(xs, ys, 'ko', 'Color', rgb('Navy'),'MarkerSize', 7)
hold on
% arrow3([startx,starty],[endx, endy]);
plot([startx, endx], [starty, endy], 'k', 'LineWidth', 1.5)
text(endx+1, endy+text_offset, [num2str(v_mag,2), ' ms^-1']);
xlim([min(xs)-10, max(xs)+10]);
ylim([min(ys)-20, max(ys)+20]);
title(['Mean spike velocity vector, unit ' num2str(unit)], 'Fontweight', 'normal');
xlabel('x (\mum)')
ylabel(('y (\mum)'))
aa=gca;
aa.FontSize=16;
box off
set(gca,'XColor', 'none','YColor','none')
% axes('Color','none','XColor','none');
