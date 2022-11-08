function [corfig] = cell_xcors(spikeStruct, params)

% This script plots the auto-correlograms and x-correlograms of selected
% units, in a way that looks similar to the plot provided in the
% PHY GUI clustering package. Allows quick reconstruction of what's seen in
% PHY to check units are well isolated. 

% Expects to be passsed (a) the spikeStruct (b) a params stuct, with the
% following fields

% params.cell_list  - a vector containing the cluster numbers to plot.
% params.binsize - binsize, in ms
% params.window - window over which to plot the xcor.
% params.baseline - 1 if we're plotting only over baseline, 0 if not.
% Updated July 2020, Anna Sales.


%% specify parameters for the auto/x corr
ac_binsize=params.binsize /1000  % binsize for the a/x cor (done for the whole recording, not just baseline)
windo=params.window/1000; %window to show in the plots.

%% specify the time window in the recording to use (here we use the baseline)
if ~isempty(params.baseline)
    if numel(params.baseline)==1
     min_t=spikeStruct.baseline_st;
     max_t=spikeStruct.baseline_end;
    else
     min_t=params.baseline(1);
     max_t=params.baseline(2);
    end
else
    min_t=0;
    max_t=spikeStruct.timeRange(2);
end

%% specify the cells 
cells2plt=params.cell_list;

%% X-correlations and auto-cors - like in PHY

sp_=length(cells2plt) ; %for calculating the number of subplots required
bincount=[];

for iU=1:length(cells2plt)
    iUnit=cells2plt(iU)
    ts_= spikeStruct.timesSorted{iUnit}; 
        
    [bincount, ~] = histcounts(ts_, min_t:ac_binsize:max_t);   %bin size is 5ms
    
    if iU==1 %if we're in the loop for the first time, define the size of bl_bincounts
        bl_bincounts=zeros(length(cells2plt), length(bincount));
    end
        
    bl_bincounts(iU,:)=bincount;        
end

% Plot xcors for the units in the candidate group (whole recording, to
% recreate what's in PHY)

fprintf('\n Calculating xcorrs')
[xc, lags]=xcorr(bl_bincounts', 'coeff'); %'lags' will run to/from +/- (centres(end)-1)
% - each col is an x-correlation. The zero lag is in the middle of each
%col, at row=number of centres (lag runs from -centres to +centres). The
%columns are in sets, with nclust cols for each cluster (the first col is
%with itself, then with each of the others in turn)

fprintf('\n Plotting xcorrs')
plotpos=1;
nGD=length(cells2plt)
corfig=figure('Name', ['Auto / X correlogram for each cluster: ' num2str(ac_binsize*1000) ' ms bins'], 'Color', 'w', 'units', 'centimeters', 'pos', [10 0.49 26.8 18.8])
plot_log=[0 0];

tbase=lags*ac_binsize;
incinds=find( tbase>=-windo & tbase <=windo);    %only plot the central stuff 
zeroind=find(~tbase(incinds));

for ttt=1:nGD
    figure(corfig)
    cols_set=((ttt-1)*nGD); %the set of columns we need.
    target=ttt;
    target_label=cells2plt(ttt);

    pp=ttt;
    
    for t=1:nGD   %get correlograms we don't have yet i.e don't calcuate 1vs 2 if we already have 2vs1
       comparator=t;
       comp_label=cells2plt(t);    
       col_=cols_set+t;
       zero=find(~lags);  
       
       if ~sum(ismember(plot_log, sort([ttt, t]), 'rows'))
           subplot(sp_,sp_, pp);
           co=xc(:, col_);
           corrs=co(incinds);
           
           if t==ttt              
               corrs(zeroind)=0; % don't plot zero for acor          
           end

           xbar=bar(1000*tbase(incinds), corrs);
           hold on
           xbar.BarWidth=0.9;
           xbar.FaceColor='b';
           xbar.FaceAlpha=0.7;
           xbar.EdgeColor='none';
           xbar.LineWidth=0.001;
           set(gca, 'FontSize', 14)
           hold on
           
           xlim(1000*[-windo, windo]);
           title(['#' int2str(target_label) ' vs #' int2str(comp_label)], 'FontSize', 10, 'FontWeight', 'normal')  ;
           aa=gca;
           aa.YLim(1)=0;
%            aa.YLim(2)=round(max(corrs(2:end)),2);
       end
       pp=pp+sp_;
       plot_log=[plot_log;sort([ttt, t])]; %keep track of which one
     
    end
end

xlabel('Time (ms)', 'FontSize', 8)
ylabel('xcor')


end

