%% ss3 - plot division events in single shift environments


% Goal: what is the lag between nutrient shift and changes in divison rate?

%  Strategy:
%
%     0. initialize complete meta data
%     0. for experiments of interest...
%           1. collect experiment date and exclude outliers (2017-10-31 and monod experiments)
%           2. initialize binning parameters
%           3. load measured data for stable condition
%           4. for stable average... 
%                  5. isolate isDrop and timestamp data
%                  6. correct time based on calculated lag in signal between junc and xy position
%                  7. isolate birth events (isDrop == 1) and corresponding timestamps
%                  8. PLOT ONE: number of tracked cells over time
%                  9. PLOT TWO: birth events over time, normalized by tracks per bin
%                 10. PLOT THREE: births normalized by cell count over period fraction
%          11. repeat analysis for fluctuating environment, plotting fluc data over stable



%           5. find average growth rate of stable average condition
%                       i. isolate data
%                      ii. remove data not in stabilized region
%                     iii. remove zeros from mu data (always bounding start and end of tracks)
%                      iv. calculate mean value for of mu and bvpr in stable
%           6. for fluctuating condition, load measured data
%                       i. isolate data of interest
%                      ii. normalize mu and bvpr data by mean of stable average condition
%                     iii. remove data not in stabilized region
%                      iv. remove zeros from mu data (always bounding start and end of tracks)
%           7. accumulate data by shifted time bin (period fraction)
%                       i. from original timestamp, subtract shift = period/4 + offset
%                      ii. re-define period to begin at start of high nutrient pulse
%                     iii. bin data by period fraction
%           8.  convert bin # to absolute time (in seconds)
%           9.  calculate average and s.e.m. per timebin
%          10.  plot, with repeated high nutrient half period at end
%    11. repeat for all fluctuating experiments



%  Last edit: Jen Nguyen, 2020 March 6
%  Commit: probability of division in single shift experiments

% OK! Lez go!


%% Part 0. initialize analysis

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
%dataIndex = find(~cellfun(@isempty,storedMetaData));
%experimentCount = length(dataIndex);


% 0. designate index of experiments in meta data
exp_upshifts = 21:22;
exp_downshifts = 26:27;
exp = [exp_upshifts, exp_downshifts];
  

% 0. initialize time binning of division events
timebin = 2; % min
binsPerHr = 60/timebin;
 

%% Part 1. plot division events over time for each experimental replicate

for ee = 1:length(exp)
    
    % 1. initialize experiment meta data
    index = exp(ee);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    shiftTime = storedMetaData{index}.shiftTime./60;
    disp(strcat(date, ': analyze!'))
    
    % 2. load measured experiment data    
    dataFolder = '/Users/jen/Documents/StockerLab/Source_data/';
    cd(dataFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    clear filename dataFolder
    
    
    % 3. for stable average and then fluctuating environment...
    environment = [1;2;4];
    ccData = cell(length(environment),1);
    
    for ii = 1:length(environment)
              
        
        % 4. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        condition = environment(ii);
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end
        
        
        % 5. determine timepoints at which births occur
        timestamps_h = getGrowthParameter(conditionData,'timestamp')./3600;  % timestamp in seconds
        birthEvents = getGrowthParameter(conditionData,'isDrop');        % isDrop, 1 marks a birth event;
        trackNum = getGrowthParameter(conditionData,'trackNum');         % track count (not ID from particle tracking)
       
        
        % 6. bin birth times and tracks
        bins = ceil(timestamps_h*binsPerHr);
        births_binned = accumarray(bins,birthEvents,[],@(x) {x});
        trackNum_binned = accumarray(bins,trackNum,[],@(x) {x});
        
        
        % 7. compute total birth events and tracks per time bin
        births = cellfun(@sum,births_binned);
        tracks_unique = cellfun(@unique,trackNum_binned,'UniformOutput',false);
        tracks = cellfun(@length,tracks_unique);
        
        
        % 8. plot probability of division!
        if condition == 1
            color = rgb('DodgerBlue');
        elseif condition == 2
            color = rgb('Indigo');
        elseif condition == 4
            color = rgb('FireBrick');
            
            figure(ee)
            title(date)
        end
        x = linspace(0,length(births)*timebin,length(births)+1);
        x = x(2:end);
        
        figure(ee)
        hold on
        plot(x,births./tracks,'o','Color',color)
        
        
    end
    clear x
    
    % 9. plot shift time
    y = linspace(0,0.12,10);
    x = ones(10,1)*shiftTime;
    
    figure(ee)
    hold on
    plot(x,y,'Color',rgb('SlateGray'))
    ylabel('fraction of tracks dividing')
    xlabel('time (min)')
    

end
    




