%% ss3 - probability of division in single shift experiments


% Goal: what is the lag between nutrient shift and changes in divison rate?

%  Strategy:
%
%     Part 0. initialize analysis
%     Part 1. plot division events over time for each experimental replicate
%             i.e. one plot per experiment

%          1. initialize experiment meta data
%          2. load measured experiment data   
%          3. loop through shift, steady low and steady high conditions to...
%          4. compile condition data matrix
%          5. determine timepoints at which births occur
%          6. bin birth times and tracks
%          7. compute total birth events and tracks per time bin
%          8. plot probability of division!
%          9. plot shift time



%  Last edit: Jen Nguyen, 2020 March 9
%  Commit: edit comments for clarity

% OK! Lez go!


%% Part 0. initialize analysis

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')


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
    
    
    % 3. loop through shift, steady low and steady high conditions to... 
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
    




