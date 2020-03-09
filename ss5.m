%% ss5 - single curves after a single downshift


% Goal: plot individual volume trajectories starting with cell that
%       experiences a single nutrient upshift


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
%  Commit: single cell volume tracks across single nutrient downshifts

% OK! Lez go!


%% Part 0. initialize analysis

% 0. initialize data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')


% 0. designate index of experiments in meta data
%exp_upshifts = 21:22;
exp_downshifts = 26:27;
%exp = [exp_upshifts, exp_downshifts];
  

% 0. initialize time binning of division events
timebin = 2; % min
binsPerHr = 60/timebin;
 

%% Part 1. plot division events over time for each experimental replicate

for ee = 1%:length(exp_upshifts)
    
    % 1. initialize experiment meta data
    index = exp_downshifts(ee);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    shiftTime = storedMetaData{index}.shiftTime;
    disp(strcat(date, ': analyze!'))
    
    
    % 2. load measured experiment data
    dataFolder = '/Users/jen/Documents/StockerLab/Source_data/';
    cd(dataFolder)
    filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    clear filename dataFolder
    
    
    % 3. for stable average and then fluctuating environment...
    environment = 1; %2;4];

    
    for ii = 1:length(environment)
        
        
        % 4. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        condition = environment(ii);
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end
        
        
        % 5. trim data to full cell cycles only
        curveIDs = getGrowthParameter(conditionData,'curveFinder');
        fullData = conditionData(curveIDs > 0,:);
        clear curveIDs
        
        
        % 6. identify cell cycles in which nutrient shift occurs
        tracks = getGrowthParameter(fullData,'trackNum');
        tracks_unique = unique(tracks);
        
        % i. eliminate tracks that start after the shift
        finalTracks = [];
        startTime = [];
        endTime = [];
        shiftTime_hours = shiftTime/3600;
        for tr = 1:length(tracks_unique)
            
            currentTR = tracks_unique(tr);
            trackData = fullData(tracks == currentTR,:);
            trackTime = getGrowthParameter(trackData,'timestamp')./3600;
            
            if trackTime(1) > shiftTime_hours || trackTime(end) < shiftTime_hours
                continue
            else
                startTime = [startTime; trackTime(1)];
                endTime = [endTime; trackTime(end)];
                finalTracks = [finalTracks; currentTR];
            end
            
        end
       
        
        
        % iii. to eliminate cell cycles prior to shift when plotting
        birthEvents = getGrowthParameter(fullData,'isDrop');
        timestamps_hours = getGrowthParameter(fullData,'timestamp')./3600;
        divisionEvents = [birthEvents(2:end);1];
        divisionTimes = timestamps_hours(divisionEvents == 1);
        divisionTimes_postshift = divisionTimes > shiftTime_hours;
        
        divisionData = fullData(divisionEvents == 1,:);
        divisionData_postshift = divisionData(divisionTimes > shiftTime_hours,:);
        postshift_cells = getGrowthParameter(divisionData_postshift,'curveFinder');
        postshift_tracks = getGrowthParameter(divisionData_postshift,'trackNum');
        clear birthEvents timestamps_hours
        
        
        
        
        % 8. for each track that hits nutrient shift... 
        plotsPerFig = 60;
        figs = ceil(length(finalTracks)/plotsPerFig);
        fig_vector = [];
        sp_vector = [];
        fig_inc = ones(plotsPerFig,1);
        sp = linspace(1,plotsPerFig,plotsPerFig);
        for ff = 1:figs
            cf = fig_inc * ff;
            fig_vector = [fig_vector; cf];
            sp_vector = [sp_vector; sp'];
        end
        
        
        for tr = 1:600%length(finalTracks)
            
            
            % 8. isolate track data
            currentTrack = finalTracks(tr);
            trackNum = getGrowthParameter(fullData,'trackNum');
            trackData = fullData(trackNum == currentTrack,:);
            
            
            % 9. trim data BEFORE cell cycle with nutrient shift
            currentCycles = postshift_cells(postshift_tracks == currentTrack);
            firstCycle = currentCycles(1);
            curveIDs_tr = getGrowthParameter(trackData,'curveFinder');
            trackData_trimmed = trackData(curveIDs_tr >= firstCycle,:);
            
            
            % 10. plot track volume vs time
            trackVols = getGrowthParameter(trackData_trimmed,'volume');
            trackTimes = getGrowthParameter(trackData_trimmed,'timestamp')./3600;
            figure(fig_vector(tr))
            subplot(plotsPerFig/3,3,sp_vector(tr))
            plot(trackTimes,trackVols)
            
            
            % 11. plot shift time
            y = linspace(0,14,10);
            x = ones(10,1)*shiftTime_hours;
            figure(fig_vector(tr))
            subplot(plotsPerFig/3,3,sp_vector(tr))
            hold on
            plot(x,y,'Color',rgb('Tomato'))
            %title(strcat('track ',num2str(currentTrack)))
            
            
        end
        
    end
    
end



