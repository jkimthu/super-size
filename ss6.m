%% ss6 - probability of division as a function of nutrient score


% Goal: from cells experiencing single nutrient downshift, plot probability
%       of division as a function of nutrient score


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
 

%% Part 1. 

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
        
        
        % 6. for each track, determine whether nutrient shift occurs
        tracks = getGrowthParameter(fullData,'trackNum');
        tracks_unique = unique(tracks);
        
        % i. eliminate tracks that start after the shift or end before the shift
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
        clear tr currentTR trackTime trackData tracks_unique
        
        
        % 7. for each remaining track, determine which cell cycle encompasses the shift
        %cnExp = cell(length(finalTracks),1);
        %ccData = nan(length(finalTracks),4); % 1. curveID, 2. birth vol,
                                             % 3. tau, 4. integrated experience 
        
        counter = 0;
        for trac = 1:length(finalTracks)
            
            % i. isolate track data
            currentTrack = finalTracks(trac);
            trackNum = getGrowthParameter(fullData,'trackNum');
            trackData = fullData(trackNum == currentTrack,:);
            trackTimes = getGrowthParameter(trackData,'timestamp')./3600;
            trackCells = getGrowthParameter(trackData,'curveFinder');
            trackVols = getGrowthParameter(trackData,'volume');
            
            % ii. determine number of cell cycles in track 
            cells = getGrowthParameter(trackData,'curveFinder');
            cells_unique = unique(cells);
            
            % iii. loop through cells until finding cell encountering shift
            for cc = 1:length(cells_unique)
                
                currentCC = cells_unique(cc);
                ccTimes = trackTimes(trackCells == currentCC);
                ccVols = trackVols(trackCells == currentCC);
                birthTime = ccTimes(1);
                divTime = ccTimes(end);
                
                if birthTime < shiftTime_hours && divTime > shiftTime_hours
                    
                    counter = counter + 1;
                    
                    % if cell experiences shift
                    % 1. define nutrient experience
                    nExperience = ccTimes < shiftTime_hours;
                    nvector = nan(length(nExperience),1);
                    nvector(nExperience == 0) = 0.01; % percent LB in low
                    nvector(nExperience == 1) = 2;    % percent LB in low
                    
                    % 2. integrate nutrient experience
                    iN = sum(nvector);
                    
                    % 3. store stats
                    currentStats = [currentCC, ccVols(1), divTime-birthTime, iN]; % 1. curveID, 2. birth vol, 3. tau, 4. integrated experience
                    ccData(counter,:) = currentStats;
                    cnExp{counter} = nvector;
                    
                    % end loop through cell cycles once shift is found
                    break
                    
                else
                    continue
                end
                
            end
            clear nExperience currentStats currentCC iN ccTimes
            clear divTime birthTime ccVols
            
        end
        
        
        % 8. bin cell cycles by integrated nutrient experience
        in = ccData(:,4);
        tau = ccData(:,3);
        inbin = 0.05;
        bins = ceil(in/inbin);
        tau_binned = accumarray(bins,tau,[],@(x) {x});
       

        % 9. plot cdf over integrated nutrient score
        counts_in = cellfun(@length,tau_binned);
        tau_in = cellfun(@mean,tau_binned);
        tau_std = cellfun(@std,tau_binned);
        in_vector = linspace(1,max(bins),max(bins)) * inbin;
        
        % 9. plot cdf over integrated nutrient score
        test = ~isnan(tau_in);
        tau_in = tau_in(test);
        tau_std= tau_std(test);
        
        in_vec = in_vector(test);
        counts_in = counts_in(test);
        total_counts = sum(counts_in);
        frac_in = counts_in./total_counts;
        
        figure(1)
        plot(in_vec,tau_in)
        hold on
        errorbar(in_vec,tau_in,tau_std,'o') %(in_vec,tau_std,'o')
        xlabel('integrated nutrient score')
        ylabel('interdivision time (h)')
        
        figure(2)
        cdf(frac_in)
        
        
    end
    
end





