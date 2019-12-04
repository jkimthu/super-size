%% Figure 1A: mean birth size vs mean growth rate



%  Goal: plot population-level birth size as a function of mean growth rate
%        across the cell cycle
%
%        steady-state populations have been demonstrated to follow a strong
%        positive correlation between the natural log of cell size and
%        growth rate.
%
%        do populations under fluctuations follow the same trend as
%        steady-state populations?
%
%        input: image processed data from each experimental replicate
%        output: plot, modeled after Figure 1C of Taheri et al., (2014)




%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell birth size and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line


%  Last edit: Jen Nguyen, 2019 December 4
%  Commit: first commit, birth size vs lambda - population mean

%  OK let's go!

%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


% 0. create array of experiments to use in analysis
exptArray = [2,3,4,5,6,7,9,10,11,12,13,14,15]; % use corresponding dataIndex values


% 0. initialize data vectors to store stats for each experiment
compiled_data = cell(length(exptArray),1);
compiled_mu = cell(length(exptArray),1);


%% Part 2. collect single cell data from all experiments
%          birth size and instantaneous growth rates


%  Strategy:
%
%       1.  loop through each experiment to collect data
%               2. initialize experiment meta data
%               3. load measured experiment data    
%               4. compile experiment data matrix
%               5. for each condition in experiment...
%                       5. isolate condition specific data
%                       6. isolate volume (Va), timestamp, drop, curve, and trackNum data
%                       7. calculate growth rate
%                       8. trim condition and growth rate data to include only full cell cycles
%                       9. isolate isDrop, timestamp, and timeSinceBirth data
%                      10. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
%                      11. remove zeros, which occur if no full track data exists at a drop
%                      12. truncate data to non-erroneous (e.g. bubbles) timestamps
%                      13. truncate data to stabilized regions
%                      14. if no div data in steady-state, skip condition
%                          else, trim outliers (those 3 std dev away from median) from final dataset
%                      15. bin growth rates by cell cycle, to match organization of birth size data
%                      16. store condition data into one variable per experiment
%               17. store experiment data into single variable for further analysis
%      18. save hard earned data


% 1. loop through each experiment to collect data
for e = 1:length(exptArray)

    
    
    % 2. initialize experiment meta data
    index = exptArray(e);
    date = storedMetaData{index}.date;
    expType = storedMetaData{index}.experimentType;
    bubbletime = storedMetaData{index}.bubbletime;
    xys = storedMetaData{index}.xys;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    
    
    
    % 3. initialize experiment specific variables
    ccData = cell(length(bubbletime),1);
    mu_instantaneous = cell(length(bubbletime),1);
    
    
    
    % 4. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    

    
    % for each condition in experiment
    for condition = 1:length(bubbletime)
            
            
        % 5. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        xy_start = min(xys(condition,:));
        xy_end = max(xys(condition,:));
        conditionData = buildDM(D5, T, xy_start, xy_end,index,expType);
        clear xy_start xy_end date
        
        
        
        % 6. calculate growth rate before trimming
        %     i) isolate required parameters
        volumes = getGrowthParameter(conditionData,'volume');            % calculated va_vals (cubic um)
        timestamps_sec = getGrowthParameter(conditionData,'timestamp');  % timestamp in seconds
        isDrop = getGrowthParameter(conditionData,'isDrop');             % isDrop, 1 marks a birth event
        curveFinder = getGrowthParameter(conditionData,'curveFinder');   % curve finder (ID of curve in condition)
        trackNum = getGrowthParameter(conditionData,'trackNum'); 
        
        %   ii) perform growht rate calculation
        growthRates_all = calculateGrowthRate(volumes,timestamps_sec,isDrop,curveFinder,trackNum);
        growthRates = growthRates_all(:,specificColumn);
        clear volumes isDrop trackNum timestamps_sec
        
        
          
        % 7. trim condition and growth rate data to include only full cell cycles
        conditionData_fullOnly = conditionData(curveFinder > 0,:);
        growthRates_fullOnly = growthRates(curveFinder > 0,:);
        curveIDs_fullOnly = curveFinder(curveFinder > 0,:);   % for trimming growth rate
        curveIDs_unique = unique(curveIDs_fullOnly);          % for assigning birth sizes
        trackNums = getGrowthParameter(conditionData_fullOnly,'trackNum');         % track number (not ID from particle tracking)
        clear curveFinder growthRates growthRates_all
        
             
        
        % 8. isolate timestamp, isDrop, width, length, surface area and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_fullOnly,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_fullOnly,'isDrop');      % isDrop, 1 marks a birth event
        volumes = getGrowthParameter(conditionData_fullOnly,'volume');     % calculated va_vals (cubic um)
        majorAxis = getGrowthParameter(conditionData_fullOnly,'length');   % length (um)
        minorAxis = getGrowthParameter(conditionData_fullOnly,'width');    % width (um) 
        sa = getGrowthParameter(conditionData_fullOnly,'surfaceArea');     % surface area (um^2)
        curveDuration = getGrowthParameter(conditionData_fullOnly,'curveDurations');  % length of cell cycle
        clear timestamps
        
        
        
        % 9. extract only final timeSinceBirth from each growth curve, this is the inter-division time!
        final_birthVolume = volumes(isDrop==1);
        final_birthLength = majorAxis(isDrop==1);
        final_birthWidth = minorAxis(isDrop==1);
        final_birthSA = sa(isDrop==1);
        final_birthSA2V = final_birthSA./final_birthVolume;
        finalTimestamps = timestamps_hr(isDrop==1); % experiment timestamp (hours) of each division event.
        final_durations = curveDuration(isDrop==1);
        final_trackNums = trackNums(isDrop==1);
        clear conditionData_fullOnly
        
        
        
        % 10. remove zeros, which occur if no full track data exists at a drop
        Vbirth = final_birthVolume(final_birthVolume > 0);
        Lbirth = final_birthLength(final_birthVolume > 0);
        Wbirth = final_birthWidth(final_birthVolume > 0);
        SA2Vbirth = final_birthSA2V(final_birthVolume > 0);
        birthTimestamps = finalTimestamps(final_birthVolume > 0);
        curveIDs = curveIDs_unique(final_birthVolume > 0);
        durations = final_durations(final_birthVolume > 0);
        tracks = final_trackNums(final_birthVolume > 0);
        clear final_birthVolume final_birthLength final_birthWidth finalTimestamps final_durations
        clear volumes majorAxis minorAxis sa final_birthSA2V final_birthSA curveDuration
        clear final_trackNums
        
        
        % 11. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        data = [Vbirth Lbirth Wbirth SA2Vbirth birthTimestamps curveIDs durations tracks];
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            data_bubbleTrimmed = data(birthTimestamps <= maxTime,:);
            birthTimestamps_bubbleTrimmed = birthTimestamps(birthTimestamps <= maxTime,:);

        else
            data_bubbleTrimmed = data;
            birthTimestamps_bubbleTrimmed = birthTimestamps;
            
        end
        clear timestamps_hr maxTime birthTimestamps data
        clear isDrop Vbirth Lbirth curveIDs SA2Vbirth Wbirth durations
        
        
        
        % 12. truncate data to stabilized regions
        minTime = 3;
        data_fullyTrimmed = data_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);       
        clear data_bubbleTrimmed birthTimestamps_bubbleTrimmed minTime
        
        
        
        % 13. isolate size from time and cell cycle information, in
        %     preparation to trim by outliers based on cell volume
        data_size = data_fullyTrimmed(:,1:4); % columns 1-4 = vol, length, width, SA:V
        data_Vbirth = data_fullyTrimmed(:,1);
        data_curves = data_fullyTrimmed(:,6);
        data_tau = data_fullyTrimmed(:,7);
        data_trackNum = data_fullyTrimmed(:,8);
        
        
        % 14. if no div data in steady-state, skip condition
        if isempty(data_curves) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            vol_median = median(data_Vbirth);
            vol_std_temp = std(data_Vbirth);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            sizes_temp = data_size(data_Vbirth <= (vol_median+vol_std_temp*3),:); % cut largest vals, over 3 std out
            IDs_temp = data_curves(data_Vbirth <= (vol_median+vol_std_temp*3));
            tau_temp = data_tau(data_Vbirth <= (vol_median+vol_std_temp*3));
            track_temp = data_trackNum(data_Vbirth <= (vol_median+vol_std_temp*3));
            vol_temp = sizes_temp(:,1);
            clear data_curves data_Vbirth 
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            sizes_final = sizes_temp(vol_temp >= (vol_median-vol_std_temp*3),:);          % cut smallest vals, over 3 std out 
            IDs_final = IDs_temp(vol_temp >= (vol_median-vol_std_temp*3));   
            tau_final = tau_temp(vol_temp >= (vol_median-vol_std_temp*3));
            trackNum_final = track_temp(vol_temp >= (vol_median-vol_std_temp*3));
            clear vol_median vol_std_temp sizes_temp IDs_temp vol_temp tau_temp track_temp
            
            % iv. remove corresponding growth rates from datasets
            trimmedIDs = setdiff(curveIDs_unique,IDs_final);    % curve IDs in growth rate dataset, NOT in final IDs trimmed by cell cycle
            toTrim = ismember(curveIDs_fullOnly,trimmedIDs);   % vector of what to trim or not in growth rate
            trimmed_curves_insta = curveIDs_fullOnly(toTrim == 0);
            trimmed_mus = growthRates_fullOnly(toTrim == 0);
            clear toTrim trimmedIDs curveIDs_fullOnly growthRates_fullOnly
            
 
            
            % 16. bin growth rates by cell cycle, to match organization of birth size data
            mus_binned = accumarray(trimmed_curves_insta,trimmed_mus,[],@(x) {x});
            mus = mus_binned(~cellfun('isempty',mus_binned));
            lambdas = cellfun(@nanmean,mus);
            clear trimmed_curves_insta trimmed_mus mus_binned
            
            
            % 17. determine birth size +1 for final full cell cycle of each track
            %  ** birth size +1 is the birth size immediately AFTER the full cycle **
            tracks_all = getGrowthParameter(conditionData,'trackNum'); % track number
            tracks_unique = unique(trackNum_final);
            
            plus1_volume = nan(length(trackNum_final),1);
            plus1_length = nan(length(trackNum_final),1);
            plus1_width = nan(length(trackNum_final),1);
            plus1_sa = nan(length(trackNum_final),1);
            
            counter = 0;
            for ut = 1:length(tracks_unique)
                
                % i. for each track, loop through all full cell cycles and
                %    determine next birth size for each full cycle
                
                %    A) identify current track and its full cell cycles
                currentTrack = tracks_unique(ut);
                cycles = IDs_final(trackNum_final == currentTrack);
                
                %    B) isolate track data including incomplete cell cycles
                %   ** will have some cell cycles not in final dataset, due
                %      to outlier removal**
                currentTrack_data = conditionData(tracks_all == currentTrack,:);
                
                %    C) loop through cell cycles of current track
                for uc = 1:length(cycles)
                    
                    counter = counter + 1;
                    currentCycle = cycles(uc);
                    %last_cycle = cycles(end);
                    
                    
                    % D) locate birth size of next cycle from not yet trimmed data
                   
                    %   i. determine final row of current cycle
                    currentTrack_cycles = getGrowthParameter(currentTrack_data,'curveFinder');
                    currentCycle_rows = find(currentTrack_cycles == currentCycle);
                    lastRow =currentCycle_rows(end);
        
                    %    ii) save first data points after final full cycle row
                    data_plus1 = currentTrack_data(lastRow+1,:);
                    
                    plus1_volume(counter) = getGrowthParameter(data_plus1,'volume');     % calculated va_vals (cubic um)
                    plus1_length(counter) = getGrowthParameter(data_plus1,'length');   % length (um)
                    plus1_width(counter) = getGrowthParameter(data_plus1,'width');    % width (um)
                    plus1_sa(counter) = getGrowthParameter(data_plus1,'surfaceArea');     % surface area (um^2)
                    
                end
            end
            plus1_sa2v = plus1_sa./plus1_volume;
            
            
            % 18. store condition data into one variable per experiment
            cc_data = [IDs_final sizes_final(:,1) plus1_volume sizes_final(:,2) plus1_length sizes_final(:,3) plus1_width sizes_final(:,4) plus1_sa2v lambdas tau_final/60 trackNum_final]; % tau here is converted from sec to min
            
            ccData{condition} = cc_data;
            mu_instantaneous{condition} = mus;
        
            
        end
        clear data_size data_fullyTrimmed data_tau
        
    end
      
    % 17. store experiment data into single variable for further analysis
    compiled_data{e} = ccData;
    compiled_mu{e} = mu_instantaneous;
    
    clear ccData mu_instantaneous
end


% 18. save hard earned data
cd('/Users/jen/super-size/')
save('A1_ccSize.mat','compiled_data','compiled_mu','exptArray')


%% Part 3. sort data by nutrient condition

% goal: plot of newborn cell size vs growth rate,
%       mean of each condition replicate

% strategy: 
%
%       i. accumulate data from each condition
%          conditions: each steady and each fluctuating timescale (7 total)
%
%      ii. calculate mean birth size and growth rate for each condition (population data)
%          plot each condition as a closed orange point and fit a line (the growth law ACROSS conditions)
%
%     iii. bin individual data by growth rate
%
%      iv. calculate mean birth size and growth rate for each bin (individual data)
%          plot each bin as an open blue point and fit a line WITHIN each condition


clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_ccSize.mat')
lamb = 10; % column in compiled_data that is lambda
sm = 2:9; % columns in compiled_data that are sizes

% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};
shape = 'o';


% 1. accumulate data from each condition
fluc = 1; % row number in data structure
low = 2; 
ave = 3; 
high = 4;


sigmas = 3;
for ee = 1:length(environment_order)
    
    condition = environment_order{ee};
    
    if ischar(condition) == 1
        
        % steady environment! concatenate data based on nutrient level
        if strcmp(condition,'low') == 1
            
            lambda_low = [];
            birthSizes_low = [];

            % loop through all experiments and store low data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{low,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{low,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{low,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_low = [lambda_low; expt_lambda];
                    birthSizes_low = [birthSizes_low; expt_sizes];
                    clear expt_lambda expt_sizes
                end
                
            end
            clear expt expt_data
            
            condition_lambda = lambda_low;
            condition_sizes = birthSizes_low;
            
            
            % isolate cycles within some error of replicate
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std * sigmas);
            upper = condition_lambda > condition_mean - (condition_std * sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            
            % store condition data
            sizes{1} = range_sizes;
            lambda{1} = range_lambda;
            
            
            
        elseif strcmp(condition,'ave') == 1
            
            lambda_ave = [];
            birthSizes_ave = [];
            
            % loop through all experiments and store ave data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{ave,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{ave,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{ave,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_ave = [lambda_ave; expt_lambda];
                    birthSizes_ave = [birthSizes_ave; expt_sizes];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data
         
            
            condition_lambda = lambda_ave;
            condition_sizes = birthSizes_ave;
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std * sigmas);
            upper = condition_lambda > condition_mean - (condition_std * sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data
            lambda{6} = range_lambda;
            sizes{6} = range_sizes;

            
        elseif strcmp(condition,'high') == 1
            
            lambda_high = [];
            birthSizes_high = [];
            
            % loop through all experiments and store high data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{high,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{high,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{high,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_high = [lambda_high; expt_lambda];
                    birthSizes_high = [birthSizes_high; expt_sizes];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data
            
            condition_lambda = lambda_high;
            condition_sizes = birthSizes_high;
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std * sigmas);
            upper = condition_lambda > condition_mean - (condition_std * sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data in cell corresponding to Condition Order
            lambda{7} = range_lambda;
            sizes{7} = range_sizes;
            
        end
    else
        
        % fluctuating environment! concatenate based on timescale
        if condition == 30
            idx = [2,3,4]; % ID of experiments with this fluc timescale
            
            lambda_30 = [];
            birthSizes_30 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{fluc,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{fluc,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_30 = [lambda_30; expt_lambda];
                    birthSizes_30 = [birthSizes_30; expt_sizes];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data
            
            condition_lambda = lambda_30;
            condition_sizes = birthSizes_30;
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std * sigmas);
            upper = condition_lambda > condition_mean - (condition_std * sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data in cell corresponding to Condition Order
            lambda{2} = range_lambda;
            sizes{2} = range_sizes;
            
            
            
        elseif condition == 300
            idx = [5,6,7]; % ID of experimennts with this fluc timescale
            
            lambda_300 = [];
            birthSizes_300 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{fluc,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{fluc,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_300 = [lambda_300; expt_lambda];
                    birthSizes_300 = [birthSizes_300; expt_sizes];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data
            
            condition_lambda = lambda_300;
            condition_sizes = birthSizes_300;
            
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std * sigmas);
            upper = condition_lambda > condition_mean - (condition_std * sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data in cell corresponding to Condition Order
            lambda{3} = range_lambda;
            sizes{3} = range_sizes;
            
            
        elseif condition == 900
            idx = [9,10,11,12]; % ID of experimennts with this fluc timescale
            
            lambda_900 = [];
            birthSizes_900 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{fluc,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{fluc,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_900 = [lambda_900; expt_lambda];
                    birthSizes_900 = [birthSizes_900; expt_sizes];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data
            
            condition_lambda = lambda_900;
            condition_sizes = birthSizes_900;
            
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std*sigmas);
            upper = condition_lambda > condition_mean - (condition_std*sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data in cell corresponding to Condition Order
            lambda{4} = range_lambda;
            sizes{4} = range_sizes;
            
            
            
        elseif condition == 3600
            idx = [13,14,15]; % ID of experimennts with this fluc timescale
            
            lambda_3600 = [];
            birthSizes_3600 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = compiled_data{expt,1}{fluc,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = compiled_data{expt,1}{fluc,1}(:,sm);
                    
                    % concanetate individual cell cycle values
                    lambda_3600 = [lambda_3600; expt_lambda];
                    birthSizes_3600 = [birthSizes_3600; expt_sizes];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data
            
            condition_lambda = lambda_3600;
            condition_sizes = birthSizes_3600;
            
            
            % isolate cycles within 3 st dev of mean
            condition_mean = nanmean(condition_lambda);
            condition_std = nanstd(condition_lambda);
            
            lower = condition_lambda < condition_mean + (condition_std*sigmas);
            upper = condition_lambda > condition_mean - (condition_std*sigmas);
            combined = lower + upper;
            
            range_lambda = condition_lambda(combined == 2);
            range_sizes = condition_sizes(combined == 2,:);
            clear condition_mean condition_std lower upper
            
            % store condition data in cell corresponding to Condition Order
            lambda{5} = range_lambda;
            sizes{5} = range_sizes;
            
        end  
    end
end
clear fluc low ave high idx condition ee expt arrayIndex
clear condition_lambda condition_sizes
clear range_lambda range_sizes combined
clear birthSizes_30 birthSizes_300 birthSizes_900 birthSizes_3600
clear birthSizes_low birthSizes_ave birthSizes_high
clear lambda_30 lambda_300 lambda_900 lambda_3600 lambda_low lambda_ave lambda_high


%% Part 4. plot cell size at birth vs mean growth rate (lambda) and best fit line from steady points

% 1. calculate mean birth size and growth rate for each condition (population data)
pop_size = cellfun(@mean,sizes,'UniformOutput',false);
pop_lambda = cellfun(@mean,lambda);


% 2. plot each condition as a closed point and fit a line (the growth law ACROSS conditions)
for metric = 1:length(sm)
    
    figure(metric)
    steady_sizes = nan(1,3);
    steady_counter = 0;
    for cond = 1:length(pop_lambda)
        
        if ischar(environment_order{cond}) == 1
            steady_counter = steady_counter + 1;
            steady_sizes(steady_counter) = pop_size{cond}(metric);
        end
        
        color = rgb(palette(cond));
        plot(pop_lambda(cond),log(pop_size{cond}(metric)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
        hold on
    end

    
    % 3. determine best fit line
    steady_lambda = [pop_lambda(1), pop_lambda(6), pop_lambda(7)];
    fit = polyfit(steady_lambda,log(steady_sizes),1);
    
    % 4. plot best fit line
    x = linspace(0,4,10);
    y = fit(1)*x + fit(2);
    
    hold on
    plot(x,y,'Color',rgb('SlateGray'))
    
    xlabel('lambda')
    xlim([0.5 3.5])
    legend('low','30','300','900','3600','ave','high')
end
clear color cc y x fit metric cond

figure(1)
ylabel('Vbirth')

figure(2)
ylabel('Vbirth+1')

figure(3)
ylabel('Lbirth')

figure(4)
ylabel('Lbirth+1')

figure(5)
ylabel('Wbirth')

figure(6)
ylabel('Wbirth+1')

figure(7)
ylabel('SA:V at birth')

figure(8)
ylabel('SA:V at birth+1')


%% Part 5. statistical analysis (std, sem, cv)

% 1. calculate standard deviation of each condition
pop_size_std = cellfun(@std,sizes,'UniformOutput',false);
pop_lambda_std = cellfun(@std,lambda);



% 2. calculate standard error of the mean (sem)
pop_size_count = cellfun(@length,sizes,'UniformOutput',false);
for condition = 1:length(lambda)
    pop_size_sem{condition} = pop_size_std{condition}./sqrt(pop_size_count{condition});
end
clear condition

pop_lambda_count = cellfun(@length,lambda);
pop_lambda_sem = pop_lambda_std./sqrt(pop_lambda_count);



% 3. calculate CV for each condition
for condi = 1:length(lambda)
    pop_size_cv{condi} = pop_size_std{condi}./pop_size{condi} * 100;
end
clear condi

pop_lambda_cv = pop_lambda_std./pop_lambda * 100;

