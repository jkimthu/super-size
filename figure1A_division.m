%% Figure 1A_division: mean division size vs mean growth rate



%  Goal: plot population-level division size as a function of mean growth rate
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


%  Last edit: Jen Nguyen, 2019 Jan 16
%  Commit: plot size at division vs growth rate

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


    
    % 4. load measured experiment data    
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
    cd(experimentFolder)
    filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    load(filename,'D5','T');
    clear filename experimentFolder


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
        clear curveFinder growthRates growthRates_all
        
           
        
        % 8. remove cell cycles that are only 2 timepoints or less long
        [instances,id_num] = hist(curveIDs_fullOnly,curveIDs_unique);
        tooShort = find(instances < 3);
        cd_temp = conditionData_fullOnly;
        curveIDs_temp = curveIDs_fullOnly;
        for ii = 1:length(tooShort)
            current_shortie = tooShort(ii); % shortie # = curveID_unique #
            cd_temp = cd_temp(curveIDs_temp ~= current_shortie,:);
            curveIDs_temp = getGrowthParameter(cd_temp,'curveFinder');
        end
        conditionData_final = cd_temp;
        clear ii cd_temp curveIDs_temp current_shortie instances id_num tooShort
       
        
        
        % 9. isolate timestamp, isDrop, width, length, surface area and volume data for cell cycle measurements
        timestamps = getGrowthParameter(conditionData_final,'timestamp');  % timestamp in seconds
        timestamps_hr = timestamps./3600;    % convert timestamp to hours
        isDrop = getGrowthParameter(conditionData_final,'isDrop');      % isDrop, 1 marks a birth event
        volumes = getGrowthParameter(conditionData_final,'volume');     % calculated va_vals (cubic um)
        majorAxis = getGrowthParameter(conditionData_final,'length');   % length (um)
        minorAxis = getGrowthParameter(conditionData_final,'width');    % width (um) 
        sa = getGrowthParameter(conditionData_final,'surfaceArea');     % surface area (um^2)
        curveDuration = getGrowthParameter(conditionData_final,'curveDurations');  % length of cell cycle
        trackNums = getGrowthParameter(conditionData_final,'trackNum');         % track number (not ID from particle tracking)
        curveIDs_final = getGrowthParameter(conditionData_final,'curveFinder');   % curve finder (ID of curve in condition)
        clear timestamps
        
        
        
        % 10. extract data associated with births, divisions, and births +1
        rows_births = find(isDrop==1);
        
        r_div_temp = rows_births - 1;
        r_div_temp2 = r_div_temp(2:end);
        rows_divisions = [r_div_temp2; length(conditionData_final)];
        clear r_div_temp r_div_temp2
        
        Vbirth = volumes(rows_births);
        Lbirth = majorAxis(rows_births);
        Wbirth = minorAxis(rows_births);
        SAbirth = sa(rows_births);
        SA2Vbirth = SAbirth./Vbirth;
        Tbirth = timestamps_hr(rows_births); % experiment timestamp (hours) of each division event.
        
        Vdiv = volumes(rows_divisions);
        Ldiv = majorAxis(rows_divisions);
        Wdiv = minorAxis(rows_divisions);
        SAdiv = sa(rows_divisions);
        SA2Vdiv = SAdiv./Vdiv;
        Tdiv = timestamps_hr(rows_divisions);
        
        durations = curveDuration(rows_births);
        tracks = trackNums(rows_births);
        curveIDs = curveIDs_final(rows_births);
        clear conditionData_fullOnly SAbirth SAdiv trackNums
        clear volumes majorAxis minorAxis sa final_birthSA2V final_birthSA curveDuration
        clear curveIDs_final
        
        
        % 11. truncate data to non-erroneous (e.g. bubbles) timestamps
        %     Note: trimming first by coursest time resolution, which is for the cell cycle.
        %           later we will trim all growth rate data that are not associated with cell cycles remaining in analysis
        data = [Vbirth Vdiv Lbirth Ldiv Wbirth Wdiv SA2Vbirth SA2Vdiv Tbirth Tdiv curveIDs durations tracks];
        maxTime = bubbletime(condition);
        
        if maxTime > 0

            data_bubbleTrimmed = data(Tbirth <= maxTime,:);
            birthTimestamps_bubbleTrimmed = Tbirth(Tbirth <= maxTime,:);

        else
            data_bubbleTrimmed = data;
            birthTimestamps_bubbleTrimmed = Tbirth;
            
        end
        clear timestamps_hr maxTime Tbirth data
        clear isDrop Vbirth Lbirth curveIDs SA2Vbirth Wbirth durations tracks
        clear Vdiv Ldiv Wdiv SA2Vdiv Tdiv
        clear rows_births rows_divisions
        
        
        
        % 12. truncate data to stabilized regions
        minTime = 3;
        data_fullyTrimmed = data_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);       
        clear data_bubbleTrimmed birthTimestamps_bubbleTrimmed minTime
        
        
        
        % 13. isolate size from time and cell cycle information, in
        %     preparation to trim by outliers based on cell volume
        data_all = data_fullyTrimmed;
        data_Vbirth = data_fullyTrimmed(:,1);
        
        
        
        % 14. if no div data in steady-state, skip condition
        if isempty(data_Vbirth) == 1
            continue
        else
            
            % 15. trim outliers (those 3 std dev away from median) from final dataset
            
            % i. determine median and standard deviation of birth size
            vol_median = median(data_Vbirth);
            vol_std_temp = std(data_Vbirth);
            
            % ii. remove cell cycles of WAY LARGE birth size, tracking IDs
            all_temp = data_all(data_Vbirth <= (vol_median+vol_std_temp*3),:);
            vol_temp = all_temp(:,1);
            clear data_curves data_Vbirth 
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            all_final = all_temp(vol_temp >= (vol_median-vol_std_temp*3),:); 
            clear vol_median vol_std_temp vol_temp all_temp
            clear data_all data_fullyTrimmed
            
            % iv. remove corresponding growth rates from datasets
            IDs_final = all_final(:,11); % new
            trimmedIDs = setdiff(curveIDs_unique,IDs_final);    % curve IDs in growth rate dataset, NOT in final IDs trimmed by cell cycle
            toTrim = ismember(curveIDs_fullOnly,trimmedIDs);   % vector of what to trim or not in growth rate
            trimmed_curves_insta = curveIDs_fullOnly(toTrim == 0);
            trimmed_mus = growthRates_fullOnly(toTrim == 0);
            clear toTrim trimmedIDs curveIDs_fullOnly growthRates_fullOnly
            clear curveIDs_unique
 
            
            % 16. bin growth rates by cell cycle, to match organization of birth size data
            mus_binned = accumarray(trimmed_curves_insta,trimmed_mus,[],@(x) {x});
            mus = mus_binned(~cellfun('isempty',mus_binned));
            lambdas = cellfun(@nanmean,mus);
            clear trimmed_curves_insta trimmed_mus mus_binned
            
            
            
            % 17. determine birth size +1 for final full cell cycle of each track
            %  ** birth size +1 is the birth size immediately AFTER the full cycle **
            tracks_total = getGrowthParameter(conditionData,'trackNum'); % track number
            tracks_all = all_final(:,13); % column 13 = track numbers
            tracks_unique = unique(tracks_all);
            
            plus1_volume = nan(length(tracks_all),1);
            plus1_length = nan(length(tracks_all),1);
            plus1_width = nan(length(tracks_all),1);
            plus1_sa = nan(length(tracks_all),1);
            
            counter = 0;
            for ut = 1:length(tracks_unique)
                
                % i. for each track, loop through all full cell cycles and
                %    determine next birth size for each full cycle
                
                %    A) identify current track and its full cell cycles
                currentTrack = tracks_unique(ut);
                cycles = IDs_final(tracks_all == currentTrack);
                
                %    B) isolate track data including incomplete cell cycles
                %   ** will have some cell cycles not in final dataset, due
                %      to outlier removal**
                currentTrack_data = conditionData(tracks_total == currentTrack,:);
                
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
            clear uc ut counter currentTrack_data currentCycle currentCycle_rows
            clear currentTrack_cycles lastRow tracks_total tracks_unique currentTrack
            clear cycles plus1_sa data_plus1
            
            
            % 18. organize fully processed data
            %     size data for next birth
            plus1_data = [plus1_volume plus1_length plus1_width plus1_sa2v]; 
            clear plus1_volume plus1_length plus1_width plus1_sa2v
            
            %     mean growth rate, interdivision time, track #, birth time, div time
            tau_all = all_final(:,12);   % column 12 = cell cycle duration in seconds
            Tbirth_all = all_final(:,9); % column 9 = time at birth in hours
            Tdiv_all = all_final(:,10);  % column 10 = time at division in hours
            meta_data = [lambdas tau_all/60 tracks_all Tbirth_all Tdiv_all];% tau here is converted from sec to min
            clear lambdas tau_all tracks_all Tbirth_all Tdiv_all
            
            %     size data for birth and division of current cell cycle
            cc_data = all_final(:,1:8);
            
            
            %     compile data for this condition into a data structure
            ccData{condition}.cc = cc_data;
            ccData{condition}.meta = meta_data;
            ccData{condition}.plus1 = plus1_data;
            ccData{condition}.mu_insta = mus;
        
 
        end
        clear mus plus1_data meta_data cc_data all_final IDs_final 
        
    end
      
    % 17. store experiment data into single variable for further analysis
    compiled_data{e} = ccData;
    clear ccData
    
end


% 18. save hard earned data
cd('/Users/jen/super-size/')
save('A1_div_ccSize.mat','compiled_data','exptArray')


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
load('A1_div_ccSize.mat')
lamb = 1; % column in compiled_data (meta) that is lambda


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
            ccSizes_low = [];
            plusSizes_low = [];

            % loop through all experiments and store low data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{low,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_low = [lambda_low; expt_lambda];
                    ccSizes_low = [ccSizes_low; expt_sizes];
                    plusSizes_low = [plusSizes_low; expt_plus]; 
                    clear expt_lambda expt_sizes expt_plus
                end
                
            end 
            clear expt expt_data expt_plus
            
            sizes{1} = ccSizes_low;
            lambda{1} = lambda_low;
            plus{1} = plusSizes_low;
            
        elseif strcmp(condition,'ave') == 1
            
            lambda_ave = [];
            ccSizes_ave = [];
            plusSizes_ave = [];
            
            % loop through all experiments and store ave data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{ave,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_ave = [lambda_ave; expt_lambda];
                    ccSizes_ave = [ccSizes_ave; expt_sizes];
                    plusSizes_ave = [plusSizes_ave; expt_plus];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data expt_plus
 
            
            % store condition data
            lambda{6} = lambda_ave; %range_lambda;
            sizes{6} = ccSizes_ave; %range_sizes;
            plus{6} = plusSizes_ave;

            
        elseif strcmp(condition,'high') == 1
            
            lambda_high = [];
            ccSizes_high = [];
            plusSizes_high = [];
            
            % loop through all experiments and store high data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{high,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_high = [lambda_high; expt_lambda];
                    ccSizes_high = [ccSizes_high; expt_sizes];
                    plusSizes_high = [plusSizes_high; expt_plus];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            lambda{7} = lambda_high; 
            sizes{7} = ccSizes_high; 
            plus{7} = plusSizes_high;
            
        end
    else
        
        % fluctuating environment! concatenate based on timescale
        if condition == 30
            idx = [2,3,4]; % ID of experiments with this fluc timescale
            
            lambda_30 = [];
            ccSizes_30 = [];
            plusSizes_30 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_30 = [lambda_30; expt_lambda];
                    ccSizes_30 = [ccSizes_30; expt_sizes];
                    plusSizes_30 = [plusSizes_30; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            lambda{2} = lambda_30; 
            sizes{2} = ccSizes_30; 
            plus{2} = plusSizes_30;
            
            
        elseif condition == 300
            idx = [5,6,7]; % ID of experimennts with this fluc timescale
            
            lambda_300 = [];
            ccSizes_300 = [];
            plusSizes_300 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_300 = [lambda_300; expt_lambda];
                    ccSizes_300 = [ccSizes_300; expt_sizes];
                    plusSizes_300 = [plusSizes_300; expt_plus];
                    
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            lambda{3} = lambda_300; %range_lambda;
            sizes{3} = ccSizes_300; %range_sizes;
            plus{3} = plusSizes_300;
            
        elseif condition == 900
            idx = [9,10,11,12]; % ID of experimennts with this fluc timescale
            
            lambda_900 = [];
            ccSizes_900 = [];
            plusSizes_900 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_900 = [lambda_900; expt_lambda];
                    ccSizes_900 = [ccSizes_900; expt_sizes];
                    plusSizes_900 = [plusSizes_900; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            lambda{4} = lambda_900; 
            sizes{4} = ccSizes_900; 
            plus{4} = plusSizes_900;
            
            
        elseif condition == 3600
            idx = [13,14,15]; % ID of experimennts with this fluc timescale
            
            lambda_3600 = [];
            ccSizes_3600 = [];
            plusSizes_3600 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_lambda = expt_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    lambda_3600 = [lambda_3600; expt_lambda];
                    ccSizes_3600 = [ccSizes_3600; expt_sizes];
                    plusSizes_3600 = [plusSizes_3600; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            lambda{5} = lambda_3600; 
            sizes{5} = ccSizes_3600;
            plus{5} = plusSizes_3600;
            
        end  
    end
end
clear fluc low ave high idx condition ee expt arrayIndex
clear condition_lambda condition_sizes
clear range_lambda range_sizes combined
clear birthSizes_30 birthSizes_300 birthSizes_900 birthSizes_3600
clear birthSizes_low birthSizes_ave birthSizes_high
clear lambda_30 lambda_300 lambda_900 lambda_3600 lambda_low lambda_ave lambda_high
clear plusSizes_30 plusSizes_300 plusSizes_900 plusSizes_3600 plusSizes_low plusSizes_ave plusSizes_high
clear ccSizes_30 ccSizes_300 ccSizes_900 ccSizes_3600 ccSizes_low ccSizes_ave ccSizes_high


%% Part 4. plot cell size at birth/div vs mean growth rate (lambda) and best fit line from steady points

% 1. calculate mean birth size and growth rate for each condition (population data)
pop_size = cellfun(@mean,sizes,'UniformOutput',false);
pop_lambda = cellfun(@mean,lambda);

% 2. plot each condition as a closed point and fit a line (the growth law ACROSS conditions)
cc = 8;
for metric = 1:cc
    
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
ylabel('Vdiv')

figure(3)
ylabel('Lbirth')

figure(4)
ylabel('Ldiv')

figure(5)
ylabel('Wbirth')

figure(6)
ylabel('Wdiv')

figure(7)
ylabel('SA:V at birth')

figure(8)
ylabel('SA:V at division')

 
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


%% Part 6. plot size of next birth over vs mean growth rate and best fit

% 1. calculate mean birth size and growth rate for each condition (population data)
pop_lambda = cellfun(@mean,lambda);
pop_plus = cellfun(@mean,plus,'UniformOutput',false);

% 2. plot each condition as a closed point and fit a line (the growth law ACROSS conditions)
pm = 4;
for metric = 1:pm
    
    figure(metric)
    steady_plusSizes = nan(1,3);
    steady_counter = 0;
    for cond = 1:length(pop_lambda)
        
        if ischar(environment_order{cond}) == 1
            steady_counter = steady_counter + 1;
            steady_plusSizes(steady_counter) = pop_plus{cond}(metric);
        end
        
        color = rgb(palette(cond));
        plot(pop_lambda(cond),log(pop_plus{cond}(metric)),'Color',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
        hold on
    end

    
    % 3. determine best fit line
    steady_lambda = [pop_lambda(1), pop_lambda(6), pop_lambda(7)];
    fit = polyfit(steady_lambda,log(steady_plusSizes),1);
    
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
ylabel('Vbirth+1')

figure(2)
ylabel('Lbirth+1')

figure(3)
ylabel('Wbirth+1')

figure(4)
ylabel('SA:V at birth+1')
