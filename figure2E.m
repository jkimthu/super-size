%% Figure 2E: fluctuating motion relative to steady relationship


%  Goal: what is the path taken in fluctuatings relative
%        to the path between two steady-states?

%        plot birth size vs mean growth rate from data,
%        binned by period fraction (like 2C fused with 2A)


%  Strategy: 
%
%       1) initialize experimental data
%       2) collect single cell data from all experiments
%       3) sort data by nutrient condition
%       4) bin data by time at birth and plot



%  Last edit: jen, 2019 December 18
%  Commit: first commit, timestep: 5 min



%  OK let's go!


%% Part 1. initialize

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));


% 0. initialize array of experiments to use in analysis, then loop through each
exptArray = 37:39; % use corresponding dataIndex values


% 0. define method of calculating growth rate
specificGrowthRate = 'log2';
specificColumn = 3;             % for selecting appropriate column in growthRates


%% Part 2. collect single cell data from all experiments
%          birth size and instantaneous growth rates

%  Same as figure1A.m
%  except: added shiftTime as parameter used to trim early data (minTime)
%          & adjusted filename to convention of single shift experiments

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
    shiftTime = storedMetaData{index}.shiftTime;
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
        xy_start = min(xys{condition,:});
        xy_end = max(xys{condition,:});
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
        
        
        
        % 12. truncate data to time after stabilization at steady high
        minTime = 3;
        data_fullyTrimmed = data_bubbleTrimmed(birthTimestamps_bubbleTrimmed >= minTime,:);       
        clear data_bubbleTrimmed birthTimestamps_bubbleTrimmed minTime
        
        
        
        % 13. isolate size from time and cell cycle information, in
        %     preparation to trim by outliers based on cell volume
        data_size = data_fullyTrimmed(:,1:4); % columns 1-4 = vol, length, width, SA:V
        data_Vbirth = data_fullyTrimmed(:,1);
        data_timestamps = data_fullyTrimmed(:,5);
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
            time_temp = data_timestamps(data_Vbirth <= (vol_median+vol_std_temp*3),:);
            IDs_temp = data_curves(data_Vbirth <= (vol_median+vol_std_temp*3));
            tau_temp = data_tau(data_Vbirth <= (vol_median+vol_std_temp*3));
            track_temp = data_trackNum(data_Vbirth <= (vol_median+vol_std_temp*3));
            vol_temp = sizes_temp(:,1);
            clear data_curves data_Vbirth 
            
            % iii. remove cell cycle of WAY SMALL birth size, tracking IDs
            sizes_final = sizes_temp(vol_temp >= (vol_median-vol_std_temp*3),:);          % cut smallest vals, over 3 std out 
            times_final = time_temp(vol_temp >= (vol_median-vol_std_temp*3),:); 
            IDs_final = IDs_temp(vol_temp >= (vol_median-vol_std_temp*3));   
            tau_final = tau_temp(vol_temp >= (vol_median-vol_std_temp*3));
            trackNum_final = track_temp(vol_temp >= (vol_median-vol_std_temp*3));
            clear vol_median vol_std_temp sizes_temp IDs_temp vol_temp tau_temp track_temp time_temp
            
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
                clear uc
            end
            clear ut
            
            plus1_sa2v = plus1_sa./plus1_volume;
            
            
            % 18. store condition data into one variable per experiment
            cc_data = [IDs_final sizes_final(:,1) plus1_volume sizes_final(:,2) plus1_length sizes_final(:,3) plus1_width sizes_final(:,4) plus1_sa2v lambdas tau_final/60 trackNum_final times_final]; % tau here is converted from sec to min
            
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
clear plus1_* track*

% 18. save hard earned data
cd('/Users/jen/super-size/')
save('B5_ccSize.mat','compiled_data','compiled_mu','exptArray')


%% Part 3. sort data by nutrient condition

%  Adapted from figure1A.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('B5_ccSize.mat')
lamb = 10; % column in compiled_data that is lambda
sm = 2:9; % columns in compiled_data that are sizes
bt = 13; % colum in compiled_ata that is timestamp at birth


% 1. accumulate data from each condition
steady2fluc = 1; % row number in data structure
low = 2; 
ave = 3; 
high = 4;
condArray = [1,2,3,4];
timePast = 3; % allowed time past shift time in hours

sigmas = 3;
for ee = 1:length(condArray)
    
    condition = condArray(ee);
    
    
    % steady environment! concatenate data based on nutrient level
    if condition == low
        
        lambda_low = [];
        birthSizes_low = [];
        birthTimes_low = [];
        
        % loop through all experiments and store low data
        for expt = 1:length(compiled_data)
            
            expt_data = compiled_data{expt}{low,1};
            if ~isempty(expt_data)
                
                % isolate data
                expt_lambda = compiled_data{expt}{low,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                expt_sizes = compiled_data{expt}{low,1}(:,sm);
                expt_times = compiled_data{expt}{low,1}(:,bt);
                
                % concanetate individual cell cycle values
                lambda_low = [lambda_low; expt_lambda];
                birthSizes_low = [birthSizes_low; expt_sizes];
                birthTimes_low = [birthTimes_low; expt_times];
                clear expt_lambda expt_sizes expt_times
            end
            
        end
        clear expt expt_data
        
        condition_lambda = lambda_low;
        condition_sizes = birthSizes_low;
        condition_times = birthTimes_low;
        
        % isolate cycles within some error of replicate
        condition_mean = nanmean(condition_lambda);
        condition_std = nanstd(condition_lambda);
        
        lower = condition_lambda < condition_mean + (condition_std * sigmas);
        upper = condition_lambda > condition_mean - (condition_std * sigmas);
        combined = lower + upper;
        
        range_lambda = condition_lambda(combined == 2);
        range_sizes = condition_sizes(combined == 2,:);
        range_times = condition_times(combined == 2);
        clear condition_mean condition_std lower upper
        
        
        % store condition data
        sizes{low} = range_sizes;
        lambda{low} = range_lambda;
        birthTimes{low} = range_times;
        
        
        
    elseif condition == ave
        
        lambda_ave = [];
        birthSizes_ave = [];
        birthTimes_ave = [];
        
        % loop through all experiments and store ave data
        for expt = 1:length(compiled_data)
            
            expt_data = compiled_data{expt}{ave,1};
            if ~isempty(expt_data)
                
                % isolate data
                expt_lambda = compiled_data{expt}{ave,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                expt_sizes = compiled_data{expt}{ave,1}(:,sm);
                expt_times = compiled_data{expt}{ave,1}(:,bt);
                
                % concanetate individual cell cycle values
                lambda_ave = [lambda_ave; expt_lambda];
                birthSizes_ave = [birthSizes_ave; expt_sizes];
                birthTimes_ave = [birthTimes_ave; expt_times];
                clear expt_lambda expt_sizes expt_times
                
            end
        end
        clear expt expt_data
        
        
        condition_lambda = lambda_ave;
        condition_sizes = birthSizes_ave;
        condition_times = birthTimes_ave;
        
        % isolate cycles within 3 st dev of mean
        condition_mean = nanmean(condition_lambda);
        condition_std = nanstd(condition_lambda);
        
        lower = condition_lambda < condition_mean + (condition_std * sigmas);
        upper = condition_lambda > condition_mean - (condition_std * sigmas);
        combined = lower + upper;
        
        range_lambda = condition_lambda(combined == 2);
        range_sizes = condition_sizes(combined == 2,:);
        range_times = condition_times(combined == 2);
        clear condition_mean condition_std lower upper
        
        % store condition data
        lambda{ave} = range_lambda;
        sizes{ave} = range_sizes;
        birthTimes{ave} = range_times;
        
        
    elseif condition == high
        
        lambda_high = [];
        birthSizes_high = [];
        birthTimes_high = [];
        
        % loop through all experiments and store high data
        for expt = 1:length(compiled_data)
            
            expt_data = compiled_data{expt}{high,1};
            if ~isempty(expt_data)
                
                % isolate data
                expt_lambda = compiled_data{expt}{high,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                expt_sizes = compiled_data{expt}{high,1}(:,sm);
                expt_times = compiled_data{expt}{high,1}(:,bt);
                
                % concanetate individual cell cycle values
                lambda_high = [lambda_high; expt_lambda];
                birthSizes_high = [birthSizes_high; expt_sizes];
                birthTimes_high = [birthTimes_high; expt_times];
                clear expt_lambda expt_sizes expt_times
                
            end
        end
        clear expt expt_data
        
        condition_lambda = lambda_high;
        condition_sizes = birthSizes_high;
        condition_times = birthTimes_high;
        
        % isolate cycles within 3 st dev of mean
        condition_mean = nanmean(condition_lambda);
        condition_std = nanstd(condition_lambda);
        
        lower = condition_lambda < condition_mean + (condition_std * sigmas);
        upper = condition_lambda > condition_mean - (condition_std * sigmas);
        combined = lower + upper;
        
        range_lambda = condition_lambda(combined == 2);
        range_sizes = condition_sizes(combined == 2,:);
        range_times = condition_times(combined == 2);
        clear condition_mean condition_std lower upper
        
        % store condition data in cell corresponding to Condition Order
        lambda{high} = range_lambda;
        sizes{high} = range_sizes;
        birthTimes{high} = range_times;
        
        
    elseif condition == steady2fluc
        
        lambda_single = [];
        birthSizes_single = [];
        birthTimes_single = [];
        
        % loop through experiments and store timescale data
        for expt = 1:length(compiled_data)
            
            expt_data = compiled_data{expt}{steady2fluc,1};
            if ~isempty(expt_data)
                
                % isolate data
                expt_lambda = compiled_data{expt}{steady2fluc,1}(:,lamb); % note: mu is all instananeous vals in each cell cycle
                expt_sizes = compiled_data{expt}{steady2fluc,1}(:,sm);
                expt_times = compiled_data{expt}{steady2fluc,1}(:,bt);
                
                % specific to steady2fluc: trim times after 6th nutrient period
                shifttime = min(expt_times);
                maxtime = shifttime + timePast;
                new_lambda = expt_lambda(expt_times < maxtime,:);
                new_sizes = expt_sizes(expt_times < maxtime,:);
                new_times = expt_times(expt_times < maxtime,:);
                
                % concanetate individual cell cycle values
                lambda_single = [lambda_single; new_lambda];
                birthSizes_single = [birthSizes_single; new_sizes];
                birthTimes_single = [birthTimes_single; new_times];
                clear expt_lambda expt_sizes
            end
        end
        clear expt expt_data
        
        
        condition_lambda = lambda_single;
        condition_sizes = birthSizes_single;
        condition_times = birthTimes_single;
        
        % isolate cycles within 3 st dev of mean
        condition_mean = nanmean(condition_lambda);
        condition_std = nanstd(condition_lambda);
        
        lower = condition_lambda < condition_mean + (condition_std * sigmas);
        upper = condition_lambda > condition_mean - (condition_std * sigmas);
        combined = lower + upper;
        
        range_lambda = condition_lambda(combined == 2);
        range_sizes = condition_sizes(combined == 2,:);
        range_times = condition_times(combined == 2);
        clear condition_mean condition_std lower upper
        
        % store condition data in cell corresponding to Condition Order
        lambda{steady2fluc} = range_lambda;
        sizes{steady2fluc} = range_sizes;
        birthTimes{steady2fluc} = range_times;
        
    end
    
end
clear single_shift low ave high condition ee expt arrayIndex
clear condition_lambda condition_sizes condition_times
clear range_lambda range_sizes range_times combined
clear birthSizes_single birthSizes_low birthSizes_ave birthSizes_high
clear birthTimes_single birthTimes_low birthTimes_ave birthTimes_high
clear lambda_single lambda_low lambda_ave lambda_high


%% Part 4. bin data by time at birth and plot


% 0. initialize parameters for plotting
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
shape = 'o';
binsPerHour = 12; % 10 min bins
shiftFactor = 10;
timescale = 3600;

% time bins of interest in this analysis 
minBin = binsPerHour*shiftFactor; 

for metric = 2%1:length(sm)
    
    figure(metric)
    
    conditions = [1,2,3,4];
    counter_steady = 0;
    for jj = 1:length(conditions) % condition
        
        
        % 1. determine current condition of interest
        c = conditions(jj);
        
        
        % 2. isolate condition data
        cond_sizes = sizes{c};
        cond_times = birthTimes{c};
        cond_lambda = lambda{c};
        
        
        % 3. bin data into bins by period fraction
        %    note: because timestamps are relative to shift time, the data
        %    must first be shifted to avoided negative values, and then
        %    shifted back after binning.
        
        % time
        cond_times_shifted = cond_times + shiftFactor;
        bins_birth = ceil(cond_times_shifted * binsPerHour);
        
        % FRACTION
        %cond_fractions = cond_times - floor(cond_times);
        %cond_fractions_shifted = cond_fractions + shiftFactor;
        %bins_birth = ceil(cond_fractions_shifted * binsPerHour);
        
        % birth volume + 1 vs time of birth
        volBirth_plus = cond_sizes(:,metric);
        volBirth_binnedByBT = accumarray(bins_birth,volBirth_plus,[],@(x) {x});
        
        % lambda vs time of birth
        lambda_binnedByBT = accumarray(bins_birth,cond_lambda,[],@(x) {x});
        
        
        
        % 4. generate time vectors for plotting
        %    note: in time vector, zero = bin immediately before shift
        timeVector_BT = ((1:length(volBirth_binnedByBT))/binsPerHour) - shiftFactor;
        
        
        
        % 6. plot what we want for reals
        if c == 1
            
            % i. isolate data to time bins of interest
            lambda_final = lambda_binnedByBT(minBin:end);
            volBirth_final = volBirth_binnedByBT(minBin:end);
            timeVector_final = timeVector_BT(minBin:end);
              
            % ii. single shift in 30 min bins
            shift_Vb = cellfun(@mean,volBirth_final);
            shift_Vb_std = cellfun(@std,volBirth_final);
            shift_Vb_count = cellfun(@length,volBirth_final);
            shift_Vb_sem = shift_Vb_std./sqrt(shift_Vb_count);
            
            shift_gr = cellfun(@mean,lambda_final);
            shift_gr_std = cellfun(@std,lambda_final);
            shift_gr_count = cellfun(@length,lambda_final);
            shift_gr_sem = shift_gr_std./sqrt(shift_gr_count);
            
            % trim nans (steady to fluc only)
            nns = ~isnan(shift_gr);
            nns_gr = shift_gr(~isnan(shift_gr));
            nns_Vb = shift_Vb(~isnan(shift_Vb));
            
            cmap = parula(length(nns_gr)); % spread parula over length of 2d vector.
            for tt = 1:length(nns_gr)
                hold on
                plot(nns_gr(tt),log(nns_Vb(tt)),'Color',cmap(tt,:),'Marker',shape,'MarkerSize',10,'LineWidth',2)
                colorbar
            end
            clear shift_Vb shift_Vb_std shift_Vb_count shift_Vb_sem shift_gr shift_gr_std shift_gr_count shift_gr_sem
            clear tt
            
        else
            
            % i. define plotting color
            color = rgb(palette(c));
            counter_steady = counter_steady + 1;
            
            % ii. isolate data from bins
            lambda_trim1 = cond_lambda(bins_birth >= minBin);
            volume_birth_trim1 = volBirth_plus(bins_birth >= minBin);
            bins_birth_trim1 = bins_birth(bins_birth >= minBin);
            
            % iii. calculate mean and sem
            lambda_steady(counter_steady) = mean(lambda_trim1);
            volume_birth_steady(counter_steady) = mean(volume_birth_trim1);
            
            % iv. plot
            hold on
            plot(lambda_steady(counter_steady),log(volume_birth_steady(counter_steady)),'Color',color,'MarkerFaceColor',color,'Marker',shape,'MarkerSize',10,'LineWidth',2)
            colorbar
            
        end
        clear bins_birth bins_birth_trim1 bt lamb sigmas
        clear cond_lambda cond_sizes cond_times cond_times_shifted condArray
        clear lambda_binnedByBT lambda_final
        clear volBirth_*  volume_birth_trim1 lambda_trim1
        
    end
    
    % 7. plot fit line between steady environments
    
    fit = polyfit(lambda_steady,log(volume_birth_steady),1);
    x = linspace(lambda_steady(1),lambda_steady(end),10);
    y = fit(1)*x + fit(2);
    
    figure(metric)
    hold on
    plot(x,y,'Color',rgb('SlateGray'))
    %axis([0 3.8 0.5 2.2])
    xlabel('lambda')
    title(strcat('metric:',num2str(metric)))
end

