%% ss15: single cell tau_i vs Vb_i across transitions


%  Goal: determine whether single cells display slope between tau_i and Vb_i.
%        plot single-cell tau_i and Vb_i over time under environmental transitions:
%         1. steady low to high
%         2. steady low to fluctuations


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. 
%  Part 2. 
%  Part 3. 
%          plot single-cell tau_i and Vb_i from trajectories by:
%            i) time of birth
%           ii) time of cell division


%  Last edit: Jen Nguyen, 2020 July 30
%  Commit: first commit, plot tau_i/Vb_i over single shifts


%  OK let's go!

%% pre-zero A. initialize analysis

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
exptArray = [21,22,37,38,39]; % 21-22 are low to high; 37-39 are low to 60 min fluc


% 0. initialize data vectors to store stats for each experiment
compiled_data = cell(length(exptArray),1);
compiled_mu = cell(length(exptArray),1);


%% pre-zero B. collect single cell data from shift experiments
%          birth size and instantaneous growth rates
%          copied from Part 2 of figure1A_division

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
    experimentFolder = strcat('/Users/jen/Documents/StockerLab/Source_data/');
    cd(experimentFolder)
    if index == 22
        filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
    else
        filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
    end
    load(filename,'D5','T');
    clear filename experimentFolder


    % for each condition in experiment
    for condition = 1:length(bubbletime)
            
            
        % 5. compile condition data matrix
        %    NOTE: compiling each condition separately restarts the curveFinder count at 1 per condition
        if index <= 22
            xy_start = min(xys(condition,:));
            xy_end = max(xys(condition,:));
        else
            xy_start = min(xys{condition,1});
            xy_end = max(xys{condition,:});
        end
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
save('ss15_div_ccSize.mat','compiled_data','exptArray')


%% Part 0. initialize analysis
%          with data organized from pre-zero parts A and B

clc
clear

cd('/Users/jen/super-size/')
load('ss15_div_ccSize.mat')
load('storedMetaData.mat')

% 0. initialize experimental conditions
expt_order = {'low','ave','high','upshift','steady2fluc'};

% 0. initialize colors for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};


%% Part 1. sort data by experimental condition, keeping replicates apart
%          adapted from ss16 Part 1, to work with dif experiment types vs fluc timescales


% 1. re-organize data by nutrient condition, keep replicates separated
upshift = 1:2; steady2fluc = 3:5;          % row number of cell in compiled_data  
shift = 1; low = 2; ave = 3; high = 4;     % row number in data structure in compiled_data{n,1}

counter = zeros(1,length(expt_order));
organized_data = [];
for exp = 1:length(compiled_data)
    
    for exp_cond = 1:4  % conditions per experiment = 4
        
        % assign to correct environment column based on condition
        if exp_cond == shift % shift
            
            if ismember(exp,upshift) == 1
                eoc = 4; % eoc = environmental order column
            elseif ismember(exp,steady2fluc) == 1
                eoc = 5; % eoc = environmental order column
            end
            
            counter(1,eoc) = counter(1,eoc) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,eoc),eoc} = data;
            
        else % steady
            
            if exp_cond == low % low
                eoc = 1;
            elseif exp_cond == ave
                eoc = 2;
            elseif exp_cond == high
                eoc = 3;
            end
            
            counter(1,eoc) = counter(1,eoc) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,eoc),eoc} = data;
            
        end
        
    end
    
end
clear eoc exp_cond data counter
clear upshift steady2fluc
clear compiled_data low ave high exp shift


%% Part 2. gather single cell data including slopes_i for each experimental replicate

%          calculate single-cell slopes (tau_i divided by Vb_i) from all

% 0. initialize step size for binning time
binsPerHour = 12; % 5 min bins


% 0. initialize variables for storing data
numexpt = length(expt_order);
slope_i = cell(numexpt,3);
Vb_i = cell(numexpt,3);
tau_i = cell(numexpt,3);
tbirth_i = cell(numexpt,3);
tdiv_i = cell(numexpt,3);


% 1. collect data from each steady condition
for sc = 1:3 % sc = steady condition
    
    % i. isolate data from current steady condition
    palette = palette_steady(1,sc);
    current_condData = organized_data(:,sc);
    
    % ii. determine replicates to loop through
    numreps = length(current_condData);
    for rr = 1:numreps
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            % iii. isolate replicate data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_tbirth = rep_meta(:,4);     % each row = cc; number = time at birth (h)
            rep_tdiv = rep_meta(:,5);       % each row = cc; number = time at division (h)
            
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % iv. calculate single-cell slopes
            rep_slopes = rep_taus./rep_birthvols;
            
            
            % v. remove items over 3 std from mean
            mm = mean(rep_slopes);
            thresh = 3 * std(rep_slopes);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = rep_slopes > upper;
            cut_down = rep_slopes < lower;
            toCut = cut_up + cut_down;
            
            slopes_final = rep_slopes(toCut == 0);
            taus_final = rep_taus(toCut == 0);
            Vb_final = rep_birthvols(toCut == 0);
            tbirth_final = rep_tbirth(toCut == 0);
            tdiv_final = rep_tdiv(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear rep_taus rep_birthvols rep_slopes rep_tdiv rep_tbirth
            
            
            % vi. plot distribution of single-cell slopes replicate
            figure(rr)
            histogram(slopes_final,100,'FaceColor',rgb(palette));
            xlabel('single-cell slope')
            ylabel('# cells')
            title(strcat('replicate',num2str(rr)))
            hold on
            
            
            % vii. store all means and standard deviations
            slope_i{rr,sc} = slopes_final;
            Vb_i{rr,sc} = Vb_final;
            tau_i{rr,sc} = taus_final;
            tbirth_i{rr,sc} = tbirth_final;
            tdiv_i{rr,sc} = tdiv_final;
            clear mf sf
            
        end
                
    end

end
clear sc rr slopes_final taus_final Vb_final tbirth_final col


% 2. collect data from each shifting condition
shf_counter = 0;
for shf = 1:2 % shf = shifting condition
    
    % i. isolate data from current steady condition
    if shf == 1
        col = 4;
    else
        col = 5;
    end
    
    palette = palette_fluc(1,shf);
    current_condData = organized_data(:,col);
    
    
    % ii. determine replicates to loop through
    numreps = length(current_condData);
    for rr = 1:numreps
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            
            shf_counter = shf_counter + 1;
            
            % iii. isolate replicate data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_tbirth = rep_meta(:,4);     % each row = cc; number = time at birth (h)
            rep_tdiv = rep_meta(:,5);       % each row = cc; number = time at division (h)
            
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % iv. calculate single-cell slopes
            rep_slopes = rep_taus./rep_birthvols;
            
            
            % v. remove items over 3 std from mean
            mm = mean(rep_slopes);
            thresh = 3 * std(rep_slopes);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = rep_slopes > upper;
            cut_down = rep_slopes < lower;
            toCut = cut_up + cut_down;
            
            slopes_final = rep_slopes(toCut == 0);
            taus_final = rep_taus(toCut == 0);
            Vb_final = rep_birthvols(toCut == 0);
            tbirth_final = rep_tbirth(toCut == 0);
            tdiv_final = rep_tdiv(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear rep_taus rep_birthvols rep_slopes rep_tdiv rep_tbirth
            
            
            % vi. plot distribution of single-cell slopes replicate
            figure(shf_counter)
            histogram(slopes_final,100,'FaceColor',rgb(palette));
            xlabel('single-cell slope')
            ylabel('# cells')
            title(strcat('replicate',num2str(rr)))
            hold on
            
            
            % vii. store all means and standard deviations
            slope_i{shf_counter,col} = slopes_final;
            Vb_i{shf_counter,col} = Vb_final;
            tau_i{shf_counter,col} = taus_final;
            tbirth_i{shf_counter,col} = tbirth_final;
            tdiv_i{shf_counter,col} = tdiv_final;
            clear mf sf
            
        end
                
    end

end
clear shf rr slopes_final taus_final Vb_final tbirth_final col
clear tdiv_final shf_counter numexpts numreps palette

%save('ss15_output.mat','slope_i','tau_i','Vb_i','tbirth_i','tdiv_i')


%% Part 3. plot slopes_i vs time for each experiment

%          go row by row (experiment)
%          plotting all columns (conditions) in same plot
%          highlight data post shift in filled circles,


cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss15/')


% 0. initialize colors for plotting shifted condition
pre_color = rgb('CornflowerBlue');
post_color = rgb('DarkTurquoise'); % different color for data after shift

% loop through experiments
for ee = 1:numexpt
    
    % 1. isolate current experiment data
    e_slopes = slope_i(ee,:);
    e_tbirths = tbirth_i(ee,:);
    e_tdivs = tdiv_i(ee,:);
    
    
    % 2. gather time of shift
    index = exptArray(ee);
    tshift = storedMetaData{index}.shiftTime/3600; % convert from seconds to hours
       
    
    % 3. loop through conditions to plot all single points
    for cond = 1:length(e_slopes)
        
        c_slopes = e_slopes{1,cond};
        if isempty(c_slopes) == 1
            continue
        else
            
            % i. isolate condition time data
            c_tbirths = e_tbirths{1,cond};
            c_tdivs = e_tdivs{1,cond};
            
            
            % continue depending on condition (steady vs shift)
            if cond < 4 % if STEADY
                
                % ii. determine plotting color by condition
                color = rgb(palette_steady{cond});
                
                % iii. plot slope_i vs tbirth_i
                figure(ee)
                scatter(c_tbirths,c_slopes,10,'MarkerEdgeColor',color)
                hold on
                
%                 % iv. plot slope_v vs tdiv_i
%                 figure(ee+10)
%                 scatter(c_tdivs,c_slopes,10,'MarkerEdgeColor',color)
%                 hold on
%                 
            else % if FLUCTUATING
                
                % ii. isolate data before and after shift
                preshift_birth = c_tbirths < tshift;
                
                pre_tbirths = c_tbirths(preshift_birth==1);
                pre_tdivs = c_tdivs(preshift_birth==1);
                pre_slopes = c_slopes(preshift_birth==1);
                
                post_tbirths = c_tbirths(preshift_birth==0);
                post_tdivs = c_tdivs(preshift_birth==0);
                post_slopes = c_slopes(preshift_birth==0);
                
                
                % iii. plot slope_i vs tbirth_i
                figure(ee)
                scatter(pre_tbirths,pre_slopes,10,'MarkerFaceColor',pre_color,'MarkerEdgeColor',pre_color)
                hold on
                scatter(post_tbirths,post_slopes,10,'MarkerFaceColor',post_color,'MarkerEdgeColor',post_color)
                ylabel('slope_i')
                xlabel('time of birth')
                
                figure(ee)
                saveas(gcf,strcat('ss15-index-',num2str(index)),'epsc')
                close(gcf)
                
%                 % iv. plot slope_v vs tdiv_i
%                 figure(ee+10)
%                 scatter(pre_tdivs,pre_slopes,10,'MarkerFaceColor',pre_color,'MarkerEdgeColor',pre_color)
%                 hold on
%                 scatter(post_tdivs,post_slopes,10,'MarkerFaceColor',post_color,'MarkerEdgeColor',post_color)

            end
            clear pre_tbirths pre_tdivs pre_slopes post_tbirths post_tdivs post_slopes
            
        end
    end
    clear e_slopes e_tbirths e_tdivs index tshift
       
end
clear ee 




