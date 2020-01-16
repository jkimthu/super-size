%% Figure 3A: plot added size as a function of initial size



%  Goal: plot individual-level added size as a function of initial size
%
%        Do cells target an added size? (constant across initial size)
%        Or is added size an artifact? (initial size and growth rate?)



%  Strategy: 
%
%  Part 1. initialize analysis
%  Part 2. collect single cell birth size and instantaneous growth rates
%  Part 3. plotting
%  Part 4. fit best line


%  Last edit: Jen Nguyen, 2019 Jan 16
%  Commit: plot added size binned by birth size

%  OK let's go!

%% Part 1. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')
dataIndex = find(~cellfun(@isempty,storedMetaData));

           
%% Part 2. sort data by nutrient condition


% 0. initialize plotting parameters
environment_order = {'low',30,300,900,3600,'ave','high'};


% 1. accumulate data from each condition
fluc = 1; % row number in data structure
low = 2; 
ave = 3; 
high = 4;
lamb = 1;


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


%% Part 3. calculate added size, store birth size in similar variable

for cond = 1:length(environment_order)
    
    Vbirth = sizes{cond}(:,1);
    Vdiv = sizes{cond}(:,2);
    added(:,1) = Vdiv - Vbirth;
    birth(:,1) = Vbirth;
    
    Lbirth = sizes{cond}(:,3);
    Ldiv = sizes{cond}(:,4);
    added(:,2) = Ldiv - Lbirth;
    birth(:,2) = Lbirth;
    
    Wbirth = sizes{cond}(:,5);
    Wdiv = sizes{cond}(:,6);
    added(:,3) = Wdiv - Wbirth;
    birth(:,3) = Wbirth;
    
    SA2Vbirth = sizes{cond}(:,7);
    SA2Vdiv = sizes{cond}(:,8);
    added(:,4) = SA2Vdiv - SA2Vbirth;
    birth(:,4) = SA2Vbirth;
    
    size_added{cond} = added;
    size_birth{cond} = birth;
    
    clear Vdiv Ldiv Wdiv SA2Vdiv Vbirth Lbirth Wbirth SA2Vbirth
    clear added birth
    
end
clear cond lamb

%% Part 4. bin added size by birth size

% 0. initialize birth size bin and size columns
binsize = 0.2; % cubic microns
vol = 1; len = 2; wid = 3; sa2v = 4; % vol is currently hard coded below!

% 0. initialize colors for plotting
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};


% 1. loop through conditions
for condition = 1:length(environment_order)
    
    % 2. isolate added and birth size of interest
    addedV = size_added{1,condition}(:,len);
    birthV = size_birth{1,condition}(:,len);
    
    % 3. assign each birth size to a bin
    birth_bin = ceil(birthV/binsize);
    
    % 4. sort added size by birth bin
    adder_binned = accumarray(birth_bin,addedV,[],@(x) {x});
    adder_mean = cellfun(@mean, adder_binned);
    adder_std = cellfun(@std, adder_binned);
    adder_count = cellfun(@length, adder_binned);
    adder_sem = adder_std./sqrt(adder_count);
    
    figure(condition)
    errorbar(adder_mean,adder_sem,'o','Color',rgb(palette{condition}))
    xlabel('birth vol x 0.2')
    ylabel('added vol')
    if isstring(environment_order)
        title(environment_order{condition})
    else
        title(num2str(environment_order{condition}))
    end
    axis([3 53 0 7])

end






