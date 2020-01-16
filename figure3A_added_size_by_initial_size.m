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
%  Commit: update to include growth rate and tau data


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


sigmas = 3;
for ee = 1:length(environment_order)
    
    condition = environment_order{ee};
    
    if ischar(condition) == 1
        
        % steady environment! concatenate data based on nutrient level
        if strcmp(condition,'low') == 1
            
            meta_low = [];
            ccSizes_low = [];
            plusSizes_low = [];

            % loop through all experiments and store low data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{low,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_low = [meta_low; expt_meta];
                    ccSizes_low = [ccSizes_low; expt_sizes];
                    plusSizes_low = [plusSizes_low; expt_plus]; 
                    clear expt_lambda expt_sizes expt_plus
                end
                
            end 
            clear expt expt_data expt_plus
            
            sizes{1} = ccSizes_low;
            metas{1} = meta_low;
            plus{1} = plusSizes_low;
            
        elseif strcmp(condition,'ave') == 1
            
            meta_ave = [];
            ccSizes_ave = [];
            plusSizes_ave = [];
            
            % loop through all experiments and store ave data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{ave,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_ave = [meta_ave; expt_meta];
                    ccSizes_ave = [ccSizes_ave; expt_sizes];
                    plusSizes_ave = [plusSizes_ave; expt_plus];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data expt_plus
 
            
            % store condition data
            metas{6} = meta_ave; %range_lambda;
            sizes{6} = ccSizes_ave; %range_sizes;
            plus{6} = plusSizes_ave;

            
        elseif strcmp(condition,'high') == 1
            
            meta_high = [];
            ccSizes_high = [];
            plusSizes_high = [];
            
            % loop through all experiments and store high data
            for expt = 1:length(compiled_data)
                
                expt_data = compiled_data{expt,1}{high,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_high = [meta_high; expt_meta];
                    ccSizes_high = [ccSizes_high; expt_sizes];
                    plusSizes_high = [plusSizes_high; expt_plus];
                    clear expt_lambda expt_sizes
                    
                end
            end
            clear expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            metas{7} = meta_high; 
            sizes{7} = ccSizes_high; 
            plus{7} = plusSizes_high;
            
        end
    else
        
        % fluctuating environment! concatenate based on timescale
        if condition == 30
            idx = [2,3,4]; % ID of experiments with this fluc timescale
            
            meta_30 = [];
            ccSizes_30 = [];
            plusSizes_30 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_30 = [meta_30; expt_meta];
                    ccSizes_30 = [ccSizes_30; expt_sizes];
                    plusSizes_30 = [plusSizes_30; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            metas{2} = meta_30; 
            sizes{2} = ccSizes_30; 
            plus{2} = plusSizes_30;
            
            
        elseif condition == 300
            idx = [5,6,7]; % ID of experimennts with this fluc timescale
            
            meta_300 = [];
            ccSizes_300 = [];
            plusSizes_300 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_300 = [meta_300; expt_meta];
                    ccSizes_300 = [ccSizes_300; expt_sizes];
                    plusSizes_300 = [plusSizes_300; expt_plus];
                    
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            metas{3} = meta_300; %range_lambda;
            sizes{3} = ccSizes_300; %range_sizes;
            plus{3} = plusSizes_300;
            
        elseif condition == 900
            idx = [9,10,11,12]; % ID of experimennts with this fluc timescale
            
            meta_900 = [];
            ccSizes_900 = [];
            plusSizes_900 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_900 = [meta_900; expt_meta];
                    ccSizes_900 = [ccSizes_900; expt_sizes];
                    plusSizes_900 = [plusSizes_900; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            metas{4} = meta_900; 
            sizes{4} = ccSizes_900; 
            plus{4} = plusSizes_900;
            
            
        elseif condition == 3600
            idx = [13,14,15]; % ID of experimennts with this fluc timescale
            
            meta_3600 = [];
            ccSizes_3600 = [];
            plusSizes_3600 = [];
            
            % loop through experiments and store timescale data
            for arrayIndex = 1:length(idx)
                
                expt = find(exptArray == idx(arrayIndex));
                expt_data = compiled_data{expt,1}{fluc,1};
                if ~isempty(expt_data)
                    
                    % isolate data
                    expt_meta = expt_data.meta; % note: mu is all instananeous vals in each cell cycle
                    expt_sizes = expt_data.cc;
                    expt_plus = expt_data.plus1;
                    
                    % concanetate individual cell cycle values
                    meta_3600 = [meta_3600; expt_meta];
                    ccSizes_3600 = [ccSizes_3600; expt_sizes];
                    plusSizes_3600 = [plusSizes_3600; expt_plus];
                    clear expt_lambda expt_sizes
                end
            end
            clear arrayIndex expt expt_data expt_plus
            
            
            % store condition data in cell corresponding to Condition Order
            metas{5} = meta_3600; 
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
lambda = 1; tau = 2;

% 0. initialize colors for plotting
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};


% 1. loop through conditions
for condition = 1:length(environment_order)
    
    % 2. isolate growth rate, tau, added and birth size of interest
    addedV = size_added{1,condition}(:,vol);
    birthV = size_birth{1,condition}(:,vol);
    gr = metas{1,condition}(:,lambda);
    interdiv = metas{1,condition}(:,tau);
    
    % 3. assign each birth size to a bin
    birth_bin = ceil(birthV/binsize);
    
    % 4. sort added size by birth bin
    adder_binned = accumarray(birth_bin,addedV,[],@(x) {x});
    adder_mean = cellfun(@mean, adder_binned);
    adder_std = cellfun(@std, adder_binned);
    adder_count = cellfun(@length, adder_binned);
    adder_sem = adder_std./sqrt(adder_count);
    
    % 5. sort growth rate by birth bin
    gr_binned = accumarray(birth_bin,gr,[],@(x) {x});
    gr_mean = cellfun(@mean, gr_binned);
    gr_std = cellfun(@std, gr_binned);
    gr_count = cellfun(@length, gr_binned);
    gr_sem = gr_std./sqrt(gr_count);
    
    % 6. sort interdivision time by birth bin
    tau_binned = accumarray(birth_bin,interdiv,[],@(x) {x});
    tau_mean = cellfun(@mean,tau_binned);
    tau_std = cellfun(@std,tau_binned);
    tau_count = cellfun(@length,tau_binned);
    tau_sem = tau_std./sqrt(tau_count);
    
    
    figure(condition)
    errorbar(gr_mean,gr_sem,'.','Color',rgb('Silver'))
    hold on
    errorbar(tau_mean/10,tau_sem/10,'.','Color',rgb('SlateGray'))
    hold on
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






