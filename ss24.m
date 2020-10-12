%% ss24: bootstrapping to determine significance of fit of slope


%  Goal: determine if slope is a statistically significant parameter with
%        bootstrap hypothesis testing


%  Strategy:
%
%       For a given experimental condition:
%        1.  Sample observed data with replacement to build many sample
%            populations
%        2.  Measure slope of tau_i vs Vb_i from each sample population
%        3.  Determine probability of getting something more extreme than
%            observed slope (p-value)


%  Last edit: jen, 2020 Oct 12
%  Commit: first commit, bootstrap determines significance of slope


%  OK let's go!

%% Part 0. initialize cell cycle data, saved from ss17.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('ss17.mat')

clear fluc_lambda_i fluc_slopes_i


%% Part 1. restructure steady slope_i and Vb_i into columns by condition
%  adapted from Part 2 of ss23.m


% 1. initialize data vectors for slope_i and Vb_i
numcells_steady = [];
steady_tau_i = [];
steady_Vb_i = [];


% 2. loop through timescale groups to concatenate data into columns by
%    nutrient conditiion
for ii = 1:length(compiled_steady)
    
    num_ii = cellfun(@length,compiled_steady{1,ii}.tau_i);
    numcells_steady = [numcells_steady; num_ii];
    
    tau_ii = compiled_steady{1,ii}.tau_i;
    steady_tau_i = [steady_tau_i; tau_ii];
    
    vb_ii = compiled_steady{1,ii}.Vb_i;
    steady_Vb_i = [steady_Vb_i; vb_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ii num_ii vb_ii tau_ii


%% Part 2. restructure data for replicate based plotting

%  re-organize data so that fluc and steady data are in the same matrix,
%  such that:

%  column order denotes nutrient condition: fluc, low, ave, high
%  row order denotes replicate of fluc timescale: 1-4 (30 s), 5-8 (5 min), 9-12 (15 min), 13-16(60 min)

%  strategy:

%     (perform both steps for tau and Vb data)
%  1. restructure fluc data from 4x4 to 16x1
%  2. concatenate 16x3 steady data such that fluc data is 1st column

%  OK let's go!

%  1. restructure fluc data from 4x4 to 16x1

fluc_Vb_i_restruct = [];
fluc_tau_i_restruct = [];

for col = 1:4
    fluc_Vb_i_restruct = [fluc_Vb_i_restruct; fluc_Vb_i(:,col)];
    fluc_tau_i_restruct = [fluc_tau_i_restruct; fluc_tau_i(:,col)];
end


%  2. concatenate 16x3 steady data such that fluc data is 1st column
Vb_i_data = [fluc_Vb_i_restruct, steady_Vb_i];
tau_i_data = [fluc_tau_i_restruct, steady_tau_i];


clear col fluc_Vb_i fluc_tau_i steady_Vb_i steady_tau_i
clear fluc_Vb_i_restruct fluc_tau_i_restruct


%% Part 3. bootstrap from each replicate condition to measure p-value

% 0. initialize num of conditions
ss = size(tau_i_data);
rows = ss(1); % replicate experiments
cols = ss(2); % nutrient conditions (fluc, low, ave, high)
clear ss

% 0. initialize number of sampled populations, k, to generate
k = 1000;


% 1. loop through each experiment
sampled_slopes = cell(rows,cols);
sampled_yints = cell(rows,cols);
measured_slopes = nan(rows,cols);
measured_yints = nan(rows,cols);
for ee = 1:rows %
    
    ee_taus = tau_i_data(ee,:);
    ee_vbs = Vb_i_data(ee,:);

    
    % 2. for each condition, sample with replacement to create 1000 populations
    for cc = 1:length(ee_taus)
        
        % i. sample tau_i & Vb_i pairs, not individually
        cc_taus = ee_taus{cc};
        if isempty(cc_taus) == 1
            continue
        end
        
        cc_vbs = ee_vbs{cc};
        pairs_measured = [cc_taus,cc_vbs];
        
        
        N = length(pairs_measured);
        cc_slopes = zeros(k,1);
        cc_yints = zeros(k,1);
        for samples = 1:k
                        
            % ii. generate k sampled populations of equal pop size, N
            pairs_sampled = datasample(pairs_measured,N); % datasample always samples the same pairs!


            % iii. measure slope
            fit_sample = polyfit(pairs_sampled(:,2),pairs_sampled(:,1),1);
  
            
            % iv. store data
            cc_slopes(samples,1) = fit_sample(1);
            cc_yints(samples,1) = fit_sample(2);
            
        end
        clear fit_sample samples
         
        sampled_slopes{ee,cc} = cc_slopes;
        sampled_yints{ee,cc} = cc_yints;
        
        fit_measured = polyfit(cc_vbs,cc_taus,1);
        measured_slopes(ee,cc) = fit_measured(1);
        measured_yints(ee,cc) = fit_measured(2);
  
    end
    clear ee_taus ee_vbs cc_taus cc_vbs pairs_measured cc fit_measured N pairs_sampled
    
end
clear ee 


%% Part 4. calculate p-value for each condition
%
%  pval = sum(s >= s0)/N, where s0 = observed slope - mean sampled slope
%                                 s  = each sampled slope farther or equal
%                                      in distance to mean of distribution

% 1. calculate mean for sampled slopes
sample_slope_mean = cellfun(@mean,sampled_slopes);


% 2. calculate s0 = observed slope - mean sampled slope
s0 = abs(measured_slopes-sample_slope_mean);


% 3. measure distance of all sampled slopes
pvals = nan(rows,cols);
for cond = 1:rows*cols
    
    currCond = sampled_slopes{cond};
    currMean = sample_slope_mean(cond);
    if isempty(currCond) == 1
        continue
    end
    currDist = abs(currCond-currMean);
    
    
    % 4. determine fraction of k that are farther from the sampled mean than s0
    currS0 = s0(cond);
    far_slopes = sum(currDist >= currS0);
    pvals(cond) = far_slopes/k;
    
end
clear cond currS0 far_slopes currDist currMean currCond 

isSig = pvals > 0.95;
     
