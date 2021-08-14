%% ss25: using residuals to fit tau vs Vb data

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

clear fluc_lambda_i


%% Part 1. restructure steady tau_i, slope_i and Vb_i into columns by condition
%  adapted from Part 2 of ss24.m


% 1. initialize data vectors for slope_i and Vb_i
numcells_steady = [];
steady_tau_i = [];
steady_slope_i = [];
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
    
    slope_ii = compiled_steady{1,ii}.slope_i;
    steady_slope_i = [steady_slope_i; slope_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ii num_ii vb_ii tau_ii slope_ii


%% Part 2. restructure data for replicate based plotting

%  re-organize data so that fluc and steady data are in the same matrix,
%  such that:

%  column order denotes nutrient condition: fluc, low, ave, high
%  row order denotes replicate of fluc timescale: 1-4 (30 s), 5-8 (5 min), 9-12 (15 min), 13-16(60 min)

%  strategy:
%  1. restructure fluc data from 4x4 to 16x1
%  2. concatenate 16x3 steady data such that fluc data is 1st column

%  OK let's go!

%  1. restructure fluc data from 4x4 to 16x1

fluc_Vb_i_restruct = [];
fluc_tau_i_restruct = [];
fluc_slope_i_restruct = [];

for col = 1:4
    fluc_Vb_i_restruct = [fluc_Vb_i_restruct; fluc_Vb_i(:,col)];
    fluc_tau_i_restruct = [fluc_tau_i_restruct; fluc_tau_i(:,col)];
    fluc_slope_i_restruct = [fluc_slope_i_restruct; fluc_slopes_i(:,col)];
end


%  2. concatenate 16x3 steady data such that fluc data is 1st column
Vb_i_data = [fluc_Vb_i_restruct, steady_Vb_i];
tau_i_data = [fluc_tau_i_restruct, steady_tau_i];
slope_i_data = [fluc_slope_i_restruct, steady_slope_i];


clear col fluc_Vb_i fluc_tau_i fluc_slope_i steady_Vb_i steady_tau_i steady_slope_i
clear fluc_Vb_i_restruct fluc_tau_i_restruct fluc_slope_i_restruct


%% Part 3. calculate residuals R(c) across range of c

%  Strategy:
%
%  1. plot scatter for all single-cells of one replicate
%  2. 

% 0. initialize num of conditions
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
ss = size(tau_i_data);
rows = ss(1); % replicate experiments
cols = ss(2); % nutrient conditions (fluc, low, ave, high)
clear ss


% 1. loop through each experiment
for ee = 1:rows %
    
    ee_taus = tau_i_data(ee,:);
    ee_vbs = Vb_i_data(ee,:);
    ee_slopes = slope_i_data(ee,:);

    
    % 2. for each condition, plot single cell data and compile into aggregate vector
    aggregate_taus = [];
    aggregate_vbs = [];
    aggregate_slopes = [];
    
    for cc = 1:length(ee_vbs)
        
        % i. isolate condition tau_i, slope_i & Vb_i
        cc_taus = ee_taus{cc};
        if isempty(cc_taus) == 1
            continue
        end
        cc_vbs = ee_vbs{cc};
        cc_slopes = ee_slopes{cc};
        cc_color = rgb(palette{cc});
        
        
        % ii. plot tau_i vs. Vb_i and slope_i vs. Vb_i
        %figure(1)
        %scatter(cc_vbs,cc_taus,'MarkerFaceColor',cc_color,'MarkerEdgeColor',cc_color)
        %hold on
        
        figure(3)
        scatter(cc_vbs,cc_slopes,'MarkerFaceColor',cc_color,'MarkerEdgeColor',cc_color)
        hold on
        
        
        % iii. compile single-cell data into one vector per experiment
        aggregate_taus = [aggregate_taus; cc_taus];
        aggregate_vbs = [aggregate_vbs; cc_vbs];
        aggregate_slopes = [aggregate_slopes; cc_slopes];
        
  
    end
    clear cc_slopes cc_taus cc_vbs cc_color cc
    
    figure(1)
    xlabel('Vb_i (cubic um)')
    ylabel('tau_i (min)')
    
    figure(2)
    xlabel('Vb_i (cubic um)')
    ylabel('slope_i (min/cubic um)')
    
    
    % 3. fit hyperbola to aggregate single-cell data
    %    hyperbola = @(b,x) b(1) + b(2)./(x);    
    %    residual norm cost func = @(b) norm(y - hyprb(b,x)); 
    %    initial parameters, B0 = [1,1]
    %    estimated parameters, B = [b(1), b(2)];
    y_tau = aggregate_taus;
    y_slope = aggregate_slopes;
    x_vb = aggregate_vbs;
    B0 = [1; 1];
    
    % tau = y parameter
    hyprb = @(b,x_vb) b(1) + b(2)./(x_vb);              
    rncf = @(b) norm(y_tau - hyprb(b,x_vb));                       
    B_tau = fminsearch(rncf, B0);                               
    
    % slope = y parameter
    hyprb = @(c,x_vb) c(1) + c(2)./(x_vb);
    rncf = @(c) norm(y_slope - hyprb(c,x_vb));   
    B_slope = fminsearch(rncf, B0);                 
    
    
    % 4. plot hyperbola over scattered data
    x = linspace(1,9,20);
    y = B_tau(1) + B_tau(2)./x;
    figure(1)
    plot(x,y,'r')
    
    y = B_slope(1) + B_slope(2)./x;
    figure(2)
    plot(x,y,'r') 
    
    
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
     
