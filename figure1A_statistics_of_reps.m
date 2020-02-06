%% Figure 1: size vs lambda (scatter)



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
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. calculate stats for each replicate


%  Last edit: Jen Nguyen, 2019 Feb 6
%  Commit: calculate mean, std, s.e.m. and cell count for each replicate

%  OK let's go!


%% Part 0. initialize data, saved from figure1A_division.m

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


%% Part 1. sort data by nutrient condition, keeping replicates apart
%          this code is the same as figure1A_scatter.m (2020 Feb 4)

% 1. re-organize data by nutrient condition, keep replicates separated
t0_5 = 1:3; t5 = 4:6; t15 = 7:10; t60 = 11:13; % row number in data structure 
fluc = 1; low = 2; ave = 3; high = 4;          % row number in data structure 

counter = zeros(1,length(environment_order));
organized_data = [];
for exp = 1:length(compiled_data)
    
    for exp_cond = 1:4
        
        % assign to correct environment column based on condition
        if exp_cond == 1 % fluc
            
            if ismember(exp,t0_5) == 1
                envr = 2;
            elseif ismember(exp,t5) == 1
                envr = 3;
            elseif ismember(exp,t15) == 1
                envr = 4;
            elseif ismember(exp,t60) == 1
                envr = 5;
            end
            
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
            
        else % steady
            
            if exp_cond == 2 % low
                envr = 1;
            elseif exp_cond == 3
                envr = 6;
            elseif exp_cond == 4
                envr = 7;
            end
            
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
            
        end
        
    end
    
end
clear envr exp_cond data counter


%% Part 2. calculate statistics for each replicate

% 0. initialize parameters for which calculate stats in organized_data
%steady = [1,6,7];
metric = 1; % 1 = volume at birth
lamb = 1;

% 0. initialize indexing of organized_data
experiments = size(organized_data);
total_idx = experiments(1)*experiments(2);


% 1. loop through indeces to calculate and store stats for each
rep_lambda_stats = cell(experiments);
rep_Vb_stats = cell(experiments);
rep_counts = NaN(experiments);

for ii = 1:total_idx
    
    i_data = organized_data{ii};
    if isempty(i_data) == 1
        continue
    end
    
    % i. gather parameters
    gr = i_data.meta(:,lamb);
    Vb = i_data.cc(:,metric);
    
    % ii. calculate and store stats
    cell_count = length(gr);

    lambda = struct;
    lambda.mean = mean(gr);
    lambda.std = std(gr);
    lambda.sem = std(gr)/sqrt(cell_count);
    lambda.cv = std(gr)/mean(gr) * 100;
    
    birthVolume = struct;
    birthVolume.mean = mean(Vb);
    birthVolume.std = std(Vb);
    birthVolume.sem = std(Vb)/sqrt(cell_count);
    birthVolume.cv = std(Vb)/mean(Vb) * 100;

    rep_lambda_stats{ii} = lambda;
    rep_Vb_stats{ii} = birthVolume;
    rep_counts(ii) = cell_count;
    
end

