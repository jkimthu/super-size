%% ss1: distributions of birth size and interdivision time per replicate


%  Goal: do these distributions appear like those in our simulations?


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. calculate stats for each replicate


%  Last edit: Jen Nguyen, 2020 Feb 27
%  Commit: histograms of Vb and tau each replicate


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


%% Part 2. plot histogram calculate statistics for each replicate

% 0. initialize parameters for which calculate stats in organized_data
%steady = [1,6,7];
metric = 1; % volume at birth = col in in cc
lamb = 1;   % mean growth rate = col 1 in meta
tau = 2;    % interdivision time = col 2 in meta


% 0. initialize indexing of organized_data
experiments = size(organized_data);
total_idx = experiments(1)*experiments(2);


% 1. loop through indeces to calculate and store stats for each
rep_lambda_stats = cell(experiments);
rep_Vb_stats = cell(experiments);
rep_tau_stats = cell(experiments);
rep_counts = NaN(experiments);

counter = 0;
for ii = 1:total_idx
    
    i_data = organized_data{ii};
    if isempty(i_data) == 1
        continue
    end
    counter = counter + 1;
    if counter < 13
        color = rgb(palette{1});
        condition = 'low';
    elseif counter < 16
        color = rgb(palette{2});
        condition = 't30';
    elseif counter < 19
        color = rgb(palette{3});
        condition = 't300';
    elseif counter < 23
        color = rgb(palette{4});
        condition = 't900';
    elseif counter < 26
        color = rgb(palette{5});
        condition = 't3600';
    elseif counter < 39
        color = rgb(palette{6});
        condition = 'ave';
    else
        color = rgb(palette{7});
        condition = 'high';
    end
    
    % i. gather parameters
    gr = i_data.meta(:,lamb);
    Vb = i_data.cc(:,metric);
    dt = i_data.meta(:,tau);
    
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
    
    divtime = struct;
    divtime.mean = mean(dt);
    divtime.std = std(dt);
    divtime.sem = std(dt)/sqrt(cell_count);
    divtime.cv = std(dt)/mean(Vb) * 100;

    rep_lambda_stats{ii} = lambda;
    rep_Vb_stats{ii} = birthVolume;
    rep_tau_stats{ii} = divtime;
    rep_counts(ii) = cell_count;
    
    
    % iii. plot histogram of all parameters
    figure(1)
    histogram(gr,'FaceColor',color)
    xlabel('mean growth rate')
    xlim([0 5])
    title(condition)
    
    figure(2)
    histogram(Vb,'FaceColor',color)
    xlabel('birth volume')
    title(condition)
    
    figure(3)
    histogram(dt,'FaceColor',color)
    xlabel('interdivision time')
    xlim([0 80])
    title(condition)
    
    
    % vi. save and close for next loop
    figure(1)
    saveas(gcf,strcat('hist-gr-',num2str(counter)),'epsc')
    
    figure(2)
    saveas(gcf,strcat('hist-Vb-',num2str(counter)),'epsc')
    
    figure(3)
    saveas(gcf,strcat('hist-tau-',num2str(counter)),'epsc')
    
end

