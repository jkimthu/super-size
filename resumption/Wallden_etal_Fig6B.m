%% Wallden et al. Fig 6B with our data


%  Goal: does our single-cell data display volume vs. generation time distributions
%        observed in Wallden et al. 2016?


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. replicate Fig 6B for each individual experimental replicate


%  Last edit: Jen Nguyen, 2021 Mar 12
%  Commit: replicate Wallden et al Fig 6B for each experimental replicate


%  OK let's go!


%% Part 0. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')

load('A1_div_ccSize.mat')
col_Vb = 1;
col_Vdiv = 2;
col_tau = 2;
col_lambda = 1;
% the code generating this data structure is "figure1A_division"
% the "cc" data matrix has 8 columns: Vbirth Vdiv Lbirth Ldiv Wbirth Wdiv SA2Vbirth SA2Vdiv
% the "meta" data matrix has 5 columns:
%       1. lambda = mean growth rate aross cell cycle (1/h) 
%       2. tau = interdivision time (min)
%       3. track number (to distinguish cells from same or different lineages)
%       4. time of birth (h), where experiment start = 0 h
%       5. time of division (h), where experiment start = 0 h


% 0. create array of experiments to use in analysis
exptArray = 7:13; % these indeces in compiled_data correspond to T = 15 and 60 min conditions


% 0. initialize data vectors to store stats for each experiment
compiled_Vb = cell(length(exptArray),1);
compiled_Vdiv = cell(length(exptArray),1);
compiled_tau = cell(length(exptArray),1);
compiled_lambda = cell(length(exptArray),1);


% 0. initialize plotting parameters
palette_steady = {'Indigo','GoldenRod','FireBrick','DarkSlateBlue','DarkGoldenRod','DarkRed'};
%palette_fluc = {'DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue'};
palette_fluc = {'DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue'};
palette_fluc_div = {'DarkCyan','MediumBlue','RoyalBlue','SteelBlue'};
shape = 'o';


%% Part 1. replicate Fig 6B for each individual experimental replicate
%
%          for each condition in each experimental replicate, plot:
%           1. volume (birth and division) vs interdivision time
%           2. volume (birth and division) vs average growth rate of cell cycle

% for each experiment of interest
for ii = 1:length(exptArray) 
    
    exp = exptArray(ii);
    expData = compiled_data{exp};
    nn = zeros(4,1);
    
    % 1. loop through fluc and steady low, ave and high conditions
    for exp_cond = 1:4 
        
        
        % 2. determine color for plotting based on condition
        if exp_cond == 1 % fluc
            if exp < 11
                color = palette_fluc{3}; % color to plot T = 15 min data
                color_div = palette_fluc_div{3};
            else
                color = palette_fluc{4}; % color to plot T = 60 min data
                color_div = palette_fluc_div{4};
            end
        else % steady    
            color = palette_steady{exp_cond-1}; % color to plot steady condition
            color_div = palette_steady{exp_cond+2};
        end
        
        
        % 3. isolate single-cell data from current condition
        condData = expData{exp_cond};
        if isempty(condData) == 1
            continue
        end
        
        Vb = condData.cc(:,col_Vb);
        Vdiv = condData.cc(:,col_Vdiv);
        tau = condData.meta(:,col_tau);
        lambda = condData.meta(:,col_lambda);
        
        
        
        % 4. remove 2 std away from mean Vb (outer most 5%, 2.5% on either side)
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
        Vdiv_trim1 = Vdiv(Vb < Vb_mean + 2*sigma);
        tau_trim1 = tau(Vb < Vb_mean + 2*sigma);
        lambda_trim1 = lambda(Vb < Vb_mean + 2*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        Vdiv_trim2 = Vdiv_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        dt_trim2 = tau_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        lambda_trim2 = lambda_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        clear tau_trim1 Vb_trim1 Vdiv_trim1 sigma Vb Vdiv tau Vb_mean 
        clear lambda lambda_trim1

        
        % 5. plot scatter of size vs. interdivision time
        figure(ii)
        scatter(dt_trim2,Vb_trim2,'MarkerEdgeColor',rgb(color))
        hold on
        scatter(dt_trim2,Vdiv_trim2,'MarkerEdgeColor',rgb(color_div))
        
        
        figure(ii+10)
        scatter(lambda_trim2,Vb_trim2,'MarkerEdgeColor',rgb(color))
        hold on
        scatter(lambda_trim2,Vdiv_trim2,'MarkerEdgeColor',rgb(color_div))
        
        
        % 6. store n = # cells for legend
        nn(exp_cond) = length(dt_trim2);
        clear dt_trim2 Vb_trim2 Vdiv_trim2
        
    end
    
    figure(ii)
    ylabel('birth volume')
    xlabel('interdivision time (min)')
    legend(num2str(nn))
    
    figure(ii+10)
    ylabel('birth volume')
    xlabel('lambda (1/h)')
    legend(num2str(nn))
    axis([0 6 0 30])
    
end
clear exp_cond 




