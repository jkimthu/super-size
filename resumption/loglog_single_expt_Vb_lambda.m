%% loglog plots for single-cells in one experiment


%  Goal: does our single-cell data display a power law relationship in loglog
%        plots of volume vs. growth rate ?


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. plot a loglog plot with single-cell data from all conditions in
%          for each experimental replicate


%  Last edit: Jen Nguyen, 2022 Jan 21
%  Commit: loglog plots of Vb_ and lambda_i for each experiment


%  OK let's go!


%% Part 0. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')

load('A1_div_ccSize.mat')
col_Vb = 1;
%col_Vdiv = 2;
%col_tau = 2;
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
exptArray = 1:13; % indeces 7-13 in compiled_data correspond to T = 15 and 60 min conditions


% 0. initialize data vectors to store stats for each experiment
compiled_Vb = cell(length(exptArray),1);
compiled_lam = cell(length(exptArray),1);


% 0. initialize plotting parameters
palette_steady = {'Indigo','GoldenRod','FireBrick','DarkSlateBlue','DarkGoldenRod','DarkRed'};
palette_fluc = {'DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue'};
shape = 'o';


%% Part 1. plot loglog plot of single-cell Vb_i vs. tau_i
%
%  - overlay data from all conditions, colored by condition
%  - make a new plot for each experiment

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
            else
                color = palette_fluc{4}; % color to plot T = 60 min data
            end
        else % steady    
            color = palette_steady{exp_cond-1}; % color to plot steady condition
        end
        
        
        % 3. isolate single-cell data from current condition
        condData = expData{exp_cond};
        if isempty(condData) == 1
            continue
        end
        
        Vb = condData.cc(:,col_Vb);
        lambda = condData.meta(:,col_lambda);
        
        
        
        % 4. remove 2 std away from mean Vb (outer most 5%, 2.5% on either side)
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        trim = 3;
        
        Vb_trim1 = Vb(Vb < Vb_mean + trim*sigma);
        lambda_trim1 = lambda(Vb < Vb_mean + trim*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        lambda_trim2 = lambda_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        
        clear tau_trim1 Vb_trim1 sigma Vb tau Vb_mean 

        
        % 5. plot scatter of size vs. interdivision time
        figure(ii)
        loglog(lambda_trim2,Vb_trim2,shape,'MarkerEdgeColor',rgb(color))
        hold on

        
        % 6. store n = # cells for legend
        nn(exp_cond) = length(lambda_trim2);
        clear Vb_trim2 lambda_trim2
        
    end
    
    figure(ii)
    ylabel('birth volume')
    xlabel('average growth rate (lambda)')
    title(strcat('expt index=',num2str(exp)))
    legend(num2str(nn))
    xlim([0.1 10])
    ylim([0 9])
    
end
clear exp_cond 
