%% measured vs. expected cell cycle parameters


%  Goal: how do various cell cycle parameters behave compared to the
%        expected values based on measurements of others?
%
%        Vdiv_i = Vb_i * 2^( lambda_i * tau_i)
%
%        i. what is lambda_i,measured vs. lambda_i,expected ?
%       ii. what is tau_i,measured vs. tau_i,expected ?
%      iii. what is Vdiv_i, measured vs. Vdiv_i, expected ?


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. plot a loglog plot with single-cell data from all conditions in
%          for each experimental replicate


%  Last edit: Jen Nguyen, 2022 Jan 21
%  Commit: comparable measured v expected values of cell cycle parameters


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
exptArray = 1:13; % indeces 7-13 in compiled_data correspond to T = 15 and 60 min conditions


% 0. initialize data vectors to store difference values for each experiment
compiled_dVdiv = cell(length(exptArray),1);
compiled_dlambda = cell(length(exptArray),1);
compiled_dtau = cell(length(exptArray),1);


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
        
        
        % 3. isolate single-cell data from current condition
        condData = expData{exp_cond};
        if isempty(condData) == 1
            continue
        end
        
        Vb = condData.cc(:,col_Vb);
        Vdiv = condData.cc(:,col_Vdiv);
        lambda = condData.meta(:,col_lambda);
        tau = condData.meta(:,col_tau);
        
        
        % 4. remove 2 std away from mean Vb (outer most 5%, 2.5% on either side)
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        trim = 2;
        
        Vb_trim1 = Vb(Vb < Vb_mean + trim*sigma);
        Vdiv_trim1 = Vdiv(Vb < Vb_mean + trim*sigma);
        lambda_trim1 = lambda(Vb < Vb_mean + trim*sigma);
        tau_trim1 = tau(Vb < Vb_mean + trim*sigma);
        
        Vb_im = Vb_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        Vdiv_im = Vdiv_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        lambda_im = lambda_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        tau_im = tau_trim1(Vb_trim1 > Vb_mean - trim*sigma);
        
        clear tau_trim1 Vb_trim1 Vdiv_trim1 lambda_trim1 sigma Vb Vdiv tau Vb_mean lambda

        
        % 5. calculate expected values
        %    from Vdiv_i = Vb_i * 2^( lambda_i * tau_i)
        
        % 5i. lambda_i, expected
        top = log((Vdiv_im./Vb_im));
        bottom = (tau_im/60 * log(2)); % convert tau from min to hours
        lambda_ie = top./bottom;
        clear bottom
        
        % 5ii. tau_i, expected
        bottom2 = (lambda_im * log(2));
        tau_ie = top./bottom2 * 60; % convert lambda from h to min
        clear top bottom2
        
        % 5iii. Vdiv_i, expected
        exponent = lambda_im .* (tau_im/60);
        Vdiv_ie = Vb_im .* 2.^exponent; % convert tau from min to hours
        clear exponent Vb_im

        
        % 6. calculate difference between measured and expected values
        %    - difference = measured - expected
        %    - store all values for easy plotting
        dlambda = lambda_im - lambda_ie;
        dVdiv = Vdiv_im - Vdiv_ie;
        dtau = tau_im - tau_ie;
        
        compiled_dVdiv{ii,exp_cond} = dVdiv;
        compiled_dlambda{ii,exp_cond} = dlambda;
        compiled_dtau{ii,exp_cond} = dtau;
        
        clear dlambda dVdiv dtau
        clear lambda_im Vdiv_im tau_im
        clear lambda_ie Vdiv_ie tau_ie
        
    end
    
end

% 7. quick check for difference between measured and expected values
cellfun(@mean,compiled_dVdiv)
cellfun(@mean,compiled_dlambda)
cellfun(@mean,compiled_dtau)

clear exp_cond
