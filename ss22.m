%% ss22: single-cell vs replicate average trend of tau vs Vb


%  Goal: compare individual and mean data in tau and Vb

%  Key question: How do single-cell data obey scaling?
%  Plot tau_i vs Vb_i with all data at the single-cell level, with: 
%           i. single-cell steady-state data in gray
%          ii. single-cell fluctuating data in color (blue)
%         iii. average of each condition replicate in color.


%  Strategy: 
%
%  Part 0. initialize cell cycle data, saved from ss17.m
%  Part 1. calculate replicate means for plotting in Part 2
%  Part 2. plot replicate means only
%  Part 3. re-structure data for replicate based plotting
%  Part 4. plot single-cell scatter overlaid by replicate means


%  Last edit: Jen Nguyen, 2020 Oct 11
%  Commit: first commit, plots spread of individual tau vs Vb with averaged
%          data to compare trends


%  OK let's go!

%% Part 0. initialize cell cycle data, saved from ss17.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('ss17.mat')


%% Part 1. calculate replicate means for plotting in part 2


% 1. calculate mean, standard deviation and coefficient of variation 
%    of each replicates from FLUCTUATING conditions
% 
% fluc_slope.mean = cellfun(@nanmean,fluc_slopes_i);
% fluc_slope.std = cellfun(@nanstd,fluc_slopes_i);
% fluc_slope.cv = (fluc_slope.std./fluc_slope.mean) * 100;

fluc_Vb.mean = cellfun(@nanmean,fluc_Vb_i);
fluc_Vb.std = cellfun(@nanstd,fluc_Vb_i);
fluc_Vb.cv = (fluc_Vb.std./fluc_Vb.mean) * 100;

fluc_tau.mean = cellfun(@nanmean,fluc_tau_i);
fluc_tau.std = cellfun(@nanstd,fluc_tau_i);
fluc_tau.cv = (fluc_tau.std./fluc_tau.mean) * 100;

% fluc_lambda.mean = cellfun(@nanmean,fluc_lambda_i);
% fluc_lambda.std = cellfun(@nanstd,fluc_lambda_i);
% fluc_lambda.cv = (fluc_lambda.std./fluc_lambda.mean) * 100;


% 2. calculate mean, standard deviation and coefficient of variation 
%    of each replicates from STEADY conditions

%    - separated by parallel fluctuation timescale

%    i. birth volume
xx = cell(1,4);
yy = cell(1,4);
zz = cell(1,4);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    x = cellfun(@nanmean,compiled_steady{1,tscale}.Vb_i);
    y = cellfun(@nanstd,compiled_steady{1,tscale}.Vb_i);
    xx{tscale} = x;
    yy{tscale} = y;
    zz{tscale} = (y./x) * 100;
end
steady_Vb.mean = xx;
steady_Vb.std = yy;
steady_Vb.cv = zz;
clear x y xx yy zz tscale


%   ii. division time (tau)
aa = cell(1,4);
bb = cell(1,4);
cc = cell(1,4);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    a = cellfun(@nanmean,compiled_steady{1,tscale}.tau_i);
    b = cellfun(@nanstd,compiled_steady{1,tscale}.tau_i);
    aa{tscale} = a;
    bb{tscale} = b;
    cc{tscale} = (b./a) * 100;
end
steady_tau.mean = aa;
steady_tau.std = bb;
steady_tau.cv = cc;
clear a b aa bb cc tscale


%  iii. growth rate (lambda)
% qq = cell(1,4);
% rr = cell(1,4);
% ss = cell(1,4);
% for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
%     q = cellfun(@nanmean,compiled_steady{1,tscale}.lambda_i);
%     r = cellfun(@nanstd,compiled_steady{1,tscale}.lambda_i);
%     qq{tscale} = q;
%     rr{tscale} = r;
%     ss{tscale} = (r./q) * 100;
% end
% steady_lambda.mean = qq;
% steady_lambda.std = rr;
% steady_lambda.cv = ss;
% clear q r qq rr ss tscale


% iv. slopes (tau_i/Vb_i)
% dd = cell(1,4);
% ee = cell(1,4);
% ff = cell(1,4);
% for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
%     d = cellfun(@nanmean,compiled_steady{1,tscale}.slope_i);
%     e = cellfun(@nanstd,compiled_steady{1,tscale}.slope_i);
%     dd{tscale} = d;
%     ee{tscale} = e;
%     ff{tscale} = (e./d) * 100;
% end
% steady_slope.mean = dd;
% steady_slope.std = ee;
% steady_slope.cv = ff;
% clear d e dd ee ff tscale



% 3. calculate mean, standard deviation and coefficient of variation 
%    from compilation of STEADY conditions

%    i. aggregate all steady data into columns by condition
gg = [];
hh = [];
%ii = [];
%jj = [];
for tscale = 1:4
    gg = [gg; compiled_steady{1,tscale}.Vb_i];     % Vb
    hh = [hh; compiled_steady{1,tscale}.tau_i];    % tau
    %ii = [ii; compiled_steady{1,tscale}.lambda_i]; % lambda
    %jj = [jj; compiled_steady{1,tscale}.slope_i];
end
clear tscale

%   ii. aggregate all data belonging to same condition
for condition = 1:3
    
    all_vb = [];
    all_tau = [];
    %all_lam = [];
    %all_slo = [];
    for repl = 1:length(gg)
        
        cr_vb = gg{repl,condition};
        if isempty(cr_vb) == 1
            continue
        else
            cr_tau = hh{repl,condition};
            %cr_lam = ii{repl,condition};
            %cr_slo = jj{repl,condition};
            
            all_vb = [all_vb; cr_vb];
            all_tau = [all_tau; cr_tau];
            %all_lam = [all_lam; cr_lam];
            %all_slo = [all_slo; cr_slo];
        end
        clear cr_vb cr_tau cr_lam cr_slo
        
    end
    aggregate_steady_Vb{1,condition} = all_vb;
    aggregate_steady_tau{1,condition} = all_tau;
    %aggregate_steady_lambda{1,condition} = all_lam;
    %aggregate_steady_slope{1,condition} = all_slo;
    
end
clear all_vb all_tau all_lam all_slo
clear cr_vb cr_tau cr_lam cr_slo
clear gg hh ii jj repl condition


%% Part 2. plot replicate means only

%          mean tau_i vs. mean Vb_i of each replicate


% 0. initialize colors and shape for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};
shape = 'o';


% 1. first, plot mean slope_i vs. mean Vb_i from all steady replicates

%   i. gather mean and std of slope_i and Vb_i into columns by condition
mean_tau_steady = [];
mean_Vb_steady = [];
std_tau_steady = [];
std_Vb_steady = [];
numcells_steady = [];
steady_tau_i = [];
steady_Vb_i = [];

for ii = 1:length(steady_tau.mean)
    
    % ii. loop through timescale groups to concatenate steady means
    ts_tau = steady_tau.mean{1,ii};
    mean_tau_steady = [mean_tau_steady; ts_tau];
    
    ts_Vb = steady_Vb.mean{1,ii};
    mean_Vb_steady = [mean_Vb_steady; ts_Vb];
    
    s_std = steady_tau.std{1,ii};
    std_tau_steady = [std_tau_steady; s_std];
    
    v_std = steady_Vb.std{1,ii};
    std_Vb_steady = [std_Vb_steady; v_std];
    
    num_ii = cellfun(@length,compiled_steady{1,ii}.tau_i);
    numcells_steady = [numcells_steady; num_ii];
    
    tau_ii = compiled_steady{1,ii}.tau_i;
    steady_tau_i = [steady_tau_i; tau_ii];
    
    vb_ii = compiled_steady{1,ii}.Vb_i;
    steady_Vb_i = [steady_Vb_i; vb_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ts_Vb ts_tau s_std v_std ii num_ii vb_ii tau_ii


% iii. calculate standard error of the mean for each steady replicate
%      s.e.m. = st dev / sqrt(n)
sem_tau_steady = std_tau_steady./sqrt(numcells_steady);
sem_Vb_steady =  std_Vb_steady./sqrt(numcells_steady);


% iv. plot data from STEADY conditions
%     color hollow points by condition

for sc = 1:3
    
    % i. define color by condition
    sc_color = rgb(palette_steady{sc});
    
    % ii. isolate data by condition column
    yvals = mean_tau_steady(:,sc);
    xvals = mean_Vb_steady(:,sc);
    
    yerr = sem_tau_steady(:,sc);
    xerr = sem_Vb_steady(:,sc);

    % iii. plot condition data
    figure(20) % slope vs Vb (mean +/- sem)
    scatter(xvals,yvals,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(xvals,yvals,xerr,'.','horizontal','Color',sc_color)
    hold on
    errorbar(xvals,yvals,yerr,'.','Color',sc_color)
    hold on
    
end
clear sc sc_color
ylabel('tau_i')
xlabel('Vb_i')



% 2. second, plot mean slope_i vs. mean Vb_i from all fluctuating replicates

%   i. gather mean and std of slope_i and Vb_i into columns by timescale
mean_tau_fluc = fluc_tau.mean;
mean_Vb_fluc = fluc_Vb.mean;
std_tau_fluc = fluc_tau.std;
std_Vb_fluc = fluc_Vb.std;

% ii. calculate standard error of the mean for each steady replicate
%     s.e.m. = st dev / sqrt(n)
sem_tau_fluc = std_tau_fluc./sqrt(numcells_fluc);
sem_Vb_fluc =  std_Vb_fluc./sqrt(numcells_fluc);


% iii. plot data from FLUCTUATING conditions
%      color hollow points by timescale

for ts = 1:4
    
    % i. define color by condition
    color = rgb(palette_fluc{ts});
    
    % ii. isolate data by condition column
    yvals = mean_tau_fluc(:,ts);
    xvals = mean_Vb_fluc(:,ts);
    
    yerr = sem_tau_fluc(:,ts);
    xerr = sem_Vb_fluc(:,ts);

    % iii. plot condition data
    figure(20) % slope vs Vb (mean +/- sem)
    scatter(xvals,yvals,100,'filled','MarkerFaceColor',color)
    hold on
    errorbar(xvals,yvals,xerr,'.','horizontal','Color',color)
    hold on
    errorbar(xvals,yvals,yerr,'.','Color',color)
    hold on
    
end
clear color ts yvals xvals yerr xerr
title('mean tau_i vs mean Vb_i of each replicate')

clear mean_tau_steady mean_Vb_steady std_tau_steady std_Vb_steady 
clear numcells_steady sem_tau_steady sem_Vb_steady
clear mean_tau_fluc mean_Vb_fluc std_tau_fluc std_Vb_fluc sem_tau_fluc sem_Vb_fluc


%% Part 3. re-stucture data for replicate based plotting

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
clear fluc_Vb_i_restruct fluc_tau_i_restruct fluc_lambda_i fluc_slopes_i
clear steady_Vb steady_tau fluc_tau fluc_Vb


%% Part 4. plot single-cell scatter overlaid by replicate means

%         in separate figures, overlay single rep condition scatter over mean data
%         filling in 'o' for replicate mean of scatter
%         other replicate means can have hollow 'o'

%  OK let's go!


% 0. initialize data
cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss22/')
counter = 0;
%palette_all = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
palette_i = {'DodgerBlue','Silver','Silver','Silver'};
palette_rep = {'RoyalBlue','Indigo','GoldenRod','FireBrick'; % hard-coded color for mean points!
               'RoyalBlue','Indigo','GoldenRod','FireBrick';
               'RoyalBlue','Indigo','GoldenRod','FireBrick';
               'RoyalBlue','Indigo','GoldenRod','FireBrick';
               'CornflowerBlue','Indigo','GoldenRod','FireBrick';
               'CornflowerBlue','Indigo','GoldenRod','FireBrick';
               'CornflowerBlue','Indigo','GoldenRod','FireBrick';
               'CornflowerBlue','Indigo','GoldenRod','FireBrick';
               'DeepSkyBlue','Indigo','GoldenRod','FireBrick';
               'DeepSkyBlue','Indigo','GoldenRod','FireBrick';
               'DeepSkyBlue','Indigo','GoldenRod','FireBrick';
               'DeepSkyBlue','Indigo','GoldenRod','FireBrick';
               'CadetBlue','Indigo','GoldenRod','FireBrick';
               'CadetBlue','Indigo','GoldenRod','FireBrick';
               'CadetBlue','Indigo','GoldenRod','FireBrick';
               'CadetBlue','Indigo','GoldenRod','FireBrick'};
timescales = {'30 s','5 min','15 min','60 min'};
timescale_index = [1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4]; %which rows correspond to each timescale


% 1. calculate mean data from all replicates for plotting at end of next step
mean_tau = cellfun(@mean,tau_i_data);
mean_Vb = cellfun(@mean,Vb_i_data);
std_tau = cellfun(@std,tau_i_data);
std_Vb = cellfun(@std,Vb_i_data);
n_Vb = cellfun(@length,Vb_i_data);
sem_tau = std_tau./sqrt(n_Vb);
sem_Vb = std_Vb./sqrt(n_Vb);


% 2. loop through experimental replicates
for rr = 1:length(tau_i_data)
    
    clear curr_tau curr_Vb
    
    rep_tau = tau_i_data(rr,:);
    if isempty(rep_tau{1}) == 1
        continue
    end
    
    counter = counter + 1;
    rep_Vb = Vb_i_data(rr,:);
    
    
    % 2. plot scatter from current rep by nutrient condition (fluc, low, ave, high)
    for cond = 1:length(rep_tau)
        
        curr_tau = rep_tau{cond};
        curr_Vb = rep_Vb{cond};
        curr_color = rgb(palette_i{cond});
        
        figure(counter)
        scatter(curr_Vb,curr_tau,'MarkerFaceColor',curr_color,'MarkerEdgeColor',curr_color)
        hold on
        
    end
    clear cond curr_tau curr_Vb curr_color rep_Vb rep_tau
    
    
    % 3. overlay mean data (iteratively by condition for coloring)
    ss = size(sem_Vb);
    for ii = 1:ss(1)*ss(2)
        
        test = mean_tau(ii);
        if isnan(test) == 1
            continue
        end
        
        color = rgb(palette_rep{ii});
        yval = mean_tau(ii);
        xval = mean_Vb(ii);
        yerr = sem_tau(ii);
        xerr = sem_Vb(ii);
        
        
        figure(counter) % slope vs Vb (mean +/- sem)
        scatter(xval,yval,100,'MarkerFaceColor',color,'MarkerEdgeColor',color)
        hold on
        errorbar(xval,yval,xerr,'.','horizontal','Color',color)
        hold on
        errorbar(xval,yval,yerr,'.','Color',color)
        hold on
        ylabel('tau (min)')
        xlabel('Vb (cubic um)')
        axis([1 9 0 120])
        
        
        
    end
    clear test xval yval yerr xerr color
end


%% Part 5. save figures

for rep = 1:13
    figure(rep)
    saveas(gcf,strcat('ss22-indi-v-rep-',num2str(rep)),'epsc')
    close(gcf)
end




