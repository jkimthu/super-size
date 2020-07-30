%% ss18: trends in single cell and replicate Vb vs lambda


%  Goal: compare trend between replicates with trend within a replicate


%  Strategy: 

%  Part 1. calculate replicate means for plotting in part 2
%  Part 2. plot replicate means only
%  Part 3. plot single-cell scatter overlaid by replicate means


%  Last edit: Jen Nguyen, 2020 July 29
%  Commit: first commit, single cell vs average trends in Vb vs lambda


%  OK let's go!


%% Part 1. calculate replicate means for plotting in part 6

clc
clear

cd('/Users/jen/super-size/')
load('ss17.mat')


% 1. calculate mean, standard deviation and coefficient of variation 
%    of each replicates from FLUCTUATING conditions

fluc_slope.mean = cellfun(@nanmean,fluc_slopes_i);
fluc_slope.std = cellfun(@nanstd,fluc_slopes_i);
fluc_slope.cv = (fluc_slope.std./fluc_slope.mean) * 100;

fluc_Vb.mean = cellfun(@nanmean,fluc_Vb_i);
fluc_Vb.std = cellfun(@nanstd,fluc_Vb_i);
fluc_Vb.cv = (fluc_Vb.std./fluc_Vb.mean) * 100;

fluc_tau.mean = cellfun(@nanmean,fluc_tau_i);
fluc_tau.std = cellfun(@nanstd,fluc_tau_i);
fluc_tau.cv = (fluc_tau.std./fluc_tau.mean) * 100;

fluc_lambda.mean = cellfun(@nanmean,fluc_lambda_i);
fluc_lambda.std = cellfun(@nanstd,fluc_lambda_i);
fluc_lambda.cv = (fluc_lambda.std./fluc_lambda.mean) * 100;


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
qq = cell(1,4);
rr = cell(1,4);
ss = cell(1,4);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    q = cellfun(@nanmean,compiled_steady{1,tscale}.lambda_i);
    r = cellfun(@nanstd,compiled_steady{1,tscale}.lambda_i);
    qq{tscale} = q;
    rr{tscale} = r;
    ss{tscale} = (r./q) * 100;
end
steady_lambda.mean = qq;
steady_lambda.std = rr;
steady_lambda.cv = ss;
clear q r qq rr ss tscale


% iv. slopes (tau_i/Vb_i)
dd = cell(1,4);
ee = cell(1,4);
ff = cell(1,4);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    d = cellfun(@nanmean,compiled_steady{1,tscale}.slope_i);
    e = cellfun(@nanstd,compiled_steady{1,tscale}.slope_i);
    dd{tscale} = d;
    ee{tscale} = e;
    ff{tscale} = (e./d) * 100;
end
steady_slope.mean = dd;
steady_slope.std = ee;
steady_slope.cv = ff;
clear d e dd ee ff tscale



% 3. calculate mean, standard deviation and coefficient of variation 
%    from compilation of STEADY conditions

%    i. aggregate all steady data into columns by condition
gg = [];
hh = [];
ii = [];
jj = [];
for tscale = 1:4
    gg = [gg; compiled_steady{1,tscale}.Vb_i];     % Vb
    hh = [hh; compiled_steady{1,tscale}.tau_i];    % tau
    ii = [ii; compiled_steady{1,tscale}.lambda_i]; % lambda
    jj = [jj; compiled_steady{1,tscale}.slope_i];
end
clear tscale

%   ii. aggregate all data belonging to same condition
for condition = 1:3
    
    all_vb = [];
    all_tau = [];
    all_lam = [];
    all_slo = [];
    for repl = 1:length(gg)
        
        cr_vb = gg{repl,condition};
        if isempty(cr_vb) == 1
            continue
        else
            cr_tau = hh{repl,condition};
            cr_lam = ii{repl,condition};
            cr_slo = jj{repl,condition};
            
            all_vb = [all_vb; cr_vb];
            all_tau = [all_tau; cr_tau];
            all_lam = [all_lam; cr_lam];
            all_slo = [all_slo; cr_slo];
        end
        clear cr_vb cr_tau cr_lam cr_slo
        
    end
    aggregate_steady_Vb{1,condition} = all_vb;
    aggregate_steady_tau{1,condition} = all_tau;
    aggregate_steady_lambda{1,condition} = all_lam;
    aggregate_steady_slope{1,condition} = all_slo;
    
end
clear all_vb all_tau all_lam all_slo
clear cr_vb cr_tau cr_lam cr_slo
clear gg hh ii jj repl condition


%% Part 2. plot replicate means only

%          mean mean Vb_i vs. mean lambda of each replicate


% 0. initialize colors and shape for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};
shape = 'o';


% 1. first, plot mean Vb_i vs. mean lambda_i from all steady replicates

%   i. gather mean and std of lambda_i and Vb_i into columns by condition
mean_lambda_steady = [];
mean_Vb_steady = [];
std_lambda_steady = [];
std_Vb_steady = [];
numcells_steady = [];
steady_lambda_i = [];
steady_Vb_i = [];

for ii = 1:length(steady_lambda.mean)
    
    % ii. loop through timescale groups to concatenate steady means
    ts_lambda = steady_lambda.mean{1,ii};
    mean_lambda_steady = [mean_lambda_steady; ts_lambda];
    
    ts_Vb = steady_Vb.mean{1,ii};
    mean_Vb_steady = [mean_Vb_steady; ts_Vb];
    
    s_std = steady_lambda.std{1,ii};
    std_lambda_steady = [std_lambda_steady; s_std];
    
    v_std = steady_Vb.std{1,ii};
    std_Vb_steady = [std_Vb_steady; v_std];
    
    num_ii = cellfun(@length,compiled_steady{1,ii}.lambda_i);
    numcells_steady = [numcells_steady; num_ii];
    
    gr_ii = compiled_steady{1,ii}.lambda_i;
    steady_lambda_i = [steady_lambda_i; gr_ii];
    
    vb_ii = compiled_steady{1,ii}.Vb_i;
    steady_Vb_i = [steady_Vb_i; vb_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ts_Vb ts_lambda s_std v_std ii num_ii vb_ii gr_ii


% iii. calculate standard error of the mean for each steady replicate
%      s.e.m. = st dev / sqrt(n)
sem_lambda_steady = std_lambda_steady./sqrt(numcells_steady);
sem_Vb_steady =  std_Vb_steady./sqrt(numcells_steady);


% iv. plot data from STEADY conditions
%     color hollow points by condition

for sc = 1:3
    
    % i. define color by condition
    sc_color = rgb(palette_steady{sc});
    
    % ii. isolate data by condition column
    xvals = mean_lambda_steady(:,sc);
    yvals = mean_Vb_steady(:,sc);
    
    xerr = sem_lambda_steady(:,sc);
    yerr = sem_Vb_steady(:,sc);

    % iii. plot condition data
    figure(1) % slope vs Vb (mean +/- sem)
    scatter(xvals,yvals,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(xvals,yvals,xerr,'.','horizontal','Color',sc_color)
    hold on
    errorbar(xvals,yvals,yerr,'.','Color',sc_color)
    hold on
    
end
clear sc sc_color
ylabel('Vb_i')
xlabel('lambda_i')



% 2. second, plot mean slope_i vs. mean Vb_i from all fluctuating replicates

%   i. gather mean and std of slope_i and Vb_i into columns by timescale
mean_lambda_fluc = fluc_lambda.mean;
mean_Vb_fluc = fluc_Vb.mean;
std_lambda_fluc = fluc_lambda.std;
std_Vb_fluc = fluc_Vb.std;

% ii. calculate standard error of the mean for each steady replicate
%     s.e.m. = st dev / sqrt(n)
sem_lambda_fluc = std_lambda_fluc./sqrt(numcells_fluc);
sem_Vb_fluc =  std_Vb_fluc./sqrt(numcells_fluc);


% iii. plot data from FLUCTUATING conditions
%      color hollow points by timescale

for ts = 1:4
    
    % i. define color by condition
    color = rgb(palette_fluc{ts});
    
    % ii. isolate data by condition column
    xvals = mean_lambda_fluc(:,ts);
    yvals = mean_Vb_fluc(:,ts);
    
    xerr = sem_lambda_fluc(:,ts);
    yerr = sem_Vb_fluc(:,ts);

    % iii. plot condition data
    figure(1) % slope vs Vb (mean +/- sem)
    scatter(xvals,yvals,100,'filled','MarkerFaceColor',color)
    hold on
    errorbar(xvals,yvals,xerr,'.','horizontal','Color',color)
    hold on
    errorbar(xvals,yvals,yerr,'.','Color',color)
    hold on
    
end
clear color ts yvals xvals yerr xerr
title('mean slope_i vs mean Vb_i of each replicate')


%% Part 3. plot single-cell scatter overlaid by replicate means

%         in separate figures, overlay single rep condition scatter over mean data
%         filling in 'o' for replicate mean of scatter
%         other replicate means can have hollow 'o'

% 1. loop through steady replicates to plot scatter from one rep per figure
%    by column nutrient condition (low, ave, high)

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss18/')
counter = 0;

for condition = 1:3
    for row = 1:length(steady_Vb_i)
         
        % 2. isolate slope_i and Vb_i of current rep condition
        Vb_i = steady_Vb_i{row,condition};
        if isempty(Vb_i) == 1
            continue
        else
            lambda_i = steady_lambda_i{row,condition};
            counter = counter + 1;
            
            % 3. plot single cell scatter
            figure(2)
            scatter(lambda_i,Vb_i,10,'MarkerEdgeColor',rgb('Silver'))
            hold on
            
            % 4. plot mean point for current rep condition (filled)
            figure(2) % slope vs Vb (mean +/- sem)
            rc_mean_vb = mean_Vb_steady(row,condition);
            rc_mean_gr = mean_lambda_steady(row,condition);
            rc_color = rgb(palette_steady{condition});
            scatter(rc_mean_gr,rc_mean_vb,100,'filled','MarkerFaceColor',rc_color)
            hold on
            
            
            % 5. plot remaining means of all replicate conditions (hollow)
            for sc = 1:3
                
                % i. define color by condition
                sc_color = rgb(palette_steady{sc});
                
                % ii. isolate data by condition column
                xvals = mean_lambda_steady(:,sc);
                yvals = mean_Vb_steady(:,sc);
                
                xerr = sem_lambda_steady(:,sc);
                yerr = sem_Vb_steady(:,sc);
                
                % iii. plot condition data
                figure(2) % slope vs Vb (mean +/- sem)
                scatter(xvals,yvals,100,'MarkerEdgeColor',sc_color)
                hold on
                errorbar(xvals,yvals,xerr,'.','horizontal','Color',sc_color)
                hold on
                errorbar(xvals,yvals,yerr,'.','Color',sc_color)
                hold on
                
            end
            clear sc sc_color xvals yvals yerr xerr
            xlabel('lambda_i')
            ylabel('Vb_i')
            axis([0 4 0 7])
            
            figure(2)
            saveas(gcf,strcat('ss18-rc-',num2str(counter)),'epsc')
            close(gcf)
        end
        
    end
end
clear condition rc_color row rc_mean_gr rc_mean_vb Vb_i lambda_i



% 6. loop through fluctuating replicates to plot scatter from one rep per figure
%    by column timescale (30 s, 5 min, 15 min, 60 min)


for timescale = 1:4
    for rw = 1:length(fluc_lambda_i)
        
        % 7. isolate slope_i and Vb_i of current rep condition
        Vb_if = fluc_Vb_i{rw,timescale};
        if isempty(Vb_if) == 1
            continue
        else
            lambda_if = fluc_lambda_i{rw,timescale};
            counter = counter + 1;
            
            % 8. plot single cell scatter
            figure(3)
            scatter(lambda_if,Vb_if,10,'MarkerEdgeColor',rgb('Silver'))
            hold on
            
            % 9. plot mean point for current rep condition (filled)
            figure(3) % slope vs Vb (mean +/- sem)
            rt_mean_vb = mean_Vb_fluc(rw,timescale);
            rt_mean_gr = mean_lambda_fluc(rw,timescale);
            rt_color = rgb(palette_fluc{timescale});
            scatter(rt_mean_gr,rt_mean_vb,100,'filled','MarkerFaceColor',rt_color)
            hold on
            
            
            % 10. plot remaining means of all replicate conditions (hollow)
            for fl = 1:4
                
                % i. define color by condition
                fl_color = rgb(palette_fluc{fl});
                
                % ii. isolate data by condition column
                xvals = mean_lambda_fluc(:,fl);
                yvals = mean_Vb_fluc(:,fl);
                
                xerr = sem_lambda_fluc(:,fl);
                yerr = sem_Vb_fluc(:,fl);
                
                % iii. plot condition data
                figure(3) % slope vs Vb (mean +/- sem)
                scatter(xvals,yvals,100,'MarkerEdgeColor',fl_color)
                hold on
                errorbar(xvals,yvals,xerr,'.','horizontal','Color',fl_color)
                hold on
                errorbar(xvals,yvals,yerr,'.','Color',fl_color)
                hold on
                
            end
            clear fl fl_color xvals yvals yerr xerr
            xlabel('lambda_i')
            ylabel('Vb_i')
            axis([0 4 0 7])
            
            
            % 10. plot means of all steady replicate conditions (hollow)
            for sc = 1:3
                
                % i. define color by condition
                sc_color = rgb(palette_steady{sc});
                
                % ii. isolate data by condition column
                xvals = mean_lambda_steady(:,sc);
                yvals = mean_Vb_steady(:,sc);
                
                xerr = sem_lambda_steady(:,sc);
                yerr = sem_Vb_steady(:,sc);
                
                % iii. plot condition data
                figure(3) % slope vs Vb (mean +/- sem)
                scatter(xvals,yvals,100,'MarkerEdgeColor',sc_color)
                hold on
                errorbar(xvals,yvals,xerr,'.','horizontal','Color',sc_color)
                hold on
                errorbar(xvals,yvals,yerr,'.','Color',sc_color)
                hold on
                
            end
            clear sc sc_color xvals yvals yerr xerr
            
            figure(3)
            saveas(gcf,strcat('ss18-rc-',num2str(counter)),'epsc')
            close(gcf)
            
        end
    end
    
end




