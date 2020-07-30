%% ss17: spread of single cell slope_i and replicate slope_i


%  Goal: compare slope_i between replicates with spread within a replicate


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. compile data from steady replicates
%  Part 3. compile data from fluctuating replicates
%  Part 4. save data!
%  Part 5. calculate replicate means for plotting in part 6
%  Part 6. plot!


%  Last edit: Jen Nguyen, 2020 July 29
%  Commit: update to plot single cell scatter for each replicate condition


%  OK let's go!

%% Part 0. initialize data, saved from figure1A_division.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')


% 0. initialize plotting parameters
environment_order = {'low',30,300,900,3600,'ave','high'};


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
clear t0_5 t5 t15 t60


%% Part 2. compile data from steady replicates

clear low ave high fluc

% 0. initialize parameters for which calculate stats in organized_data
vol_birth = 1;     % volume at birth = col in cc (compiled in figure1A_division.m)
tau = 2;    % interdivision time = col 2 in meta
lamb = 1;   % mean growth rate = col 1 in meta


% 0. group steady rows by timescale of fluctuating condition
conditions_steady = {'low','ave','high'};
Tgroup = {1:3; 4:6; 7:10; 11:13}; % T = 30 s; 5 min; 15 min; 60 min


for ts = 1:length(Tgroup) % timescale groups
    
    
    % 0. initialize data for storage
    numcells_steady = nan(4,3);
    compiled_slopes_i = cell(4,3);
    compiled_Vb_i = cell(4,3);
    compiled_tau_i = cell(4,3);
    compiled_lambda_i = cell(4,3);
    
    
    % 1. loop through steady conditions and compile replicate data
    rows_steady = Tgroup{ts};
    %color_counter = 0;
    
    for sc = 1:length(conditions_steady)
        
        %color_counter = color_counter + 1;
        
        % 2. isolate replicate data from each steady condition
        currCond = conditions_steady{sc};
        
        if strcmp(currCond,'low') == 1
            col = 1;
        elseif strcmp(currCond,'ave') == 1
            col = 6;
        elseif strcmp(currCond,'high') == 1
            col = 7;
        end
        condData = organized_data(:,col);
        clear col
        
        
        % 3. use only steady data from experiments of timescale ts_current
        currData = condData(rows_steady);
        noData = cellfun(@isempty,currData);
        currData = currData(noData == 0);
        clear noData
        
        
        % from each replicate condition...
        for rep = 1:length(currData)
            
            % 4. gather single-cell Vb, tau and lambda
            repData = currData{rep,1};
            Vb = repData.cc(:,vol_birth);
            dt = repData.meta(:,tau);
            gr = repData.meta(:,lamb);
            
            
            % 5. calculate single cell and mean of replicate slope_i 
            slope = dt./Vb;
            
            
            % 6. remove data outside of 3 st dev from mean slope (also in ss16.m)
            mm = mean(slope);
            thresh = 3 * std(slope);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = slope > upper;
            cut_down = slope < lower;
            toCut = cut_up + cut_down;

            slope_i = slope(toCut == 0);
            tau_i = dt(toCut == 0);
            Vb_i = Vb(toCut == 0);
            lambda_i = gr(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear Vb dt slope gr repData
            
            
            % 7. store data
            numcells_steady(rep,sc) = length(slope_i);
            compiled_slopes_i{rep,sc} = slope_i;
            compiled_Vb_i{rep,sc} = Vb_i;
            compiled_tau_i{rep,sc} = tau_i;
            compiled_lambda_i{rep,sc} = lambda_i;
            clear Vb_i tau_i lambda_i slope_i
            
        end 
                
    end
    clear sc rep
    
    compiled_steady{ts}.slope_i = compiled_slopes_i;
    compiled_steady{ts}.Vb_i = compiled_Vb_i;
    compiled_steady{ts}.tau_i = compiled_tau_i;
    compiled_steady{ts}.lambda_i = compiled_lambda_i;
    
    
end
clear numcells lamb tau vol_birth exp ts currCond
clear compiled_slopes_i compiled_Vb_i compiled_tau_i compiled_lambda_i


%% Part 3. compile data from fluctuating replicates


% 0. initialize parameters for which calculate stats in organized_data
metric = 1; % volume at birth = col in in cc
vol_birth = 1;     % volume at birth = col in cc (compiled in figure1A_division.m)
tau = 2;    % interdivision time = col 2 in meta
lamb = 1;   % mean growth rate = col 1 in meta


% 0. initialize data to be collected and stored
numcells_fluc = nan(4,4);
fluc_slopes_i = cell(4,4);
fluc_Vb_i = cell(4,4);
fluc_tau_i = cell(4,4);
fluc_lambda_i = cell(4,4);
    

% 0. for each condition of interest
flucdata = 2:5;


for ts = 1:length(flucdata)
 
    
    % 1. isolate replicate data from each conditions
    col = flucdata(ts);
    
    currData = organized_data(:,col);
    noData = cellfun(@isempty,currData);
    currData = currData(noData == 0);
    clear noData
    
    
    % from each replicate ...
    for rr = 1:length(currData)
        
        % 2. gather single-cell Vb, tau and lambda
        rrData = currData{rr,1};
        Vbb = rrData.cc(:,vol_birth);
        dtt = rrData.meta(:,tau);
        grr = rrData.meta(:,lamb);

        
        % 3. calculate single cell and mean of replicate slope_i
        slopee = dtt./Vbb;
        
        
        % 4. remove data outside of 3 st dev from mean slope (also in ss16.m)
        mmm = mean(slopee);
        threshh = 3 * std(slopee);
        upper = mmm + threshh;
        lower = mmm - threshh;
        cut_up = slopee > upper;
        cut_down = slopee < lower;
        toCut = cut_up + cut_down;
        
        slope_i = slopee(toCut == 0);
        tau_i = dtt(toCut == 0);
        Vb_i = Vbb(toCut == 0);
        lambda_i = grr(toCut == 0);
        clear mmm threshh upper lower cut_up cut_down toCut
        clear Vbb dtt slopee grr rrData
        
        
        % 5. store data
        numcells_fluc(rr,ts) = length(slope_i);
        fluc_slopes_i{rr,ts} = slope_i;
        fluc_Vb_i{rr,ts} = Vb_i;
        fluc_tau_i{rr,ts} = tau_i;
        fluc_lambda_i{rr,ts} = lambda_i;
        clear Vb_i tau_i lambda_i slope_i
        
      
    end
    clear rr
    
end
    
clear numcells lamb tau vol_birth exp ts currCond col


%% Part 4. save data!

save('ss17.mat', 'compiled_steady','numcells_fluc','fluc_slopes_i','fluc_Vb_i','fluc_tau_i','fluc_lambda_i')


%% Part 5. calculate replicate means for plotting in part 6

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


%% Part 6. plot replicate means only

%          mean slope_i vs. mean Vb_i of each replicate


% 0. initialize colors and shape for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};
shape = 'o';


% 1. first, plot mean slope_i vs. mean Vb_i from all steady replicates

%   i. gather mean and std of slope_i and Vb_i into columns by condition
mean_slope_steady = [];
mean_Vb_steady = [];
std_slope_steady = [];
std_Vb_steady = [];
numcells_steady = [];
steady_slopes_i = [];
steady_Vb_i = [];

for ii = 1:length(steady_slope.mean)
    
    % ii. loop through timescale groups to concatenate steady means
    ts_slope = steady_slope.mean{1,ii};
    mean_slope_steady = [mean_slope_steady; ts_slope];
    
    ts_Vb = steady_Vb.mean{1,ii};
    mean_Vb_steady = [mean_Vb_steady; ts_Vb];
    
    s_std = steady_slope.std{1,ii};
    std_slope_steady = [std_slope_steady; s_std];
    
    v_std = steady_Vb.std{1,ii};
    std_Vb_steady = [std_Vb_steady; v_std];
    
    num_ii = cellfun(@length,compiled_steady{1,ii}.slope_i);
    numcells_steady = [numcells_steady; num_ii];
    
    slo_ii = compiled_steady{1,ii}.slope_i;
    steady_slopes_i = [steady_slopes_i; slo_ii];
    
    vb_ii = compiled_steady{1,ii}.Vb_i;
    steady_Vb_i = [steady_Vb_i; vb_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ts_Vb ts_slope s_std v_std ii num_ii vb_ii slo_ii


% iii. calculate standard error of the mean for each steady replicate
%      s.e.m. = st dev / sqrt(n)
sem_slope_steady = std_slope_steady./sqrt(numcells_steady);
sem_Vb_steady =  std_Vb_steady./sqrt(numcells_steady);


% iv. plot data from STEADY conditions
%     color hollow points by condition

for sc = 1:3
    
    % i. define color by condition
    sc_color = rgb(palette_steady{sc});
    
    % ii. isolate data by condition column
    yvals = mean_slope_steady(:,sc);
    xvals = mean_Vb_steady(:,sc);
    
    yerr = sem_slope_steady(:,sc);
    xerr = sem_Vb_steady(:,sc);

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
ylabel('slope_i')
xlabel('Vb_i')



% 2. second, plot mean slope_i vs. mean Vb_i from all fluctuating replicates

%   i. gather mean and std of slope_i and Vb_i into columns by timescale
mean_slope_fluc = fluc_slope.mean;
mean_Vb_fluc = fluc_Vb.mean;
std_slope_fluc = fluc_slope.std;
std_Vb_fluc = fluc_Vb.std;

% ii. calculate standard error of the mean for each steady replicate
%     s.e.m. = st dev / sqrt(n)
sem_slope_fluc = std_slope_fluc./sqrt(numcells_fluc);
sem_Vb_fluc =  std_Vb_fluc./sqrt(numcells_fluc);


% iii. plot data from FLUCTUATING conditions
%      color hollow points by timescale

for ts = 1:4
    
    % i. define color by condition
    color = rgb(palette_fluc{ts});
    
    % ii. isolate data by condition column
    yvals = mean_slope_fluc(:,ts);
    xvals = mean_Vb_fluc(:,ts);
    
    yerr = sem_slope_fluc(:,ts);
    xerr = sem_Vb_fluc(:,ts);

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


%% Part 7. plot single-cell scatter overlaid by replicate means

%         in separate figures, overlay single rep condition scatter over mean data
%         filling in 'o' for replicate mean of scatter
%         other replicate means can have hollow 'o'

% 1. loop through steady replicates to plot scatter from one rep per figure
%    by column nutrient condition (low, ave, high)

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss17/')
counter = 0;

for condition = 1:3
    for row = 1:length(steady_Vb_i)
         
        % 2. isolate slope_i and Vb_i of current rep condition
        Vb_i = steady_Vb_i{row,condition};
        if isempty(Vb_i) == 1
            continue
        else
            slope_i = steady_slopes_i{row,condition};
            counter = counter + 1;
            
            % 3. plot single cell scatter
            figure(2)
            scatter(Vb_i,slope_i,10,'MarkerEdgeColor',rgb('Silver'))
            hold on
            
            % 4. plot mean point for current rep condition (filled)
            figure(2) % slope vs Vb (mean +/- sem)
            rc_mean_vb = mean_Vb_steady(row,condition);
            rc_mean_slo = mean_slope_steady(row,condition);
            rc_color = rgb(palette_steady{condition});
            scatter(rc_mean_vb,rc_mean_slo,100,'filled','MarkerFaceColor',rc_color)
            hold on
            
            
            % 5. plot remaining means of all replicate conditions (hollow)
            for sc = 1:3
                
                % i. define color by condition
                sc_color = rgb(palette_steady{sc});
                
                % ii. isolate data by condition column
                yvals = mean_slope_steady(:,sc);
                xvals = mean_Vb_steady(:,sc);
                
                yerr = sem_slope_steady(:,sc);
                xerr = sem_Vb_steady(:,sc);
                
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
            ylabel('slope_i')
            xlabel('Vb_i')
            axis([0 9 0 70])
            
            figure(2)
            saveas(gcf,strcat('ss17-rc-',num2str(counter)),'epsc')
            close(gcf)
        end
        
    end
end
clear condition rc_color row rc_mean_slo rc_mean_vb Vb_i slope_i



% 6. loop through fluctuating replicates to plot scatter from one rep per figure
%    by column timescale (30 s, 5 min, 15 min, 60 min)


for timescale = 1:4
    for rw = 1:length(fluc_slopes_i)
        
        % 7. isolate slope_i and Vb_i of current rep condition
        Vb_if = fluc_Vb_i{rw,timescale};
        if isempty(Vb_if) == 1
            continue
        else
            slope_if = fluc_slopes_i{rw,timescale};
            counter = counter + 1;
            
            % 8. plot single cell scatter
            figure(3)
            scatter(Vb_if,slope_if,10,'MarkerEdgeColor',rgb('Silver'))
            hold on
            
            % 9. plot mean point for current rep condition (filled)
            figure(3) % slope vs Vb (mean +/- sem)
            rt_mean_vb = mean_Vb_fluc(rw,timescale);
            rt_mean_slo = mean_slope_fluc(rw,timescale);
            rt_color = rgb(palette_fluc{timescale});
            scatter(rt_mean_vb,rt_mean_slo,100,'filled','MarkerFaceColor',rt_color)
            hold on
            
            
            % 10. plot remaining means of all replicate conditions (hollow)
            for fl = 1:4
                
                % i. define color by condition
                fl_color = rgb(palette_fluc{fl});
                
                % ii. isolate data by condition column
                yvals = mean_slope_fluc(:,fl);
                xvals = mean_Vb_fluc(:,fl);
                
                yerr = sem_slope_fluc(:,fl);
                xerr = sem_Vb_fluc(:,fl);
                
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
            ylabel('slope_i')
            xlabel('Vb_i')
            axis([0 9 0 70])
            
            
            % 10. plot means of all steady replicate conditions (hollow)
            for sc = 1:3
                
                % i. define color by condition
                sc_color = rgb(palette_steady{sc});
                
                % ii. isolate data by condition column
                yvals = mean_slope_steady(:,sc);
                xvals = mean_Vb_steady(:,sc);
                
                yerr = sem_slope_steady(:,sc);
                xerr = sem_Vb_steady(:,sc);
                
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
            saveas(gcf,strcat('ss17-rc-',num2str(counter)),'epsc')
            close(gcf)
            
        end
    end
    
end




