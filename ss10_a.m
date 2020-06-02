%% ss10_a: tau_i vs Vb_i as an individual plot per replicate


%  Goal: compile replicate data from fluctuating and steady conditions and
%        plot tau_i vs Vb_i as a scatter
%        save slopes, intercepts, Vb, tau, and lambda per replicate


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. plot tau_i vs Vb_i per condition grouped by fluc timescale 
%          determine slope and other data for each STEADY replicate
%  Part 3. plot tau_i vs Vb_i per condition grouped by fluc timescale 
%          determine slope and other data for each FLUCTUATING replicate
%  Part 4. save data!
%  Part 5. calculate mean and stdev of slopes per condition, and
%          plot slope vs. growth rate from individual replicates


%  Last edit: Jen Nguyen, 2020 June 1
%  Commit: bug fix, stats calculation and plotting slope vs growth rate


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
palette_fluc = {'DarkOrchid','DodgerBlue','Chocolate','LimeGreen'};
palette_steady = {
                  'DarkSlateBlue','Chocolate','DarkRed';
                  'SlateBlue','Peru','LightCoral';
                  'DarkMagenta','DarkGoldenrod','IndianRed';
                  'DarkOrchid','Goldenrod','FireBrick'
                  };
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
clear t0_5 t5 t15 t60


%% Part 2. compile data from each steady replicate

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
    numcells = nan(4,3);
    yint = nan(4,3);
    slopes = nan(4,3);
    Vb_mean = nan(4,3);
    Vb_std= nan(4,3);
    tau_mean = nan(4,3);
    tau_std = nan(4,3);
    lambda_mean = nan(4,3);
    lambda_std = nan(4,3);
    
    
    % 1. loop through steady conditions and compile replicate data
    rows_steady = Tgroup{ts};
    color_counter = 0;
    
    for cond_ii = 1:length(conditions_steady)
        
        color_counter = color_counter + 1;
        
        % 2. isolate replicate data from each steady condition
        currCond = conditions_steady{cond_ii};
        
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
        
        
        % 4. plot scatter for each replicate condition
        for rep = 1:length(currData)
            
            % i. gather parameters
            repData = currData{rep,1};
            Vb = repData.cc(:,vol_birth);
            dt = repData.meta(:,tau);
            gr = repData.meta(:,lamb);
            
            
            % ii. keep only data within 95% of mean Vb
            Vb_meen = mean(Vb);
            sigma = std(Vb);
            
            Vb_trim1 = Vb(Vb < Vb_meen + 2*sigma);
            dt_trim1 = dt(Vb < Vb_meen + 2*sigma);
            gr_trim1 = gr(Vb < Vb_meen + 2*sigma);
            
            Vb_i = Vb_trim1(Vb_trim1 > Vb_meen - 2*sigma);
            tau_i = dt_trim1(Vb_trim1 > Vb_meen - 2*sigma);
            lambda_i = gr_trim1(Vb_trim1 > Vb_meen - 2*sigma);
            clear dt_trim1 Vb_trim1 gr_trim1 sigma
            clear repData Vb dt gr Vb_meen
            
            
            % iii. plot
            figure(ts)
            scatter(Vb_i,tau_i,'MarkerEdgeColor',rgb(palette_steady{rep,color_counter}))
            numcells(rep,color_counter) = length(Vb_i);
            
            % overlay fit
            fit = polyfit(Vb_i,tau_i,1);
            x = linspace(1,8,10);
            y = fit(1).*x + fit(2);
            hold on
            plot(x,y,'Color',rgb(palette_steady{rep,color_counter}))
            clear x y 
            
            % iv. store data
            slopes(rep,color_counter) = fit(1);
            yint(rep,color_counter) = fit(2);
            
            Vb_mean(rep,color_counter) = mean(Vb_i);
            Vb_std(rep,color_counter) = std(Vb_i);
            
            tau_mean(rep,color_counter) = mean(tau_i);
            tau_std(rep,color_counter) = std(tau_i);
            
            lambda_mean(rep,color_counter) = mean(lambda_i);
            lambda_std(rep,color_counter) = std(lambda_i);
            clear Vb_i tau_i lambda_i
            
        end 
        
        ts_stats_steady{ts}.numcells = numcells;
        ts_stats_steady{ts}.slopes = slopes;
        ts_stats_steady{ts}.yint = yint;
        ts_stats_steady{ts}.Vb_mean = Vb_mean;
        ts_stats_steady{ts}.Vb_std = Vb_std;
        ts_stats_steady{ts}.tau_mean = tau_mean;
        ts_stats_steady{ts}.tau_std = tau_std;
        ts_stats_steady{ts}.lambda_mean = lambda_mean;
        ts_stats_steady{ts}.lambda_std = lambda_std;
        
    end
    
    figure(ts) % timescale
    axis([1 8 0 140])
    title(strcat('steady data for T = ',num2str(environment_order{ts+1})))
    ylabel('tau_i (min)')
    xlabel('Vb_i')
    
end
clear yint slopes numcells fit color_counter ts cond_ii rep vol_birth lamb
clear Vb_mean Vb_std tau_mean tau_std lambda_mean lambda_std


%% Part 3. compile data from each fluctuating replicate


% 0. initialize parameters for which calculate stats in organized_data
metric = 1; % volume at birth = col in in cc
tau = 2;    % interdivision time = col 2 in meta
lamb = 1;   % mean growth rate = col 1 in meta


% 0. initialize data to be collected and stored
numcells_fluc = nan(4,4);
yint_fluc = nan(4,4);
slopes_fluc = nan(4,4);
Vb_mean_fluc = nan(4,4);
Vb_std_fluc = nan(4,4);
tau_mean_fluc = nan(4,4);
tau_std_fluc = nan(4,4);
lambda_mean_fluc = nan(4,4);
lambda_std_fluc = nan(4,4);
    

% 0. for each condition of interest
flucdata = 2:5;


for ts = 1:length(flucdata)
 
    
    % 1. isolate replicate data from each conditions
    col = flucdata(ts);
    
    currData = organized_data(:,col);
    noData = cellfun(@isempty,currData);
    currData = currData(noData == 0);
    clear noData
    
    
    % 2. compile data from each replicate
    for rep = 1:length(currData)
        
        % i. gather parameters
        repData = currData{rep,1};
        Vb = repData.cc(:,metric);
        dt = repData.meta(:,tau);
        gr = repData.meta(:,lamb);
        
        % ii. keep only data within 95% of mean Vb
        Vb_meen = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_meen + 2*sigma);
        dt_trim1 = dt(Vb < Vb_meen + 2*sigma);
        gr_trim1 = gr(Vb < Vb_meen + 2*sigma);
        
        Vb_i = Vb_trim1(Vb_trim1 > Vb_meen - 2*sigma);
        tau_i = dt_trim1(Vb_trim1 > Vb_meen - 2*sigma);
        lambda_i = gr_trim1(Vb_trim1 > Vb_meen - 2*sigma);
        clear dt_trim1 Vb_trim1 gr_trim1 sigma
        clear repData Vb dt gr Vb_meen
        
        
        % iii. plot
        figure(ts)
        scatter(Vb_i,tau_i,'MarkerEdgeColor',rgb(palette_fluc{rep}))
        numcells_fluc(rep,ts) = length(Vb_i);

        % overlay fit
        fit = polyfit(Vb_i,tau_i,1);
        x = linspace(1,8,10);
        y = fit(1).*x + fit(2);
        hold on
        plot(x,y,'Color',rgb(palette_fluc{rep}))
        clear x y
        
        
        % iv. store slope data
        slopes_fluc(rep,ts) = fit(1);
        yint_fluc(rep,ts) = fit(2);
        
        %Vb_mean(rep,color_counter) = mean(Vb_i);
        Vb_mean_fluc(rep,ts) = mean(Vb_i);
        Vb_std_fluc(rep,ts) = std(Vb_i);
        
        tau_mean_fluc(rep,ts) = mean(tau_i);
        tau_std_fluc(rep,ts) = std(tau_i);
        
        lambda_mean_fluc(rep,ts) = mean(lambda_i);
        lambda_std_fluc(rep,ts) = std(lambda_i);
        clear Vb_i tau_i lambda_i
        
    end
    
    figure(ts)
    axis([1 8 0 140])
    ylabel('tau_i (min)')
    xlabel('Vb_i')
    title(strcat('steady data for T = ',num2str(environment_order{ts+1})))
    
end
clear fit color_counter ts col rep


%% Part 4. save data!

save('ss10_a.mat', 'ts_stats_steady','numcells_fluc','slopes_fluc','yint_fluc','Vb_mean_fluc','Vb_std_fluc','tau_mean_fluc','tau_std_fluc','lambda_mean_fluc','lambda_std_fluc')


%% Part 5. calculate stats

clear
clc

% 0. initialize data
cd('/Users/jen/super-size/')
load('ss10_a.mat')



% 1. calculate mean, standard deviation and coefficient of variation 
%    between replicates from FLUCTUATING conditions

slopes.mean = nanmean(slopes_fluc);
slopes.std = nanstd(slopes_fluc);
slopes.cv = (slopes.std./slopes.mean) * 100;

yint.mean = nanmean(yint_fluc);
yint.std = nanstd(yint_fluc);
yint.cv = (yint.std./yint.mean) * 100;

Vb.mean = nanmean(Vb_mean_fluc);
Vb.std = nanstd(Vb_mean_fluc);
Vb.cv = (Vb.std./Vb.mean) * 100;

tau.mean = nanmean(tau_mean_fluc);
tau.std = nanstd(tau_mean_fluc);
tau.cv = (tau.std./tau.mean) * 100;

lambda.mean = nanmean(lambda_mean_fluc);
lambda.std = nanstd(lambda_mean_fluc);
lambda.cv = (lambda.std./lambda.mean) * 100;


% 2. calculate mean, standard deviation and coefficient of variation 
%    between replicates from STEADY conditions

% A. separated by parallel fluctuation timescale
%    i. birth volume
xx = nan(4,3);
yy = nan(4,3);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    xx(tscale,:) = nanmean(ts_stats_steady{1,tscale}.Vb_mean);
    yy(tscale,:) = nanstd(ts_stats_steady{1,tscale}.Vb_mean);
end
Vb.mean_steady_A = xx;
Vb.std_steady_A = yy;
Vb.cv_steady_A = (yy./xx) * 100;
clear xx yy tscale

%   ii. division time (tau)
aa = nan(4,3);
bb = nan(4,3);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    aa(tscale,:) = nanmean(ts_stats_steady{1,tscale}.tau_mean);
    bb(tscale,:) = nanstd(ts_stats_steady{1,tscale}.tau_mean);
end
tau.mean_steady_A = aa;
tau.std_steady_A = bb;
tau.cv_steady_A = (bb./aa) * 100;
clear aa bb tscale

%  iii. growth rate (lambda)
qq = nan(4,3);
rr = nan(4,3);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    qq(tscale,:) = nanmean(ts_stats_steady{1,tscale}.lambda_mean);
    rr(tscale,:) = nanstd(ts_stats_steady{1,tscale}.lambda_mean);
end
lambda.mean_steady_A = qq;
lambda.std_steady_A = rr;
lambda.cv_steady_A = (rr./qq) * 100;
clear qq rr tscale

% iv. slopes of tau_i vs. Vb_i plot
cc = nan(4,3);
dd = nan(4,3);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    cc(tscale,:) = nanmean(ts_stats_steady{1,tscale}.slopes);
    dd(tscale,:) = nanstd(ts_stats_steady{1,tscale}.slopes);
end
Vb.mean_steady_A = cc;
Vb.std_steady_A = dd;
Vb.cv_steady_A = (dd./cc) * 100;
clear cc dd tscale

% v. y-intercept of tau_i vs. Vb_i plot
ee = nan(4,3);
ff = nan(4,3);
for tscale = 1:4     % 1 = 30 sec; 2 = 5 min; 3 = 15 min; 4 = 60 min
    ee(tscale,:) = nanmean(ts_stats_steady{1,tscale}.yint);
    ff(tscale,:) = nanstd(ts_stats_steady{1,tscale}.yint);
end
Vb.mean_steady_A = ee;
Vb.std_steady_A = ff;
Vb.cv_steady_A = (ff./ee) * 100;
clear ee ff tscale


% B. compiled by steady nutrient condition
ss = [];
tt = [];
uu = [];
vv = [];
ww = [];
for tscale = 1:4
    ss = [ss; ts_stats_steady{1,tscale}.Vb_mean];     % Vb
    tt = [tt; ts_stats_steady{1,tscale}.tau_mean];    % tau
    uu = [uu; ts_stats_steady{1,tscale}.lambda_mean]; % lambda
    vv = [vv; ts_stats_steady{1,tscale}.slopes];
    ww = [ww; ts_stats_steady{1,tscale}.yint];
end
clear tscale

% Vb
Vb.mean_steady_B = nanmean(ss);
Vb.std_steady_B = nanstd(ss);
Vb.cv_steady_B = (Vb.std_steady_B./Vb.mean_steady_B) * 100;

% tau
tau.mean_steady_B = nanmean(tt);
tau.std_steady_B = nanstd(tt);
tau.cv_steady_B = (tau.std_steady_B./tau.mean_steady_B) * 100;

% lambda
lambda.mean_steady_B = nanmean(uu);
lambda.std_steady_B = nanstd(uu);
lambda.cv_steady_B = (lambda.std_steady_B./lambda.mean_steady_B) * 100;

% slopes of tau_i vs. Vb_i plot
slopes.mean_steady_B = nanmean(vv);
slopes.std_steady_B = nanstd(vv);
slopes.cv_steady_B = (slopes.std_steady_B./slopes.mean_steady_B) * 100;

% y-intercepts of tau_i vs. Vb_i plot
yint.mean_steady_B = nanmean(ww);
yint.std_steady_B = nanstd(ww);
yint.cv_steady_B = (yint.std_steady_B./yint.mean_steady_B) * 100;

clear ss tt uu vv ww
%clear lambda_mean_fluc lambda_std_fluc numcells_fluc slopes_fluc yint_fluc
%clear tau_mean_fluc tau_std_fluc Vb_mean_fluc Vb_std_fluc ts_stats_steady


%% Part 6. plot slopes (per replicate) vs. lambda, tau and Vb
%          also plot y-intercepts (per replicate vs. lambda, tau and Vb)

clear
clc

% 0. initialize data
cd('/Users/jen/super-size/')
load('ss10_a.mat')

% 0. initialize colors for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};
    

% 1. concatenate slope, lambda, tau and Vb data from steady conditions
low=1; ave=2; high=3; % column in ts_stats_steady

slopes_steady = [];
yints_steady = [];
lambda_mean_steady = [];
tau_mean_steady = [];
Vb_mean_steady = [];

lambda_std_steady = [];
tau_std_steady = [];
Vb_std_steady = [];

numcells_steady = [];
for ts = 1:length(ts_stats_steady)
    
    slopes_steady = [slopes_steady; ts_stats_steady{1,ts}.slopes];
    yints_steady = [yints_steady; ts_stats_steady{1,ts}.yint];
    lambda_mean_steady = [lambda_mean_steady; ts_stats_steady{1,ts}.lambda_mean];
    tau_mean_steady = [tau_mean_steady; ts_stats_steady{1,ts}.tau_mean];
    Vb_mean_steady = [Vb_mean_steady; ts_stats_steady{1,ts}.Vb_mean];
    
    lambda_std_steady = [lambda_std_steady; ts_stats_steady{1,ts}.lambda_std];
    tau_std_steady = [tau_std_steady; ts_stats_steady{1,ts}.tau_std];
    Vb_std_steady = [Vb_std_steady; ts_stats_steady{1,ts}.Vb_std];
    
    numcells_steady = [numcells_steady; ts_stats_steady{1,ts}.numcells];
    
end
clear ts



% 2. calculate s.e.m. for lambda, tau and Vb
lambda_sem_steady = lambda_std_steady./sqrt(numcells_steady);
lambda_sem_fluc = lambda_std_fluc./sqrt(numcells_fluc);

tau_sem_steady = tau_std_steady./sqrt(numcells_steady);
tau_sem_fluc = tau_std_fluc./sqrt(numcells_fluc);

Vb_sem_steady = Vb_std_steady./sqrt(numcells_steady);
Vb_sem_fluc = Vb_std_fluc./sqrt(numcells_fluc);




% 3. plot data from STEADY conditions
for sc = 1:high
    
    % i. define color by condition
    sc_color = rgb(palette_steady{sc});

    % ii. isolate data by condition column
    sc_slope = slopes_steady(:,sc);
    sc_yint = yints_steady(:,sc);
    
    sc_lambda_means = lambda_mean_steady(:,sc);
    sc_lambda_sem = lambda_sem_steady(:,sc);
    
    sc_tau_means = tau_mean_steady(:,sc);
    sc_tau_sem = tau_sem_steady(:,sc);
    
    sc_Vb_means = Vb_mean_steady(:,sc);
    sc_Vb_sem = Vb_sem_steady(:,sc);
    
    % iii. plot data
    figure(1) % slope vs lambda (mean +/- sem)
    scatter(sc_lambda_means,sc_slope,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_lambda_means,sc_slope,sc_lambda_sem,'.','horizontal','Color',sc_color)
    
    figure(2) % slope vs tau (mean +/- sem)
    scatter(sc_tau_means,sc_slope,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_tau_means,sc_slope,sc_tau_sem,'.','horizontal','Color',sc_color)
    
    figure(3) % slope vs Vb (mean +/- sem)
    scatter(sc_Vb_means,sc_slope,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_Vb_means,sc_slope,sc_Vb_sem,'.','horizontal','Color',sc_color)
    
    figure(4) % y-int vs lambda (mean +/- sem)
    scatter(sc_lambda_means,sc_yint,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_lambda_means,sc_yint,sc_lambda_sem,'.','horizontal','Color',sc_color)
    
    figure(5) % y-int vs tau (mean +/- sem)
    scatter(sc_tau_means,sc_yint,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_tau_means,sc_yint,sc_tau_sem,'.','horizontal','Color',sc_color)
    
    figure(6) % y-int vs Vb (mean +/- sem)
    scatter(sc_Vb_means,sc_yint,100,'filled','MarkerFaceColor',sc_color)
    hold on
    errorbar(sc_Vb_means,sc_yint,sc_Vb_sem,'.','horizontal','Color',sc_color)
    
end
figure(1)
xlabel('lambda (1/h)')
ylabel('slope (tau_i vs Vb_i)')

figure(2)
xlabel('tau (min)')
ylabel('slope (tau_i vs Vb_i)')

figure(3)
xlabel('Vb (cubic um)')
ylabel('slope (tau_i vs Vb_i)')

figure(4)
xlabel('lambda (1/h)')
ylabel('y-intercept (tau_i vs Vb_i)')

figure(5)
xlabel('tau (min)')
ylabel('y-intercept (tau_i vs Vb_i)')

figure(6)
xlabel('Vb (cubic um)')
ylabel('y-intercept (tau_i vs Vb_i)')

clear sc_slope sc_lambda_means sc_lambda_sem sc_tau_means sc_yint
clear sc_tau_sem sc_Vb_means sc_Vb_sem sc sc_color



% 4. plot data from FLUCTUATING conditions


for fc = 1:4
   
    % i. define color by condition
    fc_color = rgb(palette_fluc{fc});

    % ii. isolate data by condition column
    fc_slope = slopes_fluc(:,fc);
    fc_yint = yint_fluc(:,fc);
    
    fc_lambda_means = lambda_mean_fluc(:,fc);
    fc_lambda_sem = lambda_sem_fluc(:,fc);
    
    fc_tau_means = tau_mean_fluc(:,fc);
    fc_tau_sem = tau_sem_fluc(:,fc);
    
    fc_Vb_means = Vb_mean_fluc(:,fc);
    fc_Vb_sem = Vb_sem_fluc(:,fc);
    
    % iii. plot data
    figure(1) % slope vs lambda (mean +/- sem)
    hold on
    scatter(fc_lambda_means,fc_slope,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_lambda_means,fc_slope,fc_lambda_sem,'.','horizontal','Color',fc_color)
    
    figure(2) % slope vs tau (mean +/- sem)
    hold on
    scatter(fc_tau_means,fc_slope,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_tau_means,fc_slope,fc_tau_sem,'.','horizontal','Color',fc_color)

    figure(3) % slope vs Vb (mean +/- sem)
    hold on
    scatter(fc_Vb_means,fc_slope,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_Vb_means,fc_slope,fc_Vb_sem,'.','horizontal','Color',fc_color)
    
    figure(4) % y-int vs lambda (mean +/- sem)
    hold on
    scatter(fc_lambda_means,fc_yint,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_lambda_means,fc_yint,fc_lambda_sem,'.','horizontal','Color',fc_color)
    
    figure(5) % y-int vs tau (mean +/- sem)
    hold on
    scatter(fc_tau_means,fc_yint,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_tau_means,fc_yint,fc_tau_sem,'.','horizontal','Color',fc_color)

    figure(6) % y-int vs Vb (mean +/- sem)
    hold on
    scatter(fc_Vb_means,fc_yint,100,'filled','MarkerFaceColor',fc_color)
    hold on
    errorbar(fc_Vb_means,fc_yint,fc_Vb_sem,'.','horizontal','Color',fc_color)
    
    
end
clear fc fc_color fc_slope fc_lambda_means fc_lambda_sem
clear fc_tau_means fc_tau_sem fc_Vb_means fc_Vb_sem





