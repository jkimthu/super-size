%%  figure 6 - confirm timescale of measured mu function produces shift


%  Goal: what feature of the growth rate signals produces the difference in
%        cell size between T = 15 and 60 min?

%  Result: it seems that the timescale makes a difference!
%          constraining simulations to the same timestep

%  Strategy:
%
%        Take raw mu function and slowly simplify various elements until
%        over simplification of a specific feature loses size shift


%  Last edit: jen, 2020 Feb 17
%  commit: first commit, rescaled mu functions with tstep = 5 sec


% Okie, go go let's go!

%% Part 0. initialize analysis

clear
clc

% 0. initialize stored mu functions
load('mu_functions')     % generated and saved from figure6_simulations
colorSpectrum = {'DodgerBlue','Indigo'};


% 0. initialize timescales of interest
T = [900; 3600];


% 0. initialize steady-state generations
steadies = 50:1000;


% 1. from replicates, calculate mean mu(t) for each timescale
mufun15 = mean(signals_15,length(T));
mufun60 = mean(signals_60,length(T));
mufun(1:30,1) = mufun15;
mufun(1:30,2) = mufun60;


% 2. plot for visual confirmation
% figure(1) % mean of replicate means
% plot(1:length(mufun15),mufun15,'Color',rgb(colorSpectrum(1)),'LineWidth',1,'Marker','.')
% hold on
% plot(1:length(mufun60),mufun60,'Color',rgb(colorSpectrum(2)),'LineWidth',1,'Marker','.')
% ylabel('Mu (1/h)')
% xlabel('Period fraction')
% legend('T = 15','T = 60')
% title('growth rate function over time')

clear signal_e exp_type


%% Part 1. artificially increase mu function resolution to 1 data point per 5 seconds

goalStep = 5; % seconds

for tscale = 1:length(T)
    
    % 1. define current growth rate function (T=15 or 60) and timescale
    mu_function = mufun(:,tscale);
    
    % 2. determine current timestep between values
    current_T = T(tscale);
    tstep_sec = (current_T)/length(mu_function); % timestep in sec
    
    % 3. add points between each measured step
    newRes = tstep_sec/goalStep;
    mufun_res = [];
    for pp = 1:length(mu_function)
        
        if pp == length(mu_function)
            x1 = mu_function(pp);
            x2 = mu_function(1);
        else
            x1 = mu_function(pp);
            x2 = mu_function(pp+1);
        end
        newPoints = linspace(x1,x2,newRes+1);
        mufun_res = [mufun_res, newPoints(1:newRes)];
        
    end
    clear x1 x2 pp newRes
    
    new_function{tscale} = mufun_res;
    
end
clear tscale tstep_sec current_T mufun_res


%% Part 2. simulate comparison between mu(t) using equal timestep


% 0. intialize number of cell generations simulated
gen = 1000;


% 0. initialize vectors for data storage
tau_compiled = nan(gen,length(T));
Vd_compiled = nan(gen,length(T));
Vb_compiled = nan(gen,length(T));


for tscale = 1:length(T)
    
    % 0. determine current timescale
    current_T = T(tscale);
    
    % 0. initialize mu function for current timescale
    mu_function = new_function{tscale};
    
    % 0. start generation counter
    g_counter = 1;
    
    % 0. timestep for each timescale
    %tstep = (current_T/3600)/length(mu_function); % timestep in hour
    tstep = 5/3600; % timestep in hour
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 2; %
    tau = -(17/60)*Vi + (64/60);
    period_frac = 1; % time counter for growth rate
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    
    % 1. simulate!
    counter_loop = 0; % each loop is another timestep
    while g_counter < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % i. store current cell cycle Vi
            Vb_compiled(g_counter,tscale) = Vi;
        end
        
        % ii. calculate exponent
        muT = mu_function(period_frac)*tstep; % mu(t)*t in hours
        
        % iii. calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % iv. calculate time into cell cycle
        t_cc = t_cc + tstep;
        
        % v. store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = tstep;
        else
            t(counter_loop) = t(counter_loop-1) + tstep;
        end
        
        % vi. determine whether cell divides
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_compiled(g_counter,tscale) = Vt;
            tau_compiled(g_counter,tscale) = tau;
            
            % re-set parameters for next cell cycle
            g_counter = g_counter + 1;
            Vi = Vt/2;
            tau = -(17/60)*Vi + (64/60);
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
        % vii. determine growth rate value of next timestep
        if period_frac == length(mu_function)
            period_frac = 1;
        else
            period_frac = period_frac + 1;
        end
        
    end
    
    % 2. plot volumetric growth over time!
    figure(2)
    subplot(2,1,tscale)
    hold on
    plot(t,V,'Color',rgb(colorSpectrum{tscale}))
    axis([0 40 0.5 3])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('T=',num2str(T(tscale))))
    
    
    % 3. plot histogram of birth volume
    figure(3)
    hold on
    histogram(Vb_compiled(steadies,tscale),'Facecolor',rgb(colorSpectrum{tscale}))
    ylabel('Cell count')
    xlabel('Birth volume')
    
end
clear tscale Vi Vt V t t_cc tau tstep g_counter counter_loop

figure(3)
legend('15 min','60 min')
    

%% Part 3. simulate using mu(t) from T = 15 at both timescales

% manipulates measured signal from T = 15 min experiments such that it
% stretches over a 60 min period in simulations.

% 0. intialize timescale from which measured function is drawn
mfun = 1;
oriT = T(mfun);
color15 = {'DodgerBlue','CornFlowerBlue'};

% 0. intialize number of cell generations simulated
gen = 1000;

% 0. initialize mu function for current timescale
mu_function = new_function{mfun};


% 1. prepare stretched function in which additional points increase the
%    length of function to longer T
newRes = 720/length(mu_function);
stretchedFun = [];
for pp = 1:length(mu_function)
    
    if pp == length(mu_function)
        x1 = mu_function(pp);
        x2 = mu_function(1);
    else
        x1 = mu_function(pp);
        x2 = mu_function(pp+1);
    end
    newPoints = linspace(x1,x2,newRes+1);
    stretchedFun = [stretchedFun, newPoints(1:newRes)];
    
end
clear x1 x2 pp newRes newPoints
originalFun = [mu_function, mu_function, mu_function,mu_function];



% 1. plot both original and stretched functions
figure(4)
plot(originalFun,'Color',rgb(colorSpectrum{mfun}),'LineWidth',3)
hold on
plot(stretchedFun,'Color',rgb(colorSpectrum{mfun}))
legend('original','stretched')
axis([0 720 -0.3 2.1])


% 2. loop through original and stretched functions and simulate
tstep = 5/3600; % time step in hour

tau_compiled = nan(gen,length(T));
Vd_compiled = nan(gen,length(T));
Vb_compiled = nan(gen,length(T));

for ff = 1:2
    
    % 0. determine current growth rate signal
    if ff == 1
        current_function = originalFun;
    else
        current_function = stretchedFun;
    end
    
    
    % 0. start generation counter
    g_counter = 1;
    
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 2; %
    tau = -(17/60)*Vi + (64/60);
    period_frac = 1; % time counter for growth rate
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    
    % 1. simulate!
    counter_loop = 0; % each loop is another timestep
    while g_counter < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % i. store current cell cycle Vi
            Vb_compiled(g_counter,ff) = Vi;
        end
        
        % ii. calculate exponent
        muT = current_function(period_frac)*tstep; % mu(t)*t in hours
        
        % iii. calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % iv. calculate time into cell cycle
        t_cc = t_cc + tstep;
        
        % v. store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = tstep;
        else
            t(counter_loop) = t(counter_loop-1) + tstep;
        end
        
        % vi. determine whether cell divides
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_compiled(g_counter,ff) = Vt;
            tau_compiled(g_counter,ff) = tau;
            
            % re-set parameters for next cell cycle
            g_counter = g_counter + 1;
            Vi = Vt/2;
            tau = -(17/60)*Vi + (64/60);
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
        % vii. determine growth rate value of next timestep
        if period_frac == length(current_function)
            period_frac = 1;
        else
            period_frac = period_frac + 1;
        end
        
    end
    
    % 2. plot volumetric growth over time!
    figure(5)
    subplot(2,1,ff)
    hold on
    plot(t,V,'Color',rgb(color15{ff}))
    axis([0 40 0.5 3])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('T=',num2str(T(ff))))
    
    
    % 3. plot histogram of birth volume
    figure(6)
    hold on
    histogram(Vb_compiled(steadies,ff),'Facecolor',rgb(color15{ff}))
    ylabel('Cell count')
    xlabel('Birth volume')
    
end
clear tscale Vi Vt V t t_cc tau tstep g_counter counter_loop ff muT

figure(6)
legend('15 min','60 min')
title('histogram, rescaled T=15 signal')

% 4. store mu function specific data
if mfun == 1
    ori15_15.Vb = Vb_compiled;
    ori15_15.Vd = Vd_compiled;
    ori15_15.tau = tau_compiled;
else
    ori15_60.Vb = Vb_compiled;
    ori15_60.Vd = Vd_compiled;
    ori15_60.tau = tau_compiled;
end


%% Part 4. simulate using mu(t) from T = 60 at both timescales

% manipulates measured signal from T = 60 min experiments such that it
% shrinks to a 15 min period in simulations.

clear ans color15 is15 mfun oriT originalFun stretchedFun 
clear tau_compiled Vb_compiled Vd_compiled

% 0. intialize timescale from which measured function is drawn
mfun = 2;
oriT = T(mfun);
color60 = {'MediumSlateBlue','Indigo'};

% 0. intialize number of cell generations simulated
gen = 1000;

% 0. initialize mu function for current timescale
mu_function = mufun60;

% 1. prepare scrunched and "originial" function of 60 min timescale
for tscale = 1:length(T)
    
    % 1. define current growth rate function (T=15 or 60) and timescale
    %mu_function = mufun(:,tscale);
    
    % 2. determine current timestep between values
    current_T = T(tscale);
    tstep_sec = (current_T)/length(mu_function); % timestep in sec
    
    % 3. add points between each measured step
    scale = tstep_sec/goalStep;
    scaled_fun = [];
    for pp = 1:length(mu_function)
        
        if pp == length(mu_function)
            x1 = mu_function(pp);
            x2 = mu_function(1);
        else
            x1 = mu_function(pp);
            x2 = mu_function(pp+1);
        end
        newPoints = linspace(x1,x2,scale+1);
        scaled_fun = [scaled_fun, newPoints(1:scale)];
        
    end
    clear x1 x2 pp scale newPoints
    
    scaled_function{tscale} = scaled_fun;
    
end
clear tscale tstep_sec current_T scaled_fun

originalFun = scaled_function{2};
scrunchedFun = [scaled_function{1}, scaled_function{1}, scaled_function{1}, scaled_function{1}];



% 1. plot both original and stretched functions
figure(7)
plot(scrunchedFun,'Color',rgb(colorSpectrum{mfun}))
hold on
plot(originalFun,'Color',rgb(colorSpectrum{mfun}),'LineWidth',3)
legend('scrunched','original')
axis([0 720 -0.3 2.1])

% 2. loop through original and stretched functions and simulate
tstep = 5/3600; % time step in hour

tau_compiled = nan(gen,length(T));
Vd_compiled = nan(gen,length(T));
Vb_compiled = nan(gen,length(T));

for ff = 1:2
    
    % 0. determine current growth rate signal
    if ff == 1
        current_function = scrunchedFun;
    else
        current_function = originalFun;
    end
    
    
    % 0. start generation counter
    g_counter = 1;
    
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 2; %
    tau = -(17/60)*Vi + (64/60);
    period_frac = 1; % time counter for growth rate
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    
    % 1. simulate!
    counter_loop = 0; % each loop is another timestep
    while g_counter < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % i. store current cell cycle Vi
            Vb_compiled(g_counter,ff) = Vi;
        end
        
        % ii. calculate exponent
        muT = current_function(period_frac)*tstep; % mu(t)*t in hours
        
        % iii. calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % iv. calculate time into cell cycle
        t_cc = t_cc + tstep;
        
        % v. store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = tstep;
        else
            t(counter_loop) = t(counter_loop-1) + tstep;
        end
        
        % vi. determine whether cell divides
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_compiled(g_counter,ff) = Vt;
            tau_compiled(g_counter,ff) = tau;
            tauT_compiled(g_counter,ff) = t(counter_loop);
            
            % re-set parameters for next cell cycle
            g_counter = g_counter + 1;
            Vi = Vt/2;
            tau = -(17/60)*Vi + (64/60);
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
        % vii. determine growth rate value of next timestep
        if period_frac == length(current_function)
            period_frac = 1;
        else
            period_frac = period_frac + 1;
        end
        
    end
    
    % 2. plot volumetric growth over time!
    figure(8)
    subplot(2,1,ff)
    hold on
    plot(t,V,'Color',rgb(color60{ff}))
    axis([0 40 0.5 3])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('T=',num2str(T(ff))))
    
    
    % 3. plot histogram of birth volume
    figure(9)
    hold on
    histogram(Vb_compiled(steadies,ff),'Facecolor',rgb(color60{ff}))
    ylabel('Cell count')
    xlabel('Birth volume')
    
end
clear tscale Vi Vt V t t_cc tau tstep g_counter counter_loop

figure(9)
legend('15 min','60 min')
title('histogram, rescaled T=60 signal')

% 4. store mu function specific data
if mfun == 1
    ori60_15.Vb = Vb_compiled;
    ori60_15.Vd = Vd_compiled;
    ori60_15.tau = tau_compiled;
    ori60_15.tauT = tauT_compiled;
else
    ori60_60.Vb = Vb_compiled;
    ori60_60.Vd = Vd_compiled;
    ori60_60.tau = tau_compiled;
    ori60_60.tauT = tauT_compiled;
end


mean(Vb_compiled)
std(Vb_compiled)


%% Part 5. determine periodicity of birth events

% 0. initialize binning parameters
binsPerPeriod_resolved = 12;
binsPerPeriod = [1.5,6]; % 15 and 60 min
timescale_h = [0.25, 1];


% 1. re-define period to begin at start of low nutrient pulse, by
%      subtracting quarter period from corrected timestamp
birthEvent_timestamps = ori60_60.tauT;
birthEvent_taus = ori60_60.tau;
birthEvent_Vb = ori60_60.Vb;
birthEvent_Vd = ori60_60.Vd;
shifted_birthTimestamps = birthEvent_timestamps-(timescale_h/4);


% 2. trim data to remove first 50 generations
shifted_birthTimestamps_trimmed = shifted_birthTimestamps(51:end,:);


% 3. bin birth data by period fraction
timeInPeriods_births = shifted_birthTimestamps_trimmed./timescale_h; % unit = sec/sec
timeInPeriodFraction_births = timeInPeriods_births - floor(timeInPeriods_births);
assignedBin_birthTimestamps = ceil(timeInPeriodFraction_births * binsPerPeriod_resolved);

tau_binnedByPeriodFraction_15 = accumarray(assignedBin_birthTimestamps(:,1), birthEvent_taus(51:end,1), [], @(x) {x});
tau_binnedByPeriodFraction_60 = accumarray(assignedBin_birthTimestamps(:,2), birthEvent_taus(51:end,2), [], @(x) {x});

vb_binnedByPeriodFraction_15 = accumarray(assignedBin_birthTimestamps(:,1), birthEvent_Vb(51:end,1), [], @(x) {x});
vb_binnedByPeriodFraction_60 = accumarray(assignedBin_birthTimestamps(:,2), birthEvent_Vb(51:end,2), [], @(x) {x});

vd_binnedByPeriodFraction_15 = accumarray(assignedBin_birthTimestamps(:,1), birthEvent_Vd(51:end,1), [], @(x) {x});
vd_binnedByPeriodFraction_60 = accumarray(assignedBin_birthTimestamps(:,2), birthEvent_Vd(51:end,2), [], @(x) {x});


% calculate
birthsPerPeriodFraction_15 = cellfun(@length,tau_binnedByPeriodFraction_15);
birthsPerPeriodFraction_60 = cellfun(@length,tau_binnedByPeriodFraction_60);

tausPerPeriodFraction(:,1) = cellfun(@mean,tau_binnedByPeriodFraction_15);
tausPerPeriodFraction(:,2) = cellfun(@mean,tau_binnedByPeriodFraction_60);

VbPerPeriodFraction(:,1) = cellfun(@mean, vb_binnedByPeriodFraction_15);
VbPerPeriodFraction(:,2) = cellfun(@mean, vb_binnedByPeriodFraction_60);

VdPerPeriodFraction(:,1) = cellfun(@mean, vd_binnedByPeriodFraction_15);
VdPerPeriodFraction(:,2) = cellfun(@mean, vd_binnedByPeriodFraction_60);

deltaPerPeriodFraction = VdPerPeriodFraction - VbPerPeriodFraction;



% 4.  convert bin # to absolute time (sec)
timePerBin = timescale_h./binsPerPeriod_resolved;  % in sec
binPeriod = linspace(1, binsPerPeriod_resolved, binsPerPeriod_resolved);
timePeriod(:,1) = timePerBin(1)*binPeriod';
timePeriod(:,2) = timePerBin(2)*binPeriod';


% 5. normalize binned births by total births
birthsPerPeriodFraction_normalized(:,1) = birthsPerPeriodFraction_15./950;
birthsPerPeriodFraction_normalized(:,2) = birthsPerPeriodFraction_60./950;



figure(1)
for cc = 1:2
    subplot(1,2,cc)
    stem(timePeriod(:,cc), VbPerPeriodFraction(:,cc),'Color',[0 0.7 0.7])
    axis([min(timePeriod(:,cc)) max(timePeriod(:,cc)) 0 1])
    title(strcat('fluc: ',num2str(timescale_h(cc))))
    xlabel('time (h)')
    ylabel('mean birth volume for births at period fraction')
end

figure(2)
for cc = 1:2
    subplot(1,2,cc)
    stem(timePeriod(:,cc), tausPerPeriodFraction(:,cc),'Color',[0 0.7 0.7])
    axis([min(timePeriod(:,cc)) max(timePeriod(:,cc)) 0 1])
    title(strcat('fluc: ',num2str(timescale_h(cc))))
    xlabel('time (h)')
    ylabel('mean division time for births at period fraction')
end

figure(3)
for cc = 1:2
    subplot(1,2,cc)
    stem(timePeriod(:,cc), birthsPerPeriodFraction_normalized(:,cc),'Color',[0 0.7 0.7])
    axis([min(timePeriod(:,cc)) max(timePeriod(:,cc)) 0 0.2])
    title(strcat('fluc: ',num2str(timescale_h(cc))))
    xlabel('time (h)')
    ylabel('probability of birth event, normalized by total')
end

figure(4)
for cc = 1:2
    subplot(1,2,cc)
    stem(timePeriod(:,cc), deltaPerPeriodFraction(:,cc),'Color',[0 0.7 0.7])
    axis([min(timePeriod(:,cc)) max(timePeriod(:,cc)) 0 1])
    title(strcat('fluc: ',num2str(timescale_h(cc))))
    xlabel('time (h)')
    ylabel('added v for births at period fraction')
end

grid on

%clear birthEvents birthEvent_timestamps binVector timeVector padding
%clear binnedDrops isDrops stagePositions Time condition xy_start xy_end
