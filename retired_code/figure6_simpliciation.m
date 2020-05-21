%%  figure 6 - simplification of measured mu function


%  Goal: what feature of the growth rate signals produces the difference in
%        cell size between T = 15 and 60 min?


%  Strategy:
%
%        Take raw mu function and slowly simplify various elements until
%        over simplification of a specific feature loses size shift


%  Last edit: jen, 2020 Feb 14
%  commit: first commit


% Okie, go go let's go!

%% Part 0. initialize analysis

clear
clc

% 0. initialize stored mu functions
load('mu_functions')     % generated and saved from figure6_simulations
colorSpectrum = {'DodgerBlue','Indigo'};


% 1. from replicates, calculate mean mu(t) for each timescale
mufun15 = mean(signals_15,2);
mufun60 = mean(signals_60,2);


% 2. plot for visual confirmation
figure(1) % mean of replicate means
plot(1:length(mufun15),mufun15,'Color',rgb(colorSpectrum(1)),'LineWidth',1,'Marker','.')
hold on
plot(1:length(mufun60),mufun60,'Color',rgb(colorSpectrum(2)),'LineWidth',1,'Marker','.')
ylabel('Mu (1/h)')
xlabel('Period fraction')
legend('T = 15','T = 60')
title('growth rate function over time')

clear signal_e exp_type

%% Part 1. simulate using T=15 mu(t) at both timescales


% 0. start generation counter
gen = 1000;
g_counter = 1;

% 0. initialize vectors for data storage
tau_15 = nan(gen,1);
Vd_15 = nan(gen,1);
Vi_15 = nan(gen,1);

% 0. timestep for each timescale
t15 = (15/60)/binsPerPeriod; % timestep in hour

% 0. initiation size, timestep, tau
new_cell = 1; % 0 = division did not just happen
Vi = 2; %
tau = -(17/60)*Vi + (64/60);
period_frac = 1; % time counter for growth rate
t_cc = 0; % time counter for cell cycle
V = [];
t = [];

counter_loop = 0; % each loop is another timestep
while g_counter < gen+1
    
    counter_loop = counter_loop + 1;
    
    if new_cell == 1
        % store current cell cycle Vi
        Vi_15(g_counter) = Vi;
    end
    
    % calculate exponent
    muT = mufun15(period_frac)*t15; % mu(t)*t in hours
    
    % calculate new volume after timestep
    if muT < 0
        Vt = Vi;
    else
        Vt = Vi * 2^(muT);
    end
    
    % calculate time into cell cycle
    t_cc = t_cc + t15;
    
    % store timestep data
    V(counter_loop) = Vt;
    if counter_loop == 1
        t(counter_loop) = t15;
    else
        t(counter_loop) = t(counter_loop-1) + t15;
    end
    
    % determine whether cell divides
    
    if t_cc >= tau 
        
        % division! store V_div and tau data
        Vd_15(g_counter) = Vt;
        tau_15(g_counter) = tau;
        
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
    
    if period_frac == 30
        period_frac = 1;
    else
        period_frac = period_frac + 1;
    end
    
end

figure(2)
plot(t,V)
axis([0 40 0.5 3])
xlabel('Time (h)')
ylabel('Volume')
title('visualization 15 min')


%% Part 2. simulate using T=60 mu(t) at both timescales



%% D. simulate cells in 15 min fluctuations, visualization 
%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations



%% E. simulate cells in 60 min fluctuations, visualization
%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations


% 0. start generation counter
gen = 1000;
g60 = 1;

% 0. initialize vectors for data storage
tau_60 = nan(gen,1);
Vd_60 = nan(gen,1);
Vi_60 = nan(gen,1);

% 0. timestep for each timescale
t60 = 1/binsPerPeriod;


% 0. initiation size, timestep, tau
new_cell = 1; % 0 = division did not just happen
Vi = 2; %
tau = -(17/60)*Vi + (64/60);
period_frac = 1; % time counter for growth rate
t_cc = 0; % time counter for cell cycle
V = [];
t = [];

counter_loop = 0; % each loop is another timestep
while g60 < gen+1
    
    counter_loop = counter_loop + 1;
    
    if new_cell == 1
        % store current cell cycle Vi
        Vi_60(g60) = Vi;
    end
    
    % calculate exponent
    muT = mufun60(period_frac)*t60; % mu(t)*t in hours
    
    % calculate new volume after timestep
    if muT < 0
        Vt = Vi;
    else
        Vt = Vi * 2^(muT);
    end
    
    % calculate time into cell cycle
    t_cc = t_cc + t60;
    
    % store timestep data
    V(counter_loop) = Vt;
    if counter_loop == 1
        t(counter_loop) = t60;
    else
        t(counter_loop) = t(counter_loop-1) + t60;
    end
    
    % determine whether cell divides
    
    if t_cc >= tau 
        
        % division! store V_div and tau data
        Vd_60(g60) = Vt;
        tau_60(g60) = tau;
        
        % re-set parameters for next cell cycle
        g60 = g60 + 1;
        Vi = Vt/2;
        tau = -(17/60)*Vi + (64/60);
        new_cell = 1;
        t_cc = 0;
        
    else
        Vi = Vt;
        new_cell = 0;
    end
    
    if period_frac == 30
        period_frac = 1;
    else
        period_frac = period_frac + 1;
    end
    
end

% figure(2)
% plot(t,V)
% %axis([0 40 0.5 3])
% axis([20 40 0.5 2])
% xlabel('Time (h)')
% ylabel('Volume')
% title('visualization 60 min')

clear ans counter_loop experimentFolder expType g15 g60
clear minTime new_cell period_frac specificGrowthRate
clear index growthRates_final muT t t_cc tau
clear timeInPeriodFraction V Vi Vt

%% F. simulation analysis: plot mean Vi, added V vs Vi, tau vs Vi

steadies = 50:1000;

figure(3)
histogram(Vi_15(steadies),'Facecolor',rgb(colorSpectrum{2}))
hold on
histogram(Vi_60(steadies),'Facecolor',rgb(colorSpectrum{1}))
ylabel('Cell count')
xlabel('Initial volume')
legend('15 min','60 min')


% added V vs initial volume
delta_15 = Vd_15 - Vi_15;
delta_60 = Vd_60 - Vi_60;

figure(4)
plot(Vi_15(steadies),delta_15(steadies),'o','Color',rgb(colorSpectrum{2}))
hold on
plot(Vi_60(steadies),delta_60(steadies),'o','Color',rgb(colorSpectrum{1}))
ylabel('Added volume')
xlabel('Initial volume')
legend('15 min','60 min')


% added V binned by birth size
binsize = 0.02; % cubic microns
birth_bin15 = ceil(Vi_15/binsize);
birth_bin60 = ceil(Vi_60/binsize);

% sort added size by birth bin
adder15_binned = accumarray(birth_bin15(steadies),delta_15(steadies),[],@(x) {x});
adder15_mean = cellfun(@mean, adder15_binned);
adder15_std = cellfun(@std, adder15_binned);
adder15_count = cellfun(@length, adder15_binned);
adder15_sem = adder15_std./sqrt(adder15_count);

adder60_binned = accumarray(birth_bin60(steadies),delta_60(steadies),[],@(x) {x});
adder60_mean = cellfun(@mean, adder60_binned);
adder60_std = cellfun(@std, adder60_binned);
adder60_count = cellfun(@length, adder60_binned);
adder60_sem = adder60_std./sqrt(adder60_count);

figure(5)
errorbar(adder15_mean,adder15_sem,'o','Color',rgb(colorSpectrum{2}))
hold on
errorbar(adder60_mean,adder60_sem,'o','Color',rgb(colorSpectrum{1}))
ylabel('Added volume')
xlabel('Initial volume * 0.02 (bin size)')
legend('15 min','60 min')
axis([34 52 .6 1])

% tau vs birth size
figure(6)
plot(Vi_15(steadies),tau_15(steadies)*60,'o','Color',rgb(colorSpectrum{2}))
hold on
plot(Vi_60(steadies),tau_60(steadies)*60,'o','Color',rgb(colorSpectrum{1}))
ylabel('Interdivision time (min)')
xlabel('Initial volume')
legend('15 min','60 min')

%% G. simulate cells in 15 min fluctuations, with no tau function, visualization
%     start with same Vi
%     use fixed tau
%     use different mu(t)
%     simulate 1000 generations

for tau = .5:0.1:1
    
    % 0. start generation counter
    gen = 1000;
    g_counter = 1;
    
    % 0. initialize vectors for data storage
    tau_15 = nan(gen,1);
    Vd_15 = nan(gen,1);
    Vi_15 = nan(gen,1);
    
    % 0. timestep for each timescale
    t15 = (15/60)/binsPerPeriod; % timestep in hour
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 2; %
    period_frac = 1; % time counter for growth rate
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    counter_loop = 0; % each loop is another timestep
    while g_counter < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % store current cell cycle Vi
            Vi_15(g_counter) = Vi;
        end
        
        % calculate exponent
        muT = mufun15(period_frac)*t15; % mu(t)*t in hours
        
        % calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % calculate time into cell cycle
        t_cc = t_cc + t15;
        
        % store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = t15;
        else
            t(counter_loop) = t(counter_loop-1) + t15;
        end
        
        % determine whether cell divides
        
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_15(g_counter) = Vt;
            tau_15(g_counter) = tau;
            
            % re-set parameters for next cell cycle
            g_counter = g_counter + 1;
            Vi = Vt/2;
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
        if period_frac == 30
            period_frac = 1;
        else
            period_frac = period_frac + 1;
        end
        
    end
    
    figure(6)
    hold on
    plot(t,V)
    axis([0 100 -5 100])
    xlabel('Time (h)')
    ylabel('Volume')
    title('visualization 15 min, fixed tau')
    
    clear ans counter_loop experimentFolder expType g15 g60
    clear minTime new_cell period_frac specificGrowthRate
    clear index growthRates_final muT t t_cc tau
    clear timeInPeriodFraction V Vi Vt
    
end
legend('0.5','0.6','0.7','0.8','0.9','1')

%% H. simulate cells in 60 min fluctuations, with no tau function, visualization
%     start with same Vi
%     use fixed tau
%     use different mu(t)
%     simulate 1000 generations

for tau = .5:0.1:1
    % 0. start generation counter
    gen = 1000;
    g60 = 1;
    
    % 0. initialize vectors for data storage
    tau_60 = nan(gen,1);
    Vd_60 = nan(gen,1);
    Vi_60 = nan(gen,1);
    
    % 0. timestep for each timescale
    t60 = 1/binsPerPeriod;
    
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 2; %
    period_frac = 1; % time counter for growth rate
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    counter_loop = 0; % each loop is another timestep
    while g60 < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % store current cell cycle Vi
            Vi_60(g60) = Vi;
        end
        
        % calculate exponent
        muT = mufun60(period_frac)*t60; % mu(t)*t in hours
        
        % calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % calculate time into cell cycle
        t_cc = t_cc + t60;
        
        % store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = t60;
        else
            t(counter_loop) = t(counter_loop-1) + t60;
        end
        
        % determine whether cell divides
        
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_60(g60) = Vt;
            tau_60(g60) = tau;
            
            % re-set parameters for next cell cycle
            g60 = g60 + 1;
            Vi = Vt/2;
            %tau = -(17/60)*Vi + (64/60);
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
        if period_frac == 30
            period_frac = 1;
        else
            period_frac = period_frac + 1;
        end
        
    end
    
    figure(6)
    hold on
    plot(t,V)
    axis([0 100 -5 100])
    xlabel('Time (h)')
    ylabel('Volume')
    title('visualization 60 min, fixed tau')
    
    clear ans counter_loop experimentFolder expType g15 g60
    clear minTime new_cell period_frac specificGrowthRate
    clear index growthRates_final muT t t_cc tau
    clear timeInPeriodFraction V Vi Vt
    
end
legend('0.5','0.6','0.7','0.8','0.9','1')

%% I. simulate cells with steady G_low, G_ave and G_high

% uses same tau function as T = 15 and 60 min conditions
% manually vary Vi to demonstrate that same steady state is reached
% regardless of initial condition.

clear
clc

% 0. initialize steady-state growth rates from steady environments
G_steady = [1.07, 2.31, 2.86]; % low, ave, high


% 0. initialize simulation parameters
gen = 1000;
tstep = 15/3600; % timestep for each timescale

% 0. initialize vectors for data storage
tau_G = nan(gen,length(G_steady));
Vd_G = nan(gen,length(G_steady));
Vb_G = nan(gen,length(G_steady));


for gg = 1:length(G_steady)
    
    % 0. initialize steady-state G for current condition
    currentG = G_steady(gg);
    
    % 0. start generation counter
    cycle = 1;
    
    % 0. initiation size, timestep, tau
    new_cell = 1; % 0 = division did not just happen
    Vi = 5;
    tau = -(17/60)*Vi + (64/60);
    t_cc = 0; % time counter for cell cycle
    V = [];
    t = [];
    
    counter_loop = 0; % each loop is another timestep
    while cycle < gen+1
        
        counter_loop = counter_loop + 1;
        
        if new_cell == 1
            % store current cell cycle Vi
            Vb_G(cycle,gg) = Vi;
        end
        
        % calculate exponent
        muT = currentG * tstep; % mu(t)*t in hours
        
        % calculate new volume after timestep
        if muT < 0
            Vt = Vi;
        else
            Vt = Vi * 2^(muT);
        end
        
        % calculate time into cell cycle
        t_cc = t_cc + tstep;
        
        % store timestep data
        V(counter_loop) = Vt;
        if counter_loop == 1
            t(counter_loop) = tstep;
        else
            t(counter_loop) = t(counter_loop-1) + tstep;
        end
        
        % determine whether cell divides
        
        if t_cc >= tau
            
            % division! store V_div and tau data
            Vd_G(cycle,gg) = Vt;
            tau_G(cycle,gg) = tau;
            
            % re-set parameters for next cell cycle
            cycle = cycle + 1;
            Vi = Vt/2;
            tau = -(17/60)*Vi + (64/60);
            new_cell = 1;
            t_cc = 0;
            
        else
            Vi = Vt;
            new_cell = 0;
        end
        
    end
    
    figure(6)
    hold on
    plot(t,V)
    axis([0 30 0 6])
    xlabel('Time (h)')
    ylabel('Volume')
    title('Vi = 2')
    
    clear ans counter_loop experimentFolder expType
    clear minTime new_cell period_frac specificGrowthRate
    clear index growthRates_final muT t t_cc tau
    clear timeInPeriodFraction V Vi Vt
    
end
legend('Glow','Gave','Ghigh')

