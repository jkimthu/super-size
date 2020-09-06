%% ss21: confirm that alpha and tau_o covary to fit hyperbole


%  Goal: ensure collapse of data is not a fluke


%  Strategy: 

%  edited ss17 to include quantification of average alpha and tau_o per
%  replicate condition (Parts 2 & 3).

%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. compile data from steady replicates
%  Part 3. compile data from fluctuating replicates
%  Part 4. save data!
%  Part 5. calculate replicate means for plotting in part 6
%  Part 6. plot!


%  Last edit: Jen Nguyen, 2020 Sept 6
%  Commit: visualizing covariance of alpha and tau_o


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
fluc_alpha = 1; low = 2; ave = 3; high = 4;          % row number in data structure 

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

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss21')

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
    alpha_steady = nan(4,3);
    taunot_steady = nan(4,3);
    
    
    % 0. initialize rows of steady columns associated with current timescale
    %    in organized data
    rows_steady = Tgroup{ts};
    
    
    % 1. loop through steady conditions (columns in organized data)
    %    and compile replicate data
    for sc = 1:length(conditions_steady)
        
        
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
            
            
            % 5. calculate single-cell slope_i 
            slope = dt./Vb;
            
            
            % 6. remove single-cell data outside of 3 st dev from mean slope_i (also in ss16.m)
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
            
            
            % 7. store single-cell data
            numcells_steady(rep,sc) = length(slope_i);
            compiled_slopes_i{rep,sc} = slope_i;
            compiled_Vb_i{rep,sc} = Vb_i;
            compiled_tau_i{rep,sc} = tau_i;
            compiled_lambda_i{rep,sc} = lambda_i;
            
            
            % 8. use trimmed single-cell data to quantify average data
            
            %   i. scatter tau_i vs Vb_i
            figure(rep)
            scatter(Vb_i,tau_i)
            
            %  ii. use linear fit to find alpha (slope) and tau_o (y-int)
            fit = polyfit(Vb_i,tau_i,1);
            x = linspace(min(Vb_i),max(Vb_i),10);
            y = fit(1).*x + fit(2);
            hold on
            text(3,5,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))))
            plot(x,y)
            ylabel('tau_i')
            xlabel('Vb_i')
            
            %  iii. save plot
%            saveas(gcf,strcat('ss21-taui-v-Vbi-',conditions_steady{sc},'-rep-',num2str(rep)),'epsc')
            close(gcf)
            
            %  iv. save alpha and tau_o
            alpha_steady(rep,sc) = fit(1);
            taunot_steady(rep,sc) = fit(2);
            clear Vb_i tau_i lambda_i slope_i fit x y
            
        end 
                
    end
    clear sc rep
    
    compiled_steady{ts}.slope_i = compiled_slopes_i;
    compiled_steady{ts}.Vb_i = compiled_Vb_i;
    compiled_steady{ts}.tau_i = compiled_tau_i;
    compiled_steady{ts}.lambda_i = compiled_lambda_i;
    compiled_steady{ts}.alpha = alpha_steady;
    compiled_steady{ts}.tau_o = taunot_steady;
    compiled_steady{ts}.numcells = numcells_steady;
    
end
clear numcells_steady lamb tau vol_birth exp ts currCond
clear compiled_slopes_i compiled_Vb_i compiled_tau_i compiled_lambda_i
clear taunot_steady alpha_steady currData rows_steady


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
fluc_alpha = nan(4,4);
fluc_taunot = nan(4,4);
    

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
        %clear Vb_i tau_i lambda_i slope_i
        
        
            
        
        % 6. use trimmed single-cell data to quantify average data
        
        %   i. scatter tau_i vs Vb_i
        figure(rr)
        scatter(Vb_i,tau_i)
        
        %  ii. use linear fit to find alpha (slope) and tau_o (y-int)
        fit = polyfit(Vb_i,tau_i,1);
        x = linspace(min(Vb_i),max(Vb_i),10);
        y = fit(1).*x + fit(2);
        hold on
        text(3,5,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))))
        plot(x,y)
        ylabel('tau_i')
        xlabel('Vb_i')
        
        %  iii. save plot
        %saveas(gcf,strcat('ss21-taui-v-Vbi-tscale-',num2str(ts),'-rep-',num2str(rr)),'epsc')
        close(gcf)
        
        %  iv. save alpha and tau_o
        fluc_alpha(rr,ts) = fit(1);
        fluc_taunot(rr,ts) = fit(2);
        clear Vb_i tau_i lambda_i slope_i fit x y
        
      
    end
    clear rr
    
end
    
clear numcells lamb tau vol_birth exp ts currCond col


%% Part 4. save data!

cd('/Users/jen/super-size/')
save('ss21.mat', 'compiled_steady','numcells_fluc','fluc_slopes_i','fluc_Vb_i','fluc_tau_i','fluc_lambda_i','fluc_alpha','fluc_taunot')


%% Part 5. plot hyperbolas (measured tau_o and alpha pairs)

%  Strategy:
%
%  0. load data and initialize colors

%     calculate tau (y) across a range of Vb (x) according to:
%
%           tau = tau_o + alpha * Vb
%
%     where alpha = slope of line fit to scatter (tau_i vs. Vb_i)
%           tau_o = y-intercept of line fit to scatter
%

%  1. plot tau vs Vb for measured alpha and tau_o  (should fall on hyperbola)
%  2. plot tau vs Vb for fixed alpha and all tau_o  (should fall off)
%  3. plot measured alpha vs tau_o


% SECTON 1

clc
clear

cd('/Users/jen/super-size/')
load('ss21.mat')
clear fluc_lambda_i fluc_slopes_i fluc_tau_i fluc_Vb_i


% 0. initialize colors and shape for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};


% 1. calculate tau (y) across a range of Vb (x) according to:
%
%           tau = tau_o + alpha * Vb

% 1A. with alpha and tau_o pairs measured from STEADY DATA
% 1A. i. concatenate alpha and tau_o from steady conditions of diff timescales into one matrix 
alpha_steady = [];
taunot_steady = [];
for tscale = 1:4
    
    curr_alpha = compiled_steady{1,tscale}.alpha;
    curr_taunot = compiled_steady{1,tscale}.tau_o;
    
    alpha_steady = [alpha_steady; curr_alpha];
    taunot_steady = [taunot_steady; curr_taunot];
    
end
clear curr_alpha curr_taunot tscale


% 1A. ii. for each measured alpha and tau_o pair, vary Vb to calculate tau
%      a) vary Vb
Vb_vector = 0.5:0.5:10;

%      b) prepare to loop through alpha and tau_o pairs
dim = size(alpha_steady);
numpairs = dim(1) * dim(2); 

tau_steady_calculated = cell(dim(1),dim(2));
slope_steady_calculated = cell(dim(1),dim(2));
for pair = 1:numpairs
    
    %  c) determine whether to skip index (no pair data)
    tau_o = taunot_steady(pair);
    alpha = alpha_steady(pair);
    if isnan(tau_o) == 1
        continue
    end
    
    %  d) for indeces with data, calculate tau according to:
    %     tau = tau_o + alpha * Vb
    tau = tau_o + alpha * Vb_vector';
    tau_steady_calculated{pair} = tau;
    
    %  e) plot tau/Vb vs Vb, where tau/Vb = alpha + tau_o/Vb
    slope = tau./Vb_vector';
    slope_steady_calculated{pair} = slope;
    
    if pair <= dim(1)*1
        color = palette_steady{1};
    elseif pair <= dim(1)*2
        color = palette_steady{2};
    else
        color = palette_steady{3};
    end
    
    figure(1)
    plot(Vb_vector,slope,'Color',rgb(color),'LineWidth',2)
    hold on
    
end
clear dim pair slope numpairs tau tau_o alpha color

figure(1)
xlabel('simulated Vb')
ylabel('simulated slope (tau/Vb)')
title('simulated slope vs Vb from measured alpha and tau_o')


% 1B. with alpha and tau_o pairs measured from FLUCTUATING DATA
%      a) prepare to loop through alpha and tau_o pairs
dimf = size(fluc_alpha);
numpairf = dimf(1) * dimf(2); 

tau_fluc_calculated = cell(dimf(1),dimf(2));
slope_fluc_calculated = cell(dimf(1),dimf(2));
for pf = 1:numpairf
    
    %  b) determine whether to skip index (no pair data)
    tau_o = fluc_taunot(pf);
    alpha = fluc_alpha(pf);
    if isnan(tau_o) == 1
        continue
    end
    
    %  d) for indeces with data, calculate tau according to:
    %     tau = tau_o + alpha * Vb
    tau = tau_o + alpha * Vb_vector';
    tau_fluc_calculated{pf} = tau;
    
    %  e) plot tau/Vb vs Vb, where tau/Vb = alpha + tau_o/Vb
    slope = tau./Vb_vector';
    slope_fluc_calculated{pf} = slope;
    
    if pf <= dimf(1)*1
        color = palette_fluc{1};
    elseif pf <= dimf(1)*2
        color = palette_fluc{2};
    elseif pf <= dimf(1)*3
        color = palette_fluc{3};
    else
        color = palette_fluc{4};
    end
    
    figure(1)
    plot(Vb_vector,slope,'Color',rgb(color),'LineWidth',2)
    hold on
    
end
clear dimf pf slope numpairf tau tau_o alpha color


%% Part 6. plot hyperbolas(?) for fixed alpha and measured tau_o

%  Strategy:
%
%  0. load data and initialize colors

%     calculate tau (y) across a range of Vb (x) according to:
%
%           tau = tau_o + alpha * Vb
%
%     where alpha = slope of line fit to scatter (tau_i vs. Vb_i)
%           tau_o = y-intercept of line fit to scatter
%

%  1. plot tau vs Vb for measured alpha and tau_o  (should fall on hyperbola)
%  2. plot tau vs Vb for fixed alpha and all tau_o  (should fall off)
%  3. plot measured alpha vs tau_o


% SECTON 2

% 2. calculate tau (y) across a range of Vb (x) according to:
%
%           tau = tau_o + alpha * Vb


% 2. with fixed alpha and various tau_o measured from STEADY DATA
% 2. i. prepare to loop through alphas and generate tau_o vector
%dim = size(alpha_steady);
%num_alphas = dim(1) * dim(2); 
alphas = alpha_steady(~isnan(alpha_steady));
taunot_vector = taunot_steady(~isnan(taunot_steady));


% 2. ii. for each fixed alpha/varied tau-o "pair", vary Vb to calculate tau
%      a) vary Vb
Vb_vector = 0.5:0.5:10;

%tau_steady_broken = cell(length(alphas),1);
%slope_steady_broken = cell(length(alphas),1);
for aa = 1:length(alphas)
    
    %  b) determine current alpha
    alpha = alphas(aa);
    
    
    %  c) determine current tau_o
    for tn = 1:length(taunot_vector)
        
        tau_o = taunot_vector(tn);
        
        %  d) calculate tau according to: tau = tau_o + alpha * Vb
        tau = tau_o + alpha * Vb_vector';
        %tau_steady_broken{aa} = tau;
        
        %  e) plot tau/Vb vs Vb, where tau/Vb = alpha + tau_o/Vb
        slope = tau./Vb_vector';
        %slope_steady_broken{aa} = slope;
        
        if (aa==tn) == 1
            if aa <= 12
                color = palette_steady{1};
            elseif aa <= 25
                color = palette_steady{2};
            else
                color = palette_steady{3};
            end
            lw = 2;
            
            figure(50)
            plot(alphas(aa),taunot_vector(tn),'o','Color',rgb(color),'MarkerFaceColor',rgb(color))
            hold on
        else
            color = 'Silver';
            lw = 1;
        end
        
        %figure(aa)
        %plot(Vb_vector,slope,'Color',rgb(color),'LineWidth',lw)
        %hold on
        
    end
    %figure(aa)
    %title('simulated curves with fixed alpha, varied tau_o')
    %ylabel('simulated slope (tau/Vb)')
    %xlabel('simulated Vb')
    %axis([0 10 -20 200])
    
  
end
clear aa iii bp slope numpairs tau tau_o alpha color lw tn

figure(50)
title('measured tau_o vs alpha pairs')
ylabel('measured tau_o')
xlabel('measured alpha')


%% Part 7. plot measured alpha vs tau_o

%  Strategy:
%
%  0. load data and initialize colors

%     calculate tau (y) across a range of Vb (x) according to:
%
%           tau = tau_o + alpha * Vb
%
%     where alpha = slope of line fit to scatter (tau_i vs. Vb_i)
%           tau_o = y-intercept of line fit to scatter
%

%  1. plot tau vs Vb for measured alpha and tau_o  (should fall on hyperbola)
%  2. plot tau vs Vb for fixed alpha and all tau_o  (should fall off)
%  3. plot measured alpha vs tau_o



% SECTION 3

% 3A. ii. plot alpha vs tau_o (measured pairs)
for iii = 1:length(alphas)
    if iii <= 12
        color = palette_steady{1};
    elseif iii <= 25
        color = palette_steady{2};
    else
        color = palette_steady{3};
    end
    figure(1)
    plot(alphas(iii),taunot_vector(iii),'o','Color',rgb(color),'MarkerFaceColor',rgb(color))
    hold on
end
clear iii color
figure(1)
title('measured tau_o vs alpha pairs')
ylabel('tau_o')
xlabel('alpha')


% 3B. add alpha and tau_o pairs measured from FLUCTUATING DATA
alphaf = fluc_alpha(~isnan(fluc_alpha));
taunot_f = fluc_taunot(~isnan(fluc_taunot));

% 3B. ii. plot alpha vs tau_o (measured pairs)
for iif = 1:length(alphaf)
    if iif <= 3
        color = palette_fluc{1};
    elseif iif <= 6
        color = palette_fluc{2};
    elseif iif <= 10
        color = palette_fluc{3};
    else
        color = palette_fluc{4};
    end
    figure(1)
    plot(alphaf(iif),taunot_f(iif),'o','Color',rgb(color),'MarkerFaceColor',rgb(color))
    hold on
end
clear iii color




