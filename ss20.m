%% ss20: single cell slope_i over time in fluctuating nutrient


%  Goal: determine whether average size difference between t15 and t60 can
%        be observed from the slope_i transitions in both

%        plot single-cell slope_i over time under fluctuating environments
%         1. t15
%         2. t60


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keeping replicates apart
%  Part 2. gather single cell data including slopes_i for each experimental replicate
%  Part 3. plot slopes_i vs time for each experiment


%  Last edit: Jen Nguyen, 2020 July 30
%  Commit: first commit, plot tau_i/Vb_i over time for dif timescales


%  OK let's go!


%% Part 0. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')


% 0. initialize plotting parameters
environment_order = {'low',30,300,900,3600,'ave','high'};


% 0. initialize colors for plotting
palette_steady = {'Indigo','GoldenRod','FireBrick'};
palette_fluc = {'RoyalBlue','CornflowerBlue','DeepSkyBlue','CadetBlue'};


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
                envr = 4;
            elseif ismember(exp,t5) == 1
                envr = 5;
            elseif ismember(exp,t15) == 1
                envr = 6;
            elseif ismember(exp,t60) == 1
                envr = 7;
            end
           
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
           
            
        else % steady
            
            if exp_cond == low % low
                envr = 1;
            elseif exp_cond == ave
                envr = 2;
            elseif exp_cond == high
                envr = 3;
            end
            
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
            
        end
        
    end
    
end
clear envr exp_cond data counter
clear t0_5 t5 t15 t60
clear compiled_data low ave high fluc exp


%% Part 2. gather single cell data including slopes_i for each experimental replicate

%          calculate single-cell slopes (tau_i divided by Vb_i) from all

% 0. initialize step size for binning time
binsPerHour = 12; % 5 min bins


% 0. initialize variables for storing data
numexpt = length(exptArray);
slope_i = cell(numexpt,3);
Vb_i = cell(numexpt,3);
tau_i = cell(numexpt,3);
tbirth_i = cell(numexpt,3);
tdiv_i = cell(numexpt,3);


% 1. collect data from each steady condition
for sc = 1:3 % sc = steady condition
    
    % i. isolate data from current steady condition
    palette = palette_steady(1,sc);
    current_condData = organized_data(:,sc);
    
    % ii. determine replicates to loop through
    numreps = length(current_condData);
    for rr = 1:numreps
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            % iii. isolate replicate data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_tbirth = rep_meta(:,4);     % each row = cc; number = time at birth (h)
            rep_tdiv = rep_meta(:,5);       % each row = cc; number = time at division (h)
            
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % iv. calculate single-cell slopes
            rep_slopes = rep_taus./rep_birthvols;
            
            
            % v. remove items over 3 std from mean
            mm = mean(rep_slopes);
            thresh = 3 * std(rep_slopes);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = rep_slopes > upper;
            cut_down = rep_slopes < lower;
            toCut = cut_up + cut_down;
            
            slopes_final = rep_slopes(toCut == 0);
            taus_final = rep_taus(toCut == 0);
            Vb_final = rep_birthvols(toCut == 0);
            tbirth_final = rep_tbirth(toCut == 0);
            tdiv_final = rep_tdiv(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear rep_taus rep_birthvols rep_slopes rep_tdiv rep_tbirth
            
            
            % vi. plot distribution of single-cell slopes replicate
            figure(rr)
            histogram(slopes_final,100,'FaceColor',rgb(palette));
            xlabel('single-cell slope')
            ylabel('# cells')
            title(strcat('replicate',num2str(rr)))
            hold on
            
            
            % vii. store all means and standard deviations
            slope_i{rr,sc} = slopes_final;
            Vb_i{rr,sc} = Vb_final;
            tau_i{rr,sc} = taus_final;
            tbirth_i{rr,sc} = tbirth_final;
            tdiv_i{rr,sc} = tdiv_final;
            clear mf sf
            
        end
                
    end

end
clear sc rr slopes_final taus_final Vb_final tbirth_final col


% 2. collect data from each shifting condition
ts_counter = 0;
for ts = 1:4 % ts = fluctuating timescale
    
    % i. isolate data from current steady condition
    if ts == 1
        col = 4;
    elseif ts == 2
        col = 5;
    elseif ts == 3
        col = 6;
    else
        col = 7;
    end
    
    palette = palette_fluc(1,ts);
    current_condData = organized_data(:,col);
    
    
    % ii. determine replicates to loop through
    numreps = length(current_condData);
    for rr = 1:numreps
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            
            ts_counter = ts_counter + 1;
            
            % iii. isolate replicate data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_tbirth = rep_meta(:,4);     % each row = cc; number = time at birth (h)
            rep_tdiv = rep_meta(:,5);       % each row = cc; number = time at division (h)
            
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % iv. calculate single-cell slopes
            rep_slopes = rep_taus./rep_birthvols;
            
            
            % v. remove items over 3 std from mean
            mm = mean(rep_slopes);
            thresh = 3 * std(rep_slopes);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = rep_slopes > upper;
            cut_down = rep_slopes < lower;
            toCut = cut_up + cut_down;
            
            slopes_final = rep_slopes(toCut == 0);
            taus_final = rep_taus(toCut == 0);
            Vb_final = rep_birthvols(toCut == 0);
            tbirth_final = rep_tbirth(toCut == 0);
            tdiv_final = rep_tdiv(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear rep_taus rep_birthvols rep_slopes rep_tdiv rep_tbirth
            
            
            % vi. plot distribution of single-cell slopes replicate
            figure(ts_counter)
            histogram(slopes_final,100,'FaceColor',rgb(palette));
            xlabel('single-cell slope')
            ylabel('# cells')
            title(strcat('replicate',num2str(rr)))
            hold on
            
            
            % vii. store all means and standard deviations
            slope_i{ts_counter,col} = slopes_final;
            Vb_i{ts_counter,col} = Vb_final;
            tau_i{ts_counter,col} = taus_final;
            tbirth_i{ts_counter,col} = tbirth_final;
            tdiv_i{ts_counter,col} = tdiv_final;
            clear mf sf
            
        end
                
    end

end
clear shf rr slopes_final taus_final Vb_final tbirth_final col
clear tdiv_final shf_counter numexpts numreps palette

%save('ss20_output.mat','slope_i','tau_i','Vb_i','tbirth_i','tdiv_i')


%% Part 3. plot slopes_i vs time for each experiment

%          go row by row (experiment)
%          plotting all columns (conditions) in same plot
%          highlight data post shift in filled circles,


cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss20/')


% % 0. initialize colors for plotting shifted condition
% pre_color = rgb('CornflowerBlue');
% post_color = rgb('DarkTurquoise'); % different color for data after shift

% loop through experiments
for ee = 1:numexpt
    
    % 1. isolate current experiment data
    e_slopes = slope_i(ee,:);
    e_tbirths = tbirth_i(ee,:);
    e_tdivs = tdiv_i(ee,:);
    
  
    
    % 3. loop through conditions to plot all single points
    for cond = 1:length(e_slopes)
        
        c_slopes = e_slopes{1,cond};
        if isempty(c_slopes) == 1
            continue
        else
            
            % i. isolate condition time data
            c_tbirths = e_tbirths{1,cond};
            c_tdivs = e_tdivs{1,cond};
            
            
            % ii. determine plotting color by condition
            if cond < 4
                color = rgb(palette_steady{cond});
                
                % iii. plot slope_i vs tbirth_i
                figure(ee)
                scatter(c_tbirths,c_slopes,10,'MarkerEdgeColor',color)
                hold on

            else
                color = rgb(palette_fluc{cond-3});
                
                % iii. plot slope_i vs tbirth_i
                figure(ee)
                scatter(c_tbirths,c_slopes,20,'MarkerFaceColor',color,'MarkerEdgeColor',color)
                hold on
                ylabel('slope_i')
                xlabel('time of birth')
                
                
                saveas(gcf,strcat('ss20-slopei-v-time-',num2str(ee)),'epsc')
                close(gcf)

            end
        end
        clear pre_tbirths pre_tdivs pre_slopes post_tbirths post_tdivs post_slopes
        
    end
    clear e_slopes e_tbirths e_tdivs index tshift
end
clear ee



