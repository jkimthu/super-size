%% ss13: Vdiv_i vs Vb_i compiled from all replicates


%  Goal: for comparison with single replicate data


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. calculate stats for each replicate
%  Part 3. plot predictions vs measured data

%  Last edit: Jen Nguyen, 2020 May 16
%  Commit: scatter of Vdiv_i vs Vb_i, compiled from all replicate data

%  OK let's go!

%% Part 0. initialize data, saved from figure1A_division.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')
%lamb = 1; % column in compiled_data (meta) that is lambda


% 0. initialize plotting parameters
palette = {'DarkTurquoise','Chocolate'};
environment_order = {'low',30,300,900,3600,'ave','high'};
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


%% Part 2. compile data from each replicate

% 0. initialize parameters for which calculate stats in organized_data
vol_birth = 1;     % volume at birth = col in cc (compiled in figure1A_division.m)
vol_div = 2;     % volume at division = col in cc


% 0. for each condition of interest
conditionsOI = [900; 3600];
counter = 0;
numcells = [];
Vdiv_i = [];
Vb_i = [];
for cond = 1:length(conditionsOI)
    
    % 1. isolate replicate data from each conditions
    currCond = conditionsOI(cond);
    if currCond == 900
        col = 4;
    elseif currCond == 3600
        col = 5;
    end
    currData = organized_data(:,col);
    noData = cellfun(@isempty,currData);
    currData = currData(noData == 0);
    
    
    % 2. compile data from each replicate 
    for rep = 1:length(currData)
        
        % i. tick for figure counter
        counter = counter + 1;
        
        % ii. gather parameters
        repData = currData{rep,1};
        Vb = repData.cc(:,vol_birth);
        Vdiv = repData.cc(:,vol_div);
        
        % iii. keep only data within 95% of mean Vb
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
        Vdiv_trim1 = Vdiv(Vb < Vb_mean + 2*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        Vdiv_trim2 = Vdiv_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        clear Vb_trim1 Vdiv_trim1 sigma
        clear repData Vdiv Vb sigma
        
        % iv. concatenate data
        Vdiv_i = [Vdiv_i; Vdiv_trim2];
        Vb_i = [Vb_i; Vb_trim2];
        
    end
    clear Vb_trim2 Vdiv_trim2 dt_trim2
    
    % iv. plot
    figure(cond)
    scatter(Vb_i,Vdiv_i,'MarkerEdgeColor',rgb(palette{cond}))
    xlabel('Vb_i (cubic um)')
    ylabel('Vdiv_i')
    numcells = length(Vdiv_i);
    title(strcat('compiled-t',num2str(currCond/60),'-n',num2str(numcells)))
    
    % overlay fit
    fit = polyfit(Vb_i,Vdiv_i,1);
    x = linspace(min(Vb_i),max(Vb_i),10);
    y = fit(1).*x + fit(2);
    hold on
    if currCond == 900
        text(min(Vb_i)-0.5,min(Vb_i),strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette{cond}))
    elseif currCond == 3600
        text(max(Vb_i)-1,max(Vb_i)+5,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette{cond}))
    end
    plot(x,y,'Color',rgb(palette{cond}))
        
end
clear currCond currData col rep cond
clear conditionsOI counter
clear Vb_mean x y tau fluc high low ave

