%% ss7: scatter of tau_i and birth size_i or 1/lambda_i per replicate


%  Goal: do these distributions appear like those in our simulations?


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. calculate stats for each replicate


%  Last edit: Jen Nguyen, 2020 Mar 12
%  Commit: scatter of single cell tau vs Vb or 1/lambda for each replicate


%  OK let's go!

%% Part 0. initialize data, saved from figure1A_division.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')
lamb = 1; % column in compiled_data (meta) that is lambda


% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
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


%% Part 2. plot scatter for each replicate

% 0. initialize parameters for which calculate stats in organized_data
metric = 1; % volume at birth = col in in cc
lamb = 1;   % mean growth rate = col 1 in meta
tau = 2;    % interdivision time = col 2 in meta


% 0. for each condition of interest
conditionsOI = [900; 3600];
counter = 0;
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
    
    
    % 2. for each replicate, plot data
    for rep = 1:length(currData)
        
        % i. tick for figure counter
        counter = counter + 1;
        
        % ii. gather parameters
        repData = currData{rep,1};
        gr = repData.meta(:,lamb);
        gr_inverse = 1./gr;
        Vb = repData.cc(:,metric);
        dt = repData.meta(:,tau);
        
        % iii. keep only data within 95% of mean Vb
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
        gr_inverse_trim1 = gr_inverse(Vb < Vb_mean + 2*sigma);
        dt_trim1 = dt(Vb < Vb_mean + 2*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        gr_inverse_trim2 = gr_inverse_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        dt_trim2 = dt_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        clear dt_trim1 gr_inverse_trim1 Vb_trim1 sigma
        clear repData gr gr_inverse Vb dt Vb_mean sigma
        
        % iv. plot
        figure(counter)
        subplot(1,2,1)
        scatter(gr_inverse_trim2*60,dt_trim2,'MarkerEdgeColor',rgb(palette{col}))
        axis([0 100 0 100])
        ylabel('interdivision time (min)')
        xlabel('inverse lambda (min)')
        title(strcat('T',num2str(currCond),': rep-',num2str(rep)))
        
        subplot(1,2,2)
        scatter(Vb_trim2,dt_trim2,'MarkerEdgeColor',rgb(palette{col}))
        axis([1 5 0 100])
        ylabel('interdivision time (min)')
        xlabel('birth volume')
        title(strcat('n=',num2str(length(dt_trim2))))
        
        % overlay fit
        fit = polyfit(Vb_trim2,dt_trim2,1);
        x = linspace(1,5,10);
        y = fit(1).*x + fit(2);
        subplot(1,2,2)
        hold on
        plot(x,y,'Color',rgb('SlateGray'))
        text(1.5,10,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))))
        
    end
end
clear currCond currData col clear lamb metric rep cond
clear conditionsOI Vb_trim2 dt_trim2 gr_inverse_trim2 counter

 
%% Part 2b. Alternative means of calculating lambda

% a suggested by Jonasz
% Lambda_i=(1\tau_i) log(V_division_i / V_birth_i).    
    
% 0. initialize parameters for which calculate stats in organized_data
birthV = 1; % volume at birth = col in in cc
divV = 2;   % volume at division
tau = 2;    % interdivision time = col 2 in meta


% 0. for each condition of interest
conditionsOI = [900; 3600];
counter = 0;
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
    
    
    % 2. for each replicate, plot data
    for rep = 1:length(currData)
        
        % i. tick for figure counter
        counter = counter + 1;
        
        % ii. gather parameters
        repData = currData{rep,1};
        Vb = repData.cc(:,birthV);
        Vdiv = repData.cc(:,divV);
        dt = repData.meta(:,tau);
        
        % iii. keep only data within 95% of mean Vb
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
        Vdiv_trim1 = Vdiv(Vb < Vb_mean + 2*sigma);
        dt_trim1 = dt(Vb < Vb_mean + 2*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        Vdiv_trim2 = Vdiv_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        dt_trim2 = dt_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        clear dt_trim1 Vb_trim1 Vdiv_trim1 sigma repData Vb Vdiv dt Vb_mean 

        
        % iv. calculate lambda and inverse lambda
        lambda = (log(Vdiv_trim2./Vb_trim2))./dt_trim2;
        lambda_inverse = 1./lambda;
        
        % v. plot
        figure(counter)
        subplot(1,2,1)
        scatter(lambda_inverse,dt_trim2,'MarkerEdgeColor',rgb(palette{col}))
        axis([0 100 0 100])
        ylabel('interdivision time (min)')
        xlabel('inverse lambda (js)')
        title(strcat('T',num2str(currCond),': rep-',num2str(rep)))
        
        subplot(1,2,2)
        scatter(Vb_trim2,dt_trim2,'MarkerEdgeColor',rgb(palette{col}))
        axis([1 5 0 100])
        ylabel('interdivision time (min)')
        xlabel('birth volume')
        title(strcat('n=',num2str(length(dt_trim2))))
        
         % overlay fit
        fit = polyfit(Vb_trim2,dt_trim2,1);
        x = linspace(1,5,10);
        y = fit(1).*x + fit(2);
        subplot(1,2,2)
        hold on
        plot(x,y,'Color',rgb('SlateGray'))
        text(1.5,10,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))))
        
        
    end
end
clear currCond currData col clear lamb metric rep cond
clear conditionsOI Vb_trim2 dt_trim2 gr_inverse_trim2


