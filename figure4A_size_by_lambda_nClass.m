%% Division size vs growth rate, by nutrient signal

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')
load('A9_data_volume.mat')
load('A6_data.mat','compiled_lambda')
condition = 1;


% 0. initialize plotting parameters
%palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','ForestGreen','Orange','GoldenRod'};
palette = {'SteelBlue','LightSkyBlue','DeepSkyBlue'};
shape = 'o';


% 0 . initialize plot limits
plot_tau = 90;
plot_delta = 6;
plot_divSize = 11;
max_lambda = 5;


% 0. initialize signal classifications
onlyH = 1;
onlyL = 2;
H2L = 3;
L2H = 4;
HLH = 5;
LHL = 6;
HLHL =7;
LHLH = 8;
sigArray = [onlyH; onlyL; H2L; L2H; HLH; LHL; HLHL; LHLH];
sigNames = {'onlyH','onlyL','H2L','L2H','HLH','LHL','HLHL','LHLH'};




% 1. loop through experiments to format data and plot
for ee = 11:length(exptArray)
    
    
    % 2. initialize experiment meta data
    index = exptArray(ee);
    date = storedMetaData{index}.date;
    timescale = storedMetaData{index}.timescale;
    disp(strcat(date, ': analyze!'))
    clear index
    
    
    % 3. isolate experiment data
    %eeTaus = compiled_tau{ee}{condition}./60; % convert from sec to h
    eeMus = compiled_lambda{ee}{condition};
    eeSignals = compiled_signal{ee}{condition};
    eeDeltas = compiled_addedSize{ee}{condition};
    eeDivSizes = compiled_divSize{ee}{condition};
    
    
    % 4. calculate lambda and nScores from compiled mus
    eeLambdas = cellfun(@nanmean,eeMus);
    clear eeMus
    
    
    % 5. remove data with too large growth rates
    %taus = eeTaus(eeLambdas <= max_lambda);
    nSignals = eeSignals(eeLambdas <= max_lambda);
    lambdas = eeLambdas(eeLambdas <= max_lambda);
    deltas = eeDeltas(eeLambdas <= max_lambda);
    divSizes = eeDivSizes(eeLambdas <= max_lambda);
    clear eeNscores eeLambdas eeSignals eeTaus eeDeltas eeDivSizes
    
    
    % 6. for each curve, determine signal type
    signalType = zeros(length(nSignals),1);
    if timescale == 3600
        
        for cc = 1:length(nSignals)
            
            currSignal = nSignals{cc};
            ds_dt = diff(currSignal);
            currShifts = length(find(ds_dt ~= 0));
            
            if currShifts == 0
                if currSignal(1) == 1
                    signalType(cc) = onlyH;
                else
                    signalType(cc) = onlyL;
                end
            elseif currShifts == 1
                if currSignal(1) == 1
                    signalType(cc) = H2L;
                else
                    signalType(cc) = L2H;
                end
            elseif currShifts == 2
                if currSignal(1) == 1
                    signalType(cc) = HLH;
                else
                    signalType(cc) = LHL;
                end
            elseif currShifts == 3
                if currSignal(1) == 1
                    signalType(cc) = HLHL;
                else
                    signalType(cc) = LHLH;
                end
            end
        end
        
    else
        error('wrong signal classifier for current timescale!')
    end
    clear currSignal ds_dt currShifts cc
    
    
    
    % 7. accumulate cell cycle phenotypes by signal class
    %tau_binned = accumarray(signalType,taus,[],@(x) {x});
    delta_binned = accumarray(signalType,deltas,[],@(x) {x});
    divSize_binned = accumarray(signalType,divSizes,[],@(x) {x});
    lambda_binned = accumarray(signalType,lambdas,[],@(x) {x});
    
    
    % 8. calculate mean and standard error for each signal class
    lambda_mean = cellfun(@mean,lambda_binned);
    lambda_std = cellfun(@std,lambda_binned);
    lambda_count = cellfun(@length,lambda_binned);
    lambda_sem = lambda_std./sqrt(lambda_count);
    
    delta_mean = cellfun(@mean,delta_binned);
    delta_std = cellfun(@std,delta_binned);
    delta_count = cellfun(@length,delta_binned);
    delta_sem = delta_std./sqrt(delta_count);
    
    divSize_mean = cellfun(@mean,divSize_binned);
    divSize_std = cellfun(@std,divSize_binned);
    divSize_count = cellfun(@length,divSize_binned);
    divSize_sem = divSize_std./sqrt(divSize_count);
    
    
    % 9. plot mean and standard error of each signal class
    figure(1)
    for sc = 3:5%1:length(delta_mean)
        color = rgb(palette{sc-2});
        if isempty(delta_mean(sc)) == 1
            continue
        else
            errorbar(lambda_mean(sc),delta_mean(sc),lambda_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            hold on
            errorbar(lambda_mean(sc),delta_mean(sc),delta_sem(sc),'horizontal','Color',color)
            
   
        end
    end
    xlabel('Mean growth rate')
    ylabel('Added Volume')
    axis([0.7 2 1.5 4.5])
    
    figure(2)
    for sc = 3:5%1:length(delta_mean)
        color = rgb(palette{sc-2});
        if isempty(delta_mean(sc)) == 1
            continue
        else
            errorbar(lambda_mean(sc),divSize_mean(sc),lambda_sem(sc),'Marker',shape,'MarkerSize',10,'Color',color)
            hold on
            errorbar(lambda_mean(sc),divSize_mean(sc),divSize_sem(sc),'horizontal','Color',color)
            
        end
    end
    xlabel('Mean growth rate')
    ylabel('Division Volume')
    
    figure(3)
    for sc = 3:5%1:length(delta_mean)
        color = rgb(palette{sc-2});
        if isempty(delta_mean(sc)) == 1
            continue
        else
            %errorbar(lambda_mean(sc),log(divSize_mean(sc)),log(divSize_sem(sc)),'Marker',shape,'Color',color)
            plot(lambda_mean(sc),log(divSize_mean(sc)),'Marker',shape,'Color',color)
            hold on
        end
    end
    xlabel('Mean growth rate')
    ylabel('log Division Volume')
    
    % add rough estimate of relationship between
    figure(3)
    hold on
    plot([0.5 3],[1.1 2.38],'Color',rgb('Silver'))
    clear color sc tau_mean tau_sem tau_std tau_count delta_mean delta_sem delta_std delta_count
    %clear divSize_mean divSize_sem divSize_std divSize_count
    clear signalType taus nSignals nScores delta_binned delta tau_binned
    
end
clear onlyH onlyL H2L L2H HLH LHL HLHL LHLH


% 10. plot expected tau = 1/lambda relationship
% x = linspace(0.1,plot_delta,100);
% y = 60./x; % convert from h to min
% 
% figure(1)
% hold on
% plot(x,y,'Color',rgb('SlateGray'))
