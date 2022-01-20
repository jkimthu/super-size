%% simulation1. variation in V_b


%  Question:
%  - How can we have the same growth rate and interdivision time but different slopes?

%  Strategy: 
%  - Simulate data and plot tau vs. V_b for two populations, A and B
%  - A and B have equal mean and standard deviation in tau = equal tau.
%  - However, A and B have different standard deviations in tau.
%
%  - Simulate:
%    1. A and B have equal mean V_b, B has a larger standard deviation in V_b
%    2. A and B have equal mean V_b, B has a smaller standard deviation in V_b


%  Last edit: Jen Nguyen, 2021 Oct 2
%  Commit: simulate effect of varying V_b on slope

%  OK let's go!


% 1. Define mean and s.d. for tau; generate 1000 "single-cell" values
tau_mean = 48; % min
tau_sd = 8; %min





