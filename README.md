# super-size
Scripts plotting the effect of nutrient fluctuations on cell size


## Resumption

1. loglog_single_expt.m
   - loglog plots of Vb_i vs. tau_i for each individial experiment
   - result: all fluc conditions behave like steady, with a power law except for T=60min
   - in consideration for final paper

2. loglog_single_expt_Vb_lambda.m
   - question: does T=60min preserve relationship between size and growth rate?
   - answer: no, T=60min single cells fall off the line by other conditions
   - not currently under consideration

3. measured_v_expected.m
   -


## Figures:

1A. mean birth size vs mean growth rate of each nutrient condition
    size metrics: birth volume of cell cycle, birth volume of next cell cycle
                  birth length of cell cycle, birth length of next cell cycle
                  birth width of cell cycle, birth width of next cell cycle
                  birth SA:V ratio of cell cycle, birth SA:V ratio of next cell cycle



## Supplementary Figures:





## Data processing and organization:

  All figures use input data compiled into a data structure by the following scripts:

  i. size_metrics (to be organized from figure1A_division.m)
     This script loops through data files from each experimental replicate to compile a matrix of cell size measurements, organized by cell cycle.

     currently, figure1A_division.m saves a file " A1_ccSize.mat "
     this includes a data structure "compiled_data" which includes 4 matrices:
        1) cc: 8 parameters (col)
           each row, values of each parameter for current cell cycle (cc)
        		1. Vbirth
        		2. Vdiv
        		3. Lbirth
        		4. Ldiv
        		5. Wbirth
        		6. Wdiv
        		7. SA2Vbirth
        		8. SA2Vdiv 

        2) meta_data: 5 parameters (col)
           each row, values of each parameter for current cell cycle (cc)
        		1. lambdas
        		2. tau_all/60
        		3. tracks_all
        		4. Tbirth_all
        		5. Tdiv_all
                
        3) plus1_data: 4 parameters (col)
           each row, values of each parameter at birth of next cell cycle (cc + 1)
        		1. plus1_vol
        		2. plus1_len
        		3. plus1_width
        		4. plus1_sa2v
                      
        4) mus_insta: each row = cell array of mu values for one individual e. coli



## Figures to figure stuff out

saved as ss...(#).m
1. 
2. simulate growth at different fluctuation timescales with measured mu(t)
3. probability of division over time in single shift environments
4. individual volume trajectories starting with cell experiencing a single nutrient upshift
5. individual volume trajectories starting with cell experiencing a single nutrient downshift
6. 
7. individual interdivision time vs individual growth rate (lambda) or birth volume
8. tau_i vs Vb_i for each replicate of T = 15 and 60 min, same plot with fit slopes
9. predicting V_mean from linear fit of tau_i vs Vb_i
10. tau_i vs Vb_i with data compiled from all replicates for each condition T = 15 and 60 min



## Code for on-going analyses

### ss1
1. histograms of Vb_i for each experimental replicate
2. histograms of tau for each experimental replicate

### ss2 
1. data vectors of mu(t) from T = 15 and 60 min.
- mu(t) is the mean of all mu_i calculated from all cells from all replicate conditions.
2. MATLAB figure of mu(t) from T = 15 and 60 min plotted over the nutrient period (high:low:low:high)

### ss3
1. probability of division in single shift experiments

### ss4
1. volume over time of single cell cycles that experience a single nutrient upshift

### ss5
1. single lineage volume tracks starting from cell cycle experiencing a single upshift

### ss6
1. probability of division as a function of nutrient score

### ss7
1. seven plots (4 for T = 15 min and 3 for T = 60 min) of tau_i vs Vb_i.
- n = number of individial cells (division to division event) that contribute to the scatter
- subplot 1 is tau_i vs lambda_i
- subplot 2 is tau_i vs Vb_i
- linear fit to scatter is displayed for A and B

### ss8
1. scatter of tau_i and birth size_i of each replicate overlaid
Q: are linear fit slopes different between conditions?
A: yes, T=15 min slopes are steeper! (slighter slopes suggest larger sizes)

### ss9
1. does tau_i(Vb_i) predict Vb_mean?
Q: can we predict V_mean from mu_mean and linear fit of tau_i and Vb_i? 
A: pretty ok, but no clear difference between timescales. suggests a transition?

### ss10
1. tau_i vs Vb_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)
scatter of tau_i and Vb_i 

2. ss10_a.m
data for individual replicates per condition, grouped by timescale of fluc condition
metrics: slopes and intercepts of tau_i vs. Vb_i, tau_i, lambda_i, and number of cells
looks for variability between replicates and plots slope vs lambda, tau and Vb

### ss11
1. tau_i vs delta_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)
scatter of tau_i and delta_i 
delta_i = Vdiv_i - Vb_i

### ss12
1. delta_i vs Vb_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)

### ss13
1. Vdiv_i vs Vb_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)

### ss14
single-cell determination of slope between tau_i and Vb_i.
plot single-cell tau_i and Vb_i over time in steady environments:
1. steady low
2. steady ave
3. steady high

### ss15
single-cell determination of slope between tau_i and Vb_i.
plot single-cell tau_i and Vb_i across transitions:
1. steady-state transitions (low to high)
2. steady-state transitions (high to low)
3. steady-state low to fluctuating (T = 60min) 
Slope transitions quickly between conditions

### ss16
is single-cell slope directly related to average whole dataset slope?
answer: seems so! though numbers are not a direct match
-compare distributions of single cells (tau_i divided by Vb_i) with distribution of replicate averages
-comparable distributions suggest a difference from size vs growth rate plot in Taheri which shows difference between averaged and individual level trends

### ss17
for each single cell, divide tau_i by Vb_i
plot single cell scatter of tau_i/Vb_i vs Vb_i for each replicate condition, and mean point from each (like ss10a)
Q: is individual trend different or similar to average? 
A: very similar!!

### ss18
demonstration of other metrics with mismatch between single cell and average trends
similar to ss17, except plotting Vb vs lambda
Q: is individual trend different or similar to average? 
A: no! more blobby

### ss19
identical to ss18, except plotting log_10(Vb) instead of Vb
Q: is individual trend different or similar to average? 
A: no! exaggerates flatter slope of size vs growth rate relationship

### ss20
does single-cell slope_i vs time tell us why 15 min and 60 min sizes are different?
similar to ss15 but with t15 and t60 data
Yes??! slope_i does not fluctuate with T = 15 min, but yes in T = 60min


### ss21
one might think that our collapsed data may result from the fact that we are plotting hyperbolas.
this should not be true, so let's confirm that it's really that alpha and tau_o covary to fit our hyperbole.
- quantify alpha and tau_o for each condition with a scatter of tau_i vs Vb_i
- for each average tau_o and alpha, plot hyperbola:
- for example: take one alpha, fix it and use all measured tau_o. should fall off.


### ss22
plot tau_i vs Vb_i and average tau vs average Vb (single-cell vs population per condition) as suggested by Johannes.
- similar to ss17
- used in update Oct 11


### ss23
plot slope_i vs Vb_i and average tau vs average Vb (single-cell vs population per condition) as suggested by Johannes.
- similar to ss17
- used in update Oct 11


### ss24
bootstrap hypothesis testing finds that fit slope across tau_i vs Vb_i per replicate is a signficant measure.
- builds from ss23.m
- used in update Oct 12


### ss25
Test hyperbolic fits using residuals R and varying c in y=c/x.
- used in update Oct 12




## Retired code

1. Figure 6: simulations with varied growth rate function and fixed tau function
2. Figure 7: simulations with different intialize size