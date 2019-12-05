% Performs demixed principal components analysis (dPCA) on speech initiation epoch data to identify
% a condition-invariant signal (CIS), as in Figure 4 of Stavisky et al. eLife 2019.
% Also generates additional analysis details (e.g., variance explained plots, dPC component angles) and can analyze
% additional datasets, as in Figure 4-figure supplement 1A-E.
% Also can reduce dimensionality to a different number of dimensions as in Figure 4-figure supplement 2A.
% 
% USAGE: The script can be pointed to several available datasets and use different numbers of dimensions. 
% By default it will plot the T5-words dataset as in Figure 4, left column. To examine other neural 
% recordings, uncomment/modify the dataset in the 'Dataset specification' code block below. 
% To change the number of dPCs sought, change params.numDim in the 'Analysis parameters' code block below.  
%
% Note: To find the optimal regularization (lambda) and noise covariance (C_noise), dPCA
% uses single-trial neural data, which we cannot provide out of respect for our
% participant's privacy. Thus, these variables are provided pre-computed. Trial-averaged
% neural data are provided already soft-normalized (it would not make sense to change the
% soft-normalization offset without having access to the raw data, as lamba and C_noise
% would then be mismatched to the trial-averaged data).
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
%
%   This script uses the dPCA code pack associated with Kobak*, Brendel*, et al. 
%   "Demixed principal component analysis of neural population data", eLife 2016
%   (https://elifesciences.org/articles/10989). You can download this from
%   http://github.com/machenslab/dPCA.
% 
%   This script was developed and tested in MATLAB_R2017b, but it ought to be robust across
%   versions within a couple years of this.
%
% This code is associated with "Neural ensemble dynamics in dorsal motor cortex during
% speech in people with paralysis", eLife (2019).
%
% Copyright Sergey D. Stavisky, November 2019, Stanford Neural Prosthetics Translational
% Laboratory.


%% Dataset specification
dataset = 'T5-words_dynamics_go.mat';
% dataset = 'T5-syllables_dynamics_go.mat';
% dataset = 'T5-5wordsA_dynamics_go.mat';
% dataset = 'T5-5wordsB_dynamics_go.mat';
% dataset = 'T8-words_dynamics_go.mat';
% dataset = 'T8-syllables_dynamics_go.mat';


%% Analysis Parameters
params.numDims = 8; % main analysis (Figure 4)
% params.numDims = 12;  % can try sweeping range (Figure 4-figure supplement 2A)


%% Load the data
load( dataset, 'dat' ); % variable will be called 'dat'

% firingRates are: N x S x T
N = size( dat.firingRates, 1 ); %N is the number of neurons
S = numel( dat.labels ); % S is the number of conditions (word/syllable spoken)
T = numel( dat.t ); % T is the number of time-points 

margColours = [0 0 250; 250 0 0]./255; % red for condition-indepedent, blue for condition-dependent


%% do DPCA

% this is the code that would have been run with single-trial data:
haveSingleTrialData = false;
if haveSingleTrialData % so code more readable b/c not commented
    optimalLambda = dpca_optimizeLambda( dat.firingRates, dat.firingRatesIndividualTrial, dat.trialNum, ...
        'combinedParams', dat.combinedParams, ...
        'simultaneous', true, ...
        'numRep', 10, ...  %
        'filename', 'tmp_optimalLambdas.mat' );
    Cnoise = dpca_getNoiseCovariance( dat.firingRates, ...
        dat.firingRatesIndividualTrial, dat.trialNum, 'simultaneous', true);
end

% now do DPCA with the lambda and Cnoise precomputed
[W,V,whichMarg] = dpca( dat.firingRates , params.numDims, ...
    'combinedParams', dat.combinedParams, ...
    'lambda', dat.optimalLambda, ...
    'Cnoise', dat.Cnoise);


explVar = dpca_explainedVariance( dat.firingRates, W, V, ...
    'combinedParams', dat.combinedParams, ...
        'Cnoise', dat.Cnoise); % estimate signal variance
%% Plot dPCA results using the built-in dPCA plotting tools. 

timeEvents = 0;  % we don't have multiple times of interest
dpca_plot( dat.firingRates, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', dat.margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', dat.t,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);
%          'showNonsignificantComponents', false,          ...    
    
% Note: the significance star that dpca_plot uses in the upper-right triangle (denoting principal
% axes that are significantly non-orthogonal) use p<0.001. In our paper we use a more
% conservative approach, rejecting orthogonality at p<0.01. Thus, if you want to exactly
% match the significance test stars in figure 4-figure supplement 1C, you'll need to
% change 0.001 to 0.01 in line 366 of dpca_plot [they hard-coded that threshold].
% Also note that the dpca_plot function only plots a max of 6 components. If you want it to plot more, you'll need to go in
% and remove the hard-code of max 3 pltos per row, and 4x4 subplot (i.e., expand the
% figure and # subplots).


numCDdims = nnz( whichMarg == 1 );
numCIdims = nnz( whichMarg == 2 );
fprintf('%i condition-dependent dimensions and %i condition-INDEPENDENT dims together explain %.2f%% overall variance\n', ...
    numCDdims, numCIdims, explVar.cumulativeDPCA(end) );

% ----------------------------------------------------------------  
% re-sort based on how much condition-independent activity there is
CIdims = whichMarg==2;
CIdimsInds = find( CIdims );
[~, sortIndsByCI] = sort( explVar.margVar(2,CIdims), 'descend');

eachVarExplained = explVar.margVar(:,CIdimsInds(sortIndsByCI)); % start with the CI dimensions
eachVarExplained_forPlot = [eachVarExplained, [0;0]]; % used for plotting, adds space between CI and CD components
eachVarExplained = [eachVarExplained, explVar.margVar(:,~CIdims)];
eachVarExplained_forPlot = [eachVarExplained_forPlot, explVar.margVar(:,~CIdims)];
yTickLabels = arrayfun( @(x) {mat2str(x)}, 1:size( eachVarExplained, 2 ) );
yTickLabels = [yTickLabels(1:nnz(CIdims)) {''} yTickLabels(nnz(CIdims)+1:end)];
% also grab all the W and V vectors with this order
Wreordered = [W(:,CIdimsInds(sortIndsByCI)), W(:,~CIdims)];
Vreordered = [V(:,CIdimsInds(sortIndsByCI)), V(:,~CIdims)];

% Flip CIS1 if necessary to show it as upwards-going (it's AU and thus arbirary direction anyway, but
% this keeps it from flip-flopping for different number dimensionalities and looking very different for 
% an arbtirary sign flip (this rarely happens)
acrossCondsMean = squeeze( mean( dat.firingRates, 2 ) )'; % T x E
GM = acrossCondsMean * Wreordered(:,1);
if GM(end) < GM(1)
    fprintf( 'Flipping sign of CI1 to maintain upward-going convention\n')
    Wreordered(:,1) = -Wreordered(:,1);    
    Vreordered(:,1) = -Vreordered(:,1);    
end

% What fraction of these dPC's variance does CIS 1 explain? Note I'm including both its CI and CD
% marginalization (latter ends up being tiny though)
totVarDPCA = explVar.cumulativeDPCA(end);
varCIS1 = sum( eachVarExplained(:,1) );
CIS1_CIpercent = 100*eachVarExplained(2,1)/varCIS1;
fprintf('CIS1 (which is %.1f%% condition-independent) explains %.1f%% of top %i dPCs'' variance and %.1f%% of full (%i-dim) variance\n', ...
    CIS1_CIpercent, 100*varCIS1/totVarDPCA, numCDdims+numCIdims, varCIS1, N );


% ----------------------------------------------------------------  
% How orthogonal is CIS_1 to the movement dims?
CIS1_V = Vreordered(:,1);
CDdims_V = Vreordered(:,numel(CIdimsInds)+1:end);

for i = 1: size( CDdims_V, 2 )
    angleBetween = angleBetweenVectors(CIS1_V,CDdims_V(:,i));
    if angleBetween > 90
        angleBetween = 180- angleBetween;
    end
    fprintf('  Angle between CIS1 and CD%i is %.1f deg\n', i,  angleBetween)
end



%% ----------------------------------------------------------------  
% Plot how much each component explains    
% (Presented as in Kaufman et al. 2016)
figh = figure;
figh.Color = 'w';
figh.Name = sprintf('CIS dPCA %s', dataset);
axh_var = subplot( 2, 1, 1 );
axh_var.TickDir = 'out';
hbar = barh( eachVarExplained_forPlot', 'stacked' );
axh_var.YTickLabel = yTickLabels;

axh_var.YDir = 'reverse';
hbar(1).BarWidth = 1;
hbar(2).BarWidth = 1;
hbar(2).FaceColor = margColours(2,:);
hbar(1).FaceColor = margColours(1,:);

title( dataset, 'Interpreter', 'none' )
ylabel('dPC Component');
xlabel('Variance Captured')
box off

% ----------------------------------------------------------------  
% Prep for plotting
% Define the specific colormap
colors = [];
legendLabels = {};
for iLabel = 1 : numel( dat.labels )
   colors(iLabel,1:3) = labelColors( dat.labels{iLabel} ); 
   legendLabels{iLabel} = sprintf('%s', dat.labels{iLabel} );
end


% ----------------------------------------------------------------  
% Plot each condition's CIS_1
axh_CIS = subplot( 2, 1, 2 ); hold on;
CIS1dim = Wreordered(:,1);
for iLabel = 1 : numel( dat.labels )
    myCIS1 = squeeze( dat.firingRates(:,iLabel,:) )' * CIS1dim; % T x 1 
    plot( dat.t, myCIS1, 'Color', colors(iLabel,:), ...
        'LineWidth', 1 );
end

axh_CIS.TickDir = 'out';
xlim( [dat.t(1) dat.t(end) ] );
axis tight
xlabel( sprintf( 'Time after %s (s)', dat.params.alignEvent ) );
ylabel( 'CIS_1 activity (AU)' )