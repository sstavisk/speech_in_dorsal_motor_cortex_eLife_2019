% Performs jPCA on speech production epoch data to identify
% rotatory dynamics, as in Figure 5 of Stavisky et al. eLife 2019.
% Also returns angle comparison between the top jPCA plane and the CIS_1 as in figure
% 4-figure supplement 1E.
% Also can analyze additional datasets, as in Figure 4-figure supplement 1F,G, and can be run on the post-prompt 
% epoch as in Figure 4-figure supplement 1H.
% Also can reduce dimensionality to a different number of dimensions as in Figure 4-figure supplement 2B.
% Also can reduce dimensionality to a different number of dimensions as in Figure 4-figure supplement 2B.
%
% Heads-up that the script will also plot one example surrogate dataset's jPCA neural state. This
% plot, which has a pale red background, will look very not-rotational. That's to be
% expected, since the whole point is that surrogate data doesn't show rotational
% structure.
% 
% USAGE: The script can be pointed to several available datasets and use different numbers of dimensions. 
% By default it will plot the T5-words dataset as in Figure 5, top row. To examine other neural 
% recordings, uncomment/modify the dataset in the 'Dataset specification' code block below. 
% To change the number of jPCA planes sought, change params.numDim in the 'Analysis parameters' code block below.  
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
%
%   This script uses the jPCA code pack associated with Churchland*, Cunningham*, et al. 
%   "Neural population dynamics during reaching", Nature 2012
%   (https://www.nature.com/articles/nature11129). You can download this from
%   https://churchland.zuckermaninstitute.columbia.edu/content/code.
%   Note that this version's pca.m uses the deprecated princomp function, which causes so
%   many warnings that it dramatically slows execution. I've included a jPCA_updated.m in
%   /helper_funcitons/ that swaps out princomp for pca to avoid this. 
%
%   This script also uses the Tensor-Maximum-Entropy (TME) hypothesis testing method from
%   Elsayed and Cunningham, "Structure in neural population recordings: significant or epiphenomenal?",
%   Nature Neuroscience, 2017 (http://www.nature.com/articles/nn.4617). You can download this from
%   https://github.com/gamaleldin/TME.
%
%   This script uses subspacea.m (already provided in /helper_functions), described in A. V. Knyazev 
%   and M. E. Argentati, "Principal Angles between Subspaces in an A-Based Scalar Product: 
%   Algorithms and Perturbation Estimates. SIAM Journal on Scientific Computing, 
%   23 (2002), no. 6, 2009-2041." (http://epubs.siam.org/sam-bin/dbq/article/37733)
%   which is also available at
%   https://www.mathworks.com/matlabcentral/fileexchange/55-subspacea-m.
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
dataset = 'T5-words_dynamics_speak';
% dataset = 'T5-syllables_dynamics_speak';
% dataset = 'T5-5wordsA_dynamics_speak';
% dataset = 'T5-5wordsB_dynamics_speak';
% dataset = 'T8-words_dynamics_speak';
% dataset = 'T8-syllables_dynamics_speak';

% Do analysis on audio prompt epoch (Figure 4-figure supplement 1G):
% dataset = 'T5-words_dynamics_prompt';
% dataset = 'T5-syllables_dynamics_prompt';
% dataset = 'T5-5wordsA_dynamics_prompt';
% dataset = 'T5-5wordsB_dynamics_prompt';
% dataset = 'T8-words_dynamics_prompt';
% dataset = 'T8-syllables_dynamics_prompt';


%% Analysis Parameters
params.numDims = 6; % main analysis (Figure 4)
% params.numDims = 2;  % can try sweeping range (Figure 4-figure supplement 2B). Must be a multiple of 2.

% Used for the Churchland/Cunningham jPCA code pack
params.jPCA_params.softenNorm = 10;
params.jPCA_params.suppressBWrosettes = 1;
params.jPCA_params.suppressHistograms = 1;
params.jPCA_params.meanSubtract = 1;
params.jPCA_params.numPCs = params.numDims; % inherit from above

% Used for the Elsayed/Cunningham TME-master code pack
params.numSurrogates = 1000; % how many surrogate datasets to compare to PAPER
rng(1); % consistent random seed

if rem( params.numDims, 2 ) 
    error('params.numDims must be a multiple of 2 because jPCA searches for *planes*\n')
end

%% Load the data
load( dataset ); % key variable loaded is called 'dat'
% neural data is already arranged the way jPCA.m wants it, in a struct array, one 
% element per condition. E.g., Data(1).A has the condition-averaged firing rates, arranged as
% a time x channels matrix. The corresponding times, in ms, are in  Data(1).times (these
% are  with respect to dat.params.alignEvent).


%% Do jPCA
[Projection, Summary] = jPCA_updated( dat.Data, dat.params.dataTimestampsInMS, params.jPCA_params );
fprintf('%i PCs used capture %.2f%% overall variance. Each PC''s variance is: %s\n', ...
    params.jPCA_params.numPCs, 100*sum(Summary.varCaptEachPC(1:end)), mat2str( 100.*Summary.varCaptEachPC , 4 ) )

%% Plot the neural state trajectory in the top jPCA plane
plotParams.planes2plot = [1]; % first plane
for iLabel = 1 : numel( dat.labels )   
    plotParams.circleColors{iLabel} = labelColors( dat.labels{iLabel} );
end
[figh, cmap] = makeJPCAplot(  Projection, Summary, plotParams, dat.labels );
figh.Name = sprintf('Top jPCA plane trajectories %s', dataset );



%% Significance testing against surrogate data with primary statistics matched to real data
% Uses method of Elsayed and Cunningham, Nature Neuroscience, 2017.

% Rearrange the Data matrix to a dataTensor as expected by Gamal's code. 
% dataTensor is time x neuron x condition.
T = numel( dat.params.dataTimestampsInMS );
N = size( dat.Data(1).A, 2 );
C = numel( dat.Data );
dataTensor = nan( T, N, C );
tensorTimestamps = dat.params.dataTimestampsInMS;
for iC = 1 : C
    dataTensor(:,:,iC) = dat.Data(iC).A;
end


% quantify primary features of the original data
[targetSigmaT, targetSigmaN, targetSigmaC, M] = extractFeatures(dataTensor);

GCparams.margCov{1} = targetSigmaT;
GCparams.margCov{2} = targetSigmaN;
GCparams.margCov{3} = targetSigmaC;
GCparams.meanTensor = M.TNC;
maxEntropy = fitMaxEntropy( GCparams );  % fit the maximum entropy distribution

% These are the summary statistics I'll be saving from each surrogate run
R2_Mbest_surr = nan( params.numSurrogates, 1 );
R2_Mskew_surr = R2_Mbest_surr;
R2_MbestPrimary2d_surr = R2_Mbest_surr;
R2_MskewPrimary2d_surr = R2_Mbest_surr;
surrJPCAparams = params.jPCA_params;
surrJPCAparams.suppressText = 1;
fprintf( 'Generating surrogate datasets and performing jPCA on them. This should take less than a minute... ' )
tic
for iSurr = 1 : params.numSurrogates
    surrTensor = sampleTME(maxEntropy);       % generate TME random surrogate data.
    surrData = dataTensorToDataStruct( surrTensor, tensorTimestamps );
    
    % Run same jPCA analysis on this surrogate data
    [surProjection, surSummary] = jPCA( surrData, dat.params.dataTimestampsInMS, surrJPCAparams );
    R2_Mbest_surr(iSurr) = surSummary.R2_Mbest_kD;
    R2_Mskew_surr(iSurr) = surSummary.R2_Mskew_kD;
    R2_MbestPrimary2d_surr(iSurr) = surSummary.R2_Mbest_2D;
    R2_MskewPrimary2d_surr(iSurr) = surSummary.R2_Mskew_2D;
    
    % Make the jPCA plot (for just one surrogate dataset)
    if iSurr == 1
        [figh, cmap] = makeJPCAplot(  surProjection, surSummary, plotParams, dat.labels );
        figh.Name = sprintf('SURROGATE jPCA trajectories %s', dataset );
        figh.Color = [0.9 0.8 0.8];
        axh = figh.Children;
    end
end
fprintf('DONE (took %.1f seconds)\n', toc )

% evaluate a P value
trueMetric = Summary.R2_Mskew_kD;
surrMetrics = R2_Mskew_surr;
P = 1 - nnz(trueMetric > surrMetrics ) / (params.numSurrogates + 1 );
fprintf('True R2_Mskew_kD is %g, which is > %i/%i surrogates (P value = %g). Mean of surrogotes = %g\n', ...
    trueMetric, nnz( trueMetric > surrMetrics ), numel( surrMetrics ), P, mean( surrMetrics ) )
fprintf('True R2_Mbest_kD is %g,  Mean of surrogotes = %g\n', ...
    Summary.R2_Mbest_kD, mean( R2_Mbest_surr ) )
fprintf('True R2_Mskew_2D is %g,  Mean of surrogotes = %g\n', ...
    Summary.R2_Mskew_2D, mean( R2_MskewPrimary2d_surr ) )
fprintf('True R2_Mbest_2D is %g,  Mean of surrogotes = %g\n', ...
    Summary.R2_Mbest_2D, mean( R2_MbestPrimary2d_surr ) )


%% plot null distribution
x = 0:0.03:1;
h = hist(surrMetrics, x);
figh = figure;
figh.Name = sprintf('jPCA comparison to surrogate distribution %s', dataset );
set(figh, 'color', [1 1 1]);
hold on
box on
hb = bar(x, h);
set(hb,'facecolor',[0.5000    0.3118    0.0176],'barwidth',1,'edgecolor','none')
trueline = line( [trueMetric trueMetric], [0 1.1*max( hb.YData )], ...
    'Color', [0.20 0.75 0.95], 'LineWidth', 2 );
xlabel(sprintf('M_{skew,%i-D}  R^2', params.numDims ))
ylabel('Count')
xlim([0 1])
set(gca, 'FontSize',12)
set(gca, 'xtick',[-1 0 1])
legend([trueline, hb], {'original data', 'surrogate'})
legend boxoff
box off

%% Compare jPCs to CIS1 (figure 4-figure supplement 1E)
% For convenience, the speech initiation epoch CIS1_V calculated by figure4_conditionInvariantDynamics.m is saved in
% the dataset files for all the speach onset epoch datasets (e.g.
% T5-words_dynamics_speak.mat). If one of those datasets is run, then this code cell will
% compare the CIS1 with the jPC dimensions
if exist( 'CIS1_V', 'var' ) && ~contains( dataset, 'prompt' ) % latter condition avoids mistakenly using CIS1_V from a previous dataset if you didn't clear workspace
    numJPCs = size( Summary.jPCs_highD, 2 );
    suba = rad2deg( subspacea( CIS1_V, [Summary.jPCs_highD(:,1), Summary.jPCs_highD(:,2)] ) );
    fprintf('Subspace angle between CIS1 and first 2 jPC dims is %.1f deg\n', ...
        suba );
    
     % How orthogonal is CIS_1 to the jPC dims? 
    for i = 1 : numJPCs
        myJPC = Summary.jPCs_highD(:,i);
        angleBetween = angleBetweenVectors( CIS1_V, myJPC);
        if angleBetween > 90
            angleBetween = 180- angleBetween;
        end
        % statistical test for substantial snon-orthogonality
        % taken from Kobak et al 2016 code (plot_dpca.m)
        [~, psp] = corr( CIS1_V, myJPC, 'type', 'Kendall');
        fprintf('  Angle between CIS1 and JPC%i is %.1f deg, p-value that this is significantly non-orthogonal = %f\n', ...
            i,  angleBetween, psp)
    end

    % plot the disk-and-rod
    figh = figure; 
    figh.Name = sprintf('CIS1 vs jPCs 1 and 2 %s', dataset );
    circle( [0,0],1,[.5 .5 .5], 1 )
    axh = figh.Children;
    ch = axh.Children; % get handle for the disk
    ch.Color = 'b';
    box off
    axh.XAxis.Visible = 'off';
    axh.YAxis.Visible = 'off';
    axh.ZAxis.Visible = 'off';    
    axis equal
    rayx = cos( deg2rad( suba ) );
    rayy = sin( deg2rad( suba ) );
    lh = line( [0 0], [0 rayx], [0 rayy], 'LineWidth', 2', 'Color', 'r');
    view(3);
    title( figh.Name, 'Interpreter', 'none' );
    legend('jPC_{1,2}', 'CIS_1')
end
