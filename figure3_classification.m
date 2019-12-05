% Plots single-trial neural feature vectors (using 2D tSNE) as in Figure 3A, and does single-trial classification 
% of these neural features (to predict when syllable or word was heard/spoken) as in Figure 3B,C of 
% Stavisky et al. eLife 2019.
%
% This script uses pre-assembled single-trial neural features (spikes and high-frequency LFP power in 100 ms
% bins) which were generated from the raw neural data.
%
% USAGE: The script can be pointed to several available datasets. By default it will
% analyze the participant T5, syllables dataset  (Figure 3A, 3B top-left panel). To examine other datasets,
% uncomment/modify the dataset selection in the 'Dataset specification' code block below.
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
% 
%   This script was developed and tested in MATLAB_R2017b, but it ought to be robust across
%   versions within a couple years of this.
%
%   This script uses the MATLAB statistics toolbox.
%
% This code is associated with "Neural ensemble dynamics in dorsal motor cortex during
% speech in people with paralysis", eLife (2019).
%
% Copyright Sergey D. Stavisky, November 2019, Stanford Neural Prosthetics Translational
% Laboratory.

rng(1); % consistent random seed

%% Dataset specification
% Speak epoch decoding
dataset = 'T5-syllables_classify_speak'; % panel A and  panel B top-left
% dataset = 'T8-syllables_classify_speak'; % panel B top-right
% dataset = 'T5-words_classify_speak'; % panel B bottom-left
% dataset = 'T8-words_classify_speak'; % panel B bottom-right

% Prompt epoch decoding (summarized as gray bars in panel C)
% dataset = 'T5-syllables_classify_prompt';
% dataset = 'T5-words_classify_prompt';
% dataset = 'T8-syllables_classify_prompt';
% dataset = 'T8-words_classify_prompt';


%% Analysis parameters

% tSNE parameters
params.tsne.options.MaxIter = 10000;
params.tsne.options.OutputFcn = [];
params.tsne.options.TolFun = 1e-10;
params.tsne.NumPCAComponents = 0;
params.tsne.Standardize = true; % normalizes (z-score) input data;  good given different feature types (spike counts and HLFP)
params.tsne.Perplexity = 15;

% SVM parameters
params.svn.outlierFraction = 0.05; 

% Shuffle labels (to compute chance levels)
% params.numShuffles = 101; % used for shuffle test statistics in paper; slow.
params.numShuffles = 0; % Speeds up the script to get the main (real data) result faster.

% PCA across electrodes, on trial-averaged data
params.numPCs = [];



%% Load the data
load( dataset, 'dat' ); % variable will be called 'dat'
uniqueLabels = unique( dat.labels ); % conditions 
numTrials = numel( dat.labels );
% final display order (matches paper figures)
if contains( dataset, 'syllables' )
    desiredOrder = {'silence', 'i', 'u', 'ae', 'a', 'ba', 'ga', 'da', 'k', 'p', 'sh'};
else
    desiredOrder = {'silence', 'beet', 'bat', 'bot', 'boot', 'dot', 'got', 'shot', 'keep', 'seal', 'more'};
end


% -------------------------------------------------------------------------------------------------------
%% Single trial neural feature vector tSNE plot: Panel A
% -------------------------------------------------------------------------------------------------------
fprintf('Performing tSNE (this takes a few seconds)... ')
[Y, loss] = tsne( dat.neuralFeatures, ...
    'NumPCAComponents', params.tsne.NumPCAComponents, 'Standardize', params.tsne.Standardize, ...
    'verbose', 0, 'Options', params.tsne.options, 'Algorithm', 'exact', 'Distance', 'euclidean', ...
    'Perplexity', params.tsne.Perplexity );
fprintf('Done. KL Loss = %f\n', loss);

% Plot tSNE 
figh = figure;
figh.Color = 'w';
figh.Name = sprintf( 'Single trial neural features %s', dataset );
axh = axes;
title( figh.Name );
colors = cell2mat( cellfun( @labelColors, dat.labels, 'UniformOut', false ) );
hscat = gscatter( Y(:,1), Y(:,2), dat.labels, colors );
title( figh.Name, 'Interpreter', 'none' );
axis equal;
axis off;
drawnow; % plots tSNE before doing classification; gives you something to look at while waiting.
pause( 0.1 )
% -------------------------------------------------------------------------------------------------------
%% Single trial classification: Panel B
% -------------------------------------------------------------------------------------------------------
% Do the leave-one-out SVM classification
tic
fprintf('Starting leave-one-out classification of %i trials ... this could take a few minutes...\n  ', numTrials )
realPredictions = classifySpeech( dat.neuralFeatures, dat.labels, params.svn );
fprintf( ' DONE (took %.1f minutes)\n', toc/60)

% fraction success?
result.numSuccessful = nnz( strcmp( dat.labels, realPredictions ) );
result.successRate = result.numSuccessful / numTrials;
% Report results
fprintf('Classification accuracy for %s = %.1f%%\n', dataset, 100*result.successRate  );
    
% Confusion Matrix
[result.confuseMat, result.confuseMatLabels] = confusionmat( dat.labels, realPredictions );


%% Plot Confusion Matrix
% re-order based on specified at top order
newOrder = cellfun( @(x) find( strcmp( result.confuseMatLabels, x) ), desiredOrder, 'UniformOutput', false ); % handles missing label in T5-syllables
newOrder(cellfun( @isempty, newOrder )) = [];
newOrder = cell2mat( newOrder );
for row = 1 : numel( uniqueLabels )
    for col = 1 : numel( uniqueLabels )
        orderedConfuseMat(row,col) = result.confuseMat(newOrder(row), newOrder(col) );
    end
end
% Normalize so 100% is number of true labels
orderedConfuseMat = 100*(orderedConfuseMat./ repmat( sum( orderedConfuseMat, 2 ), 1, numel( uniqueLabels ) ));

figh = figure;
figh.Color = 'w';
figh.Name = sprintf( '%s confusion matrix', dataset );
title( sprintf( '%s: %.1f%% correct', dataset, 100*result.numSuccessful/numel( dat.labels )  ), 'Interpreter', 'none' )
imagesc( orderedConfuseMat, [0 100] );
colormap('bone')
axis square

axh = figh.Children;
axh.TickLength = [0 0];
axh.XTickLabel = result.confuseMatLabels(newOrder);
axh.YTickLabel = result.confuseMatLabels(newOrder);
xlabel('Decoded syllable/word');
ylabel('Spoken syllable/word');
cbarh = colorbar;
cbarh.TickDirection = 'out';
ylabel(cbarh, '% of spoken');


%% Shuffle control 
if params.numShuffles > 0
    tic
    fprintf('Performing shuffle tests. These can take multiple hours if params.numShuffles is larger\n')
    fprintf('  shuffle     ')
    shuffledPredictions = cell( numTrials, params.numShuffles );
    shuffledAccuracies = nan( 1, params.numShuffles ); % 1 x numShuffles
    shuffledConfusedMat = nan( numel( uniqueLabels ), numel( uniqueLabels ), params.numShuffles ); % numClasses x numClasses x numShuffles

    for iShuffle = 1 : params.numShuffles
        fprintf('\b\b\b\b%4i  ', iShuffle )
        shuffledLabels = dat.labels(randperm( numTrials ) );
        shuffledPredictions(:,iShuffle) = classifySpeech( dat.neuralFeatures, shuffledLabels, params.svn );
        shuffledAccuracies(iShuffle) = nnz( strcmp( dat.labels, shuffledPredictions(:,iShuffle) ) ) ./ numTrials;
        shuffledConfusedMat(:,:,iShuffle) = confusionmat( dat.labels, shuffledPredictions(:,iShuffle) );
    end
    fprintf( '\nDONE (took %.1f minutes)\n', toc/60)    
    
    % Compute p values for overall classification accuracy based on these shuffles.
    betterThanShuffles = nnz( result.successRate > shuffledAccuracies );
    result.pValueVersusShuffle = 1- (betterThanShuffles / (params.numShuffles+1) );
    fprintf( 'Mean shuffled labels success rate is %.1f%%. Real success rate exceeds this at p = %f\n', ...
        100*mean( shuffledAccuracies ),  result.pValueVersusShuffle )

    % p values along diagonal of confusion matrix. This is effectively whether each sound
    % could be classified better than chance
    for i = 1 : size( shuffledConfusedMat, 1 )
        for j = 1 : size( shuffledConfusedMat, 2 )
            betterThanShuffles = nnz( result.confuseMat(i,j) > shuffledConfusedMat(i,j,:) );
            result.confuseMat_pValueVersusShuffle(i,j) = 1- (betterThanShuffles / (params.numShuffles+1) );
        end
    end
    fprintf( 'Individual label p-values versus shuffles:\n' )
    for iLabel = 1 : numel( result.confuseMatLabels )
        fprintf('  %8s: p=%g\n', ...
            result.confuseMatLabels{iLabel}, result.confuseMat_pValueVersusShuffle(iLabel,iLabel) )
    end
    
else
    fprintf( 'Shuffle controls were skipped. Set params.numShuffles to > 0 to run label shuffled decoding (warning: it''s very slow.)\n' )
end

% Panel C (not generated here) just plots the summary results of each of these datasets.