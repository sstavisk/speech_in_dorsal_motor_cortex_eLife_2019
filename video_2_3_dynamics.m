% Plots neural state trajectories as they unfold over the course of a speaking trial, as
% in Videos 2 and 3 of Stavisky et al. eLife 2019.
% 
% USAGE: The script can be pointed to several available datasets and either Go or Speak (voice on)
% alignment.  By default it will plot the T5-words dataset, go-aligned, as in Video 2. To examine other neural 
% recordings or alignments, uncomment/modify the dataset in the 'Dataset specification' code block below. 
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
%   This script also uses the jPCA code pack associated with Churchland*, Cunningham*, et al. 
%   "Neural population dynamics during reaching", Nature 2012
%   (https://www.nature.com/articles/nature11129). You can download this from
%   https://churchland.zuckermaninstitute.columbia.edu/content/code.
%   Note that this version's pca.m uses the deprecated princomp function, which causes so
%   many warnings that it dramatically slows execution. I've included a jPCA_updated.m in
%   /helper_funcitons/ that swaps out princomp for pca to avoid this. 
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
% Aligned to go cue
dataset = 'T5-words_video_go.mat'; % Video 2 in the paper
% dataset = 'T5-syllables_video_go.mat';
% dataset = 'T8-words_video_go.mat'; %  not much to look at
% dataset = 'T8-syllables_video_go.mat'; % not much to look at

% Aligned to voice on
% dataset = 'T5-words_video_speak.mat'; % Video 3 in the paper
% dataset = 'T5-syllables_video_speak.mat'; 
% dataset = 'T8-words_video_speak.mat'; 
% dataset = 'T8-syllables_video_speak.mat'; 


%% Load the data
load( dataset ); % key variable loaded is called 'dat'
labels = dat.labels; % shortens for frequent use

% Color scheme
viz.colorBaseline = [.4 .4 .4]; % trajectory color before the audio cue comes on (GRAY).
viz.colorDelay = [0 0 1]; % trajectory color from audio cue until go cue (BLUE).
viz.LineWidthBaseline = 2;
viz.LineWidthDelay = 3;
viz.LineWidthMove = 4;
viz.FontSize = 30;

% Video 
viz.msEachFrame = 10; % time resolution of video (in ms). Can by any multiple of 1 ms
viz.figurePixels = [720 720];

% Camera rotation scripting
viz.cam.view0 = [45, 30]; % starting view (azimuth, elevation)



%% Format for dPCA
% firingRatesAverage: N x S x D 
% Note: operates on the PSTHs from the specified dPCA analysis epoch
% 

% T is the number of time-points (note that all the trials/conditions should have the
% same length in time!)
N = size( dat.(labels{1}).psthCis, 2 ); % N is the number of neurons
S = numel( labels ); % S is the number of conditions (e.g., words)
T = numel( dat.(labels{1}).tCIS ); % T is number of samples (ms)

featureAverages = nan(N, S, T);
for iLabel = 1 : numel( labels )
    featureAverages(:,iLabel,:) = dat.(labels{iLabel}).psthCis';
end

% Soft-norm is calculated here (and stored in allChannelRange)
allChannelRange =  max( max( featureAverages, [], 3 ), [], 2 ) - min( min( featureAverages, [], 3 ), [], 2 );     % range for each channel
allChannelRange = allChannelRange + dat.params.softenNorm; % this should also be denominator for any soft-norm I want to repeat later
featureAverages = featureAverages .* repmat( 1./allChannelRange, 1, size( featureAverages, 2 ), size( featureAverages, 3 ) );


%% 1. Do dPCA on the go epoch
% single-trial derived optimalLambda and Cnoise were pre-computed, as in
% figure4_conditionInvariantDynamics.m
fprintf('dPCA on data from %s to %s\n', dat.params.CISstartEvent, dat.params.CISendEvent );
[W,V,whichMarg] = dpca(featureAverages, 8, ...
        'combinedParams', dat.combinedParams, ...
        'lambda', dat.optimalLambda, ...
        'Cnoise', dat.Cnoise); % num dims won't matter since we're only using dPC1 which is always the CIS in these data

explVar = dpca_explainedVariance(featureAverages, W, V, ...
    'combinedParams', dat.combinedParams);
numCDdims = nnz( whichMarg == 1 );
numCIdims = nnz( whichMarg == 2 );   
% re-sort based on how much condition-independent activity there is
eachVarExplained = explVar.margVar;
CIdims = whichMarg==2;
CIdimsInds = find( CIdims );
[~, sortIndsByCI] = sort( explVar.margVar(2,CIdims), 'descend');
% also grab all the W and V vectors with this order
Wreordered = [W(:,CIdimsInds(sortIndsByCI)), W(:,~CIdims)];
Vreordered = [V(:,CIdimsInds(sortIndsByCI)), V(:,~CIdims)]; % this is what we really were after


% Flip CIS1 if necessary to keep positive = more CIS activation (it's AU and thus arbitrary direction anyway, but
% this keeps it consistently left-to-right in videos).
acrossCondsMean = squeeze( mean( featureAverages,2 ) )'; % T x E
GM = acrossCondsMean * Wreordered(:,1);
if GM(end) < GM(1)
    Wreordered(:,1) = -Wreordered(:,1);    
    Vreordered(:,1) = -Vreordered(:,1);    
end

CIS1dim = Wreordered(:,1); % this too


%% 2. Do jPCA on the speak epoch
fprintf('jPCA on data from %s to %s\n', dat.params.jPCAstartEvent, dat.params.jPCAendEvent );
% apply softnorm to jPCA psth; keeps things consistent
for iLabel = 1 : numel( dat.labels )
    myLabel = dat.labels{iLabel};
    result.(myLabel).psthJPCA_softNormed = dat.(myLabel).psthJPCA .* repmat( 1./allChannelRange', size( dat.(myLabel).psthJPCA, 1 ) , 1 );
end

jPCA_params.normalize = false;  % already done to data
jPCA_params.meanSubtract =  1;
jPCA_params.softenNorm = nan;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;
jPCA_params.numPCs = 6;

% Format the data for jPCA: subsample every X ms (Mark's code won't subsample itself)
startInd = find( round( 1000.*dat.(myLabel).tjPCA ) == dat.params.dataTimestampsInMS(1), 1, 'first' ); % doesn't matter what myLabel is at this point, all have same times
endInd = find( round( 1000.*dat.(myLabel).tjPCA ) == dat.params.dataTimestampsInMS(end), 1, 'first' );

% A critical operation here is to project the data into the null space of CIS1 found by the dPCA
% operation above.
nspace = null( Vreordered(:,1)' ); 
clear('Data')
fprintf('Projecting data into %i-dim null space of CIS1, then keeping %i PCs \n', size( nspace, 2 ), jPCA_params.numPCs )
for iLabel = 1 : numel( dat.labels )
     myLabel = dat.labels{iLabel};
     Data(iLabel).A = result.(myLabel).psthJPCA_softNormed(startInd:dat.params.downSampleEveryNms:endInd,:) * nspace; % time x dPC
     Data(iLabel).times = round( 1000.* dat.(myLabel).tjPCA(startInd:dat.params.downSampleEveryNms:endInd) )';
     Data(iLabel).label = myLabel; % not used but nice to keep track of just in case.
end

[Projection, Summary] = jPCA( Data, dat.params.dataTimestampsInMS, jPCA_params ); % do the actual jPCA operation
[ indsLeftToRight, cmapRG ] = whichGroupIsWhichJpca( Projection ); % get red-green colormap

% projection matrix that does both the dPCA projection and then the jPCA projection
dim2Proj = nspace*Summary.jPCs_highD(:,1);
dim3Proj = nspace*Summary.jPCs_highD(:,2);


%% 3. Get ready to plot each condition's CIS_1 vs jPC1 vs jPC2 by projecting all trial data into this space
% apply softnorm to video psth
allVidXYZ = nan( size( dat.(myLabel).psthVid, 1 ), 3, numel( dat.labels ) ); % Time x 3 x (# Conditions)
for iLabel = 1 : numel( dat.labels )
    myLabel = dat.labels{iLabel};
    result.(myLabel).psthVid_softNormed = dat.(myLabel).psthVid .* repmat( 1./allChannelRange', size( dat.(myLabel).psthVid, 1 ) , 1 );
    myXYZ = [result.(myLabel).psthVid_softNormed * CIS1dim, ...
        result.(myLabel).psthVid_softNormed * dim2Proj, ...
        result.(myLabel).psthVid_softNormed * dim3Proj ];
    allVidXYZ(:,:,iLabel) = myXYZ;
end

% subtract across-conditions mean (standard for jPCA)
allVidXYZ(:,2,:) = allVidXYZ(:,2,:) - repmat( mean( allVidXYZ(:,2,:), 3 ), 1, 1, numel( dat.labels ) );
allVidXYZ(:,3,:) = allVidXYZ(:,3,:) - repmat( mean( allVidXYZ(:,3,:), 3 ), 1, 1, numel( dat.labels ) );
   
tVid = dat.(myLabel).tVid; % keep on hand and easy to refer to
% If this is go-aligned, replicate the trial events across conditions so that below code
% works for both go-aligned (shared events times) and voice-aligned (each condition has its own event times)
if numel( dat.sampleAudioCue ) == 1
    dat.sampleAudioCue = repmat( dat.sampleAudioCue, 1, numel( dat.labels ) );
    dat.sampleGoCue = repmat( dat.sampleGoCue, 1, numel( dat.labels ) );
end

%% Draw the axes
figh = figure;
figh.Name = sprintf('Neural state video %s', dataset );
figh.Units = 'pixels';
figh.Position(3) = viz.figurePixels(1);
figh.Position(4) = viz.figurePixels(2);

figh.Color = 'w';
axh = axes;
axis tight

hold on; 
xlabel('CIS_1', 'FontSize', viz.FontSize)
ylabel('jPC_1', 'FontSize', viz.FontSize) 
zlabel('jPC_2', 'FontSize', viz.FontSize)


% compute the extrema coordinate values so I can pre-set axis dimensions
axh.XLim = [min( min( allVidXYZ(:,1,:)  ) ), max( max( allVidXYZ(:,1,:)  ) )];
axh.YLim = [min( min( allVidXYZ(:,2,:)  ) ), max( max( allVidXYZ(:,2,:)  ) )];
axh.ZLim = [min( min( allVidXYZ(:,3,:)  ) ), max( max( allVidXYZ(:,3,:)  ) )];
view( viz.cam.view0 );

axis vis3d % locks aspect ratio so it doesn't change when rotating
totFrames = floor( numel( tVid ) / viz.msEachFrame );
axPos = axh.Position;

% Put the time text into the figure corner
axTime = axes( figh, 'OuterPosition', [0.05, 0.8, 0.1, 0.1], 'Color', 'r'); % allows for closer cropping after

axTime.Visible = 'off';
hTime = text( 0, 1, sprintf('%s%+.1fs', dat.params.vidAlignEvent, tVid(viz.msEachFrame) ), 'FontSize', viz.FontSize, ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'Parent', axTime);
axes( axh )

for iFrame = 1 : totFrames 
    mySample = iFrame * viz.msEachFrame;
    myT = tVid(mySample);
        hTime.String = sprintf('%s%+.1fs', dat.params.vidAlignEvent, myT );
    try
        delete( lh )
    catch
    end
    lh = [];
    for iLabel = 1 : numel( dat.labels )                    
        % Plot baseline (before prompt)
        endBaseline = min( dat.sampleAudioCue(iLabel)-1, mySample ); % different for each label
        x = allVidXYZ(1:endBaseline,1,iLabel);
        y = allVidXYZ(1:endBaseline,2,iLabel);
        z = allVidXYZ(1:endBaseline,3,iLabel);
        lh(end+1) = plot3( x, y, z, 'Color', viz.colorBaseline, 'LineWidth', viz.LineWidthBaseline );
        
        % Plot delay (prompt until go cue)        
        if mySample >= dat.sampleAudioCue(iLabel)
            endDelay = min( dat.sampleGoCue(iLabel)-1, mySample );
            x = allVidXYZ(dat.sampleAudioCue(iLabel):endDelay,1,iLabel);
            y = allVidXYZ(dat.sampleAudioCue(iLabel):endDelay,2,iLabel);
            z = allVidXYZ(dat.sampleAudioCue(iLabel):endDelay,3,iLabel);
            lh(end+1) = plot3( x, y, z, 'Color', viz.colorDelay, 'LineWidth', viz.LineWidthDelay );
        end

        % Plot post go
        if mySample >= dat.sampleGoCue(iLabel)
            endSpeech = min( numel( tVid ), mySample );
            x = allVidXYZ(dat.sampleGoCue(iLabel):endSpeech,1,iLabel);
            y = allVidXYZ(dat.sampleGoCue(iLabel):endSpeech,2,iLabel);
            z = allVidXYZ(dat.sampleGoCue(iLabel):endSpeech,3,iLabel);            
            lh(end+1) = plot3( x, y, z, 'Color', cmapRG(iLabel,:), 'LineWidth', viz.LineWidthMove ); % replace with Red <--> Green
        end 
    end

    drawnow;
end

