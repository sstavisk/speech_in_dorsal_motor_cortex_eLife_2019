% Analyzes neural modulation aligned to breathing during as in Figure 2-figure supplement 2 of
% Stavisky et al. eLife 2019.
%
% Note: The breath-triggered firing rate calculation requires single-trial neural data, which we are
% not providing out of respect for our participant's privacy. We instead provide the event-averaged firing rates 
% already aligned to breaths. We also pre-calculate the shuffle controls in which the same procedure was repeated
% with randomized breath times, after which these faux-breath aligned firing rates' max and min values
% are stored. These shuffles took a long time to run, so having those pre-calculated allows this
% script to be executed in a reasonable time.
%
% USAGE: The script can be pointed to several available datasets. By default it will plot
% the participant T5, unattended breathing-triggered neural modulation. To examine instructed breathing,
% uncomment/modify the dataset selection in the 'Dataset specification' code block below.
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
% 
%   This script was developed and tested in MATLAB_R2017b, but it ought to be robust across
%   versions within a couple years of this.
%
%   Script uses findpeaks.m which is in MATLAB's Signal Processing Toolbox. If you don't have
%   it, you can skip the breath finding portion of the script and move on to the neural
%   analyses.
%
% This code is associated with "Neural ensemble dynamics in dorsal motor cortex during
% speech in people with paralysis", eLife (2019).
%
% Copyright Sergey D. Stavisky, November 2019, Stanford Neural Prosthetics Translational
% Laboratory.


%% Dataset specification
% % Figure 1E:
dataset = 'T5-breathing_unattended'; % such as in panel B
% dataset = 'T5-breathing_instructed'; % such as in panel B

params.useSorted  = false; % Threshold crossings, like in panel C
% params.useSorted  = true; % Sorted single units, like in panel D.


%% Analysis parameters
params.pvalue = 0.01; % plot chance distribution boundary at this p value 

viz.numBreathsToPlot = 200; % downsample to make it manageable

% Breath detection parameters
% Note: Provided breath aligned neural data assume the parameters used below.
params.bigPeakMinDistance = 5; % in seconds, used to find 'big' breaths to calibrate to scale of data
params.minPeakDistance = 1; % in seconds
params.minPeakProminenceFractionOfBigPeak = 0.3; % look for breaths at least this big compared to median 'big' breaths

% for plotting breath-aligned stretch sensor and neural data. 
params.breath.secsBeforeBreath = 2;
params.breath.secsAfterBreath = 1.5;
params.psthSmoothingGaussianSD = 25; % in milliseconds. Matches what was used for firing rates during speaking.

% Shuffle test
params.numShuffles = 1001;


%% dataset-dependent parameters
if contains( dataset, 'instructed' )
    params.clipSecondsStart = 30; % avoids start and end of the block where breath meter is recording but task isn't running yet
    params.clipSecondsEnd = 15;
    viz.plotExampleBeltBlock = 9; % Matches Panel A, right
    viz.plotExampleBeltStart = 41.5;
    viz.plotExampleBeltEnd = 79.2;
else
    params.clipSecondsStart = 0; % can use whole block since there's no breathing task.
    params.clipSecondsEnd = 0;
    viz.plotExampleBeltBlock = 7; % Matches Panel A, left
    viz.plotExampleBeltStart = 39; 
    viz.plotExampleBeltEnd = 79.2;
end
% Examples will match panels C and D if TCs unattended dataset is used
if params.useSorted 
    viz.plotExampleChans = [6, 2, 9]; 
    chanStr = 'units';
else
    viz.plotExampleChans = [23, 186, 72 ]; %  TCs. Chan 186 is '2_90' from the figure; 
    chanStr = 'electrodes';
end

arrayMaps ={'T5_lateral', 'T5_medial'};


% --------------------------------------------------------------------------------------------------------
% BREATH MEASUREMENTS: Panels A and B
% -------------------------------------------------------------------------------------------------------

%% Find breath times
% Load breath (stretch belt) data. This has fields that are cell lists, one element per block. For each block, 
% .breathSensor_t  is  a time series of time stamps, in seconds, relative to the block start. 
% The same clock is used for the neural data. Field .breathSensor is the corresponding stretch belt measurement.
% Note: These data are already slightly pre-processed, otherwise the raw 30,000Hz data
% would be multiple gigabytes. Specifically, the data were low-pass filtered below 3 Hz
% (6th order, 0.2 passband ripple, 50 dB stop band attenuation, two-way filtfilt. These
% data were then decimated to 1000 Hz. 
breath = load( sprintf('%s_breath.mat', dataset ) );
numBlocks = numel( breath.dat.breathSensor );
FsBreath = round( 1/(breath.dat.breathSensor_t{1}(2)-breath.dat.breathSensor_t{1}(1)) ); % sampling rate of breath

breathTimesEachBlock = {}; % Breath times during each block, in seconds
allBreathTraces = []; % will fill with breathing data across all the blocks of this condition

bufferS = 0.0770; % smoothing width in seconds for neural data; used here to not get a breath where there isn't complete corresponding neural data
% Populate breath-aligned sensor times for each block
for iBlock = 1 : numBlocks
    b = breath.dat.breathSensor{iBlock}; % this block's breath
    % 1. Find 'big peaks' (with large time between them) which I use to calibrate what a large prominence is
    [pks, locs, widths, proms] = findpeaks( b, 'MinPeakDistance', params.bigPeakMinDistance*FsBreath );
    medianBigProminence = median( proms );

    % 2. Find peaks using specified min distance using the specified fraction of big prominence
    [pks, locs, widths, proms] = findpeaks( b, 'MinPeakDistance', params.minPeakDistance*FsBreath,  ...
        'MinPeakProminence', params.minPeakProminenceFractionOfBigPeak*medianBigProminence   );  
    bTimesThisBlock = locs;

    
    % 3. Don't accept breaths that happen too early or late in a block, since we won't be able to get the
    % neural data.
    startCutoff = max( (params.breath.secsBeforeBreath+bufferS), params.clipSecondsStart )*FsBreath;   % extra buffer for the Gaussian smoothing
    tooEarlyBreaths = find( bTimesThisBlock < startCutoff );
    bTimesThisBlock(tooEarlyBreaths) = [];
    endCutoff =  numel( b ) - FsBreath * max( params.breath.secsAfterBreath+bufferS, params.clipSecondsEnd );
    tooLateBreaths = find( bTimesThisBlock > (numel( b ) - (params.breath.secsAfterBreath+bufferS)*FsBreath) );
    bTimesThisBlock(tooLateBreaths) = [];
    breathTimesEachBlock{iBlock} = bTimesThisBlock ./ FsBreath; % in seconds
  
    % PLOT EXAMPLE RESPIRATORY BELT MEASUREMENTS WITH DETECTED BREATHS
    if iBlock == viz.plotExampleBeltBlock
        figh = figure;
        figh.Name = sprintf( 'Respiratory belt %s block %i', dataset, iBlock );
        title( figh.Name, 'Interpreter', 'none' );
        hold on;

        % Plot the belt measurements
        plot( breath.dat.breathSensor_t{iBlock}(startCutoff:endCutoff), b(startCutoff:endCutoff), ...
            'k', 'LineWidth', 2 )
        xlim( [viz.plotExampleBeltStart viz.plotExampleBeltEnd] )
        xlabel('Time (s)')
        ylabel('Respiratory belt (AU)')

        % Plot the detected belts
        sh = scatter( bTimesThisBlock/FsBreath, b(bTimesThisBlock), 'filled', 'MarkerFaceColor', [1 0 1], 'SizeData', 100 );
    end
    
    % 4. Save breath peak-aligned stretch belt sensor measurements (for generating the example/average
    % breath traces in panel B).
    myBreaths = nan( numel( bTimesThisBlock ), numel( -params.breath.secsBeforeBreath :1/FsBreath:params.breath.secsAfterBreath  ) ); % # breaths x time
    for i = 1 : numel( bTimesThisBlock )
        mySamples = bTimesThisBlock(i) - (params.breath.secsBeforeBreath*FsBreath) : FsBreath/1000 : bTimesThisBlock(i) + (params.breath.secsAfterBreath*FsBreath);
        myBreaths(i,:) = b(mySamples);
    end
    allBreathTraces = [allBreathTraces; myBreaths];
end

% get average breath
meanBreathTrace = mean( allBreathTraces, 1 );


%% Plot the breath peak-triggered stretch belt sensor measurements (panel B)
figh = figure;
figh.Name = sprintf( 'Aligned breaths %s', dataset );
hold on;
title( figh.Name, 'Interpreter', 'none' );
breathT = -params.breath.secsBeforeBreath : 0.001 : params.breath.secsAfterBreath;
totalNumberBreaths = size( allBreathTraces, 1 );
fprintf('Plotting %i/%i breaths\n', viz.numBreathsToPlot, totalNumberBreaths )
% pick breaths
plotTheseBreaths = indexEvenSpaced( totalNumberBreaths, viz.numBreathsToPlot);
plot( breathT, allBreathTraces(plotTheseBreaths,:), 'Color', [.7 .7 .7], 'LineWidth', 0.5 );
% plot the mean
plot( breathT, meanBreathTrace, 'Color', 'k', 'LineWidth', 2 )
xlabel('Time after breath peak (s)');
ylabel('Respiratory belt (AU)')
axh = figh.Children;
axh.TickDir = 'out';

% -------------------------------------------------------------------------------------------------------
%% NEURAL MEASUREMENTS: Panels C/D and E
% -------------------------------------------------------------------------------------------------------
% Load the neural data
if params.useSorted
    neuralFile = sprintf('%s_neural_sortedUnits.mat', dataset );
else
    neuralFile = sprintf('%s_neural_TCs.mat', dataset );
end
neural = load( neuralFile );
Nchans = size( neural.dat.trialAverageFR, 2 );
Nsamples = size( neural.dat.trialAverageFR, 1 );

% Across-channels firing rate plotted in panel G. For convenience later I load these pre-computed for all conditions and plot
grandMean = mean( neural.dat.trialAverageFR, 2 ); 

% Calculate shuffled breath-aligned max/min firing rate confidence bound for each channel
upperBoundInd = round( (1 - 0.5*2*params.pvalue) * neural.dat.numShuffles ); % multiply by 2 since two-sided test
lowerBoundInd = round( (0.5*2*params.pvalue) * neural.dat.numShuffles );

timeBoundsEachChan = nan( Nsamples, Nchans, 2 ); % will be time x chan x 2 where last dimension is (lower, upper)
significantChannels = []; % which channels significantly break above or below the chance bound.
timeBoundsEachChan(:,:,2) = repmat( neural.dat.maxFRshuffle(upperBoundInd,:), Nsamples, 1 )  ;
timeBoundsEachChan(:,:,1) = repmat( neural.dat.minFRshuffle(lowerBoundInd,:), Nsamples, 1 )  ;
for iChan = 1 : Nchans
    highCrossings = find( neural.dat.trialAverageFR(:,iChan) > timeBoundsEachChan(:,iChan,2) );
    lowCrossings = find( neural.dat.trialAverageFR(:,iChan) < timeBoundsEachChan(:,iChan,1) );
    if ~isempty( highCrossings ) || ~isempty( lowCrossings )
        significantChannels(end+1) = iChan;
    end
end
significantChannelsChanNumber = neural.dat.channelNumbers(significantChannels); % absolute electrode number, instead of indices into live channels
fprintf('%i/%i (%.1f%%) %s %s significantly modulate at p=%g (vs. %i shuffles)\n', ...
    numel( significantChannelsChanNumber ), Nchans, 100*numel( significantChannelsChanNumber )/Nchans, ...
   dataset, chanStr, params.pvalue, neural.dat.numShuffles ); 
% Note: the # tuned quoted in the panel C caption is the union of
% significantChannelsChanNumber when run for instructed and unattended datasets
% For convenenience later I load these pre-computed for both conditions and report this
% union.


%% Modulation depth
modulationDepth = nan( Nchans, 1 ); % max - minimum in the analysis window
for iChan = 1 : Nchans
    modulationDepth(iChan) = max( neural.dat.trialAverageFR(:,iChan) ) - min( neural.dat.trialAverageFR(:,iChan) );
end
fprintf('Mean %s modulation depth = %.3f Hz, median = %.3fHz\n', dataset, nanmean( modulationDepth ), nanmedian( modulationDepth ) )


%% Plot example channels
figh = figure;
figh.Name = sprintf('Example channels %s', dataset );
for iChan = 1 : numel( viz.plotExampleChans )
    myChan = find( neural.dat.channelNumbers == viz.plotExampleChans(iChan)); % index into live channels
    axh(iChan) = subplot( numel( viz.plotExampleChans ), 1, iChan );
    if params.useSorted
        myName = sprintf('neuron T5_%i.%i', neural.dat.sortInfo.unitArray(myChan), neural.dat.sortInfo.unitChannel(myChan) );
    else
        myName = sprintf('elec T5_%i', neural.dat.channelNumbers(myChan) );
    end
    hold on;
    
    % plot shuffles bounds
    plot( neural.dat.firingRate_t, squeeze( timeBoundsEachChan(:, myChan, : ) ), 'Color', [.5 .5 .5], 'LineWidth', 0.5 )
    
    % plot real data
    plot( neural.dat.firingRate_t, neural.dat.trialAverageFR(:, myChan ) - neural.dat.SEM_FR(:, myChan ), 'Color', 'k', 'LineWidth', 0.5 );
    plot( neural.dat.firingRate_t, neural.dat.trialAverageFR(:, myChan ) + neural.dat.SEM_FR(:, myChan ), 'Color', 'k', 'LineWidth', 0.5 );
    plot( neural.dat.firingRate_t, neural.dat.trialAverageFR(:, myChan ), 'Color', 'k', 'LineWidth', 2 );
    title( sprintf('%s FR, MD = %.2fHz', myName, modulationDepth(myChan) ), 'Interpreter', 'none' );
    xlim( [neural.dat.firingRate_t(1) neural.dat.firingRate_t(end) ]);
end
    

%% Modulation depth overhead plot
chanMap = channelAnatomyMap( arrayMaps );
figh = figure;
figh.Name = sprintf('%s array tuning modulation', dataset );
if params.useSorted
    figh.Name = [figh.Name ' SUA'];
else
    figh.Name = [figh.Name ' TCs'];
end
title( figh.Name, 'Interpreter', 'none' )
graymap = flipud( bone( 256 ) ); 
% start it at a light gray
graymap = graymap(26:end,:);
drawnAlready = []; % will track which electrodes were drawn as having something on them
axh = figh.Children;
hold on;
axh.XLim = chanMap.xlim;
axh.YLim = chanMap.ylim;
axis equal
disabledSize = 4;
disabledColor = graymap(1,:);

if params.useSorted
   unitElec = (neural.dat.sortInfo.unitArray-1)*96 + double( neural.dat.sortInfo.unitChannel );
end

% need to know the values range to scale the color map.
MDminMax = [floor( min( modulationDepth ) ) ceil( max( modulationDepth ) )];
for iElec = 1 : 192
    myElecNum = iElec;
    
     % Data to plot.
     if params.useSorted &&  ismember( iElec, unitElec ) % there's at least one sorted unit on this electrode
         myInd = find( unitElec == iElec ); % index into modulation depths
         myDat =  modulationDepth(myInd);
         if numel( myDat ) > 1
             fprintf('There are multiple units on elec %i. Plotting the maximum modulating one.\n', ...
                 iElec );
             myDat = max( myDat );
         end
         myColor =  graymap( floor( size(graymap,1)*(myDat-MDminMax(1))/range( MDminMax ) )+1 ,:);
         mySize = 36;
     elseif ~ params.useSorted && ismember( iElec, neural.dat.channelNumbers ) % This is an active electrode
            myInd = find( neural.dat.channelNumbers == iElec ); % index into modulation depths
            myDat = modulationDepth(myInd);
            myColor =  graymap( floor( size(graymap,1)*(myDat-MDminMax(1))/range( MDminMax ) )+1 ,:);
            mySize = 36;
     else
         % non-functional electrode
         mySize = disabledSize;
         myColor = disabledColor;
     end
   
    x = chanMap.x(myElecNum);
    y = chanMap.y(myElecNum);
    
    % draw the point
    scatter( x, y, mySize, myColor, 'filled' )
end
colormap( graymap );
cmaph =colorbar;
for i = 1 : numel( cmaph.TickLabels )
    cmaph.TickLabels{i} = [];
end
cmaph.TickLabels{1} = MDminMax(1);
cmaph.TickLabels{end} = MDminMax(end);
cmaph.Label.String = 'Modulation Depth (Hz)';
xlabel( 'Lateral      (mm)      Medial' )
ylabel( 'Posterior      (mm)      Anterior' )
 

% --------------------------------------------------------------------------------------------------------
%% UNATTENDED VS INSTRUCTED BREATHING VS SPEAKING COMPARISONS: Panels F and G
% -------------------------------------------------------------------------------------------------------
% The breathing results used here could be regenerated using the above script (run for
% both unattended breathing, and for instructed breathing). The speech results are
% precomputed, but are essentially the same as what's done in figure1_and_2_firing_rates.m,
% except that the analysis epoch is AO-2.5 to AO+1, and mean firing rates are calculated across all
% non-silent trials rather than within syllable conditions. Modulation depths are then
% taken as a given channel's max - min trial-averaged firing rate across that epoch.

if ~params.useSorted 
    % These analysis are only for TCs (with fewer sorted units and not knowing if they're the same across days,
    % it's less meaningful to compare this small pool of neurons' activities across these
    % conditions)
    
    % Load the summary results
    summary = load( 'T5-breathing_comparisons.mat' );
    
    % Report number of significantly modulating electrodes/units across both breathing
    % conditions
    fprintf('Breathing: %i/%i channels have significant breathing modulation at p=%g.\n', ...
        numel( union( summary.unattended.significantChannels , summary.instructed.significantChannels ) ), ...
        numel( summary.unattended.channelNumbers ), summary.instructed.params.pvalue );
    fprintf('Breathing: %i sorted units have significant breathing modulation at p=%g.\n', ...
        numel( union( summary.unattended.significantUnits , summary.instructed.significantUnits ) ), ...
        summary.instructed.params.pvalue );
    
    
    %% Histograms (Figure 2-figure supplement 2F)
    % Unattended breathing
    figh = figure;
    figh.Name = 'Breathing and speaking modulation depth histograms';
    h1 = histogram( summary.unattended.modulationDepth );
    h1.EdgeColor = 'none';
    h1.FaceColor = [0 0 0];
    xlabel('Modulation depth (Hz)');
    ylabel('# electrodes' );
    hold on;
    fprintf('Median modulation for unattended breathing = %.2fHz\n', median( summary.unattended.modulationDepth ) );
    
    % Instructed breathing
    h2= histogram( summary.instructed.modulationDepth );
    h2.EdgeColor = 'none';
    h2.FaceColor = [0 .3 .7];
    h2.BinWidth = h1.BinWidth;
    fprintf('Median modulation for instructed breathing = %.2fHz\n', median( summary.instructed.modulationDepth ) );
    
    
    % Speaking
    meanAcrossLabelsModDepths = mean( summary.speaking.modulationDepth, 2 );
    fprintf('Median modulation for speaking  = %.2fHz\n', median( meanAcrossLabelsModDepths ) );
    
    h3 = histogram( meanAcrossLabelsModDepths );
    h3.BinWidth = h1.BinWidth;
    h3.EdgeColor = 'none';
    h3.FaceColor = 'r';
    
    axh = gca;
    axh.TickDir = 'out';
    axh.Box = 'off';
    
    line( [nanmedian( summary.unattended.modulationDepth ) nanmedian( summary.unattended.modulationDepth)], [0 max(h1.Values)+1], 'Color', [.1 .1 .1] );
    line( [nanmedian( summary.instructed.modulationDepth ) nanmedian( summary.instructed.modulationDepth )], [0 max(h1.Values)+1], 'Color', [0 0 1] );
    line( [nanmedian( meanAcrossLabelsModDepths ) nanmedian( meanAcrossLabelsModDepths )], [0 max(h1.Values)+1], 'Color', [0.9 0 0] );
    legend({'Unattended', 'Instructed', 'Speaking'})
    
    % Compare all three distributions
    [p,h] = ranksum( summary.unattended.modulationDepth , summary.instructed.modulationDepth );
    fprintf('Unattended breathing vs. instructed distributions rank-sum test p = %g\n', p );
    [p,h] = ranksum( summary.unattended.modulationDepth , meanAcrossLabelsModDepths );
    fprintf('Unattended breathing vs. speaking distributions rank-sum test p = %g\n', p );
    [p,h] = ranksum( summary.instructed.modulationDepth , meanAcrossLabelsModDepths );
    fprintf('Instructed breathing vs. speaking distributions rank-sum test p = %g\n', p );
    
    %% Grand mean firing rates (Figure 2-figure supplement 2G)
    figh = figure;
    figh.Color = 'w';
    figh.Name = 'Breathing and speaking grand mean firing rates';
    
    axh(1) = subplot( 1, 2, 1 );
    hold on
    % Unattended
    plot( summary.unattended.t,  summary.unattended.grandMeanFR , 'Color', 'k', 'LineWidth', 2 )
    xlabel('Time after breath peak (s)');
    ylabel('\Delta population firing rate (Hz)')
    axh(1).TickDir = 'out'; axh(1).Box = 'off';
    
    % Instructed
    % align plot bottoms (accounts for inter-day firing rate shifts, focuses on depth of
    % modulation).
    plot( summary.instructed.t,  summary.instructed.grandMeanFR, 'Color', 'b', 'LineWidth', 2 )
    
    % Speaking
    axh(2) = subplot( 1, 2, 2 );
    myY = summary.speaking.grandMeanFR - (summary.speaking.grandMeanFR(1)-summary.unattended.grandMeanFR(1));
    plot( summary.speaking.t, myY, 'Color', 'r', 'LineWidth', 2 );
    axh(2).TickDir = 'out'; axh(2).Box = 'off';
    axh(2).YAxis.Visible = 'off';
    xlabel( sprintf('Time after %s', summary.speaking.alignEvent ) );
    
    linkaxes( axh, 'y' )
end
