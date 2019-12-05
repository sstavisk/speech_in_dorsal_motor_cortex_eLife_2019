% Calculates population modulation during speaking different syllables as in Figure 1E and Figure 1-figure supplement 4B of
% Stavisky et al. eLife 2019.
%
% Note 1: By 'population modulation', we mean the less-biased firing rate vector distance between a
% given spoken conditions' firing rates, and the silence conditions. See Methods: Neural
% population modulation.
%
% Note 2: The core neural distance calculation requires single-trial neural data, which we are
% not providing out of respect for our participant's privacy. We instead provide the processed neural distances. 
% The (pre-run) code is provided below, but is not run here. As an upside, this allows this script to run quickly since
% the neural distance calculation is lengthy.
% 
% USAGE: The script can be pointed to several available datasets. By default it will plot
% the participant T5, syllables dataset population modulation (Figure 1E, left panel). To examine other datasets,
% uncomment/modify the dataset selection in the 'Dataset specification' code block below.
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
% 
%   This script was developed and tested in MATLAB_R2017b, but it ought to be robust across
%   versions within a couple years of this.
%
%   The neural distance metric uses the cvDistance.m function available at
%   https://github.com/fwillett/cvVectorStats. However, since this operation isn't run in this script (see
%   above), but rather the distances are already provided, this function isn't necessary
%   here.
%
% This code is associated with "Neural ensemble dynamics in dorsal motor cortex during
% speech in people with paralysis", eLife (2019).
%
% Copyright Sergey D. Stavisky, November 2019, Stanford Neural Prosthetics Translational
% Laboratory.


%% Dataset specification
% % Figure 1E:
dataset = 'T5-syllables_popModulation'; % left panel
% dataset = 'T8-syllables_popModulation'; % right panel

% % Figure 1-figure supplement 4B:
% dataset = 'T5-words_popModulation'; % left panel
% dataset = 'T8-words_popModulation'; % right panel
%% Load the data
load( dataset, 'dat' ); % variable will be called 'dat'
uniqueLabels = dat.uniqueLabels; % conditions that will be plotted


%% Neural distance [skipped over here because single trial data is not sharable]
% Here's the code that would get the neural distance deviation from silence for each speaking condition, in each epoch
% All the single trial data is in a structure 'R' (one structure array element per trial).
% The helper function AlignedMultitrialDataMatrix.m returns a tensor consisting of data from each trial, across all
% channels, for the time window specified.
haveSingleTrialData = false;
if haveSingleTrialData % so code more readable b/c not commented
    for iEvent = 1 : numel( dat.params.alignEvent )
        myTrialInds_silence = strcmp( allLabels, 'silence' );
        silenceDat = AlignedMultitrialDataMatrix( R(myTrialInds_silence), 'featureField', params.neuralFeature, ...
            'startEvent', params.startEvent{iEvent}, 'alignEvent', params.alignEvent{iEvent}, 'endEvent', params.endEvent{iEvent} );
        
        for iLabel = 2 : numel( uniqueLabels ) % compare to silence, which is always uniqueLabels{1}
            myLabel = uniqueLabels{iLabel};
            myTrialInds = strcmp( allLabels, myLabel );
            %  get single trial data tensor
            speakDat = AlignedMultitrialDataMatrix( R(myTrialInds), 'featureField', params.neuralFeature, ...
                'startEvent', params.startEvent{iEvent}, 'alignEvent', params.alignEvent{iEvent}, 'endEvent', params.endEvent{iEvent} );
            
            speakDat.(myLabel).popDistanceFromBaseline{iEvent} = nan( speakDat.numSamples, 1 );
            fprintf('   Unbiased population distances %s (this takes a while)... ', myLabel )
            % loop over time calculating pop unbiased
            fprintf('Sample      1');
            for t = 1 : 1 : min( speakDat.numSamples, silenceDat.numSamples )
                if ~mod( t, 10)
                    fprintf('\b\b\b\b%4i', t)
                end
                myRates = squeeze( speakDat.dat(:,t,:) ); % trials x channels
                mySilence = squeeze( silenceDat.dat(:,t,:) );
                myDistance = cvDistance( myRates, mySilence ); % compare to silence
                dat.(myLabel).popDistanceFromBaseline{iEvent}(t) = myDistance;
            end
            fprintf('\n')
        end
    end
end


% ------------------------------------------
%% Get max and mean distance in each epoch and plot
% ------------------------------------------
% PLOT POP FR DISTANCE FOR EACH LABEL
figh = figure;
figh.Color = 'w';
titlestr = sprintf('Neural distance from silence %s', dataset);
figh.Name = titlestr;
axh = [];

% consistent horizontal axis between panels 
startAt = 0.1;
gapBetween = 0.05;
epochDurations = nan( numel( dat.params.alignEvent ), 1 );
epochStartPosFraction = epochDurations; % where within the figure each subplot starts. 
for iEvent = 1 : numel( dat.params.alignEvent )
    epochDurations(iEvent) = range( dat.(uniqueLabels{1}).t{iEvent} );
end
% I want to fill 0.8 of the figure with both axes, and have a 0.05 gap between subplots,
epochWidthsFraction = (1 - 2*startAt  - gapBetween*(numel( epochDurations ) - 1)) * (epochDurations ./ sum( epochDurations ));
epochStartPosFraction(1) = startAt;
for iEvent = 2 : numel( epochDurations )
    epochStartPosFraction(iEvent) = epochStartPosFraction(iEvent-1) + epochWidthsFraction(iEvent-1) + gapBetween;
end

meanDeviations = nan( numel( uniqueLabels ), numel( dat.params.alignEvent ) ); % label, event
meanDeviationsBaseline = nan( numel( uniqueLabels ), 1 ); % label x 1 (there's also a baseline epoch for event 1 alignment)
for iEvent = 1 : 2
    axh(iEvent) = subplot(1, numel( dat.params.alignEvent ), iEvent); hold on;           
    myPos =  get( axh(iEvent), 'Position');
    set( axh(iEvent), 'Position', [epochStartPosFraction(iEvent) myPos(2) epochWidthsFraction(iEvent) myPos(4)] )
    xlabel( [dat.params.alignEvent{iEvent} ' (s)']);
    set( axh(iEvent), 'TickDir', 'out' );
    
    for iLabel = 2 : numel( uniqueLabels ) % don't include silence
        myLabel = uniqueLabels{iLabel};

        % PLOT THIS POP FR DISTANCE
       myX = dat.(myLabel).t{iEvent};
       myY = dat.(myLabel).popDistanceFromBaseline{iEvent};
       plot( myX, myY, 'Color', labelColors( myLabel ), ...
           'LineWidth', 1 );
        
        % this analysis epoch start and stop
        [~, myStartInd] = FindClosest( dat.(myLabel).t{iEvent}, dat.params.compareWindow{iEvent}(1) );
        [~, myEndInd] = FindClosest( dat.(myLabel).t{iEvent}, dat.params.compareWindow{iEvent}(2) );
        
        % MEAN DEVIATION
        meanDeviations(iLabel,iEvent) = mean( dat.(myLabel).popDistanceFromBaseline{iEvent}(myStartInd:myEndInd) );
        % show it
        lh = line( dat.params.compareWindow{iEvent}, [meanDeviations(iLabel,iEvent) meanDeviations(iLabel,iEvent)], ...
            'LineWidth', 0.5, 'Color', labelColors( myLabel ) );
        
        if iEvent == 1
            % baseline
            % this analysis epoch start and stop
            [~, myStartInd] = FindClosest( dat.(myLabel).t{1}, dat.params.baselineCompareWindow(1) );
            [~, myEndInd] = FindClosest( dat.(myLabel).t{1}, dat.params.baselineCompareWindow(2) );
            meanDeviationsBaseline(iLabel) = mean( dat.(myLabel).popDistanceFromBaseline{1}(myStartInd:myEndInd) );
        end        
    end
end

% plot comparison epochs
linkaxes( axh, 'y' )
myY = get( axh(2), 'YLim' ); % for plotting analysis epoch
myY = myY(2)-0.05*range( myY );
axes( axh(1) )
line( dat.params.compareWindow{1}, [myY myY], 'LineWidth', 1.5, 'Color', [.4 .4 .4] )
line( dat.params.baselineCompareWindow, 0.95.*[myY myY], 'LineWidth', 1.5, 'Color', 'b' ) % baseline epoch
ylabel('Unbiased pop FR distance')
axes( axh(2) )
line( dat.params.compareWindow{2}, [myY myY], 'LineWidth', 1.5, 'Color', 'k' )
    set( axh(iEvent), 'TickDir', 'out' );
ya = get( axh(2), 'YAxis' );
ya.Visible = 'off';

% Report mean population modulation of each epoch, across speaking conditions
meanOfMean1 = nanmean( meanDeviations(:,1) );
meanOfMean2 = nanmean( meanDeviations(:,2) );
[p,h] = signrank( meanDeviations(:,1), meanDeviations(:,2) );

fprintf('MEAN DEVIATION: Mean across conditions 1 = %.3f, mean across conditions 2 = %.3f. ratio = %g, p = %g (sign-rank)\n', ...
    meanOfMean1, meanOfMean2, meanOfMean2/meanOfMean1, p )

% baseline to epoch 1 comparison
meanOfBaseline = nanmean( meanDeviationsBaseline (:,1) );
[p,h] = signrank(  meanDeviations(:,1), meanDeviationsBaseline );
fprintf('MEAN DEVIATION: Mean baseline epoch = %.3f, mean across conditions 1 = %.3f. p = %g (sign-rank)\n', ...
    meanOfBaseline, meanOfMean1 , p )
