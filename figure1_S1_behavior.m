% Plots acoustic spectrograms and reaction time distributions as in Figure 1-figure supplement 1 of
% Stavisky et al. eLife 2019.
%
% We do not provide raw audio (besides the one set of syllables examples in  Video 1)  out of respect for our participant's privacy.
% We instead provide the processed acoustic spectrograms. The (pre-run) code to generate this is provided below, but is not run here.
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
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
datasets_spectrograms = {...
    'T5-syllables_audioSpectrogram'; % Panel A, left
    'T8-syllables_audioSpectrogram'; % Panel A, right
    'T5-words_audioSpectrogram'; % Panel B, left
    'T8-words_audioSpectrogram'; % Panel B, right
    };

% Panel C
datasets_RTs = {...
    'T5-syllables_RTs';
    'T8-syllables_RTs';
    'T5-words_RTs';
    'T8-words_RTs';
    };

% frequency range to plot
spectrogramFreqMinMax = [100 10000];

%% Calculate acoustic spectrogram [skipped over here because single trial data is not sharable]
haveSingleTrialData = false;
if haveSingleTrialData % so code more readable b/c not commented
    for iLabel = 1 : numel( uniqueLabels )
        myLabel = uniqueLabels{iLabel};
        myTrialInds = strcmp( allLabels, uniqueLabels{iLabel} );
        
        datTensor =  AlignedMultitrialDataMatrix( R(myTrialInds), 'featureField', 'audio', ...
            'startEvent', dat.params.startEvent, 'alignEvent', dat.params.alignEvent, 'endEvent', dat.params.endEvent  );
        segmentlen = dat.params.segmentLength*audioFs;
        noverlap = dat.params.segmentOverlap*audioFs;
        NFFT = dat.params.NFFT;
        
        for iTrial = 1 : jenga.numTrials
            [s,w,t,ps] = spectrogram(jenga.dat(iTrial,:),segmentlen,noverlap,NFFT,audioFs,'yaxis','power');
            dat.(myLabel).ps(iTrial,:,:) = reshape( ps, 1, size(ps,1), size(ps,2) ); % trial x frequency x time bin
        end
        % trial-average for this label
        dat.(myLabel).psAvg = squeeze( mean( dat.(myLabel).ps, 1 ) ); % frequency x time
        % convert these faux-timestamps into actual time relative to trial alignment
        t = t + datTensor.t(1);
        
        % save spectrogram details (just once is fine, but it'll end up repeated )
        dat.spectrogram.s = s;
        dat.spectrogram.w = w;
        dat.spectrogram.t = t;
    end
end

% ------------------------------------------
%% Acoustic spectrograms
% ------------------------------------------
figh = figure;
figh.Color = 'w';
figh.Name = 'Figure 1-figure supplement 1';
numClasses = 10; % though T5-syllables is missing one.
cmap = colormap( 'bone' );


for iDS = 1 : numel( datasets_spectrograms )
    inSpec = load( datasets_spectrograms{iDS} );
    dat = inSpec.dat; % avoid having to refer to nested fields as much
    
    % keep track of each axis' natural limits, then later
    % unify these so all plots have same clim
    axh = [];
    limitExtrema = [inf -inf];
    for iLabel = 2 : numel( dat.uniqueLabels ) % start at 2 because we skip silence
        myLabel =  dat.uniqueLabels{iLabel};
        myRow = iLabel - 1;

        if isempty( dat.(myLabel) )
            axh(myRow) = nan;
            continue
            % skips 'da' in T5-syllables
        end

        mySubplotInd = numel( datasets_spectrograms )*(myRow-1)+iDS; % Subplot numbering is annoying... this will arrange things as in the figure.
        axh(myRow) = subplot( numClasses+1, numel( datasets_spectrograms ), mySubplotInd ); % last row will be for panel C
        
        % restrict to frequencies of interest
        keepTheseWind = (dat.spectrogram.w >= spectrogramFreqMinMax(1) ) & (dat.spectrogram.w <= spectrogramFreqMinMax(2) );
        myPowerSpectrum = dat.(myLabel).psAvg(keepTheseWind,:); % frequency x time
        myW = dat.spectrogram.w(keepTheseWind)/1000; % express in kHz
        imh = imagesc( dat.spectrogram.t, myW, log10( abs( myPowerSpectrum ) ).*10 ); %.*10 expresses it in decibel, /1000 for kHz
        set( axh(myRow), 'YDir', 'normal' );
        myYaxis = get( axh(myRow), 'YAxis' );
        myYaxis.Color = labelColors( myLabel );
        set( axh(myRow), 'TickDir', 'out' );
        colormap( axh(myRow), flipud( cmap ) );   
        myClim = get( axh(myRow), 'CLim' );
        limitExtrema(1) = min( [limitExtrema(1) myClim(1)] );
        limitExtrema(2) = max( [limitExtrema(2) myClim(2)] );
        box off
    end
    for i = 1 : numel( axh )
        if ~isnan( axh(i) ) % skips missing subplots
            set( axh(i), 'CLim', limitExtrema );
        end
    end
    linkaxes( axh(~isnan(axh))  )

end

% ------------------------------------------
%% Reaction times
% ------------------------------------------
% Plots distributions of reaction times for each dataset, as in
% panel C
iSubplot = 1;
for iDS = 1 : numel( datasets_RTs )
    axh_RTs(iDS) = subplot(numClasses+1, numel( datasets_spectrograms ), numClasses*numel( datasets_spectrograms )+iSubplot ); % make panel
    % load the relevant data
    inRT = load( datasets_RTs{iDS} );
    histh = histogram( inRT.RTs, 'FaceColor', 'k' );
    histh.BinWidth = 100;
    ylim([0 max( histh.Values ) + 5] );
    xlim( [500 2000] ); % known already
    line( [median( inRT.RTs ) median( inRT.RTs )], [0 max( histh.Values ) + 5], 'Color', 'b' ) ;
    title( regexprep( datasets_RTs{iDS}, '_RTs', '') );
    
    if iDS == 1
        ylabel( 'Trials' )
    end
    iSubplot = iSubplot + 1;
    box off
end
linkaxes( axh_RTs );