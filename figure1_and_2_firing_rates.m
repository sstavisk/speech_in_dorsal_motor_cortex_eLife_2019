% Generates firing rates during speaking different syllables as in Figure 1D and Figure 1-figure supplement 2 of Stavisky
% et al. eLife 2019.
% Also generates firing rates during speaking different short words (Figure 1-figure supplement 4).
% Also generates firing rates during different orofacial movements (Figure 2A).
% 
% USAGE: The script can be pointed to several available datasets and units/channels. By default it will plot
% the first two example units in Figure 1D. To examine other neural recordings,
% uncomment/modify the dataset and channel selection in the 'Dataset specification' code block below.
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
% % Figure 1D:
dataset = 'T5-syllables_psth_sortedUnits'; plotChannels = {'unit4(4)_array1_elec13(13)', 'unit15(3)_array2_elec4(100)'}; % left and center panels
% dataset = 'T8-syllables_psth_sortedUnits'; plotChannels = {'unit35(29)_array2_elec68(164)'}; % right panel

% % Figure 1-figure supplement 2:
% dataset = 'T5-syllables_psth_TCs'; plotChannels = {'chan_2.4', 'chan_1.37', 'chan_1.13', 'chan_2.91'}; % top and middle rows
% dataset = 'T8-syllables_psth_TCs'; plotChannels = {'chan_2.68', 'chan_1.35'}; % bottom row

% % Figure 1-figure supplement 4A:
% dataset = 'T5-words_psth_sortedUnits'; plotChannels = {'unit24(13)_array2_elec85(181)', 'unit5(5)_array1_elec20(20)'}; % top-left, middle-left panels
% dataset = 'T5-words_psth_TCs'; plotChannels = {'chan_2.2', 'chan_2.95'}; % top-right,  middle-right panels
% dataset = 'T8-words_psth_sortedUnits'; plotChannels = {'unit47(37)_array2_elec84(180)'}; % bottom-left panel
% dataset = 'T8-words_psth_TCs'; plotChannels = {'chan_2.39'}; % bottom-right panel

% % Figure 2A:
% dataset = 'T5-orofacial_psth_sortedUnits'; plotChannels = {'unit15(3)_array2_elec4(100)', 'unit4(4)_array1_elec13(13)', 'unit24(12)_array2_elec67(163)' };
% dataset = 'T8-orofacial_psth_sortedUnits'; plotChannels = { 'unit35(29)_array2_elec68(164)', 'unit36(30)_array2_elec71(167)', 'unit32(26)_array2_elec63(159)'};


%% Load the data
load( dataset, 'dat' ); % variable will be called 'dat'
uniqueLabels = dat.uniqueLabels; % conditions that will be plotted


%% Aesthetics
% Define the specific colormap
colors = [];
FaceAlpha = 0.3; 
legendLabels = {};
for iLabel = 1 : numel( uniqueLabels )
   colors(iLabel,1:3) = labelColors( uniqueLabels{iLabel} ); 
   legendLabels{iLabel} = sprintf('%s (n=%i)', uniqueLabels{iLabel}, dat.(uniqueLabels{iLabel}).numTrials );
end


%% Firing rates plots
% compute how long each event-aligned time window is, so that the subplots can be made of
% the right size such that time is uniformly scaled along the horizontal axis
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
    
for iCh = 1 : numel( plotChannels )
    chanInd = find( strcmp( dat.channelNames, plotChannels{iCh}) );  
    figh = figure;
    figh.Color = 'w';
    titlestr = sprintf('psth %s %s', dataset, plotChannels{iCh});
    figh.Name = titlestr;
    axh = [];
    myMax = 0; % will be used to track max FR across all conditions.

    for iEvent = 1 : numel( dat.params.alignEvent )
        % Loop through temporal events
        axh(iEvent) = subplot(1, numel( dat.params.alignEvent ), iEvent); hold on;     
        % make width proportional to this epoch's duration
        myPos =  get( axh(iEvent), 'Position');
        set( axh(iEvent), 'Position', [epochStartPosFraction(iEvent) myPos(2) epochWidthsFraction(iEvent) myPos(4)] )
        xlabel(['Time ' dat.params.alignEvent{iEvent} ' (s)']);    
        
        for iLabel = 1 : numel( uniqueLabels )
            myLabel = uniqueLabels{iLabel};
            myX = dat.(myLabel).t{iEvent};
            myY = dat.(myLabel).psthMean{iEvent}(:,chanInd);
            myMax = max([myMax, max( myY )]);
            plot( myX, myY, 'Color', colors(iLabel,:), ...
                'LineWidth', 1 );          
            mySem = dat.(myLabel).psthSem{iEvent}(:,chanInd);
            [px, py] = meanAndFlankingToPatchXY( myX, myY, mySem );
            h = patch( px, py, colors(iLabel,:), 'FaceAlpha', FaceAlpha, ...
                'EdgeColor', 'none');
            myMax = max([myMax, max( myY+mySem )]);                
        end
        
        % PRETTIFY
        % make horizontal axis nice
        xlim([myX(1), myX(end)])
        % make vertical axis nice
        if iEvent == 1
            ylabel( 'Spikes/s' );
        else
            % hide it
            yaxh = get( axh(iEvent), 'YAxis');
            yaxh.Visible = 'off';
        end
        set( axh(iEvent), 'TickDir', 'out' )
    end
    
    linkaxes(axh, 'y');
    ylim([0 ,ceil( myMax ) + 1]);
    % add legend
    axes( axh(1) );
    MakeDumbLegend( legendLabels, 'Color', colors );
end