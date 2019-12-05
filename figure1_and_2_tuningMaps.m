% Generates maps of which electrodes have significant modulation (and for how many conditions) during speaking
% different syllables as in Figure 1B insets of Stavisky et al. eLife 2019.
%   Also generates electrode array maps for each spoken sound separately, as in Figure
% 1-figure supplement 3A, and histograms showing the distribution of how many different syllables 
% evoke a significant firing rate change for electrodes' threshold crossings or sorted
% units, as in Figure 1-figure supplement 3B.
%   Also generates the analogous electrode array maps and histograms for orofacial
% movements, as in Figure 2B and Figure 2-figure supplement 1.
%
% USAGE: The script can be pointed to several available datasets. By default it will plot
% the participant T5, syllables dataset population modulation (Figure 1B, inset, top. To examine other datasets,
% uncomment/modify the dataset selection in the 'Dataset specification' code block below.
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


% % Figure 1B inset AND figure 1-figure supplement 3A
dataset = 'T5-syllables_tuningMaps_TCs'; saturateMapOneBelow = false; % top
% dataset = 'T8-syllables_tuningMaps_TCs'; saturateMapOneBelow = true; % bottom; latter option is so that the two plots in Fig 1 have same scale despite T5 missing a phoneme

% % Figure 1-figure supplement 3B sorted units 
% dataset = 'T5-syllables_tuningMaps_sortedUnits'; saturateMapOneBelow = false; % left
% dataset = 'T8-syllables_tuningMaps_sortedUnits'; saturateMapOneBelow = false; % right

% % Figure 2B AND Figure 2-figure supplement 1's threshold crossings results
% dataset = 'T5-orofacial_tuningMaps_TCs'; saturateMapOneBelow = false; % left
% dataset = 'T8-orofacial_tuningMaps_TCs'; saturateMapOneBelow = false; % right

% % Figure 2-figure supplement 1 single neuron results
% dataset = 'T5-orofacial_tuningMaps_sortedUnits'; saturateMapOneBelow = false; % left
% dataset = 'T8-orofacial_tuningMaps_sortedUnits'; saturateMapOneBelow = false; % right

% Analysis specifications
reportChannelsBelowPvalue = 0.05; % note: it gets Bonferonni corrected later


%% Load the data
load( dataset, 'dat' ); % variable will be called 'dat'
uniqueLabels = dat.uniqueLabels; % conditions that will be plotted

switch dataset(1:2)
    case 'T5'
        arrayMaps ={'T5_lateral', 'T5_medial'};
    case 'T8'
        arrayMaps = {'T8_lateral', 'T8_medial'};
end

% whether these are sorted units' SUA or threshold crossing spikes.
if strfind( dat.channelNames{1}, 'chan')
    isSorted = false;
else
    isSorted = true;
end

%% Compare each label to silence, for each channel/unit
numChans = numel( dat.channelNames );
result.pValueVsSilent = nan( numChans, numel( uniqueLabels ) );

for iChan = 1 : numChans
    % silence/stayStill rates
    myLabel = uniqueLabels{1};
    silenceDat = dat.(myLabel).meanWindowActivity(:,iChan);
    
    % loop through each of the other labels and compare them to the silent label
    for iLabel = 2 : numel( uniqueLabels )
        myLabel = uniqueLabels{iLabel};
        myDat = dat.(myLabel).meanWindowActivity(:,iChan);
        result.pValueVsSilent(iChan,iLabel) = ranksum( myDat, silenceDat );
    end  
end

% p value cutoff 
testVal = (reportChannelsBelowPvalue ./ (numel( uniqueLabels )-1) ); % Bonefonni correction
signifPvalueVsSilent = result.pValueVsSilent < testVal;
numEachTunedTo = sum( signifPvalueVsSilent, 2 ); % how many conditions is each channel/unit tuned to?
% restrict to responders, meaning having significant tuning to at least one channel 
signifResponders = find( numEachTunedTo >= 1 );

% For channels that responded to just one label, identify how often each label appeared as
% that 'sparse label';
sparseChans = find( numEachTunedTo == 1 );
sparseLabels = sum( signifPvalueVsSilent(sparseChans,:), 1 );

% ------------------------------------------------------------
%% Tuning plotted on the electrode arrays
% ------------------------------------------------------------
chanMap = channelAnatomyMap( arrayMaps );

figh = figure;
figh.Renderer = 'painters';
figh.Position = [10 10 800 1200]; % if it starts too small, its auto-layout messes up with 6 rows
titlestr = sprintf('Array tuning %s', dataset );
figh.Name = titlestr;

% grayscale colormap for the # tuned sum plot (Figure 1B, insets)
N = numel( uniqueLabels  );
if exist( 'saturateMapOneBelow', 'var' ) && saturateMapOneBelow
    N = N - 1;
    fprintf(' \nNOTE: Units tuned to %i syllables now changed to %i to saturate plot\n\n', N, N-1)
    numEachTunedTo(numEachTunedTo==N) = N-1;    
end
graymap = flipud( bone( N + 1) ); 
graymap = graymap(2:end,:); % so it doesn't start at pure white

NUM_ROWS = 6;
% will subplot as 2 columns of 6. Last one is for the # tuned plot (sum across all subplots)
for iLabel = 2 : numel( uniqueLabels )+1 % skips silence/stayStill
    drawnAlready = []; % will track which electrodes were drawn as having something on them
    % Create axis.
    axh = subplot( NUM_ROWS, 2, iLabel-1);
    hold on;
    axh.Position =  axh.OuterPosition;
    axh.XAxis.Visible = 'off';
    axh.YAxis.Visible = 'off';
    axh.XLim = chanMap.xlim;
    axh.YLim = chanMap.ylim;
    axis equal
    
    for iChan = 1 : numChans
        if iLabel == numel( uniqueLabels )+1
            % This is the "bonus" plot showing how many are tuned across conditions.
            myDat = numEachTunedTo(iChan);
            mySize = 36;
            myColor = graymap(myDat+1,:);
            
            % Disabled channels
            % here they look identical to 0, since I don't want to imply that no-activity but not
            % disabled are somehow more biologically-empty than truly disabled channels.
            disabledSize = 4;
            disabledColor = graymap(1,:);            
        else
            % This is for the single-condition (syllable/movement) plots
            myLabel = uniqueLabels{iLabel};
            myDat = result.pValueVsSilent(iChan,iLabel); % p-value
            if myDat < testVal
                mySize = 36;
                myColor = labelColors( myLabel );
            else
                mySize = 16; % distinct from disabled channels
                myColor = [.8 .8 .8];
            end
        end
      
        if isSorted
            myChanNum = dat.unitElectrodeTo192(iChan);
        else
            myChanNum = ChannelNameToNumber( dat.channelNames{iChan} );
        end
        drawnAlready(end+1) = myChanNum;
        x = chanMap.x(myChanNum);
        y = chanMap.y(myChanNum);
        
        % draw this channel
        scatter( x, y, mySize, myColor, 'filled' )
    end
    % draw the remainder of electrodes (disabled channels)
    drawThese = setdiff( 1:192, drawnAlready );
    disabledSize = 4;
    disabledColor = [.8 .8 .8];
    for iLeftover = 1 : numel( drawThese )
        x = chanMap.x(drawThese(iLeftover));
        y = chanMap.y(drawThese(iLeftover));
        % draw the point
        scatter( x, y, disabledSize, disabledColor, 'filled' )
    end
end

% This next block of code males the relationship between anatomy mm and image pixels
% consistent across participants. This makes final figure making a lot easier later.
STANDARD_MMperNormalizedUnits = 36.2076;
childs = figh.Children;
numAxes = numel( childs );

axes( childs(1 ) )
axis equal
c1Xlim = childs(1).XLim;
c1YLim = childs(1).YLim;

xRange = range( childs(1).XLim );
yRange = range( childs(1).YLim );
myPos = childs(1).Position;
MMperNormalizedUnits = yRange/myPos(4);
scaleFactor =  MMperNormalizedUnits / STANDARD_MMperNormalizedUnits;
childs(1).Position = [myPos(1) myPos(2), scaleFactor*myPos(3), scaleFactor*myPos(4)];
c1Width = scaleFactor*myPos(3);
c1Height = scaleFactor*myPos(4);
for iAxh = 1 : numAxes
    axes( childs(iAxh ) )
    axh = gca;
    axh.XLim = c1Xlim;
    axh.YLim = c1YLim;
    
    myPos = axh.Position;
    axh.Position = [myPos(1) myPos(2), c1Width, c1Height];
end

% ------------------------------------------------------------
%% Histograms of how many conditions each channel/unit responds to
% ------------------------------------------------------------

% The first bar is stacked and colored by label. The rest are just counts of #
% channels/units of that bin.
barMat = sparseLabels;
for count = 2 : max( numEachTunedTo )
    barMat(count,1) = nnz( numEachTunedTo == count);
end

figh = figure;
axh = axes;
barh = bar( barMat, 'stacked');

% color all other bars black.
barh(1).FaceColor = [0 0 0];
% color the sparse responders by their label's color
for iLabel = 2 : numel( uniqueLabels ) % skip 1 because that's silence
   barh(iLabel).FaceColor = labelColors( uniqueLabels{iLabel});
end

barh(1).BarWidth = 0.9;
xlabel( '# Labels Responding To'); 
if ~isSorted
    yStr = 'Channels';
else
    yStr = 'Units';
end
ylabel( ['# ' yStr] )
titlestr = sprintf('Labels Tuned To Histogram %s', dataset );
figh.Name = titlestr;
% draw median and report results
line( [median(  numEachTunedTo(signifResponders) ) median(  numEachTunedTo(signifResponders) )], axh.YLim, 'Color', [.5 .5 .5] );
fprintf('%i/%i %s have significant response to at least one label at p = %f (ranksum). Median = %.2f\n', ...
    numel( signifResponders ), size( result.pValueVsSilent,1 ), yStr, testVal, median(  numEachTunedTo(signifResponders) ) );


%% ------------------------------------------------------------
%% Overlaps between neural tuning to syllable and orofacial movements (Figure 2C)
% ------------------------------------------------------------
% These are just comparisons of the lists of responding channels/units from the above analyses of the 
% datasets of interest. For the reader's convenience, these lists are hard-coded below and then plotted.
figh = figure;
figh.Color = 'w';
figh.Name = 'Syllable vs orofacial overlaps';
for iPlot = 1 : 4
    switch iPlot
        case 1
            plotName = 'T5 sorted units';
            syllables_responders = {'unit2(2)_array1_elec4(4)', 'unit3(3)_array1_elec12(12)', 'unit4(4)_array1_elec13(13)', 'unit5(5)_array1_elec15(15)', 'unit11(11)_array1_elec40(40)', 'unit13(1)_array2_elec2(98)', 'unit15(3)_array2_elec4(100)', 'unit16(4)_array2_elec11(107)', 'unit17(5)_array2_elec22(118)', 'unit22(10)_array2_elec39(135)', 'unit28(16)_array2_elec76(172)', 'unit30(18)_array2_elec81(177)', 'unit31(19)_array2_elec85(181)'};
            syllables_nonResponders = {'unit7(7)_array1_elec23(23)', 'unit9(9)_array1_elec34(34)', 'unit10(10)_array1_elec37(37)', 'unit20(8)_array2_elec36(132)', 'unit23(11)_array2_elec52(148)', 'unit24(12)_array2_elec67(163)', 'unit26(14)_array2_elec74(170)', 'unit27(15)_array2_elec75(171)', 'unit32(20)_array2_elec87(183)'};
            orofacial_responders = {'unit2(2)_array1_elec4(4)', 'unit3(3)_array1_elec12(12)', 'unit4(4)_array1_elec13(13)', 'unit9(9)_array1_elec34(34)', 'unit11(11)_array1_elec40(40)', 'unit13(1)_array2_elec2(98)', 'unit15(3)_array2_elec4(100)', 'unit16(4)_array2_elec11(107)', 'unit17(5)_array2_elec22(118)', 'unit20(8)_array2_elec36(132)', 'unit22(10)_array2_elec39(135)', 'unit24(12)_array2_elec67(163)', 'unit27(15)_array2_elec75(171)', 'unit28(16)_array2_elec76(172)', 'unit30(18)_array2_elec81(177)', 'unit31(19)_array2_elec85(181)', 'unit32(20)_array2_elec87(183)'};
            orofacial_nonResponders = {'unit5(5)_array1_elec15(15)', 'unit7(7)_array1_elec23(23)', 'unit10(10)_array1_elec37(37)', 'unit23(11)_array2_elec52(148)', 'unit26(14)_array2_elec74(170)'};

        case 2
            plotName = 'T8 sorted units';
            syllables_responders = {'unit3(3)_array1_elec36(36)', 'unit4(4)_array1_elec45(45)', 'unit5(5)_array1_elec65(65)', 'unit13(7)_array2_elec19(115)', 'unit20(14)_array2_elec39(135)', 'unit21(15)_array2_elec41(137)', 'unit22(16)_array2_elec43(139)', 'unit32(26)_array2_elec63(159)', 'unit35(29)_array2_elec68(164)', 'unit36(30)_array2_elec71(167)', 'unit44(38)_array2_elec84(180)', 'unit48(42)_array2_elec93(189)'};
            syllables_nonResponders = {'unit6(6)_array1_elec80(80)', 'unit9(3)_array2_elec8(104)', 'unit10(4)_array2_elec9(105)', 'unit11(5)_array2_elec10(106)', 'unit12(6)_array2_elec12(108)', 'unit15(9)_array2_elec30(126)', 'unit23(17)_array2_elec43(139)', 'unit24(18)_array2_elec47(143)', 'unit33(27)_array2_elec65(161)', 'unit34(28)_array2_elec66(162)', 'unit40(34)_array2_elec78(174)', 'unit41(35)_array2_elec79(175)', 'unit47(41)_array2_elec90(186)'};
            orofacial_responders = {'unit3(3)_array1_elec36(36)', 'unit4(4)_array1_elec45(45)', 'unit5(5)_array1_elec65(65)', 'unit11(5)_array2_elec10(106)', 'unit13(7)_array2_elec19(115)', 'unit20(14)_array2_elec39(135)', 'unit21(15)_array2_elec41(137)', 'unit22(16)_array2_elec43(139)', 'unit23(17)_array2_elec43(139)', 'unit24(18)_array2_elec47(143)', 'unit32(26)_array2_elec63(159)', 'unit35(29)_array2_elec68(164)', 'unit36(30)_array2_elec71(167)', 'unit40(34)_array2_elec78(174)', 'unit41(35)_array2_elec79(175)', 'unit44(38)_array2_elec84(180)', 'unit48(42)_array2_elec93(189)'};
            orofacial_nonResponders = {'unit6(6)_array1_elec80(80)', 'unit9(3)_array2_elec8(104)', 'unit10(4)_array2_elec9(105)', 'unit12(6)_array2_elec12(108)', 'unit15(9)_array2_elec30(126)', 'unit33(27)_array2_elec65(161)', 'unit34(28)_array2_elec66(162)', 'unit47(41)_array2_elec90(186)'};

        case 3
            plotName = 'T5 -4.5 x RMS threshold crossings';
            syllables_responders = {'chan_1.1', 'chan_1.5', 'chan_1.6', 'chan_1.7', 'chan_1.8', 'chan_1.9', 'chan_1.10', 'chan_1.11', 'chan_1.13', 'chan_1.17', 'chan_1.19', 'chan_1.20', 'chan_1.22', 'chan_1.23', 'chan_1.28', 'chan_1.30', 'chan_1.32', 'chan_1.33', 'chan_1.34', 'chan_1.35', 'chan_1.37', 'chan_1.40', 'chan_1.48', 'chan_1.54', 'chan_1.61', 'chan_1.66', 'chan_2.1', 'chan_2.2', 'chan_2.3', 'chan_2.4', 'chan_2.6', 'chan_2.7', 'chan_2.8', 'chan_2.9', 'chan_2.10', 'chan_2.11', 'chan_2.12', 'chan_2.19', 'chan_2.21', 'chan_2.22', 'chan_2.30', 'chan_2.32', 'chan_2.34', 'chan_2.36', 'chan_2.37', 'chan_2.38', 'chan_2.41', 'chan_2.43', 'chan_2.45', 'chan_2.49', 'chan_2.55', 'chan_2.57', 'chan_2.60', 'chan_2.62', 'chan_2.63', 'chan_2.66', 'chan_2.68', 'chan_2.71', 'chan_2.74', 'chan_2.75', 'chan_2.76', 'chan_2.77', 'chan_2.79', 'chan_2.82', 'chan_2.84', 'chan_2.85', 'chan_2.86', 'chan_2.87', 'chan_2.88', 'chan_2.89', 'chan_2.91', 'chan_2.93', 'chan_2.94' };
            syllables_nonResponders  = {'chan_1.3', 'chan_1.12', 'chan_1.14', 'chan_1.15', 'chan_1.16', 'chan_1.18', 'chan_1.25', 'chan_1.36', 'chan_1.38', 'chan_1.39', 'chan_1.47', 'chan_1.56', 'chan_1.65', 'chan_1.68', 'chan_1.95', 'chan_2.13', 'chan_2.14', 'chan_2.15', 'chan_2.20', 'chan_2.26', 'chan_2.33', 'chan_2.39', 'chan_2.47', 'chan_2.52', 'chan_2.59', 'chan_2.61', 'chan_2.67', 'chan_2.69', 'chan_2.80', 'chan_2.81', 'chan_2.83' };
            orofacial_responders = {'chan_1.3', 'chan_1.6', 'chan_1.7', 'chan_1.8', 'chan_1.9', 'chan_1.10', 'chan_1.11', 'chan_1.12', 'chan_1.13', 'chan_1.17', 'chan_1.19', 'chan_1.20', 'chan_1.22', 'chan_1.25', 'chan_1.28', 'chan_1.33', 'chan_1.34', 'chan_1.35', 'chan_1.37', 'chan_1.38', 'chan_1.40', 'chan_1.47', 'chan_1.48', 'chan_1.54', 'chan_1.66', 'chan_2.1', 'chan_2.2', 'chan_2.3', 'chan_2.4', 'chan_2.6', 'chan_2.7', 'chan_2.8', 'chan_2.10', 'chan_2.11', 'chan_2.12', 'chan_2.13', 'chan_2.14', 'chan_2.20', 'chan_2.22', 'chan_2.24', 'chan_2.26', 'chan_2.30', 'chan_2.32', 'chan_2.33', 'chan_2.34', 'chan_2.36', 'chan_2.37', 'chan_2.39', 'chan_2.41', 'chan_2.43', 'chan_2.49', 'chan_2.55', 'chan_2.59', 'chan_2.60', 'chan_2.62', 'chan_2.63', 'chan_2.66', 'chan_2.68', 'chan_2.69', 'chan_2.71', 'chan_2.75', 'chan_2.77', 'chan_2.79', 'chan_2.81', 'chan_2.82', 'chan_2.83', 'chan_2.84', 'chan_2.85', 'chan_2.86', 'chan_2.87', 'chan_2.88', 'chan_2.89', 'chan_2.91', 'chan_2.94' };
            orofacial_nonResponders = {'chan_1.5', 'chan_1.14', 'chan_1.15', 'chan_1.16', 'chan_1.18', 'chan_1.23', 'chan_1.30', 'chan_1.32', 'chan_1.36', 'chan_1.39', 'chan_1.56', 'chan_1.61', 'chan_1.65', 'chan_1.95', 'chan_2.15', 'chan_2.38', 'chan_2.45', 'chan_2.47', 'chan_2.52', 'chan_2.57', 'chan_2.67', 'chan_2.74', 'chan_2.76', 'chan_2.80', 'chan_2.93'} ;

        case 4
            plotName = 'T8 -4.5 x RMS threshold crossings';            
            syllables_responders = {'chan_1.1', 'chan_1.2', 'chan_1.3', 'chan_1.6', 'chan_1.7', 'chan_1.9', 'chan_1.11', 'chan_1.15', 'chan_1.33', 'chan_1.34', 'chan_1.35', 'chan_1.36', 'chan_1.45', 'chan_1.65', 'chan_1.69', 'chan_1.71', 'chan_1.77', 'chan_1.79', 'chan_1.83', 'chan_1.90', 'chan_2.6', 'chan_2.8', 'chan_2.16', 'chan_2.28', 'chan_2.32', 'chan_2.41', 'chan_2.43', 'chan_2.44', 'chan_2.45', 'chan_2.47', 'chan_2.53', 'chan_2.60', 'chan_2.63', 'chan_2.65', 'chan_2.67', 'chan_2.68', 'chan_2.70', 'chan_2.72', 'chan_2.76', 'chan_2.79', 'chan_2.84', 'chan_2.85', 'chan_2.87', 'chan_2.88', 'chan_2.89', 'chan_2.92', 'chan_2.93'};
            syllables_nonResponders = {'chan_1.4', 'chan_1.5', 'chan_1.24', 'chan_1.28', 'chan_1.50', 'chan_1.66', 'chan_1.67', 'chan_1.68', 'chan_1.70', 'chan_1.73', 'chan_1.75', 'chan_1.78', 'chan_1.80', 'chan_1.81', 'chan_1.85', 'chan_1.86', 'chan_1.88', 'chan_2.3', 'chan_2.4', 'chan_2.5', 'chan_2.7', 'chan_2.9', 'chan_2.10', 'chan_2.12', 'chan_2.14', 'chan_2.18', 'chan_2.19', 'chan_2.23', 'chan_2.25', 'chan_2.26', 'chan_2.30', 'chan_2.31', 'chan_2.33', 'chan_2.34', 'chan_2.35', 'chan_2.36', 'chan_2.39', 'chan_2.46', 'chan_2.48', 'chan_2.51', 'chan_2.54', 'chan_2.56', 'chan_2.57', 'chan_2.64', 'chan_2.66', 'chan_2.71', 'chan_2.74', 'chan_2.78', 'chan_2.81', 'chan_2.83', 'chan_2.90', 'chan_2.94', 'chan_2.95', 'chan_2.96'};
            orofacial_responders = {'chan_1.1', 'chan_1.2', 'chan_1.3', 'chan_1.4', 'chan_1.5', 'chan_1.6', 'chan_1.7', 'chan_1.9', 'chan_1.11', 'chan_1.15', 'chan_1.28', 'chan_1.33', 'chan_1.34', 'chan_1.36', 'chan_1.45', 'chan_1.65', 'chan_1.66', 'chan_1.67', 'chan_1.68', 'chan_1.69', 'chan_1.70', 'chan_1.71', 'chan_1.73', 'chan_1.75', 'chan_1.77', 'chan_1.78', 'chan_1.79', 'chan_1.81', 'chan_1.83', 'chan_1.85', 'chan_1.86', 'chan_1.88', 'chan_1.90', 'chan_2.4', 'chan_2.5', 'chan_2.6', 'chan_2.7', 'chan_2.10', 'chan_2.14', 'chan_2.31', 'chan_2.32', 'chan_2.33', 'chan_2.34', 'chan_2.35', 'chan_2.36', 'chan_2.39', 'chan_2.41', 'chan_2.42', 'chan_2.43', 'chan_2.45', 'chan_2.47', 'chan_2.48', 'chan_2.51', 'chan_2.53', 'chan_2.54', 'chan_2.56', 'chan_2.57', 'chan_2.60', 'chan_2.63', 'chan_2.64', 'chan_2.65', 'chan_2.67', 'chan_2.68', 'chan_2.70', 'chan_2.71', 'chan_2.74', 'chan_2.76', 'chan_2.78', 'chan_2.79', 'chan_2.81', 'chan_2.83', 'chan_2.84', 'chan_2.85', 'chan_2.87', 'chan_2.88', 'chan_2.89', 'chan_2.92', 'chan_2.93', 'chan_2.94'};
            orofacial_nonResponders = {'chan_1.24', 'chan_1.35', 'chan_1.50', 'chan_1.56', 'chan_1.80', 'chan_2.3', 'chan_2.8', 'chan_2.9', 'chan_2.12', 'chan_2.16', 'chan_2.18', 'chan_2.19', 'chan_2.28', 'chan_2.30', 'chan_2.44', 'chan_2.46', 'chan_2.66', 'chan_2.72', 'chan_2.90', 'chan_2.95', 'chan_2.96'};
    end
    
    % Calculate overlaps
    syllableOnly = setdiff( syllables_responders, orofacial_responders );
    both = intersect( syllables_responders, orofacial_responders );
    orofacialOnly = setdiff( orofacial_responders, syllables_responders );
    allChans = unique( [syllables_responders, syllables_nonResponders] ); % excludes "non-functional" electrodes defined as <1 Hz during syllables and not responder during orofacial movements
    neither = setdiff( allChans, unique([syllableOnly, both, orofacialOnly]) ); 

    axh = subplot(2,2,iPlot);
    title( plotName )
    hold on;
    axh.XLim = [-1 110+1]; % hard-coded based on max channels, will keep all four plots the same
    myx = [0 numel( syllableOnly )];
    line( myx, [0.5 0.5], 'LineWidth', 5, 'Color', [1 0 0] );
    text( mean( myx ), 1, sprintf('%i', numel( syllableOnly ) ), ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
    
    myx = [numel( syllableOnly ) numel( syllableOnly )+numel(both)];
    line( myx, [0.5 0.5], 'LineWidth', 5, 'Color', [.5 0 .5] );
    text( mean( myx ), 1, sprintf('%i', numel( both ) ), ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
    
    myx = [numel( syllableOnly ) + numel( both ), numel( syllableOnly )+ numel( both ) + numel( orofacialOnly )];
    line( myx, [0.5 0.5], 'LineWidth', 5, 'Color', [0 0 1] );
    text( mean( myx ), 1, sprintf('%i', numel( orofacialOnly ) ), ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
    
    myx = [numel( syllableOnly ) + numel( both ) + numel( orofacialOnly ), numel( syllableOnly )+numel( both ) + numel( orofacialOnly ) + numel( neither )];
    line( myx, [0.5 0.5], 'LineWidth', 5, 'Color', [0.5 0.5 0.5] );
    text( mean( myx ), 1, sprintf('%i', numel( neither ) ), ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )
    
    axh.XAxis.Visible = 'off';
    axh.YAxis.Visible = 'off';
    axh.Color = 'none';
end