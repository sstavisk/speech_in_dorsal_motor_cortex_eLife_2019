% Calculates and compares neural modulation during attempted orofacial/speech and arm/hand
% movements as in Figure 2-figure supplement 3 of Stavisky et al. eLife 2019.
%
% DEPENDENCIES:
%   Add to your MATLAB path the /datasets and /helper_functions that are included in our
%   code and data release.
% 
%   This script was developed and tested in MATLAB_R2017b, but it ought to be robust across
%   versions within a couple years of this.
%
%   The neural distance metric uses the cvDistance.m function available at
%   https://github.com/fwillett/cvVectorStats. For convenience a static copy is provided with
%   this paper's code.
%
% This code is associated with "Neural ensemble dynamics in dorsal motor cortex during
% speech in people with paralysis", eLife (2019).
%
% Copyright Sergey D. Stavisky, November 2019, Stanford Neural Prosthetics Translational
% Laboratory.

%% Dataset specification
dataset = 'T5_comparisons_TCs.mat';

%% Load the data
load( dataset ); % loads variable 'dat'

faceMovements = {...
    'SayBa';
    'SayGa';
    'TongueUp';
    'TongueDown';
    'TongueLeft';
    'TongueRight';
    'MouthOpen';
    'JawClench';
    'LipsPucker';
    };

armMovements = {...
    'ShoShrug';
    'ArmRaise';
    'ElbowFlex';
    'WristExt';
    'HandClose';
    'HandOpen';
    'IndexRaise';
    'ThumbRaise';
    };

%% Calculate each conditions' modulation
allConditions = [faceMovements;armMovements];
allConditionsModulation = []; % will fill as we go 
for iCond = 1 : numel( allConditions )
    myCondName = allConditions{iCond};
    % neural distance comparison to similar time epoch of the 'do nothing' condition.
    allConditionsModulation(iCond) = cvDistance( dat.(myCondName).fr, dat.Nothing.fr );
end

%% Plot the modulations
figh = figure('Position', [10   10   300   500], 'Color', 'w');
figh.Name = 'Fig 2S3 orofacial vs arm comparison';
hold on;
% plot means with bar plot
faceMods = allConditionsModulation(1:numel( faceMovements ));
armMods = allConditionsModulation(numel( faceMovements )+1:end);
bar( 1, mean( faceMods ), 'r', 'FaceAlpha', 0.5 ) 
bar( 2, mean( armMods ), 'k', 'FaceAlpha', 0.5 )
axh = gca;
axh.XTick = [1 2];
axh.XTickLabel = {'orofacial & speech', 'arm & hand'};
ylabel( 'Population modulation (Hz)' );

% plot individual conditions
scatter( ones( size( faceMods ) ), faceMods, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'SizeData', 24 );
scatter( 2.*ones( size( armMods ) ), armMods, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'SizeData', 24 );

[p,h] = ranksum( faceMods, armMods );
fprintf('Mean orofacial & speech: %.2f, mean arm & speech: %.2f (%.2fx)\n', mean( faceMods ), mean( armMods), mean( armMods ) / mean( faceMods ) );
fprintf('  p=%g (rank-sum test)\n', p );

% Label these
for i = 1 : numel( faceMovements )
   myLabel = faceMovements{i}; 
   myY = faceMods(i);
   th = text( 0.9, myY, myLabel, 'HorizontalAlignment', 'right', 'Color', 'r' );
end
for i = 1 : numel( armMods )
   myLabel = armMovements{i}; 
   myY = armMods(i);
   th = text( 2.1, myY, myLabel, 'HorizontalAlignment', 'left' );
end