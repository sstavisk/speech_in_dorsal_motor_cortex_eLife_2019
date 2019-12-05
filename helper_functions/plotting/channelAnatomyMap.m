% channelAnatomyMap.m
%
% Returns x- and y- location lookup of each channel based on rough anatomical location.
% Useful for rendering neural activity on a more meaningful map. 
%
% USAGE: [ chanMap ] = channelAnatomyMap( arrays )
%
% EXAMPLE: chanMap = channelAnatomyMap({'T5_lateral', 'T5_medial'})
%
% INPUTS:
%     arrays                    cell list of which arrays to use
%                               
% OUTPUTS:
%     chanMap                  where to draw each chanenl, in mm
%       .x
%       .y
%       .xlim                  suggested axis limits
%       .ylim                  suggested axis limits
% Created by Sergey Stavisky on 13 Jul 2014

function [ chanMap ] = channelAnatomyMap( arrays )
    NUMCHANS = 96; % per array

    % get maps for each array
    for iArray = 1 : numel( arrays )
        [emap{iArray}, mapRotate{iArray}, mapTranslate{iArray}] = arrayMapHumans( arrays{iArray});
    end
    
    NUMARRAYS = numel( arrays );
    chanMap.x = [];
    chanMap.y = [];
    
    for iArray = 1 : NUMARRAYS
        x = zeros( NUMCHANS, 1 );
        y = zeros( NUMCHANS, 1 );
        %
        for iChan = 1 : NUMCHANS
            % get row/column of this channel
            [row, col] = find( emap{iArray} == iChan );
            % offset by 400um per row/col
            x(iChan) = col*.4;
            y(iChan) = row*-.4;
        end
        % locations based solely on the array layout
        % center it so rotation goes around its middle
        x = x - mean(x);
        y = y - mean(y);
       
        % rotate this array
        rotMatrix = [cos( deg2rad( mapRotate{iArray} ) ) -sin( deg2rad( mapRotate{iArray} ) ) ;
            sin( deg2rad( mapRotate{iArray} ) )  cos( deg2rad( mapRotate{iArray} ) ) ];
        rotated = (rotMatrix* [x y]')';
        x = rotated(:,1);
        y = rotated(:,2);
       
        % translate the array
        x = x + mapTranslate{iArray}(1);
        y = y + mapTranslate{iArray}(2);
         
        % add to final list
        chanMap.x = [chanMap.x; x];
        chanMap.y = [chanMap.y; y];
    end
    
    % Get suggested axis limits
    chanMap.xlim = [min( chanMap.x ) - 0.1, max( chanMap.x ) + 0.1];
    chanMap.ylim = [min( chanMap.y ) - 0.1, max( chanMap.y ) + 0.1];
end