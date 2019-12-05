% whichGroupIsWhichJpca.m
%
% Helps match condition indices (from Data that goes into jPCA.m) to the reuslting plots,
% by looking at how the conditions are ordered left to right at the first time sample along jPC1
% (which will correspond to brightest green to darkest red).
%
% USAGE: [ indsLeftToRight, cmap ] = whichGroupIsWhichJpca( Projection )
%
% EXAMPLE:
%
% INPUTS:
%     Projection                
%
% OUTPUTS:
%     indsLeftToRight           indices into Projection, ordered by lowest to highest jPC1
%     cmap                      colormap (already sorted according to indsLeftToRight) of
%                               each conditions' traces used in Mark's code.
%
% Created by Sergey Stavisky on 05 Oct 2017 using MATLAB version 8.5.0.197613 (R2015a)

 function [ indsLeftToRight, cmap ] = whichGroupIsWhichJpca( Projection )
    initialPC1 = arrayfun( @(x) x.proj(1,1), Projection ); % PC1 of each condition
    [~, indsLeftToRight] = sort( initialPC1, 'ascend'); 

    cmap = redgreencmap( numel(indsLeftToRight), 'interpolation', 'linear');
    cmap = cmap(indsLeftToRight,:);
end