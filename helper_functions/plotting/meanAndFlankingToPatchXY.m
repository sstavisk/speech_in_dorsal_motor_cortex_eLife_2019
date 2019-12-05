% meanAndFlankingToPatchXY.m
%
% Helper function that takes a mean and flanking offset (e.g. mean +- s.d.) and returns
% coordinates for a patch object that would plot the error bars. If series starts or ends
% with nans, these get trimmed off.
%
% USAGE: [ x, y ] = meanAndFlankingToPatchXY( t, m, s )
%
% EXAMPLE:                     [px py] = meanAndFlankingToPatchXY( t, myY, mySem );

%
% INPUTS:
%     t                      Nx1 series of horizontal coordinate (e.g. time).
%     m                      Nx1 series of mean points 
%     s                      Nx1 series of flanking points
%
% OUTPUTS:
%     x                      x coordinates for patch object
%     y                      y coordinates for patch object
%
% Created by Sergey Stavisky on 04 Jan 2018 using MATLAB version 9.3.0.713579 (R2017b)

 function [ x, y ] = meanAndFlankingToPatchXY( t, m, s )
    % find first and last non-nan in input mean.
    startInd = find( ~isnan( m ), 1, 'first');
    endInd = find( ~isnan( m ), 1, 'last' );
    m = m(startInd:endInd);
    s = s(startInd:endInd);
    t = t(startInd:endInd);

    % force columns
    if size( t, 2) ~= 1
        t = t';
    end
    if size( m, 2) ~= 1
        m = m';
    end
    if size( s, 2) ~= 1
        s = s';
    end
    x = [t; flipud( t )];
    y = [ m - s  ; flipud( m + s )];    
end