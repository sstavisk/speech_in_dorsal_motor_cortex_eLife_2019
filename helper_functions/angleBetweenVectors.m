% angleBetweenVectors.m
%
% Simple function returns the angle between two vectors.
%
% USAGE: [ degAngle ] = angleBetweenVectors( vecA, vecB )
%
% EXAMPLE:
%
% INPUTS:
%     vecA                      vector 1
%     vecB                      vector 2
%
% OUTPUTS:
%     degAngle                  angle between them, in degrees
%
% Created by Sergey Stavisky on 10 Aug 2015

function [ degAngle ] = angleBetweenVectors( vecA, vecB )
    costheta = dot(vecA,vecB)/(norm(vecA)*norm(vecB));
    degAngle = rad2deg( real( acos(costheta) ) ); % real so with 0 angle it doesnt give weird complex tiny numbers
end