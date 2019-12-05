% indexEvenSpaced.m
%
% Given a number of indices and how many evenly sampled samples you want from it,
% returns the index of these samples. Guarantees first and last sample are included.
%
% E.g. if 100 elemnts and want 4 samples, will return [1, 34, 67, 100]
% USAGE: [ sampleInds ] = indexEvenSpaced( numElements, numSamples )
%
% EXAMPLE:   sliceInds = indexEvenSpaced( numel( datasets ), numSlices );
%
% INPUTS:
%     numElements               How many total elements there are in the list
%     numSamples                How many samples to grab
%
% OUTPUTS:
%     sampleInds                
%
% Created by Sergey Stavisky on 11 May 2015

function [ sampleInds ] = indexEvenSpaced( numElements, numSamples )
    sampleInds =  0  : numElements/(numSamples-1) : numElements;
    sampleInds = ceil( sampleInds );
    sampleInds(1) = 1;
end