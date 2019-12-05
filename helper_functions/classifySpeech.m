% classifySpeech.m
%
% Does a leave-one-out classification on neural features + corresponding labels, using a provided
% set of classification  parameters.
%
% USAGE: [ result ] = classifySpeech( data, labels, params )
%
% EXAMPLE: classifySpeech( dat.neuralFeatures, dat.labels, params );
%
% INPUTS:
%     data                      Trials x features matrix of data.
%     labels                    Trials x 1 cell list of labels of each trial.
%     params                    Classification parameters.
%
% OUTPUTS:
%     predictedLabels           Trials x 1 cell list of leave-one-out predicted labels
%
% Created by Sergey Stavisky on 23 Sep 2017 using MATLAB version 8.5.0.197613 (R2015a)

 function [ predictedLabels ] = classifySpeech( data, labels, params )

    kernelFunction = 'linear'; % default

    uniqueLabelsStr = unique( labels );
    numTrials = numel( labels );

   
    % record how many classification features there are total
    result.datMatSize = size( data );
    
    %% Ready it for leave one out
    % Parameters for each SVM
    tSVM = templateSVM( 'Standardize', true, 'KernelFunction', kernelFunction, 'OutlierFraction', ...
        params.outlierFraction );
    

    clear('LOOpredictedLabels');
    fprintf('LOO Trial    ');
    for iTrial = 1 : numTrials % recommended if your machine has few cores (e.g., laptop).
%     parfor iTrial = 1 : numTrials % recommended if you have many cores (e.g, beefy workstation).
        fprintf('\b\b\b\b%4i', iTrial)

        trainDat = data;
        trainLabels = labels;
        % leave this trial out of training
        trainDat(iTrial,:) = []; 
        trainLabels(iTrial) = [];
        testDat = data(iTrial,:);
        testLabel = labels{iTrial};
        
        % Multi-class model of binary support vector machines
        Mdl = fitcecoc( trainDat, trainLabels, 'Learners', tSVM); 
        LOOpredictedLabels(iTrial,1) = predict( Mdl, testDat );
    end
   
    predictedLabels = LOOpredictedLabels;
    
end