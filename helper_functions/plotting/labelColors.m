% labelColors.m
%
% Lookup table of colors used for each condition. These come from the
% 10-category 'qualitative' maps from http://colorbrewer2.org/.
%
% USAGE: [ c ] = speechColors( label )
%
% EXAMPLE:  color = speechColors( 'ga' ); 
%
% INPUTS:
%     label                     string
%
% OUTPUTS:
%     c                         1x3 RGB color
%
% Created by Sergey Stavisky on 14 Dec 2017 using MATLAB version 9.3.0.713579 (R2017b)

 function [ c ] = labelColors( label )

    switch label
        case 'silence'
            c = [0 0 0];            
        case 'stayStill'
            c = [0 0 0];
            
        % Phonemes
        case 'i'
            c = [166 206 227];
        case 'u'
            c = [31 120 180];
        case 'ae'
            c = [178 223 138];
        case 'a'
            c = [51 160 44];
        case 'ba'
            c = [251 154 153];
        case 'ga'
            c = [227 26 28];
        case 'da'
            c = [253 191 111];
        case 'k'
            c = [255 127 0];
        case 'p'
            c = [202 178 214];
        case 'sh'
            c = [106 21 154];
            
        % Words     
        case 'beet'
            c = [141 211 199];
        case 'bat'
            c = [229 222 91];
        case 'bot'
            c = [190 186 218];
        case 'boot'
            c = [251 128 114];
        case 'dot'
            c = [128 177 211];
        case 'got'
            c = [253 180 98];
        case 'shot'
            c = [179 222 105];
        case 'keep'
            c = [252 205 229];
        case 'seal'
            c = [150 150 150];
        case 'more'
            c =[188 128 189];
            
        % Movements
        case 'tongueLeft'
            c = [230 171 2];
        case 'tongueRight'
            c = [166 118 29];
        case 'tongueDown'
            c = [27 158 119];
        case 'tongueUp'
            c = [102 166 30];
        case 'lipsForward'
           c = [117 112 179];
        case 'lipsBack'
            c = [231 41 138];
        case 'mouthOpen'
            c = [217 95 2];
        otherwise
            error('Color %s is not in speechColors lookup', label )
    end
    
    c = c ./ 255;
end