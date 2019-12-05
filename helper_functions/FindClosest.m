% FindClosest.m
%
% Returns the index in vector u that is closest to value val, along with
% the difference ebtween val and u(ind). Useful when searching in a dataseries
% when an exact numerical match may not be expected (for instance a specific
% time in a time series).
% Optional arguments allow you to enforce that the matching element in u
% be {>, >=, <, <=} than val.
%
% USAGE: [ closestVal, ind, diff ] = FindClosest( u, val, relationOp )
%
% EXAMPLE:
%
% INPUTS:
%     u                         vector of values to check
%     val                       Value you want to (most closely) find in val
%     (relationOp)              (optional)  Can be '<=', '>=', '<', '>' to specify that you
%                               require the value in <u> to be, for example,
%                               greater than <val>.
% OUTPUTS:
%     closestVal                closest matching value found
%     ind                       the index in <u> where the match was found
%     diff                      closestVal - val
%                       
% Created by Sergey Stavisky on 16 Mar 2012
function [closestVal, ind, diff] = FindClosest( u, val, relationOp )

    % Input processing / defaults
    if ~isvector( u )
        error( 'input <u> must be a vector' )
    end
    if ~isscalar( val )
        error( 'input <val> must be a vector' )
    end
    
    if nargin < 3
        relationOp = [];
    else
        if ~any( strcmp( relationOp, {'<=', '>=', '<', '>'}) )
            error( 'relationOp input must be one of ''>'', ''<'', ''>='', ''<=''. You entered ''%s''', ...
                relationOp )
        end 
    end 

    % Ignore certain values
    if ~isempty( relationOp )
        checkInds = true( numel(u) );
        
        switch relationOp
            case '>'
                checkInds = u > val;
            case '>='
                checkInds = u >= val;
            case '<'
                checkInds = u < val;
            case '<='
                checkInds = u <= val;
        end
        u(~checkInds) = NaN;
    end % if ~isempty( relationOp )
    %%
    [~, ind] =  min( abs( val - u ) ); 
    closestVal = u(ind);
    diff = closestVal - val;
end