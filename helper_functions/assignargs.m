function par = assignargs(varargin)
% par = assignargs(defaults, varargin)
%
%   Like structargs except additionally assigns values individually by their names in
%   caller.
%
%   Overwrites fields in struct defaults with those specified by:
%    - if arg(1) is a structure, the values therein
%    - the values specified in 'name', value pairs in the arguments list
%      (these values take precedence over arg(1)
%   Assigns new values or old defaults in caller workspace
%
% par = structargs(varargin)
%
%   Same functionality as above, except uses all existing variables in the 
%   calling workspace as defaults
%
% Author: Dan O'Shea (dan@djoshea.com), (c) 2008


par = structargs(varargin{:});

fields = fieldnames(par);
for i = 1:length(fields)
    assignin('caller', fields{i}, par.(fields{i}))
end

end

