% ChannelNameToNumber.m
%
% Converts from string specifying array/channel to a unique number (which can
% be over 96 if on the second array). Its inverse is chanNumToName.m
%
% USAGE: [ chanNum ] = ChannelNameToNumber( channelname )
%
% EXAMPLE: ChannelNameToNumber('chan_2.3')  returns  99
%
% INPUTS:
%     channelname               string of format 'chan_ARRAY.CHANNEL'
%
% OUTPUTS:
%     chanNum                   channel number
%
% Created by Sergey Stavisky on 08 May 2012

function [ chanNum ] = ChannelNameToNumber( channelname, varargin )

    def.MAXCHANSPERARRAY = 96;
    assignargs( def, varargin );

    if ~iscell( channelname )
        channelname = {channelname};
    end
    
    for i = 1 : numel( channelname )
        if strfind( channelname{i}, 'chan_' )
            % assumes no more than 9 arrays
            dotind = regexp(channelname{i}, '[0-9].', 'end');
            arrayN = str2double( channelname{i}(dotind-1:dotind-1) );
            chanN  = str2double( channelname{i}(dotind+1:end) );
            chanNum(i) = (arrayN-1)*MAXCHANSPERARRAY + chanN;
            
        else
            error('input %s did not have keyword ''chan_'' in it.', ...
                channelname{i} )
        end
    end


end