% MakeDumbLegend.m
%
% I've found the default matlab legend function not well-suited to
% complicated plots where there are many lineseries. So I wanted one that
% was much simpler; just give it strings and colors and it puts it up 
% somewhere.
% I'm going to keep expanding this out to meet growing feature needs.
% Current features:
%      Scales text size to fit within the window
%      Scalar expansion of color, or use the rows of the 'Color' input
%      to color each line of text accordingly.
%
% Part of the Generalized Figures Function Set
%
%
% USAGE: [ texth ] = makeDumbLegend( labels )
%
% EXAMPLE: texth = makeDumbLegend( 'LFP', 'spiking' )
%
% INPUTS:
%     labels                    cell array of string labels for the categories
%     OPTIONAL ARGUMENTS SET IN PARAMETER-VALUE PAIR FORMAT
%     (Color)                   RGB color matrix
%     (Position)                2x1 vector: distance from left/bottom in normalized (0 to 1)
%                               units
%     (BackgroundColor)         Color of the box
%     (LineStyles)              Cell list of line styles for each entry.
%     (FontSize)                Font size for text.
%     (Marker)                  Cell list of the marker style (if you want a marker) 
%                               for each entry. 
%
% OUTPUTS:
%     texth                     
%
% Created by Sergey Stavisky on 10 Mar 2012
% Last modified by Sergey Stavisky on 10 Mar 2012

function [ texth ] = MakeDumbLegend( labels, varargin )

    % Argument Processing / Defaults

    if ~iscell( labels )
        % if labels isn't a cell array, make it one
        labels = {labels};
    end
    
    N = numel( labels ); % number of labels

    def.Color = ones( N, 3 );     % Silly to use default colors but I it can be used
                                  % just for annotation this way.
                                  

    def.axish = gca;
    def.Position = [0.05 0.95];
    def.FontSize = 14;
    def.BackgroundColor = [];
    def.LineStyles = []; % by default lines won't be plotted.
    def.Marker = []; % by default markers won't be plotted.
    def.MarkerSize = 15^2;
    def.EdgeColor = 'none';
    def.HorizontalAlignment = 'left';
    def.VerticalAlignment = 'top';
    
    assignargs( def, varargin );
    
    
    % Input processing.
    if ~isempty( LineStyles )
        if ~iscell( LineStyles )
            LineStyles = {LineStyles};
        end
        if numel( LineStyles ) == 1
            LineStyles = repmat( LineStyles, numel( labels ), 1 );
        end
    end
   
    % Because I'm using multiline text, I need to enter the RGB as a stirng
    % into each cell. to get the colors to show up. If only one color was entered,
    % then repmat it to so each row gets the same color.
    if numel( Color ) == 3 %#ok<NODEF>
        Color = repmat( forceRow( Color ), N, 1 );
    end
    for i = 1 : N
        labels{i} = sprintf('\\color[rgb]{%f %f %f}%s', ...
            Color(i,1), Color(i,2), Color(i,3), labels{i});
        % Since I'm using tex, make underscores not do subscript by 
        % putting an escape character before them
        labels{i} = regexprep( labels{i}, '_', '\\_' );
    end
    
    
    % Write the text 
    texth = text( Position(1), Position(2), labels, ...
        'HorizontalAlignment', HorizontalAlignment, 'VerticalAlignment', VerticalAlignment', ...
        'Units', 'normalized', 'FontSize', FontSize, 'EdgeColor', EdgeColor);
    
    % Line Styles, if so specified.
    if ~isempty( LineStyles )
        % Aim the lines evenly along the extent of the text
        set(texth, 'Units', 'data'); % need this or line will be out in nowhere
        myPos = get( texth, 'Extent' );
        top = myPos(2) + myPos(4);
        spacing = myPos(4) / (2*N);
        xends = [myPos(1)-0.25*myPos(3) myPos(1)-0.05*myPos(3)]; % scale it to total size of textbox
        for i = 1 : N
            myY =  top - spacing*(1+2*(i-1));            
            lineh = line( xends, [myY myY], 'Color', Color(i,:), 'LineStyle', LineStyles{i}, ...
                'LineWidth', 3);    
        end
    end %if ~isempty( LineStyles )
    
    % Marker, if so specified.
    if ~isempty( Marker )
        set(texth, 'Units', 'data'); % need this or line will be out in nowhere
        myPos = get( texth, 'Extent' );
        top = myPos(2) + myPos(4);
        spacing = myPos(4) / (2*N);
        myX = myPos(1)-0.33*myPos(3); % scale it to total size of textbox
        for i = 1 : N
            myY =  top - spacing*(1+2*(i-1));
            markerh = scatter( myX, myY, MarkerSize, ...
                'MarkerEdgeColor', Color(i,:), 'Marker', Marker{i}); %#ok<NASGU>
        end        
    end %if ~isempty( Marker )
    
    
    
    % Background color if so specified
    if ~isempty( BackgroundColor )
        set( texth, 'BackgroundColor', BackgroundColor );
    end
    
    % Parameter-value setting of other properties specified as optional args
    for i = 1 : 2 : numel( varargin )
        switch varargin{i}
            % some properties can't be set or they'll breka it; enumerate them here
            case 'Position'
                continue
            otherwise
                try
                    set( texth, varargin{i}, varargin{i+1} );
                catch %#ok<CTCH>
                    
                end
        end
    end
    % make text fit - use try catch since it sometimes fails
    try
        texth = DecreaseFontUntilTextFits( texth );
    catch
        fprintf( 'DecreseFontUntilTextFits failed. Moving on.\n')
    end
end


function texth = DecreaseFontUntilTextFits( texth )
    % Decreases the fontsize of text object <texth> until its extent does not
    % exceed the x:[0 1], y:[0 1] boundaries.
    itFits = false;
    % the code relies on units being normalized
    set( texth, 'Units', 'normalized');
    while ~itFits
        currExtent = get( texth, 'Extent' );
        currFont = get( texth, 'FontSize' );
        itFits = ~( currExtent(1) < 0 | currExtent(1)+currExtent(3) > 1 ...
            | currExtent(2) < 0 | currExtent(2)+currExtent(4) > 1 );
        
        if currFont <= 1
            cprintf('red', '[%s] Could not make text fit even with size 1 font.', ...
                mfilename)
            return
        end
        
        if ~itFits 
            set( texth, 'FontSize', currFont - 1 );
        end 
    end %while ~itFits
end  %function texth = DecreaseFontUntilTextFits( texth )