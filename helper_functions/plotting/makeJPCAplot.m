% Makes the jPCA plot using Mark's code, and labels each condition with its name.
% Sergey Stavisky, January 5 2018
% Stanford Neural Prosthetics Translational Laboratory   
function [figh, cmap] = makeJPCAplot( Projection, Summary, plotParams, uniqueLabels )
    [colorStruct, haxP, vaxP] = phaseSpace( Projection, Summary, plotParams );

    % identify the groups
    [ indsLeftToRight, cmap ] = whichGroupIsWhichJpca( Projection );
    
    % Label the end points of each condition with its label. Note this won't work if more than
    % one plane was plotted (because it'll try to plot on the last plane using data from first
    % plane. So make plotParams.planes2plot = 1 to use this .
    figh = gcf;
    axh = figh.Children;
    axh.Children;
    

    for iLabel = 1 : numel( uniqueLabels )
        myH = Projection(iLabel).proj(end,1);
        myV = Projection(iLabel).proj(end,2);
        th(iLabel) = text( myH, myV, uniqueLabels{iLabel} );        
    end

end