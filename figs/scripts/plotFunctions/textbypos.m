function [xpos, ypos] = textbypos(posx, posy, txtstr, varargin)
% -------------------------------------------------------------------------
%TEXTBYPOS makes text box based off normalized axes position
% -------------------------------------------------------------------------
% Syntax: 
%   [xpos, ypos] = textbypos(posx, posy, txtstr, varargin)
% 
% Inputs:
%   posx
%   posy        y position in normalized figure coordinates
%   txtstr      text string
%   varargin    passed to text varargin
% 
% Sample:
%   figure(5); clf;
%   x = [1:100]; y = exp([1:100]*0.05 + randn(size(x))./4);
%   plot(x,y, 'ko-')
%   textbypos(0.05, 0.8, 'textbypos(0.05, 0.8)', 'Color', 'b', 'BackgroundColor', 'w', 'EdgeColor', 'k');
% -------------------------------------------------------------------------

pos = get(gca, 'position');
xlims = get(gca, 'xlim' );
ylims = get(gca, 'ylim');

xsc = get(gca, 'xscale' );
ysc = get(gca, 'yscale');

[xmin xmax] = feval(@(x) x{:}, num2cell(xlims));
[ymin ymax] = feval(@(x) x{:}, num2cell(ylims));

if strcmp(ysc, 'linear') & strcmp(ysc, 'linear')

    xpos = (posx-pos(1))/pos(3) * (range(xlims)) + xmin;
    ypos = (posy-pos(2))/pos(4) * (range(ylims)) + ymin;

elseif strcmp(ysc, 'log') & strcmp(ysc, 'log')

    xpos = 10.^((posx-pos(1))/pos(3) * range(log10(xlims)) + log10(xmin));
    ypos = 10.^((posy-pos(2))/pos(4) * range(log10(ylims)) + log10(ymin));

end
text(xpos, ypos, txtstr, varargin{:})


end

%%%



 