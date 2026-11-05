function style_panel(cb, cticks)
% STYLE_PANEL: Common geometry for every 2D field panel in the paper.
%
%   style_panel()              % no colorbar
%   style_panel(cb)            % colorbar handle, automatic ticks
%   style_panel(cb, 0:4)       % colorbar handle, whole-number ticks
%
% Produces a 990 x 910 px canvas at 300 dpi, so panels need little or no 
% trim in LaTeX and all figures share an identical plot rectangle.
%


W = 3.30;  H = 3.0333;                 % inches -> 990 x 910 px at 300 dpi

set(gcf, 'Color','w', 'Units','inches', 'Position',[1 1 W H]);
set(gcf, 'PaperPositionMode','auto');

set(gca, 'FontSize', 16);
set(gca, 'XTick', 0:2, 'YTick', 0:2);
set(gca, 'TickLabelInterpreter', 'latex');      
set(cb, 'TickLabelInterpreter', 'latex');

xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
daspect([1 1 1]);
set(gca, 'Units','normalized', 'Position',[0.165 0.205 0.610 0.665]);

if nargin >= 1 && ~isempty(cb) && isvalid(cb)
    if nargin >= 2 && ~isempty(cticks)
        cb.Ticks = cticks;
    end
    set(cb, 'Units','normalized', 'Position',[0.815 0.205 0.038 0.665]);
end
end
