function style_line_panel()
%STYLE_LINE_PANEL  Geometry for wide line/scatter panels.
%
% Companion to STYLE_PANEL.
%


W = 6.60;  H = 3.0333;                 % inches -> 1980 x 910 px at 300 dpi

set(gcf, 'Color','w', 'Units','inches', 'Position',[1 1 W H]);
set(gcf, 'PaperPositionMode','auto');

set(gca, 'FontSize', 16, 'Box', 'on');
set(gca, 'Units','normalized', 'Position',[0.165 0.205 0.790 0.665]);
end
