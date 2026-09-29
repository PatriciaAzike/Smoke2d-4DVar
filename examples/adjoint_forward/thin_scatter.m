function n = thin_scatter(h)
%THIN_SCATTER  Collapse a ClawPack PlotType==4 scatter to its distinct points.
%
%   hclaw = findobj(gca, '-property', 'XData');   % capture BEFORE adding more
%   ... plot q_F, the gauges, the legend ...
%   thin_scatter(hclaw);
%
% ClawPack draws one graphics object per AMR patch, so frame 4 of the pseudo-1D
% run puts 188,416 markers on the axes. That run is constant in y, so only 416
% of them are at distinct (x,q) positions; the other ~188,000 sit exactly on top
% of each other. Painters still rasterizes every one, and with
% ScatterStyle = {'*'} each marker is three line segments, so print at -r300 has
% to emit on the order of 5e5 vector primitives. That is the "Busy" hang.
%
% This keeps the distinct points and blanks the duplicates, so the picture is
% pixel-for-pixel identical and the export is instant. Objects are emptied,
% never deleted, so the legend and any handle ClawPack is holding stay valid.
%
% Objects are grouped by colour and marker before thinning, so if ClawPack ever
% colours refinement levels differently that distinction survives.

n = 0;
if nargin < 1 || isempty(h), return; end
h = h(ishghandle(h));
h = h(arrayfun(@(g) isprop(g,'XData') && isprop(g,'YData'), h));
if isempty(h), return; end

key = cell(numel(h), 1);
for k = 1:numel(h)
    c = [0 0 0];
    if isprop(h(k), 'Color'), c = get(h(k), 'Color'); end
    m = get(h(k), 'Marker');
    if ~ischar(m), m = '?'; end
    key{k} = sprintf('%.4f_%.4f_%.4f_%s', c, m);
end
[~, ~, g] = unique(key);

before = 0;
for gi = 1:max(g)
    hg = h(g == gi);
    x = []; y = [];
    for k = 1:numel(hg)
        x = [x; get(hg(k),'XData').'];  y = [y; get(hg(k),'YData').'];   %#ok<AGROW>
    end
    if isempty(x), continue; end
    before = before + numel(x);

    [~, ia] = unique([x y], 'rows', 'stable');
    set(hg(1), 'XData', x(ia), 'YData', y(ia));
    for k = 2:numel(hg)
        set(hg(k), 'XData', [], 'YData', []);
    end
    n = n + numel(ia);
end

fprintf('thin_scatter: %d markers -> %d distinct (%d objects)\n', ...
        before, n, numel(h));
end
