function n = thin_scatter(h)
% THIN_SCATTER: Collapse a ClawPack PlotType==4 scatter to its distinct points by
% removing duplicate (x,q) positions to prevent a  "Busy" hang during export.
%
% The duplicate objects are emptied without deleting handles, preserving legends.  
% This also reduces export time without changing the plotted points
%
%   hclaw = findobj(gca, '-property', 'XData');   % capture handles adding overlays
%   thin_scatter(hclaw);


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
