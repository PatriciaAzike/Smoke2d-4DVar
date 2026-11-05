function color_gauges(hdl, values)
% COLOR_GAUGES  Fill each gauge marker with the colour of its value.
%
%    color_gauges(hdl, values)
%
% hdl are the gauge markers (from add_gauges) and values(k) the field at gauge
% k (from get_gauge_values). The fill uses the current colormap and colour
% limits, so a gauge can be read against the field around it; the black edge
% marks its location. Values outside the limits take the end colours.

cmap = colormap(gca);
cl = clim(gca);
n = size(cmap, 1);

for k = 1:numel(hdl)
    % The colormap row Matlab uses for a value in a scaled CData plot,
    % so the marker has exactly the colour the field would have there
    idx = fix((values(k) - cl(1))/(cl(2) - cl(1))*n) + 1;
    idx = min(max(idx, 1), n);
    c = cmap(idx, :);
    set(hdl(k), 'Marker', 'o', 'MarkerFaceColor', c, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1);
end

end
