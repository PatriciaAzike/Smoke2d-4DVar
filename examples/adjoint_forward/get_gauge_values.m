function values = get_gauge_values(amrdata, xg, yg, mq)
% GET_GAUGE_VALUES : q at the gauges (xg(k), yg(k)) in the frame being plotted.
%
%    values = get_gauge_values(amrdata, xg, yg)
%    values = get_gauge_values(amrdata, xg, yg, mq)
%
% amrdata is the frame plotclaw2 has already read. Uses ForestClaw's gauge
% interpolation (fclaw2d_clawpatch46_fort_gauges_update). That is, in the
% finest patch containing the gauge, take the cell (i,j) holding it and
% interpolate bilinearly with cells i+1 and j+1, using the offsets of the
% gauge from the centre of (i,j). mq is the component of q, which is 1 by 
% default.

if nargin < 4
    mq = 1;
end
amr = amrdata;

values = nan(numel(xg), 1);
for k = 1:numel(xg)
    values(k) = interpolate_at_gauge(amr, xg(k), yg(k), mq);
end

end


function v = interpolate_at_gauge(amr, x, y, mq)
p = find_leaf(amr, x, y);
P = amr(p);

% Cell (i,j) holding the gauge and the offsets of the gauge from its centre,
% in units of the cell size
i = floor((x - P.xlow)/P.dx);
j = floor((y - P.ylow)/P.dy);
xo = (x - (P.xlow + (i + 0.5)*P.dx))/P.dx;
yo = (y - (P.ylow + (j + 0.5)*P.dy))/P.dy;

v = (1 - xo)*(1 - yo)*cell_value(amr, p, i,     j,     mq) ...
  +      xo *(1 - yo)*cell_value(amr, p, i + 1, j,     mq) ...
  + (1 - xo)*     yo *cell_value(amr, p, i,     j + 1, mq) ...
  +      xo *     yo *cell_value(amr, p, i + 1, j + 1, mq);
end


function p = find_leaf(amr, x, y)
% The finest patch whose area contains (x,y)
p = 0;
dxmin = Inf;
for n = 1:numel(amr)
    P = amr(n);
    inside = x >= P.xlow && x < P.xlow + P.mx*P.dx && ...
             y >= P.ylow && y < P.ylow + P.my*P.dy;
    if inside && P.dx < dxmin
        p = n;
        dxmin = P.dx;
    end
end
end


function v = cell_value(amr, p, i, j, mq)
% q in cell (i,j) (0-based) of patch p. A cell one past the patch edge stands
% in for ForestClaw's ghost cell. Use the leaf cell that contains its centre.
P = amr(p);
if i >= 0 && i < P.mx && j >= 0 && j < P.my
    v = P.data(mq, 1 + i + j*P.mx);
    return
end
L = 2;                                  % the domain is [0,L]^2, periodic
x = mod(P.xlow + (i + 0.5)*P.dx, L);
y = mod(P.ylow + (j + 0.5)*P.dy, L);
n = find_leaf(amr, x, y);
N = amr(n);
ii = floor((x - N.xlow)/N.dx);
jj = floor((y - N.ylow)/N.dy);
v = N.data(mq, 1 + ii + jj*N.mx);
end
