import configparser
import fclaw_analysis

def parse_array(value):
    return [float(v) for v in value.split()]

config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
config.read("adjoint_options.ini")   # use your actual ini filename

mdata = config.getint("adjoint-user", "mdata")

xm = parse_array(config.get("adjoint-user", "xm"))
ym = parse_array(config.get("adjoint-user", "ym"))
tm = parse_array(config.get("adjoint-user", "tm"))

dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)



# Offset the gauge from (xm, ym) = (1.5, 1.0). That point is a patch corner 
# at every refinement level used here (2 and up). As a result, ForestClaw's 
# point-location search assigns the gauge to no patch and it records nothing.
# The original gauges, e.g. (0.67, 1.0), were not corners. The delta source
# itself stays at the (xm, ym) from the .ini. It is evaluated at cell centres
# and needs no point search. We choose EPS = 0.001, which is under 1% of the
# kernel width sigma = 0.1414.
EPS = 0.001

# Time window [tm - dt, tm], not zero-width as in adjoint_forward. With [tm, tm]
# the gauge buffer was empty. The upper end misses the step at tm by round-off
# (the clock reaches 1.5000000000000002), so the last sample is about one step
# before tm and the model records none.  That's fine here: the delta_impulse
# build never reads the gauges; the verification uses the output frames.
dt = 0.02  # time window

for i in range(mdata):

    gaugedata.gauges.append([
        i,
        xm[i] - EPS,
        ym[i] + EPS,
        tm[i] - dt,
        tm[i]
    ])



gaugedata.write(data_source='write_gauges.py')

# Each run reads its own copy from its run directory; keep them identical.
import os, shutil
for d in ("adjoint", "forward", "model"):
    os.makedirs(d, exist_ok=True)
    shutil.copyfile("gauges.data", os.path.join(d, "gauges.data"))
