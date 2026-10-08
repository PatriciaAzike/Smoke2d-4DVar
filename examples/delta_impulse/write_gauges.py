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



# The gauge location must not sit on a patch boundary. (1.5, 1.0) is a patch
# corner at every refinement level. 1.5 and 1.0 are both exact multiples of the
# patch width 2/2^L for every L and ForestClaw's point-location search then
# fails to assign the gauge to any patch. The gauge records nothing and
# get_gauges() aborts with "No gauge buffer available".
#
# None of the original three gauges hit this: (0.67, 1.0) is interior in x,
# (0.80, 1.2) and (0.55, 1.2) are interior in both. This only appeared when the
# impulse moved to a dyadic point.
#
# So offset the gauge by a non-dyadic amount. The delta source stays at the
# (xm, ym) from the .ini. It is evaluated at cell centres and needs no
# point location. 0.001 is 1% of sigma = 0.1414, far below anything that matters.
EPS = 0.001

# The window also ends at tm rather than being zero-width, because get_gauges()
# reads gauge_buffer[kmax-1], the last sample. A window of [tm-dt, tm] gives it
# several samples and makes the last one the sample closest to tm (at or before it).
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
