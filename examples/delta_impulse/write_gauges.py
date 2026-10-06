import configparser
import fclaw_analysis

def parse_array(value):
    return [float(v) for v in value.split()]

config = configparser.ConfigParser()
config.read("adjoint_options.ini")   # use your actual ini filename

mdata = config.getint("adjoint-user", "mdata")

xm = parse_array(config.get("adjoint-user", "xm"))
ym = parse_array(config.get("adjoint-user", "ym"))
tm = parse_array(config.get("adjoint-user", "tm"))

dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)



# The gauge location must NOT sit on a patch boundary. (1.5, 1.0) is a patch
# CORNER at every refinement level -- 1.5 and 1.0 are both exact multiples of the
# patch width 2/2^L for every L -- and ForestClaw's point-location search then
# fails to assign the gauge to any patch. The gauge records nothing and
# get_gauges() aborts with "No gauge buffer available".
#
# None of the original three gauges hit this: (0.67, 1.0) is interior in x,
# (0.80, 1.2) and (0.55, 1.2) are interior in both. This only appeared when the
# impulse moved to a dyadic point.
#
# So offset the GAUGE by a non-dyadic amount. The delta SOURCE stays at the
# clean (xm, ym) from the .ini -- it is evaluated at cell centres and needs no
# point location. 0.001 is 1% of sigma = 0.1414, far below anything that matters.
EPS = 0.001

# The window also ends AT tm rather than being zero-width, because get_gauges()
# reads gauge_buffer[kmax-1], the last sample. A window of [tm-dt, tm] gives it
# several samples and makes the last one the closest at or before tm.
TWIN = 0.02 #time window

for i in range(mdata):

    gaugedata.gauges.append([
        i,
        xm[i] - EPS,
        ym[i] + EPS,
        tm[i] - TWIN,
        tm[i]
    ])



gaugedata.write(data_source='write_gauges.py')
