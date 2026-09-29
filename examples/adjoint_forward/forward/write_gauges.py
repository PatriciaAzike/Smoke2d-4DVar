import configparser
import fclaw_analysis

# Define observations here — single source of truth
observations = [
    {"xm": 0.67, "ym": 1.0, "tm": 0.4, "dm": 3.0},
    {"xm": 0.80, "ym": 1.2, "tm": 0.4, "dm": 3.0},
    {"xm": 0.55, "ym": 1.2, "tm": 0.4, "dm": 2.0},
    #{"xm": 0.6, "ym": 0.35, "tm": 0.4, "dm": 1.0},
]

mdata = len(observations)
xm = [o["xm"] for o in observations]
ym = [o["ym"] for o in observations]
tm = [o["tm"] for o in observations]
dm = [o["dm"] for o in observations]

# Write to adjoint_options.ini
config = configparser.ConfigParser()
config.read("../adjoint_options.ini")

config["adjoint-user"]["mdata"] = str(mdata)
config["adjoint-user"]["xm"]    = " ".join(str(v) for v in xm)
config["adjoint-user"]["ym"]    = " ".join(str(v) for v in ym)
config["adjoint-user"]["tm"]    = " ".join(str(v) for v in tm)
config["adjoint-user"]["dm"]    = " ".join(str(v) for v in dm)

with open("../adjoint_options.ini", "w") as f:
    config.write(f)

# Set up gauge data
dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)

# Periodic domain : [-1,1]x[-1,1]
#
# Format : [id, x, y, t0, t1]
#   id   : integer, identifying the gauge
#  x,y   : Location of the gauge
#  t0,t1 : (t0,t1) interval over which to monitor the gauge.


# xm = 0.67
# ym = 1.0
# zm = 0
# t0 = 0.4
# t1 = 0.4
#gaugedata.gauges.append([  0, xm, ym, t0, t1])

#gaugedata.gauges.append([  1, xm, ym, t0,  t1])

for i in range(mdata):

    gaugedata.gauges.append([
        i,
        xm[i],
        ym[i],
        tm[i],
        tm[i]
    ])



gaugedata.write(data_source='write_gauges.py')
