import numpy as np
import matplotlib.pyplot as plt
import stress
import plotting

Sv = 88.2
SH = 90
Sh = 51.5
SH_azi = 90

# stress state
ss = stress.StressState(Sv, SH, Sh, SH_azi)

# hole geometry
inc = 40 # inclination
azi = 0 # azimuth
Pp = 31.5
Pw = 31.5
pr = 0.25
mu = 1
ucs = 45

# create well instance
well1 = stress.Wellbore(ss, inc=inc, azi=azi)

# donut plot (required rock strength)
plotting.wall(well1, Pp, Pw, pr, mu, 'Required Co (MPa)', contour=[ucs], cmap='jet')
plt.show()

