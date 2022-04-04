import numpy as np 
import matplotlib.pyplot as plt 
import stress 

# stress state
Sv = 88.2
SH = 90
Sh = 51.5
SH_azi = 90

# hole geometry
inc = 30
azi = 45
Pp = 31.5
Pw = 36.5
pr = 0.25
mu = 0.6
ucs = 40 

r = 1
theta = np.linspace(0, 360, 360)
ss = stress.StressState(Sv, SH, Sh, SH_azi)
well = stress.Wellbore(ss, inc=inc, azi=azi)
s1, s2, s3 = well.wall_stress(r, theta, Pp, Pw, pr, principal=True)
srr, stt, sz, trt, ttz, trz = well.wall_stress(r, theta, Pp, Pw, pr, principal=False)

# plt.plot(theta, s1, label='s1')
plt.plot(theta, stt, label='stt')
plt.plot(theta, srr, label='srr')
plt.plot(theta, sz, label='sz')
plt.plot(theta, ttz, label='ttz')
plt.plot(theta, trt, label='trt')
plt.plot(theta, trz, label='trz')
# plt.plot(theta, s1, '*-', label='s1')
plt.legend()
plt.show()
