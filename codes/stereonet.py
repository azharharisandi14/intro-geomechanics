import numpy as np
import matplotlib.pyplot as plt
import stress
import plotting

# Sv = 88.5
# SH = 97.3
# Sh = 59.5
# SH_azi = 168

Sv = 88.2
SH = 90
Sh = 51.5
SH_azi = 90

# Sv = 5483
# SH = 6425
# Sh = 4681
# SH_azi = 255

# stress state
ss = stress.StressState(Sv, SH, Sh, SH_azi)

# data
Pp = 31.5
Pw = 31.5
pr = 0.25
mu = 1
ucs = 45

# Pp = 2558
# Pw = 2670
# pr = 0.177
# mu = 1
# ucs = 6298

# incs = 90
# azis = np.linspace(0, 360, 180)
# azis = [90, 180]
# # debug
# plt.figure()
# for azi in azis:
#     t = np.linspace(0, 360, 500)
#     well = stress.Wellbore(ss, incs,azi)
#     s1, s3 = well.wall_stress(1, t, Pp, Pw, pr, principal = True)
#     plt.axhline(ucs)
#     wbo, _ = well.breakout_width(Pp, Pw, pr, mu, ucs)
#     print(f'wbo = {wbo}:.2f degree')
#     plt.plot(t, s1, label=f's1 azi = {azi}')
#     plt.plot(t, s3, label=f's3 azi = {azi}')
#     plt.legend()
# plt.grid()

# bikin matriks untuk azimuth dan inklinasi
azi = np.linspace(0, 360, 180)
inc = np.linspace(0, 90, 60)

# bikin 2d grid
azis, incs = np.meshgrid(azi, inc)


# bikin 2d grid kosong
co = np.zeros_like(azis)
wbo = np.zeros_like(azis)

for i in range(len(inc)):
    for j in range(len(azi)):
        well = stress.Wellbore(ss, incs[i, j], azis[i, j])
        wbo[i, j], co[i, j] = well.breakout_width(Pp, Pw, pr, mu, ucs)


plotting._plot_polar(azis, incs, wbo, 'breakout width (degree)', donut=False, vmin=0, vmax=100)
plotting._plot_polar(azis, incs, co, 'Required Co (MPa)', donut=False, vmin=0, vmax=200)

print(np.amin(wbo), np.amax(wbo))
print(np.amin(co), np.amax(co))

plt.show()


