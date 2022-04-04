import numpy as np 
import matplotlib.pyplot as plt 
import stress 
import plotting 

def dsin(t):
    t = np.deg2rad(t)
    return np.sin(t)

def dcos(t):
    t = np.deg2rad(t)
    return np.cos(t)


def rgw(azi, inc):
    r = np.array([[-dcos(azi)*dcos(inc), -dsin(azi)*dcos(inc), dsin(inc)],
                   [dsin(azi), -dcos(azi), 0],
                   [dcos(azi)*dsin(inc), dsin(azi)*dsin(inc), dcos(inc)]])
    return r 

def cylindrical2cartesian(r, theta):
    vecs = np.zeros((3, len(theta)))
    vecs[0,:] = r*dcos(theta)
    vecs[1,:] = r*dsin(theta)
    return vecs

# stress state
Sv = 88.5 
SH = 97.3
Sh = 59.5
SH_azi = np.linspace(0, 180, 180)

# hole geometry
inc = 45
azi = 180
Pp = 35.5
Pw = 60
pr = 0.25
mu = 0.6 
ucs = 40 

n_theta = 360 # resolusi sampel di dinding bor (azimuthal)
theta = np.linspace(0, 360, n_theta)
hoop_stress_grid = np.zeros((len(SH_azi), n_theta))

for idx, sh_azimuth in enumerate(SH_azi):
    # stress state
    ss = stress.StressState(Sv, SH, Sh, sh_azimuth)

    # create well 
    well1 = stress.Wellbore(ss, inc=inc, azi=azi)

    # get nearfield stress around borehole
    r = 1
    _, stt, _, _, _, _ = well1.wall_stress(r, theta, Pp, Pw, pr, principal=False)

    # create rotation matrix
    rotmat = rgw(azi, inc)

    # convert angle to vector components in wellbore coordinate
    vecs = cylindrical2cartesian(1, theta)

    # convert vector from wellbore coordinate to NED geographical coordinate
    vecs_geo = rotmat @ vecs 
    azimuth_geo = np.rad2deg(np.arctan2(vecs_geo[1,:],vecs_geo[0,:]))
    negs = np.where(azimuth_geo < 0)[0]
    azimuth_geo[negs] = 360 + azimuth_geo[negs]

    # sort based on azimuth
    stt_geo = stt[azimuth_geo.argsort()] 
    azimuth_geo_sorted = azimuth_geo[azimuth_geo.argsort()]

    hoop_stress_grid[idx, :] = stt_geo

plt.title(f'Well inclination : {inc} degree, well azimuth : {azi} (degree from North)', pad=20)
plt.pcolormesh(hoop_stress_grid, cmap='jet')
plt.gca().invert_yaxis()

plt.xlabel('Azimuth (degree from North)')
plt.ylabel('SHmax azimuth (degree from North)')
plt.colorbar(label='Hoop Stress (MPa)')
plt.show()



