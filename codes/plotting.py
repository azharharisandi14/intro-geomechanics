import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import numpy as np

from stress import Wellbore

def wall(well, Pp, Pw, pr, mu, label, contour=None, **kwargs):
    q = (np.sqrt(mu ** 2 + 1) + mu) ** 2
    r = np.linspace(1, 1.5, 40)
    t = np.linspace(0, 360, 180)
    r, t = np.meshgrid(r, t)
    s1, _, s3 = well.wall_stress(r, t, Pp, Pw, pr, principal=True)

    t = np.deg2rad(t)
    ucs = s1 - q * s3

    _plot_polar(t, r, ucs, label, contour, **kwargs)
    #plt.title('Required Co')
    #plt.show()

def _plot_polar(t, r, z, label, contour=None, donut=True, **kwargs):
    f, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(4, 4))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    if donut:
        ax.set_ylim(0, np.amax(r))
        cont = ax.contourf(t, r, z)
        a = ax.pcolormesh(t, r, z, shading='auto', 
                          **kwargs)
        ax.plot((np.pi, np.pi), (0.85, 1), c='k', linewidth=0.8)
        ax.scatter(0, 0, marker='x', s=40, c='k')
        ax.grid(color='k', linewidth=0.5)

        # new add
        ax.contour(cont, levels=contour, colors='black', linestyles='dashed')
        ax.set_rgrids((0.99, 1.5))
        ax.set_thetagrids([])
        # no radial label
        ax.set_yticklabels([])
        # ax.set_xticks([np.pi])
        # ax.set_xticklabels(['Bottom'])
    else:
        cont = ax.contourf(np.deg2rad(t), r, z)
        a = ax.pcolormesh(np.deg2rad(t), r, z, shading='auto', cmap='jet', **kwargs)
        ax.grid(color='k', linewidth=0.5)
        ax.contour(cont, levels=contour, colors='black', linestyles='dashed')
        ax.set_rgrids((30, 60), angle=90)
        ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2])
        ax.set_xticklabels(['N', 'E', 'S', 'W'])
        ax.set_thetagrids(range(0, 360, 30))

    # ticks location in colorbar
    ticks_cbar = [np.amin(z) + (np.amax(z) - np.amin(z)) / 4 * (i) for i in range(5)]

    cb = f.colorbar(a, fraction=0.04, pad=0.1,
                    orientation="horizontal", ticks=ticks_cbar)
    cb.set_label(label)
