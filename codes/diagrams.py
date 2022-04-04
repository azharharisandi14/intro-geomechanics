import matplotlib.pyplot as plt 

def _get_center_and_radius(smax:float, smin:float) -> tuple:
    """get center and radius of a Mohr's circle
    given two points that are separated by 180 degree
    along a circle

    Args:
        smax (float): maximum principal stress
        smin (float): minimum principal stress

    Returns:
        tuple: (center, radius)
    """
    c = (smax + smin)/2
    r = (smax - smin)/2
    return c, r

def mohr2d(sigma1:float, sigma3:float, ax, color='k', **kwargs):
    """create mohr diagram for 2d stress state

    Args:
        sigma1 (float): max principal stress
        sigma3 (float): min principal stress
        ax (matplotlib axes): axes to draw circle. 
        color (str, optional): color of the line. Defaults to 'k' (black).

    Returns:
        figure and axes (optional): if not given ax, returns a figure and an axes
    """
    c, r = _get_center_and_radius(sigma1, sigma3)
    circle1 = plt.Circle((c, 0), r, fill=False, color=color, **kwargs)

    ymin, ymax = 0, sigma1-sigma3
    xmin, xmax = 0, sigma1+sigma3

    ax.add_patch(circle1)
    ax.set_aspect('equal')
    ax.axhline(0, c='k')
    ax.axvline(0, c='k')
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)


def mohr3d(sigma1:float, sigma2:float, sigma3:float, ax, color='k', **kwargs):
    """create mohr diagram for 3d stress state

    Args:
        sigma1 (float): max principal stress
        sigma2 (float): intermediate principal stress
        sigma3 (float): min principal stress
        ax (matplotlib axes): axes to draw circle. 
        color (str, optional): color of the line. Defaults to 'k' (black).

    Returns:
        figure and axes (optional): if not given ax, returns a figure and an axes
    """
    
    c1, r1 = _get_center_and_radius(sigma1, sigma3)
    c2, r2 = _get_center_and_radius(sigma1, sigma2)
    c3, r3 = _get_center_and_radius(sigma2, sigma3)

    ymin, ymax = 0, sigma1 - sigma3
    xmin, xmax = 0, sigma1 + sigma3

    # largest circle (s1 and s3)
    circle1 = plt.Circle((c1, 0), r1, fill=False, color=color, **kwargs)
    # (s1 and s2)
    circle2 = plt.Circle((c2, 0), r2, fill=False, color=color, **kwargs)
    # (s2 and s3)
    circle3 = plt.Circle((c3, 0), r3, fill=False, color=color, **kwargs)

    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.add_patch(circle3)
    ax.axhline(0, c='k')
    ax.axvline(0, c='k')
    ax.set_aspect('equal')
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)


## zobackogram
## formulas for frictional faulting theory 
from math import sqrt

def get_q(mu:float) -> float:
    """calculate ratio of s1/s3 as a function of mu (friction coefficient)

    Args:
        mu (float): friction coefficient

    Returns:
        float: q
    """
    return (sqrt(mu**2 + 1) + mu)**2

def get_s1(s3:float, pp:float, q:float) -> float:
    return q * (s3-pp) + pp 

def get_s3(s1:float, pp:float, q:float) -> float:
    return pp + ((s1-pp)/q)


def zobackogram(sv:float, pp:float, ax=None, mu=0.6, color='k', **kwargs):
    """construct stress polygon (zoback-o-gram) for a constant depth


    Args:
        sv (float): vertical/overburden stress
        pp (float): pore pressure
        mu (float, optional): friction coefficient. Defaults to 0.6.
    """
    q = get_q(mu)
    s1 = get_s1(sv, pp, q)
    s3 = get_s3(sv, pp, q)

    # params = {'color':'k'}
    if ax is None:
        f, ax = plt.subplots()

    # sh=sH line
    ax.plot((s3, s1), (s3, s1), color=color, **kwargs)
    # normal fault constraint
    ax.plot((s3, s3), (s3, sv), color=color, **kwargs)
    # normal - strike-slip limit
    ax.plot((s3, sv), (sv, sv), color=color, **kwargs)
    # strike-slip fault constraint
    ax.plot((s3, sv), (sv, s1), color=color, **kwargs)
    # strike-slip -  reverse limit
    ax.plot((sv, sv), (sv, s1), color=color, **kwargs)
    # reverse fault limit
    ax.plot((sv, s1), (s1, s1), color=color, **kwargs)

if __name__ == "__main__":
    from failure_criterion import MohrCoulomb

    # normal, shear = MohrCoulomb(.9, 0).mohr(0, 25)

    # # zobackogram(10, 3, mu=0.6, c='k', label='0.6')
    # # zobackogram(10, 3, mu=0.7, c='r')
    # f, ax = plt.subplots(1, 2, figsize=(10, 5))
    # mohr3d(25, 10, 5, ax=ax[0])
    # ax[0].set_ylim(0, 15)
    # ax[0].plot(normal, shear)
    import matplotlib.pyplot as plt

    Sv = 11000
    Pp = 4400
    zobackogram(Sv, Pp)
    plt.show()