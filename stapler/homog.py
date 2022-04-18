"""
Directly ported from: https://github.com/willsheffler/homog
Under the following license: http://www.apache.org/licenses/LICENSE-2.0
"""
import numpy as np


def hnorm(a):
    a = np.asanyarray(a)
    return np.sqrt(np.sum(a[..., :3] * a[..., :3], axis=-1))


def hnormalized(a):
    a = np.asanyarray(a)
    if (not a.shape and len(a) == 3) or (a.shape and a.shape[-1] == 3):
        a, tmp = np.zeros(a.shape[:-1] + (4, )), a
        a[..., :3] = tmp
    return a / hnorm(a)[..., None]


def hcross(a, b):
    a = np.asanyarray(a)
    b = np.asanyarray(b)
    c = np.zeros(np.broadcast(a, b).shape, dtype=a.dtype)
    c[..., :3] = np.cross(a[..., :3], b[..., :3])
    return c


def hpoint(point):
    point = np.asanyarray(point)
    if point.shape[-1] == 4: return point
    elif point.shape[-1] == 3:
        r = np.ones(point.shape[:-1] + (4, ))
        r[..., :3] = point
        return r
    else:
        raise ValueError('point must len 3 or 4')


def hstub(u, v, w, cen=None):
    u, v, w = hpoint(u), hpoint(v), hpoint(w)
    assert u.shape == v.shape == w.shape
    if not cen: cen = u
    cen = hpoint(cen)
    assert cen.shape == u.shape
    stubs = np.empty(u.shape[:-1] + (4, 4))
    stubs[..., :, 0] = hnormalized(u - v)
    stubs[..., :, 2] = hnormalized(hcross(stubs[..., :, 0], w - v))
    stubs[..., :, 1] = hcross(stubs[..., :, 2], stubs[..., :, 0])
    stubs[..., :, 3] = hpoint(cen[..., :])
    return stubs
