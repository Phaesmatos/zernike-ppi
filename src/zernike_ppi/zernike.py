"""Small, dependency-light 2D Zernike backend."""
import math
import numpy as np

def radial(n, m, r):
    out = np.zeros_like(r, dtype=float)
    for s in range((n-m)//2 + 1):
        out += ((-1)**s * math.factorial(n-s) /
                (math.factorial(s)*math.factorial((n+m)//2-s)*math.factorial((n-m)//2-s))) * r**(n-2*s)
    return out

class ZernikeBackend:
    def __init__(self, order=20, grid_size=64):
        self.order, self.grid_size = int(order), int(grid_size)
        # One invariant per (n, |m|) pair; this is the 121-value order-20
        # representation used by the Zernike2D/zepyros workflow.
        self.modes = [(n,m) for n in range(self.order+1) for m in range(0,n+1,2)]
    @property
    def dimension(self): return len(self.modes)
    def compute(self, field):
        a = np.asarray(field, float); g = a.shape[0]
        y,x = np.mgrid[-1:1:complex(g), -1:1:complex(g)]
        r=np.hypot(x,y); t=np.arctan2(y,x); mask=r<=1
        vals=np.nan_to_num(a[mask]); rr=r[mask]; tt=t[mask]
        raw_mean=vals.mean(); vals=(vals-raw_mean)/(vals.std()+1e-12); result=[]
        for n,m in self.modes:
            basis=radial(n,abs(m),rr)*np.exp(1j*m*tt)
            # Keep the signed DC term: it makes sign-complement channels
            # (shape/charge) distinguishable while all higher modes remain
            # rotation-invariant magnitudes, matching the verso convention.
            coeff=np.mean(vals*np.conj(basis))
            result.append(float(raw_mean) if (n,m)==(0,0) else abs(coeff))
        return np.asarray(result, dtype=float)

def zernike_descriptor(field, order=20, grid_size=64):
    return ZernikeBackend(order, grid_size).compute(field)
