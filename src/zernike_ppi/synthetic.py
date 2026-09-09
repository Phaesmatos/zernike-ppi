import numpy as np
def patch(shape="dome",charge="positive",hydro="center",size=64):
    y,x=np.mgrid[-1:1:complex(size),-1:1:complex(size)]; r2=x*x+y*y
    z=np.exp(-5*r2)*(1 if shape in ("dome","convex") else -1); q=np.exp(-8*r2)*(1 if charge=="positive" else -1)
    h=np.exp(-12*r2) if hydro=="center" else np.exp(-40*(np.sqrt(r2)-.6)**2)
    return {"shape":z,"charge":q,"hydrophobicity":h}
