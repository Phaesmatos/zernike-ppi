import numpy as np
from zernike_ppi.interface_map import interface_masks, run_interface_map, deterministic_basis

def opposing(n=9, gap=2.0):
 x,y=np.meshgrid(np.linspace(-2,2,n),np.linspace(-2,2,n)); xy=np.c_[x.ravel(),y.ravel()]
 # A outward normal points up into the gap; B outward normal points down.
 a={'xyz':np.c_[xy,np.zeros(len(xy))],'normal':np.tile([0.,0.,1.],(len(xy),1))}
 b={'xyz':np.c_[xy,np.full(len(xy),gap)],'normal':np.tile([0.,0.,-1.],(len(xy),1))}
 return a,b
def test_interface_and_axis():
 a,b=opposing(); ma,mb=interface_masks(a,b,3.0); assert ma.all() and mb.all()
 rows,m=run_interface_map(a,b,sample_every=10,order=4,grid_size=24,grid_spacing=.5)
 assert len(rows)>0 and all(r['axial_parameter_t']>0 for r in rows)
 assert all(r['normal_opposition']>.99 for r in rows)
 assert m['axis'][2]>0 and abs(m['axis']@m['e1'])<1e-10 and np.isfinite(m['matrix']).any()
def test_tilted_interface_is_detected():
 a,b=opposing(); rot=np.array([[1,0,0],[0,0,-1],[0,1,0.]])
 for s in (a,b): s['xyz']=s['xyz']@rot.T; s['normal']=s['normal']@rot.T
 rows,m=run_interface_map(a,b,sample_every=20,order=4,grid_size=20)
 assert len(rows)>0 and m['axis'][1]<0 and all(r['normal_opposition']>.99 for r in rows)
