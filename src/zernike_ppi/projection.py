import numpy as np
def local_frame(normal):
 z=np.asarray(normal,float); z/=np.linalg.norm(z)+1e-12; ref=np.array([0.,0.,1.]) if abs(z[2])<.9 else np.array([1.,0.,0.]); x=np.cross(ref,z); x/=np.linalg.norm(x)+1e-12; return x,np.cross(z,x),z
def project_patch(surface, index, radius=6., grid_size=64, sample_every=1):
 c=surface["xyz"][index]; n=surface.get("normal",np.tile([0,0,1.],(len(surface["xyz"]),1)))[index]; x,y,z=local_frame(n); d=surface["xyz"]-c; u=d@x; v=d@y; w=d@z; keep=(u*u+v*v<=radius*radius)
 yy,xx=np.mgrid[-1:1:complex(grid_size),-1:1:complex(grid_size)]; fields={}
 for name,vals in surface.get("fields",{}).items():
  arr=np.zeros_like(xx); sel=np.where(keep)[0][::sample_every]; gx=((u[sel]/radius+1)*(grid_size-1)/2).astype(int).clip(0,grid_size-1); gy=((v[sel]/radius+1)*(grid_size-1)/2).astype(int).clip(0,grid_size-1); arr[gy,gx]=vals[sel]; fields[name]=arr
 fields.setdefault("shape",np.zeros_like(xx)); return {"patch_id":index,"center":c,"normal":n,"fields":fields}
