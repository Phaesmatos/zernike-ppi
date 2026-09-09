import numpy as np
from scipy.spatial import cKDTree
def approximate_surface(atoms, probe_radius=1.4, points_per_atom=24):
 radii={"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80}; xyz=atoms["xyz"]; elems=np.array(["C"]*len(xyz))
 # PDB reader does not retain element in the minimal format; carbon is a safe fallback.
 dirs=[]; phi=(1+5**.5)/2
 for k in range(points_per_atom):
  z=1-2*(k+.5)/points_per_atom; t=2*np.pi*k/phi; dirs.append([np.sqrt(1-z*z)*np.cos(t),np.sqrt(1-z*z)*np.sin(t),z])
 dirs=np.asarray(dirs); pts=[]; normals=[]; rn=[]; rid=[]; ch=[]
 allxyz=xyz
 for i,c in enumerate(xyz):
  r=radii.get(elems[i],1.7)+probe_radius; cand=c+r*dirs
  tree=cKDTree(allxyz); buried=tree.query(cand,k=1)[0] < r-0.25
  pts.extend(cand[~buried]); normals.extend(dirs[~buried]); rn.extend([atoms.get("residue_name",["UNK"])[i]]*int(sum(~buried))); rid.extend([atoms.get("residue_id",[""])[i]]*int(sum(~buried))); ch.extend([atoms.get("chain",[""])[i]]*int(sum(~buried)))
 return {"xyz":np.asarray(pts),"normal":np.asarray(normals),"residue_name":np.asarray(rn),"residue_id":np.asarray(rid),"chain":np.asarray(ch)}
