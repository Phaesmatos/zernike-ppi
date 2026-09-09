"""Minimal PDB reader and DMS-compatible surface CSV reader."""
import csv, numpy as np
def read_pdb(path):
 xyz=[]; names=[]; residues=[]; chains=[]
 for line in open(path, encoding="utf-8", errors="ignore"):
  if line.startswith(("ATOM  ","HETATM")) and line[12:16].strip() not in ("H",):
   try: xyz.append([float(line[30:38]),float(line[38:46]),float(line[46:54])]); names.append(line[17:20].strip()); residues.append(line[22:26].strip()); chains.append(line[21].strip())
   except ValueError: pass
 return {"xyz":np.asarray(xyz,float),"residue_name":np.asarray(names),"residue_id":np.asarray(residues),"chain":np.asarray(chains)}
def read_surface_csv(path):
 rows=list(csv.DictReader(open(path,newline=""))); xyz=np.array([[float(r["x"]),float(r["y"]),float(r["z"])] for r in rows])
 out={"xyz":xyz}
 for k in ("nx","ny","nz","residue_name","residue_id","chain"): out[k]=np.array([r.get(k,0) for r in rows])
 return out
