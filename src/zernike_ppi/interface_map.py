"""Bound-complex geometric interface complementarity maps from surface points."""
from __future__ import annotations
import csv, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree, QhullError
from scipy.interpolate import griddata
from .projection import local_frame
from .zernike import zernike_descriptor

def arrays(surface):
    xyz=np.asarray(surface['xyz'],float)
    if 'normal' in surface: normal=np.asarray(surface['normal'],float)
    else: normal=np.column_stack([surface['nx'],surface['ny'],surface['nz']]).astype(float)
    normal/=np.linalg.norm(normal,axis=1,keepdims=True)+1e-12
    return xyz,normal

def interface_masks(a,b,distance=3.0):
    xa,_=arrays(a); xb,_=arrays(b)
    da=cKDTree(xb).query(xa)[0]; db=cKDTree(xa).query(xb)[0]
    return da < distance, db < distance

def mean_normal(normals, central, eps=1e-8):
    v=np.sum(normals,axis=0); mag=np.linalg.norm(v)
    if mag>eps: return v/mag,'mean'
    c=np.asarray(central,float); return c/(np.linalg.norm(c)+eps),'central-fallback'

def height_field(xyz,normals,indices,center,normal,radius,grid_size):
    e1,e2,n=local_frame(normal); d=xyz[indices]-center; u=d@e1; v=d@e2; w=d@n
    field=np.zeros((grid_size,grid_size)); count=np.zeros_like(field)
    gx=np.clip(((u/radius+1)*(grid_size-1)/2).astype(int),0,grid_size-1)
    gy=np.clip(((v/radius+1)*(grid_size-1)/2).astype(int),0,grid_size-1)
    np.add.at(field,(gy,gx),w); np.add.at(count,(gy,gx),1); nz=count>0; field[nz]/=count[nz]
    return field

def deterministic_basis(n):
    n=np.asarray(n,float); n/=np.linalg.norm(n)
    axes=np.eye(3); ref=axes[np.argmin(np.abs(axes@n))]
    e1=ref-(ref@n)*n; e1/=np.linalg.norm(e1); return e1,np.cross(n,e1)

def run_interface_map(surface_a, surface_b, interface_distance=3., patch_radius=6., sample_every=10, order=20, grid_size=64, grid_spacing=.5, axis_max_distance=None):
    xa,na=arrays(surface_a); xb,nb=arrays(surface_b); ma,mb=interface_masks(surface_a,surface_b,interface_distance)
    if not ma.any() or not mb.any(): raise ValueError('No interface points: increase --interface-distance or supply a bound complex.')
    ia=np.flatnonzero(ma); ib=np.flatnonzero(mb); ta,tb=cKDTree(xa),cKDTree(xb)
    pa, pb = xa[ma].mean(0), xb[mb].mean(0); P=np.vstack([xa[ma],xb[mb]]).mean(0)
    NA,_=mean_normal(na[ma],na[ia[0]]); NB,_=mean_normal(nb[mb],nb[ib[0]])
    axis=NA-NB if NA@NB<0 else NA
    axis/=np.linalg.norm(axis)+1e-12
    if axis@(pb-pa)<0: axis=-axis
    e1,e2=deterministic_basis(axis); selected=ia[::max(1,int(sample_every))]; rows=[]
    candidates=ib
    axis_max_distance=axis_max_distance or (patch_radius+interface_distance)
    for pid,i in enumerate(selected):
        aidx=np.asarray(ta.query_ball_point(xa[i],patch_radius),int)
        nA,methodA=mean_normal(na[aidx],na[i]); delta=xb[candidates]-xa[i]; t=delta@nA; perp=np.linalg.norm(delta-t[:,None]*nA,axis=1)
        valid=(t>0)&(t<=axis_max_distance)
        if valid.any(): j=candidates[np.where(valid)[0][np.argmin(perp[valid])]]
        else: j=candidates[np.argmin(np.linalg.norm(delta,axis=1))]
        bidx=np.asarray(tb.query_ball_point(xb[j],patch_radius),int); nB,methodB=mean_normal(nb[bidx],nb[j])
        fa=height_field(xa,na,aidx,xa[i],nA,patch_radius,grid_size); fb=height_field(xb,nb,bidx,xb[j],nB,patch_radius,grid_size)
        za=zernike_descriptor(fa,order,grid_size); zb=zernike_descriptor(fb,order,grid_size); dist=float(np.linalg.norm(za+zb)/np.sqrt(len(za))); comp=1/(1+dist)
        mid=(xa[i]+xb[j])/2; q=mid-P
        rows.append(dict(patch_id=pid,patch_a=int(i),patch_b=int(j),center_a=xa[i],center_b=xb[j],mean_normal_a=nA,mean_normal_b=nB,axis_distance=float(np.linalg.norm((xb[j]-xa[i])-(xb[j]-xa[i])@nA*nA)),axial_parameter_t=float((xb[j]-xa[i])@nA),normal_opposition=float(-(nA@nB)),suspicious_normals=bool(nA@nB>0),normal_method_a=methodA,normal_method_b=methodB,zernike_distance=dist,complementarity_score=comp,x=float(q@e1),y=float(q@e2)))
    pts=np.array([[r['x'],r['y']] for r in rows]); vals=np.array([r['complementarity_score'] for r in rows]); lo=pts.min(0)-grid_spacing; hi=pts.max(0)+grid_spacing
    gx=np.arange(lo[0],hi[0]+grid_spacing,grid_spacing); gy=np.arange(lo[1],hi[1]+grid_spacing,grid_spacing); XX,YY=np.meshgrid(gx,gy)
    try:
        if len(rows)>=3:
            # Keep cells outside the convex hull as NaN: values must be
            # supported by measured patch-pair samples, not extrapolated.
            mat=griddata(pts,vals,(XX,YY),method='linear')
        else: mat=griddata(pts,vals,(XX,YY),method='nearest')
    except QhullError:
        # Sparse/collinear patch samples cannot support 2D triangulation.
        # Nearest-neighbour is an explicit conservative fallback.
        mat=griddata(pts,vals,(XX,YY),method='nearest')
    return rows,dict(interface_center=P,interface_center_a=pa,interface_center_b=pb,axis=axis,e1=e1,e2=e2,interface_count_a=int(ma.sum()),interface_count_b=int(mb.sum()),normal_dot=float(NA@NB),grid_x=gx,grid_y=gy,matrix=mat,mask_a=ma,mask_b=mb)

def write_interface_map(output, rows, meta, config, surface_a, surface_b):
    import yaml, matplotlib.pyplot as plt
    out=Path(output); out.mkdir(parents=True,exist_ok=True)
    np.savez(out/'interface_points_A.npz',mask=meta['mask_a']); np.savez(out/'interface_points_B.npz',mask=meta['mask_b'])
    simple=[]
    for r in rows:
        d={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in r.items()}; simple.append(d)
    keys=['patch_id','patch_a','patch_b','x','y','complementarity_score','zernike_distance','normal_opposition','axis_distance','axial_parameter_t','suspicious_normals','normal_method_a','normal_method_b']
    with (out/'patch_pairs.csv').open('w',newline='') as f: w=csv.DictWriter(f,keys); w.writeheader(); w.writerows([{k:r[k] for k in keys} for r in simple])
    with (out/'projected_complementarity.csv').open('w',newline='') as f: w=csv.DictWriter(f,['patch_id','x','y','complementarity_score','zernike_distance','normal_opposition','axis_distance']); w.writeheader(); w.writerows([{k:r[k] for k in w.fieldnames} for r in simple])
    np.save(out/'complementarity_matrix.npy',meta['matrix']); np.savetxt(out/'complementarity_matrix.csv',meta['matrix'],delimiter=',')
    geom={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in meta.items() if k not in ('matrix','mask_a','mask_b')}; (out/'interface_geometry.json').write_text(json.dumps(geom,indent=2)); (out/'config.yaml').write_text(yaml.safe_dump(config,sort_keys=True))
    fig,ax=plt.subplots(figsize=(6,5)); im=ax.imshow(meta['matrix'],origin='lower',extent=[meta['grid_x'][0],meta['grid_x'][-1],meta['grid_y'][0],meta['grid_y'][-1]],aspect='auto'); ax.scatter([r['x'] for r in rows],[r['y'] for r in rows],c='k',s=8); ax.set(xlabel='plane e1 (Å)',ylabel='plane e2 (Å)',title='Geometric complementarity'); fig.colorbar(im,ax=ax,label='complementarity'); fig.tight_layout(); fig.savefig(out/'complementarity_map.png',dpi=160); plt.close(fig)
    xa,_=arrays(surface_a); xb,_=arrays(surface_b); fig=plt.figure(figsize=(7,6)); ax=fig.add_subplot(111,projection='3d'); ax.scatter(*xa.T,s=1,alpha=.18); ax.scatter(*xb.T,s=1,alpha=.18); ax.scatter(*xa[meta['mask_a']].T,s=4); ax.scatter(*xb[meta['mask_b']].T,s=4)
    for r in rows:
        pair=np.vstack([r['center_a'],r['center_b']]); ax.plot(*pair.T,color='k',alpha=.35,lw=.6)
    P=meta['interface_center']; N=meta['axis']; ax.quiver(*P,*N,length=8,color='k'); ax.set_title('Interface points, facing pairs, and A→B axis'); fig.tight_layout(); fig.savefig(out/'interface_geometry.png',dpi=160); plt.close(fig)
