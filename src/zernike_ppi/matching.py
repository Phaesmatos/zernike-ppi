import numpy as np
from scipy.spatial.distance import cdist
def match_descriptors(a,b,channels,weights=None,top_k=None):
    weights=weights or {}; rows=[]
    for name in channels:
        target="acceptor" if name=="donor" else "donor" if name=="acceptor" else name
        x=np.stack([p.channels[name] for p in a]); y=np.stack([p.channels[target] for p in b])
        if name in ("shape","charge"): y=-y
        d=cdist(x,y,"euclidean")/np.sqrt(x.shape[1]); rows.append((name,d))
    total=sum(weights.get(n,1.0)*d for n,d in rows); ia,ib=np.unravel_index(np.argsort(total,axis=None),total.shape)
    result=[]
    for rank,(i,j) in enumerate(zip(ia,ib),1):
        r={"rank":rank,"total_score":float(total[i,j]),"patch_a":a[i].patch_id,"patch_b":b[j].patch_id}; r.update({f"score_{n}":float(d[i,j]) for n,d in rows}); result.append(r)
    return result if top_k is None else result[:top_k]
