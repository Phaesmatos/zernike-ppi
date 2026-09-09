import numpy as np
from zernike_ppi.zernike import ZernikeBackend
from zernike_ppi.synthetic import patch
from zernike_ppi.descriptors import PatchDescriptor
from zernike_ppi.matching import match_descriptors

def pd(fields): return PatchDescriptor(0,np.zeros(3),np.array([0,0,1]),{k:ZernikeBackend(order=4,grid_size=32).compute(v) for k,v in fields.items()})
def test_order20_dimension(): assert ZernikeBackend().dimension==121
def test_rotation_invariant():
 a=patch()["shape"]; b=np.rot90(a); assert np.allclose(ZernikeBackend(8,64).compute(a),ZernikeBackend(8,64).compute(b),atol=.08)
def test_shape_complement():
 a=pd(patch("convex")); b=pd(patch("concave")); wrong=pd(patch("convex"));
 assert match_descriptors([a],[b],["shape"])[0]["total_score"] < match_descriptors([a],[wrong],["shape"])[0]["total_score"]
def test_channel_weight_and_cross_mapping():
 a=pd({"donor":np.ones((32,32))}); b=pd({"acceptor":np.ones((32,32))})
 assert match_descriptors([a],[b],["donor"])[0]["score_donor"] < 1e-8
def test_selections_change_score():
 a=pd(patch()); b=pd(patch("concave","negative")); r1=match_descriptors([a],[b],["shape"])[0]; r2=match_descriptors([a],[b],["shape","charge"])[0]; assert r1["total_score"] != r2["total_score"]
