from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Channel:
    name: str; relation: str; description: str; complement: str = "similarity"
    def values(self, surface): raise NotImplementedError

class PropertyChannel(Channel):
    def __init__(self, name, relation, description, mapping, complement="similarity"):
        super().__init__(name, relation, description, complement); object.__setattr__(self, "mapping", mapping)
    def values(self, surface):
        names=np.asarray(surface.get("residue_name", ["UNK"]*len(surface["xyz"])))
        return np.array([self.mapping.get(str(x).upper(), 0.0) for x in names], float)

HYDRO={"ILE":4.5,"VAL":4.2,"LEU":3.8,"PHE":2.8,"CYS":2.5,"MET":1.9,"ALA":1.8,"GLY":-0.4,"THR":-0.7,"SER":-0.8,"TRP":-0.9,"TYR":-1.3,"PRO":-1.6,"HIS":-3.2,"GLU":-3.5,"GLN":-3.5,"ASP":-3.5,"ASN":-3.5,"LYS":-3.9,"ARG":-4.5}
CHARGE={"ARG":1.0,"LYS":1.0,"HIS":0.25,"ASP":-1.0,"GLU":-1.0}
DONOR={"ARG":1.,"LYS":1.,"HIS":.7,"ASN":.6,"GLN":.6,"SER":.5,"THR":.5,"TYR":.4,"TRP":.3}
ACCEPTOR={"ASP":1.,"GLU":1.,"ASN":.7,"GLN":.7,"SER":.5,"THR":.5,"TYR":.4,"HIS":.4}
CHANNELS={
 "hydrophobicity":PropertyChannel("hydrophobicity","similar ↔ similar","Kyte-Doolittle residue scale",HYDRO),
 "charge":PropertyChannel("charge","positive ↔ negative","residue charge proxy; HIS=+0.25",CHARGE,"negate"),
 "donor":PropertyChannel("donor","donor ↔ acceptor","annotated donor propensity",DONOR,"cross:acceptor"),
 "acceptor":PropertyChannel("acceptor","acceptor ↔ donor","annotated acceptor propensity",ACCEPTOR,"cross:donor")}
class ShapeChannel(Channel):
    def __init__(self): super().__init__("shape","convex ↔ concave","local surface height; orientation transformed","negate")
    def values(self,surface): return np.asarray(surface.get("height",surface.get("z",np.zeros(len(surface["xyz"])))),float)
CHANNELS["shape"]=ShapeChannel()
def available_channels(): return list(CHANNELS)
def get_channels(names):
    if names==["all"]: names=available_channels()
    bad=set(names)-set(CHANNELS)
    if bad: raise ValueError(f"Unknown channel(s): {', '.join(sorted(bad))}")
    return {n:CHANNELS[n] for n in names}
