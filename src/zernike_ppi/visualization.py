import csv
import numpy as np
import matplotlib.pyplot as plt
from .pipeline import load_surface, build_descriptors
def plot_match(protein_a, protein_b, matches_csv, rank, output):
    row=next(r for r in csv.DictReader(open(matches_csv,newline='')) if int(r['rank'])==rank)
    sa=load_surface(protein_a); sb=load_surface(protein_b); ia=int(row['patch_a']); ib=int(row['patch_b'])
    fig,ax=plt.subplots(1,2,figsize=(10,4));
    for a,s,i,title in ((ax[0],sa,ia,'Protein A'),(ax[1],sb,ib,'Protein B')):
        a.scatter(s['xyz'][:,0],s['xyz'][:,1],s=4,alpha=.35); a.scatter(*s['xyz'][i,:2],c='red',s=50); a.set_title(title); a.set_aspect('equal')
    fig.tight_layout(); fig.savefig(output,dpi=150); plt.close(fig)
