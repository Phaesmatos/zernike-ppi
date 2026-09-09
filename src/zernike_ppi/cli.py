import argparse, json, csv, time
from pathlib import Path
from .channels import available_channels, get_channels
from .zernike import zernike_descriptor
from .synthetic import patch
def main(argv=None):
    p=argparse.ArgumentParser(prog="zernike-ppi"); sub=p.add_subparsers(dest="cmd",required=True)
    sub.add_parser("channels"); d=sub.add_parser("demo"); d.add_argument("--channels",default="shape")
    def common(q, pair=False):
        q.add_argument('--protein-a' if pair else '--protein',required=True); 
        if pair: q.add_argument('--protein-b',required=True)
        q.add_argument('--surface-a'); q.add_argument('--surface-b') if pair else None
        q.add_argument('--channels',default='shape'); q.add_argument('--patch-radius',type=float,default=6.); q.add_argument('--sample-every',type=int,default=10); q.add_argument('--zernike-order',type=int,default=20); q.add_argument('--grid-size',type=int,default=64); q.add_argument('--output',required=True)
    s=sub.add_parser('surface'); common(s); s.add_argument('--probe-radius',type=float,default=1.4)
    dsc=sub.add_parser('describe'); common(dsc)
    m=sub.add_parser('match'); common(m,True); m.add_argument('--top-k',type=int,default=20); m.add_argument('--channel-weight',action='append',default=[]); m.add_argument('--descriptors-a'); m.add_argument('--descriptors-b')
    v=sub.add_parser('visualize'); v.add_argument('--protein-a',required=True); v.add_argument('--protein-b',required=True); v.add_argument('--matches',required=True); v.add_argument('--rank',type=int,default=1); v.add_argument('--output',default='match.png')
    b=sub.add_parser('benchmark'); b.add_argument('--dataset'); b.add_argument('--config'); b.add_argument('--output',required=True)
    args=p.parse_args(argv)
    if args.cmd=="channels":
        for n in available_channels():
            c=get_channels([n])[n]; print(f"{n}\n    {c.description}\n    relation: {c.relation.replace('↔','<->')}")
        return
    if args.cmd=='surface':
        from .pipeline import surface_from_input, save_surface
        t=time.time(); s=surface_from_input(args.protein); save_surface(args.output,s); print(f"wrote {len(s['xyz'])} surface points in {time.time()-t:.3f}s"); return
    if args.cmd=='describe':
        from .pipeline import build_descriptors, save_descriptor_bundle
        names=[x.strip() for x in args.channels.split(',')]; t=time.time(); ds,_=build_descriptors(args.protein,names=names,radius=args.patch_radius,sample_every=args.sample_every,order=args.zernike_order,grid_size=args.grid_size); save_descriptor_bundle(args.output,ds,vars(args)); print(f"wrote {len(ds)} patches in {time.time()-t:.3f}s"); return
    if args.cmd=='match':
        from .pipeline import build_descriptors, load_descriptor_bundle
        from .matching import match_descriptors
        names=available_channels() if args.channels=='all' else [x.strip() for x in args.channels.split(',')]; get_channels(names); weights={x.split('=',1)[0]:float(x.split('=',1)[1]) for x in args.channel_weight}
        out=Path(args.output); out.mkdir(parents=True,exist_ok=True); t=time.time();
        if args.descriptors_a and args.descriptors_b:
            from .pipeline import load_descriptor_bundle
            da=load_descriptor_bundle(args.descriptors_a,names); db=load_descriptor_bundle(args.descriptors_b,names)
        else:
            da,_=build_descriptors(args.protein_a,names=names,radius=args.patch_radius,sample_every=args.sample_every,order=args.zernike_order,grid_size=args.grid_size); db,_=build_descriptors(args.protein_b,names=names,radius=args.patch_radius,sample_every=args.sample_every,order=args.zernike_order,grid_size=args.grid_size)
        rows=match_descriptors(da,db,names,weights,args.top_k)
        amap={p.patch_id:p for p in da}; bmap={p.patch_id:p for p in db}
        for r in rows:
            pa=amap[r['patch_a']]; pb=bmap[r['patch_b']]; r.update({'center_a_x':float(pa.center[0]),'center_a_y':float(pa.center[1]),'center_a_z':float(pa.center[2]),'center_b_x':float(pb.center[0]),'center_b_y':float(pb.center[1]),'center_b_z':float(pb.center[2]),'nearest_residue_a':pa.metadata.get('residue_name',''),'nearest_residue_b':pb.metadata.get('residue_name','')})
        keys=list(rows[0]) if rows else ['rank']; (out/'matches.json').write_text(json.dumps(rows,indent=2));
        with (out/'matches.csv').open('w',newline='') as f: w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)
        (out/'config.yaml').write_text(json.dumps(vars(args),indent=2)); print(f"wrote {len(rows)} matches in {time.time()-t:.3f}s"); return
    if args.cmd=='visualize':
        from .visualization import plot_match; plot_match(args.protein_a,args.protein_b,args.matches,args.rank,args.output); print(f'wrote {args.output}'); return
    if args.cmd=='benchmark':
        from .benchmark import write_summary
        if not args.dataset: raise SystemExit('--dataset is required for benchmark ranking diagnostics')
        rows=[]
        for r in csv.DictReader(open(args.dataset,newline='')): rows.append({'complex_id':r.get('complex_id',''),'status':'ranking-only','channels':'configured','ROC_AUC':'nan','PR_AUC':'nan'})
        write_summary(rows,args.output); (Path(args.output)/'config.yaml').write_text(Path(args.config).read_text() if args.config else 'ranking-only: true\n'); print(f'wrote benchmark results for {len(rows)} complexes'); return
    names=[x.strip() for x in args.channels.split(",")]; get_channels(names); a=patch(); b=patch("concave","negative")
    from .descriptors import PatchDescriptor
    pa=PatchDescriptor(0,[0,0,0],[0,0,1],{n:zernike_descriptor(a[n]) for n in names}); pb=PatchDescriptor(0,[0,0,0],[0,0,1],{n:zernike_descriptor(b[n]) for n in names})
    from .matching import match_descriptors
    print(json.dumps(match_descriptors([pa],[pb],names)[0],indent=2))

if __name__ == "__main__":
    main()
