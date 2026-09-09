import argparse, json
from .channels import available_channels, get_channels
from .zernike import zernike_descriptor
from .synthetic import patch
def main(argv=None):
    p=argparse.ArgumentParser(prog="zernike-ppi"); sub=p.add_subparsers(dest="cmd",required=True)
    sub.add_parser("channels"); d=sub.add_parser("demo"); d.add_argument("--channels",default="shape")
    args=p.parse_args(argv)
    if args.cmd=="channels":
        for n in available_channels():
            c=get_channels([n])[n]; print(f"{n}\n    {c.description}\n    relation: {c.relation.replace('↔','<->')}")
        return
    names=[x.strip() for x in args.channels.split(",")]; get_channels(names); a=patch(); b=patch("concave","negative")
    from .descriptors import PatchDescriptor
    pa=PatchDescriptor(0,[0,0,0],[0,0,1],{n:zernike_descriptor(a[n]) for n in names}); pb=PatchDescriptor(0,[0,0,0],[0,0,1],{n:zernike_descriptor(b[n]) for n in names})
    from .matching import match_descriptors
    print(json.dumps(match_descriptors([pa],[pb],names)[0],indent=2))

if __name__ == "__main__":
    main()
