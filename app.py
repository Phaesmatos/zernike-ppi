"""Local Streamlit interface for the existing ZernikePPI interface-map API."""
import tempfile, traceback
from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np
from zernike_ppi.ui import InterfaceMapConfig, read_pdb_upload, chains, analyse_atoms, result_files

st.set_page_config(page_title='ZernikePPI',layout='wide')
st.title('ZernikePPI')
st.caption('Protein–Protein Interface Complementarity using 2D Zernike Descriptors')

def save_upload(upload, folder, label):
    if upload is None: return None
    p=Path(folder)/f'{label}.pdb'; p.write_bytes(upload.getvalue()); return p

mode=st.radio('Input mode',['One PDB complex','Two PDB files'],horizontal=True)
with tempfile.TemporaryDirectory() as input_dir:
    if mode=='One PDB complex':
        uploaded=st.file_uploader('Upload complex PDB',type='pdb',key='complex')
        pa=save_upload(uploaded,input_dir,'complex') if uploaded else None
        atoms_a=atoms_b=read_pdb_upload(pa) if pa else None
        if atoms_a is not None:
            available=chains(atoms_a)
            if len(available)<2: st.error('Choose a complex PDB containing at least two non-empty chains.'); ca=cb=None
            else:
                c1,c2=st.columns(2); ca=c1.selectbox('Protein A chain',available,index=0); cb=c2.selectbox('Protein B chain',[c for c in available if c!=ca],index=0)
        else: ca=cb=None
    else:
        c1,c2=st.columns(2); ua=c1.file_uploader('Protein A — Upload PDB',type='pdb',key='a'); ub=c2.file_uploader('Protein B — Upload PDB',type='pdb',key='b')
        pa=save_upload(ua,input_dir,'protein_a') if ua else None; pb=save_upload(ub,input_dir,'protein_b') if ub else None
        atoms_a=read_pdb_upload(pa) if pa else None; atoms_b=read_pdb_upload(pb) if pb else None
        ca=cb=None
        if atoms_a is not None: ca=c1.selectbox('Protein A chain',chains(atoms_a))
        if atoms_b is not None: cb=c2.selectbox('Protein B chain',chains(atoms_b))
    if atoms_a is not None:
        st.info(f"A: {len(atoms_a['xyz'])} atoms; chains: {', '.join(chains(atoms_a))}")
    if atoms_b is not None and mode!='One PDB complex': st.info(f"B: {len(atoms_b['xyz'])} atoms; chains: {', '.join(chains(atoms_b))}")

with st.expander('Advanced settings'):
    source=st.selectbox('Surface source',['Automatic from PDB','DMS, if available'],index=0)
    if source!='Automatic from PDB': st.warning('The local UI currently uses the validated automatic PDB surface generator. DMS surfaces remain available through the CLI interface-map command.')
    c1,c2,c3=st.columns(3); interface_distance=c1.number_input('Interface distance (Å)',value=3.0,min_value=.1); patch_radius=c2.number_input('Patch radius (Å)',value=6.0,min_value=.5); sample_every=c3.number_input('Sample every N points',value=10,min_value=1,step=1)
    c1,c2=st.columns(2); order=c1.number_input('Zernike order',value=20,min_value=2,step=1); spacing=c2.number_input('Grid spacing (Å)',value=.5,min_value=.1)
st.subheader('Zernike channels')
st.checkbox('Shape',value=True,disabled=True,help='Geometric shape is the validated interface-plane channel.')
st.checkbox('Hydrophobicity',disabled=True,help='Experimental / not yet included in the interface plane.')
st.checkbox('Charge',disabled=True,help='Experimental / not yet included in the interface plane.')
st.checkbox('H-bond donor/acceptor',disabled=True,help='Experimental / not yet included in the interface plane.')

with st.expander('How does this work?'):
    st.code('PDB A        PDB B\n  │            │\nSurface      Surface\n  └── Interface ──┘\n       │\nFacing patches\n       │\n2D Zernike comparison\n       │\nInterface projection plane\n       │\nComplementarity heatmap')
    st.write('Surface points closer than 3 Å define the observed interface. Around sampled points, local 6 Å patches are built; A normals choose facing B patches. Rotationally invariant 2D Zernike descriptors quantify local geometry. The projection plane is an analysis coordinate system, not a claim that the real molecular interface is flat.')

if st.button('Run analysis',type='primary',disabled=not(atoms_a is not None and atoms_b is not None and ca and cb)):
    cfg=InterfaceMapConfig(float(interface_distance),float(patch_radius),int(sample_every),int(order),64,float(spacing))
    try:
        with st.status('Running geometric interface analysis...',expanded=True) as status:
            st.write('Generating surfaces...')
            with tempfile.TemporaryDirectory() as out:
                st.write('Detecting interface, building facing patches, and computing Zernike descriptors...')
                rows,meta,_,_=analyse_atoms(atoms_a,ca,atoms_b,cb,cfg,out)
                st.write('Building complementarity plane and figures...')
                files=result_files(out)
            status.update(label='Done',state='complete')
        st.session_state['zppi_result']={'rows':rows,'meta':meta,'files':files,'config':cfg,'chain_a':ca,'chain_b':cb,'source':'Automatic from PDB'}
    except ValueError as e:
        st.error(f'{e} Verify that structures are in their bound pose, selected chains are correct, or increase interface distance for exploration.')
    except Exception:
        st.error('Analysis failed during the scientific pipeline.')
        with st.expander('Technical details'): st.code(traceback.format_exc())

r=st.session_state.get('zppi_result')
if r:
    rows=r['rows']; meta=r['meta']; df=pd.DataFrame(rows)
    tabs=st.tabs(['Overview','Complementarity Map','3D Interface','Patch Pairs','Diagnostics'])
    with tabs[0]:
        vals=[meta['interface_count_a'],meta['interface_count_b'],len(rows),float(df.complementarity_score.mean()),float(df.complementarity_score.max())]
        for col,label,val in zip(st.columns(5),['Interface points A','Interface points B','Facing patch pairs','Mean complementarity','Best complementarity'],vals): col.metric(label,f'{val:.3f}' if isinstance(val,float) else val)
        st.caption(f"Chains {r['chain_a']} / {r['chain_b']} · interface distance {r['config'].interface_distance} Å · patch radius {r['config'].patch_radius} Å · order {r['config'].zernike_order}")
    with tabs[1]:
        st.image(r['files']['complementarity_map.png'],use_container_width=True)
    with tabs[2]: st.image(r['files']['interface_geometry.png'],use_container_width=True)
    with tabs[3]:
        cols=['patch_id','zernike_distance','complementarity_score','normal_opposition','axis_distance','suspicious_normals','x','y']; st.dataframe(df[cols].sort_values('complementarity_score',ascending=False),use_container_width=True)
    with tabs[4]:
        st.json({'interface_center_P':meta['interface_center'].tolist(),'global_axis_A_to_B':meta['axis'].tolist(),'mean_normal_dot':meta['normal_dot'],'suspicious_normal_pairs':int(df.suspicious_normals.sum()),'NaN_grid_cells':int(np.isnan(meta['matrix']).sum()),'surface_generation_method':r['source']})
    st.subheader('Download results')
    for name,data in r['files'].items(): st.download_button(name,data,file_name=name,key='download_'+name)
