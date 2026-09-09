from pathlib import Path
import numpy as np
import pytest
from zernike_ppi.ui import InterfaceMapConfig, chains, extract_chain, read_pdb_upload, analyse_atoms

ROOT=Path(__file__).parent/'data'
def test_two_pdb_input_and_chain_extraction():
 a=read_pdb_upload(ROOT/'fixture_a.pdb'); b=read_pdb_upload(ROOT/'fixture_b.pdb')
 assert chains(a)==['A'] and chains(b)==['B']
 assert len(extract_chain(a,'A')['xyz'])==4
 with pytest.raises(ValueError): extract_chain(a,'Z')
def test_one_complex_pdb_chain_extraction():
 complex_atoms=read_pdb_upload(ROOT/'fixture_complex.pdb')
 assert chains(complex_atoms)==['A','C']
 assert len(extract_chain(complex_atoms,'A')['xyz'])==2
 assert len(extract_chain(complex_atoms,'C')['xyz'])==2
def test_config_and_no_interface():
 cfg=InterfaceMapConfig(); assert cfg.interface_distance==3 and cfg.zernike_order==20
 a=read_pdb_upload(ROOT/'fixture_a.pdb'); b=read_pdb_upload(ROOT/'fixture_b.pdb')
 b={k:(v+100 if k=='xyz' else v) for k,v in b.items()}
 with pytest.raises(ValueError,match='No interface'):
  analyse_atoms(a,'A',b,'B',cfg,Path('tests/_tmp_ui'))
