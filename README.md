# PDBvisualizer
<b>Generate 3D- Ribbon visualization of PDB File with VTK which can be incorporated into your Softwares.</b>

<img src=".\protman.png">

## Details
<li>Beta Sheet - Coloured to Yellow</li>
<li>Alpha Helix - Coloured to Purple</li>
<li>Heteroatoms other than IONS, is shown as sticks with O,N,S and Halogens coloured</li>
<li>IONS are Spheres in GREEN colour can be changed accordingly</li>

## Dependency
The 'visualize.py' Python File requires VTK and Biopython <br>
<b>- pip install vtk</b><br>
<b>- pip install biopython</b>

## NOTE
Please create folder named 'lig' in the current directory . This folder will intake the HETATM files
