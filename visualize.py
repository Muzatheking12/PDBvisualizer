import vtk
import os
from Bio.PDB import PDBParser, PDBIO

# Coded by Muzammil Kabier
# Directory for extracted ligands
hetlig = os.path.join(os.path.dirname(__file__), r'lig')
# Change r'2v5z.pdb to your PDB Filename'
protein = os.path.join(os.path.dirname(__file__), r'2v5z.pdb')

# Function to check if a residue is a heteroatom
def is_het(residue):
    return residue.id[0] not in (" ", "W")

# Function to extract ligands from PDB (Reside with 4 and less characters will be considered as IONS and more will be treated as other HETATM)
def extract_ligands_from_pdb(pdb_path, output_directory):
    pdb_file = os.path.basename(pdb_path)
    pdb_id = os.path.splitext(pdb_file)[0].upper()
    pdb = PDBParser().get_structure(pdb_id, pdb_path)
    
    io = PDBIO()
    
    for model in pdb:
        for chain in model:
            for residue in chain:
                if is_het(residue):
                    chain_id = chain.get_id()
                    ligand_resname = residue.get_resname()
                    ligand_filename = f"{ligand_resname}_{chain_id}.pdb"
                    ligand_path = os.path.join(output_directory, ligand_filename)

                    io.set_structure(residue)
                    io.save(ligand_path)
                    print(f"Ligand extracted: {ligand_resname} in Chain {chain_id}")

def pdb_to_vtk_actor(pdb_file):
            """Determine if the molecule is an ion and render accordingly."""
            filename = os.path.basename(pdb_file)
            ligand_name = os.path.splitext(filename)[0]  # Remove extension

            if len(ligand_name) <= 4:  # IONS usually has 2 string character like for eg. CL,MG 
                return create_sphere_actor()
            else:
                return create_tube_actor(pdb_file)

def create_sphere_actor():
            """Create a VTK Sphere for ions."""
            sphere_source = vtk.vtkSphereSource()
            sphere_source.SetRadius(3.0)  # Adjust size as needed
            sphere_source.SetThetaResolution(20)
            sphere_source.SetPhiResolution(20)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere_source.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(0, 1, 0)  # Green Colour (CHANGE AS YOU WISH FOR IONS)
            return actor

def create_tube_actor(pdb_file):
            """Create a VTK Tube for normal ligands."""
            reader = vtk.vtkPDBReader()
            reader.SetFileName(pdb_file)
            reader.Update()

            tube_filter = vtk.vtkTubeFilter()
            tube_filter.SetInputConnection(reader.GetOutputPort())
            tube_filter.SetRadius(0.3)
            tube_filter.SetNumberOfSides(20)
            tube_filter.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(tube_filter.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            return actor

# Function to load ligands into the scene (HETATM PDB IS STORED IN FOLDER 'lig')
def load_pdb_files(directory):
    actors = []
    for file in os.listdir(directory):
        if file.endswith(".pdb"):
            pdb_path = os.path.join(directory, file)
            actor= pdb_to_vtk_actor(pdb_path)
            actors.append(actor)
            
    return actors



    


# Load initial PDB file (Change as You Wish)
initial_pdb = protein
mainreader = vtk.vtkPDBReader()
mainreader.SetFileName(protein)
ribbonfilter = vtk.vtkProteinRibbonFilter()
ribbonfilter.SetInputConnection(mainreader.GetOutputPort())
ribbonfilter.SetDrawSmallMoleculesAsSpheres(False) # Hides HETATM visualized as SPHERES
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(ribbonfilter.GetOutputPort())
renderer = vtk.vtkRenderer()
actorx = vtk.vtkActor()
actorx.SetMapper(mapper)
renderer.AddActor(actorx)
for file in os.listdir(hetlig):
      pathway = os.path.join(hetlig, file)
      if os.path.exists(pathway):
            os.remove(pathway) # REMOVE HETATM BEFORE GENERATION
extract_ligands_from_pdb(initial_pdb, hetlig)
actors = load_pdb_files(hetlig)
for ligand_actor in actors:
    renderer.AddActor(ligand_actor)


renderer.SetBackground(0.1, 0.1, 0.1)

# Create the render window
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
render_window.SetSize(800, 600)

# Create interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# Start the interaction loop
interactor.Initialize()
render_window.Render()
interactor.Start()
