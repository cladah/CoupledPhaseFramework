from dolfinx.io import XDMFFile, gmshio
from dolfinx import io
from mpi4py import MPI
import gmsh
def convert_mesh():
    msh, cell_markers, facet_markers = gmshio.read_from_msh("Resultfiles/Mesh.msh", MPI.COMM_WORLD, 0)
    msh.name = 'Sphere'
    cell_markers.name = f"{msh.name}_cells"
    facet_markers.name = f"{msh.name}_facets"
    with io.XDMFFile(MPI.COMM_WORLD, "Resultfiles/Mesh.xdmf", "w") as file:
        msh.topology.create_connectivity(1, 2)
        file.write_mesh(msh)
        file.write_meshtags(cell_markers)
        file.write_meshtags(facet_markers)
        # bottom_facets = mesh.locate_entities_boundary(domain, domain.topology.dim-1, bottom)