"""Generate a strongly spatially-heterogeneous active-stress directional
distribution VTU for the directionally_distributed_active_stress_heterogeneous
test case.

The eta fractions are bound to the element centroid along the z (fiber) axis,
producing a sharp two-region split:

  - centroid z below the mesh midplane: (eta_f, eta_s, eta_n) = (0.1, 0.2, 0.7)
  - centroid z above the mesh midplane: (eta_f, eta_s, eta_n) = (0.7, 0.2, 0.1)

Both rows sum to 1.0 and are valid for the Holzapfel-Ogden modified-anisotropy
model used here. Because the field is sharply heterogeneous in space, any
incorrect per-element mapping (e.g. applying eta in the wrong post-partition
element order under MPI) moves stress onto the wrong cells and changes the
deformation. The accompanying regression test therefore fails if the per-element
distribution is mis-mapped in parallel.

Output CellData arrays: eta_f, eta_s, eta_n (Float64) and GlobalElementID
(Int32), in the same cell order as the mesh, keyed by the mesh GlobalElementID.

Run in an environment with meshio + numpy (e.g. the `svenv` conda env).
"""

import meshio
import numpy as np

mesh = meshio.read("mesh/mesh-complete.mesh.vtu")

# Single tetra cell block.
cells = mesh.cells_dict["tetra"]              # (nEl, 4)
centroids = mesh.points[cells].mean(axis=1)   # (nEl, 3)
z = centroids[:, 2]
z_mid = 0.5 * (mesh.points[:, 2].min() + mesh.points[:, 2].max())

low = z < z_mid
n_el = len(cells)

eta_f = np.where(low, 0.1, 0.7)
eta_s = np.full(n_el, 0.2)
eta_n = np.where(low, 0.7, 0.1)

# Use the mesh's own GlobalElementID so the values are keyed consistently with
# the solver's global element numbering (the loader indexes by GlobalElementID-1).
if "GlobalElementID" in mesh.cell_data:
    gid = np.concatenate(mesh.cell_data["GlobalElementID"]).astype(np.int32)
else:
    gid = np.arange(1, n_el + 1, dtype=np.int32)

print(f"nEl={n_el}  low={int(low.sum())}  high={int((~low).sum())}  z_mid={z_mid}")

out = meshio.Mesh(
    mesh.points,
    [("tetra", cells)],
    cell_data={
        "eta_f": [eta_f],
        "eta_s": [eta_s],
        "eta_n": [eta_n],
        "GlobalElementID": [gid],
    },
)
out.write("eta_heterogeneous.vtu", binary=True)
print("wrote eta_heterogeneous.vtu")
