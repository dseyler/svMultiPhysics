"""Generate a strongly spatially-heterogeneous per-element magnitude (scale) VTU
for the directionally_distributed_active_stress_scale test case.

`scale` is a dimensionless per-element multiplier on the active stress
(Ta_e = scale_e * Ta). It is bound to the element centroid along the z (fiber)
axis for a sharp two-region contrast:

  - centroid z below the mesh midplane: scale = 0.5  (attenuated)
  - centroid z above the mesh midplane: scale = 2.0  (amplified)

Both halves are active but with very different magnitude, producing an asymmetric
deformation. A wrong per-element scale mapping under MPI moves the strong region
onto the wrong cells and diverges from the serial reference.

Output CellData: Float64 `scale` + Int32 `GlobalElementID` (keyed consistently
with the mesh's global element numbering; the loader indexes by GlobalElementID-1).
The VTU contains only `scale` (no eta/delay arrays), exercising the
independent/optional path where eta falls back to the uniform domain values and
delay defaults to 0.

Run in an environment with meshio + numpy (e.g. the `svenv` conda env).
"""

import meshio
import numpy as np

mesh = meshio.read("mesh/mesh-complete.mesh.vtu")

cells = mesh.cells_dict["tetra"]              # (nEl, 4)
centroids = mesh.points[cells].mean(axis=1)   # (nEl, 3)
z = centroids[:, 2]
z_mid = 0.5 * (mesh.points[:, 2].min() + mesh.points[:, 2].max())

n_el = len(cells)
scale = np.where(z < z_mid, 0.5, 2.0).astype(np.float64)

if "GlobalElementID" in mesh.cell_data:
    gid = np.concatenate(mesh.cell_data["GlobalElementID"]).astype(np.int32)
else:
    gid = np.arange(1, n_el + 1, dtype=np.int32)

print(f"nEl={n_el}  scale=0.5:{int((scale == 0.5).sum())}  scale=2.0:{int((scale == 2.0).sum())}  z_mid={z_mid}")

out = meshio.Mesh(
    mesh.points,
    [("tetra", cells)],
    cell_data={"scale": [scale], "GlobalElementID": [gid]},
)
out.write("scale_heterogeneous.vtu", binary=True)
print("wrote scale_heterogeneous.vtu")
