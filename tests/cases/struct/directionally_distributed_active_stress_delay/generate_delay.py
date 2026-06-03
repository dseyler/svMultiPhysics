"""Generate a strongly spatially-heterogeneous per-element activation-delay VTU
for the directionally_distributed_active_stress_delay test case.

The `delay` is bound to the element centroid along the z (fiber) axis:

  - centroid z below the mesh midplane: delay = 0.0   (activates immediately)
  - centroid z above the mesh midplane: delay = 10.0  (>> final sim time, so these
                                                       elements stay at Ta = 0)

This isolates the delay feature (uniform eta = domain default) and makes the
spatial mapping observable: at the final compared step the lower half is actively
contracting while the upper half is passive (Ta = 0, via the zero-until-activation
override). Any incorrect per-element delay mapping under MPI moves the passive
region onto the wrong cells and diverges from the serial reference.

Output CellData: Float64 `delay` + Int32 `GlobalElementID` (keyed consistently
with the mesh's global element numbering; the loader indexes by GlobalElementID-1).
The VTU contains only `delay` (no eta arrays), exercising the independent/optional
path where eta falls back to the uniform domain values.

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
DELAY_BIG = 10.0                               # >> final sim time -> stays Ta == 0
delay = np.where(z < z_mid, 0.0, DELAY_BIG).astype(np.float64)

if "GlobalElementID" in mesh.cell_data:
    gid = np.concatenate(mesh.cell_data["GlobalElementID"]).astype(np.int32)
else:
    gid = np.arange(1, n_el + 1, dtype=np.int32)

print(f"nEl={n_el}  delay0={int((delay == 0.0).sum())}  delayed={int((delay > 0.0).sum())}  z_mid={z_mid}")

out = meshio.Mesh(
    mesh.points,
    [("tetra", cells)],
    cell_data={"delay": [delay], "GlobalElementID": [gid]},
)
out.write("delay_heterogeneous.vtu", binary=True)
print("wrote delay_heterogeneous.vtu")
