"""Generate the mesh domain split and two distinct per-domain active-stress
directional distributions for the directionally_distributed_active_stress_multidomain
test case.

The 459-element block mesh is split into two material domains along the z
(fiber) axis:

  - domain 1: element centroid z below the mesh midplane  (DOMAIN_ID = 1)
  - domain 2: element centroid z at/above the mesh midplane (DOMAIN_ID = 2)

Each domain is given its OWN spatial directional distribution, written as a
WHOLE-MESH VTU (per the loader convention: keyed by the mesh GlobalElementID,
covering all gnEl elements; only that domain's elements consume its file). The
two fields are deliberately DIFFERENT and both vary in space (with x), so:

  - a multi-domain *routing* bug (one domain's distribution clobbering the
    other, the bug this test guards) sends the wrong eta to a domain's elements,
    and
  - a per-element *ordering* bug under MPI moves eta onto the wrong cells,

either of which changes the deformation and fails the regression test against
the single serial reference.

eta_s is held at 0.2; the remaining 0.8 is split between fiber and sheet-normal
by a position factor t in [0,1] (eta_f = 0.8 t, eta_n = 0.8 (1-t)), so every row
sums to 1.0 and is valid for the Holzapfel-Ogden modified-anisotropy model:

  - domain 1 file: t = x_frac        (fiber fraction grows with x)
  - domain 2 file: t = 1 - x_frac    (fiber fraction shrinks with x)

Run in an environment with meshio + numpy (e.g. the `svenv` conda env).
"""

import meshio
import numpy as np

mesh = meshio.read("mesh/mesh-complete.mesh.vtu")
# Work in float64 throughout: the mesh points are Float32, so float32 arithmetic
# would leave eta_f + eta_s + eta_n off 1.0 by ~1e-7 and fail the solver's
# 1e-10 sum-to-1 check.
pts = mesh.points.astype(np.float64)
cells = mesh.cells_dict["tetra"]                 # (nEl, 4)
centroids = pts[cells].mean(axis=1)              # (nEl, 3)
n_el = len(cells)

if "GlobalElementID" in mesh.cell_data:
    gid = np.concatenate(mesh.cell_data["GlobalElementID"]).astype(np.int32)
else:
    gid = np.arange(1, n_el + 1, dtype=np.int32)

# --- domain split along z ---
z = centroids[:, 2]
z_mid = 0.5 * (pts[:, 2].min() + pts[:, 2].max())
domain_id = np.where(z < z_mid, 1, 2).astype(np.int32)
print(f"nEl={n_el}  domain1={(domain_id==1).sum()}  domain2={(domain_id==2).sum()}  z_mid={z_mid}")

# DOMAIN_ID file (whole mesh; value = domain id = bit position matching <Domain id>).
meshio.Mesh(pts, [("tetra", cells)],
            cell_data={"DOMAIN_ID": [domain_id]}).write("mesh/domain_ids.vtu", binary=True)
print("wrote mesh/domain_ids.vtu")

# --- per-domain directional distributions (whole-mesh, distinct, x-varying) ---
x = centroids[:, 0]
x_frac = (x - x.min()) / (x.max() - x.min())     # in [0, 1]

def write_eta(path, t):
    # Float64 (the solver's copy_cell_data silently skips non-Float64 CellData),
    # with eta_n the exact remainder so each row sums to 1.0 to machine precision.
    t = np.asarray(t, dtype=np.float64)
    eta_s = np.full(n_el, 0.2, dtype=np.float64)
    eta_f = 0.8 * t
    eta_n = 1.0 - eta_f - eta_s
    meshio.Mesh(pts, [("tetra", cells)],
        cell_data={
            "eta_f": [eta_f], "eta_s": [eta_s], "eta_n": [eta_n],
            "GlobalElementID": [gid],
        }).write(path, binary=True)
    print(f"wrote {path}  (eta_f range {eta_f.min():.3f}..{eta_f.max():.3f})")

# Keep t in [0.1, 0.9] so every fraction is comfortably positive (avoids tiny
# negative round-off at the endpoints that the validator would reject).
write_eta("eta_domain1.vtu", 0.1 + 0.8 * x_frac)   # fiber fraction grows with x
write_eta("eta_domain2.vtu", 0.9 - 0.8 * x_frac)   # fiber fraction shrinks with x
