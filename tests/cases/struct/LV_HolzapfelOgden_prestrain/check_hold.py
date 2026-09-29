"""Equilibrium hold test for the imprinted-deformation-gradient prestrain.

After prestrain.xml has imprinted the 10 mmHg state and forward.xml has run
at the same 10 mmHg, the forward run must not move: its displacement stays at
zero and its stress equals the stress the prestrain run ended with.

Usage, from this directory:
    mpirun -np 1 <path>/svmultiphysics prestrain.xml
    mpirun -np 1 <path>/svmultiphysics forward.xml
    python check_hold.py [prestrain_last.vtu] [forward_last.vtu]
"""
import glob
import os
import sys

import numpy as np
import pyvista as pv

# Tolerances, relative. Displacements in the VTU are in the mesh file's units
# (metres here; the solver works in cm via Mesh_scale_factor), so the hold
# displacement is measured against the mesh's bounding-box diagonal.
DISP_REL_TOL = 1.0e-4
STRESS_REL_TOL = 1.0e-3


def last_result(folder):
    files = sorted(glob.glob(os.path.join(folder, "result_*.vtu")))
    if not files:
        sys.exit(f"no result_*.vtu in {folder}")
    return files[-1]


def main():
    pre_file = sys.argv[1] if len(sys.argv) > 1 else last_result("prestrain")
    fwd_file = sys.argv[2] if len(sys.argv) > 2 else last_result("forward")

    pre = pv.read(pre_file)
    fwd = pv.read(fwd_file)

    # The prestrain run's last displacement is that step's increment, which
    # should already be small if the imprint converged.
    b = np.array(fwd.bounds)
    size = np.linalg.norm(b[1::2] - b[0::2])   # bounding-box diagonal, mesh units
    pre_disp = np.linalg.norm(pre.point_data["Displacement"], axis=1).max() / size
    fwd_disp = np.linalg.norm(fwd.point_data["Displacement"], axis=1).max() / size

    s_pre = pre.point_data["Stress"]
    s_fwd = fwd.point_data["Stress"]
    stress_scale = np.abs(s_pre).max()
    stress_diff = np.abs(s_fwd - s_pre).max() / stress_scale

    has_U = "Prestrain_displacement" in pre.point_data
    U_max = np.linalg.norm(pre.point_data["Prestrain_displacement"], axis=1).max() if has_U else float("nan")

    print(f"prestrain run : {pre_file}")
    print(f"  mesh bounding-box diagonal       : {size:.4e} (mesh units)")
    print(f"  last-step max nodal displacement : {pre_disp:.3e} of the mesh size")
    print(f"  max |Prestrain_displacement|     : {U_max:.3e} (mesh units)")
    print(f"forward run   : {fwd_file}")
    print(f"  max nodal displacement           : {fwd_disp:.3e} of the mesh size   (tol {DISP_REL_TOL:.0e})")
    print(f"  max |Stress - Stress_prestrain|  : {stress_diff:.3e} relative   (tol {STRESS_REL_TOL:.0e})")
    print(f"  max |Stress|                     : {stress_scale:.4e} dyne/cm^2")

    ok = True
    if not has_U:
        print("FAIL: the prestrain result carries no Prestrain_displacement point array")
        ok = False
    if fwd_disp > DISP_REL_TOL:
        print("FAIL: the prestrained geometry moved under the load it was prestrained at")
        ok = False
    if stress_diff > STRESS_REL_TOL:
        print("FAIL: the forward stress does not reproduce the prestrain stress")
        ok = False

    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
