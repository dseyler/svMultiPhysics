# Prestraining the LV by an imprinted deformation gradient

Two-phase example of `<Prestrain>`, the modified updated
Lagrangian prestraining of Gee, Förster & Wall (Int. J. Numer. Meth. Biomed.
Engng. 26, 2010, section 3), on the `LV_HolzapfelOgden_passive` geometry at an
endocardial pressure of 10 mmHg (13332.2 dyne/cm²).

The imaged geometry is taken to be the loaded configuration. Its coordinates
never change; instead each converged load step imprints the deformation
gradient it reached into every Gauss point, `F0`, and the next step starts from
`F = F0 + Grad(u)` with the displacement reset to zero. The state converges to
one that carries the load without deforming.

```
mpirun -np 1 <build>/bin/svmultiphysics prestrain.xml    # writes prestrain/result_NNN.vtu
mpirun -np 1 <build>/bin/svmultiphysics forward.xml      # reads prestrain/result_012.vtu
python check_hold.py                                     # equilibrium hold test
```

`prestrain.xml` applies the full load from the first step as a dead load
(`Follower_pressure_load false`, as in the paper) and prints
`Prestrain: max nodal displacement` every step. With
`Prestrain_adaptive_time_step` the pseudo time step grows as the state
approaches equilibrium and the displacement reaches round-off in 6
steps; at a fixed dt = 1e-2 it takes about 160. The last result carries the
imprint as cell arrays `Prestrain_F_g<g>`, one per Gauss point.

`forward.xml` loads that imprint through `<Prestrain_file_path>`
and holds the same 10 mmHg as a follower load. `check_hold.py` asserts that the
geometry does not move and that the reported stress reproduces the prestrain
stress. Raise the pressure in `forward.xml` to simulate the prestrained LV; run
it without the file path to see the naive inflation for contrast.

The forward run may be started on any number of ranks regardless of how many
the prestrain run used: the imprint is written in the mesh's original element
order and partitioned on input like the fibre data.
