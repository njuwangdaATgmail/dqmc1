I am designing a DQMC code with only one effective component of fermions. All the other components can be related to the working component through duplicate or symmetry transformations, e.g. particle-hole transformation. In fact in principle, this single component program can also be applied to multi component cases as long as we view all fermions as one "big component".

Phonon or gauge fields should be able to be added, so the auxiliary fields can be both descrete and continuous. Besides, the auxiliary fields can couple to both onsite and offsite fermionic billinear operators whose dimension can be arbitrary.

Interactions already been added on square lattice:
  + U, V, J
  + onsite Holstein phonon
  + onbond breathing phonon
  + onbond buckling phonon coupling to density
  + onbond buckling phonon coupling to hoppinp

How to check the code?
  + 2-site model without phonon, by exact diagonalization
  + 1-site model (or t=0 in 2-site model), by Lang-Firsov transformation
  + 2-site model, see Alexandrov, 1994'? (only T=0)
  + Johnston, 2013'.

Todo:
  + remove unnecessary parts outside of the core. For example: init_phi, init_ising, generate_new_ising, generate_new_phi.
