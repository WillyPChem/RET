# RET
Simulator for resonant energy transfer in photosynthetic light harvesting complexes

CoupledDensityMatrix.c simulates resonant energy transfer (via dipole-dipole coupling) between two systems described by 
two independent hamiltonians.  The system designated the Nanoparticle can have an arbitrary number of states...
the user will be prompted for the number of states at runtime, and this should be consistent with
the Hamiltonian matrices provided in the MATRICES directory).
Currently, the system designated the molecule (nicknamed MG for malachite green) is hardcoded to have three states.
This can be changed fairly easily in the code if needed.

There is a subdirectory called UNCOUPLED with a file called NP.c that contains code for a single system (no coupling).  
It is a good check on the Hamiltonian parameters to run the uncoupled code first.  

- Compling

-- Both codes come with a makefile that will properly link the code to the Fastest Fourier Transform in the West library
-- To comple simply type 'make' 
-- The coupled code will compile to an executable called 'Coupled.x'
-- the uncoupled code will compile to an executable called 'NP.x'
