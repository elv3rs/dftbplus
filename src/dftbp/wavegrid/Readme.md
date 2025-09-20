PR: waveplot molorb kernel rewrite
- Separate Molorb Kernel from waveplot
- Move molorb kernel to 
- Move molorb kernel from app/waveplot to  
- Rewrite Molorb Kernel to employ OMP CPU parallelisation
Allow inplace total charge accumulation
Add optional gpu offloaded (cuda) kernel implementation
Basis:
    Allow supplyying LUT for the radial part instead of sto parameters
    Expand tesserals to l=4



# Wavegrid
This is a molecular orbital calculator.
How does it work?
## GridParams
You declare the calculation grid, specifying grid points, basis vectors, and origin.
## System
Specify a number of atoms, their positions, and species type (H, He,...)
Specify a mapping of type to a range of basis functions. (e.g. H->1s)
## Periodic
is the systemperiodic?
If so, here go the latice vectors. They are used to fold coords into the unit cell.
If you have complex eigenvectors, you must also specify the kPoints and kIndices.
## Basis
You provide the basis set as STOs.
sto = radial(r) * Y_lm(theta, phi)
The radial part is calculated as a linear combination of polynomials multiplied by exponentials.
radial = \sum_i^{n_alpha} [ \sum_j^{n_pow} a_{ij} r^{l + j -1} ] exp(-alpha_i r)

## Calculation Parameters
You provide the eigenvector coefficients for each orbital and state, and are faced with three options:

Regular calculation:
Calculate the wavefunction summed over all orbitals at each grid point, for each state.
Psi(x,y,z,iEig) = \sum_atoms \sum_species orbitals c_{iOrbital, iEig} * sto_iOrbital(x,y,z)
                 = sum_{iOrbital} c_{iOrbital, iEig} * b_{iOrbital}(x,y,z) 

AtomicDensities Calculation:
Calculate the wavefunction squared at each grid point, for each state.
Psi(x,y,z,iEig) = sum_{iOrbital} c_{iOrbital, iEig} * [b_{iOrbital}(x,y,z)]^2

TotalCharge Calculation:
Calculate the total charge density at each grid point, summing over all states.
You provide the occupation vector f_i.
Psi(x,y,z) = \sum_{iEig} f_i |Psi(x,y,z,iEig)|^2






