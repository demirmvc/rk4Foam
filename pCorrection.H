U.correctBoundaryConditions();
solve(pEqn == fvc::div(U)/runTime.deltaT());
#include "continuityErrs.H"
U = U - fvc::grad(p)*runTime.deltaT();
U.correctBoundaryConditions();
