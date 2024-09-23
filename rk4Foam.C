/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow,"
        " using the RK4 algorithm."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Time loop
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "CourantNo.H"

        piso.read();
        Info<< "Time = " << runTime.timeName() << nl << endl;
        //tensor A(tensor::zero);    
        std::array<std::array<double, 4>, 4> A = {{
            {0.5, 0.0, 0.0, 0.0},      
            {0., 0.5, 0.0, 0.0},       
            {0.0, 0., 1., 0.0},       
            {0.16666666666, 0.33333333333, 0.33333333333, 0.16666666666}  
        }};
        std::array<double, 4> ct = {0.5, 0.5, 1.0, 1.0};

        Uold = U;

        // 1st Stage
        phi = fvc::interpolate(U) & mesh.Sf();
        const volVectorField dU1 = runTime.deltaT() * (fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi, U)); //intermediate velocity k1
        // use this velocity while calculating k2. y +0.5*k1 
        U = Uold + (A[0][0]) * dU1; 
        U.correctBoundaryConditions();
        solve(fvm::laplacian(p) == (1.0 / (ct[0] * runTime.deltaT())) * fvc::div(U)); 
        #include "continuityErrs.H"
        U = U - (ct[0]) * runTime.deltaT() * fvc::grad(p); // div free velocity
        U.correctBoundaryConditions();

        // 2nd Stage
        phi = fvc::interpolate(U) & mesh.Sf();
        const volVectorField dU2 = runTime.deltaT() * (fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi, U)); //k2 
         // use this velocity while calculating k2. y +0.5*k2
        U = Uold + (A[1][1]) * dU2; 
        U.correctBoundaryConditions();
        solve(fvm::laplacian(p) == (1.0 / (ct[1] * runTime.deltaT())) * fvc::div(U));
        #include "continuityErrs.H"
        U = U - ct[1] * runTime.deltaT() * fvc::grad(p);
        U.correctBoundaryConditions();

        // 3rd Stage
        phi = fvc::interpolate(U) & mesh.Sf();
        const volVectorField dU3 = runTime.deltaT() * (fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi, U)); // k3
        // use this velocity while calculating k2. y + k3
        U = Uold + (A[2][2]) * dU3; //k4
        U.correctBoundaryConditions();
        solve(fvm::laplacian(p) == (1.0 / ((ct[2]) * runTime.deltaT())) * fvc::div(U));
        #include "continuityErrs.H"
        U = U - ct[2] * runTime.deltaT() * fvc::grad(p);
        U.correctBoundaryConditions();

        // 4th Stage
        phi = fvc::interpolate(U) & mesh.Sf();
        const volVectorField dU4 = runTime.deltaT() * (fvc::laplacian(turbulence->nuEff(), U) - fvc::div(phi, U)); //k4
        // U_(n+1) = U_n + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4
        U = Uold + (A[3][0]) * dU1 + (A[3][1]) * dU2 + (A[3][2]) * dU3 + (A[3][3]) * dU4;
        U.correctBoundaryConditions();
        solve(fvm::laplacian(p) == (1.0 / ((ct[3]) * runTime.deltaT())) * fvc::div(U));
        #include "continuityErrs.H"
        U = U - ct[3] * runTime.deltaT() * fvc::grad(p);
        U.correctBoundaryConditions();

        phi = fvc::interpolate(U) & mesh.Sf();

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //