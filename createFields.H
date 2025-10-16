#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<twoPhaseSystem> fluidPtr
(
    twoPhaseSystem::New(mesh)
);
twoPhaseSystem& fluid = fluidPtr();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = fluid.phase1().thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());

dimensionedScalar delay("delay", dimTime, 0.0); 

FlotationModel flotation (mesh, fluid, delay);


// volScalarField Q
// (
//     IOobject
//     (
//         "Q",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Q", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0) 
// );




// Info << "Reading flotationProperties" << endl;
// IOdictionary flotationProperties
// (
//     IOobject
//     (
//         "flotationProperties",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ_IF_MODIFIED,
//         IOobject::NO_WRITE
//     )
// );

// dimensionedScalar rhoP
// (
//     "rho0", flotationProperties
// );

// dimensionedScalar C_1
// (
//     "C_1", flotationProperties
// );

// dimensionedScalar d_p
// (
//     "d_p", flotationProperties
// );

// dimensionedScalar d_b
// (
//     "d_b", flotationProperties
// );

// dimensionedScalar nu_F
// (
//     "nu_F", flotationProperties
// );

// dimensionedScalar sigma
// (
//     "sigma", flotationProperties
// );

// dimensionedScalar theta_c
// (
//     "theta_c", flotationProperties
// );

// dimensionedScalar Z_0 = Foam::sqrt((8.0 * M_PI) / (15.0 * nu_F)) * (((d_p + d_b) / 2) * ((d_p + d_b) / 2) * ((d_p + d_b) / 2));

// dimensionedScalar Z_prime_0 = Foam::sqrt(C_1) * Foam::pow(d_p + d_b, -2.0/3.0);

// dimensionedScalar Nmax = ((d_b / d_p) * (d_b / d_p)) / 2;

// dimensionedScalar E_c = (3 / 2) * ((d_p / d_b) * (d_p / d_b));

// dimensionedScalar t_i("t_i", dimensionSet(0, 0, 1, 0, 0, 0, 0), (75.0 / theta_c.value()) * Foam::pow(d_p.value(), 0.6)); 

// theta_c = M_PI * theta_c / 180.0;

// // dimensionedScalar my_epsilon("my_epsilon", dimensionSet(0, 2, -3, 0, 0, 0, 0), 1.5e-4);

// dimensionedScalar A_s = 0.5;

// volScalarField KA
// (
//     IOobject
//     (
//         "KA",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("KA", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0) 
// );

// volScalarField KD
// (
//     IOobject
//     (
//         "KD",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("KD", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0) 
// );

// volScalarField E_s
// (
//     IOobject
//     (
//         "E_s",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("E_s", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0) 
// );

// volScalarField oneDivBo_prime
// (
//     IOobject
//     (
//         "oneDivBo_prime",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("oneDivBo_prime", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0) 
// );

// volScalarField arY1
// (
//     IOobject
//     (
//         "arY1",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("arY1", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0) 
// );

// volScalarField arY2
// (
//     IOobject
//     (
//         "arY2",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("arY2", dimensionSet(1, -3, 0, 0, 0, 0, 0), 0.0) 
// );
