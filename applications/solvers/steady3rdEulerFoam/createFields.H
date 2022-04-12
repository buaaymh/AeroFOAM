// Gas information
IOdictionary testDict
(
    IOobject
    (
        "fluidProperties", // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::MUST_READ // indicate that reading this dictionary is compulsory
    )
);

word simulationType;
testDict.lookup("simulationType") >> simulationType;

const scalar gamma(testDict.lookupOrDefault<scalar>("gamma", 1.4));
const scalar PrDyna(testDict.lookupOrDefault<scalar>("PrDyna", 0.0));
const scalar PrTurb(testDict.lookupOrDefault<scalar>("PrTurb", 0.0));

Info << "=================Fluid Properties===================" << nl
     << "simulationType: "  << simulationType << nl
     << "gamma:          "  << gamma << nl
     << "PrDyna:         "  << PrDyna << nl
     << "PrTurb:         "  << PrTurb << nl
     << "====================================================" << endl;

const dictionary& freeStream = testDict.subDict("freeStream");
const scalar Mach = freeStream.lookup<scalar>("Mach");
const scalar T = freeStream.lookup<scalar>("T");
const scalar AOA = freeStream.lookup<scalar>("AOA");
const scalar Pressure = freeStream.lookupOrDefault<scalar>("Pressure", 101325.0);
const scalar Re = freeStream.lookupOrDefault<scalar>("Re", 1e6);
const scalar length = freeStream.lookupOrDefault<scalar>("length", 1.0);

Info << "===============Free-steam Definition================" << nl
     << "Temperature    : "  << T << " K" << nl
     << "Mach number    : "  << Mach << nl
     << "Angle of attack: "  << AOA  << " deg" << nl
     << "Pressure       : "  << Pressure  << " N/m2" << nl
     << "Reynolds number: "  << Re << nl
     << "Reynolds length: "  << length << " m" << nl
     << "====================================================" << endl;

Info << "==================Reading Fields====================" << endl;
// Primitive variables

Info<< "Reading field rho" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "====================================================" << endl;

EulerSolver solver
(
    mesh,
    rho,
    U,
    p
);

// Resident

volScalarField resRho
(
    IOobject
    (
        "resRho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("resRho", rho.dimensions(), 1.0)
);

volVectorField resRhoU
(
    IOobject
    (
        "resRhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("resRhoU", solver.rhoU().dimensions(), vector::one)
);

volScalarField resRhoE
(
    IOobject
    (
        "resRhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("resRhoE", solver.rhoE().dimensions(), 1.0)
);