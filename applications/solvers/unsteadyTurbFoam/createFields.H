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

fluidProperties fluidProps;

testDict.lookup("simulationType") >> fluidProps.simulationType;
fluidProps.withSourceTerm = testDict.lookupOrDefault<bool>("withSourceTerm", false);

Info << "=================Fluid Properties===================" << nl
     << "simulationType: "  << fluidProps.simulationType << nl
     << "====================================================" << endl;

const dictionary& freeStream = testDict.subDict("freeStream");
fluidProps.AOA = freeStream.lookup<scalar>("AOA");
fluidProps.Mach_inf = freeStream.lookupOrDefault<scalar>("Mach", 1.0);
fluidProps.Re_inf = freeStream.lookupOrDefault<scalar>("Re", 10000);
fluidProps.T_inf = freeStream.lookupOrDefault<scalar>("T", 288);
fluidProps.p_inf = freeStream.lookupOrDefault<scalar>("Pressure", 101325.0);

Info << "===============Free-stream Definition================" << nl
     << "Temperature    : "  << fluidProps.T_inf << " K" << nl
     << "Mach number    : "  << fluidProps.Mach_inf << nl
     << "Angle of attack: "  << fluidProps.AOA  << " deg" << nl
     << "Pressure       : "  << fluidProps.p_inf  << " N/m2" << nl
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

Info<< "Reading field nuTilda" << endl;
volScalarField nuTilda
(
    IOobject
    (
        "nuTilda",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info << "====================================================" << endl;

std::unique_ptr<Source> source;
if (fluidProps.withSourceTerm) source = std::make_unique<Source>(rho, U);

std::unique_ptr<turbulenceSolver> solver;
label order = mesh.schemesDict().lookup<label>("order");

if (order == 2) solver = std::make_unique<turbulence2ndSolver>(fluidProps, rho, U, p, nuTilda);
else if (order == 3) solver = std::make_unique<turbulence3rdSolver>(fluidProps, rho, U, p, nuTilda);

solver->correctFields();
// Resident
scalarField resRho_1(mesh.nCells());
vectorField resRhoU_1(mesh.nCells());
scalarField resRhoE_1(mesh.nCells());
scalarField resNuTilda_1(mesh.nCells());