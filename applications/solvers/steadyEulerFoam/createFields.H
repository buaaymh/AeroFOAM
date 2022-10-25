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

fluidProps.gamma = testDict.lookupOrDefault<scalar>("gamma", 1.4);

Info << "=================Fluid Properties===================" << nl
     << "simulationType: "  << fluidProps.simulationType << nl
     << "gamma:          "  << fluidProps.gamma << nl
     << "====================================================" << endl;

const dictionary& freeStream = testDict.subDict("freeStream");
fluidProps.AOA = freeStream.lookup<scalar>("AOA");
fluidProps.Mach_inf = freeStream.lookup<scalar>("Mach");
fluidProps.T_inf = freeStream.lookup<scalar>("T");
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
rho.correctBoundaryConditions();
U.correctBoundaryConditions();
p.correctBoundaryConditions();

// File pointer to direct the output to
autoPtr<OFstream> outputFilePtr;
if (Pstream::master())
{
    // Create the output path directory
    fileName outputDir = mesh.time().path()/"postProcessing";
    // Creathe the directory
    mkDir(outputDir);
    // Open the file in the newly created directory
    outputFilePtr.reset(new OFstream(outputDir/"residuals.dat"));
    outputFilePtr() << "# Residuals" << endl;
    outputFilePtr() << "# Step" << tab << "CPU" << tab << "L2(resRho)" << endl;
}

Info << "====================================================" << endl;

std::unique_ptr<solver> solver;
label order = mesh.schemesDict().lookup<label>("order");
if (fluidProps.simulationType == "Euler")
{
    if (order == 2) solver = std::make_unique<euler2ndSolver>(fluidProps, rho, U, p);
    else if (order == 3) solver = std::make_unique<eulerPrimVar3rdSolver>(fluidProps, rho, U, p);
}
else
{
    Info << "Error in simulationType!" << endl;
}

solver->correctFields();
// Resident
scalarField resRho(mesh.nCells());
vectorField resRhoU(mesh.nCells());
scalarField resRhoE(mesh.nCells());