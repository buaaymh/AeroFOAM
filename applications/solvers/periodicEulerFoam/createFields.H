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

// Initialize Condition
// rho = 1.0 + 0.2*sin((x+y)*pi)
// forAll(mesh.C(), cellI)
// {
//     std::vector<scalar> weights;
//     std::vector<vector> quadPoints;
//     build2ndCell(mesh, cellI, weights, quadPoints);
//     scalar rhoTemp = 0;
//     for (label i = 0; i != weights.size(); ++i)
//     {
//         scalar x = quadPoints[i].x();
//         scalar y = quadPoints[i].y();
//         scalar temp = 1.0 + 0.2*Foam::sin((x+y)*constant::mathematical::pi);
//         rhoTemp  += temp * weights[i];
//     }
//     rhoTemp  /= mesh.V()[cellI];
//     rho[cellI] = rhoTemp;
// }
// scalarField rhoExact(rho.primitiveField());

// vortex evolution problem
scalar beta = 5.0;
scalar gamma = fluidProps.gamma;
scalar pi = constant::mathematical::pi;
forAll(mesh.C(), cellI)
{
    std::vector<scalar> weights;
    std::vector<vector> quadPoints;
    build2ndCell(mesh, cellI, weights, quadPoints);
    scalar rhoTemp  = 0;
    vector rhoUTemp = vector::zero;
    scalar rhoETemp = 0;
    for (label i = 0; i != weights.size(); ++i)
    {
        scalar x = quadPoints[i].x();
        scalar y = quadPoints[i].y();
        scalar r2 = sqr(x-5.0) + sqr(y-5.0);
        scalar TLocal   = 1.0 - (gamma-1.0)*sqr(beta)*Foam::exp(1.0-r2)/(8.0*gamma*sqr(pi));
        scalar rhoLocal = Foam::pow(TLocal, 1.0/(gamma-1.0));
        vector ULocal   = vector(1.0 - beta*Foam::exp((1.0-r2)/2.0)*(y-5.0)/(2.0*pi),
                                 1.0 + beta*Foam::exp((1.0-r2)/2.0)*(x-5.0)/(2.0*pi), 0);
        scalar pLocal   = TLocal*rhoLocal;
        rhoTemp  += rhoLocal * weights[i];
        rhoUTemp += rhoLocal * ULocal * weights[i];
        rhoETemp += (pLocal/(gamma-1.0) + 0.5*rhoLocal*magSqr(ULocal)) * weights[i];
    }
    rhoTemp  /= mesh.V()[cellI];
    rhoUTemp /= mesh.V()[cellI];
    rhoETemp /= mesh.V()[cellI];
    rho[cellI] = rhoTemp;
    U[cellI]   = rhoUTemp/rhoTemp;
    p[cellI]   = (rhoETemp-0.5*magSqr(rhoUTemp)/rhoTemp)*(gamma-1.0);
}
scalarField rhoExact(mesh.nCells());
forAll(mesh.C(), cellI)
{
    std::vector<scalar> weights;
    std::vector<vector> quadPoints;
    build2ndCell(mesh, cellI, weights, quadPoints);
    scalar rhoTemp  = 0;
    for (label i = 0; i != weights.size(); ++i)
    {
        scalar x = quadPoints[i].x();
        scalar y = quadPoints[i].y();
        scalar r2 = sqr(x-7.0) + sqr(y-7.0);
        scalar TLocal   = 1.0 - (gamma-1.0)*sqr(beta)*Foam::exp(1.0-r2)/(8.0*gamma*sqr(pi));
        scalar rhoLocal = Foam::pow(TLocal, 1.0/(gamma-1.0));
        rhoTemp  += rhoLocal * weights[i];
    }
    rhoTemp  /= mesh.V()[cellI];
    rhoExact[cellI] = rhoTemp;
}


rho.correctBoundaryConditions();
U.correctBoundaryConditions();
p.correctBoundaryConditions();

Info << "====================================================" << endl;

std::unique_ptr<solver> solver;
label order = mesh.schemesDict().lookup<label>("order");
if (fluidProps.simulationType == "Euler")
{
    if (order == 2) solver = std::make_unique<euler2ndSolver>(fluidProps, rho, U, p);
    else if (order == 3)
    {
        solver = std::make_unique<euler3rdSolver>(fluidProps, rho, U, p);
    }
}
else
{
    Info << "Error in simulationType!" << endl;
}
// Resident
scalarField resRho_1(mesh.nCells());
vectorField resRhoU_1(mesh.nCells());
scalarField resRhoE_1(mesh.nCells());
