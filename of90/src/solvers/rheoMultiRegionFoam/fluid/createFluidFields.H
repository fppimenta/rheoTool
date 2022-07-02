// Initialise fluid field pointer lists

// Fields
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> CFluid(fluidRegions.size());
 
// Models
PtrList<constitutiveModel> constEqFluid(fluidRegions.size());
PtrList<EDFModel> elecMFluid(fluidRegions.size());
PtrList<fluidThermoModel> thermoFluid(fluidRegions.size());

// Solvers
PtrList<coupledSolver> cpsFluid(fluidRegions.size());
PtrList<sparseSolver<scalar> > spSolverp_rghFluid(fluidRegions.size());
PtrList<sparseSolver<vector> > spSolverUFluid(fluidRegions.size());
PtrList<pressureReference> pRefFluid(fluidRegions.size());

// Options
List<bool> solveCoupledFluid(fluidRegions.size());
List<bool> sPSFluid(fluidRegions.size());   
List<bool> frozenFlowList(fluidRegions.size());
List<bool> simplecFluidList(fluidRegions.size()); 
List<int> nInIterFluid(fluidRegions.size()); 

// Extras
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<Foam::fvModels> fluidFvModels(fluidRegions.size());
PtrList<fvConstraints> fluidFvConstraints(fluidRegions.size());
PtrList<ppUtilInterface> postProcFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    
    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(UFluid[i])  
        )
    );
    
    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set
    (
        i,
        fluidThermoModel::New
        (
          word("thermo"+fluidRegions[i].name()),
          fluidRegions[i]
        ).ptr()
    ); 
    
    Info<< "    Adding to p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    
    Info<< "    Adding to pFluid\n" << endl;
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].hasGravity() 
            ?
            (p_rghFluid[i] + thermoFluid[i].rhok()*thermoFluid[i].gh())()
            :
            p_rghFluid[i]
        )
    );
    
    pRefFluid.set
    (
        i,
        new pressureReference(pFluid[i], p_rghFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"))
    );

    if (p_rghFluid[i].needReference())
    {
      pFluid[i] += dimensionedScalar
      (
        "p",
        pFluid[i].dimensions(),
        pRefFluid[i].refValue()
      - getRefCellValue(pFluid[i], pRefFluid[i].refCell())
      );
 
      p_rghFluid[i] = pFluid[i] - thermoFluid[i].rhok()*thermoFluid[i].gh();
    }

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name()); 
        
    // Create matrix solver
    #include "createSolver.H" 
    
    
    Info<< "    Adding to constEqFluid\n" << endl;
    constEqFluid.set
    (
        i,
        new constitutiveModel
        (
          UFluid[i],
          phiFluid[i]
        )
    ); 
    
    Info<< "    Adding to elecMFluid\n" << endl;
    elecMFluid.set
    (
        i,
        new EDFModel 
        (
          phiFluid[i]
        )
    );  
    
    // C
    IOdictionary cttProperties
    (
      IOobject
      (
        "constitutiveProperties",
        runTime.constant(),
        fluidRegions[i],
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
      )
    );  
    
    sPSFluid[i] = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false);
    
       
    Info<< "    Adding to CFluid\n" << endl;
    CFluid.set
    (
        i,
        sPSFluid[i] != true
        ?
        NULL
        :
        new volScalarField
        (
          IOobject
          (
             "C",
             runTime.timeName(),
             fluidRegions[i],
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
          ),
          fluidRegions[i]
        )
    );
    
    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );
    
    Info<< "    Adding fvModelsFluid\n" << endl;
    fluidFvModels.set
    (
        i,
        new fvModels(fluidRegions[i])
    );

    Info<< "    Adding fvConstraintsFluid\n" << endl;
    fluidFvConstraints.set
    (
        i,
        new fvConstraints(fluidRegions[i])
    );

    
    Info<< "    Adding ppUtil to Fluid\n" << endl;
    const dictionary* ppUDict(fluidRegions[i].solutionDict().subDictPtr("PostProcessing"));   
    postProcFluid.set
    (
        i,
        new ppUtilInterface(ppUDict, fluidRegions[i])
    );
    
    // More options    
    frozenFlowList[i] = fluidRegions[i].solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("frozenFlow", false); 
    simplecFluidList[i] = fluidRegions[i].solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("simplec", true); 
    nInIterFluid[i] = fluidRegions[i].solutionDict().subDict("PIMPLE").lookupOrDefault<int>("nInIter", 1); 
}
