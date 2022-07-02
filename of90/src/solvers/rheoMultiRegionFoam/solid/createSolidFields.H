// Initialise solid field pointer lists

//- Models
PtrList<solidThermoModel> thermoSolid(solidRegions.size());

//- Options
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());

//- Extras
PtrList<ppUtilInterface> postProcSolid(solidRegions.size());

 
// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;
        
        
    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set
    (
        i,
        solidThermoModel::New
        (
          word("thermo"+solidRegions[i].name()),
          solidRegions[i]
        ).ptr()
    ); 
    
    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );
    
    Info<< "    Adding ppUtil to Solid\n" << endl;
    const dictionary* ppUDict(solidRegions[i].solutionDict().subDictPtr("PostProcessing"));   
    postProcSolid.set
    (
        i,
        new ppUtilInterface(ppUDict, solidRegions[i])
    ); 
}
