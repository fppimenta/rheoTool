// Initialise solid field pointer lists

//- Models
PtrList<solidThermoModel> thermoSolid(solidRegions.size());

//- Options
PtrList<fv::options> solidHeatSources(solidRegions.size());

//- Extras
PtrList<ppUtilInterface> postProcSolid(solidRegions.size());

 
// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding fvOptions\n" << endl;
    solidHeatSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );
    
    Info<< "    Adding ppUtil to Fluid\n" << endl;
    const dictionary* ppUDict(solidRegions[i].solutionDict().subDictPtr("PostProcessing"));   
    postProcSolid.set
    (
        i,
        new ppUtilInterface(ppUDict, solidRegions[i])
    );
    
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
    
}
