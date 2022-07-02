    const fvMesh& mesh = fluidRegions[i];
    
    constitutiveModel& constEq = constEqFluid[i];
    EDFModel& elecM = elecMFluid[i];
    fluidThermoModel& thermo = thermoFluid[i];
     
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& p_rgh = p_rghFluid[i];
    volScalarField& p = pFluid[i];
      
    ppUtilInterface& postProc = postProcFluid[i];
 
    IOMRFZoneList& MRF = MRFfluid[i];
    fvModels& fvModels = fluidFvModels[i];
    fvConstraints& fvConstraints = fluidFvConstraints[i];
    pressureReference& pressureReference = pRefFluid[i]; 
 
    pimpleNoLoopControl& pimple = pimples.pimple(i);
 
    scalar cumulativeContErr = cumulativeContErrs[i];
    
    // Options
    bool frozenFlow = frozenFlowList[i];
    bool sPS = sPSFluid[i];
    bool solveCoupled = solveCoupledFluid[i];
    bool simplec = simplecFluidList[i];
    int nInIter = nInIterFluid[i];
    
