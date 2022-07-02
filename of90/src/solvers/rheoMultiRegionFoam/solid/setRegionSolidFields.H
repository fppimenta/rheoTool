solidThermoModel& thermo = thermoSolid[i];

const fvModels& fvModels = fvModelsSolid[i];
fvConstraints& fvConstraints = fvConstraintsSolid[i];

solidNoLoopControl& pimple = pimples.solid(i);

ppUtilInterface& postProc = postProcSolid[i];
 
