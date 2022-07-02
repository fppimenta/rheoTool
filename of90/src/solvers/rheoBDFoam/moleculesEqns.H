int nSubCycles = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nSubCycles", 1);

Info<< "Moving molecules." << endl;

scalar t0(runTime.elapsedCpuTime());

bool cont;
for (subCycleTime molcSubCycle(runTime, nSubCycles); !(++molcSubCycle).end();)
{ 
  cont = molecules.update(); 
  // Exit subcycle  
  if (!cont)
   break;
}

execTimeLagrang += (runTime.elapsedCpuTime() - t0);
Info<< "ExecutionTime Lagrangian = " <<  execTimeLagrang << " s" << endl;

// Exit simple loop 
if (!cont)
{
   Info << nl << "No more valid molecules inside the computational domain." << nl << "Exiting time loop." << nl << endl;
   break;
}
   

