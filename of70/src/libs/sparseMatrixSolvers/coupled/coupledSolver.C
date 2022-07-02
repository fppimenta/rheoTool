/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM. 

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coupledSolver.H"
#include <chrono>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(coupledSolver, 0);
  
  // This is the time index until which we force the preconditioner to be updated
  // when saveSystem is enabled. Should be at least 3 because of the matSum check
  // in getResiduals() (fisrt branch of the if will not execute if nEvalInit_ < 3)  
  int Foam::coupledSolver::nEvalInit_ = 3;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledSolver::coupledSolver
(
  const word& solverName,
  const word& meshName,
  const Time& time,
  const dictionary& dict
)
:
IOList<label>
(
  IOobject
  (
    word(solverName + "." + meshName),  
    time.timeName(),  
    time,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
  ),
  0
),
isSet(false),
isSysSized(false),
resetX(true),
times_(0),
saveSystem_(dict.subDict("coupledSolvers").subDict(solverName).lookupOrDefault<Switch>("saveSystem", false)),
name_(solverName),
prefix_(word(solverName + "_")),
updatePrecondFreq_(saveSystem_? readInt(dict.subDict("coupledSolvers").subDict(solverName).lookup("updatePrecondFrequency")) : 1e4),
updateA_(saveSystem_?readBool(dict.subDict("coupledSolvers").subDict(solverName).lookup("updateMatrixCoeffs")):false),
isRobustSumCheck(saveSystem_?dict.subDict("coupledSolvers").subDict(solverName).lookupOrDefault<bool>("robustSumCheck", true):false),
sumCheckDone_(false),
initTimeFlag(true),
initTimeIndex(time.timeIndex()),
autoPrecond(false),
isThereCyclicAMI_(false)
{

// Detect auto mode for update of preconditioner
if (saveSystem_ && updateA_ && updatePrecondFreq_ < 1)
 autoPrecond = true;
  
// Init MPI if running in singleton.
if ( this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
  MPI_Init(NULL,NULL);  
 
// Change prefix for multi-region cases
label nRegions(time.db().names("polyMesh").size());
if (nRegions>1)    
  prefix_ = word(meshName + "." + solverName + "_");
     
// Initialize PETSC
if (this->counterPetsc_ == 0)
{
  // Global database
  std::string petscOptFile
  (
    Pstream::parRun()
    ?
    time.path()/".."/"system"/"petscDict"   
    :
    time.path()/"system"/"petscDict"
  );  
     
  int argc = 0;
  char** argv = new char*[1];
  argv[0] = NULL;
 
  PetscInitialize(&argc,&argv,petscOptFile.c_str(),NULL);
    
  delete [] argv; 
}

this->counterPetsc_++; 
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::coupledSolver::~coupledSolver()
{      
  // Note: if resetX = false then this is an indicator
  // that createSystem() has been called at least once
  // and so there are things to clean up. This is to 
  // prevent Destroy() on uncreated object, which could 
  // happen for an unused coupledSolver object.      
  
  // Destroy
  if (saveSystem_ && !resetX)
  {     
    KSPDestroy(&ksp);
    VecDestroy(&b); VecDestroy(&x);  
    MatDestroy(&A);           
  }
  else if (!saveSystem_ && !resetX)
  {
    VecDestroy(&x); 
  }
  
  // It may happen that the number of solve() calls is less 
  // than nEvalInit_, so in this case the robust sum check was not done,
  // hence A0 is still alive and we need to destoy it. Only applies if
  // robustsumcheck flag is on, otherwise the autoPtr is always NULL and
  // nothing is done
  if (!A0.empty())
   MatDestroy(&A0()); 
            
  this->counterPetsc_--;
  
  // Finalize petsc for the last object of this type
  if (this->counterPetsc_==0)
   ierr = PetscFinalize();
  
  // Since we start MPI, we are also responsible to finalize it 
  if (this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
    MPI_Finalize();           
}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledSolver::insertMesh
(
  const fvMesh& mesh
)
{ 
  if (meshList.size() == 0)
  {  
    meshList.append(meshListS());
    meshList.last().mesh = &mesh;
        
    // Build inter-processor data if needed               
    meshList.last().ID = sparseSolverBase::buildSharedDataOnDemand(mesh); 
    
    // Verify limitations of the interface
    checkLimitations(mesh);   
  }
  else
  {
    bool exists(false);
    forAll(meshList, i)
    {
      if ( meshList[i].mesh->name() == mesh.name())
      {
        exists = true;
      }
    } 
    
    if (!exists) 
    {
      meshList.append(meshListS());
      meshList.last().mesh = &mesh;
        
      // Build inter-processor data if needed               
      meshList.last().ID = sparseSolverBase::buildSharedDataOnDemand(mesh); 
      
      // Verify limitations of the interface
      checkLimitations(mesh);
    }   
  }
}

void Foam::coupledSolver::checkLimitations
(
  const fvMesh& mesh
)
{
  
 forAll(mesh.boundary(), patchI)
 { 
   const fvPatch& pfvPatch = mesh.boundary()[patchI];
 
   if (pfvPatch.coupled() && pfvPatch.size() > 0)
   {     
     if (isType<cyclicFvPatch>(pfvPatch))
     {  
       // No support for rotational transform in non-scalar variables  
       if (
            refCast<const cyclicFvPatch>(pfvPatch).cyclicPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           FatalErrorInFunction
           << "No support in Petsc interface for patches of type cyclic rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << exit(FatalError);
       }
     }
     else if (isType<processorCyclicFvPatch>(pfvPatch))
     {
       // No support for rotational transform in non-scalar variables     
       if (
            refCast<const processorCyclicFvPatch>(pfvPatch).procPolyPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           FatalErrorInFunction
           << "No support in Petsc interface for patches of type cyclic rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << exit(FatalError);
       }      
     }
     else if (isType<cyclicAMIFvPatch>(pfvPatch))
     {
       // No support for rotational transform in non-scalar variables     
       if (
            refCast<const cyclicAMIFvPatch>(pfvPatch).cyclicAMIPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           FatalErrorInFunction
           << "No support in Petsc interface for patches of type cyclicAMI rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << exit(FatalError);
       } 
       
       isThereCyclicAMI_ = true;      
     }
     else if (isType<processorFvPatch>(pfvPatch))
     {
       // Full support
     }
     else
     {
       FatalErrorInFunction
       << "Sorry, but Petsc interface cannot support patches of type: "
       << pfvPatch.type() << ", which is currently set for "
       << "patch " << pfvPatch.name() << "."
       << exit(FatalError);  
     }     
   }
 }


}


void Foam::coupledSolver::resetA()
{
  ierr = MatZeroEntries(A); CHKERRV(ierr);
}

void Foam::coupledSolver::resetb()
{
  ierr = VecSet(b,0.); CHKERRV(ierr);
  // Cannot reset x because we need it to compute the initial residuals 
}

void Foam::coupledSolver::createSystem()
{
 if (isSysSized)
  return;
 
 // If topo changes, update number of cells and shared data, destroy
 // x and create it again with correct dimensions.  
 bool topoChanging(false);
 bool changing(false);
 forAll(meshList, i)
 {
  topoChanging = (topoChanging || meshList[i].mesh->topoChanging()); 
  changing = (changing || meshList[i].mesh->changing()); 
 }
 
 if (topoChanging)
 {
  forAll(meshList, i)
    sparseSolverBase::setSharedData(*meshList[i].mesh);
    
  // Since we destroy x , the initial residual
  // displayed each time the mesh is updated will be 1. Ignore it. 
  resetX = true;
  
  // Detroy x every time it enters here, except upon the first call,
  // since x would still not exist.
  if (times_>0)
    VecDestroy(&x);
 }
 
 // Create
 ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRV(ierr);
 ierr = VecCreate(PETSC_COMM_WORLD, &b);CHKERRV(ierr);
 ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRV(ierr);

 // Update mesh-dependent members of varInfo
 forAll(varInfo, i)
 {
   label meshID = varInfo[i].meshID;
   label ncI = meshList[meshID].mesh->nCells();
   reduce(ncI, sumOp<label>());
   varInfo[i].nCells = ncI;
 }
 
 forAll(varInfo, i)
 {
   // First field always starts at 0 (no need to change)
   if (i != 0)
   {
     varInfo[i].firstElem = varInfo[i-1].firstElem + varInfo[i-1].nValidCmp*varInfo[i-1].nCells;
   }
 }

 // Size. If the number of rows is not divisible by nProcs,
 // the remainder is added to the last process.
 
 label nglb(0);
 forAll(meshList, i)
 {
   label nMeshI = meshList[i].mesh->nCells();
   reduce(nMeshI, sumOp<label>());
   nglb += nMeshI*meshList[i].nValidCmp;
 }
 label nloc(nglb);
  
 if (Pstream::parRun())
 {
   nloc = nglb/Pstream::nProcs();
   label dif = nglb - nloc*Pstream::nProcs();  
   if (dif > 0)
   {
     if (Pstream::myProcNo() == Pstream::nProcs()-1)
       nloc += dif;
   }   
 }

 //- A
 ierr = MatSetOptionsPrefix(A,prefix_.c_str());CHKERRV(ierr);
 ierr = MatSetSizes(A,nloc,nloc,nglb,nglb);CHKERRV(ierr);
 ierr = MatSetFromOptions(A);CHKERRV(ierr);

 if (topoChanging || maxInProcBlocks_.size() == 0 || (isThereCyclicAMI_ && changing))
  computeAllocationPetsc(nloc, nglb);
 
 ierr = MatSeqAIJSetPreallocation(A,0,maxInProcBlocks_.begin());
 ierr = MatMPIAIJSetPreallocation(A,0,maxInProcBlocks_.begin(),0,maxOutProcBlocks_.begin());

// ierr = MatSetUp(A);CHKERRV(ierr);
 
 //- b
 ierr = VecSetOptionsPrefix(b,prefix_.c_str());CHKERRV(ierr);
 ierr = VecSetSizes(b,nloc,nglb);CHKERRV(ierr);
 ierr = VecSetFromOptions(b);CHKERRV(ierr);
 
 //- x is saved even if savesSystem_ = false. Hence, create at the begining only or
 // create/destroy always if mesh topo changes.
 if (resetX)
 {
  ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRV(ierr); 
  ierr = VecSetOptionsPrefix(x,prefix_.c_str());CHKERRV(ierr);
  ierr = VecSetSizes(x,nloc,nglb);CHKERRV(ierr);
  ierr = VecSetFromOptions(x);CHKERRV(ierr);
  ierr = VecSet(x,0.); CHKERRV(ierr);
 }
 
 //- Set flag 
 isSysSized = true; 
 resetX = false;
}

void Foam::coupledSolver::printSystem
(
  Mat& A,
  Vec& b
)
{
 std::string timeOut = std::to_string(meshList[0].mesh->time().timeIndex());
 std::string Alab = timeOut + "_A.txt";
 std::string blab = timeOut + "_b.txt";
  
 // Matrix of coeffs
 {
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, Alab.c_str(), &viewer);
   PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT); //PETSC_VIEWER_ASCII_DENSE   
   MatView(A,viewer);
   PetscViewerPopFormat(viewer);
   PetscViewerDestroy(&viewer);
 }
    
 // RHS vector
 {    
   PetscViewer viewer;
   PetscViewerASCIIOpen(PETSC_COMM_WORLD, blab.c_str(), &viewer);
   PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT); //PETSC_VIEWER_ASCII_DENSE 
   VecView(b,viewer);
   PetscViewerPopFormat(viewer);
   PetscViewerDestroy(&viewer);
 }      
}


void Foam::coupledSolver::solve
()
{
  bool topoChanging(false);
  bool changing(false);
  forAll(meshList, i)
  {
   topoChanging = (topoChanging || meshList[i].mesh->topoChanging()); 
   changing = (changing || meshList[i].mesh->changing()); 
  }

  // Do some checks and updates if the mesh topology changes
  if (topoChanging && saveSystem_)
  {
    // The system cannot be saved for a mesh changing its topo. Abort.
    if (Pstream::master)
    {
      FatalErrorInFunction
      << "Coupled matrix " << name_ << " is set to be saved "
      << "over time, but the mesh is changing and has AMIs. Option saveSystem should be set as false."
      << exit(FatalError); 
    }
  }
  
  // Do some checks and updates if the mesh topology changes
  if (changing && saveSystem_ && isThereCyclicAMI_)
  {
    // The system cannot be saved for a mesh changing and having AMIs. Abort.
    if (Pstream::master)
    {
      FatalErrorInFunction
      << "Coupled matrix " << name_ << " is set to be saved "
      << "over time, but the mesh is changing and has AMIs. Option saveSystem should be set as false."
      << exit(FatalError);
    }
  }
  
  if (!isSet)
  {
    FatalErrorInFunction
    << "Cannot solve for an empty system. "
    << exit(FatalError);
  }
 
  solvePetsc();
  
  times_++;
}

void Foam::coupledSolver::solvePetsc()
{
 
 // Assemble matrix/vector
 if (!saveSystem_ || updateA_ || initTimeFlag) 
 {
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
 }
  
 ierr = VecAssemblyBegin(b);CHKERRV(ierr);
 ierr = VecAssemblyEnd(b);CHKERRV(ierr);
  
 if (this->debug_)
  printSystem(A, b);
    
 ierr = KSPSetOptionsPrefix(ksp,prefix_.c_str());CHKERRV(ierr);
 
 // Set operators   
 if (initTimeFlag || times_ > updatePrecondFreq_ || !saveSystem_)
 {
   ierr = KSPSetReusePreconditioner(ksp, PETSC_FALSE); CHKERRV(ierr);
   ierr = KSPSetOperators(ksp,A,A); CHKERRV(ierr);
 
   // Reset counter
   times_ = 1;    
 }
 else
 {
   ierr = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERRV(ierr);   
 }
  
 // Once we reach time index nEvalInit_, for the first time in an inner 
 // iteration loop, the system should be ready and is no longer updated 
 // in case saveSystem_ is true and updateA is false.
 if ( (meshList[0].mesh->time().timeIndex() - initTimeIndex) == nEvalInit_)
  initTimeFlag = false;
    
 ierr = KSPSetFromOptions(ksp);CHKERRV(ierr);
 
 KSPType  ksptype; 
 ierr = KSPGetType(ksp, &ksptype);CHKERRV(ierr);
  
 // Prevent users from using direct solvers when A is changing and ksp is
 // frozen 
 if (strcmp(ksptype,"preonly") == 0 && saveSystem_ && updateA_)
 {
   if (updatePrecondFreq_ != 1)
   {
     if (Pstream::master)
     {
      FatalErrorIn
      (
        "Foam::coupledSolver::solvePetsc"
      )   << nl << nl << "You cannot freeze the sparse matrix solver "
      << name_ << " while the matrix is being"
      << " updated, for a direct solver (preonly). Either set 'updatePrecondFrequency' to 1 or"
      << " choose an iterative matrix solver." << nl << nl << endl
      << exit(FatalError);
     }
   }
 } 
  
 // Get initial residuals
 int nBlocks = varInfo.last().firstCmp + varInfo.last().nValidCmp;
 scalarList initResidual(nBlocks, 0.);
  
 getResiduals(A,b,x,initResidual);

 auto start = std::chrono::high_resolution_clock::now();
  
 // Solve
 ierr = KSPSolve(ksp,b,x);CHKERRV(ierr);
 
 // Adjust precond update frequency if auto mode.
 // If parallel run get the average (could be maxOp either
 // to be more conservative) cpu time across processors
 // because we need the same time on all processors to 
 // ensure same updatePrecondFreq_ across all of them and so
 // a synchronized entry in precond update. 
 auto elapsed = std::chrono::high_resolution_clock::now() - start;
 scalar solveTime = scalar(std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count());
 
 if (Pstream::parRun())
 {
   reduce(solveTime, sumOp<double>());
   solveTime /= Pstream::nProcs();
 }
   
 if (autoPrecond)
 {  
   if (times_ == 2)
     solveTime2 = solveTime;  
   
   if (times_ == 1)
   {
     solveTime1 = solveTime;
     updatePrecondFreq_ = 1e5;
   }
   else
   {
     updatePrecondFreq_ = min(1e5, max(1, int( (solveTime1/solveTime) * (1./mag(1.-solveTime2/solveTime + 1e-6)) )));
   }    
 }
 
 // Print info
 PetscInt its;
 PCType   pctype;
 PC pc;
   
 ierr = KSPGetIterationNumber(ksp,&its);CHKERRV(ierr);
 ierr = KSPGetPC(ksp, &pc);
 ierr = PCGetType(pc, &pctype);CHKERRV(ierr);
  
 // Get final residuals
 scalarList finalResidual(nBlocks, 0.);
 getResiduals(A,b,x,finalResidual);
   
 Info << "Coupled system: " << word(varNames.substr(0, varNames.size()-2)) << endl;
 if (meshList.size() == 1)
 {
  forAll(varNamesCmpList, ni)
  { 
    Info << "Petsc:" << ksptype << ":" << pctype << ": Solving coupled for "<< varNamesCmpList[ni]  << ", Initial residual = " 
         << initResidual[ni] << ", Final residual = " << finalResidual[ni] << ", No Iterations " << its << endl;   
  }
 }
 else
 {
  // If there are fields from different regions in the same coupled system, then prepend region's name
  // to var's name to distinguish 
  List<word> meshName(nBlocks," ");
  int ii(0);
  forAll(varInfo, i)
  {
    for (int j=0; j<varInfo[i].nValidCmp; j++)
    {
     meshName[ii] = meshList[varInfo[i].meshID].mesh->name();
     ii++;
    }
  }
  
  forAll(varNamesCmpList, ni)
  { 
    Info << "Petsc:" << ksptype << ":" << pctype << ": Solving coupled for "<< meshName[ni] << "." 
         << varNamesCmpList[ni]  << ", Initial residual = " << initResidual[ni] << ", Final residual = " 
         << finalResidual[ni] << ", No Iterations " << its << endl;   
  }
 }
  
 // Collect and transfer solution for all field types.
 // If there are no fields for a particular type, it will do nothing.
 getSolution<scalar>();
 getSolution<vector>();
 getSolution<symmTensor>();
 getSolution<tensor>();
 getSolution<sphericalTensor>();
 
 if (saveSystem_)
 { 
   // Since we always ADD_VALUES in matrix/vector assembly, we 
   // need to reset the coefficients before a new assembly
   // routine takes place (except if !updateA_ for A) 
   resetb();
   
   if (initTimeFlag || updateA_)
     resetA();
 }
 else
 {
   VecDestroy(&b); 
   MatDestroy(&A); 
   KSPDestroy(&ksp);
   isSysSized = false;
 }
 
}

// The residuals are computed using openfoam's definition, which
// is applied to each individual sub-equation of the coupled system (rectangular matrix; block row).
// Filters are used to compute xmean and the norms individually for each sub-equation
// (component to be solved). Building the filter (xpart) can be optimized. 
void Foam::coupledSolver::getResiduals 
( 
  const Mat& A,
  const Vec& b,
  const Vec& x, 
  scalarList& residual
)
{
 // Number of valid cmps inside the matrix/vectors
 label nCmp(0);
 forAll(varInfo, i)
  nCmp += varInfo[i].nValidCmp;
 
 // Aux list {on/off, firstElem, nCells, meshID} for each solvable component
 List<List<label>> aux(nCmp, {0,-1,-1,-1});
 forAll(varInfo, i)
 {
   aux[varInfo[i].firstCmp][0] = 1; // 1 means that the given component is first one of the field
   aux[varInfo[i].firstCmp][1] = varInfo[i].firstElem; // first element position
   aux[varInfo[i].firstCmp][2] = varInfo[i].nCells;  // ncells
   aux[varInfo[i].firstCmp][3] = varInfo[i].meshID;  // meshID
 }
 
 // Fill-in aux entries corresponding to all non-first components (aux[*][0]=0)
 // Note that the components of a multi-component var are contiguous and that is 
 // why we can use (i-1) below.
 forAll(aux, i)
 {
   if (aux[i][0] == 0)
   {     
     aux[i][1] = aux[i-1][1]+aux[i-1][2]; // first element is set as the previous first element + previous ncells
     aux[i][2] = aux[i-1][2]; // n cells is the same 
     aux[i][3] = aux[i-1][3]; // mesh ID is the same
     aux[i][0] = 1;           // flag it on 
   }
 }
 
 /*
 forAll(aux, i)
 {
    Info << nl << "AUX" << tab << i << nl << endl;
    Info << aux[i][0] << nl
         << aux[i][1] << nl
         << aux[i][2] << nl
         << aux[i][3] << nl
         << endl;
 }
 */
 
 // Note: VecNorm() automatically scatters->gathers the norm.
 // Output of Info and Pout of norm is the same.  
 
 Vec Ax, xpart, tmpV, AxMb, xCum;
   
 // Create & size 
 ierr = VecDuplicate(b,&Ax);CHKERRV(ierr);
 ierr = VecDuplicate(x,&xpart);CHKERRV(ierr);
 ierr = VecDuplicate(x,&tmpV);CHKERRV(ierr); 
 ierr = VecDuplicate(b,&AxMb);CHKERRV(ierr);
 ierr = VecDuplicate(x,&xCum);CHKERRV(ierr);
 
 ierr = VecSet(xpart,0.); CHKERRV(ierr);
 ierr = VecSet(xCum,0.); CHKERRV(ierr);
 ierr = VecCopy(b,AxMb); CHKERRV(ierr); // AxMb is b
 
 //- Numerator: |Ax - b|
 ierr = MatMult(A,x,Ax); CHKERRV(ierr);  
 ierr = VecAXPY(AxMb, -1., Ax); CHKERRV(ierr);  
 
 scalarList numNormW(nCmp, 0.);
 scalarList den1W(nCmp, 0.);
 scalarList den2W(nCmp, 0.);
 
 scalar val(1.); 
 scalar xAvW(0.);
 forAll(numNormW, i)
 { 
   label ilower = this->sharedData[meshList[aux[i][3]].ID].ilower;
   label nLocal = meshList[aux[i][3]].mesh->nCells();
   label nGlobalCells = aux[i][2]; // global number of cells (equivalent to reduceOp(nLocal))
 
   for (label j=0; j<nLocal; j++)
   {
     label row = j + ilower + aux[i][1]; 
     ierr = VecSetValues(xpart,1,&row,&val,INSERT_VALUES);CHKERRV(ierr);
   } 
   
   ierr = VecAssemblyBegin(xpart);CHKERRV(ierr);
   ierr = VecAssemblyEnd(xpart);CHKERRV(ierr);
   
   // |Ax-b|
   ierr = VecPointwiseMult(tmpV, xpart, AxMb); CHKERRV(ierr);
   ierr = VecNorm(tmpV, NORM_1, &numNormW[i]); CHKERRV(ierr);  
   
   // xmean
   ierr = VecPointwiseMult(tmpV, xpart, x);CHKERRV(ierr);
   ierr = VecSum(tmpV, &xAvW); CHKERRV(ierr); 
   xAvW /= nGlobalCells;  
   
   // xCum +=
   ierr = VecCopy(xpart, tmpV);CHKERRV(ierr);  
   ierr = VecScale(tmpV, xAvW); CHKERRV(ierr); 
   ierr = VecAXPY(xCum, 1., tmpV); CHKERRV(ierr); 
     
   ierr = VecSet(xpart, 0.); CHKERRV(ierr);
 }

 //- Denominator (norm scale): |A*xav - b| + |Ax - A*xav| + SMALL
 
 // A*xav
 ierr = MatGetRowSum(A, AxMb); CHKERRV(ierr);  
 ierr = VecPointwiseMult(tmpV, AxMb, xCum);CHKERRV(ierr); 
 
 
 
   // Here we check if the coefficients matrix A is changing over time.
   // If this is the case and updateA_ is disabled, then retrieve error.
   // The code is here to re-use AxMb. 
   if (!updateA_ && saveSystem_ && !sumCheckDone_)
   {
    // Robust check through matrix subtraction
    if (isRobustSumCheck)
     {
      // Udpdate A0
      if (initTimeFlag)   
      {  
        if (A0.empty())
        {
          A0.set(new Mat);    
          ierr = MatDuplicate(A,MAT_COPY_VALUES,&A0());CHKERRV(ierr);         
        }
        else
        {
          ierr = MatCopy(A,A0(),SAME_NONZERO_PATTERN);CHKERRV(ierr);          
        }          
      }
      // Perform mat subtraction and decide
      else
      {        
        Vec xc;
        PetscScalar norm;
        ierr = VecDuplicate(b,&xc);CHKERRV(ierr);
        ierr = MatAXPY(A0(),-1,A,SAME_NONZERO_PATTERN);CHKERRV(ierr);
        ierr = MatGetRowMaxAbs(A0(),xc,NULL); CHKERRV(ierr);
        ierr = VecNorm(xc, NORM_1, &norm); CHKERRV(ierr);
        VecDestroy(&xc);
        MatDestroy(&A0()); A0.clear(); // Make autoPtr NULL as flag in dtor
        if (mag(norm) > 1e-20)
        {
          if (Pstream::master)
          {
            FatalErrorIn
            (
            "Foam::coupledSolver::getResiduals"
            )   << nl << nl << "Norm of coupled matrix " << name_ << " is changing by " << norm
            << " between consecutive time-steps and option updateMatrixCoeffs is disabled."
            << " You should enable this option or set option saveSystem to false." << nl << nl << endl
            << exit(FatalError);
          }
        }
        
        sumCheckDone_ = true;
      }
     }
    // Weak check through matrix coeffs sum
     {
      // Update Asum
      if (initTimeFlag)
      {
        ierr = VecNorm(AxMb, NORM_1, &Asum); CHKERRV(ierr);
      }
      // Evaluate norm variation and decide
      else 
      {
        PetscScalar Asumti;
        ierr = VecNorm(AxMb, NORM_1, &Asumti); CHKERRV(ierr);
        PetscScalar dif = mag(Asumti-Asum); 
        if (dif > 1e-20)
        {
          if (Pstream::master)
          {
            FatalErrorIn
            (
            "Foam::coupledSolver::getResiduals"
            )   << nl << nl << "Sum of coupled matrix " << name_ << " is changing by " << dif
            << " between consecutive time-steps and option updateMatrixCoeffs is disabled."
            << " You should enable this option or set option saveSystem to false." << nl << nl << endl
            << exit(FatalError);
          }
        }
        
        sumCheckDone_ = true;
      } 
     }
   } 
 
 
 
 // First term
 ierr = VecCopy(b,AxMb);CHKERRV(ierr);  
 ierr = VecAXPY(AxMb, -1., tmpV); CHKERRV(ierr);
 
 // Second term
 ierr = VecAXPY(tmpV, -1., Ax); CHKERRV(ierr);
 
 // Filter per block equation
 forAll(numNormW, i)
 {    
   label ilower = this->sharedData[meshList[aux[i][3]].ID].ilower;
   label nLocal = meshList[aux[i][3]].mesh->nCells();
 
   for (label j=0; j<nLocal; j++)
   {
     label row = j + ilower + aux[i][1];  
     ierr = VecSetValues(xpart,1,&row,&val,INSERT_VALUES);CHKERRV(ierr);
   }

   ierr = VecAssemblyBegin(xpart);CHKERRV(ierr);
   ierr = VecAssemblyEnd(xpart);CHKERRV(ierr);
   
   // |A*xav - b| = AxMb 
   ierr = VecPointwiseMult(xCum, xpart, AxMb);CHKERRV(ierr);
   ierr = VecNorm(xCum, NORM_1, &den1W[i]); CHKERRV(ierr);  
   
   // |Ax - A*xav|  = tmpV
   ierr = VecPointwiseMult(xCum, xpart, tmpV);CHKERRV(ierr);
   ierr = VecNorm(xCum, NORM_1, &den2W[i]); CHKERRV(ierr);    

   // Reset filter   
   ierr = VecSet(xpart,0.); CHKERRV(ierr);
 }
  
 // Compute residual per block equation
 forAll(numNormW, i)
 {
   scalar n = numNormW[i];
   scalar d1 = den1W[i];
   scalar d2 = den2W[i];
   if ( n + d1 + d2 < 1e-30)                     
     residual[i] = 1;                           
   else                                         
     residual[i] = n/(d1+d2+1e-20);
 }
  
 // Cleanup
 VecDestroy(&Ax); 
 VecDestroy(&xpart); 
 VecDestroy(&tmpV); 
 VecDestroy(&AxMb); 
 VecDestroy(&xCum); 
}


// Algorithm: first we check how many non-zero values per row we have in a
// matrix composed of a single scalar field. We do this for all different
// meshes inserted in meshList. We multiply the value in each row by the 
// number of valid components, as if the scalar-matrix was repeating in
// each block (horizontally). We then make a list by appending vertically
// the list (allocation) expected for that scalar-field matrix. The list appended
// is the one corresponding to mesh of the field in that position.   

// Note 1: the prealloc we do is clearly oversized, since we assume
// that each block/field communicates with all the other ones. Hence, it 
// is like if we had nb x nb fvMatrix<>. In practice, each field only 
// couples directly to few others. In parallel, it is even worse, because
// we use the total number of expected rows in the FINAL matrix for both
// diagonal and off-diagonal parts. That is, we replace each part size by
// their sum. Needs improvement.

void Foam::coupledSolver::computeAllocationPetsc
(
 label nloc,
 label nglb
)
{ 
 
 //-- PART I (as if single field)
 
 // Resize arrays to 0 because setSize() might be conservative
 maxInProcBlocks_.clear();
 maxOutProcBlocks_.clear();
 
 // Each element of this list corresponds to a different mesh in meshList.
 // Each element is the number of non-zero values expected for each row of 
 // a matrix having a single field (size nxn where n is the total number of
 // cells of the given mesh). Both diagonal and non-diagonal elements are 
 // accounted for (summed).
 List<labelList> meshIAlloc;  
 
 int nBlocks = varInfo.last().firstCmp + varInfo.last().nValidCmp;
 
 forAll(meshList, mI)
 { 
   labelList maxInProcFaces(meshList[mI].mesh->nCells(), 1); // (n,n) diagonal elements should always exist 
   labelList maxOutProcFaces(meshList[mI].mesh->nCells(), 0);
 
   const labelUList& owner = meshList[mI].mesh->owner();
   const labelUList& neig = meshList[mI].mesh->neighbour();

   //- Diagonal of MPI matrix
   forAll(owner, facei)
   {
     maxInProcFaces[owner[facei]] += 1;
     maxInProcFaces[neig[facei]] += 1;
   }
 
   const fvBoundaryMesh& bMesh(meshList[mI].mesh->boundary());
 
   //- Off-diagonal of MPI matrix 
   forAll(bMesh, patchI)
   { 
    #include "specialBoundariesAlloc.H"  
     
    if (bMesh[patchI].coupled())
    {
      const labelUList& addr = meshList[mI].mesh->lduAddr().patchAddr(patchI); 
      const fvPatch& pfvPatch = bMesh[patchI];
     
      // cyclic
      if (isType<cyclicFvPatch>(bMesh[patchI]))
      {
        forAll(addr, faceI)
        {
          maxInProcFaces[addr[faceI]] += 1;
        }
      }
      // cyclicAMI
      else if (isType<cyclicAMIFvPatch>(pfvPatch))
      {       
       const cyclicAMIPolyPatch& camipp = refCast<const cyclicAMIFvPatch>(pfvPatch).cyclicAMIPatch();
       const cyclicAMIPolyPatch& neicamipp = camipp.neighbPatch();
       
       const labelList& ownFC = camipp.faceCells();   
     
       if (camipp.owner())
       {      
         forAll(camipp.AMIs(), i)
         {
           const labelListList& srcAd = camipp.AMIs()[i].srcAddress();
 
           // Each patch of the AMI interface in a different processor
           if (camipp.AMIs()[i].singlePatchProc() == -1) 
           {            
             forAll(srcAd, facei)
             {
               forAll(srcAd[facei], kk)
               {
                 if (camipp.AMIs()[i].applyLowWeightCorrection())
                 {            
                   if (camipp.AMIs()[i].srcWeightsSum()[facei] < camipp.AMIs()[i].lowWeightCorrection())
                   {
                     maxInProcFaces[ownFC[facei]] += 1;
                   }
                   else
                   {
                     maxOutProcFaces[ownFC[facei]] += 1;                 
                   }
                 }
                 else
                 {
                   maxOutProcFaces[ownFC[facei]] += 1;                 
                 }
               }            
             }
             
           }
           // Both patches of the AMI interface in same processor
           else
           {
             forAll(srcAd, facei)
             {
               forAll(srcAd[facei], kk)
               {
                 maxInProcFaces[ownFC[facei]] += 1;
               }            
             }
           }
   
        }   
           
       }  
       else
       {      
         forAll(neicamipp.AMIs(), i)
         {
           const labelListList& tgtAd = neicamipp.AMIs()[i].tgtAddress();
      
           // Each patch of the AMI interface in a different processor
           if (neicamipp.AMIs()[i].singlePatchProc() == -1) 
           {             
             forAll(tgtAd, facei)
             {
               forAll(tgtAd[facei], kk)
               {
                 if (neicamipp.AMIs()[i].applyLowWeightCorrection())
                 {           
                   if (neicamipp.AMIs()[i].tgtWeightsSum()[facei] < neicamipp.AMIs()[i].lowWeightCorrection())
                   {
                     maxInProcFaces[ownFC[facei]] += 1;
                   }
                   else
                   {                 
                     maxOutProcFaces[ownFC[facei]] += 1;
                   }
                 }
                 else
                 {                 
                   maxOutProcFaces[ownFC[facei]] += 1;
                 }
               }            
             }
             
           }
           // Both patches of the AMI interface in same processor
           else
           {
             forAll(tgtAd, facei)
             {
               forAll(tgtAd[facei], kk)
               {
                 maxInProcFaces[ownFC[facei]] += 1;
               }           
             }
           }
   
        }   
           
       }   
      }
      // processor and processorCyclic
      else
      {
        forAll(addr, faceI)
        {
          maxOutProcFaces[addr[faceI]] += 1;
        }
      }
    }
   }
  
 
   if (!Pstream::parRun())
   {     
     // Consider other fields in the horizontal  
     forAll(maxInProcFaces, i)
       maxInProcFaces[i] *= nBlocks;
   
     meshIAlloc.append(maxInProcFaces);
   }
   else
   {
     forAll(maxInProcFaces, i)
        maxInProcFaces[i] += maxOutProcFaces[i]; 
 
     // Now maxInProcFaces has all entries expected in a row for a single field,
     // as if the case was run in serial. We will now send this allocation to all
     // processors and, inside each processor, we will combine all the lists to get
     // a single (and the same) list with the alocation equivalent to a serial run.
     List< labelList > allMeshFacesInP(Pstream::nProcs());

     allMeshFacesInP[Pstream::myProcNo()] = maxInProcFaces;
 
     Pstream::gatherList(allMeshFacesInP);
   
     Pstream::scatterList(allMeshFacesInP);

     labelList allMeshFacesIn
     (
       ListListOps::combine<List<label> >
       (
        allMeshFacesInP,
        accessOp<List<label> >()
       )
     ); 

     // Consider other fields in the horizontal
     forAll(allMeshFacesIn, i)
       allMeshFacesIn[i] *= nBlocks;
     
     meshIAlloc.append(allMeshFacesIn);
  
   }
 }
   
 //-- PART II (multi-fields) 
 
 if (!Pstream::parRun())
 { 
   // We assume that all fields contribute to all the others (oversized)
   forAll(varInfo, i)
   { 
     for (int j=0; j<varInfo[i].nValidCmp; j++)
     {
       label mID(varInfo[i].meshID);
       maxInProcBlocks_.append(meshIAlloc[mID]);
     }
   }
 }
 else
 {
   // We assume that all fields contribute to all the others (oversized)
   labelList listTmp;
   forAll(varInfo, i)
   { 
     for (int j=0; j<varInfo[i].nValidCmp; j++)
     {
       label mID(varInfo[i].meshID);
       listTmp.append(meshIAlloc[mID]);
     }
   }

   maxInProcBlocks_.setSize(nloc, 0);
   maxOutProcBlocks_.setSize(nloc, 0);
   
   int myProcNo = Pstream::myProcNo();
   
   int nlocRound = nglb/Pstream::nProcs();   
   
   // All procs have a complet list of the allocation needed for the complet
   // super-matrix (as if serial run). Now each processor just needs to grab
   // the allocation corresponding to the rows it owns.          
   forAll(maxInProcBlocks_, i)
   {
     label p = myProcNo*nlocRound + i; // Global index
     label total = listTmp[p];
     maxInProcBlocks_[i] = min(nloc,total); // min: we cannot allocate more nnz than the local size of diagonal matrix 
     maxOutProcBlocks_[i] = min(nglb-nloc,total); // min: we cannot allocate more nnz than the local size of off-diagonal matrix    
   } 
 }
  
}



// Explicit instantiation 
makeTypeLMatrix(scalar); makeTypefvMatrix(scalar);
makeTypeLMatrix(vector); makeTypefvMatrix(vector);
makeTypeLMatrix(symmTensor); makeTypefvMatrix(symmTensor);
makeTypeLMatrix(tensor); makeTypefvMatrix(tensor);
makeTypeLMatrix(sphericalTensor); makeTypefvMatrix(sphericalTensor);

 
// ************************************************************************* //
