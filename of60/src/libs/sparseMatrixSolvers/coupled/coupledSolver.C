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
  const word& name,
  const fvMesh& mesh
)
:
IOList<label>
(
  IOobject
  (
    word(name),  
    mesh.time().timeName(),  
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE,
    true
  ),
  0
),
mesh_(mesh),
nLocalCells(mesh.nCells()),
nGlobalCells(mesh.nCells()),
isSet(false),
isSysSized(false),
resetX(true),
times_(0),
saveSystem_(mesh.solutionDict().subDict("coupledSolvers").subDict(name).lookupOrDefault<Switch>("saveSystem", false)),
name_(name),
prefix_(word(name + "_")),
updatePrecondFreq_(saveSystem_? readInt(mesh.solutionDict().subDict("coupledSolvers").subDict(name).lookup("updatePrecondFrequency")) : 1e4),
updateA_(saveSystem_?readBool(mesh.solutionDict().subDict("coupledSolvers").subDict(name).lookup("updateMatrixCoeffs")):false),
isRobustSumCheck(saveSystem_?mesh.solutionDict().subDict("coupledSolvers").subDict(name).lookupOrDefault<bool>("robustSumCheck", true):false),
sumCheckDone_(false),
initTimeFlag(true),
initTimeIndex(mesh.time().timeIndex()),
autoPrecond(false)
{

reduce(nGlobalCells, sumOp<int>());

// For the first variable, its first component will
// always correspond to block 0, independent of his type 
firstCmpList.append(0);

// Build inter-processor data if needed               
sparseSolverBase::buildSharedDataOnDemand(mesh); 

// Detect auto mode for update of preconditioner
if (saveSystem_ && updateA_ && updatePrecondFreq_ < 1)
 autoPrecond = true;
  
// Init MPI if running in singleton.
if ( this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
  MPI_Init(NULL,NULL);  
  
// Initialize PETSC
if (this->counterPetsc_ == 0)
{
  // Global database
  std::string petscOptFile(mesh.time().path()/"system"/"petscDict"); 
  if (Pstream::parRun())
    petscOptFile = std::string(mesh.time().path()/".."/"system"/"petscDict");
   
  int argc = 1;
  char** argv = new char*[argc+1];
  argv[1] = NULL;
 
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
            
  this->counterPetsc_--;
  
  // Finalize petsc for the last object of this type
  if (this->counterPetsc_==0)
   ierr = PetscFinalize();
  
  // Since we start MPI, we are also responsible to finalize it 
  if (this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
    MPI_Finalize();           
}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
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
 if (mesh().topoChanging())
 {
  sparseSolverBase::setSharedData(mesh());
  nLocalCells = mesh().nCells();
  nGlobalCells = mesh().nCells();
  reduce(nGlobalCells, sumOp<int>());
   
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

 // Size. If the number of rows is not divisible by nProcs,
 // the remainder is added to the last process.
 int nglb = nGlobalCells*firstCmpList.last();
 int nloc = nglb;
 if (Pstream::parRun())
 {
   nloc = nglb/Pstream::nProcs();
   int dif = nglb - nloc*Pstream::nProcs();  
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
 
 if (mesh().topoChanging() || maxInProcBlocks_.size() == 0)
  computeAllocationPetsc(nloc, nglb);
    
 ierr = MatSeqAIJSetPreallocation(A,0,maxInProcBlocks_.begin());
 ierr = MatMPIAIJSetPreallocation(A,0,maxInProcBlocks_.begin(),0,maxOutProcBlocks_.begin());
 // ierr = MatSetUp(A);CHKERRV(ierr);
 
 //- b
 ierr = VecSetOptionsPrefix(b,prefix_.c_str());CHKERRV(ierr);
 ierr = VecSetSizes(b,nloc,nglb);CHKERRV(ierr);
 ierr = VecSetFromOptions(b);CHKERRV(ierr);
  
 //- x is saved even is savesSystem_ = false. Hence, create at the begining only or
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
 std::string timeOut = std::to_string(mesh().time().timeIndex());
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
  // Do some checks and updates if the mesh topology changes
  if (mesh().topoChanging() && saveSystem_)
  {
    // The system cannot be saved for a mesh changing its topo. Abort.
    if (Pstream::master)
    {
      FatalErrorIn
      (
        "Foam::petscSolver<Type>::solve"
      )   << nl << nl << "coupled matrix " << name_ << " is set to be saved "
          << "over time, but the mesh topology is changing. "
          << exit(FatalError);
    }
  }
  
  if (!isSet)
  {
    FatalErrorIn
    (
       "Foam::coupledSolver::solve"
    )   << nl << nl << "Cannot solve for an empty system. "
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
 if ( (mesh().time().timeIndex() - initTimeIndex) == nEvalInit_)
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
 scalarList initResidual(firstCmpList.last(), 0.);
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
 scalarList finalResidual(firstCmpList.last(), 0.);
 getResiduals(A,b,x,finalResidual);
   
 Info << "Coupled system: " << word(varNames.substr(0, varNames.size()-2)) << endl;
 forAll(varNamesCmpList, ni)
 { 
  Info << "Petsc:" << ksptype << ":" << pctype << ": Solving coupled for "<< varNamesCmpList[ni]  << ", Initial residual = " 
       << initResidual[ni] << ", Final residual = " << finalResidual[ni] << ", No Iterations " << its << endl;   
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
 ierr = VecCopy(b,AxMb);CHKERRV(ierr); // AxMb is b
 
 //- Numerator: |Ax - b|
 ierr = MatMult(A,x,Ax); CHKERRV(ierr);  
 ierr = VecAXPY(AxMb, -1., Ax); CHKERRV(ierr);  
 
 scalarList numNormW(firstCmpList.last(), 0.);
 scalarList den1W(firstCmpList.last(), 0.);
 scalarList den2W(firstCmpList.last(), 0.);
  
 int ilower = this->sharedData.ilower;
 scalar val(1.); 
 scalar xAvW(0.);
 forAll(numNormW, i)
 { 
   for (int j = 0; j<nLocalCells; j++)
   {
     int row = j + ilower + i*nGlobalCells; 
     ierr = VecSetValues(xpart,1,&row,&val,INSERT_VALUES);CHKERRV(ierr);
   } 
   
   ierr = VecAssemblyBegin(xpart);CHKERRV(ierr);
   ierr = VecAssemblyEnd(xpart);CHKERRV(ierr);
   
   // |Ax-b|
   ierr = VecPointwiseMult(tmpV, xpart, AxMb);CHKERRV(ierr);
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
        ierr = MatDuplicate(A,MAT_COPY_VALUES,&A0);CHKERRV(ierr); 
      }
      // Perform mat subtraction and decide
      else
      {        
        Vec xc;
        PetscScalar norm;
        ierr = VecDuplicate(b,&xc);CHKERRV(ierr);
        ierr = MatAXPY(A0,-1,A,SAME_NONZERO_PATTERN);CHKERRV(ierr);
        ierr = MatGetRowMaxAbs(A0,xc,NULL);
        ierr = VecNorm(xc, NORM_1, &norm); CHKERRV(ierr);
        VecDestroy(&xc);
        MatDestroy(&A0); 
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
   for (int j = 0; j<nLocalCells; j++)
   {
     int row = j + ilower + i*nGlobalCells; 
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

// Note: the prealloc we do is clearly oversized, since we assume
// that each block/field communicates with all the other ones. Hence, it 
// is like if we had nb x nb fvMatrix<>. In practice, each field only 
// couples directly to few others. In parallel, it is even worse, because
// we use the total number of expected rows in the FINAL matrix for both
// diagonal and off-diagonal parts. That is, we replace each part size by
// their sum. Needs improvement.
void Foam::coupledSolver::computeAllocationPetsc
(
 int nloc,
 int nglb
)
{ 
 // Resize arrays to 0 because setSize() might be conservative
 if (mesh().topoChanging())
 {
   maxInProcBlocks_.clear();
   maxOutProcBlocks_.clear();
 }
 
 labelList maxInProcFaces(mesh().nCells(), 1); // (n,n) diagonal elements should always exist 
 labelList maxOutProcFaces(mesh().nCells(), 0);
 
 const labelUList& owner = mesh().owner();
 const labelUList& neig = mesh().neighbour();

 //- Diagonal of MPI matrix
 forAll(owner, facei)
 {
   maxInProcFaces[owner[facei]] += 1;
   maxInProcFaces[neig[facei]] += 1;
 }
 
 const fvBoundaryMesh& bMesh(mesh().boundary());
 
 //- Off-diagonal of MPI matrix
 if (Pstream::parRun())
 {
  forAll(bMesh, patchI)
  { 
   if (bMesh[patchI].coupled())
   {
     const labelUList& addr = mesh().lduAddr().patchAddr(patchI); 
     
     if (isType<cyclicFvPatch>(bMesh[patchI]))
     {
       forAll(addr, faceI)
       {
         maxInProcFaces[addr[faceI]] += 1;
       }
     }
     else
     {
       forAll(addr, faceI)
       {
         maxOutProcFaces[addr[faceI]] += 1;
       }
     }
   }
  }
 }
 else
 {
   forAll(bMesh, patchI)
   { 
     if (bMesh[patchI].coupled())
     {
       const labelUList& addr = mesh().lduAddr().patchAddr(patchI); 
    
       if (isType<cyclicFvPatch>(bMesh[patchI]))
       {
         forAll(addr, faceI)
         {
           maxInProcFaces[addr[faceI]] += 1;
         }
       }
     }    
   }
 }
  
 int nBlocks = firstCmpList.last();
 
 if (!Pstream::parRun())
 { 
   forAll(maxInProcFaces, i)
    maxInProcFaces[i] *= nBlocks;
  
   // We assume that all fields contribute to all the others (oversized)
   for (int i = 0; i<nBlocks; i++)
    maxInProcBlocks_.append(maxInProcFaces);
 }
 else
 {
   forAll(maxInProcFaces, i)
     maxInProcFaces[i] += maxOutProcFaces[i]; 
 
   // Now maxInProcFaces has all entries expected in a row for a single field,
   // as if the case was run in serial 
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

   forAll(allMeshFacesIn, i)
     allMeshFacesIn[i] *= nBlocks;
  
   // We assume that all fields contribute to all the others (oversized)
   labelList listTmp;
   for (int i = 0; i<nBlocks; i++)
     listTmp.append(allMeshFacesIn);
  
   maxInProcBlocks_.setSize(nloc, 0);
   maxOutProcBlocks_.setSize(nloc, 0);
   
   int myProcNo = Pstream::myProcNo();
   int nlocRound = nGlobalCells*firstCmpList.last()/Pstream::nProcs();
   forAll(maxInProcBlocks_, i)
   {
     int p = myProcNo*nlocRound + i;
     int total = listTmp[p];
     maxInProcBlocks_[i] = total;  
     maxOutProcBlocks_[i] = total;    
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
