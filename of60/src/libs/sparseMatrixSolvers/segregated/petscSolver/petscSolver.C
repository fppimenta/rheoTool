/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "petscSolver.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"
#include <chrono>

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

// This is the time index until which we force the preconditioner to be updated
// when saveSystem is enabled.  
template<class Type>
int Foam::petscSolver<Type>::nEvalInit_ = 3;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::petscSolver<Type>::petscSolver
(
    const GeometricField<Type, fvPatchField, volMesh>& T,
    const fvMesh& mesh,
    const dictionary& fvSolution
)
:
sparseSolver<Type>(T, mesh, fvSolution),
saveSystem_(fvSolution.subDict("solvers").subDict(word(T.name())).lookupOrDefault<Switch>("saveSystem", false)),
times_(0),
updatePrecondFreq_(saveSystem_ ? readInt(fvSolution.subDict("solvers").subDict(word(T.name())).lookup("updatePrecondFrequency")) : 1000),
updateA_(saveSystem_ ? readBool(fvSolution.subDict("solvers").subDict(word(T.name())).lookup("updateMatrixCoeffs")) : false),
nValidCmp_(0),
sumofA_(0.),
prefix_(word( word(T.name() + "_"))),
initTimeIndex(mesh.time().timeIndex()),
autoPrecond(false)
{   
  
  // Build inter-processor data if needed               
  sparseSolver<Type>::buildSharedDataOnDemand(mesh); 
  
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
   
    int argc = 0;
    char** argv = new char*[argc];
    argv[0] = NULL;
 
    PetscInitialize(&argc,&argv,petscOptFile.c_str(),NULL);
    
    delete [] argv; 
  }
  
  // Build/set structures if the system is to be saved  
  if (this->saveSystem_)
  {   
    typename pTraits<Type>::labelType validComponents
    (
      mesh.template validComponents<Type>()
    );

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
      if (component(validComponents, cmpt) == -1) continue; 
      
      listA_.append(new Mat);
      listSolvers_.append(new KSP); 
      nValidCmp_++;
    }
    
    listb_.append(new Vec);
    listx_.append(new Vec); 
    
    forAll(listA_, i)
    {     
      ierr = MatCreate(PETSC_COMM_WORLD, &listA_[i]);CHKERRV(ierr);
      ierr = KSPCreate(PETSC_COMM_WORLD, &listSolvers_[i]); CHKERRV(ierr);
    }
    
    ierr = VecCreate(PETSC_COMM_WORLD,&listb_[0]);CHKERRV(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&listx_[0]);CHKERRV(ierr);
  } 
  
  this->counterPetsc_++;   
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::petscSolver<Type>::solve
(
    fvMatrix<Type>& matrix
)
{
  // Do some checks and updates if the mesh topology changes
  if (matrix.psi().mesh().topoChanging())
  {
     // The system cannot be saved for a mesh changing its topo. Abort.
     if (this->saveSystem_)
     {
       if (Pstream::master)
       {
         FatalErrorIn
         (
            "Foam::petscSolver<Type>::solve"
         )   << nl << nl << "fvMatrix(" << matrix.psi().name() << ") is set to be saved "
            << "over time, but the mesh topology is changing. Set saveSystem option to false"
            << " in dictionary fvSolution.solvers(" << matrix.psi().name() <<")." << nl << nl << endl
            << exit(FatalError);
       }
     }
     else
     {
       // If the mesh topo is changing we need to update the inter-processor
       // struct each time solve() is called. Actually, we only need to update
       // each time mesh changes (we could skip udpdate after it is called by
       // the first eqn to be solved).
       sparseSolver<Type>::setSharedData(matrix.psi().mesh());       
     } 
  }
  
  if (saveSystem_)
  {
    solveReuse(matrix);
  }
  else
  {
    solveNoReuse(matrix);
  }
  
  times_++;
}


template<class Type>
void Foam::petscSolver<Type>::solveReuse
(
    fvMatrix<Type>& matrix
)
{
 
 GeometricField<Type, fvPatchField, volMesh>& T =
 const_cast< GeometricField<Type, fvPatchField, volMesh>& >
   (
      matrix.psi()
   ); 

 typename pTraits<Type>::labelType validComponents
 (
    matrix.psi().mesh().template validComponents<Type>()
 );
 
 // Cmp-independent re-used PETSC stuff
 PetscInt  its;
 KSPType   ksptype;
 PCType    pctype;
 PC        pc;
 
 scalarField saveDiag(matrix.lduMatrix::diag());

 Field<Type> source(matrix.source());
 
 this->addBoundarySource(source, matrix, T);
 
 int i = 0;
 for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
 {
   int cmpI(cmpt);
   if (component(validComponents, cmpt) == -1) continue;
   
   // The steps below up to the call to this->getFoamResiduals() are only needed
   // to compute the residuals using the Foam machinery and definition. They
   // are not needed for the solve step because the boundary contributions
   // in this class are explicitly added in the assemble() rountines.     
   
   scalarField psiCmpt(T.primitiveField().component(cmpt));
   scalarField sourceCmpt(source.component(cmpt));
   
   FieldField<Field, scalar> bouCoeffsCmpt
   (
     matrix.boundaryCoeffs().component(cmpt)
   );

   FieldField<Field, scalar> intCoeffsCmpt
   (
     matrix.internalCoeffs().component(cmpt)
   );

   lduInterfaceFieldPtrsList interfaces =
      T.boundaryField().scalarInterfaces();
      
   matrix.initMatrixInterfaces
   (
     bouCoeffsCmpt,
     interfaces,
     psiCmpt,
     sourceCmpt,
     cmpt
   );

   matrix.updateMatrixInterfaces
   (
     bouCoeffsCmpt,
     interfaces,
     psiCmpt,
     sourceCmpt,
     cmpt
   );
    
   scalar initResidual = 
   this->getFoamResiduals
   (
     T,
     matrix,
     sourceCmpt,
     psiCmpt, 
     saveDiag,
     bouCoeffsCmpt,
     interfaces,
     this->nEvalInit_,
     this->saveSystem_,
     cmpt,
     i   
   ); 
  
   // Now we start PETSC related stuff 
   label timeID = T.time().timeIndex()- initTimeIndex;
   
   if ( timeID < this->nEvalInit_ || updateA_)
   {  
     // Last arg: only need to size/allocate for first time-step.
     // Updating matrix coeffs does not change allocation.
     assemblePetscAbx
     (
       listA_[i],
       listb_[0],
       listx_[0],
       matrix,
       T,
       cmpI,
       timeID < this->nEvalInit_ ? 0 : -1
     );  
   }
   else
   {
     assemblePetscbx
     (
       listb_[0],
       listx_[0],
       matrix,
       T,
       cmpI,
       i 
     ); 
   }
   
   ierr = KSPSetOptionsPrefix(listSolvers_[i],prefix_.c_str());CHKERRV(ierr);
   
   // Set operators and solve    
   if (timeID < this->nEvalInit_ || times_ > updatePrecondFreq_)
    {
      ierr = KSPSetReusePreconditioner(listSolvers_[i], PETSC_FALSE); CHKERRV(ierr);
      ierr = KSPSetOperators(listSolvers_[i], listA_[i], listA_[i]); CHKERRV(ierr);
   
      // Reset counter
      if (nValidCmp_ == i+1)
       times_ = 1;    
    }
   else
    {
      ierr = KSPSetReusePreconditioner(listSolvers_[i], PETSC_TRUE); CHKERRV(ierr);
    }
    
   ierr = KSPSetFromOptions(listSolvers_[i]);CHKERRV(ierr);
   ierr = KSPGetType(listSolvers_[i], &ksptype);CHKERRV(ierr);
   
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
          "Foam::petscSolver<Type>::solveReuse"
        )   << nl << nl << "You cannot freeze the sparse matrix solver for field " 
        <<  T.name() << " while the matrix is being"
        << " updated, for a direct solver (preonly). Either set 'updatePrecondFrequency' to 1 or"
        << " choose an iterative matrix solver." << nl << nl << endl
        << exit(FatalError);
       }
     }
   } 
    
   auto start = std::chrono::high_resolution_clock::now();
   
   // Solve
   ierr = KSPSolve(listSolvers_[i],listb_[0],listx_[0]);CHKERRV(ierr);
   
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
   
   if (autoPrecond && nValidCmp_ == i+1)
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
  
   ierr = KSPGetIterationNumber(listSolvers_[i],&its);CHKERRV(ierr);
  
   ierr = KSPGetPC(listSolvers_[i], &pc);
   ierr = PCGetType(pc, &pctype);CHKERRV(ierr);
   
   // Get PETSC
   transferPetscSolution(listx_[0], T, cmpI);   
   
   // Compute final residuals using Foam definition
   scalar finalResidual = 
   this->getFoamResiduals
   (
     T,
     matrix,
     sourceCmpt,
     T.primitiveField().component(cmpt), 
     saveDiag,
     bouCoeffsCmpt,
     interfaces,
     this->nEvalInit_,
     this->saveSystem_,
     cmpt,
     i    
   ); 
   
   // Print solver info (solver/PC name, iters, residuals) 
   word cmpName(T.name() + pTraits<Type>::componentNames[cmpt]);  
   Info << "Petsc:" << ksptype << ":" << pctype << ": Solving for "<< cmpName  << ", Initial residual = " << initResidual  
         << ", Final residual = " << finalResidual << ", No Iterations " << its << endl; 
       
   i++;   
 }
 
 T.correctBoundaryConditions(); 
}

template<class Type>
void Foam::petscSolver<Type>::solveNoReuse
(
    fvMatrix<Type>& matrix
)
{

 GeometricField<Type, fvPatchField, volMesh>& T =
 const_cast< GeometricField<Type, fvPatchField, volMesh>& >
 (
    matrix.psi()
 ); 

 typename pTraits<Type>::labelType validComponents
 (
    matrix.psi().mesh().template validComponents<Type>()
 );
 
 // Create A,b,x and solver/PC outside the components loop
 Mat A;
 Vec b,x;
 KSP ksp;
 PC pc;
   
 ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRV(ierr);
 ierr = VecCreate(PETSC_COMM_WORLD, &b);CHKERRV(ierr);
 ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRV(ierr);
 ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRV(ierr);
 
 scalarField saveDiag(matrix.lduMatrix::diag());

 Field<Type> source(matrix.source());
 
 this->addBoundarySource(source, matrix, T);
 
 int i = 0;
 for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
 {
   int cmpI(cmpt);
   if (component(validComponents, cmpt) == -1) continue; 
   
   // The steps below up to the call to this->getFoamResiduals() are only needed
   // to compute the residuals using the Foam machinery and definition. They
   // are not needed for the solve step because the boundary contributions
   // in this class are explicitly added in the assemble() rountines.     
   
   scalarField psiCmpt(T.primitiveField().component(cmpt));
   scalarField sourceCmpt(source.component(cmpt));
   
   FieldField<Field, scalar> bouCoeffsCmpt
   (
     matrix.boundaryCoeffs().component(cmpt)
   );

   FieldField<Field, scalar> intCoeffsCmpt
   (
     matrix.internalCoeffs().component(cmpt)
   );

   lduInterfaceFieldPtrsList interfaces =
      T.boundaryField().scalarInterfaces();
      
   matrix.initMatrixInterfaces
   (
     bouCoeffsCmpt,
     interfaces,
     psiCmpt,
     sourceCmpt,
     cmpt
   );

   matrix.updateMatrixInterfaces
   (
     bouCoeffsCmpt,
     interfaces,
     psiCmpt,
     sourceCmpt,
     cmpt
   );
    
   scalar initResidual = 
   this->getFoamResiduals
   (
     T,
     matrix,
     sourceCmpt,
     psiCmpt, 
     saveDiag,
     bouCoeffsCmpt,
     interfaces,
     this->nEvalInit_,
     this->saveSystem_,
     cmpt,
     i
   ); 
   
   // Now we start PETSC related stuff      
   PetscInt  its;
   KSPType  ksptype;
   PCType  pctype;
   
   assemblePetscAbx
   (
     A,
     b,
     x,
     matrix,
     T,
     cmpI, 
     i
   );  
    
   ierr = KSPSetOptionsPrefix(ksp,prefix_.c_str());CHKERRV(ierr);
  
   ierr = KSPSetOperators(ksp,A,A); CHKERRV(ierr);
   
   ierr = KSPSetFromOptions(ksp);CHKERRV(ierr);
 
   ierr = KSPSolve(ksp,b,x);CHKERRV(ierr);
  
   ierr = KSPGetIterationNumber(ksp,&its);CHKERRV(ierr);
  
   ierr = KSPGetType(ksp, &ksptype);CHKERRV(ierr);
   ierr = KSPGetPC(ksp, &pc);
   ierr = PCGetType(pc, &pctype);CHKERRV(ierr);
   
   transferPetscSolution(x, T, cmpI);   
   
   // Compute final residuals using Foam definition
   scalar finalResidual = 
   this->getFoamResiduals
   (
     T,
     matrix,
     sourceCmpt,
     T.primitiveField().component(cmpt), 
     saveDiag,
     bouCoeffsCmpt,
     interfaces,
     this->nEvalInit_,
     this->saveSystem_,
     cmpt,
     i    
   ); 
   
   // Print solver info (solver/PC name, iters, residuals) 
   word cmpName(T.name() + pTraits<Type>::componentNames[cmpt]);  
   Info << "Petsc:" << ksptype << ":" << pctype << ": Solving for "<< cmpName  << ", Initial residual = " << initResidual  
         << ", Final residual = " << finalResidual << ", No Iterations " << its << endl; 
  
   i++;
 }
 
 // Destroy
 KSPDestroy(&ksp); 
 VecDestroy(&b);  VecDestroy(&x);
 MatDestroy(&A);
 
 T.correctBoundaryConditions(); 
}

template<class Type>
void Foam::petscSolver<Type>::computeAllocationPetsc
(
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T
)
{
 // Set Allocation. Note: this will only work for direct-neig based
 // stencil, ie, each cell only 'depends' on its direct neigs.
 
 int n = T.size();
 int nglb = T.size();
 reduce(nglb, sumOp<int>()); 
 
 // Resize arrays to 0 because setSize() might be conservative
 if (T.mesh().topoChanging())
 {
   maxInProcFaces_.clear();
   maxOutProcFaces_.clear();
 }
 
 maxInProcFaces_.setSize(n, 1); // (n,n) diagonal elements should always exist 
 maxOutProcFaces_.setSize(n, 0);
 
 const labelUList& owner = T.mesh().owner();
 const labelUList& neig = T.mesh().neighbour();

 //- Diagonal of MPI matrix
 forAll(owner, facei)
 {
   maxInProcFaces_[owner[facei]] += 1;
   maxInProcFaces_[neig[facei]] += 1;
 }
 
 //- Off-diagonal of MPI matrix
 if (Pstream::parRun())
 {
  forAll(T.boundaryField(), patchI)
  { 
   if (T.boundaryField()[patchI].coupled())
   {
     const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI); 
     
     if (isType<cyclicFvPatch>(T.mesh().boundary()[patchI]))
     {
       forAll(addr, faceI)
       {
         maxInProcFaces_[addr[faceI]] += 1;
       }
     }
     else
     {
       forAll(addr, faceI)
       {
         maxOutProcFaces_[addr[faceI]] += 1;
       }
     }
   }
  }
 }
 else
 {
   forAll(T.boundaryField(), patchI)
   { 
     if (T.boundaryField()[patchI].coupled())
     {
       const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI); 
    
       if (isType<cyclicFvPatch>(T.mesh().boundary()[patchI]))
       {
         forAll(addr, faceI)
         {
           maxInProcFaces_[addr[faceI]] += 1;
         }
       }
     }    
    }
 }
  
}

// Note: in functions assemblePetscAbx() and assemblePetscbx() we mix 
// INSERT_VALUES and ADD_VALUES without intermediate MAT_FLUSH_ASSEMBLY.
// However, the code is structured such that for a given matrix position:
//  - we first call INSERT_VALUES at that position and this is done once
//    by the owning process;
//  - ADD_VALUES is only called after INSERT_VALUES has been called;
//  - INSERT_VALUES is never called if ADD_VALUES has already ben called 
//    at that position;
//  - each process only sets values on its own matrix elements. There is 
//    never inter-processor exchange of data.
// Therefore, we can break that PETSC's rule.    

template<class Type>
void Foam::petscSolver<Type>::assemblePetscAbx
(
  Mat& A,
  Vec& b,
  Vec& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI,
  int cmpV
)
{
 int ilower = this->sharedData.ilower;

 int n = T.size();
 int nglb = T.size();
 reduce(nglb, sumOp<int>()); 
 
 // Matrix/Vectors can be reused for other components. Only need
 // to size/preallocate upon call from first valid component.
 if (cmpV == 0)
 {
  ierr = MatSetOptionsPrefix(A,prefix_.c_str());CHKERRV(ierr);
  ierr = MatSetSizes(A,n,n,nglb,nglb);CHKERRV(ierr);
  ierr = MatSetFromOptions(A);CHKERRV(ierr);
  
  ierr = VecSetOptionsPrefix(b,prefix_.c_str());CHKERRV(ierr);
  ierr = VecSetSizes(b,n,nglb);CHKERRV(ierr);
  ierr = VecSetFromOptions(b);CHKERRV(ierr);
  
  ierr = VecSetOptionsPrefix(x,prefix_.c_str());CHKERRV(ierr);
  ierr = VecSetSizes(x,n,nglb);CHKERRV(ierr);
  ierr = VecSetFromOptions(x);CHKERRV(ierr);  
 
  //- Allocate. Compute allocation once for static meshes or 
  // every solve() for moving meshes.
  if (T.mesh().topoChanging() || maxInProcFaces_.size() == 0)
    computeAllocationPetsc(eqn, T);
 
  ierr = MatSeqAIJSetPreallocation(A,0,maxInProcFaces_.begin());
  ierr = MatMPIAIJSetPreallocation(A,0,maxInProcFaces_.begin(),0,maxOutProcFaces_.begin());CHKERRV(ierr);  
 
  // Use this (let PETSC allocate on its own) if the above is not working
  // ierr = MatSetUp(A);CHKERRV(ierr);
 }
 else
 {
  // Reuse matrix structure. Simply reset all coeffs.
  ierr = MatZeroEntries(A); CHKERRV(ierr);
  ierr = VecSet(b,0.); CHKERRV(ierr);
  // x is zeroed automatically before kspsolve(), unless -ksp_initial_guess_nonzero 1
 }
 
 // Start filling the matrix/vector
 
 //- Off diagonal elements   
 const lduAddressing& offDiag = eqn.lduMatrix::lduAddr();
 int col; int row;
 int nFaces = offDiag.lowerAddr().size();
 
 // Symmetric matrices only have upper(). No need to force creation of lower().
 if (eqn.symmetric())
 {
   for (register label face=0; face<nFaces; face++)
   { 
     row = ilower + offDiag.upperAddr()[face];
     col = ilower + offDiag.lowerAddr()[face];
   
     // ij Off-diagonal
     ierr = MatSetValues(A,1,&row,1,&col,&eqn.upper()[face],INSERT_VALUES);CHKERRV(ierr);
     // ji Off-diagonal
     ierr = MatSetValues(A,1,&col,1,&row,&eqn.upper()[face],INSERT_VALUES);CHKERRV(ierr);   
   }
 }
 else
 {
   for (register label face=0; face<nFaces; face++)
   {
     row = ilower + offDiag.upperAddr()[face];
     col = ilower + offDiag.lowerAddr()[face];
   
     // ij Off-diagonal
     ierr = MatSetValues(A,1,&row,1,&col,&eqn.lower()[face],INSERT_VALUES);CHKERRV(ierr);
     // ji Off-diagonal
     ierr = MatSetValues(A,1,&col,1,&row,&eqn.upper()[face],INSERT_VALUES);CHKERRV(ierr);   
   }
 }
 
 // Diagonal and source  
 scalarField source = eqn.source().component(cmpI);   
 for (int cellI=0; cellI<n; cellI++) 
 {  
   // Diagonal
   row = ilower + cellI;  
   ierr = MatSetValues(A,1,&row,1,&row,&eqn.diag()[cellI],INSERT_VALUES);CHKERRV(ierr);
   
   // Source vector   
   ierr = VecSetValues(b,1,&row,&source[cellI],INSERT_VALUES);CHKERRV(ierr);
 }
 
 //- Contribution from BCs
 forAll(T.boundaryField(), patchI)
 {
   scalarField bC = eqn.boundaryCoeffs()[patchI].component(cmpI);  
   scalarField iC = eqn.internalCoeffs()[patchI].component(cmpI);  
   const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI);      
    
   // Non-coupled
   if (!T.boundaryField()[patchI].coupled())
   {
     forAll(addr, facei)
     {
       // Matrix of coefs - Diagonal
       row = ilower + addr[facei];  
       ierr = MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);CHKERRV(ierr);
       
       // Source vector
       ierr = VecSetValues(b,1,&row,&bC[facei],ADD_VALUES);CHKERRV(ierr);       
     }
   }
   // Processor or processorCyclic.
   // Note: from exchange ops point of view both patch types are equivalent.
   // Only difference is when building the sharedData struct, where the names
   // used to find neig patch in neig proc are different.
   else if (   isType<processorFvPatch>(T.mesh().boundary()[patchI])
            || isType<processorCyclicFvPatch>(T.mesh().boundary()[patchI])
           )
   { 
     // Matrix of coefs - Diagonal
     forAll(addr, facei)
     {      
       row = ilower + addr[facei];  
       ierr = MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);CHKERRV(ierr);       
     } 
     
     // Matrix of coefs - off-diagonal (row -> this processor; col -> other processors)
     // Unique face contribution to off-diag (set instead of add). If use add, then the
     // inter-processor off-diag coefs need to be reset to 0 before.
     forAll(this->sharedData.procInfo, pI)
     {
       if (this->sharedData.procInfo[pI][2] == patchI)
       {
         forAll(bC, facei)
         {        
           double v = -bC[facei];
           ierr = MatSetValues
           (
             A,
             1,
             &this->sharedData.fCo[pI][facei],
             1,
             &this->sharedData.fCn[pI][facei],
             &v,
             INSERT_VALUES
           ); 
           CHKERRV(ierr);     
         } 
       }
     }  
   }
   // Cyclic. 
   // Note: cyclic patches (if not empty after decomposition) always contain
   // both own and neig cells in the same processor. Therefore, no data exchange 
   // is needed, and biasing with ilower is enough to account for global indexing.  
   else if (isType<cyclicFvPatch>(T.mesh().boundary()[patchI]))
   { 
     // Row
     const labelList& owfC = T.mesh().boundaryMesh()[patchI].faceCells();
     
     const fvPatch& cyclicPatch = T.mesh().boundary()[patchI];
    
     const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();
      
     // Col   
     const labelList& nbFC = nbrPatch.patch().faceCells();
     
     forAll(owfC, facei)
     {           
       row = ilower + owfC[facei];  
       col = ilower + nbFC[facei]; 
        
       // Matrix of coefs - Diagonal  
       ierr = MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);CHKERRV(ierr);  
       // Matrix of coefs - off-diagonal
       double v = -bC[facei];
       ierr = MatSetValues(A,1,&row,1,&col,&v,INSERT_VALUES);CHKERRV(ierr);       
     } 
       
   }
    
 } 
 
 ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
 ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
 
 ierr = VecAssemblyBegin(b);CHKERRV(ierr);
 ierr = VecAssemblyEnd(b);CHKERRV(ierr);

 if (sparseSolver<Type>::debug_) 
 {
   std::string timeOut = std::to_string(T.mesh().time().timeIndex());
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
} 

template<class Type>
void Foam::petscSolver<Type>::assemblePetscbx
(
  Vec& b,
  Vec& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI,
  int cmpV
)
{
 int ilower = this->sharedData.ilower;

 int n = T.size();
 int nglb = T.size();
 reduce(nglb, sumOp<int>()); 
 
 // Matrix/Vectors can be reused for other components. Only need
 // to size/preallocate upon call from first valid component.
 if (cmpV == 0)
 {  
   ierr = VecSetOptionsPrefix(b,prefix_.c_str());CHKERRV(ierr);
   ierr = VecSetSizes(b,n,nglb);CHKERRV(ierr);
   ierr = VecSetFromOptions(b);CHKERRV(ierr);
 
   ierr = VecSetOptionsPrefix(x,prefix_.c_str());CHKERRV(ierr);
   ierr = VecSetSizes(x,n,nglb);CHKERRV(ierr);
   ierr = VecSetFromOptions(x);CHKERRV(ierr);
 }
 else
 {
   // Reuse matrix structure. Simply reset all coeffs.
   ierr = VecSet(b,0.); CHKERRV(ierr);
   // x is zeroed automatically before kspsolve(), unless -ksp_initial_guess_nonzero 1
 }
 
 // Start filling the vector
 int row;
  
 // Source  
 scalarField source = eqn.source().component(cmpI);   
 for (int cellI=0; cellI<n; cellI++) 
 {  
   row = ilower + cellI;   
   ierr = VecSetValues(b,1,&row,&source[cellI],INSERT_VALUES);CHKERRV(ierr);
 }
 
 //- Contribution from BCs
 forAll(T.boundaryField(), patchI)
 {
   scalarField bC = eqn.boundaryCoeffs()[patchI].component(cmpI);   
   const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI);      
    
   // Non-coupled
   if (!T.boundaryField()[patchI].coupled())
   {
     forAll(addr, facei)
     {
       row = ilower + addr[facei];  
       ierr = VecSetValues(b,1,&row,&bC[facei],ADD_VALUES);CHKERRV(ierr);       
     }
   }    
 } 
  
 ierr = VecAssemblyBegin(b);CHKERRV(ierr);
 ierr = VecAssemblyEnd(b);CHKERRV(ierr);

 if (sparseSolver<Type>::debug_) 
 {    
   // RHS vector
   {    
     PetscViewer viewer;
     PetscViewerASCIIOpen(PETSC_COMM_WORLD, "b.txt", &viewer);
     PetscViewerPushFormat(viewer,PETSC_VIEWER_DEFAULT);
     VecView(b,viewer);
     PetscViewerPopFormat(viewer);
     PetscViewerDestroy(&viewer);
   }    
 } 
} 

template<class Type>
void Foam::petscSolver<Type>::transferPetscSolution
(
  Vec& x,
  GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{
   int ilower = this->sharedData.ilower;
   
   /* get the local solution */
   int nvalues = T.size();
   std::vector<int> rows(nvalues);
        
   for (int i = 0; i < nvalues; i++)
    rows[i] = ilower + i;
    
   scalarField tt(T.size(), 0.);
   VecGetValues(x, nvalues, &rows[0], tt.begin());  
   
   T.primitiveFieldRef().replace(cmpI, tt);
} 
 
template<class Type>
void Foam::petscSolver<Type>::checkMatrixSum
(
    const scalarField& rowSum,
    const word name,
    const int tindex, 
    const int vcmpt   
) 
{ 
  if (tindex < this->nEvalInit_ && vcmpt == 0)
  {
    sumofA_ = gSumMag(rowSum);
  }
  else if (tindex == this->nEvalInit_ && vcmpt == 0)
  {
   scalar dif = mag(sumofA_ - gSumMag(rowSum) );
   if ( dif > 1e-20 && !updateA_ )
    {
     if (Pstream::master)
     {
      FatalErrorIn
        (
            "Foam::petscSolver<Type>::checkMatrixSum"
        )   << nl << nl << "Sum of fvMatrix(" << name << ") is changing by " << dif
            << " between consecutive time-steps and option updateMatrixCoeffs is disabled."
            << " You should enable this option or set option saveSystem to false." << nl << nl << endl
            << exit(FatalError);
     }
    }
  }
    
}
// ************************************************************************* //
