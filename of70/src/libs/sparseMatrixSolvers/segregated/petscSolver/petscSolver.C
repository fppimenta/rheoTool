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
#include "cyclicAMIFvPatch.H"   
#include "cyclicFvPatchField.H" 
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
autoPrecond(false),
isThereCyclicAMI_(false)
{   
  // Verify limitations of the interface
  checkLimitations(T);
  
  // Build inter-processor data if needed               
  meshID_ = sparseSolver<Type>::buildSharedDataOnDemand(mesh); 
  
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
void Foam::petscSolver<Type>::checkLimitations
(
  const GeometricField<Type, fvPatchField, volMesh>& T
)
{

 bool isScalar(pTraits<Type>::rank == 0 ? true : false);

 forAll(T.boundaryField(), patchI)
 { 
   const fvPatchField<Type>& pfield = T.boundaryField()[patchI];
   const fvPatch& pfvPatch = T.mesh().boundary()[patchI];
 
   if (pfield.coupled() && pfvPatch.size() > 0)
   {     
     if (isType<cyclicFvPatch>(pfvPatch))
     {  
       // Partial support for rotational transform in non-scalar variables (coupling is explicit)     
       if (!isScalar &&
            refCast<const cyclicFvPatch>(pfvPatch).cyclicPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           WarningInFunction
           << nl << "Because patch " << pfvPatch.name() 
           << " is rotational cyclic, the matrix connection between the cyclic patches is explicit."
           << " This is innacurate for transient computations and may generate local instabilities."
           << " Problem raised while solving for field: " << T.name()   
           << nl << endl;
       }
     }
     else if (isType<processorCyclicFvPatch>(pfvPatch))
     {
       // No support for rotational transform in non-scalar variables     
       if (!isScalar &&
            refCast<const processorCyclicFvPatch>(pfvPatch).procPolyPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           FatalErrorInFunction
           << "No support in Petsc interface for patches of type cyclic rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << " Problem raised while solving for field: " << T.name() 
           << exit(FatalError);
       }     
     }
     else if (isType<cyclicAMIFvPatch>(pfvPatch))
     {
       // No support for rotational transform in non-scalar variables     
       if (!isScalar &&
            refCast<const cyclicAMIFvPatch>(pfvPatch).cyclicAMIPatch().transform()
            ==
            coupledPolyPatch::transformType::ROTATIONAL
          )
       {
           FatalErrorInFunction
           << "No support in Petsc interface for patches of type cyclicAMI rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << " Problem raised while solving for field: " << T.name() 
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
       << " Problem raised while solving for field: " << T.name() 
       << exit(FatalError);  
     }     
   }
 }


}



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
   // in this class are directly added in the assemble() rountines.     
   
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
 maxInProcFaces_.clear();
 maxOutProcFaces_.clear();
 
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
 forAll(T.boundaryField(), patchI)
 {   
   if (T.boundaryField()[patchI].coupled())
   {
     const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI); 
     const fvPatch& pfvPatch = T.mesh().boundary()[patchI];
      
     // cyclic
     if (isType<cyclicFvPatch>(pfvPatch))
     {
       forAll(addr, faceI)
       {
         maxInProcFaces_[addr[faceI]] += 1;
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
                     maxInProcFaces_[ownFC[facei]] += 1;
                   }
                   else
                   {
                     maxOutProcFaces_[ownFC[facei]] += 1;                 
                   }
                 }
                 else
                 {
                   maxOutProcFaces_[ownFC[facei]] += 1;                 
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
                 maxInProcFaces_[ownFC[facei]] += 1;
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
                     maxInProcFaces_[ownFC[facei]] += 1;
                   }
                   else
                   {                 
                     maxOutProcFaces_[ownFC[facei]] += 1;
                   }
                 }
                 else
                 {                 
                   maxOutProcFaces_[ownFC[facei]] += 1;
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
                 maxInProcFaces_[ownFC[facei]] += 1;
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
         maxOutProcFaces_[addr[faceI]] += 1;
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
 int ilower = this->sharedData[meshID_].ilower;

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
  
  // Each process must only assemble its own values (error otherwise) and the allocation
  // must be exact.
  ierr = MatSetOption(A,MAT_NO_OFF_PROC_ENTRIES,PETSC_TRUE);CHKERRV(ierr);
  
  ierr = VecSetOptionsPrefix(b,prefix_.c_str());CHKERRV(ierr);
  ierr = VecSetSizes(b,n,nglb);CHKERRV(ierr);
  ierr = VecSetFromOptions(b);CHKERRV(ierr);
  
  ierr = VecSetOptionsPrefix(x,prefix_.c_str());CHKERRV(ierr);
  ierr = VecSetSizes(x,n,nglb);CHKERRV(ierr);
  ierr = VecSetFromOptions(x);CHKERRV(ierr);  
  
  //- Allocate. Compute allocation once for static meshes or 
  // every solve() for moving meshes.
  if (T.mesh().topoChanging() || maxInProcFaces_.size() == 0 || (isThereCyclicAMI_ && T.mesh().changing()))
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
 scalarField source(eqn.source().component(cmpI));   
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
   scalarField bC(eqn.boundaryCoeffs()[patchI].component(cmpI));  
   scalarField iC(eqn.internalCoeffs()[patchI].component(cmpI));  
   const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI);  
   const fvPatch& pfvPatch = T.mesh().boundary()[patchI];    
    
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
   else if (   isType<processorFvPatch>(pfvPatch)
            || isType<processorCyclicFvPatch>(pfvPatch)
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
     forAll(this->sharedData[meshID_].procInfo, pI)
     {
       if (this->sharedData[meshID_].procInfo[pI][2] == patchI)
       {
         forAll(bC, facei)
         {        
           double v = -bC[facei];
           ierr = MatSetValues
           (
             A,
             1,
             &this->sharedData[meshID_].fCo[pI][facei],
             1,
             &this->sharedData[meshID_].fCn[pI][facei],
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
   else if (isType<cyclicFvPatch>(pfvPatch))
   { 
     // Row
     const labelList& ownFC = T.mesh().boundaryMesh()[patchI].faceCells();
     
     const fvPatch& cyclicPatch = T.mesh().boundary()[patchI];
    
     const fvPatch& nbrPatch = refCast<const cyclicFvPatch>(cyclicPatch).neighbFvPatch();
      
     // Col   
     const labelList& nbFC = nbrPatch.patch().faceCells();
         
     // Now we need to differentiate the treatment. If the variable is not a scalar and its transform type
     // is rotational, the coupling is explicit, ie, we multiply matrix coefficients by field values 
     // and send them to the source vector... 
     if (pTraits<Type>::rank != 0 &&
         refCast<const cyclicFvPatch>(pfvPatch).cyclicPatch().transform()
         ==
         coupledPolyPatch::transformType::ROTATIONAL
        )
     {           
       const labelUList& nbrFaceCells = refCast<const cyclicFvPatch>( cyclicPatch ).neighbFvPatch().faceCells();

       Field<Type> pnf(T.internalField(), nbrFaceCells);
       
       // Apply the rotation
       refCast<const cyclicFvPatchField<Type> >(T.boundaryField()[patchI]).transformCoupleField(pnf);

       scalarField pnfc(pnf.component(cmpI));
      
       forAll(ownFC, facei)
       {           
         row = ilower + ownFC[facei];  
        
         // Matrix of coefs - Diagonal  
         ierr = MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);CHKERRV(ierr);  
         // Source vector
         double v = bC[facei]*pnfc[facei];       
         ierr = VecSetValues(b,1,&row,&v,ADD_VALUES);CHKERRV(ierr); 
       }       
      }
      //... Otherwise, the coupling between cyclic patches is fully implicit and there are only contributions 
      // to the matrix of coefficients.
      else
      {
        forAll(ownFC, facei)
        {           
          row = ilower + ownFC[facei];  
          col = ilower + nbFC[facei]; 
        
          // Matrix of coefs - Diagonal  
          ierr = MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);CHKERRV(ierr);  
          // Matrix of coefs - off-diagonal
          double v = -bC[facei];
          ierr = MatSetValues(A,1,&row,1,&col,&v,INSERT_VALUES);CHKERRV(ierr);       
        } 
      }
     
   }
   // cyclicAMIFvPatch
   // Note: contrary to cyclic, cyclicAMI patches do not get transformed in processorCyclicAMI
   // when they are split between processors. The AMI machinery allows mapping between 
   // processors. 
   else if (isType<cyclicAMIFvPatch>(pfvPatch))
   {           
     // camipp and neicampi always lay in the same processor, but one of them can be empty
     // if the other half of a given AMI was sent to another processor. In pratice, nothing
     // happens bellow if camipp is zero-sized.
     const cyclicAMIPolyPatch& camipp = refCast<const cyclicAMIFvPatch>(pfvPatch).cyclicAMIPatch();
     const cyclicAMIPolyPatch& neicamipp = camipp.neighbPatch();
     
     const labelList& ownFC = camipp.faceCells();   
     const labelList& neiFC = neicamipp.faceCells();
     
     // The neiFC faceCells already get the proc's ilower. No further correction needed.
     labelList neiFCproc(neiFC+ilower);
     
     // AMIs are shared by at least 2 patches, but the AMI interpolator is only accessible
     // from one of them, which is considered the "owner" patch.
     if (camipp.owner())
     {      
       forAll(camipp.AMIs(), i)
       {
         // Distribute faceCells of neighbour (tgt) patch if parallel. We get
         // the tgt faceCells transfered to the proc where the src patch (camipp) is.
         if (camipp.AMIs()[i].singlePatchProc() == -1) 
           camipp.AMIs()[i].tgtMap().distribute(neiFCproc);   
     
         const labelListList& srcAd = camipp.AMIs()[i].srcAddress();
         const scalarListList& srcW = camipp.AMIs()[i].srcWeights();
         forAll(srcAd, k)
         {             
           // Matrix of coefs - Diagonal  
           row = ilower + ownFC[k];  
           ierr = MatSetValues(A,1,&row,1,&row,&iC[k],ADD_VALUES);CHKERRV(ierr); 
         
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (camipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (camipp.AMIs()[i].srcWeightsSum()[k] < camipp.AMIs()[i].lowWeightCorrection())
             {
               col = row;
               double v = -bC[k];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr); 
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(srcAd[k], kk)
               {                        
                col = neiFCproc[srcAd[k][kk]];
                double v = -bC[k]*srcW[k][kk];
                ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);          
               }
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(srcAd[k], kk)
              {                        
               col = neiFCproc[srcAd[k][kk]];
               double v = -bC[k]*srcW[k][kk];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);          
              }   
           }
         }
         
       }
  
     }
     else
     { 
      forAll(neicamipp.AMIs(), i)
       {
         // Distribute faceCells of neighbour (src) patch if parallel. We get
         // the src faceCells transfered to the proc where the tgt patch (camipp) is
         if (neicamipp.AMIs()[i].singlePatchProc() == -1) 
           neicamipp.AMIs()[i].srcMap().distribute(neiFCproc);   
       
         const labelListList& tgtAd = neicamipp.AMIs()[i].tgtAddress();
         const scalarListList& tgtW = neicamipp.AMIs()[i].tgtWeights();
         forAll(tgtAd, k)
         {                         
           // Matrix of coefs - Diagonal 
           row = ilower + ownFC[k];   
           ierr = MatSetValues(A,1,&row,1,&row,&iC[k],ADD_VALUES);CHKERRV(ierr);            
                     
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (neicamipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (neicamipp.AMIs()[i].tgtWeightsSum()[k] < neicamipp.AMIs()[i].lowWeightCorrection())
             {
               col = row;
               double v = -bC[k];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr); 
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(tgtAd[k], kk)
               {                
                col = neiFCproc[tgtAd[k][kk]];
                double v = -bC[k]*tgtW[k][kk];
                ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);          
               }  
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(tgtAd[k], kk)
              {                
               col = neiFCproc[tgtAd[k][kk]];
               double v = -bC[k]*tgtW[k][kk];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);          
              }   
           }
              
         }
       }
    }
    
   }
   // NOT IMPLEMENTED 
   else
   {   
     FatalErrorInFunction
     << "Sorry, but Petsc interface cannot support patches of type: "
     << pfvPatch.type() << ", which is currently set for "
     << "patch " << pfvPatch.name() << "."
     << exit(FatalError);   
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
 int ilower = this->sharedData[meshID_].ilower;

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
 scalarField source(eqn.source().component(cmpI));   
 for (int cellI=0; cellI<n; cellI++) 
 {  
   row = ilower + cellI;   
   ierr = VecSetValues(b,1,&row,&source[cellI],INSERT_VALUES);CHKERRV(ierr);
 }
 
 //- Contribution from BCs
 forAll(T.boundaryField(), patchI)
 {
   scalarField bC(eqn.boundaryCoeffs()[patchI].component(cmpI));   
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
   // Cylic Bcs only contribute to source vector if transform is rotational
   // and variable has rank > 0
   else if (isType<cyclicFvPatch>(T.mesh().boundary()[patchI]) &&
           pTraits<Type>::rank != 0 &&
           refCast<const cyclicFvPatch>(T.mesh().boundary()[patchI]).cyclicPatch().transform()
           ==
           coupledPolyPatch::transformType::ROTATIONAL
           )
   {   
     const fvPatch& cyclicPatch = T.mesh().boundary()[patchI];
       
     const labelUList& nbrFaceCells = refCast<const cyclicFvPatch>( cyclicPatch ).neighbFvPatch().faceCells();

     Field<Type> pnf(T.internalField(), nbrFaceCells);
         
     refCast<const cyclicFvPatchField<Type> >(T.boundaryField()[patchI]).transformCoupleField(pnf);

     scalarField pnfc(pnf.component(cmpI));
      
     const labelList& ownFC = T.mesh().boundaryMesh()[patchI].faceCells(); 
     forAll(ownFC, facei)
     {           
       row = ilower + ownFC[facei];  
       double v = bC[facei]*pnfc[facei];       
       ierr = VecSetValues(b,1,&row,&v,ADD_VALUES);CHKERRV(ierr); 
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
   int ilower = this->sharedData[meshID_].ilower;
   
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
