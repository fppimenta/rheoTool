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

#include "hypreSolver.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"
#include "cyclicAMIFvPatch.H"   
#include "cyclicFvPatchField.H" 
#include <chrono>

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

// This is the time index until which we force the preconditioner to be updated
// when saveSystem is enabled.  
template<class Type>
int Foam::hypreSolver<Type>::nEvalInit_ = 3;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::hypreSolver<Type>::hypreSolver
(
    const GeometricField<Type, fvPatchField, volMesh>& T,
    const fvMesh& mesh,
    const dictionary& fvSolution
)
:
sparseSolver<Type>(T, mesh, fvSolution),
saveSystem_(fvSolution.subDict("solvers").subDict(word(T.name())).lookupOrDefault<Switch>("saveSystem", false)),
HPsolver_(NULL),
HPprecond_(NULL),
times_(0),
updatePrecondFreq_(saveSystem_ ? readInt(fvSolution.subDict("solvers").subDict(word(T.name())).lookup("updatePrecondFrequency")) : 1000),
updateA_(saveSystem_ ? readBool(fvSolution.subDict("solvers").subDict(word(T.name())).lookup("updateMatrixCoeffs")) : false),
nValidCmp_(0),
sumofA_(0.),
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
  // Direct MPI_Init instead of Foam::UPstream::init()
  // to keep Foam structure as if it was a serial run
  // but runs Hypre in mpi singleton (needed).
  if ( this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
    MPI_Init(NULL,NULL); 
    
  // Build/set structures if the system is to be saved  
  if (this->saveSystem_)
  {
    HPsolver_.reset
    (
       HypreIJSolver::New(fvSolution.subDict("solvers").subDict(word(T.name()))).ptr()
    );
    
    HPprecond_.reset
    (
       HypreIJPreconditioner::New(fvSolution.subDict("solvers").subDict(word(T.name()))).ptr()
    );
      
    label ilower = this->sharedData[meshID_].ilower;
    label iupper = this->sharedData[meshID_].iupper;
 
    typename pTraits<Type>::labelType validComponents
    (
      mesh.template validComponents<Type>()
    );

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
      if (component(validComponents, cmpt) == -1) continue; 
      
      listA_.append(new HYPRE_IJMatrix);
      listSolvers_.append(new HYPRE_Solver); 
      listPrecond_.append(new HYPRE_Solver); 
      nValidCmp_++;
    }
    
    listb_.append(new HYPRE_IJVector);
    listx_.append(new HYPRE_IJVector); 
    
    forAll(listA_, i)
    {
      HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &listA_[i]);
      HYPRE_IJMatrixSetObjectType(listA_[i], HYPRE_PARCSR);
         
      // Create solver and preconditioner
      HPsolver_->initialize(listSolvers_[i]);
      HPprecond_->initialize(listPrecond_[i]);
    }
    
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&listb_[0]);
    HYPRE_IJVectorSetObjectType(listb_[0], HYPRE_PARCSR);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&listx_[0]);
    HYPRE_IJVectorSetObjectType(listx_[0], HYPRE_PARCSR);
  }
  
  this->counterHypre_++;
  
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::hypreSolver<Type>::checkLimitations
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
           << "No support in Hypre interface for patches of type cyclic rotational, which is the type"
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
           << "No support in Hypre interface for patches of type cyclicAMI rotational, which is the type"
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
       << "Sorry, but Hypre interface cannot support patches of type: "
       << pfvPatch.type() << ", which is currently set for "
       << "patch " << pfvPatch.name() << "."
       << " Problem raised while solving for field: " << T.name() 
       << exit(FatalError);  
     }     
   }
 }


}

template<class Type>
void Foam::hypreSolver<Type>::solve
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
            "Foam::hypreSolver<Type>::solve"
         )   << nl << nl << "fvMatrix(" << matrix.psi().name() << ") is set to be saved "
            << "over time, but the mesh topology is changing. Set saveSystem option to false"
            << " in dictionary fvSolution.solvers(" << matrix.psi().name() <<")." << nl << nl << endl
            << exit(FatalError);
       }
     }
     else
     {
       // If the mesh topo is changing we need to update the inter-processor
       // struct each time solve() is called (only applies for parRun).
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
void Foam::hypreSolver<Type>::solveReuse
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
  
   // Now we start HYPRE related stuff
   label timeID = T.time().timeIndex()- initTimeIndex;
   
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_ParVector par_b;
   HYPRE_ParVector par_x; 
   
   if (timeID < this->nEvalInit_ || updateA_)
   {
     assembleHypreAbx
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
     assembleHypreBx
     (
       listb_[0],
       listx_[0],
       matrix,
       T,
       cmpI 
     );
   }
    
   // Get the parcsr matrix object to use (ugly syntax with casts) 
   HYPRE_IJMatrixGetObject(listA_[i], static_cast<void**>(static_cast<void*>(&parcsr_A)));  
   HYPRE_IJVectorGetObject(listb_[0], static_cast<void**>(static_cast<void*>(&par_b))); 
   HYPRE_IJVectorGetObject(listx_[0], static_cast<void**>(static_cast<void*>(&par_x)));  
   
   // Now setup... 
   auto start = std::chrono::high_resolution_clock::now();
    
   if (timeID < this->nEvalInit_ || times_ > updatePrecondFreq_)
    {
      HPsolver_->setup
      (
        listSolvers_[i],
        listPrecond_[i], 
        parcsr_A, 
        par_b, 
        par_x, 
        HPprecond_->solvePtr(), 
        HPprecond_->setupPtr(),
        matrix.symmetric()
      ); 
      
      // Reset counter
      if (nValidCmp_ == i+1)
       times_ = 1;    
    }
   
   // ...and solve          
   word cmpName(T.name() + pTraits<Type>::componentNames[cmpt]);                              
   HPsolver_->solve(listSolvers_[i], parcsr_A, par_b, par_x); 
   
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

   // Transfer solution to OF 
   transferHypreSolution(listx_[0], T, cmpI);
   
   //  Print the system for debug 
   if (sparseSolver<Type>::debug_) 
    this->printSystem(listA_[i], listb_[0], cmpName);
    
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
   HPsolver_->printInfo
   (
     cmpName,
     HPprecond_->precondName(),
     initResidual,
     finalResidual
   );
       
   i++;   
 }
 
 T.correctBoundaryConditions(); 
}

template<class Type>
void Foam::hypreSolver<Type>::solveNoReuse
(
    fvMatrix<Type>& matrix
)
{

 label ilower = this->sharedData[meshID_].ilower;
 label iupper = this->sharedData[meshID_].iupper;
 
 GeometricField<Type, fvPatchField, volMesh>& T =
 const_cast< GeometricField<Type, fvPatchField, volMesh>& >
 (
    matrix.psi()
 ); 

 typename pTraits<Type>::labelType validComponents
 (
    matrix.psi().mesh().template validComponents<Type>()
 );
 
 scalarField saveDiag(matrix.lduMatrix::diag());

 Field<Type> source(matrix.source());
 
 this->addBoundarySource(source, matrix, T);
 
 // Sparsity patterns do not change among components. Thus
 // create matrix/vectors once and simply change values
 HYPRE_IJMatrix A;
 HYPRE_IJVector b;
 HYPRE_IJVector x;
 
 HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
 HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  
 HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
 HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);

 HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
 HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
 
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
   
   // Now we start HYPRE related stuff   
  
   HYPRE_ParCSRMatrix parcsr_A;    
   HYPRE_ParVector par_b;  
   HYPRE_ParVector par_x;
   HYPRE_Solver solver, precond; 
   
   autoPtr<HypreIJSolver> HPsolver(HypreIJSolver::New(this->solDict_.subDict(word(T.name()))));
   autoPtr<HypreIJPreconditioner> HPprecond(HypreIJPreconditioner::New(this->solDict_.subDict(word(T.name()))));
   
   assembleHypreAbx
   (
     A,
     b,
     x,
     matrix,
     T,
     cmpI,
     i
   ); 

   // Get the parcsr matrix object to use (ugly syntax with casts) 
   HYPRE_IJMatrixGetObject(A, static_cast<void**>(static_cast<void*>(&parcsr_A)));  
   HYPRE_IJVectorGetObject(b, static_cast<void**>(static_cast<void*>(&par_b))); 
   HYPRE_IJVectorGetObject(x, static_cast<void**>(static_cast<void*>(&par_x))); 
      
   // Create solver and preconditioner 
   HPsolver->initialize(solver);
   HPprecond->initialize(precond);
    
   // Now setup and solve!   
   HPsolver->setup
   (
     solver,
     precond, 
     parcsr_A, 
     par_b, 
     par_x, 
     HPprecond->solvePtr(), 
     HPprecond->setupPtr(),
     matrix.symmetric()
   ); 
                      
   word cmpName(T.name() + pTraits<Type>::componentNames[cmpt]);                              
   HPsolver->solve(solver, parcsr_A, par_b, par_x);  

   // Transfer solution to OF 
   transferHypreSolution(x, T, cmpI);
   
   //  Print the system for debug 
   if (sparseSolver<Type>::debug_) 
    this->printSystem(A, b, cmpName);  
    
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
   HPsolver->printInfo
   (
     cmpName,
     HPprecond->precondName(),
     initResidual,
     finalResidual
   );
    
   // Destroy solver and preconditioner 
   HPsolver->destroy(solver); 
   HPprecond->destroy(precond);
   
   i++;
 }
 
 T.correctBoundaryConditions();
 
 // Clean up 
 HYPRE_IJMatrixDestroy(A);
 HYPRE_IJVectorDestroy(b);
 HYPRE_IJVectorDestroy(x);
 
}

template<class Type>
void Foam::hypreSolver<Type>::computeAllocationHypre
(
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T
)
{
 // Set Allocation. Note: this will only work for direct-neig based
 // stencil, ie, each cell only 'depends' on its direct neigs.
 
 label n = T.size();
 label nglb = T.size();
 reduce(nglb, sumOp<label>()); 
 
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

template<class Type>
void Foam::hypreSolver<Type>::printSystem
(
    const HYPRE_IJMatrix& A,
    const HYPRE_IJVector& b,
    const word& name
) const
{
   std::string outputAfile("IJ.out.A." + name); 
   std::string outputBfile("IJ.out.b." + name);
   
   // Print out the system  - files names will be IJ.out.A.var.XXXXX
   // and IJ.out.b.var.XXXXX, where XXXXX = processor id     
   HYPRE_IJMatrixPrint(A, outputAfile.c_str());
   HYPRE_IJVectorPrint(b, outputBfile.c_str());
} 

template<class Type>
void Foam::hypreSolver<Type>::assembleHypreAbx
(
  HYPRE_IJMatrix& A,
  HYPRE_IJVector& b,
  HYPRE_IJVector& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI,
  int cmpV
)
{

 // Matrix can be reused for other components. Only need
 // to size/preallocate upon call from first valid component.
 if (cmpV == 0)
 { 
  //- Allocate. Compute allocation once for static meshes or 
  // every solve() for moving meshes.
  if (T.mesh().topoChanging() || maxInProcFaces_.size() == 0 || (isThereCyclicAMI_ && T.mesh().changing()))
    computeAllocationHypre(eqn, T);
 
  if (Pstream::parRun())
  {
    HYPRE_IJMatrixSetDiagOffdSizes(A, maxInProcFaces_.begin(), maxOutProcFaces_.begin());
  }
  else
  {  
    // This is retrieving alloc warnings. It might be due to Petsc.
    // HYPRE_IJMatrixSetRowSizes(A, maxInProcFaces_.begin());
  }
 }
 
 
 label ilower = this->sharedData[meshID_].ilower;
 
 // Initialize before setting coefficients 
 HYPRE_IJMatrixInitialize(A);
 HYPRE_IJVectorInitialize(b);
 HYPRE_IJVectorInitialize(x);
 
 label nc = T.size();  
  
 //- Off diagonal elements   
 const lduAddressing& offDiag = eqn.lduMatrix::lduAddr();
 label nnz = 1; label col; label row;
 label nFaces = offDiag.lowerAddr().size();
 
 // Symmetric matrices only have upper(). No need to force creation of lower().
 if (eqn.symmetric())
 {
   for (register label face=0; face<nFaces; face++)
   {
     row = ilower + offDiag.upperAddr()[face];
     col = ilower + offDiag.lowerAddr()[face];
   
     // ij Off-diagonal
     HYPRE_IJMatrixSetValues(A, 1, &nnz, &row, &col, &eqn.upper()[face]); 
     // ji Off-diagonal
     HYPRE_IJMatrixSetValues(A, 1, &nnz, &col, &row, &eqn.upper()[face]);   
   }
 }
 else
 {
   for (register label face=0; face<nFaces; face++)
   {
     row = ilower + offDiag.upperAddr()[face];
     col = ilower + offDiag.lowerAddr()[face];
   
     // ij Off-diagonal
     HYPRE_IJMatrixSetValues(A, 1, &nnz, &row, &col, &eqn.lower()[face]); 
     // ji Off-diagonal
     HYPRE_IJMatrixSetValues(A, 1, &nnz, &col, &row, &eqn.upper()[face]);   
   }
 }
 
 //- Diagonal elements and source vector
 std::vector<double> rhs_values(nc);
 std::vector<double> x_values(nc);
 std::vector<label> rows(nc);

 scalarField source = eqn.source().component(cmpI);   
 for (label cellI=0; cellI<nc; cellI++) 
 {  
   // Diagonal
   row = ilower + cellI;  
   HYPRE_IJMatrixSetValues(A, 1, &nnz, &row, &row, &eqn.diag()[cellI]); 
   
   // Source vector   
   rhs_values[cellI] = source[cellI];
   x_values[cellI] = 0.;
   rows[cellI] = row;
 }
  
 //- Contribution from BCs
 forAll(T.boundaryField(), patchI)
 {
   scalarField bC = eqn.boundaryCoeffs()[patchI].component(cmpI);  
   scalarField iC = eqn.internalCoeffs()[patchI].component(cmpI);  
   const labelUList& addr = eqn.lduMatrix::lduAddr().patchAddr(patchI);      
   const fvPatch& pfvPatch = T.mesh().boundary()[patchI]; 
   
   // Non-coupled
   if (!T.boundaryField()[patchI].coupled())
   {
     forAll(addr, facei)
     {
       // Matrix of coefs - Diagonal
       row = ilower + addr[facei];  
       HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[facei]);
       
       // Source vector
       rhs_values[addr[facei]] += bC[facei];        
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
       HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[facei]);       
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
           HYPRE_IJMatrixSetValues
           (
             A,
             1,
             &nnz,
             &this->sharedData[meshID_].fCo[pI][facei],
             &this->sharedData[meshID_].fCn[pI][facei],
             &v
           );      
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

       scalarField pnfc = pnf.component(cmpI);
      
       forAll(ownFC, facei)
       {           
         row = ilower + ownFC[facei];  
        
         // Matrix of coefs - Diagonal   
         HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[facei]); 
         // Source vector
         rhs_values[ownFC[facei]] += bC[facei]*pnfc[facei]; 
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
          HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[facei]);   
       
          // Matrix of coefs - off-diagonal
          double v = -bC[facei];
          HYPRE_IJMatrixSetValues(A, 1, &nnz, &row, &col, &v);       
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
           HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[k]);  
         
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (camipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (camipp.AMIs()[i].srcWeightsSum()[k] < camipp.AMIs()[i].lowWeightCorrection())
             {
               col = row;
               double v = -bC[k];
               HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);  
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(srcAd[k], kk)
               {                        
                col = neiFCproc[srcAd[k][kk]];
                double v = -bC[k]*srcW[k][kk];
                HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);           
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
               HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);                    
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
           HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &row, &iC[k]);                     
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (neicamipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (neicamipp.AMIs()[i].tgtWeightsSum()[k] < neicamipp.AMIs()[i].lowWeightCorrection())
             {
               col = row;
               double v = -bC[k];
               HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);   
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(tgtAd[k], kk)
               {                
                col = neiFCproc[tgtAd[k][kk]];
                double v = -bC[k]*tgtW[k][kk];
                HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);           
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
               HYPRE_IJMatrixAddToValues(A, 1, &nnz, &row, &col, &v);            
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
     << "Sorry, but Hypre interface cannot support patches of type: "
     << pfvPatch.type() << ", which is currently set for "
     << "patch " << pfvPatch.name() << "."
     << exit(FatalError);   
   } 
      
 }
 
 //- Assemble the matrix
 HYPRE_IJMatrixAssemble(A);
 
 //- Set the source and solution vectors from the arrays. Assemble and clean.
 HYPRE_IJVectorSetValues(b, nc, &rows[0], &rhs_values[0]); 
 HYPRE_IJVectorAssemble(b);
 
 HYPRE_IJVectorSetValues(x, nc, &rows[0], &x_values[0]);
 HYPRE_IJVectorAssemble(x);
} 

template<class Type>
void Foam::hypreSolver<Type>::assembleHypreBx
(
  HYPRE_IJVector& b,
  HYPRE_IJVector& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{
 
 label ilower = this->sharedData[meshID_].ilower;
 
 // Initialize before setting coefficients 
 HYPRE_IJVectorInitialize(b);
 HYPRE_IJVectorInitialize(x);
 
 label nc = T.size();  
  
 //- Diagonal elements and source vector
 std::vector<double> rhs_values(nc);
 std::vector<double> x_values(nc);
 std::vector<label> rows(nc);

 scalarField source = eqn.source().component(cmpI);   
 for (label cellI=0; cellI<nc; cellI++) 
 {  
   // Source vector   
   rhs_values[cellI] = source[cellI];
   x_values[cellI] = 0.;
   rows[cellI] = ilower + cellI; 
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
       // Source vector
       rhs_values[addr[facei]] += bC[facei];        
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

     scalarField pnfc = pnf.component(cmpI);
      
     const labelList& ownFC = T.mesh().boundaryMesh()[patchI].faceCells(); 
     forAll(ownFC, facei)
     {           
       rhs_values[ownFC[facei]] += bC[facei]*pnfc[facei];    
     }       
   }    
 }
 
 //- Set the source and solution vectors from the arrays. Assemble and clean.
 HYPRE_IJVectorSetValues(b, nc, &rows[0], &rhs_values[0]); 
 HYPRE_IJVectorAssemble(b);
 
 HYPRE_IJVectorSetValues(x, nc, &rows[0], &x_values[0]);
 HYPRE_IJVectorAssemble(x);
} 


template<class Type>
void Foam::hypreSolver<Type>::transferHypreSolution
(
  HYPRE_IJVector& x,
  GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{
   label ilower = this->sharedData[meshID_].ilower;
   
   // get the local solution 
   label nvalues = T.size();
   std::vector<label> rows(nvalues);
        
   for (label i = 0; i < nvalues; i++)
    rows[i] = ilower + i;
    
   scalarField tt(T.size(), 0.);
   HYPRE_IJVectorGetValues(x, nvalues, &rows[0], tt.begin());   
   T.primitiveFieldRef().replace(cmpI, tt);
} 
 

template<class Type>
void Foam::hypreSolver<Type>::checkMatrixSum
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
   if ( dif > 1e-20 && !updateA_)
    {
     if (Pstream::master)
     {
      FatalErrorIn
        (
            "Foam::hypreSolver<Type>::checkMatrixSum"
        )   << nl << nl << "Sum of fvMatrix(" << name << ") is changing by " << dif
            << " between consecutive time-steps and option updateMatrixCoeffs is disabled."
            << " You should enable this option or set option saveSystem to false." << nl << nl << endl
            << exit(FatalError);
     }
    }
  }
    
}
// ************************************************************************* //
