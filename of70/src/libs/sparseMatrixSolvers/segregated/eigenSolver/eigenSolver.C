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

#include "eigenSolver.H"
#include "cyclicFvPatch.H"
#include "cyclicAMIFvPatch.H"   
#include "cyclicFvPatchField.H" 
#include <chrono>

// * * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * //

// This is the time index until which we force the preconditioner to be updated
// when saveSystem is enabled.  
template<class Type>
int Foam::eigenSolver<Type>::nEvalInit_ = 3;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::eigenSolver<Type>::eigenSolver
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
initTimeIndex(mesh.time().timeIndex()),
autoPrecond(false)
{  
  // Verify limitations of the interface
  checkLimitations(T);
   
  // Detect auto mode for update of preconditioner
  if (saveSystem_ && updateA_ && updatePrecondFreq_ < 1)
    autoPrecond = true;
 
  if (this->saveSystem_)
  { 
     
    typename pTraits<Type>::labelType validComponents
    (
      mesh.template validComponents<Type>()
    );

    label nc = T.size();
    
    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
      if (component(validComponents, cmpt) == -1) continue; 
      
      listA_.append(new spmat(nc,nc));
      listSolvers_.append(EigenIterDirSolver::New(fvSolution.subDict("solvers").subDict(word(T.name()))));
      nValidCmp_++; 
    }
    
    listb_.append(new VectorXd(nc));
    listx_.append(new VectorXd(nc)); 
    
    forAll(listSolvers_, i)
    {
      // Initialize options
      listSolvers_[i]->initialize();
      
      // Check options compatibility for direct solver.
      // Direct solvers need to factorize all times the matrix coeffs change.
      if (!listSolvers_[i]->isIterative())
      {
        if (updateA_ && updatePrecondFreq_ != 1)
        {
           FatalErrorIn
           (
             "Foam::eigenSolver<Type>::eigenSolver"
           )   << nl << nl << "You cannot freeze the sparse matrix solver for field " 
           <<  T.name() << " while the matrix is being"
           << " updated, for a direct solver. Either set 'updatePrecondFrequency' to 1 or"
           << " choose an iterative matrix solver." << nl << nl << endl
           << exit(FatalError);
        }
      }
      
    }  
  }
 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::eigenSolver<Type>::checkLimitations
(
  const GeometricField<Type, fvPatchField, volMesh>& T
)
{

 // Check-abort if parallel run
 if (Pstream::parRun())
 {
   if (Pstream::master())
   {
       FatalErrorInFunction
           << "Eigen solvers cannot be run in parallel."
           << "Either run in serial or choose other solver for field: " << T.name() << nl << nl 
           << exit(FatalError);
   }
 }

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
           << "No support in Eigen interface for patches of type cyclicAMI rotational, which is the type"
           << " of patch " << pfvPatch.name() << "."
           << " Problem raised while solving for field: " << T.name() 
           << exit(FatalError);
       }      
     }
     else
     {
       FatalErrorInFunction
       << "Sorry, but Eigen interface cannot support patches of type: "
       << pfvPatch.type() << ", which is currently set for "
       << "patch " << pfvPatch.name() << "."
       << " Problem raised while solving for field: " << T.name() 
       << exit(FatalError);  
     }     
   }
 }


}

template<class Type>
void Foam::eigenSolver<Type>::solve
(
    fvMatrix<Type>& matrix
)
{
  // The system cannot be saved for a mesh changing its topo. Abort.
  if (matrix.psi().mesh().topoChanging() && this->saveSystem_)
  {      
    FatalErrorIn
    (
      "Foam::eigenSolver<Type>::solve"
    )   << nl << nl << "fvMatrix(" << matrix.psi().name() << ") is set to be saved "
        << "over time, but the mesh topology is changing. Set saveSystem option to false"
        << " in dictionary fvSolution.solvers(" << matrix.psi().name() <<")." << nl << nl << endl
        << exit(FatalError);  
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
void Foam::eigenSolver<Type>::solveReuse
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
  
   // Note: here we allow matrix update and preconditioner update
   // to happen at possibly different times. The two operations
   // are allowed to be decoupled. This works with Eigen 3.2.x but
   // will retrieve error for Eigen 3.3.x, because matrix A becomes
   // invalid upon changing its values if the preconditioner is not
   // updated before calling solve().
    
   // Now we start Eigen related stuff
   label timeID = T.time().timeIndex()- initTimeIndex;
      
   if (timeID < this->nEvalInit_ || updateA_)
   {
     assembleEigenAbx
     (
       listA_[i],
       listb_[0],
       listx_[0],
       matrix,
       T,
       cmpI 
     ); 
   }
   else
   {
     assembleEigenBx
     (
       listb_[0],
       listx_[0],
       matrix,
       T,
       cmpI 
     );
   }
   auto start = std::chrono::high_resolution_clock::now();

   // Now setup 
   if (timeID < this->nEvalInit_ || times_ > updatePrecondFreq_)
   {
    listSolvers_[i]->setup(listA_[i], matrix.symmetric(), false);    

    // Reset counter
    if (nValidCmp_ == i+1)
     times_ = 1;   
   }
   
   //... and solve                                          
   listSolvers_[i]->solve(listb_[0], listx_[0]); 
   
   auto elapsed = std::chrono::high_resolution_clock::now() - start;
   scalar solveTime = scalar(std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count());
  
   // Adjust precond update frequency if auto mode
  
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
   transferEigenSolution(listx_[0], T, cmpI);
    
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
   listSolvers_[i]->printInfo
   (
     cmpName,
     initResidual,
     finalResidual
   );
      
   i++;   
 }
 
 T.correctBoundaryConditions(); 
}

template<class Type>
void Foam::eigenSolver<Type>::solveNoReuse
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
 
 // Create and Initialize the solver outside the cmp loop
 // because we can reuse the sparsity pattern
 autoPtr<EigenIterDirSolver> solver(EigenIterDirSolver::New(this->solDict_.subDict(word(T.name())))); 
 solver->initialize(); 
 
 // Matrix of coefficients
 spmat A(T.size(), T.size());
 
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
  
   // Now we start Eigen related stuff
   VectorXd b(T.size());
   VectorXd x(T.size());
  
   assembleEigenAbx
   (
     A,
     b,
     x,
     matrix,
     T,
     cmpI 
   );  
   
   // Now setup. Analyze pattern for first cmp,
   // but reuse it for remaining components.
   if (i>0)
   {
     solver->setup(A, matrix.symmetric(), true); 
   }
   else
   {
     solver->setup(A, matrix.symmetric(), false);
   }   

   //... and solve                                          
   solver->solve(b, x);    

   // Transfer solution to OF 
   transferEigenSolution(x, T, cmpI);
    
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
   solver->printInfo
   (
     cmpName,
     initResidual,
     finalResidual
   );
      
   i++;   
 }
 
 T.correctBoundaryConditions(); 
}
 

template<class Type>
void Foam::eigenSolver<Type>::assembleEigenAbx
(
  spmat& A,
  VectorXd& b,
  VectorXd& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{
 label nc = T.size();  
  
 std::vector<trip> tripList;
  
 //- Off diagonal elements   
 const lduAddressing& offDiag = eqn.lduMatrix::lduAddr();
 label nIFaces = offDiag.lowerAddr().size();
 
 tripList.reserve(nc+nIFaces+T.mesh().nFaces());
 
 // Symmetric matrices only have upper(). No need to force creation of lower().
 if (eqn.symmetric())
 {
   for (register label face=0; face<nIFaces; face++)
   { 
    // ij Off-diagonal 
    tripList.push_back( trip(offDiag.upperAddr()[face], offDiag.lowerAddr()[face], eqn.upper()[face]) );
   
    // ji Off-diagonal
    tripList.push_back( trip(offDiag.lowerAddr()[face], offDiag.upperAddr()[face], eqn.upper()[face]) );
   }
 }
 else
 { 
   for (register label face=0; face<nIFaces; face++)
   { 
    // ij Off-diagonal
    tripList.push_back( trip(offDiag.upperAddr()[face], offDiag.lowerAddr()[face], eqn.lower()[face]) );
    
    // ji Off-diagonal
    tripList.push_back( trip(offDiag.lowerAddr()[face], offDiag.upperAddr()[face], eqn.upper()[face]) );
   }
 }
 
 //- Diagonal elements and source vector
 scalarField source = eqn.source().component(cmpI);   
 for (label cellI=0; cellI<nc; cellI++) 
 {  
   // Diagonal 
   tripList.push_back( trip(cellI, cellI, eqn.diag()[cellI]) );
   
   // Source vector   
   b(cellI) = source[cellI];
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
       tripList.push_back( trip(addr[facei], addr[facei], iC[facei]) );
       
       // Source vector
       b(addr[facei]) += bC[facei];        
     }
   }
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
         // Matrix of coefs - Diagonal  
         tripList.push_back( trip(ownFC[facei], ownFC[facei], iC[facei]) );  
         
         // Source vector
         b(ownFC[facei]) += bC[facei]*pnfc[facei];  
       }       
      }
      //... Otherwise, the coupling between cyclic patches is fully implicit and there are only contributions 
      // to the matrix of coefficients.
      else
      {
        forAll(ownFC, facei)
        {                  
          // Matrix of coefs - Diagonal
          tripList.push_back( trip(ownFC[facei], ownFC[facei], iC[facei]) );
       
          // Matrix of coefs - off-diagonal
          tripList.push_back( trip(ownFC[facei], nbFC[facei], -bC[facei]) );      
        } 
      }
     
   }
   else if (isType<cyclicAMIFvPatch>(pfvPatch))
   {           
     // camipp and neicampi always lay in the same processor, but one of them can be empty
     // if the other half of a given AMI was sent to another processor. In pratice, nothing
     // happens bellow if camipp is zero-sized.
     const cyclicAMIPolyPatch& camipp = refCast<const cyclicAMIFvPatch>(pfvPatch).cyclicAMIPatch();
     const cyclicAMIPolyPatch& neicamipp = camipp.neighbPatch();
     
     const labelList& ownFC = camipp.faceCells();   
     const labelList& neiFC = neicamipp.faceCells();
     
     labelList neiFCproc(neiFC);
     
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
           tripList.push_back( trip(ownFC[k], ownFC[k], iC[k]) ); 
         
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (camipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (camipp.AMIs()[i].srcWeightsSum()[k] < camipp.AMIs()[i].lowWeightCorrection())
             { 
               tripList.push_back( trip(ownFC[k], ownFC[k], -bC[k]) );
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(srcAd[k], kk)
               {                        
                 tripList.push_back( trip(ownFC[k], neiFCproc[srcAd[k][kk]], -bC[k]*srcW[k][kk]) );         
               }
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(srcAd[k], kk)
              {                        
                tripList.push_back( trip(ownFC[k], neiFCproc[srcAd[k][kk]], -bC[k]*srcW[k][kk]) );          
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
           tripList.push_back( trip(ownFC[k], ownFC[k], iC[k]) );            
                     
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (neicamipp.AMIs()[i].applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (neicamipp.AMIs()[i].tgtWeightsSum()[k] < neicamipp.AMIs()[i].lowWeightCorrection())
             {
                tripList.push_back( trip(ownFC[k], ownFC[k], -bC[k]) );  
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(tgtAd[k], kk)
               {                
                 tripList.push_back( trip(ownFC[k], neiFCproc[tgtAd[k][kk]], -bC[k]*tgtW[k][kk]) );          
               }  
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(tgtAd[k], kk)
              {                
                tripList.push_back( trip(ownFC[k], neiFCproc[tgtAd[k][kk]], -bC[k]*tgtW[k][kk]) );          
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
 
 // Assemble matrix
 A.setFromTriplets(tripList.begin(), tripList.end());
 
 // Simply reset x to 0
 x *= 0.;
 
} 

template<class Type>
void Foam::eigenSolver<Type>::assembleEigenBx
(
  VectorXd& b,
  VectorXd& x,
  fvMatrix<Type>& eqn,
  const GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{
 label nc = T.size();  
 
 //- Diagonal elements and source vector
 scalarField source = eqn.source().component(cmpI);   
 for (label cellI=0; cellI<nc; cellI++) 
 { 
   // Source vector   
   b(cellI) = source[cellI];
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
       b(addr[facei]) += bC[facei];        
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
        b(ownFC[facei]) += bC[facei]*pnfc[facei];
     }       
   } 
 }
 
 // Simply reset x to 0
 x *= 0.;
} 


template<class Type>
void Foam::eigenSolver<Type>::transferEigenSolution
(
  VectorXd& x,
  GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI
)
{  
   label nc = T.size();
   scalarField tt(nc, 0.);
         
   for (label i = 0; i < nc; i++)
     tt[i] = x(i);
    
   T.primitiveFieldRef().replace(cmpI, tt);
} 
 

template<class Type>
void Foam::eigenSolver<Type>::checkMatrixSum
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
      FatalErrorIn
        (
            "Foam::eigenSolver<Type>::checkMatrixSum"
        )   << nl << "Sum of fvMatrix(" << name << ") is changing by " << dif
            << " between consecutive time-steps and option updateMatrixCoeffs is disabled."
            << " You should enable this option or set option saveSystem to false." << endl
            << exit(FatalError);
    }
  }
    
}

// ************************************************************************* //
