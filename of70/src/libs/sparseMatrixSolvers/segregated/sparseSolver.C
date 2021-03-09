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
 
#include "objectRegistry.H"
#include "dictionary.H"
#include <vector>
 
// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::sparseSolver<Type>> Foam::sparseSolver<Type>::New
(
    const GeometricField<Type, fvPatchField, volMesh>& T,
    const fvMesh& mesh,
    const dictionary& fvSolution
)
{
  // Using lookupOrDefault for backward compatibility
  word modelType(fvSolution.subDict("solvers").subDict(word(T.name())).lookupOrDefault<word>("solverType", "openFoamSolver"));

  Info<< "Using "<< modelType << " solver for " << word(T.name()) << endl;

  typename dictionaryConstructorTable::iterator cstrIter =
      dictionaryConstructorTablePtr_->find(modelType);

  if (cstrIter == dictionaryConstructorTablePtr_->end())
  {
    FatalErrorInFunction
     << "Unknown sparseSolver type "
     << modelType << nl << nl
     << "Valid sparseSolvers are : " << endl
     << dictionaryConstructorTablePtr_->sortedToc()
     << exit(FatalError);
  }

  return autoPtr<sparseSolver<Type> > (cstrIter()(T, mesh, fvSolution));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::sparseSolver<Type>::~sparseSolver()
{
   
}

// * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * //

template<class Type>
void Foam::sparseSolver<Type>::addBoundarySource
(
    Field<Type>& source,
    fvMatrix<Type>& matrix,
    const GeometricField<Type, fvPatchField, volMesh>& T,
    const bool couples
) const
{
  forAll(T.boundaryField(), patchi)
  {
    const fvPatchField<Type>& ptf = T.boundaryField()[patchi];
    Field<Type>& pbc = matrix.boundaryCoeffs()[patchi];

    if (!ptf.coupled())
    {
      const labelUList& addr(matrix.lduMatrix::lduAddr().patchAddr(patchi));
      forAll(addr, facei)
      {
        source[addr[facei]] += pbc[facei];
      }
    }
    else if (couples)
    {
      const tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
      const Field<Type>& pnf = tpnf();

      const labelUList& addr = matrix.lduMatrix::lduAddr().patchAddr(patchi);

      forAll(addr, facei)
      {
        source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
      }
    }
 }
}

template<class Type>
void Foam::sparseSolver<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt,
    fvMatrix<Type>& matrix
) const
{
  forAll(matrix.internalCoeffs(), patchi)
  {
    const labelUList& addr(matrix.lduMatrix::lduAddr().patchAddr(patchi));
    scalarField iF(matrix.internalCoeffs()[patchi].component(solveCmpt));
        
    forAll(addr, facei)
    {
      diag[addr[facei]] += iF[facei];
    }
  }
}
 
template<class Type>
Foam::scalar Foam::sparseSolver<Type>::getFoamResiduals
(
    const GeometricField<Type, fvPatchField, volMesh>& T,
    fvMatrix<Type>& matrix,
    const scalarField& sourceCmpt,
    const scalarField& psiCmpt, 
    const scalarField& saveDiag,
    const FieldField<Field, scalar>& bouCoeffsCmpt,
    const lduInterfaceFieldPtrsList& interfaces,
    const int nEvalInit,
    const bool saveSystem,
    const direction cmpt,
    const int vcmpt   
) 
{   
   // We need the diagonal to be complete for the computation
   // of A.psi
   addBoundaryDiag(matrix.diag(), cmpt, matrix);
   
   scalarField wA(T.size());
   scalarField pA(T.size());
   
   // --- Calculate A.psi
   matrix.lduMatrix::Amul(wA, psiCmpt, bouCoeffsCmpt, interfaces, cmpt);

   // --- Calculate initial residual field
   scalarField rA(sourceCmpt - wA);

   // --- Calculate normalisation factor
   // ---- Calculate A dot reference value of psi
   matrix.lduMatrix::sumA(pA, bouCoeffsCmpt, interfaces);
   
   // Check matrix constness
   if (T.time().timeIndex() < nEvalInit + 1 && saveSystem)
     checkMatrixSum(pA, T.name(), T.time().timeIndex(),  vcmpt);  
 
   pA *= gAverage(psiCmpt, T.mesh().comm());

   scalar normF = 
   gSum
   (
      (mag(wA - pA) + mag(sourceCmpt - pA))(),
      T.mesh().comm()
   )
   + 1e-20;
   
   // Restore diagonal
   matrix.diag() = saveDiag; 
   
   return  gSumMag(rA, T.mesh().comm())/normF;
} 

template<class Type>
void Foam::sparseSolver<Type>::checkMatrixSum
(
    const scalarField& rowSum,
    const word name,
    const int tindex,
    const int vcmpt   
) 
{
  // Re-implemented in derived
}

template<class Type>
void Foam::sparseSolver<Type>::solve
(
  fvMatrix<Type>& matrix,
  const dictionary&
)
{
  solve(matrix);
}

template<class Type>
void Foam::sparseSolver<Type>::solve
(
  const tmp<fvMatrix<Type>>& tmatrix
)
{
  solve(const_cast<fvMatrix<Type>&>(tmatrix()));
  
  tmatrix.clear();
}

template<class Type>
void Foam::sparseSolver<Type>::solve
(
  const tmp<fvMatrix<Type>>& tmatrix,
  const dictionary& dict
)
{
  solve(const_cast<fvMatrix<Type>&>(tmatrix()));
  
  tmatrix.clear();
}

// ************************************************************************* //
