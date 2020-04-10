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

#include "grad.H"
 
namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
Note: List LMatrix.rel contains the row-col-cmp indices. Ie, LMatrix.field(cmp) 
is to be inserted in (row, col) of the block matrix. Both row/col assume that the
first valid component of the fields involved is 0 (fields ordering in the block 
matrix is accounted for later inside the coupled solver). The coupled solver will
only insert the equations specified in LMatrix.rel, even though the operators 
are built for all the 3 geometric dimensions.
*/

namespace fvmb
{
 
// vector-scalar
tmp<LMatrix<vector> >  grad
(
   const dimensionedScalar c,
   const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
  return grad(c.value(), vf);
}

tmp<LMatrix<vector> >  grad
(
   const scalar c,
   const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
  volScalarField cVolField
  (
    IOobject
    (
      "c",
      vf.time().constant(),
      vf.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
    ),
    vf.mesh(),
    dimensionedScalar("cs", dimless, c)
  );
  
  return grad(cVolField, vf);    
}

// Multiplication by c is such that it reproduces the behavior
// that we would have for c*fvc::grad(vf)*V. Thus, c is not interpolated 
// on faces.
// vector-scalar
tmp<LMatrix<vector> >  grad
(
   const GeometricField<scalar, fvPatchField, volMesh>& c,
   const GeometricField<scalar, fvPatchField, volMesh>& vf
)
{
  const fvMesh& mesh = vf.mesh();
   
  tmp< LMatrix<vector> > tfvm
  (
    LMatrix<vector>::New(mesh)
  );
  LMatrix<vector>& fvm = tfvm.ref();
 
  vectorField& upper = fvm.upper();
  vectorField& lower = fvm.lower();
  vectorField& diag = fvm.diag();
  List<labelList>& rel = fvm.rowColList();
   
  const labelUList& own = mesh.owner();
  const labelUList& neig = mesh.neighbour();
  vectorField Sf = mesh.faceAreas();  
  
  surfaceScalarField weights = vf.mesh().surfaceInterpolation::weights();

  forAll(own, facei)
  {    
    scalar w = weights[facei];
    
    // For the owner[facei] cell 
    diag[own[facei]] += w*Sf[facei]*c[own[facei]];
    upper[facei] = (1. - w)*Sf[facei]*c[own[facei]];
    
    // For the neig[facei] cell 
    diag[neig[facei]] -= (1. - w)*Sf[facei]*c[neig[facei]];
    lower[facei] = -w*Sf[facei]*c[neig[facei]];
  } 
  
  forAll(vf.boundaryField(), patchi)
  {
    const fvPatchField<scalar>& pvf = vf.boundaryField()[patchi];
    const fvsPatchField<scalar>& bw = weights.boundaryField()[patchi];   
    vectorField Sfb = pvf.patch().Sf(); 
    const scalarField cpI = c.boundaryField()[patchi].patchInternalField();
            
    // bw refer to the f-C fractional distance on the other processor. The
    // cell on each processor is always the owner, so bw still corresponds 
    // to the fractional distance of the neighbour.
    if (pvf.coupled())
    {      
       // To procowner
       fvm.internalCoeffs()[patchi] =  bw * Sfb * cpI;
        
       // To procneigh
       fvm.boundaryCoeffs()[patchi] = -(1. - bw) * Sfb * cpI; // the (-) is because it would go to source
    }
    else
    {
      const scalarField internalCoeffs(pvf.valueInternalCoeffs(bw));
      const scalarField boundaryCoeffs(pvf.valueBoundaryCoeffs(bw));
      
      fvm.internalCoeffs()[patchi] = internalCoeffs * Sfb * cpI;
        
      fvm.boundaryCoeffs()[patchi] = -boundaryCoeffs * Sfb * cpI;
    }
  }
  
  // Fill the row-col-cmpt list
  typename pTraits<vector>::labelType validComponents
  (
    mesh.template validComponents<vector>()
  );
 
  int i = 0;
  for (direction cmpt=0; cmpt<pTraits<vector>::nComponents; cmpt++)
  {
    if (component(validComponents, cmpt) == -1) continue; 
    
    // 0 is for the scalar var
    rel.append({i,0,cmpt});   
    i++;
  }

  return tfvm;
    
}

} // namespace fvmb
} // namespace Foam
