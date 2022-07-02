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

#include "twoSymmGrad.H"
 
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
 
// symmTensor-vector
// Operator: c*[grad(U) + T(grad(U))]
tmp< LMatrix<tensor> >  twoSymmGrad
(
   const dimensionedScalar c,
   const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
  return twoSymmGrad(c.value(), vf);
}


tmp< LMatrix<tensor> >  twoSymmGrad
(
   const scalar c,
   const GeometricField<vector, fvPatchField, volMesh>& vf
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
  
  return twoSymmGrad(cVolField, vf);    
}

tmp< LMatrix<tensor> >  twoSymmGrad
(
   const GeometricField<scalar, fvPatchField, volMesh>& c,
   const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
  const fvMesh& mesh = vf.mesh();
  
  tmp< LMatrix<tensor> > tfvm
  (
    LMatrix<tensor>::New(mesh)
  );
  LMatrix<tensor>& fvm = tfvm.ref();
 
  tensorField& upper = fvm.upper();
  tensorField& lower = fvm.lower();
  tensorField& diag = fvm.diag();
  List<labelList>& rel = fvm.rowColList();
   
  const labelUList& own = mesh.owner();
  const labelUList& neig = mesh.neighbour();
  const vectorField Sf = mesh.faceAreas(); 
  
  surfaceScalarField weights = vf.mesh().surfaceInterpolation::weights();
  
  // Here, the diag components of tensor directly account for the sum of gradA + gradA.T()
  // with the multiplication by 2 of these elements (we can do it because the diagonal elements
  // of both gradA and gradA.T() are the same, and involve the same pair SiTj). Off-diagonal elements are
  // computed as if it was for gradA alone and the sum gradA + gradA.T() for these coefficients
  // is accounted for later via adressing in fvm.rowColList(), since both elements ij and ji
  // contribute to ij, but with different coefficients each.  
  forAll(own, facei)
  {    
    scalar w = weights[facei];
    
    tensor tt = tensor(vector::one*Sf[facei].x(),vector::one*Sf[facei].y(),vector::one*Sf[facei].z());
    tt.xx() *= 2.;
    tt.yy() *= 2.;
    tt.zz() *= 2.;
    tensor wo = w*tt;
    tensor wn = (1.-w)*tt;
    
    // For owner
    diag[own[facei]] += wo*c[own[facei]];
    upper[facei] = wn*c[own[facei]];
    
    // For neig
    diag[neig[facei]] -= wn*c[neig[facei]];  
    lower[facei] = -wo*c[neig[facei]];
  } 
  
  forAll(vf.boundaryField(), patchi)
  {
    const fvPatchField<vector>& pvf = vf.boundaryField()[patchi];
    const fvsPatchField<scalar>& bw = weights.boundaryField()[patchi];
    const vectorField& Sfb = vf.mesh().Sf().boundaryField()[patchi];
    tensorField& internalCoeffsFvm(fvm.internalCoeffs()[patchi]);
    tensorField& boundaryCoeffsFvm(fvm.boundaryCoeffs()[patchi]);
    const scalarField cpI(c.boundaryField()[patchi].patchInternalField());
        
    // bw refer to the f-C fractional distance on the other processor. The
    // cell on each processor is always the owner, so bw still corresponds 
    // to the fractional distance of the neighbour.
    if (pvf.coupled())
    {         
      forAll(pvf, facei)
      {
        tensor SfbTensor = tensor
        (
          vector::one*Sfb[facei].x(),
          vector::one*Sfb[facei].y(),
          vector::one*Sfb[facei].z()
        );
        SfbTensor.xx() *= 2.;
        SfbTensor.yy() *= 2.;
        SfbTensor.zz() *= 2.;
        
        // To procowner
        internalCoeffsFvm[facei] = bw[facei] * SfbTensor * cpI[facei]; 
        
        // To procneigh
        boundaryCoeffsFvm[facei] = -(1. - bw[facei]) * SfbTensor * cpI[facei]; // the (-) is because it would go to source
      }
    }
    else
    {
      const vectorField internalCoeffs(pvf.valueInternalCoeffs(bw));
      const vectorField boundaryCoeffs(pvf.valueBoundaryCoeffs(bw));
     
      forAll(pvf, facei)
      {
        internalCoeffsFvm[facei] = Sfb[facei] * internalCoeffs[facei] * cpI[facei];
        boundaryCoeffsFvm[facei] = -Sfb[facei] * boundaryCoeffs[facei] * cpI[facei];
        internalCoeffsFvm[facei].xx() *= 2.; internalCoeffsFvm[facei].yy() *= 2.; internalCoeffsFvm[facei].zz() *= 2.;
        boundaryCoeffsFvm[facei].xx() *= 2.; boundaryCoeffsFvm[facei].yy() *= 2.; boundaryCoeffsFvm[facei].zz() *= 2.;
      }          
    }
  }
  
  // Fill the list
  
  //-- Tensor
  typename pTraits<symmTensor>::labelType validComponentsT
  (
    mesh.template validComponents<symmTensor>()
  );
  
  int i = 0;
  labelList cT(6,-1);
  forAll(cT, cmpt)
  {
    if (component(validComponentsT, cmpt) == -1) continue; 
     
    cT[cmpt]= i;
    i++;
  }
  
  //-- Vector
  typename pTraits<vector>::labelType validComponentsV
  (
    mesh.template validComponents<vector>()
  );
  
  i = 0;
  labelList cV(3,-1);
  forAll(cV, cmpt)
  {
    if (component(validComponentsV, cmpt) == -1) continue; 
    
    cV[cmpt] = i;     
    i++;
  }
 
  // Fill the row-col-cmpt list
  List<labelList> VT({
  {cV[0],cT[0]}, {cV[1],cT[1]}, {cV[2],cT[2]},
  {cV[0],cT[1]}, {cV[1],cT[3]}, {cV[2],cT[4]},
  {cV[0],cT[2]}, {cV[1],cT[4]}, {cV[2],cT[5]}
  });
   
  forAll(VT, i)
  {
    if (VT[i][0] == -1 || VT[i][1] == -1 ) continue;
    
    // Eq for tensor, contrib from vector
    rel.append( {VT[i][1], VT[i][0], i} );
  }
   
  return tfvm;    
}

tmp< LMatrix<tensor> >  twoSymmGrad
(
   const tmp<GeometricField<scalar, fvPatchField, volMesh>>& tc,
   const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
  tmp< LMatrix<tensor> >  twoSymmGradOp(twoSymmGrad(tc(), vf));
  tc.clear();
  return twoSymmGradOp;
}

} // namespace fvmb
} // namespace Foam
