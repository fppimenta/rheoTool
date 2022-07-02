/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    
\*---------------------------------------------------------------------------*/

#include "Hookean.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace springModels
{
    defineTypeNameAndDebug(Hookean, 0);
    addToRunTimeSelectionTable(springModel, Hookean, dictFS);
}
}
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::springModels::Hookean::Hookean
(
    const dictionary& dict,
    const volVectorField& U,
    sPCloudInterface& sPCI
)
:
springModel(dict, U, sPCI)
{  
 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 
Foam::tmp<Foam::vectorField> Foam::springModels::Hookean::fSpringI
(
  vectorField& p, 
  label mi,
  bool hasDimensions,
  bool diagonalAlso
)
{

 label gI(mIds_[mi][0][2]);  
 label nb(p.size());
 scalar Ls(Ls_[gI]);
 
 // LsNorm is only used to normalize the spring lengths. If the
 // positions p are already dimensionless, then LsNorm = 1;
 scalar LsNorm(hasDimensions == true ? Ls : 1.);
  
 tmp<Field<vector> > tmUI
  (
    new Field<vector>
     (
       nb,
       vector::zero
     )
  );
 Field<vector>& mUI = tmUI.ref();
 
 //- First compute the force 
 
 scalar fu(3. * Nks_[gI]); 
  
 vector F(vector::zero);
  
 forAll(mSpr_[mi], bi)
  {
     label& b0 = mSpr_[mi][bi][0];
     label& b1 = mSpr_[mi][bi][1];
     
     vector rj( (p[b0] - p[b1])/LsNorm );
    
     F = fu * rj;
       
     mUI[b1] += F;
     mUI[b0] -= F;   
  }
  
 // Remove force from first bead if tethered
 if (isTethered_) 
  {
     mUI[0] *= 0.;   
  }

 //- ... then dot it with D 
 
 if (isHI_)
 {
   tmp<Field<vector> > tUout
   (
     new Field<vector>
     (
       nb,
       vector::zero
     )
   );
   Field<vector>& Uout = tUout.ref();
 
   forAll(mUI, i)
    {
      vector U(vector::zero);
    
      forAll(mUI, j)
       {
          U += (mD_[mi][i][j] & mUI[j])/Ls; 
       } 
    
      // Remove diagonal ii if not needed 
      if (!diagonalAlso)
       {
          U -= (mD_[mi][i][i] & mUI[i])/Ls;
       }
     
      Uout[i] = U; 
    } 
    
   return tUout;
 }
 else
 {
   mUI *= D_[gI]/Ls;
   return tmUI;   
 }
 
}

void Foam::springModels::Hookean::checkSpringsLength
(
  const Foam::PtrList<Field<vector > >& mxStar,
  const Foam::PtrList<Field<vector > >& mx0
) 
{ 
  // Only explicit Euler scheme available 
  if (timeSch_ == "semiImplicit")  
  { 
    FatalErrorIn("Foam::springModels::Hookean::checkSpringsLength()")
    << "\nOnly the 'explicit' time integration scheme is available for the Hookean spring model.\n"  
    << exit(FatalError); 
  } 
}



void Foam::springModels::Hookean::fSIM
(
  label mi,
  label cmpi, 
  const Field<scalar>& xStar, 
  const Field<vector>& x,
  Field<scalar>& fm
)
{
  FatalErrorIn("Foam::springModels::Hookean::fSIM()")
  << "\nOnly the 'explicit' time integration scheme is available for the Hookean spring model.\n"  
  << exit(FatalError); 
} 

void Foam::springModels::Hookean::jacobianSIM
(
  label mi,
  label cmpi,
  const Field<vector>& x,
  scalarSquareMatrix& J
)    
{  
  FatalErrorIn("Foam::springModels::Hookean::jacobianSIM()")
  << "\nOnly the 'explicit' time integration scheme is available for the Hookean spring model.\n"  
  << exit(FatalError);
} 
// ************************************************************************* //
