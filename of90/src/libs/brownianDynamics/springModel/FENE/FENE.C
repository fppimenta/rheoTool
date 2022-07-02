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

#include "FENE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace springModels
{
    defineTypeNameAndDebug(FENE, 0);
    addToRunTimeSelectionTable(springModel, FENE, dictFS);
}
}
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::springModels::FENE::FENE
(
    const dictionary& dict,
    const volVectorField& U,
    sPCloudInterface& sPCI
)
:
springModel(dict, U, sPCI),
MAXFRACL_(dict.subDict("springModelProperties").lookupOrDefault<scalar>("cutOffStretch", .999))   
{  
 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 
Foam::tmp<Foam::vectorField> Foam::springModels::FENE::fSpringI
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
     scalar lrj( min(mag(rj),MAXFRACL_) );
     
     F =  
        fu
       *(
          1./(1.-lrj*lrj) 
        )
       * rj;
       
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


void Foam::springModels::FENE::fSIM
(
  label mi,
  label cmpi, 
  const Field<scalar>& xStar, 
  const Field<vector>& x,
  Field<scalar>& fm
)
{
 
 label  gI(mIds_[mi][0][2]); 
 scalar dt(U().mesh().time().deltaTValue());
 scalar f(3. * dt * D_[gI] * Nks_[gI] / (Ls_[gI] * Ls_[gI]) ); // dimless
  
 // x and x0 terms
 forAll(x, bi)
  {
    fm[bi] = x[bi].component(cmpi) - xStar[bi];
  }
 
 // Spring terms
 scalar F(0.);
 forAll(mSpr_[mi], bi)
  {
     label& b0 = mSpr_[mi][bi][0];
     label& b1 = mSpr_[mi][bi][1];
     
     scalar rj( x[b0].component(cmpi) - x[b1].component(cmpi) );
     scalar mrj( mag(x[b0] - x[b1]) );
     scalar lmrj( min(mrj,MAXFRACL_) );
     
     F =  -f * ( 1./(1.-lmrj*lmrj) ) * rj; 
     
     fm[b1]  += F;
     fm[b0]  -= F;   
  }
   
} 
 
void Foam::springModels::FENE::jacobianSIM
(
  label mi,
  label cmpi,
  const Field<vector>& x,
  scalarSquareMatrix& J
)    
{  

  label  gI(mIds_[mi][0][2]); 
  scalar dt(U().mesh().time().deltaTValue());
  scalar f(3. * dt * D_[gI] * Nks_[gI] / (Ls_[gI] * Ls_[gI]) ); // dimless
   
  // d(xa)/dxa = 1 goes to diagonal
  forAll(x, bi)
  {
    J[bi][bi] = 1.;
  }
   
  // Spring terms
  scalar Rab(0.);
  forAll(mSpr_[mi], bi)
  {
     label& b0 = mSpr_[mi][bi][0];
     label& b1 = mSpr_[mi][bi][1];
        
     scalar mab( mag(x[b0] - x[b1])  );
     scalar lmab(min(mab,MAXFRACL_));
     scalar rab( Foam::pow(x[b0].component(cmpi) - x[b1].component(cmpi), 2) );
     
     Rab =
        f 
      * ( 
           1./(lmab*lmab - 1.)
         - 2.*rab/Foam::pow( (lmab*lmab - 1.), 2)          
        ); 
      
     // Assign to matrix
      
     J[b1][b0] = Rab;
     J[b0][b1] = Rab;
     J[b1][b1] -= Rab;
     J[b0][b0] -= Rab;
  }
 
} 
// ************************************************************************* //
