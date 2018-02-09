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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void gaussDefCmpwConvectionScheme<Type>::lims
(
   scalarField& alpha,
   scalarField& beta,
   scalarField& bounds,
   const surfaceScalarField& faceFlux,
   const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (scheme_ == "upwind")
     {
       alpha[0]=1.; 
       beta[0]=0.; 
       bounds[0]=1.; 
     } 
     else if (scheme_ == "cubista")
     {
       alpha[0]=7./4.;   alpha[1]=3./4.;   alpha[2]=1./4.; 
       beta[0]=0.;       beta[1]=3./8.;    beta[2]=3./4.; 
       bounds[0]=3./8.;  bounds[1]=3./4.;
     } 
     else if (scheme_ == "minmod")
     {
       alpha[0]=1.5;     alpha[1]=.5;      alpha[2]=.5; 
       beta[0]=0.;       beta[1]=.5;       beta[2]=.5;   
       bounds[0]=.5;     bounds[1]=1.; 
     }  
     else if (scheme_ == "smart")
     {
       alpha[0]=3.;      alpha[1]=3./4.;   alpha[2]=0.; 
       beta[0]=0.;       beta[1]=3./8.;    beta[2]=1.; 
       bounds[0]=1./6.;  bounds[1]=5./6.;
     } 
     else if (scheme_ == "waceb")
     {
       alpha[0]=2.;      alpha[1]=3./4.;   alpha[2]=0.; 
       beta[0]=0.;       beta[1]=3./8.;    beta[2]=1.; 
       bounds[0]=3./10.; bounds[1]=5./6.;
     } 
     else if (scheme_ == "superbee")
     {
       alpha[0]=0.5;     alpha[1]=1.5;     alpha[2]=0.; 
       beta[0]=0.5;      beta[1]=0.;       beta[2]=1.; 
       bounds[0]=1./2.;  bounds[1]=2./3.;
     } 
     else if (scheme_ == "none")
     {
       // Do nothing
     } 
     else
     {
         FatalErrorIn("gaussDefCmpwConvectionScheme<Type>::lims()\n")
            << "\nError in div("<< faceFlux.name() << ", " << vf.name() << ")\n"
            << "\nThe deferred limited scheme is not specified or does not exist.\n"
            << "\nValid schemes are:\n"
            << "\n. upwind" <<"\n. cubista" << "\n. minmod" << "\n. smart"
            << "\n. waceb" << "\n. superbee" << "\n. none (no convection)" 
            << abort(FatalError);
 
     } 

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
