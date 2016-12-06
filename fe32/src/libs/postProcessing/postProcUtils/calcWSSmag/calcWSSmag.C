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

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"

#include "calcWSSmag.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcWSSmag, 0);
    addToRunTimeSelectionTable(ppUtil, calcWSSmag, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcWSSmag::calcWSSmag
(
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(dict, U),
WSSmag_
    (
        IOobject
        (
            "WSSmag",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        (
          "zero",
          dimensionSet( 1, -1, -2, 0, 0, 0, 0),
          0.
        )
    )  
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcWSSmag::update()
{

 if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {
	   
	if (WSSmag_.time().outputTime()) // Only spends time when it is outputed
         {
         
           const volScalarField& eta_ = U().mesh().lookupObject<volScalarField>("eta");
           
          // Compute WSSmag on boundaries
 
           volTensorField L = fvc::grad(U());
           volTensorField tauNewt = eta_ * ( L + L.T() ) ;

           forAll(WSSmag_.boundaryField(), patchI)
	    {
		WSSmag_.boundaryField()[patchI] =
                mag(
			 (
			    U().mesh().Sf().boundaryField()[patchI]
			   /U().mesh().magSf().boundaryField()[patchI]
			 ) & tauNewt.boundaryField()[patchI]
                   );
	    }
	    
          } 
          
         counter_ = 0;
    }

 }
  
}


// ************************************************************************* //
