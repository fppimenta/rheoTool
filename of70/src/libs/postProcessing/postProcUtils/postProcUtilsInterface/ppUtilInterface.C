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

#include "ppUtilInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ppUtilInterface, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ppUtilInterface::ppUtilInterface
(
    const dictionary* dict,
    const volVectorField& U
) 
:
 ppUPtr_()
{ 
  if (dict!=NULL)
  {
    PtrList<entry> l((*dict).lookup("functions"));
    ppUPtr_.setSize(l.size());

    
    forAll (ppUPtr_, i)
    {
        ppUPtr_.set
        (
            i,
            ppUtil::New(l[i].keyword(), l[i].dict(), U)
        );
    }
  }
 else
  {
    ppUPtr_.setSize(0);
  }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void ppUtilInterface::update() 
{
    forAll (ppUPtr_, i)
    {
        ppUPtr_[i].update();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
