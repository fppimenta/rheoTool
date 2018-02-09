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

Class
    calcBalance
    
Description
    This util computes and writes the total mass (divided by the 
    domain's volume) and sums-up the flux through all the domain boundaries, 
    for each specie.
    
    This class is part of rheoTool.

SourceFiles
    calcBalance.C

\*---------------------------------------------------------------------------*/

#ifndef calcBalance_H
#define calcBalance_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcBalance Declaration
\*---------------------------------------------------------------------------*/

class calcBalance
:
    public ppUtil
{
    // Private data

        List<scalar> D_;
	List<scalar> z_;
	List<word> specL_;	 
	
	scalar T_;
	
	List<word> patchL_;
	 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcBalance(const calcBalance&);

        //- Disallow default bitwise assignment
        void operator=(const calcBalance&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcBalance");


    // Constructors

        //- Construct from components
        calcBalance
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcBalance()
        {}


    // Member Functions

        //- Correct the viscoelastic stress
        virtual void update();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
