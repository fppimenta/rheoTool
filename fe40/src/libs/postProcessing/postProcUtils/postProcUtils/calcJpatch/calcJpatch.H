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
    calcJpatch
    
Description
    This util computes and writes the dimensional surface-averaged current
    density over a list of patches (user-defined). A value is retrieved for
    each specie.
    
    J_ci = (1/|Sf_patch|) * sum[ F*zi*phi*ci - F*zi*D*{ grad(ci)&n + (e*z/k*T)*ci*grad(psi)&n } * |Sf|]
    
    In case the log-version is used, the diffusive flux is ci*grad(lci)&n,
    with lci = log(ci). 
     
    It is used in the selecMembrane tutorial of rheoTool.  
    
    This class is part of rheoTool.

SourceFiles
    calcJpatch.C

\*---------------------------------------------------------------------------*/

#ifndef calcJpatch_H
#define calcJpatch_H

#include "ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{

/*---------------------------------------------------------------------------*\
                           Class calcJpatch Declaration
\*---------------------------------------------------------------------------*/

class calcJpatch
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
        calcJpatch(const calcJpatch&);

        //- Disallow default bitwise assignment
        void operator=(const calcJpatch&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcJpatch");


    // Constructors

        //- Construct from components
        calcJpatch
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcJpatch()
        {}


    // Member Functions

        //- Correct the viscoelastic stress
        virtual void update();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ppUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
