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
    ppUtil
    
Description
    This is the abstract base class for the post-processing utilities used in
    the tutorials provided with rheoTool for foam-extend versions. Those are 
    used to compute and extract some usefull results for the set of tutorials 
    included in rheoTool.
    
    The controls for this class are expected to be defined in fvSolution, in a 
    subDict named 'postProcessing'. Although the controls of these utils where 
    moved to the solver itself, here is an example of such a subDict:
    
    PostProcessing
    {
      enabled            true; // Is the util enabled? (true/false)
      evaluateInterval   20; // When is supposed to evaluate the util? (# of time-steps) 
      funcType           calcKineticE; // TypeName of the util
    }
    
    In order to not use any of the utils set 'enabled' as false and 'funcType' as none. 
    
    This class is part of rheoTool.

SourceFiles
    ppUtil.C

\*---------------------------------------------------------------------------*/

#ifndef ppUtil_H
#define ppUtil_H

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class ppUtil Declaration
\*---------------------------------------------------------------------------*/

class ppUtil
{
private:

    // Private data


        //- Reference to velocity field
        const volVectorField& U_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ppUtil(const ppUtil&);

        //- Disallow default bitwise assignment
        void operator=(const ppUtil&);
        
protected:

    // Protected Member Functions
    
      // OStream
      autoPtr<OFstream> outS_;

      // Is it enabled
      bool enabled_;
      
      // Evaluate frequency
      int nEval_; 

      // Counter
      int counter_;

public:

    //- Runtime type information
    TypeName("ppUtil");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ppUtil,
            dictFS,
            (
                const dictionary& dict,
                const volVectorField& U
            ),
            (dict, U)
        );


    // Constructors

        //- Construct from components
        ppUtil
        (
            const dictionary& dict,
            const volVectorField& U
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<ppUtil> New
        (
            const dictionary& dict,
            const volVectorField& U
        );


    // Destructor

        virtual ~ppUtil()
        {}


    // Member Functions

      
        //- Return the velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Correct the viscoelastic stress
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
