/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::fv::gaussDefCmpwConvectionScheme

Description
    NVA convective schemes implemented in a deferred and component-wise approach,
    ie, only the upwind part of the scheme is introduced implicitly (remaining goes
    to the source term) and the limiters are individually computed for each component, 
    when the convected variable is of rank >= 1. The available schemes can be found
    at the end of source file gaussDefCmpwConvectionScheme.C.
    
    This class is part of rheoTool.

SourceFiles
    gaussDefCmpwConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussDefCmpwConvectionScheme_H
#define gaussDefCmpwConvectionScheme_H

#include "convectionScheme.H"
#include "etcFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

//- Temporary debug switch to provide warning about backward-compatibility
//  issue with setting div schemes for steady-state
extern int warnUnboundedgaussDefCmpw;

/*---------------------------------------------------------------------------*\
                       Class gaussDefCmpwConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussDefCmpwConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private data
 
        word scheme_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        gaussDefCmpwConvectionScheme(const gaussDefCmpwConvectionScheme&);

        //- Disallow default bitwise assignment
        void operator=(const gaussDefCmpwConvectionScheme&);

public:

    //- Runtime type information
    TypeName("GaussDefCmpw");


    // Constructors

        //- Construct from flux and interpolation scheme
        gaussDefCmpwConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const word& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_(scheme)
        {}

        //- Construct from flux and Istream
        gaussDefCmpwConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is,
            const word& scheme
        )
        :
            convectionScheme<Type>(mesh, faceFlux),
            scheme_(scheme)
        {
            is.rewind();
            word bounded(is);

            if
            (
                warnUnboundedgaussDefCmpw
             && word(mesh.ddtScheme("default")) == "steadyState"
             && bounded != "bounded"
            )
            {
                fileNameList controlDictFiles(findEtcFiles("controlDict"));

                IOWarningInFunction(is)
                    << "Unbounded 'gaussDefCmpw' div scheme used in "
                       "steady-state solver, use 'bounded gaussDefCmpw' "
                       "to ensure boundedness.\n"
                    << "    To remove this warning switch off "
                    << "'boundedgaussDefCmpw' in "
                    << controlDictFiles[controlDictFiles.size()-1]
                    << endl;
            }
        }

    // Member Functions


        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type> > fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > phifDefC
        (
            const surfaceScalarField&,
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&,
            const bool& 
        ) const;

       
        void lims
        (
            scalarList&, scalarList&, scalarList&,
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&    
        ) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "gaussDefCmpwConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
