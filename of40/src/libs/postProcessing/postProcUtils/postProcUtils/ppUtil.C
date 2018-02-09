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

#include "ppUtil.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ppUtil, 0);
    defineRunTimeSelectionTable(ppUtil, dictFS);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ppUtil::ppUtil
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
U_(U),
outS_(NULL),
ppDir_(word::null),
enabled_(dict.lookupOrDefault<bool>("enabled", false)),
nEval_(dict.lookupOrDefault<int>("evaluateInterval", 10000)),
counter_(1)
{
  //- Create folder with postProcessing data
   
  if (Pstream::parRun()) 
   {
     ppDir_ = U.time().path()/".."/"rheoToolPP"/U.time().timeName()/name;
   }
  else
   {
     ppDir_ = U.time().path()/"rheoToolPP"/U.time().timeName()/name;
   }  
    
   if (enabled_)
      if (Pstream::master()) 
         mkDir(ppDir_);
}


autoPtr<ppUtil> ppUtil::New
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
{
    word typeName = dict.lookup("funcType");

    Info<< "Selecting post-processing utility: " << typeName << endl;

    dictFSConstructorTable::iterator cstrIter =
        dictFSConstructorTablePtr_->find(typeName);

    if (cstrIter == dictFSConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ppUtil::New(const dictionary& dict, const volVectorField&)"
        )   << "Unknown ppUtil type " << typeName
            << endl << endl
            << "Valid ppUtil types are :" << endl
            << dictFSConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<ppUtil>(cstrIter()(name, dict, U));
}


} //End namespace

// ************************************************************************* //
