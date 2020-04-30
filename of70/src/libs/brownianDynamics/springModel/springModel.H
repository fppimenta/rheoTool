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
    springModel
    
Description
    Base spring model.
    This class is part of rheoTool.     

SourceFiles
    springModel.C

\*---------------------------------------------------------------------------*/

#ifndef springModel_H
#define springModel_H

 
#include "sPCloudInterface.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class springModel Declaration
\*---------------------------------------------------------------------------*/

class springModel
{
private:

    // Private data

        //- Reference to continuum velocity field
        const volVectorField& U_;
        
        //- Reference to the sPCloudInterface object 
        sPCloudInterface& sPCref_;
         
        
        //- Newton-Raphson parameters
        
        //-- Max number of iterations to exit
        int maxIter_;
        
        //-- Relative tolerance to exit (normalized by Ls)
        scalar relTol_;
        
        //-- Ls fraction to trigger the semi-implicit scheme (below that is remains explicit Euler)
        scalar tresholdF_;  
        
        //-- Matrix solver
        word solverType_;      

    // Private Member Functions

        //- Disallow default bitwise copy construct
        springModel(const springModel&);

        //- Disallow default bitwise assignment
        void operator=(const springModel&);
        
        //- Adds spring force contribution implicitly
        void implicitfSpring(label, const Field<vector>&, const Field<vector>&); 
  
        //- Deletes a molecule if a spring is overstretched (call to delete function of interface)       
        void deleteMolecule(label);
        
protected:

  // Protected data
  
    //- Time integration scheme
    word  timeSch_;   

    //- References from sPCloudInterface (definitions therein)   
    List<label>& linkM_;
    int& nMolc_;
         
    PtrList<Field<vector > >& mx_;       
    PtrList<Field<vector > >& mU_;              
    PtrList<List<label > >& mAct_;      
    PtrList<List<List<label > > >& mIds_;
    PtrList<List<List<symmTensor > > >& mD_;
    PtrList<List<List<tensor > > >& mSigma_;
    PtrList<List<List<label > > >& mSpr_;
	   
    const  List<scalar>& Nks_;
    const  List<scalar>& D_;
    const  List<scalar>& Ls_;
    
    const bool& isTethered_;
    const bool& isHI_;
    
 // Protected member functions   

    //- Implementation of the TDMA method
    void TDMA(scalarSquareMatrix&, scalarField&);
    
    //- Model-dependent implementation of the spring force (call per molecule)
    virtual tmp<vectorField> fSpringI(vectorField&, label, bool, bool) = 0;
        
    //- Function f of Newton-Raphson method
    virtual void fSIM
    (
      label mi,
      label cmpi, 
      const Field<scalar>& xStar, 
      const Field<vector>& x,
      Field<scalar>& fm
    ) = 0;

    //- Jacobian of f in Newton-Raphson method
    virtual void jacobianSIM
    (
      label mi,
      label cmpi,
      const Field<vector>& x,
      scalarSquareMatrix& J
    ) = 0; 

public:

    //- Runtime type information
    TypeName("springModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            springModel,
            dictFS,
            (
                const dictionary& dict,
                const volVectorField& U,
                sPCloudInterface& sPCI
            ),
            (dict, U, sPCI)
        );


    // Constructors

        //- Construct from components
        springModel
        (
            const dictionary& dict,
            const volVectorField& U,
            sPCloudInterface& sPCI
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<springModel> New
        (
            const dictionary& dict,
            const volVectorField& U,
            sPCloudInterface& sPCI
        );


    // Destructor

        virtual ~springModel()
        {}


    // Member Functions
        
        //- Return reference to continuum velocity field
        const volVectorField& U() const
        {
            return U_;
        }
        
        //- Adds contribution of the spirng force to the beads velocity and 
        // position. Auto loop over all valid molecules.
        void fSpring(); 
        
        //- For the explicit scheme, check for overstretch and delete if needed.
        //  For the semi-implicit scheme, check is tresholdF_ is exceed and call
        // the implicit force if needed. Delete the molecules if the implicit force
        // was not able to avoid overstretch.            
        virtual void checkSpringsLength
        (         
          const PtrList<Field<vector > >& mxStar,
          const PtrList<Field<vector > >& mx0
        );
	
	//- Return true if the model limits the springs length
	virtual bool isLimitedSpring() = 0;      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
