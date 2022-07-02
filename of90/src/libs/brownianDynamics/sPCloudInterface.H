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
    sPCloudInterface

Description
    Interface class used for Brownian dynamics simulations.   
      
    This class is part of rheoTool.

SourceFiles
    sPCloudInterface.C

\*---------------------------------------------------------------------------*/

#ifndef sPCloudInterface_H
#define sPCloudInterface_H

#include "solidParticleCloud.H"
#include <random>
#include "springModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class springModel;

/*---------------------------------------------------------------------------*\
                        Class sPCloudInterface Declaration
\*---------------------------------------------------------------------------*/
 
class sPCloudInterface
{

private:
        
    // Private data
        
         //- Reference to continuum U field 
         const volVectorField& U_;
         
         //- Mersenne Twister generator for random numbers 
         std::mt19937  gen_;
         
         //- Random number uniform distribution 
         std::uniform_real_distribution<> randUN_;
         
         //- Dictionary with molecule controls
         IOdictionary molcDict_;
         
         //- True if found molecule data structure in the case dir
         bool exist_;
        
         //- solidParticleCloud class charged of particle creation and tracking 
         solidParticleCloud spc_;
         
         //- To write or not continuum fields
         bool writeContFields_;
         
         //- Dictionary containing physical properties, active molecules and springs
         IOdictionary runTimeInfoDict_;
         
         //- Number of active molecules (the ones which are still tracked)
         int nMolc_;
         
         //- (x y z) position of each bead, for each molecule. [molecule][bead] = vector 
         PtrList<Field<vector > > mx_;
         
         //- This is mx_ saved at the beginning of a new loop (all springs are Ls-compliant) 
         PtrList<Field<vector > > mx0_;
         
         //- This is mx_ saved at the beginning of a new simulation
         PtrList<Field<vector > > mxFix_;
         
         //- (x y z) position of each bead, for each molecule, after accounting for drag, EV and Brownian forces.
         PtrList<Field<vector > > mxStar_;
         
         //- (u v w) velocity of each bead, for each molecule. [molecule][bead] = vector 
         PtrList<Field<vector > > mU_;
         
         //- List of active molecules. [i] = labelList = (moleculeID number_of_beads groupsID)
         // On starting, all [i][0] = moleculeID =! -1, but as the simulation runs
         // untracked molecules are flaged with moleculeID = -1, which are not written
         // when mAct_ gets written in runTimeInfoDict_.   
         PtrList<List<label > > mAct_;
         
         //- List of indexes for each bead, of each molecule. [molecule][bead] = labelList = (globalID localID groupID)  
         PtrList<List<List<label > > > mIds_;
         
         //- N x N list of diffusion tensors, for each molecule. [molecule][beadi][beadj] = symmTensor  
         PtrList<List<List<symmTensor > > > mD_;
         
         //- N x N list of sigma tensors (Brownian force term), for each molecule. [molecule][beadi][beadj] = symmTensor  
         PtrList<List<List<tensor > > > mSigma_;
         
         //- List of springs for each molecule. [molecule][spring] = labelList = (beadi beadi+1 molecule) 
         PtrList<List<List<label > > > mSpr_;
         
         //- 1D list for linking. The ith element of linkM_ contains the index in all m#_ PtrList<>, which corresponds to 
         // moleculeID (name) i. Once a molecule with moleculeID (name) i gets untracked, then linkM_[i] is set to -1.
         // The size of linkM_ correponds to max(moleculeID), but only the elements at a position i corresponding to any valid
         // moleculeID are different from -1. To be removed after code cleanup.         
         List<label> linkM_;  
         
         //- 1D list of groups names. (all: active and not ative)
         List<word>  names_;
         
         //- 1D list of Nks (number of Kuhn steps per spring) for each group. (all: active and not ative)
         List<scalar> Nks_;
         
         //- 1D list of nuEV (nu^{EV} parameter in the EV force) for each group. (all: active and not ative)
         List<scalar> nuEV_;
         
         //- 1D list of D (isotropic diffusion coefficient = kT/(6a\pi\eta) ) for each group. (all: active and not ative)
         List<scalar> D_;
         
         //- 1D list of a (beads diameter) for each group. (all: active and not ative)
         List<scalar> a_;
         
         //- 1D list of Ls (individual spring length upon full extension) for each group. (all: active and not ative)
         List<scalar> Ls_;
         
         
         //- Data output
         //-- Should the molecules ID/position/stretch be outputed over time?   
         Switch writeStats_;
         
         //-- Should the molecules VTK representation be outputed over time?   
         Switch writeVTK_;
         
         //-- Directory to output VTK
         fileName VTKDir_;
         
         //-- Number of time-steps between output calls
         int outputStatsN_;
         
         //-- Ouput streams. One per group per field.
         PtrList<PtrList<autoPtr<OFstream > > > outS_;
         
         //-- Counter to output data 
         int ppCount_;
         
         
         //- Algorithm controls
         //-- Include exclusion volume forces?
         bool isExclusionVolumeF_;
         
         //-- Account for hydrodynamic interactions?
         bool isHI_;
         
         //-- Vector whose components indicate the molecules center of mass components to fix. Example: (1 0 0) fixes the x-component,
         // (1 1 1) fixes all the components, (0 0 0) fixes none, ... All the vector components should be either 1 or 0.
         // Not used if the molecule is tethered. 
         vector pBackV_; 
         
         //-- Push back frequency 
         int pBackFreq_;
         
         //-- Push back counter 
         int pBackCount_;
         
         //-- Are the molecules tethered (fixed position of their first bead; ALWAYS first bead ONLY)?  
         bool isTethered_;
         
         //- Spring model object 
         autoPtr<springModel>   spModel_;
         
    // Private Member Functions

        //- Disallow default bitwise copy construct
        sPCloudInterface(const sPCloudInterface&);

        //- Disallow default bitwise assignment
        void operator=(const sPCloudInterface&);
        
        //- Checks if the molecules data strutuctre exists in the case directory (does not check if it is valid).
        bool checkIfExist();
        
        //- Writes runTimeInfoDict_ when time == outputTime. 
        void writeRunTimeInfoDict();
        
        //- Writes the molecules in VTK format
        void writeVTK();
        
        //- Calls the particle tracker and updates mx_ with the resulting beads positions (conversion from barycentric). 
        //  Deletes any particle belonging to an untracked molecule. The argument indicates if the drag force is to
        //  be added or not to mU_ during the tracking. 
        void moveAndReceive(bool);
        
        //- Copies mU_ (beads velocity) to spc
        void sendU(); 
        
        //- Computes tensor mD_ and uses Cholesky decomposition to compute mSigma_ 
        void computeDSigma();
        
        //- Brownian force contribution to mU_ (first resets mU_; important for the calling order)
        void fBrownian();
        
        //- Exclusion volume force contribution to mU_ (adds to mU_; important for the calling order)
        void fEV();  
        
        //- Writes ID/position/stretch of each molecule (lost molecules get stretch = 0 from the moment they are lost; this is also a flag)
        // Splite the write considering groups.
        void writeStatistics(); 
        
        //- Resets the molecules center of mass
        void pushToX0();
                 

public:

    //- Runtime type information
    TypeName("sPCloudInterface");
    

    // Constructors

        //- Construct from U
        sPCloudInterface
        (
  	    const volVectorField& U
        );


    // Destructor
       virtual ~sPCloudInterface()
        {}


    // Member Functions
    
       //- Deletes a lost molecule, which becomes not active and untracked. Set -1 flags in linkM_ and actM_, and clears all m#_
       // (keeps pointer, but calls element dtor)  
       void deleteMolecule(label);
       
       //- Updates the molecules positions and configurations
       bool update();
       
       //- Writes the molecules data structure to the case directory (this emulates the AUTO_WRITE option for the class)
       void writeM();
       
       //- Return the continuum velocity field
       const volVectorField& U() const
       {
            return U_;
       }
        
       //- Access and Modify
       
       PtrList<Field<vector > >& mx()
       {
          return mx_;
       }
       
       PtrList<Field<vector > >& mU()
       {
          return mU_;
       }
       
       PtrList<List<label > >& mAct()
       {
          return mAct_;
       }
       
       PtrList<List<List<label > > >& mIds()
       {
          return mIds_;
       }
       
       PtrList<List<List<symmTensor > > >& mD()
       {
          return mD_;
       }
       
       PtrList<List<List<tensor > > >& mSigma()
       {
          return mSigma_;
       }
       
       PtrList<List<List<label > > >& mSpr()
       {
          return mSpr_;
       }
       
       List<label>& linkM()
       {
          return linkM_;
       } 
       
       int& nMolc()
       {
          return nMolc_;
       }
       
       // Const-access    
       const List<scalar>& Nks() const
       {
          return  Nks_;
       }
       
       const List<scalar>& nuEV() const
       {
          return  nuEV_;
       }
       
       const List<scalar>& D() const
       {
          return  D_;
       }
       
       const List<scalar>& Ls() const
       {
          return  Ls_;
       }   
       
       const bool& isTethered() const
       {
          return isTethered_;
       }  
       
       const bool& isHI() const
       {
          return isHI_;
       }      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
