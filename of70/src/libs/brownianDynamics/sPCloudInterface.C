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

#include "sPCloudInterface.H"
#include "sPCloudInterfaceIO.H"
#include <string> 
#include "IFstream.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sPCloudInterface, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sPCloudInterface::sPCloudInterface
(
    const volVectorField& U
) 
:
U_(U),
gen_(std::random_device()()),
randUN_(-1, 1),
molcDict_
(
  IOobject
  (
    "moleculesControls",
    U.time().constant(),
    U.mesh(),
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
),
exist_(checkIfExist()),
spc_(U.mesh(), "molecules"),
writeContFields_(molcDict_.subDict("externalFlow").lookupOrDefault<Switch>("writeFields", true)),
runTimeInfoDict_
(
  IOobject
  (
    "MoleculesInfo",
    U.time().constant(),
    "runTimeInfo"/U.time().timeName(),
    U.mesh(),
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
  )
),
nMolc_(readInt(runTimeInfoDict_.subDict("ActiveMolecules").lookup("nActive"))),
mx_(nMolc_),
mx0_(nMolc_),
mxFix_(nMolc_),
mxStar_(nMolc_),
mU_(nMolc_),
mAct_(nMolc_),
mIds_(nMolc_),
mD_(nMolc_),
mSigma_(nMolc_),
mSpr_(nMolc_),
linkM_(),
names_(runTimeInfoDict_.subDict("GroupProperties").lookup("names")),
Nks_(runTimeInfoDict_.subDict("GroupProperties").lookup("Nks")),
nuEV_(runTimeInfoDict_.subDict("GroupProperties").lookup("nuEV")),
D_(runTimeInfoDict_.subDict("GroupProperties").lookup("D")),
a_(runTimeInfoDict_.subDict("GroupProperties").lookup("a")),
Ls_(runTimeInfoDict_.subDict("GroupProperties").lookup("Ls")),
writeStats_(molcDict_.subDict("outputOptions").lookupOrDefault<Switch>("writeStats", false)),
writeVTK_(molcDict_.subDict("outputOptions").lookupOrDefault<Switch>("writeVTK", false)),
VTKDir_(U.time().path()/"VTKMolecules"/U.time().timeName()),
outputStatsN_(0),
outS_(names_.size()),
ppCount_(1),
isExclusionVolumeF_(readBool(molcDict_.subDict("exclusionVolumeProperties").lookup("activeExclusionVolume"))),
isHI_(readBool(molcDict_.subDict("HIProperties").lookup("activeHI"))),
pBackV_(0.,0.,0.),
pBackFreq_(0),
pBackCount_(1),
isTethered_(readBool(molcDict_.subDict("externalFlow").lookup("tethered"))),
spModel_()
{ 

 // Push back vector read and normalize (components should be either 0 or 1)
 if (!isTethered_)
  {
    molcDict_.subDict("externalFlow").lookup("pushBackCmp") >> pBackV_;
    molcDict_.subDict("externalFlow").lookup("pushBackFreq") >> pBackFreq_;
    if (mag(pBackV_.x())>0)
     pBackV_.x() = 1.;
    if (mag(pBackV_.y())>0)
     pBackV_.y() = 1.;
    if (mag(pBackV_.z())>0)
     pBackV_.z() = 1.;
  }

 // Active molecules
 List<List<label> > actMolc_(runTimeInfoDict_.subDict("ActiveMolecules").lookup("idList"));
 
 int mMax(0); // Max name id
 for (int i=0; i<nMolc_; i++)
  { 
    // Act molec
    mAct_.set
     (
        i,
        new List<label>(3,0) 
     );
       
    mAct_[i] = actMolc_[i]; 
    int nBeads(mAct_[i][1]);
    
    // Create beads id
    mIds_.set
     (
        i,
        new List<List<label> >(nBeads, List<label>(3,0)) 
     );
     
    // Create positions
    mx_.set
     (
        i,
        new Field<vector>(nBeads, vector::zero) 
     );
     
    mx0_.set
     (
        i,
        new Field<vector>(nBeads, vector::zero) 
     );
    
    mxFix_.set
     (
        i,
        new Field<vector>(nBeads, vector::zero) 
     );
     
    mxStar_.set
     (
        i,
        new Field<vector>(nBeads, vector::zero) 
     );
     
    // Create velocity
    mU_.set
     (
        i,
        new Field<vector>(nBeads, vector::zero) 
     );
     
    // Create beads D and Sigma if HI
    if (isHI_) 
     {
       mD_.set
        (
           i,
           new List<List<symmTensor> >(nBeads, List<symmTensor>(nBeads, symmTensor::zero)) 
        );
     
       mSigma_.set
        (
           i,
           new List<List<tensor> >(nBeads, List<tensor>(nBeads, tensor::zero)) 
        );
     }  
     
    // Create springs
    mSpr_.set
     (
        i,
        new List<List<label> >(nBeads-1, List<label>(3,0)) 
     );
     
    // Find max active id for the linking matrix        
    if (mAct_[i][0] > mMax)
     {
       mMax = mAct_[i][0];    
     }
  } 
  
 linkM_.clear(); 
 linkM_.setSize(mMax+1,-1);
 
 // Warn: assuming the same order of initilization for actMolc and linkM
 forAll(actMolc_, i)
  {
    linkM_[actMolc_[i][0]] = i; 
  }
  
 // Fill beads index and positions 
 forAllIter(Cloud<solidParticle>, spc_, iter)
  {
    solidParticle& p = iter();
    
    label molcN(linkM_[p.molcID_]);
    label beadN(p.ids_[1]);
        
    mIds_[molcN][beadN] = p.ids_;
    mx_[molcN][beadN] = p.position();	
    mU_[molcN][beadN] = p.U_;	
  }
 
 // Set or clear mxFix
 if (!isTethered_ && mag(pBackV_)>SMALL)
  {
    mxFix_ = mx_;
  }
 else
  {
    mxFix_.clear(); 
  }
     
 // Fill the springs
 IOList<List<label> > spr
  (
    IOobject
     (
       "springs",
       U.time().constant(),
       "runTimeInfo"/U.time().timeName(),
       U.mesh(),
       IOobject::MUST_READ,
       IOobject::NO_WRITE,
       false
     )
  );
  
 List<label>  posmSpr(mSpr_.size(), 0);
 forAll(spr, si)
 { 
    List<label> sprI = spr[si];     
    label mi = linkM_[sprI[2]];      
    label bi = posmSpr[mi];
   
    // If the branches are not closed the size is fine and there is no
    // redimensioning. If not, one additional spring exists per closed
    // branch, so append (redim+assign)
    if (bi < mSpr_[mi].size() )
     {
       mSpr_[mi][bi] = sprI;
     }
    else
     {
       mSpr_[mi].append(sprI);
     }
    posmSpr[mi]++;
 }
   
 spr.clear();
 posmSpr.clear();
   
 // Spring model
 spModel_ = springModel::New(molcDict_, U, *this);
 
 // Stats post-processing  
 if (writeStats_)
  {
   fileName ppDirRoot(U.time().path()/"rheoToolPP"/U.time().timeName()/"moleculesStats");
    
   // Create one [ X | Stretch | SumSprings | IDs ] file for each group
   forAll(names_, gi)
    {
       fileName ppDir(ppDirRoot/names_[gi]);
      
       mkDir(ppDir);
       
       outS_.set(gi, new PtrList<autoPtr<OFstream > >(3));
       outS_[gi].set(0, new autoPtr<OFstream >); outS_[gi][0].reset(new OFstream(ppDir/"X.txt"));
       outS_[gi].set(1, new autoPtr<OFstream >); outS_[gi][1].reset(new OFstream(ppDir/"Stretch.txt"));
       outS_[gi].set(2, new autoPtr<OFstream >); outS_[gi][2].reset(new OFstream(ppDir/"IDs.txt"));
       
       // Make header for indexes files
       outS_[gi][2]() << "Name" << tab << "Nb_of_beads" << tab << "Group" << nl;      
    }
 
   // Fill indexes for each group     
   forAll(actMolc_, i)
    {
       int fi = actMolc_[i][2]; // GroupID = 1st index of PtrList outS_
       outS_[fi][2]() << actMolc_[i][0] << tab << actMolc_[i][1] << tab << actMolc_[i][2] << nl; 
    } 
    
   // Flush stream
   forAll(outS_, gi)
    {     
      outS_[gi][2]() << endl;
    } 
    
   molcDict_.subDict("outputOptions").lookup("outputStatsInterval") >> outputStatsN_; 
  } 

 // VTK files
 if (writeVTK_)
  {
   mkDir(VTKDir_);
   // Call the writer on first dt
   writeVTK();
  }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool sPCloudInterface::update() 
{
  // Compute Brownian and Exclusion Volume force terms. Save mx0_ in case it is 
  // needed in the implicit corrector for spring force. In that case mx0_ is used
  // as the starting guess for x. We could use mxStar_ saved next for that, but mxStar_
  // may already have overstretched springs and this would be a bad starting 
  // guess. Using mx0_ we guarantee that there is no overstretch.
  // Order is important: 1st Fb and then fEV because velocity is incremented. 
  
  mx0_ = mx_;
  
  if (isHI_) 
     computeDSigma();
     
  fBrownian();
    
  if (isExclusionVolumeF_)   
     fEV(); 
   
  sendU(); 
  
  // Move particles: Drag + Brownian + Exclusion Volume
  
  moveAndReceive(true);
  
  // Compute Spring force term explicitly and update local mU and mx. We save mxStar_ 
  // in case it is needed in the implicit corrector for spring force. In that case 
  // mxStar_ represents the indermediary beads positions after Drag+B+EV and it is the
  // actual value of p.position_.  
  mxStar_ = mx_;
  spModel_->fSpring();
  
  // Check for violations in spring lengths and add spring force implicitly if needed
  spModel_->checkSpringsLength(mxStar_, mx0_);
  
  // Push back the molecules if not tethered and if pushback vector is not negligible 
  if (!isTethered_ && mag(pBackV_)>SMALL)
    pushToX0();
  
  sendU();
   
  // Move particles: spring force only  
  moveAndReceive(false);
  
  // Write data (controlled by output time)
  writeM();
 
  // Write statistics (has its own control for output)
  if (writeStats_)
    writeStatistics(); 
    
  if (nMolc_>0)
   {
     return true;
   }
  else
   {
     return false;
   } 
}

void sPCloudInterface::moveAndReceive(bool includeDrag) 
{
  
  // Get the name of the molc to remove. We remove a molc if at least one 
  // of its beads is lost. In general, molcToDelete may have repeated molecules
  // to remove, so this is why we check if molcN != -1 is already not set. 
  List<label> molcToDelete(spc_.move(includeDrag));
   
  forAll(molcToDelete, i)
   {
     label molcN(linkM_[molcToDelete[i]]);
     
     if (molcN != -1)
      {
        deleteMolecule(molcN);     
      }     
   }   
  
  // Fill local beads index and positions 
  forAllIter(Cloud<solidParticle>, spc_, iter)
    {
      solidParticle& p = iter();
      
      label molcN(linkM_[p.molcID_]);
      label beadN(p.ids_[1]);
      
      if (molcN != -1)
       {          
         mIds_[molcN][beadN] = p.ids_; // optional (constant since ctor)
         mx_[molcN][beadN] = p.position(); // tensorial operation involved (!cput)    
       }	
      else
       {
        // This will remove all the particles which are still tracked, but
        // that belong to a molc that has been deleted, because one of its
        // beads get lost.
        
         spc_.deleteParticle(p);
       }
    }
}

void sPCloudInterface::sendU() 
{

  // Transfer U to spc
  forAllIter(Cloud<solidParticle>, spc_, iter)
    {
      solidParticle& p = iter();
     
      label molcN(linkM_[p.molcID_]);
      label beadN(p.ids_[1]);
      
      if (molcN != -1)
       {
         p.U_ = mU_[molcN][beadN];    
       }	
    }   
}

void sPCloudInterface::deleteMolecule(label mi)
{
    
  //- Better to remove the molecule to avoid unphysical conformations due to overstretch 
  // Flag                        
  linkM_[mAct_[mi][0]] = -1;
  mAct_[mi][0] = -1;

  // Clear the list but keep the pointer to it
  mx_[mi].clear();  
  mU_[mi].clear();    
  mIds_[mi].clear();
  mSpr_[mi].clear(); 
  if (isHI_)
  {
    mD_[mi].clear();
    mSigma_[mi].clear();
  }

  // Update counter of active molecules
  nMolc_--;   																
}

void sPCloudInterface::computeDSigma() 
{
  symmTensor I_(symmTensor::I);
   
  forAll(mD_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
          label gI(mIds_[mi][0][2]); // All beads belong to the same group, thus check for the first bead and save         
          scalar D = D_[gI];
          scalar a = a_[gI];
          
          int sz = mD_[mi].size();
          int n = 3*sz;
                
          scalarField A(n*n, 0.);
          int row = 0;
          int col = 0;
         
          // Loop over all beads
          forAll(mD_[mi], i)
           {
             // Loop over the lower triangular
             for (int j=0; j<=i; j++)
              {
                 if (i == j)
                  {
                    mD_[mi][i][j] = D*I_;
                  }
                 else
                  {
                    vector rij(mx_[mi][j] - mx_[mi][i]); 
                    scalar mrij(mag(rij));
                    
                    if (mrij>=2.*a)
                     {
                       mD_[mi][i][j] =
                        ( 3.*D*a/(4*mrij) )
                       *(
                            ( 1. + (2./3.)*(a/mrij)*(a/mrij) )*I_
                          + ( 1. - 2.*(a/mrij)*(a/mrij) ) * symm(rij*rij)/(mrij*mrij)
                        );
                     }
                    else
                     {
                       mD_[mi][i][j] =
                        D
                       *(
                            ( 1. - 9.*mrij/(32.*a) )*I_
                          + 3.*symm(rij*rij)/(32.*mrij*a)
                        );
                     }
                    
                    // Fill also the upper triangular (symmetry). #Wasting memory. 
                    mD_[mi][j][i] = mD_[mi][i][j];
                                      
                  }
                 
                 // Populate A to be used in Cholesky decomposition. A is equivalent 
                 // to mD; the individual tensors in mD become concatenated in A for
                 // sequential access. Only the lower triangular is copied, because
                 // mD is symmetric and only the lower triangular is used in Cholesky
                 // decomposition. The upper triangular remains as initialized (0).
                 // Although only lower is filled, A is full size. Accessor is row
                 // first, col then.  
                 symmTensor& tt = mD_[mi][i][j]; 
                 
                 A[n*row + col] = tt.xx();
                 A[n*row + col+1] = tt.xy(); 
                 A[n*row + col+2] = tt.xz(); 
                 
                 A[n*(row+1) + col] = tt.xy();
                 A[n*(row+1) + col+1] = tt.yy(); 
                 A[n*(row+1) + col+2] = tt.yz();
                 
                 A[n*(row+2) + col] = tt.xz();
                 A[n*(row+2) + col+1] = tt.yz(); 
                 A[n*(row+2) + col+2] = tt.zz();  
                 
                 col = col + 3;                 
              }
              
              col = 0;
              row = row + 3;
           }
           
        
          //-* * * * * Cholesky Decomposition (from: C Rosetta.org) * * * * *-//
          // Note: one could use directly mD and mSigma instead of relying on 
          // L and A. However, it is slower due to row/col mapping and nested
          // accessors (at least 4) in the ptrLists of tensors.
         
          scalarField L(n*n, 0.);
                
          for (int i = 0; i < n; i++)
          for (int j = 0; j < (i+1); j++) {
             scalar s = 0;
               for (int k = 0; k < j; k++)
                 s += L[i * n + k] * L[j * n + k];
               L[i * n + j] = (i == j) ?
                                        Foam::sqrt(A[i * n + i] - s) :
                                        (1.0 / L[j * n + j] * (A[i * n + j] - s));
          }
          
         row = 0;
         col = 0;
         
         //- Transfer solution to mSigma
         forAll(mD_[mi], i)
         {
           for (int j=0; j<=i; j++)
            {
              tensor& tt = mSigma_[mi][i][j]; 
              
              tt.xx() = L[n*row + col];
              tt.xy() = L[n*row + col+1]; 
              tt.xz() = L[n*row + col+2]; 
                 
              tt.yx() = L[n*(row+1) + col];
              tt.yy() = L[n*(row+1) + col+1]; 
              tt.yz() = L[n*(row+1) + col+2];
                 
              tt.zx() = L[n*(row+2) + col];
              tt.zy() = L[n*(row+2) + col+1]; 
              tt.zz() = L[n*(row+2) + col+2];  
                 
              col = col + 3;
            }
           col = 0;
           row = row + 3;
         }
            
       }       
   }   
}
 
void sPCloudInterface::fBrownian() 
{

  scalar dt(U().mesh().time().deltaTValue());
  scalar f( Foam::sqrt(6./dt) );
  
  forAll(mU_, mi)
   {
      if (mAct_[mi][0] != -1)
       {   
         if (isHI_)
          {                        
            forAll(mU_[mi], bi)
             {    
               // Need to reset because of the loop next
               vector& U = mU_[mi][bi]; 
               U *= 0.;
            
               for (int bj = 0; bj<=bi; bj++)
                {        
                  vector Rnd(randUN_(gen_),randUN_(gen_),randUN_(gen_));           
                  U += f * (mSigma_[mi][bi][bj] & Rnd);            
                }
             }            
          }
         else
          {
            label gI(mIds_[mi][0][2]); // All beads belong to the same group, thus check for the first bead and save
            scalar fac( D_[gI]/dt );
            
            forAll(mU_[mi], bj)
             {         
               vector Rnd(randUN_(gen_),randUN_(gen_),randUN_(gen_)); 
            
               mU_[mi][bj] = Foam::sqrt(6.*fac) * Rnd; 
             }
          }
         
        // Correct tethered case 
        if (isTethered_)
         {
            mU_[mi][0] *= 0;
         } 
                 
       }     
   }
}

void sPCloudInterface::fEV() 
{
  // The matrix inside the sum is anti-symmetric, thus run upper-triang only
  // and fill both triangs.
 
   forAll(mU_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         label gI(mIds_[mi][0][2]); // All beads belong to the same group, thus check for the first bead and save
         
         scalar Ls(Ls_[gI]);
         scalar Nks(Nks_[gI]);
         
         scalar fac = (9./2.)*nuEV_[gI]*Foam::pow(3./(Ls*4.*Foam::sqrt(M_PI)), 3)*Foam::pow(2.*Nks, 4.5);
       
         vector vi(0.,0.,0.);
         
         int nb(mU_[mi].size());
         
         vectorField mUI(nb, vector::zero);
       
         // First compute the EV force
         forAll(mU_[mi], bj)
          {
            vector& xj = mx_[mi][bj];
            
            // No worries for bk==bj, because rjk = 0
            for (int bk=bj; bk<nb ; bk++)
             {
                vector rjk( (mx_[mi][bk] - xj)/ Ls );
                
                // This is one bottleneck of the code. Computing the exp() is expensive.
                // The CPU time is OS and compiler dependent.
                  
             // vi = fac * Foam::exp( (-9./2.) * Nks * Foam::magSqr(rjk) ) * rjk;               
                vi = fac * Foam::pow(2.718281828459045, (-9./2.) * Nks * Foam::magSqr(rjk) ) * rjk;
                               
                mUI[bj] -= vi;
                mUI[bk] += vi;                 
             }     
          }
  
         if (isHI_)
          { 
            // ... then dot it with D and save to mU         
            forAll(mU_[mi], i)
             {
               vector& U = mU_[mi][i];
               forAll(mU_[mi], j)
               { 
                  // The multiplication is for dimensions (kT/L)/kT
                  U += (mD_[mi][i][j] & mUI[j]) / Ls;  
               } 
             } 
          }
         else
          {
            mU_[mi] += mUI * D_[gI] / Ls;
          }
          
         // Correct tethered case 
         if (isTethered_)
          {
             mU_[mi][0] *= 0;
          }
        
       }     
   }    
}

void sPCloudInterface::pushToX0() 
{

// Update counter
pBackCount_++;

if(pBackCount_>pBackFreq_)
 {
  // Change internal mx and mU
  forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {          
         // Find old center of mass
         vector cM0(0.,0.,0.);         
         forAll(mxFix_[mi], bj)
          {            
            cM0 += mxFix_[mi][bj];      
          }
         cM0 /= mxFix_[mi].size();
         
         // Find new center of mass
         vector cMf(0.,0.,0.);         
         forAll(mx_[mi], bj)
          {            
            cMf += mx_[mi][bj];      
          }
         cMf /= mx_[mi].size();
         
         // Push back the molecule
         vector dX( cM0-cMf );
         dX.x() *= pBackV_.x(); dX.y() *= pBackV_.y(); dX.z() *= pBackV_.z();
         forAll(mx_[mi], bj)
          {            
            mx_[mi][bj] += dX;
            mU_[mi][bj] += dX/U().mesh().time().deltaTValue();      
          }       
       }     
   }
   
 // Reset counter  
 pBackCount_ = 1;
 }
 
}
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
