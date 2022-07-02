// Optimized for TDMA method, since no data transfer
// is needed. However, TDMA can only be used for
// tridiagonal matrices, ie, linear molecules.

for (int cmpI=0; cmpI<3 ; cmpI++)
{
  if (solverType_ == "TDMA")
   {

     TDMA(J[cmpI], f[cmpI]);  // f[cmpI] is dX on return.   
     xn.replace(cmpI, xn.component(cmpI) + f[cmpI]);
	 
     // Reset matrix     
     // No need due to direct assignment and invariance of non-zero indexes.  
   }
  else if (solverType_ == "Gaussian")
   {
     scalarDiagonalMatrix rhs(f[cmpI].size(), 0.);
     forAll(rhs, lk)
      {
        rhs[lk] = -f[cmpI][lk]; 
      }
	  
     solve(J[cmpI], rhs);
	 
     forAll(rhs, lk)
      {
	xn[lk].component(cmpI) += rhs[lk]; 
      }  
     
      // Need to reset J, because errors may accumulate in 0 elements
      for (int i = 0; i<J[cmpI].n(); i++)
       {            
        for (int j = 0; j<J[cmpI].m(); j++)
         {
           J[cmpI][i][j] = 0.;
         }
       }    

   }
  else if (solverType_ == "LLT" || solverType_ == "QR")
   {    
     Eigen::VectorXd b(n);
     Eigen::MatrixXd A(Eigen::MatrixXd::Zero(n,n));

     // Fill b and diagonal of A
     forAll(f[cmpI], i)
      {
        b(i) = -f[cmpI][i]; 
        A(i,i) = J[cmpI][i][i];
      }

     // Fill off-diagonal of A (only non-zero for spring elements)	 
     forAll(mSpr_[mi], bi)
      {
        label& b0 = mSpr_[mi][bi][0];
        label& b1 = mSpr_[mi][bi][1];
	   
        A(b0,b1) = J[cmpI][b0][b1]; 
        A(b1,b0) = J[cmpI][b1][b0]; 
      }

     if (solverType_ == "LLT")
      {
        b = A.llt().solve(b);
      }
     else
      {
        b = A.colPivHouseholderQr().solve(b);
      }
     
     forAll(b, ii)
      {
        xn[ii].component(cmpI) += b(ii); 
      }
   }
  else
   {
     FatalErrorIn("springModel::implicitfSpring()")
  	<< "\nUnknow solver type:" << solverType_ 
	<< "\nValid solvers are: \n .TDMA \n .Gaussian \n .LLT \n .QR \n"
	<< exit(FatalError);
   }
}
