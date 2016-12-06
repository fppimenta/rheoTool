#ifndef jacobi_H 
#define jacobi_H 

namespace Foam
{

inline int jacobi(const Foam::symmTensor& At, int& N, Foam::tensor& Dt, Foam::tensor& Vt, int& NROT) {

/*
 Compute eigenvalues/vectors by the jacobi method. Converted from fortran rountine:
http://www.public.iastate.edu/~jba/Fortran/jacobi.f
Modified to return the exponential of the eigenvalues.
*/


Foam::scalar  c,g,h,s,sm,t,tau,theta,tresh;



int     i,j,ip,iq;


  Eigen::Vector3d B(0.0, 0.0, 0.0);
  Eigen::Vector3d Z(0.0, 0.0, 0.0);
  Eigen::Vector3d D(0.0, 0.0, 0.0);
  Eigen::MatrixXd V=Eigen::MatrixXd::Identity(3, 3);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);


  for (i=0; i<9; i++) {


       A(0,0)=At.xx();
       A(1,1)=At.yy();
       A(2,2)=At.zz();

       A(0,1)=At.xy();
       A(1,0)=At.xy();

       A(0,2)=At.xz();
       A(2,0)=At.xz();

       A(1,2)=At.yz();
       A(2,1)=At.yz();

  }



  for (ip=0; ip<=N-1; ip++) {
    B(ip)=A(ip,ip);
    D(ip)=B(ip);
    Z(ip)=0;    
  }
  NROT=0;

  for (i=1; i<=50; i++) {
    sm=0;
    for (ip=0; ip<N-1; ip++)    //sum off-diagonal elements
      for (iq=ip+1; iq<=N-1; iq++)
	    sm=sm+Foam::mag(A(ip,iq));
	
    if (i < 4)
      tresh=0.2*sm*sm;
    else
      tresh=0;
    for (ip=0; ip<N-1; ip++) {
	  for (iq=ip+1; iq<=N-1; iq++) {
        g=100*Foam::mag(A(ip,iq));
// after 4 sweeps, skip the rotation if the off-diagonal element is small
        if ((i > 4) && (Foam::mag(D(ip))+g == Foam::mag(D(ip))) && (Foam::mag(D(iq))+g == Foam::mag(D(iq))))
		  A(ip,iq)=0;
		else if (Foam::mag(A(ip,iq)) > tresh) {
	      h=D(iq)-D(ip);
	      if (Foam::mag(h)+g == Foam::mag(h))
	        t=A(ip,iq)/h;
		  else {
	        theta=0.5*h/A(ip,iq);  
            t=1/(Foam::mag(theta)+Foam::sqrt(1.0+theta*theta));
	        if (theta < 0)  t=-t;
		  }
	      c=1.0/Foam::sqrt(1.0+t*t);
	      s=t*c;
              tau=s/(1.0+c);
	      h=t*A(ip,iq);
	      Z(ip) -= h;
	      Z(iq) += h;
	      D(ip) -= h;
	      D(iq) += h;
	      A(ip,iq)=0;
		  for (j=0; j<ip; j++) {
	        g=A(j,ip);
	        h=A(j,iq);
	        A(j,ip) = g-s*(h+g*tau);
	        A(j,iq) = h+s*(g-h*tau);
		  }
		  for (j=ip+1; j<iq; j++) {
	        g=A(ip,j);
	        h=A(j,iq);
	        A(ip,j) = g-s*(h+g*tau);
	        A(j,iq) = h+s*(g-h*tau);
		  }		      
		  for (j=iq+1; j<=N-1; j++) {
	        g=A(ip,j);
	        h=A(iq,j);
	        A(ip,j) = g-s*(h+g*tau);
		    A(iq,j) = h+s*(g-h*tau);
		  }		  
		  for (j=0; j<=N-1; j++) {
	        g=V(j,ip);
	        h=V(j,iq);
	        V(j,ip) = g-s*(h+g*tau);
	        V(j,iq) = h+s*(g-h*tau);
		  }		  
          NROT=NROT+1;
		} //end ((i.gt.4)...else if
      } // main iq loop
    } // main ip loop
    for (ip=0; ip<=N-1; ip++) {
      B(ip) += Z(ip);
      D(ip)=B(ip);
      Z(ip)=0;
    }
  } //end of main i loop
  //printf("\n 50 iterations !\n");
  
  
  for (i=0; i<9; i++) {


       Vt.xx()=V(0,0);
       Vt.yy()=V(1,1);
       Vt.zz()=V(2,2);

       Vt.xy()=V(0,1);
       Vt.yx()=V(1,0);

       Vt.xz()=V(0,2);
       Vt.zx()=V(2,0);

       Vt.yz()=V(1,2);
       Vt.zy()=V(2,1);

       Dt.xx()=Foam::exp(D(0));
       Dt.yy()=Foam::exp(D(1)); 
       Dt.zz()=Foam::exp(D(2));

      // Dt.xx()=D(0);
       //Dt.yy()=D(1); 
       //Dt.zz()=D(2);

  }


  return(0);  //too many iterations
  
 
}

} //end nameSpace


#endif 

