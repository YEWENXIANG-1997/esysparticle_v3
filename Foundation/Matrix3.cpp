/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#include "Foundation/Matrix3.h"
#include "Foundation/Matrix3.hpp"
#include "Foundation/cube_eq.h" // for eigenvalues

#include <iostream>

// --- STL includes ---
#include <set>

using std::set;

// Solve 3x3 equation system via determinants
// throws exception MatSingularError() in case of a singular matrix
// needs 59 Flops ( 36 mult , 20 add/sub , 3 div )  
Vec3 Matrix3::solve(const Vec3& rhs) const
{
  double det1,detx,dety,detz;

  det1=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if(det1 == 0){
    throw  MatSingularError();
  }
  detx=rhs.data[0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*rhs.data[2]-rhs.data[1]*m[2][2])+m[0][2]*(rhs.data[1]*m[2][1]-m[1][1]*rhs.data[2]);
  dety=m[0][0]*(rhs.data[1]*m[2][2]-m[1][2]*rhs.data[2])+rhs.data[0]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*rhs.data[2]-rhs.data[1]*m[2][0]);
  detz=m[0][0]*(m[1][1]*rhs.data[2]-rhs.data[1]*m[2][1])+m[0][1]*(rhs.data[1]*m[2][0]-m[1][0]*rhs.data[2])+rhs.data[0]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  double x=detx/det1;
  double y=dety/det1;
  double z=detz/det1;
  return Vec3(x,y,z);
}    

// in situ inversion
// based on a`[ij]=det(A[ij])/det(M)
// where A[ij] is the adjoined submatrix to a[ij]
// 50 Flops (27 mult, 14 add/sub, 9 div)
void Matrix3::invert()
{
  double det1,temp[3][3];

  det1=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if(det1 == 0){ // det(M)==0 -> singular matrix
    throw  MatSingularError();
  } 
  temp[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det1;
  temp[0][1]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det1;
  temp[0][2]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det1;
  temp[1][0]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])/det1;
  temp[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det1;
  temp[1][2]=(m[0][1]*m[2][0]-m[0][0]*m[2][1])/det1;
  temp[2][0]=(m[1][2]*m[0][1]-m[1][1]*m[0][2])/det1;
  temp[2][1]=(m[1][0]*m[0][2]-m[1][2]*m[0][0])/det1;
  temp[2][2]=(m[1][1]*m[0][0]-m[1][0]*m[0][1])/det1;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      m[i][j]=temp[i][j];
    }
  }
}

/*!
  Solve the homogeneous system Mx=(0,0,0). It is assumed that rank(M)=2.
  A unit vector is returned;
*/
Vec3 Matrix3::solve_homogeneous() const
{
  Vec3 res;
  
  // check if any colums vector is all zeros
  int idx = 0;
  bool zerocol=false;
  Vec3 cv1=Vec3(m[0][0],m[1][0],m[2][0]);
  Vec3 cv2=Vec3(m[0][1],m[1][1],m[2][1]);
  Vec3 cv3=Vec3(m[0][2],m[1][2],m[2][2]);
  if(cv1.norm()<1e-7){
    idx=0;
    zerocol=true;
  }
  if(cv2.norm()<1e-7){
    idx=1;
    zerocol=true;
  }
  if(cv3.norm()<1e-7){
    idx=2;
    zerocol=true;
  }
  if(zerocol){
    res=Vec3(0.0,0.0,0.0);
    res[idx]=1.0;
  } else {
    // find non-zero element in 1st column
    int i=0;
    while((m[i][0]==0) && (i<3)) i++;
    // the other 2 rows indices
    int j=(i+1)%3;
    int k=(i+2)%3; 

    // copy data
    Matrix3 temp=*this;
    // zero out the 1st column in the other 2 rows
    double f=temp.m[j][0]/temp.m[i][0];
    temp.m[j][0]=0.0; 
    temp.m[j][1]-=f*temp.m[i][1];
    temp.m[j][2]-=f*temp.m[i][2];
    f=temp.m[k][0]/temp.m[i][0];
    temp.m[k][0]=0.0; 
    temp.m[k][1]-=f*temp.m[i][1];
    temp.m[k][2]-=f*temp.m[i][2];
//     cout << "step 1 : " << temp << endl;
    // check if row j is zero (can't have produced a zero column)
    double nj=Vec3(temp.m[j][0],temp.m[j][1],temp.m[j][2]).norm();
    if(nj<1e-7){
      // check if m[k][1] is zero 
      if(fabs(temp.m[k][1])<10e-7){ // diagonalize m[i][0],m[k][2]
	// zero m[i][2]
	f=temp.m[i][2]/temp.m[k][2];
	temp.m[i][2]=0.0;
	temp.m[i][1]-=f*temp.m[k][1];
	//      cout << "step 2a : " << temp << endl;
	// solve for res[0],res[1]
	res[0]=-1.0*temp.m[i][1]/temp.m[i][0];
	res[1]=1.0;
	res[2]=-1.0*temp.m[k][1]/temp.m[k][2];
      } else { // m[k][1]!=0.0 -> diagonalize m[i][0],m[k][1] 
	// zero m[i][1]
	f=temp.m[i][1]/temp.m[k][1];
	temp.m[i][1]=0.0;
	temp.m[i][2]-=f*temp.m[k][2];
	//       cout << "step 2b : " << temp << endl;
	// solve for res[0],res[1]
	res[0]=-1.0*temp.m[i][2]/temp.m[i][0];
	res[1]=1.0;
	res[2]=-1.0*temp.m[k][2]/temp.m[k][1];
      }
    } else {
      // check if m[j][1] is zero 
      if(fabs(temp.m[j][1])<10e-7){ // diagonalize m[i][0],m[j][2]
	// zero m[i][2]
	f=temp.m[i][2]/temp.m[j][2];
	temp.m[i][2]=0.0;
	temp.m[i][1]-=f*temp.m[j][1];
	//      cout << "step 2a : " << temp << endl;
	// solve for res[0],res[1]
	res[0]=-1.0*temp.m[i][1]/temp.m[i][0];
	res[1]=-1.0*temp.m[j][1]/temp.m[j][2];
	res[2]=1.0;
      } else { // m[j][1]!=0.0 -> diagonalize m[i][0],m[j][1] 
	// zero m[i][1]
	f=temp.m[i][1]/temp.m[j][1];
	temp.m[i][1]=0.0;
	temp.m[i][2]-=f*temp.m[j][2];
	//       cout << "step 2b : " << temp << endl;
	// solve for res[0],res[1]
	res[0]=-1.0*temp.m[i][2]/temp.m[i][0];
	res[1]=-1.0*temp.m[j][2]/temp.m[j][1];
	res[2]=1.0;
      }
    }
  }

//   cout << " res: " << res << endl;
  return res.unit();
}

/*!
  calculate eigenvectors & eigenvalues

  \returns v1 eigenvector with largest eigenvalue
  \returns v2 eigenvector with middle eigenvalue
  \returns v3 eigenvector with smalles eigenvalue
  \returns e1 largest eigenvalue
  \returns e2 middle eigenvalue
  \returns e3 smalles eigenvalue
*/
void Matrix3::eigen(Vec3& v1,Vec3& v2,Vec3& v3,double& e1,double& e2,double& e3)
{
  double te1,te2;
  // --- eigenvalues first ---
  // get determinant
  double I1=-1.0*(m[0][0]+m[1][1]+m[2][2]);
  double I2=m[0][0]*m[1][1]-m[0][1]*m[1][0]+
    m[0][0]*m[2][2]-m[0][2]*m[2][0]+
    m[1][1]*m[2][2]-m[1][2]*m[2][1];
  double I3=-1.0*det();
  //  std::cout << "I1,I2,I3 " << I1 << " , " << I2 << " , " << I3 << std::endl; 
  if(I3==0){ // simple case, matrix singluar -> 1 eigenvalue 0
    // solve quadratic eqn.
    double p=0.5*I1;
    double q=sqrt(p*p-I2);
    te1=p+q;
    te2=p-q;
    // calc eigenvectors
    Matrix3 m=*this;
    m+=Matrix3::Unit()*-1.0*te1;
    v1=(m.solve(Vec3(0.0,0.0,0.0))).unit();
    m=*this;
    m+=Matrix3::Unit()*-1.0*te2;
    v2=(m.solve(Vec3(0.0,0.0,0.0))).unit();
    v3=cross(v1,v2); 
  } else { // matrix not singular -> solve full cubic equation
    CubicEquation CE(I1,I2,I3);
    set<double> roots=CE.getRealRoots(10e-8);
    if(roots.size()==3){
      set<double>::iterator iter=roots.begin();
      e1=*iter;
      iter++;
      e2=*iter;
      iter++;
      e3=*iter;
      Matrix3 m=*this;
      m+=Matrix3::Unit()*-1.0*e1;
      v1=(m.solve_homogeneous()).unit();
      m=*this;
      m+=Matrix3::Unit()*-1.0*e2;
      v2=(m.solve_homogeneous()).unit();
      m=*this;
      m+=Matrix3::Unit()*-1.0*e3;
      v3=(m.solve_homogeneous()).unit();
    }else{ //can't happen
      // throw some error
    }
  }
}


