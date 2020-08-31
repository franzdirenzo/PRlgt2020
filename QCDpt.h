/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: QCDpt.h
Copyright (C) 2020

Author: Michele Brambilla <michele.brambilla@psi.ch>
Author: Francesco Di Renzo <francesco.direnzo@unipr.it>
Author: Gianluca Filaci <gianluca.filaci@monozukuri.eu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL NOTICE */

#ifndef _MY_QCD_H_
#define _MY_QCD_H_

#include"MyQCD.h"

class ptSU3;
class ptCVector;
class ptGluon;

class ptSU2;
class ptCSU2Vector;
class ptSU2Gluon;

template<class, int> class serie;


class ptSU3{
 private:

 public:
    
  Cplx flag;    // 0 = matrice nell'algebra, 1 = vuoto banale, z3 ...  
  SU3 ptU[allocORD];

  ptSU3(Cplx i = VACUUM){
    flag = i;
  }
    
  ~ptSU3(){};
 
  SU3* handle();

  int write(FILE *filept) {
    if(fwrite(&flag, SZ_DB*2, 1, filept) 
       && fwrite(&ptU, SZ_DB*2*NC*NC, allocORD, filept) ) return 0;
    return 1;
  }

  int read(FILE *filept) {
    if(fread(&flag, SZ_DB*2, 1, filept) 
       && fread(&ptU, SZ_DB*2*NC*NC, allocORD, filept) ) return 0;
    return 1;
  }

  void zero(){
    memset(ptU, 0, SZ_DB*2*NC*NC*PTORD);
    flag  = 0;
  };

  void id(){
    memset(ptU, 0, SZ_DB*2*NC*NC*PTORD);
    flag  = VACUUM;
  };

  void Tr(Cplx *tt);  

  ptSU3& Trless();
  
  SU3& operator[] (int& ord) { return ptU[ord]; }

  ptSU3& operator=(const ptSU3&);

  ptSU3 operator+(const ptSU3&) const ;

  ptSU3 operator-(const ptSU3&) const;

  ptSU3 operator*(const ptSU3& A) const{
    ptSU3 B;

    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	B.ptU[i+j+1] += ptU[j]*A.ptU[i];
      }
    }
  
    if( (flag == 0) && (A.flag == 0) ) {
      B.flag = 0;
      return B;
    }
    else {

      B.flag = flag * A.flag;
      for(int i = 0; i < PTORD; i++){
	B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
      }
          
      return B;
    }
  };
  

  ptSU3& operator+=(const ptSU3&);

  ptSU3& operator-=(const ptSU3&);

  ptSU3& operator*=(const ptSU3& A);

  ptSU3 operator*(const Cplx& z) const{
    ptSU3 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]*z;
    }
    B.flag = z*flag;
    return (B);
  };
  
  ptSU3 operator/(const Cplx &z) const{
    ptSU3 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]/z;
    }
    B.flag = flag/z;
    return B;
  };
  
  ptSU3& operator*=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] *= z;
    }
    flag *= z;
    return *this;
  };
  
  ptSU3& operator/=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] /= z;
    }
    flag /= z;
    return *this;
  };

  inline friend ptSU3 operator*(const Cplx&, const ptSU3&);

  ptSU3 operator*(const double&) const ;
  ptSU3 operator/(const double&) const ;

  ptSU3& operator*=(const double&);
  ptSU3& operator/=(const double&);

  inline friend ptSU3 operator*(const double& x, const ptSU3& U){
    ptSU3 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = x*U.ptU[i];
    }
    B.flag = x*U.flag;
    return B;
  }

  friend ptSU3 dag(const ptSU3&);

  inline friend ptSU3 log(const ptSU3& U){
    ptSU3 B = U, Bnew, res = U;
    double segno = -1, aux;
    B.flag = 0;

    for(int i = 2; i <= PTORD; i++){
     Bnew.zero();
     for(int iU = 1; iU < PTORD; iU++){
       for(int iB = i-1; iB <= PTORD - iU; iB++){
 	Bnew.ptU[iU+iB-1] += B.ptU[iB-1]*U.ptU[iU-1];
       }
     }
     for(int eq = 0; eq < PTORD; eq++){
       B.ptU[eq] = Bnew.ptU[eq];
     }
     aux = segno/(double)i;
     res += aux*B;
     segno = -segno;
    }
    res.flag = 0;
    
    return res;
  }
  
  inline friend ptSU3 exp(const ptSU3& A){
    ptSU3 B = A, Bnew, res = A;
    double den;

    for(int i = 2; i <= PTORD; i++){
      
      Bnew.zero();
      for(int iA = 1; iA < PTORD; iA++){
        for(int iB = i-1; iB <= PTORD - iA; iB++){
	  Bnew.ptU[iA+iB-1] += B.ptU[iB-1]*A.ptU[iA-1];
        }
      }
      
      den = 1./(double)i;
      
      B = den*Bnew;
      res += B;
    }
    res.flag = 1;
    return res;
 }

  ptSU3& reH();

  void prout();    

};



SU3 multiply_fixed_order(const ptSU3&,const ptSU3&, const int);
  
SU3 multiply_fixed_order(const ptSU3&,const ptSU3&,
                          const ptSU3&, const int);

SU3 multiply_fixed_order( const ptSU3&, const ptSU3&,
                          const ptSU3&, const ptSU3&, const int);


#define ptSU3U ptSU3(1.0)
#define ptSU3A ptSU3(0.0)



// --- end of ptSU3 declarations ---




class ptGluon {

 public:
  ptSU3 U[dim];

  ptGluon(){
    /* for (int mu = 0; mu < dim; mu++) */
    /*   U[mu].flag = 1; */
  };

  ptGluon(const ptGluon& A);

  int write(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].write(filept)) return 0;
    return 1;
  }
  
  int read(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].read(filept)) return 0;
    return 1;
  }

  ptGluon operator*(const Cplx& z) const
    {
      ptGluon res;
      for( int mu = 0; mu < dim; ++mu)
	res.U[mu] = U[mu]*z;
      return res;
    }

  ptGluon& operator*=(const Cplx& z) 
    {
      for( int mu = 0; mu < dim; ++mu)
	U[mu] *= z;
      return *this;
    }


  ptGluon& operator=(const ptGluon& A);
  
  ptSU3 operator*(const ptGluon& A) const;

  void prout();
  void zero() {
    for( int mu = 0; mu < dim; mu++)
      U[mu].zero();
  }

  friend ptGluon dag(const ptGluon& A);

};

// --- end of ptGluon declarations ---



class ptCVector {

public:

  CVector ptCV[allocORD + 1];

  int write(FILE *filept) {
    if(fwrite(&ptCV, SZ_DB*2*NC, PTORD+1, filept) ) return 0;
    return 1;
  }

  int read(FILE *filept) {
    if(fread(&ptCV, SZ_DB*2*NC, PTORD+1, filept) ) return 0;
    return 1;
  }
  
  void prout(){
    for (int i = 0; i <= PTORD; i++){
      ptCV[i].prout();
      printf("\n");
    }
  }

  ptCVector& operator=(const ptCVector &ptcv) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] = ptcv.ptCV[i];
    return *this;
  }
  
  ptCVector operator+(const ptCVector &cv0) const{
    ptCVector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = ptCV[i] + cv0.ptCV[i];
    return cv1;
  }
  
  ptCVector operator-(const ptCVector &cv0) const{
    ptCVector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = ptCV[i] - cv0.ptCV[i];
    return cv1;
  }
  
  ptCVector& operator+=(const ptCVector &cv0){
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] += cv0.ptCV[i];
    return *this;
  }
  
  ptCVector& operator-=(const ptCVector &cv0) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] += cv0.ptCV[i];
    return *this;
  }
  
  ptCVector operator*(const Cplx& z) const{
    ptCVector cv;    
    for (int i = 0; i <= PTORD; i++)
      cv.ptCV[i] = ptCV[i] * z;
    return cv;
  }

  ptCVector operator/(const Cplx& z) const{
    ptCVector cv;
    for (int i = 0; i <= PTORD; i++)
      cv.ptCV[i] = ptCV[i] / z;
    return cv;
  }

  ptCVector& operator*=(const Cplx& z) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] *= z;
    return *this;
  }
  
  ptCVector& operator/=(const Cplx& z) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] = ptCV[i] / z;
    return *this;
  }
  
  friend ptCVector operator*(const Cplx& z, const ptCVector& cv0) {
    ptCVector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = cv0.ptCV[i] * z;
    return cv1;
  }
  


  ptCVector operator*(ptSU3& A) const{
    ptCVector cv;
    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	cv.ptCV[i+j+1] += ptCV[j]*A.ptU[i];
      }
    }
    if(A.flag != 0) {
      for(int i = 0; i < PTORD; i++){
	cv.ptCV[i] += A.flag*ptCV[i];
      }
    }
    return cv;
  }

  ptCVector operator-() const{
    ptCVector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = -ptCV[i];
    return cv1;
  }
  

  friend ptCVector operator*(ptSU3 &A, ptCVector& cv0){
    ptCVector cv1;
    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	cv1.ptCV[i+j+1] += A.ptU[j]*cv0.ptCV[i];
      }
    }
    if(A.flag != 0) {
      for(int i = 0; i < PTORD; i++){
	cv1.ptCV[i] += A.flag*cv0.ptCV[i];
      }
    }
    
    return cv1;
  };
  
};


// ---- end ptCVector declarations -------


/***********   generic series   *************/



template<class T, int ptord>
class serie{

 private:
  
  
 public:
  
  T ord[1+ptord];

  serie() {};
  
  ~serie() {};

  serie& operator=(const serie& M) 
  { 
    for( int i = 0; i <= ptord; ++i) 
      ord[i] = M.ord[i]; 
    return *this;
  }
  
  T operator[](int n)    {  return ord[n];  }

  serie& set(T& X, int oo)
  {
    ord[oo] = X;
    return *this;
  }

  serie operator+(const serie& M) 
  { 
    serie<T,ptord> res; 
    for( int i = 0; i <= ptord; ++i) 
      res.ord[i] = ord[i] + M.ord[i]; 
    return res;
  }
  
  serie operator-(const serie& M) 
  { 
    serie<T,ptord> res; 
    for( int i = 0; i <= ptord; ++i) 
      res.ord[i] = ord[i] - M.ord[i]; 
    return res;
  }
  
  serie operator*(const serie& M) 
  { 
    serie<T,ptord> res; 
    
    for(int i = 0; i <= ptord; ++i)
      {
	res.ord[i] += ord[0]*M.ord[i];
	for(int j = 1; j <= ptord-i; ++j){
	  res.ord[i+j] += ord[j]*M.ord[i];
	}
      }

    return res;
  }
  
  serie& operator+=(const serie& M) 
  { 
    for( int i = 0; i <= ptord; ++i) 
      ord[i] += M.ord[i]; 
    return *this;
  }
  
  serie& operator-=(const serie& M) 
  { 
    for( int i = 0; i <= ptord; ++i) 
      ord[i] -= M.ord[i]; 
    return *this;
  }
  
  serie& operator*=(const serie& M) 
  { 
    serie<T,ptord> res; 

    for(int i = 0; i <= ptord; ++i)
      {
	res.ord[i] += ord[0]*M.ord[i];
	for(int j = 1; j <= ptord-i; ++j){
	  res.ord[i+j] += ord[j]*M.ord[i];
	}
      }
    for(int i = 0; i <= ptord; ++i) ord[i] = res.ord[i];

    return *this;
  }

  serie& operator%(const serie& M) 
  { 
    for(int i = 0; i <= ptord; ++i)
      {
	ord[i] *= M.ord[i];
      }
    return *this;
  }

  serie operator%=(const serie& M) 
  { 
    serie<T,ptord> res;
    for(int i = 0; i <= ptord; ++i)
      {
	res.ord[i] = ord[i]*M.ord[i];
      }
    return res;
  }


  template<class X> serie& operator=(const X& M)   {   ord[0] = M;    return *this;  }
  
  template<class X> serie operator*(const X& M) 
  { 
    serie<T,ptord> res; 
    for(int i = 0; i <= ptord; ++i)   res.ord[i] = ord[i]*M;
    return res;
  }
  
  template<class X> serie& operator*=(const X& M) 
  { 
    serie<T,ptord> res; 
    for(int i = 0; i <= ptord; ++i)      ord[i] *= M;
    return *this;
  }

  template<class X> serie operator/(const X& M) 
  { 
    X tmp = 1.0/M;
    serie<T,ptord> res; 
    for(int i = 0; i <= ptord; ++i)   res.ord[i] = ord[i]*tmp;
    return res;
  }
  
  template<class X> serie& operator/=(const X& M) 
  { 
    X tmp = 1.0/M;
    serie<T,ptord> res; 
    for(int i = 0; i <= ptord; ++i)      ord[i] *= tmp;
    return *this;
  }

  // perturbative inversion
  serie& invert(serie& M, const T& I0)
    {
      ord[0] = I0;
      
      for(int iO = 1; iO <= ptord; ++iO)
	{
	  for(int jO = 0; jO < iO; ++jO)
	    {
	      ord[iO] += M.ord[iO-jO]*ord[jO];
	    }
	  ord[iO] = -I0*ord[iO];
	}

      return *this;
    }


  // perturbative inversion
  serie& invert2(serie& M, const T& I0)
    {
      ord[0] = I0;
      
      for(int iO = 2; iO <= ptord; iO+=2)
	{
	  for(int jO = 0; jO < iO; jO+=2)
	    {
	      ord[iO] += M.ord[iO-jO]*ord[jO];
	    }
	  ord[iO] = -I0*ord[iO];
	}
      for(int iO = 1; iO <= ptord/2; ++iO)  ord[iO] = ord[2*iO];

      return *this;
    }
  

  template<class X, int N> inline friend std::ostream& operator<<(std::ostream&, serie<X,N>);

};




template<class T,int ptord> inline std::ostream& operator<<(std::ostream& out, serie<T,ptord> p)
{
  for( int oo = 0; oo <= ptord; ++oo)
    {
      out << "ord = "  << oo         << std::endl << std::endl;
      out << p.ord[oo] << std::endl  << std::endl;
    }
  return out;
}


/*********** end generic series *************/


/***SU(2)****/

class ptSU2{
 private:

 public:

  Cplx flag;    // 0 = matrice nell'algebra, 1 = vuoto banale, z3 ...  
  SU2 ptU[allocORD];

  ptSU2(Cplx i = VACUUM){
    flag = i;
  }
  
  ~ptSU2(){};
 
  SU2* handle();

  int write(FILE *filept) {
    if(fwrite(&flag, SZ_DB*2, 1, filept) 
       && fwrite(&ptU, SZ_DB*2*2*2, allocORD, filept) ) return 0;
    return 1;
  }

  int read(FILE *filept) {
    if(fread(&flag, SZ_DB*2, 1, filept) 
       && fread(&ptU, SZ_DB*2*2*2, allocORD, filept) ) return 0;
    return 1;
  }

  void zero(){
    memset(ptU, 0, SZ_DB*2*2*2*PTORD);
    flag  = 0;
  };

  void id(){
    memset(ptU, 0, SZ_DB*2*2*2*PTORD);
    flag  = VACUUM;
  };

  void Tr(Cplx *tt);  

  ptSU2& Trless();
  
  SU2& operator[] (int& ord) { return ptU[ord]; }

  ptSU2& operator=(const ptSU2&);

  ptSU2 operator+(const ptSU2&) const ;

  ptSU2 operator-(const ptSU2&) const;

  ptSU2 operator*(const ptSU2& A) const{
    ptSU2 B;
    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	B.ptU[i+j+1] += ptU[j]*A.ptU[i];
      }
    }
  
    if( (flag == 0) && (A.flag == 0) ) {
      B.flag = 0;
      return B;
    }
    else {

      B.flag = flag * A.flag;
      for(int i = 0; i < PTORD; i++){
	B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
      }
          
      return B;
    }
  };
  

  ptSU2& operator+=(const ptSU2&);

  ptSU2& operator-=(const ptSU2&);

  ptSU2& operator*=(const ptSU2& A);

  ptSU2 operator*(const Cplx& z) const{
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]*z;
    }
    B.flag = z*flag;
    return (B);
  };
  
  ptSU2 operator/(const Cplx &z) const{
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = ptU[i]/z;
    }
    B.flag = flag/z;
    return B;
  };
  
  ptSU2& operator*=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] *= z;
    }
    flag *= z;
    return *this;
  };
  
  ptSU2& operator/=(const Cplx &z){
    for(int i = 0; i < PTORD; i++){
      ptU[i] /= z;
    }
    flag /= z;
    return *this;
  };

  inline friend ptSU2 operator*(const Cplx&, const ptSU2&);

  ptSU2 operator*(const double&) const ;
  ptSU2 operator/(const double&) const ;

  ptSU2& operator*=(const double&);
  ptSU2& operator/=(const double&);

  inline friend ptSU2 operator*(const double& x, const ptSU2& U){
    ptSU2 B;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] = x*U.ptU[i];
    }
    B.flag = x*U.flag;
    return B;
  }

  friend ptSU2 dag(const ptSU2&);

  inline friend ptSU2 log(const ptSU2& U){
    ptSU2 B = U, Bnew, res = U;
    double segno = -1, aux;
    B.flag = 0;

    for(int i = 2; i <= PTORD; i++){
     Bnew.zero();
     for(int iU = 1; iU < PTORD; iU++){
       for(int iB = i-1; iB <= PTORD - iU; iB++){
 	Bnew.ptU[iU+iB-1] += B.ptU[iB-1]*U.ptU[iU-1];
       }
     }
     for(int eq = 0; eq < PTORD; eq++){
       B.ptU[eq] = Bnew.ptU[eq];
     }
     aux = segno/(double)i;
     res += aux*B;
     segno = -segno;
    }
    res.flag = 0;
    
    return res;
  }
  
  inline friend ptSU2 exp(const ptSU2& A){
    ptSU2 B = A, Bnew, res = A;
    double den;

    for(int i = 2; i <= PTORD; i++){
      
      Bnew.zero();
      for(int iA = 1; iA < PTORD; iA++){
        for(int iB = i-1; iB <= PTORD - iA; iB++){
	  Bnew.ptU[iA+iB-1] += B.ptU[iB-1]*A.ptU[iA-1];
        }
      }
      
      den = 1./(double)i;
      
      B = den*Bnew;
      res += B;
    }
    res.flag = 1;
    return res;
 }

  ptSU2& reH();

  void prout();    

};


SU2 multiply_fixed_order(const ptSU2&,const ptSU2&, const int);
  
SU2 multiply_fixed_order(const ptSU2&,const ptSU2&,
                          const ptSU2&, const int);

SU2 multiply_fixed_order( const ptSU2&, const ptSU2&,
                          const ptSU2&, const ptSU2&, const int);


#define ptSU2U ptSU2(1.0)
#define ptSU2A ptSU2(0.0)



// --- end of ptSU2 declarations ---




class ptSU2Gluon {

 public:
  ptSU2 U[dim];  

  ptSU2Gluon(){
    /* for (int mu = 0; mu < dim; mu++) */
    /*   U[mu].flag = 1; */
  };

  ptSU2Gluon(const ptSU2Gluon& A);

  int write(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].write(filept)) return 0;
    return 1;
  }
  
  int read(FILE *filept) {
    for (int mu = 0; mu < dim; mu++)
      if (U[mu].read(filept)) return 0;
    return 1;
  }

  ptSU2Gluon operator*(const Cplx& z) const
    {
      ptSU2Gluon res;
      for( int mu = 0; mu < dim; ++mu)
	res.U[mu] = U[mu]*z;
      return res;
    }

  ptSU2Gluon& operator*=(const Cplx& z) 
    {
      for( int mu = 0; mu < dim; ++mu)
	U[mu] *= z;
      return *this;
    }


  ptSU2Gluon& operator=(const ptSU2Gluon& A);
  
  ptSU2 operator*(const ptSU2Gluon& A) const;

  void prout();
  void zero() {
    for( int mu = 0; mu < dim; mu++)
      U[mu].zero();
  }

  friend ptSU2Gluon dag(const ptSU2Gluon& A);

};

// --- end of ptSU2Gluon declarations ---



class ptCSU2Vector {

public:

  CSU2Vector ptCV[allocORD + 1];

  int write(FILE *filept) {
    if(fwrite(&ptCV, SZ_DB*2*2, PTORD+1, filept) ) return 0;
    return 1;
  }

  int read(FILE *filept) {
    if(fread(&ptCV, SZ_DB*2*2, PTORD+1, filept) ) return 0;
    return 1;
  }
  
  void prout(){
    for (int i = 0; i <= PTORD; i++){
      ptCV[i].prout();
      printf("\n");
    }
  }

  ptCSU2Vector& operator=(const ptCSU2Vector &ptcv) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] = ptcv.ptCV[i];
    return *this;
  }
  
  ptCSU2Vector operator+(const ptCSU2Vector &cv0) const{
    ptCSU2Vector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = ptCV[i] + cv0.ptCV[i];
    return cv1;
  }
  
  ptCSU2Vector operator-(const ptCSU2Vector &cv0) const{
    ptCSU2Vector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = ptCV[i] - cv0.ptCV[i];
    return cv1;
  }
  
  ptCSU2Vector& operator+=(const ptCSU2Vector &cv0){
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] += cv0.ptCV[i];
    return *this;
  }
  
  ptCSU2Vector& operator-=(const ptCSU2Vector &cv0) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] += cv0.ptCV[i];
    return *this;
  }
  
  ptCSU2Vector operator*(const Cplx& z) const{
    ptCSU2Vector cv;    
    for (int i = 0; i <= PTORD; i++)
      cv.ptCV[i] = ptCV[i] * z;
    return cv;
  }

  ptCSU2Vector operator/(const Cplx& z) const{
    ptCSU2Vector cv;
    for (int i = 0; i <= PTORD; i++)
      cv.ptCV[i] = ptCV[i] / z;
    return cv;
  }

  ptCSU2Vector& operator*=(const Cplx& z) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] *= z;
    return *this;
  }
  
  ptCSU2Vector& operator/=(const Cplx& z) {
    for (int i = 0; i <= PTORD; i++)
      ptCV[i] = ptCV[i] / z;
    return *this;
  }
  
  friend ptCSU2Vector operator*(const Cplx& z, const ptCSU2Vector& cv0) {
    ptCSU2Vector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = cv0.ptCV[i] * z;
    return cv1;
  }
  


  ptCSU2Vector operator*(ptSU2& A) const{
    ptCSU2Vector cv;
    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	cv.ptCV[i+j+1] += ptCV[j]*A.ptU[i];
      }
    }
    if(A.flag != 0) {
      for(int i = 0; i < PTORD; i++){
	cv.ptCV[i] += A.flag*ptCV[i];
      }
    }
    return cv;
  }

  ptCSU2Vector operator-() const{
    ptCSU2Vector cv1;
    for (int i = 0; i <= PTORD; i++)
      cv1.ptCV[i] = -ptCV[i];
    return cv1;
  }
  

  friend ptCSU2Vector operator*(ptSU2 &A, ptCSU2Vector& cv0){
    ptCSU2Vector cv1;
    for(int i = 0; i < (PTORD-1); i++){
      for(int j = 0; j < (PTORD-1-i); j++){
	cv1.ptCV[i+j+1] += A.ptU[j]*cv0.ptCV[i];
      }
    }
    if(A.flag != 0) {
      for(int i = 0; i < PTORD; i++){
	cv1.ptCV[i] += A.flag*cv0.ptCV[i];
      }
    }
    
    return cv1;
  };
  
};


// ---- end ptCSU2Vector declarations -------




inline std::ostream& operator<<(std::ostream& out, ptSU3 U)
{
  out << U.flag << std::endl << std::endl;
  for( int oo = 0; oo < PTORD; ++oo)
    {
      out << "ord = "  << oo         << std::endl << std::endl;
      out << U.ptU[oo] << std::endl  << std::endl;
    }
  return out;
}


inline std::ostream& operator<<(std::ostream& out, ptGluon G)
{

  for( int mu = 0; mu < dim; ++mu)
    {
      out << "\t\t[ mu = "  << mu << "]" << std::endl << std::endl;
      out << G.U[mu] << std::endl  << std::endl;
    }
  return out;
}


inline std::ostream& operator<<(std::ostream& out, ptCVector P)
{
  for( int oo = 0; oo <= PTORD; ++oo)
    {
      out << "ord = "  << oo         << std::endl << std::endl;
      out << P.ptCV[oo] << std::endl  << std::endl;
    }
  return out;
}



inline std::ostream& operator<<(std::ostream& out, ptSU2 U)
{
  out << U.flag << std::endl << std::endl;
  for( int oo = 0; oo < PTORD; ++oo)
    {
      out << "ord = "  << oo         << std::endl << std::endl;
      out << U.ptU[oo] << std::endl  << std::endl;
    }
  return out;
}


inline std::ostream& operator<<(std::ostream& out, ptSU2Gluon G)
{

  for( int mu = 0; mu < dim; ++mu)
    {
      out << "\t\t[ mu = "  << mu << "]" << std::endl << std::endl;
      out << G.U[mu] << std::endl  << std::endl;
    }
  return out;
}


inline std::ostream& operator<<(std::ostream& out, ptCSU2Vector P)
{
  for( int oo = 0; oo <= PTORD; ++oo)
    {
      out << "ord = "  << oo         << std::endl << std::endl;
      out << P.ptCV[oo] << std::endl  << std::endl;
    }
  return out;
}



#endif
