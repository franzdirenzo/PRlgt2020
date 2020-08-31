/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: QCDpt.cc
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

#include<stdio.h>
#include"QCDpt.h"

#define SWAP(a,b) { tmp = a; a = b; b = tmp; }

int PTORD = allocORD;

SU3* ptSU3::handle(){
  return (SU3*)&ptU;
}

ptSU3& ptSU3::operator=(const ptSU3& A){
  flag  = A.flag;
    for(int i = 0; i < PTORD; i++){
      for(int j = 0; j < 9; j++){
        ptU[i].whr[j] = A.ptU[i].whr[j];
      }
    }
  return *this;
}


ptSU3 ptSU3::operator+(const ptSU3& A) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] + A.ptU[i];
  }
  B.flag = flag + A.flag;

  return B;
}


ptSU3 ptSU3::operator-(const ptSU3& A) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] - A.ptU[i];
  }
  
  B.flag = flag - A.flag;

  return B;
}

ptSU3& ptSU3::operator+=(const ptSU3& A) {

  for(int i = 0; i < PTORD; i++){
    ptU[i] += A.ptU[i];
  }

  flag += A.flag;
  return *this;
}


ptSU3& ptSU3::operator-=(const ptSU3& A){

  for(int i = 0; i < PTORD; i++){
    ptU[i] -= A.ptU[i];
  }
  flag -= A.flag;
  return *this;
}


ptSU3& ptSU3::operator*=(const ptSU3 &A){
  ptSU3 B;

  for(int i = 0; i < (PTORD-1); i++){
    for(int j = 0; j < (PTORD-1-i); j++){
      B.ptU[i+j+1] += ptU[j]*A.ptU[i];
    }
  }

  if( (flag == 0) && (A.flag == 0) ) {
    B.flag = 0;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
    }
  }
  else{
    B.flag = flag * A.flag;
  }

  *this = B;
  return *this;
}


void ptSU3::Tr(Cplx *tt){
  for(int i = 0; i < PTORD; i++){
    tt[i+1] = (ptU[i].whr[0] +
	       ptU[i].whr[4] +
	       ptU[i].whr[8]);
  }
  tt[0] = 3.*flag;
}

ptSU3& ptSU3::Trless(){
  Cplx z;
  for(int i = 0; i < PTORD; i++){
    z = D3*(ptU[i].whr[0] +
	    ptU[i].whr[4] +
	    ptU[i].whr[8]);
    ptU[i].whr[0] -= z;
    ptU[i].whr[4] -= z;
    ptU[i].whr[8] -= z;
  }
  flag = 0;
  return *this;
}
  

ptSU3 operator*(const Cplx& z, const ptSU3 &U) {
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = z*U.ptU[i];
  }
  B.flag = z*U.flag;
  return B;
}


ptSU3 ptSU3::operator*(const double& x) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]*x;
  }
  B.flag = x*flag;
  return B;
}

ptSU3 ptSU3::operator/(const double& x) const{
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]/x;
  }
  B.flag = flag/x;
  return B;
}


ptSU3& ptSU3::operator*=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] *= x;
  }
  flag *= x;
  return *this;
}

ptSU3& ptSU3::operator/=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] /= x;
  }
  flag /= x;
  return *this;
}


ptSU3 dag(const ptSU3& A){
  ptSU3 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = dag(A.ptU[i]);
  }
  
  B.flag.re =  A.flag.re;
  B.flag.im = -A.flag.im;
  return B;
}


ptSU3& ptSU3::reH(){
  Cplx tr;
  for(int i = 0; i < PTORD; i++){
    ptU[i] -= dag(ptU[i]);
    ptU[i] = .5*ptU[i];
    tr = ptU[i].Tr()*D3;
    ptU[i].whr[0] -= tr;
    ptU[i].whr[4] -= tr;
    ptU[i].whr[8] -= tr;
  }
  flag = 0;
  return *this;
}


void ptSU3::prout(){
  flag.prout();
  printf("\n");
  for(int i = 0; i < PTORD; i++){
    ptU[i].prout();
  }
}

  SU3 multiply_fixed_order( ptSU3& A, ptSU3& B, const int ord) {
    SU3 result;
    for(int i = 0; i < ord; i++) result += A.ptU[i]*B.ptU[ord-1-i];
    if( (B.flag != 0) || (A.flag != 0) )
      result += A.flag*B.ptU[ord] + A.ptU[ord]*B.flag;
    
    return result;
  }
  
  SU3 multiply_fixed_order( ptSU3& A, ptSU3& B,
                            ptSU3& C, const int ord) {
    SU3 result;
    
    for(int i = 0; i < ord; i++) result += A[i]*multiply_fixed_order(B,C,ord-1-i);
    if( A.flag != 0 )
      result += A.flag*multiply_fixed_order(B,C,ord)+A.ptU[ord]*B.flag*C.flag;

    return result;
    //    return std::move(result);
  }

  SU3 multiply_fixed_order( ptSU3& A, ptSU3& B,
                            ptSU3& C, ptSU3& D, const int ord) {
    SU3 result;
 
    for(int i = 0; i < ord; i++) result += A[i]*multiply_fixed_order(B,C,D,ord-1-i);
    if( A.flag != 0 )
      result += A.flag*multiply_fixed_order(B,C,D,ord)+A.ptU[ord]*B.flag*C.flag*D.flag;

    return result;
  }



// ---- end of ptSU3 --------------- //




ptGluon& ptGluon::operator=(const ptGluon& A) {
#if dim == 4
  U[0] = A.U[0];
  U[1] = A.U[1];
  U[2] = A.U[2];
  U[3] = A.U[3];
#else
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
#endif
  return *this;
}


ptGluon::ptGluon(const ptGluon& A) {
  for (int i=0; i < dim; i++)
    U[i] = A.U[i];
}



ptSU3 ptGluon::operator*(const ptGluon &A) const{
  ptSU3 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++){
    res += U[i]*A.U[i];
  }
#endif
  return res;
}


void ptGluon::prout(){
  for(int i = 0; i < dim; i++)
    U[i].prout();
}


ptGluon dag(const ptGluon &A){  
  ptGluon B;
#if dim == 4
    B.U[0] = dag(A.U[0]);
    B.U[1] = dag(A.U[1]);
    B.U[2] = dag(A.U[2]);
    B.U[3] = dag(A.U[3]);
#else
  for(int i = 0; i < dim; i++){
    B.U[i] = dag(A.U[i]);
  } 
#endif
  return B;
}


// ------------- end of ptGluon  -----------------//


// SU2

SU2* ptSU2::handle(){
  return (SU2*)&ptU;
}

ptSU2& ptSU2::operator=(const ptSU2& A){
  flag  = A.flag;
    for(int i = 0; i < PTORD; i++){
      for(int j = 0; j < 4; j++){
        ptU[i].whr[j] = A.ptU[i].whr[j];
      }
    }
  return *this;
}


ptSU2 ptSU2::operator+(const ptSU2& A) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] + A.ptU[i];
  }
  B.flag = flag + A.flag;

  return B;
}


ptSU2 ptSU2::operator-(const ptSU2& A) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i] - A.ptU[i];
  }
  
  B.flag = flag - A.flag;

  return B;
}

ptSU2& ptSU2::operator+=(const ptSU2& A) {

  for(int i = 0; i < PTORD; i++){
    ptU[i] += A.ptU[i];
  }

  flag += A.flag;
  return *this;
}


ptSU2& ptSU2::operator-=(const ptSU2& A){

  for(int i = 0; i < PTORD; i++){
    ptU[i] -= A.ptU[i];
  }
  flag -= A.flag;
  return *this;
}


ptSU2& ptSU2::operator*=(const ptSU2 &A){
  ptSU2 B;

  for(int i = 0; i < (PTORD-1); i++){
    for(int j = 0; j < (PTORD-1-i); j++){
      B.ptU[i+j+1] += ptU[j]*A.ptU[i];
    }
  }

  if( (flag == 0) && (A.flag == 0) ) {
    B.flag = 0;
    for(int i = 0; i < PTORD; i++){
      B.ptU[i] += A.flag*ptU[i] + flag*A.ptU[i];
    }
  }
  else{
    B.flag = flag * A.flag;
  }

  *this = B;
  return *this;
}


void ptSU2::Tr(Cplx *tt){
  for(int i = 0; i < PTORD; i++){
    tt[i+1] = (ptU[i].whr[0] +
	       ptU[i].whr[3]);
  }
  tt[0] = 2.*flag;
}

ptSU2& ptSU2::Trless(){
  Cplx z;
  for(int i = 0; i < PTORD; i++){
    z = 0.5*(ptU[i].whr[0] +
	     ptU[i].whr[3] );
    ptU[i].whr[0] -= z;
    ptU[i].whr[3] -= z;
  }
  flag = 0;
  return *this;
}
  

ptSU2 operator*(const Cplx& z, const ptSU2 &U) {
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = z*U.ptU[i];
  }
  B.flag = z*U.flag;
  return B;
}


ptSU2 ptSU2::operator*(const double& x) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]*x;
  }
  B.flag = x*flag;
  return B;
}

ptSU2 ptSU2::operator/(const double& x) const{
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = ptU[i]/x;
  }
  B.flag = flag/x;
  return B;
}


ptSU2& ptSU2::operator*=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] *= x;
  }
  flag *= x;
  return *this;
}

ptSU2& ptSU2::operator/=(const double& x){
  for(int i = 0; i < PTORD; i++){
    ptU[i] /= x;
  }
  flag /= x;
  return *this;
}


ptSU2 dag(const ptSU2& A){
  ptSU2 B;
  for(int i = 0; i < PTORD; i++){
    B.ptU[i] = dag(A.ptU[i]);
  }
  
  B.flag.re =  A.flag.re;
  B.flag.im = -A.flag.im;
  return B;
}


ptSU2& ptSU2::reH(){
  Cplx tr;
  for(int i = 0; i < PTORD; i++){
    ptU[i] -= dag(ptU[i]);
    ptU[i] = .5*ptU[i];
    tr = ptU[i].Tr()*0.5;
    ptU[i].whr[0] -= tr;
    ptU[i].whr[3] -= tr;
  }
  flag = 0;
  return *this;
}


void ptSU2::prout(){
  flag.prout();
  printf("\n");
  for(int i = 0; i < PTORD; i++){
    ptU[i].prout();
  }
}

  SU2 multiply_fixed_order( ptSU2& A, ptSU2& B, const int ord) {
    SU2 result;
    for(int i = 0; i < ord; i++) result += A.ptU[i]*B.ptU[ord-1-i];
    if( (B.flag != 0) || (A.flag != 0) )
      result += A.flag*B.ptU[ord] + A.ptU[ord]*B.flag;
    
    return result;
  }
  
  SU2 multiply_fixed_order( ptSU2& A, ptSU2& B,
                            ptSU2& C, const int ord) {
    SU2 result;
    
    for(int i = 0; i < ord; i++) result += A[i]*multiply_fixed_order(B,C,ord-1-i);
    if( A.flag != 0 )
      result += A.flag*multiply_fixed_order(B,C,ord)+A.ptU[ord]*B.flag*C.flag;

    return result;
    //    return std::move(result);
  }

  SU2 multiply_fixed_order( ptSU2& A, ptSU2& B,
                            ptSU2& C, ptSU2& D, const int ord) {
    SU2 result;
 
    for(int i = 0; i < ord; i++) result += A[i]*multiply_fixed_order(B,C,D,ord-1-i);
    if( A.flag != 0 )
      result += A.flag*multiply_fixed_order(B,C,D,ord)+A.ptU[ord]*B.flag*C.flag*D.flag;

    return result;
  }



// ---- end of ptSU2 --------------- //




ptSU2Gluon& ptSU2Gluon::operator=(const ptSU2Gluon& A) {
#if dim == 4
  U[0] = A.U[0];
  U[1] = A.U[1];
  U[2] = A.U[2];
  U[3] = A.U[3];
#else
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
#endif
  return *this;
}


ptSU2Gluon::ptSU2Gluon(const ptSU2Gluon& A) {
  for (int i=0; i < dim; i++)
    U[i] = A.U[i];
}



ptSU2 ptSU2Gluon::operator*(const ptSU2Gluon &A) const{
  ptSU2 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++){
    res += U[i]*A.U[i];
  }
#endif
  return res;
}


void ptSU2Gluon::prout(){
  for(int i = 0; i < dim; i++)
    U[i].prout();
}


ptSU2Gluon dag(const ptSU2Gluon &A){  
  ptSU2Gluon B;
#if dim == 4
    B.U[0] = dag(A.U[0]);
    B.U[1] = dag(A.U[1]);
    B.U[2] = dag(A.U[2]);
    B.U[3] = dag(A.U[3]);
#else
  for(int i = 0; i < dim; i++){
    B.U[i] = dag(A.U[i]);
  } 
#endif
  return B;
}


// ------------- end of ptSU2Gluon  -----------------//


