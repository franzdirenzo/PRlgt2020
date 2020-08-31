/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: MyQCD.cc
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

#include"MyQCD.h"
#define SWAP(a,b) { CVector tmp = a; a = b; b = tmp; }

SU3 Gluon::operator*(const Gluon &A) const{
  SU3 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++) res += U[i]*A.U[i];
#endif
  return res;
}


Gluon Gluon::operator*(const Cplx& z) const{
  Gluon res;
#if dim == 4
  res.U[0] = U[0]*z;
  res.U[1] = U[1]*z;
  res.U[2] = U[2]*z;
  res.U[3] = U[3]*z;
#else
  for(int i = 0; i < dim; i++) res.psi[i] = U[i]*z;
#endif
  return res;
}


Gluon dag(const Gluon &A) {  
  Gluon res;
#if dim == 4
  res.U[0] = dag(A.U[0]); 
  res.U[1] = dag(A.U[1]); 
  res.U[2] = dag(A.U[2]); 
  res.U[3] = dag(A.U[3]); 
#else
  for(int i = 0; i < dim; i++) res.U[i] = dag(A.U[i]);
#endif
  return res;
}



SU2 SU2Gluon::operator*(const SU2Gluon &A) const{
  SU2 res;
#if dim == 4
  res = U[0]*A.U[0] + U[1]*A.U[1] + U[2]*A.U[2] + U[3]*A.U[3] ;
#else
  for(int i = 0; i < dim; i++) res += U[i]*A.U[i];
#endif
  return res;
}


SU2Gluon SU2Gluon::operator*(const Cplx& z) const{
  SU2Gluon res;
#if dim == 4
  res.U[0] = U[0]*z;
  res.U[1] = U[1]*z;
  res.U[2] = U[2]*z;
  res.U[3] = U[3]*z;
#else
  for(int i = 0; i < dim; i++) res.psi[i] = U[i]*z;
#endif
  return res;
}


SU2Gluon dag(const SU2Gluon &A) {  
  SU2Gluon res;
#if dim == 4
  res.U[0] = dag(A.U[0]); 
  res.U[1] = dag(A.U[1]); 
  res.U[2] = dag(A.U[2]); 
  res.U[3] = dag(A.U[3]); 
#else
  for(int i = 0; i < dim; i++) res.U[i] = dag(A.U[i]);
#endif
  return res;
}


