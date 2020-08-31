/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: QCDenvNODE.h
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

#include "QCDpt.h"
#include "lattice.h"
#include<iostream>
#include<fstream>

#ifdef _OPENMP
#include <omp.h>
#else
static const int omp_get_max_threads() { return 1; }
static const int omp_get_num_threads() { return 1; }
static const int omp_get_thread_num() { return 0; }
#endif


#define DIFF(a, b) (( fabs((a) - (b)) < 1e-15) ? 0 : 1)  
#define DIFF13(a, b) (( fabs((a) - (b)) < 1e-13) ? 0 : 1)  
#define ISDIFF(a, b, p) (( fabs((a) - (b)) < (p)) ? 0 : 1)

using namespace std;


class SU3_fld;
class Gluon_fld;
  
class SU2_fld;
class SU2Gluon_fld;
  
  
class SU3_fld{

 public:

  latt *Z;
  SU3  *W;

  SU3_fld(latt* z)  {
    Z = z;
    W = new SU3[Z->TSize];
  }

  ~SU3_fld(){
    delete [] W;
  }
  
  SU3_fld(FILE*,int); // SU3_fld(&input_file,read_mode)

  SU3* handle(){ return W; }

  SU3 get(int *);  
  
  SU3 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(SU3_fld *, int*);
  friend int get(SU3_fld *, int, int, int);
};  
  

  
inline int get(SU3_fld *W,int n,int step,int dir){
  return W->Z->get(n,step,dir);
}

  

  
class Gluon_fld{
 private:
  latt   *Z;

 public:
  Gluon  *W;
  Gluon_fld(latt *z){
    Z = z;
    W = new Gluon[z->TSize];
  }
  
  ~Gluon_fld(){
    delete [] W;
  }
  
  Gluon_fld(FILE*,int); // Gluon_fld(&input_file,read_mode)
  

  void operator=(const Gluon_fld& A) {
    Z = A.Z;
    for(int i = 0; i < Z->Size; i++)
      W[i] = A.W[i];
  }
  
  
  Gluon get(int n, int step, int dir) {return W[Z->get(n, step, dir)];}
  
  

  SU3 staple(int n, int mu, int nu){

  /*
    b ->
    ^    ^
    |    |
    n == a
    ^    ^
    |    |
    d -> c

   */
  
    return( W[Z->get(n, 1, mu)].U[nu]*
	    dag(W[Z->L[n][4]].U[nu]*
		W[Z->get(n, 1, nu)].U[mu]) 
	    +
	    dag(W[Z->get(n, -1, nu)].U[mu]*
		W[Z->get(n, 1, mu, -1, nu)].U[nu])*
	    W[Z->get(n, -1, nu)].U[nu]);
  }
  
};



// SU(2)

class SU2_fld{

 public:

  latt *Z;
  SU2  *W;

  SU2_fld(latt* z)  {
    Z = z;
    W = new SU2[Z->TSize];
  }

  ~SU2_fld(){
    delete [] W;
  }
  
  SU2_fld(FILE*,int); // SU2_fld(&input_file,read_mode)

  SU2* handle(){ return W; }

  SU2 get(int *);  
  
  SU2 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(SU2_fld *, int*);
  friend int get(SU2_fld *, int, int, int);
};  
  




inline int get(SU2_fld *W,int n,int step,int dir){
  return W->Z->get(n,step,dir);
}

  

  
class SU2Gluon_fld{
 private:
  latt   *Z;

 public:
  SU2Gluon  *W;
  SU2Gluon_fld(latt *z){
    Z = z;
    W = new SU2Gluon[z->TSize];
  }
  
  ~SU2Gluon_fld(){
    delete [] W;
  }
  
  SU2Gluon_fld(FILE*,int); // SU2Gluon_fld(&input_file,read_mode)
  

  void operator=(const SU2Gluon_fld& A) {
    Z = A.Z;
    for(int i = 0; i < Z->Size; i++)
      W[i] = A.W[i];
  }
  
  



  SU2Gluon get(int n, int step, int dir) {return W[Z->get(n, step, dir)];}
  
  SU2 staple(int n, int mu, int nu){

  /*
    b ->
    ^    ^
    |    |
    n == a
    ^    ^
    |    |
    d -> c

   */
  
    return( W[Z->get(n, 1, mu)].U[nu]*
	    dag(W[Z->L[n][4]].U[nu]*
		W[Z->get(n, 1, nu)].U[mu]) 
	    +
	    dag(W[Z->get(n, -1, nu)].U[mu]*
		W[Z->get(n, 1, mu, -1, nu)].U[nu])*
	    W[Z->get(n, -1, nu)].U[nu]);
  }
  
};





