/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: MyQCD.h
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

#ifndef _MY_QCD_H
#define _MY_QCD_H

#include"MyMath.h"
#include"choices.h"

class Gluon;
class SU2Gluon;

class Gluon {

 public:
  SU3 U[dim];  

  Gluon(){};

  Gluon(const Gluon& A) {
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
  }

  int write(FILE *filept){
    if(fwrite(&U, SZ_DB*NC*NC*2, dim, filept)) return 0;
    return 1;
  }

  int read(FILE *filept){
    if(fread(&U, SZ_DB*9*2, dim, filept)) return 0;
    return 1;
  }

  Gluon& operator=(const Gluon& A) {
    for (int i=0; i < dim; i++){
      U[i] = A.U[i];
    }
    return *this;
  };

  Gluon operator-() {
    Gluon res;
    for (int i=0; i < dim; i++){
      res.U[i] = -U[i];
    }
    return res;
  };
  
  SU3 operator*(const Gluon&) const;

  Gluon operator*(const Cplx&) const;

  void dag() {  for(int i = 0; i < dim; i++) U[i].dag(); }

  void prout(){
    for(int i = 0; i < dim; i++)
      U[i].prout();
  };

};

Gluon dag(const Gluon &U);


/************* end class gluon *****************/
/***********************************************/


class SU2Gluon {

 public:
  SU2 U[dim];  

  SU2Gluon(){};

  SU2Gluon(const SU2Gluon& A) {
    for (int i=0; i < dim; i++)
      U[i] = A.U[i];
  }

  int write(FILE *filept){
    if(fwrite(&U, SZ_DB*4*2, dim, filept)) return 0;
    return 1;
  }

  int read(FILE *filept){
    if(fread(&U, SZ_DB*4*2, dim, filept)) return 0;
    return 1;
  }

  SU2Gluon& operator=(const SU2Gluon& A) {
    for (int i=0; i < dim; i++){
      U[i] = A.U[i];
    }
    return *this;
  };

  SU2Gluon operator-() {
    SU2Gluon res;
    for (int i=0; i < dim; i++){
      res.U[i] = -U[i];
    }
    return res;
  };
  
  SU2 operator*(const SU2Gluon&) const;

  SU2Gluon operator*(const Cplx&) const;

  void dag() {  for(int i = 0; i < dim; i++) U[i].dag(); }

  void prout(){
    for(int i = 0; i < dim; i++)
      U[i].prout();
  };

};

SU2Gluon dag(const SU2Gluon &U);


/************* end class gluon *****************/
/***********************************************/




inline std::ostream& operator<<(std::ostream& out, Gluon G)
{

  for( int mu = 0; mu < dim; ++mu)
    {
      out << "\t\t[ mu = "  << mu << "]" << std::endl << std::endl;
      out << G.U[mu] << std::endl  << std::endl;
    }
  return out;
}


inline std::ostream& operator<<(std::ostream& out, SU2Gluon G)
{

  for( int mu = 0; mu < dim; ++mu)
    {
      out << "\t\t[ mu = "  << mu << "]" << std::endl << std::endl;
      out << G.U[mu] << std::endl  << std::endl;
    }
  return out;
}



#endif
