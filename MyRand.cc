/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: MyRand.cc
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

#include "MyRand.h"
#include <stdio.h>
#include <string.h>

long Seed = 161803398L;

MyRand::MyRand(long idum)
{
  gen = 0;
#ifdef __GEN_RAND55__
    int i,j;
    long tmp,aux;

    for(i=0;i<54;++i) {
      next[i] = i+1;
    }
    next[54]= 0;
    tmp = Seed + idum;
    ring[54] = tmp;
    aux = 1;
    j = 20;
    for(i=0;i<54;++i,j+=21){
	j %= 55;
	ring[j] = aux;
	aux += tmp;
	tmp = ring[j];
    }

    n1 = 0;
    n2 = 31;

    for(j=0;j<55*4;++j) {
      Ranf();
    }

    gss_flag = 0;
    gss_value = 0.0;
    gss_sigma = 0.0;
    gss_x0 = 0.0;
#elif defined __GEN_RAND__
    srand(idum);
#elif defined __GEN_DSFMT__
    dsfmt_init_gen_rand ( &mydsfmt, idum );
#endif
}


void MyRand::init(long idum)
{
#ifdef __GEN_RAND55__
    int i,j;
    long tmp,aux;

    for(i=0;i<54;++i) {
      next[i] = i+1;
    }
    next[54]= 0;
    tmp = Seed + idum;
    ring[54] = tmp;
    aux = 1;
    j = 20;
    for(i=0;i<54;++i,j+=21){
	j %= 55;
	ring[j] = aux;
	aux += tmp;
	tmp = ring[j];
    }

    n1 = 0;
    n2 = 31;

    for(j=0;j<55*4;++j) {
      Ranf();
    }

    gss_flag = 0;
    gss_value = 0.0;
    gss_sigma = 0.0;
    gss_x0 = 0.0;
#elif defined __GEN_RAND__   
    srand(idum);
#elif defined __GEN_DSFMT__
    dsfmt_init_gen_rand ( &mydsfmt, idum );
#endif
}




double MyRand::generate_gauss(){
  double x;

  if(gen == 0){
    t = MY_2pi*Ranf();
    r = sqrt( -log((1. - Ranf())) );
    x = r*cos(t);
    gen = 1;
  }
  else{
    x = r*sin(t);
    gen = 0;
  }
  
  return(x);
}






double* MyRand::generate_gauss(double* v, int n){
  register double r,t;

  Ranf(v,n);

  for( int i = 0; i < n; i+=2)
    {
      t = MY_2pi*v[i];
      r = sqrt( -log((1. - v[i+1])) );
      v[i]   = r*cos(t);
      v[i+1] = r*sin(t);
  }
  
  return v;
}

