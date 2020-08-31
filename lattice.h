/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: lattice.h
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

#include <iostream>

#define dim 4
#define x0 1
#define x1 2
#define x2 3
#define x3 4

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MAX4(x, y, z, t) (MAX(MAX((x), (y)), MAX((z), (t))))

#define GETINDEX(x, y, z, t) ((x)*(sizeYZT) + (y)*(sizeZT) + (z)*(sizeT) + (t))
#define MOD(i, m) (((m) + (i)) % (m))

#define ERRORVALUE -999999999

#define SGN(a) ((a) >= 0 ? +1 : -1)

//#define ABC

const int lat_offset = dim+1;

class latt;

class SU3_fld;
class Gluon_fld;
class CVector_fld;

class ptSU3_fld;
class ptGluon_fld;
class ptCVector_fld;


// An important notice on the structure of lattice class: the main structure
// is the array of array LL. This is a vector of size Size, and each element
// is an array of size 9. The ordering of the array is the "phisical" ordering
// of space-time. The 5th(#4) element of the array refers to where a given degree
// of freedom attached to a site is stored within a field in computer memory.
// The elements from 0 to 3 are the location of the backward nearest neighbours,
// while the elements from 5 to 8 are the location of forward nearest neighbours.
// (See below)
//
//    #0    #1    #2    #3    #4    #5  #6   #7   #8
//  .     .                 .     .                   .
//  .......                 -------                   .
//  |     |                 | n-1 |                   |
//  |-----------------------|-----|--------------------
//  | -x0 | -x1 | -x2 | -x3 | n   | x0 | x1 | x2 | x3 |
//  |-----------------------|-----|--------------------
//  |                       | n+1 |                   |
//  .                       |-----|                   .
//  .                       .     .                   .
//


class latt{

 public:
  int *Sz;

#if (dim == 4)

  int sizeYZT; 
  int sizeZT;
  int sizeT;
  int Size;   

  int TSize;
#ifdef TWISTED_BC
    int *SizeBoundary,*SizeBoundary_f,*SizeBoundary_b;
    int *copyoffset_f,*copyoffset_b;
    int **justtwist;
    int *sizejusttwist;
#ifdef _TWIST_X
    int **copyandtwistF_x;
    int sizecopyandtwistF_x;
    int **copyandtwistB_x;
    int sizecopyandtwistB_x;
#endif
#ifdef _TWIST_Y
    int **copyandtwistF_y;
    int sizecopyandtwistF_y;
    int **copyandtwistB_y;
    int sizecopyandtwistB_y;
#endif
#ifdef _TWIST_Z
    int **copyandtwistF_z;
    int sizecopyandtwistF_z;
    int **copyandtwistB_z;
    int sizecopyandtwistB_z;
#endif
#ifdef _TWIST_T
    int **copyandtwistF_t;
    int sizecopyandtwistF_t;
    int **copyandtwistB_t;
    int sizecopyandtwistB_t;
#endif
#endif
    
  double **p;
  double **pbar;
  double **p2bar;
  double **p2hat;

  double ***tw_pperp;
  double ****tw_pbar;
  double ***tw_sumpbar2;
  double ***tw_halfsumphat2;
    
#endif


  int **L;  
  latt(int *size) {
    int c = 0;

#if (dim == 4)
    Sz = new int[4];
    Sz[0] = size[0];
    Sz[1] = size[1];
    Sz[2] = size[2];
    Sz[3] = size[3];
    sizeT = Sz[3];
    sizeZT= Sz[2]*Sz[3];
    sizeYZT  = Sz[1]*Sz[2]*Sz[3];
    Size    = Sz[0]*Sz[1]*Sz[2]*Sz[3];
    
    
    
    TSize = Size;
#ifdef TWISTED_BC
      // allocation of forward and backward copies
      SizeBoundary = new int[4];
      SizeBoundary_f = new int[4];
      SizeBoundary_b = new int[4];
#ifdef _TWIST_X
      SizeBoundary_f[0] = Sz[1]*Sz[2]*Sz[3];
      SizeBoundary_b[0] = (Sz[1]+1)*(Sz[2]+1)*(Sz[3]+1);
      SizeBoundary[0] = SizeBoundary_f[0] + SizeBoundary_b[0];
      TSize += SizeBoundary[0];
#else
      SizeBoundary_f[0] = 0;
      SizeBoundary_b[0] = 0;
      SizeBoundary[0] = 0;
#endif
#ifdef _TWIST_Y
      SizeBoundary_f[1] = Sz[0]*Sz[2]*Sz[3];
      SizeBoundary_b[1] = (Sz[0]+1)*(Sz[2]+1)*(Sz[3]+1);
      SizeBoundary[1] = SizeBoundary_f[1] + SizeBoundary_b[1];
      TSize += SizeBoundary[1];
#else
      SizeBoundary_f[1] = 0;
      SizeBoundary_b[1] = 0;
      SizeBoundary[1] = 0;
#endif
#ifdef _TWIST_Z
      SizeBoundary_f[2] = Sz[0]*Sz[1]*Sz[3];
      SizeBoundary_b[2] = (Sz[0]+1)*(Sz[1]+1)*(Sz[3]+1);
      SizeBoundary[2] = SizeBoundary_f[2] + SizeBoundary_b[2];
      TSize += SizeBoundary[2];
#else
      SizeBoundary_f[2] = 0;
      SizeBoundary_b[2] = 0;
      SizeBoundary[2] = 0;
#endif
#ifdef _TWIST_T
      SizeBoundary_f[3] = Sz[0]*Sz[1]*Sz[2];
      SizeBoundary_b[3] = (Sz[0]+1)*(Sz[1]+1)*(Sz[2]+1);
      SizeBoundary[3] = SizeBoundary_f[3] + SizeBoundary_b[3];
      TSize += SizeBoundary[3];
#else
      SizeBoundary_f[3] = 0;
      SizeBoundary_b[3] = 0;
      SizeBoundary[3] = 0;
#endif
#endif
      
      
    L = new int*[TSize];
    
    for(int ii=0; ii<TSize; ii++){
      L[ii] = new int[1+2*dim];      
    }    
  
    for (int i = 0; i < TSize; i++){
      L[i][dim] = c++;
      //printf("%d  ->%d \n", i,L[i][4]);
    }  


    int i;
    for (int x = 0; x < Sz[0]; x++)
	for (int y = 0; y < Sz[1]; y++)
	    for (int z = 0; z < Sz[2]; z++)
		for (int t = 0; t < Sz[3]; t++) {
		    i = GETINDEX(x,y,z,t);

		    L[i][0] = GETINDEX(MOD(x-1,Sz[0]), y, z, t);
		    L[i][1] = GETINDEX(x, MOD(y-1,Sz[1]), z, t);
		    L[i][2] = GETINDEX(x, y, MOD(z-1,Sz[2]), t);
		    L[i][3] = GETINDEX(x, y, z, MOD(t-1,Sz[3]));

		    L[i][5] = GETINDEX((x+1)%Sz[0], y, z, t);
		    L[i][6] = GETINDEX(x, (y+1)%Sz[1], z, t);
		    L[i][7] = GETINDEX(x, y, (z+1)%Sz[2], t);
		    L[i][8] = GETINDEX(x, y, z, (t+1)%Sz[3]);
		    
		}
      
      
#ifdef TWISTED_BC
      // from position Size on there are all the copies of the boundary.
      // in the copyandtwist array there is the correspondence given by
      // periodic boundary conditions, both for forward and backward sites.
      // 0 is the original
      // 1 is the copy
      // the justtwist[mu][i] array collects additional sites that
      // undergo a mu twist.
      copyoffset_f = new int[4];
      copyoffset_b = new int[4];
      justtwist = new int*[4];
      sizejusttwist = new int[4];
      
      sizejusttwist[0] = 0;
      if(SizeBoundary[1]>0) sizejusttwist[0] += Sz[2]*Sz[3];
      if(SizeBoundary[2]>0) sizejusttwist[0] += Sz[1]*Sz[3];
      if(SizeBoundary[3]>0) sizejusttwist[0] += Sz[1]*Sz[2];
      
      sizejusttwist[1] = 0;
      if(SizeBoundary[0]>0) sizejusttwist[1] += Sz[2]*Sz[3];
      if(SizeBoundary[2]>0) sizejusttwist[1] += Sz[0]*Sz[3];
      if(SizeBoundary[3]>0) sizejusttwist[1] += Sz[0]*Sz[2];
      
      sizejusttwist[2] = 0;
      if(SizeBoundary[0]>0) sizejusttwist[2] += Sz[1]*Sz[3];
      if(SizeBoundary[1]>0) sizejusttwist[2] += Sz[0]*Sz[3];
      if(SizeBoundary[3]>0) sizejusttwist[2] += Sz[0]*Sz[1];
      
      sizejusttwist[3] = 0;
      if(SizeBoundary[0]>0) sizejusttwist[3] += Sz[1]*Sz[2];
      if(SizeBoundary[1]>0) sizejusttwist[3] += Sz[0]*Sz[2];
      if(SizeBoundary[2]>0) sizejusttwist[3] += Sz[0]*Sz[1];
      
#ifdef _TWIST_X
      copyoffset_f[0] = Size;
      copyoffset_b[0] = copyoffset_f[0] + SizeBoundary_f[0];
      
      sizecopyandtwistF_x = Sz[1]*Sz[2]*Sz[3];
      copyandtwistF_x = new int*[sizecopyandtwistF_x];
      for (int i=0; i<sizecopyandtwistF_x; i++) {
          copyandtwistF_x[i] = new int[2];
      }
      sizecopyandtwistB_x = Sz[1]*Sz[2]*Sz[3] + Sz[1]*Sz[2] + Sz[1]*Sz[3] + Sz[3]*Sz[2];
      copyandtwistB_x = new int*[sizecopyandtwistB_x];
      for (int i=0; i<sizecopyandtwistB_x; i++) {
          copyandtwistB_x[i] = new int[2];
      }
      
      if(sizejusttwist[0]>0) justtwist[0] = new int[sizejusttwist[0]];
#else
      copyoffset_f[0] = 0;
      copyoffset_b[0] = 0;
#endif
#ifdef _TWIST_Y
      copyoffset_f[1] = Size + SizeBoundary[0];
      copyoffset_b[1] = copyoffset_f[1] + SizeBoundary_f[1];
      
      sizecopyandtwistF_y = Sz[0]*Sz[2]*Sz[3];
      copyandtwistF_y = new int*[sizecopyandtwistF_y];
      for (int i=0; i<sizecopyandtwistF_y; i++) {
          copyandtwistF_y[i] = new int[2];
      }
      sizecopyandtwistB_y = Sz[0]*Sz[2]*Sz[3] + Sz[0]*Sz[2] + Sz[0]*Sz[3] + Sz[3]*Sz[2];
      copyandtwistB_y = new int*[sizecopyandtwistB_y];
      for (int i=0; i<sizecopyandtwistB_y; i++) {
          copyandtwistB_y[i] = new int[2];
      }
      
      if(sizejusttwist[1]>0) justtwist[1] = new int[sizejusttwist[1]];
#else
      copyoffset_f[1] = 0;
      copyoffset_b[1] = 0;
#endif
#ifdef _TWIST_Z
      copyoffset_f[2] = Size + SizeBoundary[0] + SizeBoundary[1];
      copyoffset_b[2] = copyoffset_f[2] + SizeBoundary_f[2];
      
      sizecopyandtwistF_z = Sz[1]*Sz[0]*Sz[3];
      copyandtwistF_z = new int*[sizecopyandtwistF_z];
      for (int i=0; i<sizecopyandtwistF_z; i++) {
          copyandtwistF_z[i] = new int[2];
      }
      sizecopyandtwistB_z = Sz[1]*Sz[0]*Sz[3] + Sz[1]*Sz[0] + Sz[1]*Sz[3] + Sz[0]*Sz[3];
      copyandtwistB_z = new int*[sizecopyandtwistB_z];
      for (int i=0; i<sizecopyandtwistB_z; i++) {
          copyandtwistB_z[i] = new int[2];
      }
      
      if(sizejusttwist[2]>0) justtwist[2] = new int[sizejusttwist[2]];
#else
      copyoffset_f[2] = 0;
      copyoffset_b[2] = 0;
#endif
#ifdef _TWIST_T
      copyoffset_f[3] = Size + SizeBoundary[0] + SizeBoundary[1] + SizeBoundary[2];
      copyoffset_b[3] = copyoffset_f[3] + SizeBoundary_f[3];
      
      sizecopyandtwistF_t = Sz[1]*Sz[2]*Sz[0];
      copyandtwistF_t = new int*[sizecopyandtwistF_t];
      for (int i=0; i<sizecopyandtwistF_t; i++) {
          copyandtwistF_t[i] = new int[2];
      }
      sizecopyandtwistB_t = Sz[1]*Sz[2]*Sz[0] + Sz[1]*Sz[2] + Sz[1]*Sz[0] + Sz[0]*Sz[2];
      copyandtwistB_t = new int*[sizecopyandtwistB_t];
      for (int i=0; i<sizecopyandtwistB_t; i++) {
          copyandtwistB_t[i] = new int[2];
      }
      
      if(sizejusttwist[3]>0) justtwist[3] = new int[sizejusttwist[3]];
#else
      copyoffset_f[3] = 0;
      copyoffset_b[3] = 0;
#endif
      
      // Now neighbours are defined,
      // the initialization is redundant but safe.
      // Neighbours assigned to ERRORVALUE=-999 must not participate to anything.
      // If they are, it is an error.
      // There are some true neighbourhood relations that are missing,
      // but they should not play any role in the dynamics.
      
      int x[4],counterF, counterB;
      
#ifdef _TWIST_X
      // FORWARD COPIES
      x[0] = Sz[0]-1;
      for (x[1]=0; x[1]<Sz[1]; x[1]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  int  site_copy = getxcopy_f(x,copyoffset_f[0],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=0){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getxcopy_f(x,copyoffset_f[0],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]-1) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getxcopy_f(x,copyoffset_f[0],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          app = getreal(x,Sz);
                          L[site_copy][mu] = app;
                          L[site_copy][lat_offset+mu] = ERRORVALUE;
                          L[app][lat_offset+mu] = site_copy;
                      }
                  } // mu loop
              }
      // BACKWARD COPIES
      x[0] = 0;
      for (x[1]=0; x[1]<=Sz[1]; x[1]++)
          for (x[2]=0; x[2]<=Sz[2]; x[2]++)
              for (x[3]=0; x[3]<=Sz[3]; x[3]++){
                  int  site_copy = getxcopy_b(x,copyoffset_b[0],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=0){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getxcopy_b(x,copyoffset_b[0],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getxcopy_b(x,copyoffset_b[0],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          if(x[1]!=Sz[1] && x[2]!=Sz[2] && x[3]!=Sz[3]){
                              app = getreal(x,Sz);
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                              x[mu] = Sz[mu]-1;
                              app = getreal(x,Sz);
                              x[mu] = 0;
                              L[site_copy][mu] = app; // useful for stochastic gauge fixing
                              
                          } else {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                              L[site_copy][mu] = ERRORVALUE;
                          }
                      }
                  } // mu loop
              }
#endif
#ifdef _TWIST_Y
      // FORWARD COPIES
      x[1] = Sz[1]-1;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  int  site_copy = getycopy_f(x,copyoffset_f[1],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=1){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getycopy_f(x,copyoffset_f[1],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]-1) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getycopy_f(x,copyoffset_f[1],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          app = getreal(x,Sz);
                          L[site_copy][mu] = app;
                          L[site_copy][lat_offset+mu] = ERRORVALUE;
                          L[app][lat_offset+mu] = site_copy;
                      }
                  } // mu loop
              }
      // BACKWARD COPIES
      x[1] = 0;
      for (x[0]=0; x[0]<=Sz[0]; x[0]++)
          for (x[2]=0; x[2]<=Sz[2]; x[2]++)
              for (x[3]=0; x[3]<=Sz[3]; x[3]++){
                  int  site_copy = getycopy_b(x,copyoffset_b[1],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=1){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getycopy_b(x,copyoffset_b[1],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getycopy_b(x,copyoffset_b[1],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          if(x[0]!=Sz[0] && x[2]!=Sz[2] && x[3]!=Sz[3]){
                              app = getreal(x,Sz);
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                              x[mu] = Sz[mu]-1;
                              app = getreal(x,Sz);
                              x[mu] = 0;
                              L[site_copy][mu] = app; // useful for stochastic gauge fixing
                          } else {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                              L[site_copy][mu] = ERRORVALUE;
                          }
                      }
                  } // mu loop
              }
#endif
#ifdef _TWIST_Z
      // FORWARD COPIES
      x[2] = Sz[2]-1;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  int  site_copy = getzcopy_f(x,copyoffset_f[2],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=2){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getzcopy_f(x,copyoffset_f[2],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]-1) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getzcopy_f(x,copyoffset_f[2],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          app = getreal(x,Sz);
                          L[site_copy][mu] = app;
                          L[site_copy][lat_offset+mu] = ERRORVALUE;
                          L[app][lat_offset+mu] = site_copy;
                      }
                  } // mu loop
              }
      // BACKWARD COPIES
      x[2] = 0;
      for (x[0]=0; x[0]<=Sz[0]; x[0]++)
          for (x[1]=0; x[1]<=Sz[1]; x[1]++)
              for (x[3]=0; x[3]<=Sz[3]; x[3]++){
                  int  site_copy = getzcopy_b(x,copyoffset_b[2],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=2){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = getzcopy_b(x,copyoffset_b[2],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = getzcopy_b(x,copyoffset_b[2],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          if(x[0]!=Sz[0] && x[1]!=Sz[1] && x[3]!=Sz[3]){
                              app = getreal(x,Sz);
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                              x[mu] = Sz[mu]-1;
                              app = getreal(x,Sz);
                              x[mu] = 0;
                              L[site_copy][mu] = app; // useful for stochastic gauge fixing
                          } else {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                              L[site_copy][mu] = ERRORVALUE;
                          }
                      }
                  } // mu loop
              }
#endif
#ifdef _TWIST_T
      // FORWARD COPIES
      x[3] = Sz[3]-1;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++)
              for (x[2]=0; x[2]<Sz[2]; x[2]++){
                  int  site_copy = gettcopy_f(x,copyoffset_f[3],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=3){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = gettcopy_f(x,copyoffset_f[3],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]-1) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = gettcopy_f(x,copyoffset_f[3],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          app = getreal(x,Sz);
                          L[site_copy][mu] = app;
                          L[site_copy][lat_offset+mu] = ERRORVALUE;
                          L[app][lat_offset+mu] = site_copy;
                      }
                  } // mu loop
              }
      // BACKWARD COPIES
      x[3] = 0;
      for (x[0]=0; x[0]<=Sz[0]; x[0]++)
          for (x[1]=0; x[1]<=Sz[1]; x[1]++)
              for (x[2]=0; x[2]<=Sz[2]; x[2]++){
                  int  site_copy = gettcopy_b(x,copyoffset_b[3],Sz);
                  int app;
                  
                  for (int mu=0; mu<dim; mu++) {
                      if(mu!=3){
                          // --- //
                          if (x[mu]==0) {
                              L[site_copy][mu] = ERRORVALUE;
                          } else {
                              x[mu]--;
                              app = gettcopy_b(x,copyoffset_b[3],Sz);
                              x[mu]++;
                              L[site_copy][mu] = app;
                              L[app][lat_offset+mu] = site_copy;
                          }
                          if (x[mu]==Sz[mu]) {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                          } else {
                              x[mu]++;
                              app = gettcopy_b(x,copyoffset_b[3],Sz);
                              x[mu]--;
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                          }
                          // --- //
                      } else {
                          if(x[0]!=Sz[0] && x[1]!=Sz[1] && x[2]!=Sz[2]){
                              app = getreal(x,Sz);
                              L[site_copy][lat_offset+mu] = app;
                              L[app][mu] = site_copy;
                              x[mu] = Sz[mu]-1;
                              app = getreal(x,Sz);
                              x[mu] = 0;
                              L[site_copy][mu] = app; // useful for stochastic gauge fixing
                          } else {
                              L[site_copy][lat_offset+mu] = ERRORVALUE;
                              L[site_copy][mu] = ERRORVALUE;
                          }
                      }
                  } // mu loop
              }
#endif
      
      
      
      // SET UP CORRESPONDENCE BETWEEN ORIGINALS AND COPIES
      int additionalmu, counterarray[4]={0};
#ifdef _TWIST_X
      counterF = 0;
      counterB = 0;
      for (x[1]=0; x[1]<Sz[1]; x[1]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  
                  // FORWARD
                  x[0] = Sz[0]-1;
                  int site_copy = getxcopy_f(x,copyoffset_f[0],Sz);
                  x[0] = 0;
                  int app = getreal(x,Sz);
                  
                  copyandtwistF_x[counterF][0] = app;
                  copyandtwistF_x[counterF++][1] = site_copy;
                  
                  // BACKWARD
                  site_copy = getxcopy_b(x,copyoffset_b[0],Sz);
                  x[0] = Sz[0]-1;
                  app = getreal(x,Sz);
                  
                  copyandtwistB_x[counterB][0] = app;
                  copyandtwistB_x[counterB++][1] = site_copy;
              }
      
      // additional backwards
      additionalmu = 1;
      for (x[2]=0; x[2]<Sz[2]; x[2]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[0] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getxcopy_b(x,copyoffset_b[0],Sz);
              x[0] = Sz[0]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_x[counterB][0] = app;
              copyandtwistB_x[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 2;
      for (x[1]=0; x[1]<Sz[1]; x[1]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[0] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getxcopy_b(x,copyoffset_b[0],Sz);
              x[0] = Sz[0]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_x[counterB][0] = app;
              copyandtwistB_x[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 3;
      for (x[2]=0; x[2]<Sz[2]; x[2]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++){
              
              x[0] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getxcopy_b(x,copyoffset_b[0],Sz);
              x[0] = Sz[0]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_x[counterB][0] = app;
              copyandtwistB_x[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
#endif
#ifdef _TWIST_Y
      counterF = 0;
      counterB = 0;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  
                  // FORWARD
                  x[1] = Sz[1]-1;
                  int site_copy = getycopy_f(x,copyoffset_f[1],Sz);
                  x[1] = 0;
                  int app = getreal(x,Sz);
                  
                  copyandtwistF_y[counterF][0] = app;
                  copyandtwistF_y[counterF++][1] = site_copy;
                  
                  // BACKWARD
                  site_copy = getycopy_b(x,copyoffset_b[1],Sz);
                  x[1] = Sz[1]-1;
                  app = getreal(x,Sz);
                  
                  copyandtwistB_y[counterB][0] = app;
                  copyandtwistB_y[counterB++][1] = site_copy;
              }
      
      // additional backwards
      additionalmu = 0;
      for (x[2]=0; x[2]<Sz[2]; x[2]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[1] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getycopy_b(x,copyoffset_b[1],Sz);
              x[1] = Sz[1]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_y[counterB][0] = app;
              copyandtwistB_y[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 2;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[1] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getycopy_b(x,copyoffset_b[1],Sz);
              x[1] = Sz[1]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_y[counterB][0] = app;
              copyandtwistB_y[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 3;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++){
              
              x[1] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getycopy_b(x,copyoffset_b[1],Sz);
              x[1] = Sz[1]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_y[counterB][0] = app;
              copyandtwistB_y[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
#endif
#ifdef _TWIST_Z
      counterF = 0;
      counterB = 0;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++)
              for (x[3]=0; x[3]<Sz[3]; x[3]++){
                  
                  // FORWARD
                  x[2] = Sz[2]-1;
                  int site_copy = getzcopy_f(x,copyoffset_f[2],Sz);
                  x[2] = 0;
                  int app = getreal(x,Sz);
                  
                  copyandtwistF_z[counterF][0] = app;
                  copyandtwistF_z[counterF++][1] = site_copy;
                  
                  // BACKWARD
                  site_copy = getzcopy_b(x,copyoffset_b[2],Sz);
                  x[2] = Sz[2]-1;
                  app = getreal(x,Sz);
                  
                  copyandtwistB_z[counterB][0] = app;
                  copyandtwistB_z[counterB++][1] = site_copy;
              }
      
      // additional backwards
      additionalmu = 0;
      for (x[1]=0; x[1]<Sz[1]; x[1]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[2] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getzcopy_b(x,copyoffset_b[2],Sz);
              x[2] = Sz[2]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_z[counterB][0] = app;
              copyandtwistB_z[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 1;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[3]=0; x[3]<Sz[3]; x[3]++){
              
              x[2] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getzcopy_b(x,copyoffset_b[2],Sz);
              x[2] = Sz[2]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_z[counterB][0] = app;
              copyandtwistB_z[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 3;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++){
              
              x[2] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = getzcopy_b(x,copyoffset_b[2],Sz);
              x[2] = Sz[2]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_z[counterB][0] = app;
              copyandtwistB_z[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
#endif
#ifdef _TWIST_T
      counterF = 0;
      counterB = 0;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++)
              for (x[2]=0; x[2]<Sz[2]; x[2]++){
                  
                  // FORWARD
                  x[3] = Sz[3]-1;
                  int site_copy = gettcopy_f(x,copyoffset_f[3],Sz);
                  x[3] = 0;
                  int app = getreal(x,Sz);
                  
                  copyandtwistF_t[counterF][0] = app;
                  copyandtwistF_t[counterF++][1] = site_copy;
                  
                  // BACKWARD
                  site_copy = gettcopy_b(x,copyoffset_b[3],Sz);
                  x[3] = Sz[3]-1;
                  app = getreal(x,Sz);
                  
                  copyandtwistB_t[counterB][0] = app;
                  copyandtwistB_t[counterB++][1] = site_copy;
              }
      
      // additional backwards
      additionalmu = 0;
      for (x[1]=0; x[1]<Sz[1]; x[1]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++){
              
              x[3] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = gettcopy_b(x,copyoffset_b[3],Sz);
              x[3] = Sz[3]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_t[counterB][0] = app;
              copyandtwistB_t[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 1;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[2]=0; x[2]<Sz[2]; x[2]++){
              
              x[3] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = gettcopy_b(x,copyoffset_b[3],Sz);
              x[3] = Sz[3]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_t[counterB][0] = app;
              copyandtwistB_t[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
      
      additionalmu = 2;
      for (x[0]=0; x[0]<Sz[0]; x[0]++)
          for (x[1]=0; x[1]<Sz[1]; x[1]++){
              
              x[3] = 0;
              x[additionalmu] = Sz[additionalmu];
              int site_copy = gettcopy_b(x,copyoffset_b[3],Sz);
              x[3] = Sz[3]-1;
              x[additionalmu] = 0;
              int app = getreal(x,Sz);
              
              copyandtwistB_t[counterB][0] = app;
              copyandtwistB_t[counterB++][1] = site_copy;
              
              if (SizeBoundary[additionalmu]>0) {
                  justtwist[additionalmu][counterarray[additionalmu]++] = site_copy;
              }
              
          }
#endif
#endif // #ifdef TWISTED_BC
  }

  

  void p_init(){
    p     = new double* [4];
    pbar  = new double* [4];
    p2bar = new double* [4];
    p2hat = new double* [4];

    for (int i = 0; i < dim; i++) {
	p[i]     = new double [Sz[i]];
	pbar[i]  = new double [Sz[i]];
	p2bar[i] = new double [Sz[i]];
	p2hat[i] = new double [Sz[i]];
    }
      
      
#ifdef TWISTED_BC
      // tw_pperp[n1][n2][mu]
      // tw_pbar[x][n1][n2][mu]
      // tw_sumpbar2[x][n1][n2]
      // tw_halfsumphat2[x][n1][n2]
      
      tw_pperp = new double** [NC];
      tw_pbar = new double*** [Size];
      tw_sumpbar2 = new double** [Size];
      tw_halfsumphat2 = new double** [Size];
      
      for (int n1=0; n1<NC; n1++) {
          tw_pperp[n1] = new double* [NC];
          for (int n2=0; n2<NC; n2++) {
              tw_pperp[n1][n2] = new double [dim];
          }
      }
      for (int i=0; i<Size; i++) {
          tw_pbar[i] = new double** [NC];
          tw_sumpbar2[i] = new double* [NC];
          tw_halfsumphat2[i] = new double* [NC];
          for (int n1=0; n1<NC; n1++) {
              tw_pbar[i][n1] = new double* [NC];
              tw_sumpbar2[i][n1] = new double [NC];
              tw_halfsumphat2[i][n1] = new double [NC];
              for (int n2=0; n2<NC; n2++) {
                  tw_pbar[i][n1][n2] = new double [dim];
              }
          }
      }
      
      for (int n1=0; n1<NC; n1++) {
          for (int n2=0; n2<NC; n2++) {
              tw_pperp[n1][n2][0] = 0;
              tw_pperp[n1][n2][1] = 2.*M_PI*(double)n1/(Sz[1]*NC);
              tw_pperp[n1][n2][2] = 2.*M_PI*(double)n2/(Sz[2]*NC);
              tw_pperp[n1][n2][3] = 0;
          }
      }
      
      for (int x = 0; x < Sz[0]; x++)
          for (int y = 0; y < Sz[1]; y++)
              for (int z = 0; z < Sz[2]; z++)
                  for (int t = 0; t < Sz[3]; t++)
                      for (int n1 = 0; n1 < NC; n1++)
                          for (int n2 = 0; n2 < NC; n2++){
                              tw_pbar[GETINDEX(x, y, z, t)][n1][n2][0] = sin(2. * M_PI * (x + .5) / Sz[0]);
                              tw_pbar[GETINDEX(x, y, z, t)][n1][n2][1] = sin(2. * M_PI * (y + (double)n1/NC) / Sz[1]);
                              tw_pbar[GETINDEX(x, y, z, t)][n1][n2][2] = sin(2. * M_PI * (z + (double)n2/NC) / Sz[2]);
                              tw_pbar[GETINDEX(x, y, z, t)][n1][n2][3] = sin(2. * M_PI * t / Sz[3]);
                              
                              tw_halfsumphat2[GETINDEX(x, y, z, t)][n1][n2] =
                                    2*sin(M_PI*(x+.5)/Sz[0])*sin(M_PI*(x+.5)/Sz[0])+
                                    2*sin(M_PI*(y+(double)n1/NC)/Sz[1])*sin(M_PI*(y+(double)n1/NC)/Sz[1])+
                                    2*sin(M_PI*(z+(double)n2/NC)/Sz[2])*sin(M_PI*(z+(double)n2/NC)/Sz[2])+
                                    2*sin(M_PI*t/Sz[3])*sin(M_PI*t/Sz[3]);
                              
                              tw_sumpbar2[GETINDEX(x, y, z, t)][n1][n2] = 0;
                              for (int mu = 0; mu < dim; mu++) {
                                  tw_sumpbar2[GETINDEX(x, y, z, t)][n1][n2] +=
                                        tw_pbar[GETINDEX(x, y, z, t)][n1][n2][mu]*tw_pbar[GETINDEX(x, y, z, t)][n1][n2][mu];
                              }
                          }
#endif

#ifdef ABC
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * (i+.5) / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
#elif defined PBC
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * i / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
#endif
    for (int i = 0; i < Sz[1]; i++) {
	p[1][i] = 2 * M_PI * i / Sz[1];
	pbar[1][i] = sin(p[1][i]);
	p2bar[1][i] = pbar[1][i]*pbar[1][i];
	p2hat[1][i] = sin(p[1][i] * 0.5); p2hat[1][i] *= 4*p2hat[1][i];
    }
    for (int i = 0; i < Sz[2]; i++) {
	p[2][i] = 2 * M_PI * i / Sz[2];
	pbar[2][i] = sin(p[2][i]);
	p2bar[2][i] = pbar[2][i]*pbar[2][i];
	p2hat[2][i] = sin(p[2][i] * 0.5); p2hat[2][i] *= 4*p2hat[2][i];
    }
    for (int i = 0; i < Sz[3]; i++) {
	p[3][i] = 2 * M_PI * i / Sz[3];
	pbar[3][i] = sin(p[3][i]);
	p2bar[3][i] = pbar[3][i]*pbar[3][i];
	p2hat[3][i] = sin(p[3][i] * 0.5); p2hat[3][i] *= 4*p2hat[3][i];
    }

#endif
  }


  void p_pbc(){
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * i / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
  }

  void p_abc(){
    for (int i = 0; i < Sz[0]; i++) {
      p[0][i] = 2 * M_PI * (i+.5) / Sz[0];
      pbar[0][i] = sin(p[0][i]);
      p2bar[0][i] = pbar[0][i]*pbar[0][i];
      p2hat[0][i] = sin(p[0][i] * 0.5); p2hat[0][i] *= 4*p2hat[0][i];
    }
  }


  int get(int *n) {
      return (GETINDEX(n[0], n[1], n[2], n[3])); 
  }

  void get(int i, int *coord) {
    int modYZT, modZT;
    
    coord[0] = i / sizeYZT;
    coord[1] = (modYZT = i % sizeYZT) / sizeZT;
    coord[2] = (modZT = modYZT % sizeZT) / sizeT;
    coord[3] = modZT % sizeT;
  }
  
  const int& get(const int n) { return L[n][dim]; }
  
  const int& get(int n, const int step, const int dir) {
    if (step < 0)
      for (int i = 0; i < -step; ++i) n = L[n][dir];
    else
      for (int i = 0; i <  step; ++i) n = L[n][lat_offset+dir];
    return L[n][dim];
  }
  
  const int& get(int n, 
          const int step1, const int dir1,
          const int step2,const int dir2) {
    n=get(n,step1,dir1);
    return get(n,step2,dir2);
  }
  
  const int near(int n, int sign, int dir){
    if (sign < 0)
      return L[n][dir]+1;
    return L[n][lat_offset+dir]+1;
  }

  const int near2(int n, int sign, int dir){
    if ((dir == 0 && sign == +1 && (n / sizeYZT) == Sz[0]-1) || (dir == 0 && sign == -1 && (n / sizeYZT) == 0)) {
      if (sign < 0)
        return -(L[n][dir]+1);
      else
        return -(L[n][lat_offset+dir]+1);
    }
    else {
      if (sign < 0)
        return +(L[n][dir]+1);
      else
        return +(L[n][lat_offset+dir]+1);
    }
  }
    
#ifdef TWISTED_BC
    int getreal(int *x,int *Sz){
        return x[3] + x[2]*Sz[3] + x[1]*Sz[3]*Sz[2] + x[0]*Sz[3]*Sz[2]*Sz[1];
    }
#ifdef _TWIST_X
    int getxcopy_f(int *x,int offset,int *Sz){
        return x[3] + x[2]*Sz[3] + x[1]*Sz[3]*Sz[2] + offset;
    }
    int getxcopy_b(int *x,int offset,int *Sz){
        return x[3] + x[2]*(Sz[3]+1) + x[1]*(Sz[3]+1)*(Sz[2]+1) + offset;
    }
#endif
#ifdef _TWIST_Y
    int getycopy_f(int *x,int offset,int *Sz){
        return x[3] + x[2]*Sz[3] + x[0]*Sz[3]*Sz[2] + offset;
    }
    int getycopy_b(int *x,int offset,int *Sz){
        return x[3] + x[2]*(Sz[3]+1) + x[0]*(Sz[3]+1)*(Sz[2]+1) + offset;
    }
#endif
#ifdef _TWIST_Z
    int getzcopy_f(int *x,int offset,int *Sz){
        return x[3] + x[1]*Sz[3] + x[0]*Sz[2]*Sz[3] + offset;
    }
    int getzcopy_b(int *x,int offset,int *Sz){
        return x[3] + x[1]*(Sz[3]+1) + x[0]*(Sz[2]+1)*(Sz[3]+1) + offset;
    }
#endif
#ifdef _TWIST_T
    int gettcopy_f(int *x,int offset,int *Sz){
        return x[2] + x[1]*Sz[2] + x[0]*Sz[1]*Sz[2] + offset;
    }
    int gettcopy_b(int *x,int offset,int *Sz){
        return x[2] + x[1]*(Sz[2]+1) + x[0]*(Sz[1]+1)*(Sz[2]+1) + offset;
    }
#endif
#endif
};    // end class latt

