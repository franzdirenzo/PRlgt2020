/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: QCDenvNODEpt_test.h
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

#ifndef _QCD_ENV_NODE_PT_
#define _QCD_ENV_NODE_PT_

#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

#include "QCDenvNODE.h"

#define ADD 1
#define NOADD 0
#define DAG 1
#define NODAG 0


#define SGN(a) ((a) >= 0 ? +1 : -1)

#define SWAP(a, b, s) {(s) = (a); (a) = (b); (b) = (s);}


namespace io {
  
  template<class T, bool contains_plaquette=false>
    int load(const char *filename, T* data) {
    std::ifstream is (filename, std::ifstream::binary);
    int length;
    if (is) {
      // get length of file:
      is.seekg (0, is.end);
      length = is.tellg();
      if(contains_plaquette)
        length-=(allocORD+1)*sizeof(Cplx);// we are not interested in the plaquette values
      is.seekg (0, is.beg);
      is.read (reinterpret_cast<char*>(data),length);
    }
    if (is)
      std::cout << "data read successfully." << std::endl;
    else
      std::cout << "error: only " << (length = is.gcount()) << " could be read";
    is.close();
    if(contains_plaquette)
      length+=(allocORD+1)*sizeof(Cplx);// we are wow interested in having read the plaquette for cross-check sizes
    return length;
  }


  
  template<class T>
    int load(const char *filename, T* data, Cplx* plaquette) {
    std::ifstream is (filename, std::ifstream::binary);
    int length;
    if (is) {
      // get length of file:
      is.seekg (0, is.end);
      length = is.tellg();
      length-=(allocORD+1)*sizeof(Cplx);// we are not interested in the plaquette values
      is.seekg (0, is.beg);
      is.read (reinterpret_cast<char*>(data),length);
      is.read (reinterpret_cast<char*>(plaquette),(allocORD+1)*sizeof(Cplx));
    }
    if (is)
      std::cout << "data read successfully." << std::endl;
    else
      std::cout << "error: only " << (length = is.gcount()) << " could be read";
    is.close();
    return length;
  }


};


class ptSU3_fld{
 private:
  latt *Z;
 public:
  ptSU3  *W;

  ptSU3_fld(latt* z) { Z=z; W = new ptSU3[Z->TSize]; }
  
  ptSU3_fld(FILE*,int); // SU3_fld(&input_file,read_mode)

  ~ptSU3_fld(){
    delete[] W;
  }

  int save(char *filename) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptSU3)*Z->Size);
    return 0;
  }


  int load(const char *filename) {
    if( io::load(filename,W) != Z->Size*sizeof(ptSU3) )
      return -1;
    return 0;

  }
  

  int load_ape(const char *in){
    
    FILE *filept;
    int *xx, *YY;
    xx = new int[dim];
    YY = new int[2];
    
    char *nome;
    nome = new char[100];

    sprintf(nome,"%s.X0",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    


    
    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		


	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X1",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    


    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+8;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X2",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    
    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+16;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    /***********************************/

    sprintf(nome,"%s.X3",in);
    if(!(filept = fopen(nome, "r"))){
      std::cout << "no configurazione " << nome << " !\n";
      exit(0);
    }
    std::cout << "File " << nome << " trovato, ora carico i dati\n";    

    for(YY[0] = 0; YY[0] < 8 ; YY[0]++){
      for(YY[1] = 0; YY[1] < 8 ; YY[1]++){
	
	for (int l = 0; l < 32; l++){ // t
	  for (int k = 0; k < 4; k++){ // z
	    for (int j = 0; j < 4; j++){ // y 
	      for (int i = 0; i < 8; i++){ // per noi e' x
		
		xx[3] = i+24;
		xx[2] = j + YY[0]*4;
		xx[1] = k + YY[1]*4;
		xx[0] = l;
		
		fread(W+Z->get(xx), sizeof(double), 2*(1+PTORD*NC*NC),filept);
		
	      }
	    }
	  }
	}
	
      }
    }

    fclose(filept);
    delete [] nome;

    return 1;
  }




  ptSU3* handle(){ return W; }

  ptSU3 get(int *);

  ptSU3 get(int n){
    return W[Z->L[n][4]];
  }
  
  ptSU3 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(ptSU3_fld *, int*);
  friend int get(ptSU3_fld *, int, int, int);
  friend int get(ptSU3_fld *, int);

};
  




inline ptSU3 ptSU3_fld::get(int *n) {
  return W[Z->get(n)];
}



inline int get(ptSU3_fld *W, int *n) {
  return (W->Z)->get(n);

}


inline int get(ptSU3_fld *W,int n,int step,int dir){
  return (W->Z)->get(n,step,dir);
}

inline int get(ptSU3_fld *W,int n){
  return (W->Z)->get(n);
}


// ------------- end class pt SU3_fld  --------------



class ptGluon_fld{
 private:

 public:
  latt   *Z;
  ptGluon  *W;
  
  ptGluon_fld(latt *z) { Z=z; W = new ptGluon[z->TSize];  }
  
  ~ptGluon_fld(){  delete [] W; }

  int save(char *filename) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptGluon)*Z->Size);
    if ( os.rdstate() != 0 )
      std::cerr << "Error writing " << filename << "\n";
    return 0;
  }


  int load(const char *filename) {
    if( io::load<ptGluon,false>(filename,W) != Z->Size*sizeof(ptGluon) )
      return -1;
    return 0;
  }
    
  // plaquette value is saves as well
  // as a check for correct configuration loading
  int save(char *filename, Cplx *placchetta) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptGluon)*Z->Size);
    os.write(reinterpret_cast<char*>(placchetta),sizeof(Cplx)*(allocORD+1));
    if ( os.rdstate() != 0 )
      std::cerr << "Error writing " << filename << "\n";
    return 0;
  }

  // plaquette value is read as well
  // as a check for correct configuration loading
  int load(const char *filename, Cplx *placchetta) {
    if( io::load<ptGluon>(filename,W,placchetta) != (Z->Size*sizeof(ptGluon)/*+(allocORD+1)*sizeof(Cplx)*/) )
      return -1;
    return 0;
  }

    int loadlowerorder(const char *filename, int loworder, Cplx *plaquette) {
      ptSU3 myapp;
        if(loworder>allocORD){
            cout<<"ERROR: trying to load more orders than the ones allocated."<<endl;
            return -1;
        }
        // ios::ate set the initial position at the end of the file, for a check on dimensions
        ifstream handle(filename, ios::binary | ios::ate);
        if (!handle.is_open()){
            cout<<"ERROR: file "<<filename<<" not found."<<endl;
            return -1;
        }
        if(handle.tellg()!=sizeof(double)*(2*(NC*NC*loworder+1)*dim*Z->Size+2*(loworder+1))){
            cout<<"ERROR: not possible to read the configuration, unexpected file dimension."<<endl;
            return -1;
        }
        handle.close();
        // open the file to read
        handle.open(filename, ios::binary);
        for (int x=0; x<Z->Size; x++){
            for (int mu=0; mu<dim; mu++){
                // read flag
                handle.read(reinterpret_cast<char*>( &(W[x].U[mu].flag.re) ), sizeof(double));
                handle.read(reinterpret_cast<char*>( &(W[x].U[mu].flag.im) ), sizeof(double));
                // read orders from file
                for (int k=0; k<loworder; k++){
                    for (int a=0; a<NC*NC; a++){
                        handle.read(reinterpret_cast<char*>( &(W[x].U[mu].ptU[k].whr[a].re) ), sizeof(double));
                        handle.read(reinterpret_cast<char*>( &(W[x].U[mu].ptU[k].whr[a].im) ), sizeof(double));
                    }
                }
                // set higher orders to zero
                for (int k=loworder+1; k<allocORD; k++){
                    for (int a=0; a<NC*NC; a++){
                        W[x].U[mu].ptU[k].whr[a].re = 0.;
                        W[x].U[mu].ptU[k].whr[a].im = 0.;
                    }
                }
		// the field has to be unitary
		myapp = log(W[x].U[mu]);
		myapp = myapp.reH();
		W[x].U[mu] = exp(myapp);
            }
        }
        // read plaquette
         for (int k=0; k<allocORD+1; k++){
             handle.read(reinterpret_cast<char*>( &(plaquette[k].re) ), sizeof(double));
             handle.read(reinterpret_cast<char*>( &(plaquette[k].im) ), sizeof(double));
         }
        return 0;
    }
    

  int load(const char *filename, Cplx *placchetta, int rank) {
    FILE *filept;
    if( (filept = fopen(filename, "r") ) == NULL ) return 1;

    long offset = Z->Sz[0]*rank*sizeof(ptGluon);
    fseek(filept, offset, SEEK_SET);
    fread(W,sizeof(ptGluon), Z->Size,filept);

    offset = -(allocORD+1)*sizeof(Cplx);
    fseek(filept, offset, SEEK_END);
    
    fread(placchetta, sizeof(Cplx), allocORD+1 , filept);

    fclose(filept);
    return 0;
  }
  


  int load_ape(const char* in){
    
    FILE *filept;
    int *xx, *YY;
    xx = new int[dim];
    YY = new int[2];
    
    char *nome;
    nome = new char[100];
    const int off=8;

    /* I files di ape sono spezzati in 4 fette. Qui ci giro sopra e li
       ricostruisco */
    for(int slice=0;slice<4;++slice) {

      sprintf(nome,"%s.X%d",in,slice);
      if(!(filept = fopen(nome, "r"))) {
        std::cout << "no configurazione " << nome << " !\n";
        exit(0);
      }
      std::cout << "File " << nome << " trovato, ora carico i dati\n";

      for(YY[0] = 0; YY[0] < 8 ; YY[0]++)
        for(YY[1] = 0; YY[1] < 8 ; YY[1]++)
          
          for (int l = 0; l < 32; l++) // t
            for (int k = 0; k < 4; k++) // z
              for (int j = 0; j < 4; j++) // y 
                for (int i = 0; i < 8; i++) { // per noi e' x
                  
                  xx[3] = l;
                  xx[2] = k + YY[1]*4;
                  xx[1] = j + YY[0]*4;
                  xx[0] = i+slice*off;
                  fread(W+Z->get(xx), sizeof(double), 2*(1+allocORD*NC*NC)*4,filept);
                }
      fclose(filept);
      
    } // slice
    delete [] nome;

    return 1;
  }
  
  
  void prout(){
    for(int i = 0; i < Z->Size; i++){
      printf("i = %d\n",i);
      W[i].prout();
      printf("\n\n");
    }
  }
  
  
  void operator=(const ptGluon_fld& A) {
    Z = A.Z;
    std::copy(A.W,A.W+Z->Size,W);
    /* for(int i = 0; i < Z->Size; i++) */
    /*   W[i] = A.W[i]; */
  }
  
  ptSU3* handle(){ return (ptSU3*)(W->U); }
  
  ptGluon& get(int* n) { return W[Z->get(n)]; }

  ptGluon& get(int n, int step, int dir) { return W[Z->get(n, step, dir)]; }

  ptSU3& get(int* n, int mu) { return W[Z->get(n)].U[mu]; }

  ptSU3& get(int n, int mu){ return W[Z->get(n)].U[mu]; }

  ptSU3& get(int n, int step, int dir, int mu){ return W[Z->get(n, step, dir)].U[mu]; }
  
  friend int get(ptGluon_fld *W, int n, int mu) { return ( 4*((W->Z)->get(n))+mu ); }

  friend int get(ptGluon_fld *W, int *n) { return W->Z->get(n); }

  friend int get(ptGluon_fld *W, int n) { return ( 4*((W->Z)->get(n)) ); }

  friend int get(ptGluon_fld *W, int n, int step , int dir, int mu) {
    return (4*((W->Z)->get(n,step,dir))+mu); 
  }

  friend int get(ptGluon_fld *W, int n, int step1 , int dir1,
		 int step2, int dir2, int mu) {
    return (4*((W->Z)->get( ((W->Z)->get(n, step1, dir1)), step2, dir2))+mu);
  }

  ptGluon& near(int n, int sign, int dir) {
    return W[Z->near(n, sign, dir)-1];
  }    


  /*
    b ->
    ^    ^
    |    |
    n == a

    n == 
    ^    ^
    |    |
    b -> a
  */
  ptSU3 staple(int*, int, int); //(x[dim], mu, nu)
  
  ptSU3 staple(int n, int mu, int nu) {
    return( W[Z->get(n, 1, mu)].U[nu]*dag(W[Z->L[n][4]].U[nu]*W[Z->get(n, 1, nu)].U[mu]) +
            dag(W[Z->get(n, -1, nu)].U[mu]*W[Z->get(n, -1, nu, 1, mu)].U[nu])*W[Z->get(n, -1, nu)].U[nu]);
  }

  ptSU3 singlestaple(int n, int mu, int nu) {
      return( W[Z->get(n, 1, mu)].U[nu]*dag(W[Z->L[n][4]].U[nu]*W[Z->get(n, 1, nu)].U[mu]) );
  }
  /*
         d -> 
         ^    ^
         |    |
    g -> c -> b ->
    ^    ^    ^    ^
    |    |    |    |
    f -> n == a -> e

    h -> n == a -> 
    ^    ^    ^    ^
    |    |    |    |
    g -> c -> b -> f
         ^    ^
         |    |
         d -> e
  */
  ptSU3 staple2x1(int n, int mu, int nu) {

    int curr = Z->L[n][4];
    int a,b,c,d,e,f,g,h;
    ptSU3 tmp;

    a = Z->L[curr][5+mu];
    b = Z->L[a][5+nu];
    c = Z->L[curr][5+nu];
    d = Z->L[c][5+nu];
    e = Z->L[a][5+mu];
    f = Z->L[curr][mu  ];
    g = Z->L[f][5+nu];

    tmp = ( W[a].U[nu] * W[b].U[nu] * dag(W[curr].U[nu] * W[c].U[nu] * W[d].U[mu] ) +
            W[a].U[mu] * W[e].U[nu] * dag(W[curr].U[nu] * W[c].U[mu] * W[b].U[mu] ) +
            W[a].U[nu] * dag(W[f].U[nu] * W[g].U[mu] * W[c].U[mu]) * W[f].U[mu] );
  
    b = Z->L[a][nu  ];
    c = Z->L[curr][nu  ];
    d = Z->L[c][nu  ];
    e = Z->L[b][nu  ];
    f = Z->L[b][5+mu];
    g = Z->L[c][mu  ];
    h = Z->L[curr][mu   ];

    tmp += ( dag(W[d].U[mu] * W[e].U[nu] * W[b].U[nu]) * W[d].U[nu] * W[c].U[nu]  +
             W[a].U[mu] * dag(W[c].U[mu] * W[b].U[mu] * W[f].U[nu]) * W[c].U[nu]  +
             dag(W[g].U[mu] * W[c].U[mu] * W[b].U[nu]) * W[g].U[nu] * W[h].U[mu] );

    return tmp;

  }



};


// ********* end class Gluon_fld **************




// SU2

class ptSU2_fld{
 private:
  latt *Z;
 public:
  ptSU2  *W;

  ptSU2_fld(latt* z) { Z=z; W = new ptSU2[Z->TSize]; }
  
  ptSU2_fld(FILE*,int); // SU2_fld(&input_file,read_mode)

  ~ptSU2_fld(){
    delete[] W;
  }

  int save(char *filename) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptSU2)*Z->Size);
    return 0;
  }


  int load(const char *filename) {
    if( io::load(filename,W) != Z->Size*sizeof(ptSU2) )
      return -1;
    return 0;

  }
  



  ptSU2* handle(){ return W; }

  ptSU2 get(int *);

  ptSU2 get(int n){
    return W[Z->L[n][4]];
  }
  
  ptSU2 get(int n,int step,int dir){
    return W[Z->get(n, step, dir)];
  }


  friend int get(ptSU2_fld *, int*);
  friend int get(ptSU2_fld *, int, int, int);
  friend int get(ptSU2_fld *, int);

};
  




inline ptSU2 ptSU2_fld::get(int *n) {
  return W[Z->get(n)];
}



inline int get(ptSU2_fld *W, int *n) {
  return (W->Z)->get(n);

}


inline int get(ptSU2_fld *W,int n,int step,int dir){
  return (W->Z)->get(n,step,dir);
}

inline int get(ptSU2_fld *W,int n){
  return (W->Z)->get(n);
}



// ------------- end class pt SU2_fld  --------------



class ptSU2Gluon_fld{
 private:

 public:
  latt   *Z;
  ptSU2Gluon  *W;
  
  ptSU2Gluon_fld(latt *z) { Z=z; W = new ptSU2Gluon[z->TSize];  }
  
  ~ptSU2Gluon_fld(){  delete [] W; }

  int save(char *filename) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptSU2Gluon)*Z->Size);
    if ( os.rdstate() != 0 )
      std::cerr << "Error writing " << filename << "\n";
    return 0;
  }


  int load(const char *filename) {
    if( io::load<ptSU2Gluon,false>(filename,W) != Z->Size*sizeof(ptSU2Gluon) )
      return -1;
    return 0;
  }
  
  // plaquette value is saves as well
  // as a check for correct configuration loading
  int save(char *filename, Cplx *placchetta) {
    std::ofstream os (filename, std::ofstream::binary);
    os.write(reinterpret_cast<char*>(W),sizeof(ptSU2Gluon)*Z->Size);
    os.write(reinterpret_cast<char*>(placchetta),sizeof(Cplx)*(allocORD+1));
    if ( os.rdstate() != 0 )
      std::cerr << "Error writing " << filename << "\n";
    return 0;
  }

  // plaquette value is read as well
  // as a check for correct configuration loading
  int load(const char *filename, Cplx *placchetta) {
    if( io::load<ptSU2Gluon>(filename,W,placchetta) != (Z->Size*sizeof(ptSU2Gluon)/*+(allocORD+1)*sizeof(Cplx)*/) )
      return -1;
    return 0;
  }

    int loadlowerorder(const char *filename, int loworder, Cplx *plaquette) {
      ptSU2 myapp;
        if(loworder>allocORD){
            cout<<"ERROR: trying to load more orders than the ones allocated."<<endl;
            return -1;
        }
        // ios::ate set the initial position at the end of the file, for a check on dimensions
        ifstream handle(filename, ios::binary | ios::ate);
        if (!handle.is_open()){
            cout<<"ERROR: file "<<filename<<" not found."<<endl;
            return -1;
        }
        if(handle.tellg()!=sizeof(double)*(2*(NC*NC*loworder+1)*dim*Z->Size+2*(loworder+1))){
            cout<<"ERROR: not possible to read the configuration, unexpected file dimension."<<endl;
            return -1;
        }
        handle.close();
        // open the file to read
        handle.open(filename, ios::binary);
        for (int x=0; x<Z->Size; x++){
            for (int mu=0; mu<dim; mu++){
                // read flag
                handle.read(reinterpret_cast<char*>( &(W[x].U[mu].flag.re) ), sizeof(double));
                handle.read(reinterpret_cast<char*>( &(W[x].U[mu].flag.im) ), sizeof(double));
                // read orders from file
                for (int k=0; k<loworder; k++){
                    for (int a=0; a<NC*NC; a++){
                        handle.read(reinterpret_cast<char*>( &(W[x].U[mu].ptU[k].whr[a].re) ), sizeof(double));
                        handle.read(reinterpret_cast<char*>( &(W[x].U[mu].ptU[k].whr[a].im) ), sizeof(double));
                    }
                }
                // set higher orders to zero
                for (int k=loworder+1; k<allocORD; k++){
                    for (int a=0; a<NC*NC; a++){
                        W[x].U[mu].ptU[k].whr[a].re = 0.;
                        W[x].U[mu].ptU[k].whr[a].im = 0.;
                    }
                }
                // the field has to be unitary
		myapp = log(W[x].U[mu]);
                myapp = myapp.reH();
                W[x].U[mu] = exp(myapp);
            }
        }
        // read plaquette
        for (int k=0; k<allocORD+1; k++){
            handle.read(reinterpret_cast<char*>( &(plaquette[k].re) ), sizeof(double));
            handle.read(reinterpret_cast<char*>( &(plaquette[k].im) ), sizeof(double));
        }
        return 0;
    }
    
  int load(const char *filename, Cplx *placchetta, int rank) {
    FILE *filept;
    if( (filept = fopen(filename, "r") ) == NULL ) return 1;

    long offset = Z->Sz[0]*rank*sizeof(ptSU2Gluon);
    fseek(filept, offset, SEEK_SET);
    fread(W,sizeof(ptSU2Gluon), Z->Size,filept);

    offset = -(allocORD+1)*sizeof(Cplx);
    fseek(filept, offset, SEEK_END);
    
    fread(placchetta, sizeof(Cplx), allocORD+1 , filept);

    fclose(filept);
    return 0;
  }
  


  
  void prout(){
    for(int i = 0; i < Z->Size; i++){
      printf("i = %d\n",i);
      W[i].prout();
      printf("\n\n");
    }
  }
  
  
  void operator=(const ptSU2Gluon_fld& A) {
    Z = A.Z;
    std::copy(A.W,A.W+Z->Size,W);
    /* for(int i = 0; i < Z->Size; i++) */
    /*   W[i] = A.W[i]; */
  }
  
  ptSU2* handle(){ return (ptSU2*)(W->U); }
  
  ptSU2Gluon& get(int* n) { return W[Z->get(n)]; }

  ptSU2Gluon& get(int n, int step, int dir) { return W[Z->get(n, step, dir)]; }

  ptSU2& get(int* n, int mu) { return W[Z->get(n)].U[mu]; }

  ptSU2& get(int n, int mu){ return W[Z->get(n)].U[mu]; }

  ptSU2& get(int n, int step, int dir, int mu){ return W[Z->get(n, step, dir)].U[mu]; }
  
  friend int get(ptSU2Gluon_fld *W, int n, int mu) { return ( 4*((W->Z)->get(n))+mu ); }

  friend int get(ptSU2Gluon_fld *W, int *n) { return W->Z->get(n); }

  friend int get(ptSU2Gluon_fld *W, int n) { return ( 4*((W->Z)->get(n)) ); }

  friend int get(ptSU2Gluon_fld *W, int n, int step , int dir, int mu) {
    return (4*((W->Z)->get(n,step,dir))+mu); 
  }

  friend int get(ptSU2Gluon_fld *W, int n, int step1 , int dir1,
		 int step2, int dir2, int mu) {
    return (4*((W->Z)->get( ((W->Z)->get(n, step1, dir1)), step2, dir2))+mu);
  }

  ptSU2Gluon& near(int n, int sign, int dir) {
    return W[Z->near(n, sign, dir)-1];
  }    


  /*
    b ->
    ^    ^
    |    |
    n == a

    n == 
    ^    ^
    |    |
    b -> a
  */
  ptSU2 staple(int*, int, int); //(x[dim], mu, nu)
  
  ptSU2 staple(int n, int mu, int nu) {
    return( W[Z->get(n, 1, mu)].U[nu]*dag(W[Z->L[n][4]].U[nu]*W[Z->get(n, 1, nu)].U[mu]) +
            dag(W[Z->get(n, -1, nu)].U[mu]*W[Z->get(n, -1, nu, 1, mu)].U[nu])*W[Z->get(n, -1, nu)].U[nu]);
  }

  ptSU2 singlestaple(int n, int mu, int nu) {
    return( W[Z->get(n, 1, mu)].U[nu]*dag(W[Z->L[n][4]].U[nu]*W[Z->get(n, 1, nu)].U[mu]) );
  }
  
  /*
         d -> 
         ^    ^
         |    |
    g -> c -> b ->
    ^    ^    ^    ^
    |    |    |    |
    f -> n == a -> e

    h -> n == a -> 
    ^    ^    ^    ^
    |    |    |    |
    g -> c -> b -> f
         ^    ^
         |    |
         d -> e
  */
  ptSU2 staple2x1(int n, int mu, int nu) {

    int curr = Z->L[n][4];
    int a,b,c,d,e,f,g,h;
    ptSU2 tmp;

    a = Z->L[curr][5+mu];
    b = Z->L[a][5+nu];
    c = Z->L[curr][5+nu];
    d = Z->L[c][5+nu];
    e = Z->L[a][5+mu];
    f = Z->L[curr][mu  ];
    g = Z->L[f][5+nu];

    tmp = ( W[a].U[nu] * W[b].U[nu] * dag(W[curr].U[nu] * W[c].U[nu] * W[d].U[mu] ) +
            W[a].U[mu] * W[e].U[nu] * dag(W[curr].U[nu] * W[c].U[mu] * W[b].U[mu] ) +
            W[a].U[nu] * dag(W[f].U[nu] * W[g].U[mu] * W[c].U[mu]) * W[f].U[mu] );
  
    b = Z->L[a][nu  ];
    c = Z->L[curr][nu  ];
    d = Z->L[c][nu  ];
    e = Z->L[b][nu  ];
    f = Z->L[b][5+mu];
    g = Z->L[c][mu  ];
    h = Z->L[curr][mu   ];

    tmp += ( dag(W[d].U[mu] * W[e].U[nu] * W[b].U[nu]) * W[d].U[nu] * W[c].U[nu]  +
             W[a].U[mu] * dag(W[c].U[mu] * W[b].U[mu] * W[f].U[nu]) * W[c].U[nu]  +
             dag(W[g].U[mu] * W[c].U[mu] * W[b].U[nu]) * W[g].U[nu] * W[h].U[mu] );

    return tmp;

  }


};


// ********* end class SU2Gluon_fld **************





#endif
