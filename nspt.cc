/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: nspt.cc
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
/*  END LEGAL NOTICE*/

#include "nspt.h"
#include "twistmatrices.h"

ptSUX W, W1, W2, Wsym;
#ifdef REMOVE_GLUON_ZEROMOM
ptSUX *MM = new ptSUX[dim];
#endif

int link_c, site_c, site_up, xx[4], curr;
const int offset = dim+1;

extern Cplx* w;  // Will be useful ...
extern Cplx* w1; // Useful for 1x1 plaquette
extern Cplx* w2; // Useful for 2x1 plaquette
extern double* norm;
extern double* norm1;
extern Cplx** trU;
extern int ratio; 
extern ptSUX* U;
#ifndef UPDATE_ONTHEFLY
extern ptSUX* F;
#endif
extern nspt_params_t nspt_pars;
extern act_params_t  act_pars;
extern thread_params_t thr_pars[NTHR];

extern std::ofstream plaqfile, normfile, logfile, trufile;

#ifdef __WILLOOP_2x2__
extern std::ofstream loop2x2;
#endif
#ifdef __MANY_LOOP__
extern std::ofstream loop_nxm;
#endif

#ifdef USE_RAND55
extern MyRand Rand[NTHR];
#else
std::vector<ranlxd::Rand> Rand;
#endif
#ifdef __PARALLEL_OMP__
ptSUX Ww[NTHR], Ww1[NTHR], Ww2[NTHR], Wwsym[NTHR];
#ifdef REMOVE_GLUON_ZEROMOM
ptSUX MMm[NTHR][dim];
#endif
SUX P[NTHR];
Cplx ww[NTHR][allocORD+1], ww1[NTHR][allocORD+1], ww2[NTHR][allocORD+1];
double normm[NTHR][allocORD+1];
Cplx ttrU[NTHR][dim][allocORD];
int tid;

#endif

// for (single step of) Fourier Accelerated gauge fixing
int coord[dim];
ptSUX_fld* Wgauge;
fftw_plan *planFA;

// for momentum space analysis of the norm
fftw_plan *planGluon;

#ifdef __K_MOM_ANALYSIS__
double *knorm[allocORD];
extern FILE *FPkn;
#endif

#ifdef __TIMING__
extern PRlgtTime Time;
#endif


#ifdef TWISTED_BC
#ifdef _TWIST_X
ptSUX twist_x(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_x*A.ptU[i]*omega_xdag;
    }
    return B;
}
ptSUX invtwist_x(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_xdag*A.ptU[i]*omega_x;
    }
    return B;
}
#endif
#ifdef _TWIST_Y
ptSUX twist_y(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_y*A.ptU[i]*omega_ydag;
    }
    return B;
}
ptSUX invtwist_y(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_ydag*A.ptU[i]*omega_y;
    }
    return B;
}
#endif
#ifdef _TWIST_Z
ptSUX twist_z(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_z*A.ptU[i]*omega_zdag;
    }
    return B;
}
ptSUX invtwist_z(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_zdag*A.ptU[i]*omega_z;
    }
    return B;
}
#endif
#ifdef _TWIST_T
ptSUX twist_t(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_t*A.ptU[i]*omega_tdag;
    }
    return B;
}
ptSUX invtwist_t(const ptSUX &A){
    ptSUX B;
    for (int i=0; i < PTORD; i++) {
        B.ptU[i] = omega_tdag*A.ptU[i]*omega_t;
    }
    return B;
}
#endif

void twist_boundary(ptSUXGluon_fld& Umu){
#ifdef _TWIST_X
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistF_x; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistF_x[i][1]].U[mu] = twist_x(Umu.W[Umu.Z->copyandtwistF_x[i][0]].U[mu]);
        }
    }
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistB_x; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistB_x[i][1]].U[mu] = invtwist_x(Umu.W[Umu.Z->copyandtwistB_x[i][0]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_Y
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistF_y; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistF_y[i][1]].U[mu] = twist_y(Umu.W[Umu.Z->copyandtwistF_y[i][0]].U[mu]);
        }
    }
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistB_y; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistB_y[i][1]].U[mu] = invtwist_y(Umu.W[Umu.Z->copyandtwistB_y[i][0]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_Z
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistF_z; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistF_z[i][1]].U[mu] = twist_z(Umu.W[Umu.Z->copyandtwistF_z[i][0]].U[mu]);
        }
    }
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistB_z; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistB_z[i][1]].U[mu] = invtwist_z(Umu.W[Umu.Z->copyandtwistB_z[i][0]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_T
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistF_t; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistF_t[i][1]].U[mu] = twist_t(Umu.W[Umu.Z->copyandtwistF_t[i][0]].U[mu]);
        }
    }
#pragma omp parallel for
    for (int i=0; i<Umu.Z->sizecopyandtwistB_t; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->copyandtwistB_t[i][1]].U[mu] = invtwist_t(Umu.W[Umu.Z->copyandtwistB_t[i][0]].U[mu]);
        }
    }
#endif
    
    // additional twist
    
#ifdef _TWIST_X
    #pragma omp parallel for
    for (int i=0; i<Umu.Z->sizejusttwist[0]; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->justtwist[0][i]].U[mu] = twist_x(Umu.W[Umu.Z->justtwist[0][i]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_Y
    #pragma omp parallel for
    for (int i=0; i<Umu.Z->sizejusttwist[1]; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->justtwist[1][i]].U[mu] = twist_y(Umu.W[Umu.Z->justtwist[1][i]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_Z
    #pragma omp parallel for
    for (int i=0; i<Umu.Z->sizejusttwist[2]; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->justtwist[2][i]].U[mu] = twist_z(Umu.W[Umu.Z->justtwist[2][i]].U[mu]);
        }
    }
#endif
#ifdef _TWIST_T
    #pragma omp parallel for
    for (int i=0; i<Umu.Z->sizejusttwist[3]; i++) {
        for (int mu=0; mu<4; mu++) {
            Umu.W[Umu.Z->justtwist[3][i]].U[mu] = twist_t(Umu.W[Umu.Z->justtwist[3][i]].U[mu]);
        }
    }
#endif
}
#endif

#ifdef UPDATE_ONTHEFLY
void gauge_wilson(ptSUXGluon_fld& Umu){
#else
void gauge_wilson(ptSUXGluon_fld& Umu, ptSUXGluon_fld& Fmu){
#endif
    
#ifndef __PARALLEL_OMP__
  for(curr = 0; curr < act_pars.iVol; curr++){
    
    for(int mu = 0; mu < dim; mu++){
        
      // set up current link
      link_c = get(&Umu, curr, mu);
      W1.zero();	
#if GAUGE_ACTION != WIL
      W2.zero();
#endif				
      // Compute staple, possible 2x1 as well
      // for each relevant plane
      for(int nu = 0; nu < dim; nu++){		
	if(nu != mu){					
	  W1 += Umu.staple(curr, mu, nu);			
#if GAUGE_ACTION != WIL
	  W2 += Umu.staple2x1(curr, mu, nu);			
#endif
	}						
      }						
    
        
      // From staple to plaquette, then compute trace
      W1 = U[link_c]*W1;
      W1.Tr(w);
      // Sum plaq value and compute norm
      w1[0] += w[0];
      for (int i1=1; i1 <= PTORD; i1++){		
	w1[i1] += w[i1];
	norm[i1-1] += ((U[link_c].ptU[i1-1]*dag(U[link_c].ptU[i1-1])).Tr()).re;
	trU[mu][i1-1] += U[link_c].ptU[i1-1].Tr();
      }


      // idem for improved actions
#if GAUGE_ACTION != WIL
      W2 = U[link_c]*W2;			
      W2.Tr(w);
      for (int i1=0; i1 <= PTORD; i1++){		
	w2[i1] += w[i1];
      }

      W1 = act_pars.c0*W1 + act_pars.c1*W2;
#endif
        
      // Stay in the algebra!
      W1 = W1.reH();					
      W1 *= act_pars.tau_g;

        
      // Add (Langevin) fluctuation 
#ifdef _runSU2
      W1.ptU[0] -= act_pars.stau*SU2rand(Rand[0]);
#else
      W1.ptU[0] -= act_pars.stau*SU3rand(Rand[0]);
#endif

      U[link_c] = exp(W1)*U[link_c];
        
      // Zero momentum contribution (in case...)
#ifdef REMOVE_GLUON_ZEROMOM
#ifdef UPDATE_ONTHEFLY
        MM[mu] += log(U[link_c]);
#else
        MM[mu] += log(F[link_c]);
#endif
#endif
    } //end mu
    
  } // all sites done


#else
  
  for(tid = 0;tid < NTHR; tid++) {
    for(int ord = 0; ord < PTORD; ord++){
      ww[tid][ord]  = 0;
      ww1[tid][ord] = 0;
      ww2[tid][ord] = 0;
      normm[tid][ord] = 0;
      for( int mu = 0; mu < dim; mu++)
	{
	  ttrU[tid][mu][ord] = Cplx(0,0);
	}
    }
    ww[tid][PTORD]  = 0;
    ww1[tid][PTORD] = 0;
    ww2[tid][PTORD] = 0;
    
#ifdef REMOVE_GLUON_ZEROMOM
    for(int mu = 0; mu < dim; mu++){
      MMm[tid][mu].zero();
    }
#endif
  }
  
  
#pragma omp parallel private(curr,tid,link_c) num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );


	    for(int mu = 0; mu < dim; mu++){
	      // set up current link
	      link_c = get(&Umu, curr, mu);
      
	      Ww1[tid].zero();	
#if GAUGE_ACTION != WIL
	      Ww2[tid].zero();
#endif				
	      // Compute staple, possible 2x1 as well
	      // for each relevant plane
	      for(int nu = 0; nu < dim; nu++){		
		if(nu != mu){
		  Ww1[tid] += Umu.staple(curr, mu, nu);
#if GAUGE_ACTION != WIL
		  Ww2[tid] += Umu.staple2x1(curr, mu, nu);			
#endif
		}						
	      }						
	      // From staple to plaquette, then compute trace
	      Ww1[tid] = U[link_c]*Ww1[tid];
	      Ww1[tid].Tr(ww[tid]);
	      // Sum plaq value and compute norm
	      ww1[tid][0] += ww[tid][0];
	      for (int i1=1; i1 <= PTORD; i1++){
		ww1[tid][i1] += ww[tid][i1];

		normm[tid][i1-1] += ((U[link_c].ptU[i1-1]*dag(U[link_c].ptU[i1-1])).Tr()).re;

		ttrU[tid][mu][i1-1] += U[link_c].ptU[i1-1].Tr();
	      }

	      // idem for improved actions
#if GAUGE_ACTION != WIL
	      Ww2[tid] = U[link_c]*Ww2[tid];		
	      Ww2[tid].Tr(ww[tid]);
	      for (int i1=0; i1 <= PTORD; i1++){		
		ww2[tid][i1] += ww[tid][i1];
	      }
      
	      Ww1[tid] = act_pars.c0*Ww1[tid] + act_pars.c1*Ww2[tid];
#endif
            
	      // Stay in teh algebra!
	      Ww1[tid]  = Ww1[tid].reH();
	      Ww1[tid] *= act_pars.tau_g;
      
	      // Add (Langevin) fluctuation
#ifdef _runSU2
	      Ww1[tid].ptU[0] -= act_pars.stau*SU2rand(Rand[tid]);
#else
          Ww1[tid].ptU[0] -= act_pars.stau*SU3rand(Rand[tid]);
#endif
            
#ifdef UPDATE_ONTHEFLY
	      U[link_c] = exp(Ww1[tid])*U[link_c];
#else
          F[link_c] = exp(Ww1[tid])*U[link_c];
#endif
            
#ifdef REMOVE_GLUON_ZEROMOM
	      // Zero momentum contribution (in case...)
#ifdef UPDATE_ONTHEFLY
            MMm[tid][mu] += log(U[link_c]);
#else
            MMm[tid][mu] += log(F[link_c]);
#endif
#endif
      
	    } //end mu
	    
#if ntz > 1
#pragma omp barrier
#endif
	  } // end z
	  
#if nty > 1
#pragma omp barrier
#endif
	} // end y
	
#if ntx > 1
#pragma omp barrier
#endif
      } // end x
      
#if ntt > 1
#pragma omp barrier
#endif
    }  // end t

#pragma omp barrier  
  } // end parallel

  // Reduce from parallel
  for(tid = 0; tid < NTHR; tid++){
  
    for(int ord = 0; ord < PTORD; ord++){
      w1[ord] += ww1[tid][ord];
      w2[ord] += ww2[tid][ord];
      norm[ord] += normm[tid][ord];
      for(int mu = 0; mu < dim; mu++){
	trU[mu][ord] += ttrU[tid][mu][ord];
      }
    }
    w1[PTORD] += ww1[tid][PTORD];
    w2[PTORD] += ww2[tid][PTORD];
      
#ifdef REMOVE_GLUON_ZEROMOM
    for(int mu = 0; mu < dim; mu++){
      MM[mu] += MMm[tid][mu];
    }
#endif
  }

#endif
    

#ifdef SELECTED_PLAQ_MEASURE
    int plaqcounter = 0;
    for (int i1=0; i1 <= PTORD; i1++) w1[i1] = 0.0;
    for (int y0=0; y0<Umu.Z->Sz[0]; y0++)
        for (int y1=0; y1<Umu.Z->Sz[1]; y1++)
            for (int y2=0; y2<Umu.Z->Sz[2]; y2++)
                for (int y3=0; y3<Umu.Z->Sz[3]; y3++){
                    
                    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );
/////////////////////
                    W1.zero();
                    
                        int mu = 1;
                        int nu = 2;
                        W1 += Umu.singlestaple(curr, mu, nu);
                        plaqcounter++;
                    
                        mu = 1;
                        nu = 3;
                        W1 += Umu.singlestaple(curr, mu, nu);
                        plaqcounter++;
                    
                    link_c = get(&Umu, curr, mu);
                    W1 = U[link_c]*W1;
                    W1.Tr(w);
                    
                    for (int i1=0; i1 <= PTORD; i1++) w1[i1] += w[i1];
                    
                    W1.zero();
                    
                        mu = 2;
                        nu = 3;
                        W1 += Umu.singlestaple(curr, mu, nu);
                        plaqcounter++;
                    
                    link_c = get(&Umu, curr, mu);
                    W1 = U[link_c]*W1;
                    W1.Tr(w);
                    
                    for (int i1=0; i1 <= PTORD; i1++) w1[i1] += w[i1];
/////////////////////
    }
#ifdef _runSU2
    for (int i1=0; i1 <= PTORD; i1++) w1[i1] *= 48.*act_pars.iVol/(2.*plaqcounter);
#else
    for (int i1=0; i1 <= PTORD; i1++) w1[i1] *= 72.*act_pars.iVol/(3.*plaqcounter);
#endif
#endif
    
    
#ifdef _runSU2
  for (int i1=0; i1 < PTORD; i1++){
#if GAUGE_ACTION == WIL
    plaqfile << w1[i1].re/(double)(48*act_pars.iVol) << "\t"
	     << w1[i1].im/(double)(48*act_pars.iVol) << std::endl;
#else
    plaqfile << w1[i1].re/(double)(48*act_pars.iVol) << "\t"
	     << w1[i1].im/(double)(48*act_pars.iVol) << "\t"
	     << w2[i1].re/(double)(48*act_pars.iVol) << "\t"
	     << w2[i1].im/(double)(48*act_pars.iVol) << std::endl;
#endif
    normfile << norm[i1]/(double)(8*act_pars.iVol) << std::endl;
    for( int mu = 0; mu < dim; mu++)
      {
	trufile  << trU[mu][i1].re /(double)(2*act_pars.iVol) << "\t"
		 << trU[mu][i1].im /(double)(2*act_pars.iVol) << "\t";
      }
    trufile  << std::endl;
  }

#if GAUGE_ACTION == WIL
  plaqfile << w1[PTORD].re/(double)(48*act_pars.iVol) << "\t"
	   << w1[PTORD].im/(double)(48*act_pars.iVol) << std::endl;
#else
  plaqfile << w1[PTORD].re/(double)(48*act_pars.iVol) << "\t"
	   << w1[PTORD].im/(double)(48*act_pars.iVol) << "\t"
	   << w2[PTORD].re/(double)(48*act_pars.iVol) << "\t"
	   << w2[PTORD].im/(double)(48*act_pars.iVol) << std::endl;
#endif
#else // else ... SU3
    for (int i1=0; i1 < PTORD; i1++){
#if GAUGE_ACTION == WIL
        plaqfile << w1[i1].re/(double)(72*act_pars.iVol) << "\t"
        << w1[i1].im/(double)(72*act_pars.iVol) << std::endl;
#else
        plaqfile << w1[i1].re/(double)(72*act_pars.iVol) << "\t"
        << w1[i1].im/(double)(72*act_pars.iVol) << "\t"
        << w2[i1].re/(double)(72*act_pars.iVol) << "\t"
        << w2[i1].im/(double)(72*act_pars.iVol) << std::endl;
#endif
        normfile << norm[i1]/(double)(12*act_pars.iVol) << std::endl;
        for( int mu = 0; mu < dim; mu++)
        {
            trufile  << trU[mu][i1].re /(double)(3*act_pars.iVol) << "\t"
            << trU[mu][i1].im /(double)(3*act_pars.iVol) << "\t";
        }
        trufile  << std::endl;
    }
    
#if GAUGE_ACTION == WIL
    plaqfile << w1[PTORD].re/(double)(72*act_pars.iVol) << "\t"
	   << w1[PTORD].im/(double)(72*act_pars.iVol) << std::endl;
#else
    plaqfile << w1[PTORD].re/(double)(72*act_pars.iVol) << "\t"
	   << w1[PTORD].im/(double)(72*act_pars.iVol) << "\t"
	   << w2[PTORD].re/(double)(72*act_pars.iVol) << "\t"
	   << w2[PTORD].im/(double)(72*act_pars.iVol) << std::endl;
#endif
#endif
    
} // gluonic evolution done!




#ifdef REMOVE_GLUON_ZEROMOM
void zero_modes_subtraction(ptSUXGluon_fld& Umu){
  // Normalizzo il modo nullo
  for (int i1=0; i1 < dim; i1++){
    MM[i1] *= act_pars.rVol;
  }

  // We will deal with the algebra field to subtract zero modes
  // The way we visit the lattice is pointless

  #ifndef __PARALLEL_OMP__
  for(int i = 0; i < act_pars.iVol; i++){
    for(int mu = 0; mu < dim; mu++){
      W = log(Umu.W[i].U[mu]);
      W -= MM[mu];
      W.reH();
      Umu.W[i].U[mu]  = exp(W);
    }
  }
  #else
    
  int curr;
#pragma omp parallel private(curr, tid) num_threads(NTHR) 
  {
    tid = omp_get_thread_num();
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );

	    for(int mu = 0; mu < dim; mu++){
	      Ww[tid] = log(Umu.W[curr].U[mu]);
	      Ww[tid] -= MM[mu];
	      Ww[tid].reH();
	      Umu.W[curr].U[mu] = exp(Ww[tid]);
	    } // mu

	  } // y3
	} // y2
      } // y1
    } // y0
    
  } // parallel

#endif

}
#endif





void stochastic_gauge_fixing(ptSUXGluon_fld& Umu){

#ifndef __PARALLEL_OMP__
  for(curr = 0; curr < act_pars.iVol; curr++){

    // set up current link
    site_c = Umu.Z->get(curr);
    W.zero();
    
    // compute DmuUmu and constrain it in the algebra
    for(int mu = 0; mu < dim; mu++){
      W += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
    }
    W -= dag(W);
    
//     Yet another version...
//     for(int mu = 0; mu < dim; mu++){
//       W += log(Umu.W[site_c].U[mu])- log(Umu.W[Umu.Z->L[curr][mu]].U[mu]);
//     }

    W.Trless();
      
    // Set up gauge transformations
    W1 = exp( act_pars.alpha*W);
    W2 = exp(-act_pars.alpha*W);
    
    // Action of GT on relevant links 
      int mu=0;
      Umu.W[site_c].U[mu] = W1*Umu.W[site_c].U[mu];
#ifdef _TWIST_X
      if ((int)site_c/Umu.Z->sizeYZT==0) {
          Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu] = Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu]*twist_x(W2);
      }
      else
#endif
          Umu.W[Umu.Z->L[site_c][mu]].U[mu] = Umu.W[Umu.Z->L[site_c][mu]].U[mu]*W2;
      
      mu=1;
      Umu.W[site_c].U[mu] = W1*Umu.W[site_c].U[mu];
#ifdef _TWIST_Y
      if(((int)site_c/Umu.Z->sizeZT)%Umu.Z->Sz[1]==0){
          Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu] = Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu]*twist_y(W2);
      }
      else
#endif
          Umu.W[Umu.Z->L[site_c][mu]].U[mu] = Umu.W[Umu.Z->L[site_c][mu]].U[mu]*W2;
      
      mu=2;
      Umu.W[site_c].U[mu] = W1*Umu.W[site_c].U[mu];
#ifdef _TWIST_Z
      if(((int)site_c/Umu.Z->sizeT)%Umu.Z->Sz[2]==0){
          Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu] = Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu]*twist_z(W2);
      }
      else
#endif
          Umu.W[Umu.Z->L[site_c][mu]].U[mu] = Umu.W[Umu.Z->L[site_c][mu]].U[mu]*W2;
      
      mu=3;
      Umu.W[site_c].U[mu] = W1*Umu.W[site_c].U[mu];
#ifdef _TWIST_T
      if (site_c%Umu.Z->Sz[3]==0) {
          Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu] = Umu.W[Umu.Z->get(site_c,-2,mu)].U[mu]*twist_t(W2);
      }
      else
#endif
          Umu.W[Umu.Z->L[site_c][mu]].U[mu] = Umu.W[Umu.Z->L[site_c][mu]].U[mu]*W2;
    
  } // loop on sites
  
#else

#pragma omp parallel private(curr,tid,site_c)  num_threads(NTHR)
  {
    tid = omp_get_thread_num();
    
    for(int y0 = thr_pars[tid].xi[0]; y0 < thr_pars[tid].xf[0]; y0++){
      for(int y1 = thr_pars[tid].xi[1]; y1 < thr_pars[tid].xf[1]; y1++){
	for(int y2 = thr_pars[tid].xi[2]; y2 < thr_pars[tid].xf[2]; y2++){
	  for(int y3 = thr_pars[tid].xi[3]; y3 < thr_pars[tid].xf[3]; y3++){

	    curr = y3 + act_pars.sz[3]*(y2 + act_pars.sz[2]*(y1 + act_pars.sz[1]*y0) );
	    
	    site_c = Umu.Z->L[curr][dim];
	    
	    // compute DmuUmu and constrain it in the algebra
	    Ww[tid].zero();
	    for(int mu = 0; mu < dim; mu++){
	      Ww[tid] += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
	    }

	    // yet another version ...
	    // Ww[tid] = Umu.W[site_c].U[0] - Umu.W[Umu.Z->L[curr][0]].U[0];
	    // for(int mu = 1; mu < dim; mu++){
	    //   Ww[tid] += Umu.W[site_c].U[mu] - Umu.W[Umu.Z->L[curr][mu]].U[mu];
	    // }

	    Ww[tid] -= dag(Ww[tid]);
	    Ww[tid].Trless();
	    
	    // Set up gauge transformations
	    Ww1[tid] = exp( act_pars.alpha*Ww[tid]);
	    Ww2[tid] = exp(-act_pars.alpha*Ww[tid]);
	    
	    // Action of GT on relevant links
          int mu=0;
          Umu.W[site_c].U[mu] = Ww1[tid]*Umu.W[site_c].U[mu];
#ifdef _TWIST_X
          if (y0==0) {
              Umu.W[Umu.Z->get(curr,-2,mu)].U[mu] = Umu.W[Umu.Z->get(curr,-2,mu)].U[mu]*twist_x(Ww2[tid]);
          }
          else
#endif
              Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*Ww2[tid];
          
          
          mu=1;
          Umu.W[site_c].U[mu] = Ww1[tid]*Umu.W[site_c].U[mu];
#ifdef _TWIST_Y
          if (y1==0) {
              Umu.W[Umu.Z->get(curr,-2,mu)].U[mu] = Umu.W[Umu.Z->get(curr,-2,mu)].U[mu]*twist_y(Ww2[tid]);
          }
          else
#endif
              Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*Ww2[tid];
          
              
          mu=2;
          Umu.W[site_c].U[mu] = Ww1[tid]*Umu.W[site_c].U[mu];
#ifdef _TWIST_Z
          if (y2==0) {
              Umu.W[Umu.Z->get(curr,-2,mu)].U[mu] = Umu.W[Umu.Z->get(curr,-2,mu)].U[mu]*twist_z(Ww2[tid]);
          }
          else
#endif
              Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*Ww2[tid];

          
          mu=3;
          Umu.W[site_c].U[mu] = Ww1[tid]*Umu.W[site_c].U[mu];
#ifdef _TWIST_T
          if (y3==0) {
              Umu.W[Umu.Z->get(curr,-2,mu)].U[mu] = Umu.W[Umu.Z->get(curr,-2,mu)].U[mu]*twist_t(Ww2[tid]);
          }
          else
#endif
              Umu.W[Umu.Z->L[curr][mu]].U[mu] = Umu.W[Umu.Z->L[curr][mu]].U[mu]*Ww2[tid];

#if ntz > 1
#pragma omp barrier
#endif
	  } // end z
	  
#if nty > 1
#pragma omp barrier
#endif
	} // end y
	
#if ntx > 1
#pragma omp barrier
#endif
      } // end x
      
#if ntt > 1
#pragma omp barrier
#endif
    }  // end t


#pragma omp barrier
  } // end parallel
#endif


#ifdef __K_MOM_ANALYSIS__

  fftw_execute(planGluon[0]);
  for(curr = 0; curr < act_pars.iVol; curr++){
    for (int i1=0; i1 < allocORD; i1++){		
      knorm[i1][curr] = 0;
    }
    for(int mu = 0; mu < dim; mu++){
      for (int i1=0; i1 < allocORD; i1++){		
	knorm[i1][curr] += ( Umu.W[curr].U[mu].ptU[i1]*
			     dag(Umu.W[curr].U[mu].ptU[i1]) 
			     ).Tr().re;
      }
    }
  }

  fftw_execute(planGluon[1]);
  
  for(curr = 0; curr < act_pars.iVol; curr++){
    for(int mu = 0; mu < dim; mu++){
      Umu.W[curr].U[mu] *= act_pars.rVol;
    }
  }

  fwrite(knorm[0],sizeof(double),act_pars.iVol*allocORD ,FPkn);

#endif


}



void FAstochastic_gauge_fixing(ptSUXGluon_fld& Umu){

#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    Wgauge->W[Umu.Z->get(x)].zero();
    for(int mu = 0; mu < dim; mu++){
      Wgauge->W[Umu.Z->get(x)] += ( log(U[get(&Umu,x,mu)]) - 
				    log(U[get(&Umu,x,-1,mu,mu)]) );
    }
  }
  
  fftw_execute(planFA[0]);
  
#ifdef PBC
  Wgauge->W[0] *= 0.0;
#endif  
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 1; x < act_pars.iVol; x++){
    Umu.Z->get(x,coord);
    Wgauge->W[Umu.Z->get(x)] *= -(act_pars.alpha/(Umu.Z->p2hat[0][coord[0]]+
						  Umu.Z->p2hat[1][coord[1]]+
						  Umu.Z->p2hat[2][coord[2]]+
						  Umu.Z->p2hat[3][coord[3]]));
  }
  
  
  fftw_execute(planFA[1]);
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    Wgauge->W[x] /= (double)act_pars.iVol;
    Wgauge->W[x].reH();
    Wgauge->W[x] = exp(Wgauge->W[x]);
  }
    
  
#ifdef  __PARALLEL_OMP__
#pragma omp for schedule(static)
#endif
  for(int x = 0; x < act_pars.iVol; x++){
    for(int mu = 0; mu < dim; mu++){
      U[get(&Umu,x,mu)] = Wgauge->get(x)*Umu.get(x,mu)*dag(Wgauge->get(x,1,mu));
    }
  }
  

#ifdef __K_MOM_ANALYSIS__

  fftw_execute(planGluon[0]);

  for(curr = 0; curr < act_pars.iVol; curr++){
    for (int i1=0; i1 < allocORD; i1++){		
      knorm[i1][curr] = 0;
    }
    for(int mu = 0; mu < dim; mu++){
      for (int i1=0; i1 < allocORD; i1++){		
	knorm[i1][curr] += ( Umu.W[curr].U[mu].ptU[i1]*
			     dag(Umu.W[curr].U[mu].ptU[i1]) 
			     ).Tr().re;
      }
    }
  }

  fftw_execute(planGluon[1]);
  
  for(curr = 0; curr < act_pars.iVol; curr++){
    for(int mu = 0; mu < dim; mu++){
      Umu.W[curr].U[mu] *= act_pars.rVol;
    }
  }

  fwrite(knorm[0],sizeof(double),act_pars.iVol*allocORD ,FPkn);

#endif
  
}



// Beat step (quenched evolution); a funny legacy name ...
#ifdef UPDATE_ONTHEFLY
void NsptEvolve(ptSUXGluon_fld& Umu){
#else
void NsptEvolve(ptSUXGluon_fld& Umu, ptSUXGluon_fld& Fmu){
#endif
    
#ifdef __TIMING__
  Time.reset();
#endif

  for( int t2 = 0; t2 < nspt_pars.Beat; t2++ ){
    
    // Set a few variables to zero
#ifdef REMOVE_GLUON_ZEROMOM
    for (int i1=0; i1 < dim; i1++){
      MM[i1].zero();
    }
#endif
      for (int i1=0; i1 < PTORD; i1++){
          w1[i1]    = 0.0;
          norm[i1]  = 0.0;
          for (int mu=0; mu < dim; mu++){
              trU[mu][i1]   = 0.0;
          }
#if GAUGE_ACTION != WIL
          w2[i1]    = 0.0;
#endif
      }
      w1[PTORD] = 0.0;

    // gluonic evolution
#ifdef __TIMING__
    Time.tic_g();
#endif

#ifdef UPDATE_ONTHEFLY
    gauge_wilson(Umu);
#else
      gauge_wilson(Umu,Fmu);
      Umu = Fmu; // Very inefficient...
#ifdef TWISTED_BC
      twist_boundary(Umu);
#endif
#endif
      
#ifdef __TIMING__
    Time.toc_g();
#endif

    // Sottrazione dei momenti nulli
#ifdef __TIMING__
    Time.tic_zm();
#endif
#ifdef REMOVE_GLUON_ZEROMOM
    zero_modes_subtraction(Umu);
#endif
#ifdef __TIMING__
    Time.toc_zm();
#endif

    // Stochastic Gauge fixing 
#ifdef __TIMING__
    Time.tic_gf();
#endif
#ifndef __FA_GAUGE_FIXING__
    stochastic_gauge_fixing(Umu);
#else
    // Fourier acceleration
    FAstochastic_gauge_fixing(Umu);
#endif
      
#ifdef TWISTED_BC
      twist_boundary(Umu);
#endif
      
#ifdef __TIMING__
    Time.toc_gf();
    Time.reduce();
#endif

  }

#ifdef __TIMING__
  Time.out();
#endif

}






// Compare two plaquette values
// Zero is returned if they are equal; otherwise we return the difference
int plaquette_check(Cplx* w1, Cplx* w2){
  // return ( memcmp(w1, w2, (PTORD+1)*sizeof(Cplx)) );
#ifdef LOADLOWERORDER
  for(int i1 = 0; i1 <= LOADLOWERORDER;i1++){
#else
  for(int i1 = 0; i1 <= PTORD;i1++){
#endif
    if( (w1[i1]-w2[i1]).mod() > 1e-7 ) return 1;
  }
  return 0;
}





#ifdef __WILLOOP_2x2__

void WL2x2(ptSUXGluon_fld& Umu) {
  
#ifndef __PARALLEL_OMP__
  ptSUX tmp, tmp1;
#else
  ptSUX tmp[NTHR], tmp1[NTHR];
#endif
  
  int a,b,d,e;
  

#ifndef __PARALLEL_OMP__
  for(int site_x = 0; site_x < act_pars.iVol; site_x++) {

    for(int mu = 3; mu > 0; --mu){
      a = Umu.Z->L[site_x][5+mu];
      b = Umu.Z->L[a][5+mu];
      tmp = Umu.W[site_x].U[mu] * Umu.W[a].U[mu];

      for( int nu = 0; nu < mu; ++nu){
	d = Umu.Z->L[site_x][5+nu];
	e = Umu.Z->L[d][5+nu];
	tmp1 += tmp * ( Umu.W[b].U[nu] * Umu.W[ Umu.Z->L[b][5+nu] ].U[nu] * 
			dag( Umu.W[site_x].U[nu] * Umu.W[d].U[nu] * Umu.W[e].U[mu] * 
			     Umu.W[ Umu.Z->L[e][5+mu] ].U[mu] ));

      } // nu
    } // mu
  } // site

  for( int i1 = 1; i1 < PTORD; i1+=2){
    loop2x2 << (tmp1.ptU[i1].Tr()).re*act_pars.rVol*D18 << "\t";
  }
  loop2x2 << std::endl;

#else
#pragma omp parallel private(tid,a,b,d,e)  num_threads(NTHR) 
  {
    int chunk = act_pars.iVol/NTHR;
    tid = omp_get_thread_num();
    for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
      
      for(int mu = 3; mu > 0; --mu){
	a = Umu.Z->L[site_x][5+mu];
	b = Umu.Z->L[a][5+mu];
	tmp[tid] = Umu.W[site_x].U[mu] * Umu.W[a].U[mu];
	
	for( int nu = 0; nu < mu; ++nu){
	  d = Umu.Z->L[site_x][5+nu];
	  e = Umu.Z->L[d][5+nu];
	  tmp1[tid] += tmp[tid] * ( Umu.W[b].U[nu] * Umu.W[ Umu.Z->L[b][5+nu] ].U[nu] * 
				    dag( Umu.W[site_x].U[nu] * Umu.W[d].U[nu] * Umu.W[e].U[mu] * 
					 Umu.W[ Umu.Z->L[e][5+mu] ].U[mu] ));
	} // nu
      } // mu
      
    } // end site_x
    
  } // end parallel
  

  // Reduce
  for(int tid = 1; tid < NTHR; tid++){
    tmp1[0] += tmp1[tid];
  }

  // std::cout << std::endl;
  for( int i1 = 1; i1 <= PTORD; i1+=2){
    // std::cout << (tmp1[0].ptU[i1].Tr()).re*act_pars.rVol*D18 << "\t";
    loop2x2 << (tmp1[0].ptU[i1].Tr()).re*act_pars.rVol*D18 << "\t";
  }

  // std::cout << std::endl;
  loop2x2 << std::endl;

#endif

}

#endif


#ifdef __MANY_LOOP__

Cplx l21[allocORD+1];
Cplx l22[allocORD+1];
Cplx l42[allocORD+1];
Cplx l44[allocORD+1];


void ComputeLoops(ptSUXGluon_fld &Umu){

  ptSUX br, bl, mr, ml, tr, tl;
  ptSUX ld, lu, cd, cu, rd, ru;


  ptSUX w22,w21,w44,w42;
#ifdef __PARALLEL_OMP__
  ptSUX t_w22[NTHR],t_w21[NTHR],t_w44[NTHR],t_w42[NTHR];
#endif

  ptSUX tmp;

  int a,b,c,d,e,f,g,h,j,k;
  int i,l,m,n,o,p,q,r,s,t;

  w22.zero();
  w21.zero();
  w42.zero();
  w44.zero();

#ifdef __PARALLEL_OMP__
  for( int thr = 0; thr < NTHR; thr++)
    {
      t_w22[thr].zero();
      t_w21[thr].zero();
      t_w42[thr].zero();
      t_w44[thr].zero();    
    }
#endif

  for( int mu = 0; mu < dim-1; mu++)
    {
      for( int nu = mu+1; nu < dim; nu++)
	{

#ifdef __PARALLEL_OMP__
#pragma omp parallel for private(a,b,c,d,e,f,g,h,j,k, i,l,m,n,o,p,q,r,s,t,br, bl, mr, ml, tr, tl, ld, lu, cd, cu, rd, ru) num_threads(NTHR)
#endif
	  for( int x = 0; x < Umu.Z->Size; x++)
	    {
#ifdef __PARALLEL_OMP__
	      int thr = omp_get_thread_num();
#endif
	      a = Umu.Z->get(x, 1, mu);
	      b = Umu.Z->get(a, 1, mu);
	      c = Umu.Z->get(b, 1, nu);
	      e = Umu.Z->get(x, 1, nu);
	      d = Umu.Z->get(e, 1, mu);

  	      f = Umu.Z->get(c, 1, nu);
	      g = Umu.Z->get(d, 1, nu);
	      h = Umu.Z->get(e, 1, nu);


	      bl = Umu.W[x].U[mu]*Umu.W[a].U[mu];
	      ml = Umu.W[h].U[mu]*Umu.W[g].U[mu];
	      cd = Umu.W[b].U[nu]*Umu.W[c].U[nu];
	      ld = Umu.W[x].U[nu]*Umu.W[e].U[nu];

#ifdef  __PARALLEL_OMP__
	      t_w22[thr] += bl*cd*dag(ld*ml);
	      
	      t_w21[thr] += ( bl*Umu.W[b].U[nu]*
			      dag(Umu.W[x].U[nu]*Umu.W[e].U[mu]*Umu.W[d].U[mu]) +
			      Umu.W[x].U[mu]*Umu.W[a].U[nu]*Umu.W[d].U[nu]*
			      dag( ld*Umu.W[h].U[mu]) );	      
#else
	      w22 += bl*cd*dag(ld*ml);

	      w21 += ( bl*Umu.W[b].U[nu]*
		       dag(Umu.W[x].U[nu]*Umu.W[e].U[mu]*Umu.W[d].U[mu]) +
		       Umu.W[x].U[mu]*Umu.W[a].U[nu]*Umu.W[d].U[nu]*
		       dag( ld*Umu.W[h].U[mu]) );	      
#endif

              i = Umu.Z->get(h, 1, nu);
              j = Umu.Z->get(f, 1, nu);
              k = Umu.Z->get(j, 1, nu);
              n = Umu.Z->get(i, 1, nu);
              l = Umu.Z->get(n, 1, mu);
              m = Umu.Z->get(b, 1, mu);
              o = Umu.Z->get(m, 1, mu);
              p = Umu.Z->get(o, 1, nu);
              q = Umu.Z->get(p, 1, nu);
              r = Umu.Z->get(f, 1, mu);
              s = Umu.Z->get(q, 1, nu);
              t = Umu.Z->get(k, 1, mu);

              br = Umu.W[b].U[mu]*Umu.W[m].U[mu];
              mr = Umu.W[f].U[mu]*Umu.W[r].U[mu];
              tr = Umu.W[k].U[mu]*Umu.W[t].U[mu];
              tl = Umu.W[n].U[mu]*Umu.W[l].U[mu];
              lu = Umu.W[h].U[nu]*Umu.W[i].U[nu];             
              cu = Umu.W[f].U[nu]*Umu.W[j].U[nu];
              ru = Umu.W[q].U[nu]*Umu.W[s].U[nu];
              rd = Umu.W[o].U[nu]*Umu.W[p].U[nu];
              

#ifdef  __PARALLEL_OMP__
              t_w42[thr] += bl*(br*rd*dag(ml*mr) 
				+ cd*cu*dag(lu*tl))*dag(ld);
	      
              t_w44[thr] += bl*br*rd*ru*dag(ld*lu*tl*tr);	      
#else
              w42 += bl*(br*rd*dag(ml*mr) 
                         + cd*cu*dag(lu*tl))*dag(ld);

              w44 += bl*br*rd*ru*dag(ld*lu*tl*tr);
#endif
	    } // i
	} // nu
    } // mu



#ifdef  __PARALLEL_OMP__
  for( int thr = 0; thr < NTHR; thr++)
    {
      w21+=t_w21[thr];
      w22+=t_w22[thr];
      w42+=t_w42[thr];
      w44+=t_w44[thr];
    }
#endif


  w21.Tr(l21);
  w22.Tr(l22);
  w42.Tr(l42);
  w44.Tr(l44);

  const double loopNorm = 1./(3*2*3*(Umu.Z->Size));
  // loop_nxm << "Plaquette 1x2" << endl;
  for(int i1 = 0; i1 <= PTORD; i1+=2){  
    loop_nxm << l21[i1].re*.5*loopNorm << "\t";
  }					
  
  // loop_nxm << "Plaquette 2x2" << endl;
  for(int i1 = 0; i1 <= PTORD; i1+=2){  
    loop_nxm << l22[i1].re*loopNorm << "\t";
  }					
  
  // loop_nxm << "Plaquette 2x4" << endl;
  for(int i1 = 0; i1 <= PTORD; i1+=2){  
    loop_nxm << l42[i1].re*.5*loopNorm << "\t";
  }					
  
  // loop_nxm << "Plaquette 4x4" << endl;
  for(int i1 = 0; i1 <= PTORD; i1+=2){  
    loop_nxm << l44[i1].re*loopNorm << "\t";
  }					
  loop_nxm << std::endl;

}

#endif


// Static plaquette measurement
void plaquette_measure(ptSUXGluon_fld& Umu, act_params_t& act_pars){
#ifndef __PARALLEL_OMP__
  for(int i1 = 0; i1 < PTORD+1; i1++){
    w[i1] = 0.0;
  }					
  
  Cplx* app = new Cplx[PTORD+1];
  
  for(int i = 0; i < act_pars.iVol;i++){
    W1.zero();				
    link_c = Umu.Z->L[i][dim];
    for(int mu = 0; mu < dim-1; mu++){
      W.zero();				
      for(int nu = mu+1; nu < dim; nu++){
	//if(nu != mu ){
	  W += Umu.staple(i, mu, nu);
	//}
      }
      W1 = Umu.W[link_c].U[mu]*W;
      W1.Tr(app);
      for(int i1 = 0; i1 <= PTORD; i1++){
	w[i1] += app[i1];
      }
    }	
  }
  for(int i1 = 0; i1 <= PTORD; i1++){
#ifdef _runSU2
    w[i1] /= (act_pars.iVol*24);
#else
    w[i1] /= (act_pars.iVol*36);
#endif
  }

  W.zero();
  W1.zero();
  delete [] app;

#else

  for(int nt = 0; nt < NTHR; nt++){
    for(int i1 = 0; i1 < PTORD+1; i1++){
      ww[nt][i1] = 0.0;
    }				
  }
  
#pragma omp parallel private(tid,link_c) num_threads(NTHR)
  {
    int chunk = act_pars.iVol/NTHR;

    tid = omp_get_thread_num();    
    for(int site_x = chunk*tid; site_x < chunk*(tid+1);site_x++){
      Ww1[tid].zero();
      for(int mu = 0; mu < dim; mu++){
 	Ww[tid].zero();
 	for(int nu = 0; nu < dim; nu++){
 	  if(nu != mu ){
 	    Ww[tid] += Umu.staple(site_x, mu, nu);
 	  }
 	}
 	Ww1[tid] = Umu.W[site_x].U[mu]*Ww[tid];
 	Ww1[tid].Tr(ww1[tid]);
	
 	for(int i1 = 0; i1 <= PTORD; i1++){
 	  ww[tid][i1] += ww1[tid][i1];
 	}
	
      }// a given link at a given site
    } // end of loop on sites
    
    Ww[tid].zero();
    Ww1[tid].zero();
  } // end parallel
#pragma omp barrier

  for(int i1 = 0; i1 <= PTORD; i1++){
    w[i1] = 0;
    for(int nt = 0; nt < NTHR; nt++){
      //      std::cout << ww[nt][i1].re << "\t";
      w[i1] += ww[nt][i1];
    }
    //    std::cout << std::endl;
#ifdef _runSU2
    w[i1] /= (act_pars.iVol*48);
#else
    w[i1] /= (act_pars.iVol*72);
#endif
  }
  //  std::cout << std::endl;
#endif

}


