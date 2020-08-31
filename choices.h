/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: choices.h
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


// Please notice: dim set to 4 is de facto the only meaningful choice...
#define dim 4

#define allocORD 6

extern int PTORD;

// You can choose to run SU(2); otherwise it is SU(3)
// Generic SU(N) via templating on its way...
//#define _runSU2

//////////////////////////////////////////////////////
////////////////      GAUGE FIELD      ///////////////
//////////////////////////////////////////////////////

#define PERIODIC_BC
//#define TWISTED_BC

#ifdef TWISTED_BC
//#define _TWIST_X
#define _TWIST_Y
#define _TWIST_Z
//#define _TWIST_T
#endif

// ------------------------ //
// this option sets all the twist matrices to the identity
//#define IDENTITY_TWIST
// ------------------------ //

// ------------------------ //
// When this option is active, the gluonic evolution is performed "on the fly",
// that is every link is updated in succession.
// Otherwise the new links are computed freezing the current configuration.
#define UPDATE_ONTHEFLY
// ------------------------ //

// ------------------------ //
// this option allows for a dedicated (and slower) measure of the plaquette
// (only in some planes, avoiding boundaries, etc)
// everything has to be specified in nspt.cc
//#define SELECTED_PLAQ_MEASURE
// ------------------------ //

// ------------------------ //
// with PBC you can compute larger loops beyond
// the basic plaquette which is always measured;
// either only 2x2 or a bit more (as an example)
#ifdef PERIODIC_BC
#define __MANY_LOOP__
//#define __WILLOOP_2x2__
#endif
// ------------------------ //


#ifdef PERIODIC_BC
#define REMOVE_GLUON_ZEROMOM
#endif

#ifdef _runSU2
#define NC 2
#else
#define NC 3
#endif

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


#define _SQRT_BETA_EXPANSION_

#ifdef _runSU2

#define WIL 0
#define GAUGE_ACTION WIL

#else

#ifdef PERIODIC_BC
#define WIL 0
#define TLSYM 1
#define IWA 2
#define DBW2 3
#define GAUGE_ACTION WIL //Alternatives for SU(3) with PBC: WIL, TLSYM, IWA, DBW2
#else
#define WIL 0
#define GAUGE_ACTION WIL
#endif

#endif


//////////////////////////////////////////////////////
// Please notice: you typically want ntt set to 1
//////////////////////////////////////////////////////

#ifdef __PARALLEL_OMP__
#define ntt  1
#define ntx  1
#define nty  1
#define ntz  1
#define NTHR  (ntx*nty*ntz*ntt)

#else
#define NTHR 1
#endif

// you could in principle want a
// Z3 (or Z2) element as your vacuum
#ifdef __NONTRIVIAL_VACUUM__

// --------- SU3 --------

/* Identity */
/* #define VACUUM Cplx(1,0) */
/* #define VACUUM_CONJ Cplx(1,0) */

/* exp(2\pi/3 \im) */
#define VACUUM Cplx(-.5,.5*sqrt(3))
#define VACUUM_CONJ Cplx(-.5,-.5*sqrt(3))

/* exp(-2\pi/3 \im) */
/* #define VACUUM Cplx(-.5,.5*sqrt(3)) */
/* #define VACUUM_CONJ Cplx(-.5,-.5*sqrt(3)) */

// --------- SU2 -----------

/* Identity */
/* #define SU2_VACUUM Cplx(1,0) */
/* #define SU2_VACUUM_CONJ Cplx(1,0) */

/* exp(2\pi / 2 \im) */
#define SU2_VACUUM Cplx(-1, 0)
#define SU2_VACUUM_CONJ Cplx(-1,0)

#else
#define VACUUM 1
#define SU2_VACUUM 1
#endif

