/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: nspt.h
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

#include "input.h"
#include "MyTime.h"

#include <vector>

#ifdef _runSU2
typedef SU2_fld SUX_fld;
typedef SU2Gluon_fld SUXGluon_fld;
typedef ptSU2Gluon_fld ptSUXGluon_fld;
typedef ptSU2_fld ptSUX_fld;
typedef SU2Gluon SUXGluon;
typedef ptSU2Gluon ptSUXGluon;
typedef ptSU2 ptSUX;
typedef SU2 SUX;
typedef ptCSU2Vector ptCSUXVector;
typedef CSU2Vector CSUXVector;
#else
typedef SU3_fld SUX_fld;
typedef Gluon_fld SUXGluon_fld;
typedef ptGluon_fld ptSUXGluon_fld;
typedef ptSU3_fld ptSUX_fld;
typedef Gluon SUXGluon;
typedef ptGluon ptSUXGluon;
typedef ptSU3 ptSUX;
typedef SU3 SUX;
typedef ptCVector ptCSUXVector;
typedef CVector CSUXVector;
#endif

#ifdef UPDATE_ONTHEFLY
void gauge_wilson(ptSUXGluon_fld&);
#else
void gauge_wilson(ptSUXGluon_fld&,ptSUXGluon_fld&);
#endif

#ifdef TWISTED_BC
ptSUX twist_x(const ptSUX &);
ptSUX invtwist_x(const ptSUX &);
ptSUX twist_y(const ptSUX &);
ptSUX invtwist_y(const ptSUX &);
ptSUX twist_z(const ptSUX &);
ptSUX invtwist_z(const ptSUX &);
ptSUX twist_t(const ptSUX &);
ptSUX invtwist_t(const ptSUX &);
void twist_boundary(ptSUXGluon_fld&);
#endif


#ifdef REMOVE_GLUON_ZEROMOM
void zero_modes_subtraction(ptSUXGluon_fld&);
#endif

void stochastic_gauge_fixing(ptSUXGluon_fld&);
void FAstochastic_gauge_fixing(ptSUXGluon_fld&);

// Quenched evolution
#ifdef UPDATE_ONTHEFLY
void NsptEvolve(ptSUXGluon_fld&);
#else
void NsptEvolve(ptSUXGluon_fld&, ptSUXGluon_fld&);
#endif


void plaquette_measure(ptSUXGluon_fld&, act_params_t&);

int plaquette_check(Cplx*, Cplx*);

void WL2x2(ptSUXGluon_fld&);
void ComputeLoops(ptSUXGluon_fld&);


