/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: twistmatrices.h
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

#include"MyMath.h"

#ifdef _runSU2

#ifdef _TWIST_X
#ifdef IDENTITY_TWIST
// identity
Cplx omega_x_list0(1,0);
Cplx omega_x_list1(0,0);
Cplx omega_x_list2(0,0);
Cplx omega_x_list3(1,0);
#else
// twist matrix...
Cplx omega_x_list0(1,0);
Cplx omega_x_list1(0,0);
Cplx omega_x_list2(0,0);
Cplx omega_x_list3(1,0);
#endif

Cplx omega_x_list[4]={omega_x_list0,omega_x_list1,omega_x_list2,omega_x_list3};
SU2 omega_x(omega_x_list);
Cplx omega_x_listdag[4]={~omega_x_list0,~omega_x_list2,~omega_x_list1,~omega_x_list3};
SU2 omega_xdag(omega_x_listdag);
#endif

#ifdef _TWIST_Y
#ifdef IDENTITY_TWIST
// identity
Cplx omega_y_list0(1,0);
Cplx omega_y_list1(0,0);
Cplx omega_y_list2(0,0);
Cplx omega_y_list3(1,0);
#else
Cplx omega_y_list0(0,-1);
Cplx omega_y_list1(0,0);
Cplx omega_y_list2(0,0);
Cplx omega_y_list3(0,1);
// maximal twist

#endif

Cplx omega_y_list[4]={omega_y_list0,omega_y_list1,omega_y_list2,omega_y_list3};
SU2 omega_y(omega_y_list);
Cplx omega_y_listdag[4]={~omega_y_list0,~omega_y_list2,~omega_y_list1,~omega_y_list3};
SU2 omega_ydag(omega_y_listdag);
#endif

#ifdef _TWIST_Z
#ifdef IDENTITY_TWIST
// identity
Cplx omega_z_list0(1,0);
Cplx omega_z_list1(0,0);
Cplx omega_z_list2(0,0);
Cplx omega_z_list3(1,0);
#else
// maximal twist
Cplx omega_z_list0(0,0);
Cplx omega_z_list1(1,0);
Cplx omega_z_list2(-1,0);
Cplx omega_z_list3(0,0);
#endif

Cplx omega_z_list[4]={omega_z_list0,omega_z_list1,omega_z_list2,omega_z_list3};
SU2 omega_z(omega_z_list);
Cplx omega_z_listdag[4]={~omega_z_list0,~omega_z_list2,~omega_z_list1,~omega_z_list3};
SU2 omega_zdag(omega_z_listdag);
#endif

#ifdef _TWIST_T
#ifdef IDENTITY_TWIST
// identity
Cplx omega_t_list0(1,0);
Cplx omega_t_list1(0,0);
Cplx omega_t_list2(0,0);
Cplx omega_t_list3(1,0);
#else
// maximal twist
Cplx omega_t_list0(0,0);
Cplx omega_t_list1(0,-1);
Cplx omega_t_list2(0,-1);
Cplx omega_t_list3(0,0);
#endif

Cplx omega_t_list[4]={omega_t_list0,omega_t_list1,omega_t_list2,omega_t_list3};
SU2 omega_t(omega_t_list);
Cplx omega_t_listdag[4]={~omega_t_list0,~omega_t_list2,~omega_t_list1,~omega_t_list3};
SU2 omega_tdag(omega_t_listdag);
#endif

#else

#ifdef _TWIST_X
#ifdef IDENTITY_TWIST
// identity
Cplx omega_x_list0(1,0);
Cplx omega_x_list1(0,0);
Cplx omega_x_list2(0,0);
Cplx omega_x_list3(0,0);
Cplx omega_x_list4(1,0);
Cplx omega_x_list5(0,0);
Cplx omega_x_list6(0,0);
Cplx omega_x_list7(0,0);
Cplx omega_x_list8(1,0);
#else
// mio twist
//Cplx omega_x_list0(0,0);
//Cplx omega_x_list1(1,0);
//Cplx omega_x_list2(0,0);
//Cplx omega_x_list3(0,0);
//Cplx omega_x_list4(0,0);
//Cplx omega_x_list5(1,0);
//Cplx omega_x_list6(1,0);
//Cplx omega_x_list7(0,0);
//Cplx omega_x_list8(0,0);
//Cplx omega_x_list0(1,0);
//Cplx omega_x_list1(0,0);
//Cplx omega_x_list2(0,0);
//Cplx omega_x_list3(0,0);
//Cplx omega_x_list4(1,0);
//Cplx omega_x_list5(0,0);
//Cplx omega_x_list6(0,0);
//Cplx omega_x_list7(0,0);
//Cplx omega_x_list8(1,0);
#endif

Cplx omega_x_list[9]={omega_x_list0,omega_x_list1,omega_x_list2,omega_x_list3,omega_x_list4,omega_x_list5,omega_x_list6,omega_x_list7,omega_x_list8};
SU3 omega_x(omega_x_list);

Cplx omega_x_listdag[9]={~omega_x_list0,~omega_x_list3,~omega_x_list6,~omega_x_list1,~omega_x_list4,~omega_x_list7,~omega_x_list2,~omega_x_list5,~omega_x_list8};
SU3 omega_xdag(omega_x_listdag);
#endif

#ifdef _TWIST_Y
#ifdef IDENTITY_TWIST
// identity
Cplx omega_y_list0(1,0);
Cplx omega_y_list1(0,0);
Cplx omega_y_list2(0,0);
Cplx omega_y_list3(0,0);
Cplx omega_y_list4(1,0);
Cplx omega_y_list5(0,0);
Cplx omega_y_list6(0,0);
Cplx omega_y_list7(0,0);
Cplx omega_y_list8(1,0);
#else
Cplx omega_y_list0(cos(-2.*M_PI/3.),sin(-2.*M_PI/3.));
Cplx omega_y_list1(0,0);
Cplx omega_y_list2(0,0);
Cplx omega_y_list3(0,0);
Cplx omega_y_list4(1,0);
Cplx omega_y_list5(0,0);
Cplx omega_y_list6(0,0);
Cplx omega_y_list7(0,0);
Cplx omega_y_list8(cos(2.*M_PI/3.),sin(2.*M_PI/3.));
#endif

Cplx omega_y_list[9]={omega_y_list0,omega_y_list1,omega_y_list2,omega_y_list3,omega_y_list4,omega_y_list5,omega_y_list6,omega_y_list7,omega_y_list8};
SU3 omega_y(omega_y_list);

Cplx omega_y_listdag[9]={~omega_y_list0,~omega_y_list3,~omega_y_list6,~omega_y_list1,~omega_y_list4,~omega_y_list7,~omega_y_list2,~omega_y_list5,~omega_y_list8};
SU3 omega_ydag(omega_y_listdag);
#endif

#ifdef _TWIST_Z
#ifdef IDENTITY_TWIST
// identity
Cplx omega_z_list0(1,0);
Cplx omega_z_list1(0,0);
Cplx omega_z_list2(0,0);
Cplx omega_z_list3(0,0);
Cplx omega_z_list4(1,0);
Cplx omega_z_list5(0,0);
Cplx omega_z_list6(0,0);
Cplx omega_z_list7(0,0);
Cplx omega_z_list8(1,0);
#else
Cplx omega_z_list0(0,0);
Cplx omega_z_list1(1,0);
Cplx omega_z_list2(0,0);
Cplx omega_z_list3(0,0);
Cplx omega_z_list4(0,0);
Cplx omega_z_list5(1,0);
Cplx omega_z_list6(1,0);
Cplx omega_z_list7(0,0);
Cplx omega_z_list8(0,0);
#endif

Cplx omega_z_list[9]={omega_z_list0,omega_z_list1,omega_z_list2,omega_z_list3,omega_z_list4,omega_z_list5,omega_z_list6,omega_z_list7,omega_z_list8};
SU3 omega_z(omega_z_list);
Cplx omega_z_listdag[9]={~omega_z_list0,~omega_z_list3,~omega_z_list6,~omega_z_list1,~omega_z_list4,~omega_z_list7,~omega_z_list2,~omega_z_list5,~omega_z_list8};
SU3 omega_zdag(omega_z_listdag);
#endif

#ifdef _TWIST_T
#ifdef IDENTITY_TWIST
// identity
Cplx omega_t_list0(1,0);
Cplx omega_t_list1(0,0);
Cplx omega_t_list2(0,0);
Cplx omega_t_list3(0,0);
Cplx omega_t_list4(1,0);
Cplx omega_t_list5(0,0);
Cplx omega_t_list6(0,0);
Cplx omega_t_list7(0,0);
Cplx omega_t_list8(1,0);
#else
Cplx omega_t_list0(0,0);
Cplx omega_t_list1(cos(-2.*M_PI/3.),sin(-2.*M_PI/3.));
Cplx omega_t_list2(0,0);
Cplx omega_t_list3(0,0);
Cplx omega_t_list4(0,0);
Cplx omega_t_list5(1,0);
Cplx omega_t_list6(cos(2.*M_PI/3.),sin(2.*M_PI/3.));
Cplx omega_t_list7(0,0);
Cplx omega_t_list8(0,0);
#endif

Cplx omega_t_list[9]={omega_t_list0,omega_t_list1,omega_t_list2,omega_t_list3,omega_t_list4,omega_t_list5,omega_t_list6,omega_t_list7,omega_t_list8};
SU3 omega_t(omega_t_list);
Cplx omega_t_listdag[9]={~omega_t_list0,~omega_t_list3,~omega_t_list6,~omega_t_list1,~omega_t_list4,~omega_t_list7,~omega_t_list2,~omega_t_list5,~omega_t_list8};
SU3 omega_tdag(omega_t_listdag);
#endif


#endif
