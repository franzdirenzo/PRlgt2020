/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: input.h
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

#ifndef _INPUT_H_
#define _INPUT_H_

#include<iostream>
#include<fstream>

#include"QCDenvNODEpt_test.h"

#ifdef MACOS
#include <mach/mach_time.h>
#endif


typedef struct{
  int Sweep;
  int Beat;
  int Kill;
  int Save;
  int Init;
  char *plaqn;
  char *confn;
  char *damon;
  char *normn;
  char *trun;
  char *logn;
  char *name2x2;
  char *name_nxm;
} nspt_params_t;


typedef struct{
  int iVol;
  int *sz;
  double rVol;
  double tau_f;
  double tau_g;
  double stau;
  double alpha;
  double c0;
  double c1;
} act_params_t;

typedef struct{
  int *xi;
  int *xf;
} thread_params_t;

#endif
