/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: MyTime.h
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

#ifndef _MYTIME_H
#define _MYTIME_H
#include<iostream>
#include<fstream>

#include<map>

#ifdef HAVE_STDCXX_0X
#include<chrono>
#else
#include<ctime>
#endif

#include "input.h"

extern nspt_params_t nspt_pars;

namespace mytime {

  struct Timer {
    Timer () : t(0.0) { };
    void tic() { 
#ifdef HAVE_STDCXX_0X
      t1 = std::chrono::high_resolution_clock::now();
#else
      t1 = clock();
#endif
    }
    void toc() { 
#ifdef HAVE_STDCXX_0X
      t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = 
	std::chrono::duration_cast<std::chrono::seconds> (t2 - t1);
#else
      t2 = clock();
      clock_t elapsed = t2-t1;
#endif
      t += elapsed;
    }
    double out() { 
#ifdef HAVE_STDCXX_0X
      return t.count(); 
#else
      return static_cast<double>(t)/CLOCKS_PER_SEC;
#endif
    }
    void reset() { t *= 0.0; }
  private:
#ifdef HAVE_STDCXX_0X
    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<double> t;
#else
    clock_t t1, t2, t;
#endif
  };


}

typedef std::map<std::string, mytime::Timer> PRlgtTime;

#endif //MYTIME_H
