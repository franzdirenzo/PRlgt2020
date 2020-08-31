/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: MyMath.cc
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

#include<stdlib.h>
#include<string.h>
#include<vector>

#include"MyMath.h"

#define SWAP(a,b) { Cplx tmp = a; a = b; b = tmp; }

#define C0MAX_CONST 2./sqrt(27.)

double mod(Cplx& z) {
  return sqrt(z.re*z.re + z.im*z.im);;
}


Cplx sqrt(Cplx b) {
  double r = sqrt(b.re*b.re + b.im*b.im);
  return Cplx(sqrt(.5*(r+b.re)), b.im/sqrt(2.*(r+b.re)));
}

Cplx cbrt(Cplx b) {
  double theta = atan(b.im/b.re);
  double r = sqrt(b.re*b.re + b.im*b.im);
  return Cplx(cbrt(r)*cos(D3*theta), cbrt(r)*sin(D3*theta));
}


Cplx sin(Cplx b) {
  return Cplx(cosh(b.im)*sin(b.re), sinh(b.im)*cos(b.re));
}

Cplx cos(Cplx b) {
  return Cplx(cosh(b.im)*cos(b.re), -sinh(b.im)*sin(b.re));
}

Cplx tan(Cplx b) {
  double den;
  b.re *= 2.;
  b.im *= 2.;
  den   = 1./(cos(b.re) + cosh(b.im));
  return Cplx(sin(b.re)*den, sinh(b.re)*den);
}

Cplx exp(Cplx b){
  double rho = exp(b.re);
  return   Cplx(rho*cos(b.im), rho*sin(b.im));
}


Cplx atan(Cplx b) {
  Cplx z;
  double r2, i2;
  r2 = b.re*b.re;
  i2 = b.im*b.im;
  z.re = .5*atan2(2.*b.re,(1. - r2 - i2));
  z.im = .25*(log((1 + b.im)*(1 + b.im) - r2) -
	      log((1 - b.im)*(1 - b.im) + r2) ); 
  return(z);
}

//--------------------------------------------//


void SU3::prout() {
  for(int i = 0; i < 9; i++){
    if(i%3 == 0)    std::cout << std::endl;
    std::cout << " ";
    whr[i].prout();
  }
  std::cout << std::endl;
}


void SU3::dag() {
  double app;
  app = whr[1].re; whr[1].re =  whr[3].re; whr[3].re =  app;
  app = whr[1].im; whr[1].im = -whr[3].im; whr[3].im = -app;
  app = whr[2].re; whr[2].re =  whr[6].re; whr[6].re =  app;
  app = whr[2].im; whr[2].im = -whr[6].im; whr[6].im = -app;
  app = whr[5].re; whr[5].re =  whr[7].re; whr[7].re =  app;
  app = whr[5].im; whr[5].im = -whr[7].im; whr[7].im = -app;
  whr[0].im = -whr[0].im;
  whr[4].im = -whr[4].im;
  whr[8].im = -whr[8].im;
}

   


bool SU3::operator==(SU3 A){
  if(memcmp(whr, A.whr, 144)) return 1;
  else return 0;}

bool SU3::operator!=(SU3 A){
  if(memcmp(whr, A.whr, 144)) return 0;
  else return 1;}



SU3 SU3::operator~() const{  
  SU3 res;

  res.whr[0].re =  whr[0].re;
  res.whr[0].im = -whr[0].im;
  res.whr[1].re =  whr[1].re;
  res.whr[1].im = -whr[1].im;
  res.whr[2].re =  whr[2].re;
  res.whr[2].im = -whr[2].im;
  res.whr[3].re =  whr[3].re;
  res.whr[3].im = -whr[3].im;
  res.whr[4].re =  whr[4].re;
  res.whr[4].im = -whr[4].im;
  res.whr[5].re =  whr[5].re;
  res.whr[5].im = -whr[5].im;
  res.whr[6].re =  whr[6].re;
  res.whr[6].im = -whr[6].im;
  res.whr[7].re =  whr[7].re;
  res.whr[7].im = -whr[7].im;
  res.whr[8].re =  whr[8].re;
  res.whr[8].im = -whr[8].im;

  return res; 
}



void SU3::reU(){

  double den_A, den_B;
  Cplx num_A, tmp_C;
  
  den_A = mod2(whr[0]) + mod2(whr[1]) + mod2(whr[2]);
  num_A = ( whr[3] * ~(whr[0]) + whr[4] * ~(whr[1]) + 
	    whr[5] * ~(whr[2]) );
  tmp_C = num_A/den_A;

  whr[3] = whr[3] - tmp_C*whr[0];
  whr[4] = whr[4] - tmp_C*whr[1];
  whr[5] = whr[5] - tmp_C*whr[2];

  den_B = mod2(whr[3]) + mod2(whr[4]) + mod2(whr[5]);
  den_A = sqrt(den_A);
  den_B= sqrt(den_B);

  whr[0] /= den_A;
  whr[1] /= den_A;
  whr[2] /= den_A;
  whr[3] /= den_B;
  whr[4] /= den_B;
  whr[5] /= den_B;
  whr[6] = ~(whr[1] * whr[5] - whr[2] * whr[4]);
  whr[7] = ~(whr[2] * whr[3] - whr[0] * whr[5]);
  whr[8] = ~(whr[0] * whr[4] - whr[1] * whr[3]);
}


void SU3::reH(){
  SU3 a = *this, b = *this;
  a.dag();
  b -= a;
  (*this) = b*.5;
}







SU3 exp(SU3& Q){
  Cplx tra, det, h0, h1, h2, im(0,1);
  double c0, c1, c0max, u, w, th, q3, q1, q2;
  SU3 id, QQ;
  bool app = 0;
  QQ = Q*Q;

  id.whr[0] = 1;
  id.whr[4] = 1;
  id.whr[8] = 1;
  
  c0 = Q.Det().re;

  c1 = (.5*(QQ).Tr()).re;

  c0max = C0MAX_CONST*sqrt(c1*c1*c1);
  
  th = acos(c0/c0max);
  u  = sqrt(c1*D3)*cos(D3*th);
  w  = sqrt(c1)*sin(D3*th);

  //  printf("%lf    %lf     %lf\n",th, u , w);
  
  q1 = 2.*u;
  q2 = w-u;
  q3 = -u-w;
  //  printf("%lf    %lf     %lf\n",q1, q2 , q3);
  double x0, uv, u3v;
  if( fabs(w) > .5 )
    x0 = sin(w)/w;
  else
    x0 = 1-w*w*(1-.05*w*w*(1-w*w/42.))/6.;
  
  uv  = (u*u-w*w);

  Cplx e1 = exp(2.0*u*im);
  Cplx e2 = exp(-u*im);

  h0 = uv*e1+e2*(8.0*u*u*cos(w)+2.0*im*u*(3.0*u*u+w*w)*x0);
  h1 = 2*u*e1 - e2*((Cplx)2*u*cos(w)-im*(3*u*u-w*w)*x0);  
  h2 = e1-e2*((Cplx)cos(w)+3.0*im*u*x0);

  u3v = 1./(9.*u*u - w*w);
  if(app == 0){
    h0 *= u3v;    
    h1 *= u3v;
    h2 *= u3v;
  }
  else{
    h0 =   (~h0)*u3v;
    h1 = (~h1)*(-u3v);
    h2 =   (~h2)*u3v;
  }

  Q = h0*id + h1*Q + h2*QQ;

  return Q;
}



template<>
SU3 SU3rand(MyRand &Rand){
  SU3 P;
  register double g,gg;

  g = Rand.generate_gauss();
  P.whr[1].im = g;
  P.whr[3].im = g;

  g = Rand.generate_gauss();
  P.whr[1].re =  g;
  P.whr[3].re = -g;

  g = Rand.generate_gauss();
  P.whr[2].im =  g;
  P.whr[6].im =  g;  

  g = Rand.generate_gauss();
  P.whr[2].re =  g;
  P.whr[6].re = -g;  

  g = Rand.generate_gauss();
  P.whr[5].im =  g;
  P.whr[7].im =  g;  

  g = Rand.generate_gauss();
  P.whr[5].re =  g;
  P.whr[7].re = -g;  

  g = Rand.generate_gauss();
  gg = Rand.generate_gauss()*sD3;
  P.whr[0].im =     gg+g;
  P.whr[4].im =     gg-g;
  P.whr[8].im = -2.*gg;
    
  return P;
}


  template <>
  SU3 SU3rand(ranlxd::Rand& Rand){
    SU3 result;
    static const double twopi = std::atan(1.0) * 8.0;
    static const double soneo3 = std::sqrt( 1./ 3.);
    double g[8], r, t;
    // get flat distribution
    Rand.ranlxd(g, g+8);
    // make gaussian
    for (int i = 0; i < 8; i+=2){
      t = twopi*g[i];
      r = std::sqrt( -std::log((1. - g[i+1])) );
      g[i]   = r*std::cos(t);
      g[i+1] = r*std::sin(t);
    }
    g[7] *= soneo3;

    result.whr[0] = Cplx(0,g[7]+g[6]);
    result.whr[1] = Cplx(g[1],g[0]);
    result.whr[2] = Cplx(g[3],g[2]);
    result.whr[3] = Cplx(-g[1],g[0]);
    result.whr[4] = Cplx(0,g[7]-g[6]);
    result.whr[5] = Cplx(g[5],g[4]);
    result.whr[6] = Cplx(-g[3],g[2]);
    result.whr[7] = Cplx(-g[5],g[4]);
    result.whr[8] = Cplx(0,-g[7]*2);
    return result;
  }



template<>
SU2 SU2rand(MyRand &Rand){
  SU2 P;
  register double g,gg;

  g = Rand.generate_gauss();
  P.whr[1].im = g;
  P.whr[2].im = g;

  g = Rand.generate_gauss();
  P.whr[1].re =  g;
  P.whr[2].re = -g;

  g = Rand.generate_gauss();
  P.whr[0].im =  g;
  P.whr[3].im = -g;  

  return P;
}


    template <>
  SU2 SU2rand(ranlxd::Rand& Rand){
    SU2 result;
    static const double twopi = std::atan(1.0) * 8.0;
    static const double soneo3 = std::sqrt( 1./ 3.);
    double g[3], r, t;
    // get flat distribution
    Rand.ranlxd(g, g+3);
    // make gaussian
    for (int i = 0; i < 3; i+=2){
      t = twopi*g[i];
      r = std::sqrt( -std::log((1. - g[i+1])) );
      g[i]   = r*std::cos(t);
      g[i+1] = r*std::sin(t);
    }

    result.whr[0] = Cplx(0,g[2]);
    result.whr[1] = Cplx(g[1],g[0]);
    result.whr[2] = Cplx(-g[1],g[0]);
    result.whr[3] = Cplx(0,-g[2]);
    return result;
  }

