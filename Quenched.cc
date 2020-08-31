/*************************************************************************************
PRlgt2020 NSPT library, www.github.com/franzdirenzo/PRlgt2020
Source file: Quenched.cc
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



// APE configurations legacy; not useful, but it was to sad to cut it 

#include<fstream>
#include<string>
#include<vector>

#include"nspt.h"

using namespace std;

#ifdef __MINIMAL_MPI__
#include "mpi.h"
#endif

Cplx* w       = new Cplx[PTORD+1];
Cplx* w1      = new Cplx[PTORD+1];
Cplx* w2      = new Cplx[PTORD+1];
double* norm  = new double[PTORD];
double* norm1 = new double[PTORD];
Cplx** trU;
ptSUX* U;
#ifndef UPDATE_ONTHEFLY
ptSUX* F;
#endif


nspt_params_t   nspt_pars;
act_params_t    act_pars;
thread_params_t thr_pars[NTHR];

int nInt = 0;
int config_counter=0;
char config_name[100];
char append[10];

#ifdef USE_RAND55
MyRand Rand[NTHR];
#else
extern std::vector<ranlxd::Rand> Rand;
#endif

std::ofstream plaqfile, normfile, logfile, trufile;

#ifdef __WILLOOP_2x2__
std::ofstream loop2x2;
#endif

#ifdef __MANY_LOOP__
std::ofstream loop_nxm;
#endif

#ifdef __TIMING__
PRlgtTime Time;
#endif


#ifdef __K_MOM_ANALYSIS__
extern ptSUX_fld* Wgauge;
extern fftw_plan *planFA;

extern fftw_plan *planGluon;

extern double *knorm[allocORD+1];
FILE *FPkn;
#endif

// Safely read parameters
// check matching string in front of a parameter value

static int get_val(FILE* fp, const char *str, const char* fmt,  void* val)
{
    char c[128];

    if(1!=fscanf(fp,"%s",c)){
      fprintf(stderr,"Error reading input file at %s\n",str);
      exit(1);
    }
    
    if(strcmp(str,c)!=0){
      fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
      exit(1);
    }

    if(1!=fscanf(fp,fmt,val)){
      fprintf(stderr,"Error reading input file at %s\n",str);
      fprintf(stderr,"Cannot read value format %s\n",fmt);
      exit(1);
    }
    
    return 0;
}




// Parameters from file or from command line
// relevant files are opened

int initialize(int argc, char** argv){

  long unsigned seed[NTHR];
  act_pars.sz       = new int[dim];
  nspt_pars.plaqn   = new char[100];
  nspt_pars.confn   = new char[100];
  nspt_pars.damon   = new char[100];
  nspt_pars.normn   = new char[100];
  nspt_pars.logn    = new char[100];
  nspt_pars.trun    = new char[100];
#ifdef __WILLOOP_2x2__
  nspt_pars.name2x2 = new char[100];
#endif
#ifdef __MANY_LOOP__
  nspt_pars.name_nxm = new char[100];
#endif

  trU = new Cplx* [dim];
  for( int i = 0; i < dim; i++)
    {
      trU[i] = new Cplx[PTORD];
    }


  // Read parameters from file
  // error returned in case of failure
  // Messages and errors in Italian! Not changed ...
  FILE *fp;
  if( (fp = fopen("Quench.cfg","r")) == NULL){
    std::cout << "Errore: impossibile leggere il file di configurazione." 
	      << std::endl;
    return 1;
  }

  fscanf( fp,"taglia %d %d %d %d\n", &act_pars.sz[0], &act_pars.sz[1], 
	  &act_pars.sz[2], &act_pars.sz[3] );

  get_val(fp, "SWEEP",       "%d" ,&(nspt_pars.Sweep) );
  get_val(fp, "BEAT",        "%d" ,&(nspt_pars.Beat)  );
  get_val(fp, "tau_g",       "%lf",&(act_pars.tau_g)  );
  get_val(fp, "alpha",       "%lf",&(act_pars.alpha)  );
  get_val(fp, "init_status", "%d" ,&(nspt_pars.Init)  );
  get_val(fp, "plaq_out",    "%s" ,nspt_pars.plaqn    );
  get_val(fp, "last_conf",   "%s" ,nspt_pars.confn    );
  get_val(fp, "PTORD",       "%d" ,&(PTORD)           );
  for( int tid = 0; tid < NTHR; tid++)
    {
      get_val(fp, "seed",        "%ld",&seed[tid]     );
    }

  fclose(fp);

  // File damocle
  strcpy(nspt_pars.damon,"damocle.dag");


  // What is read from file can be supersided
  // by command line options

  int opt;
  extern char *optarg;
  while ((opt = getopt(argc, argv, "l:s:b:d:r:g:a:t:p:c:n:h")) != -1) {
    switch (opt) {
    case 'l':
      act_pars.sz[0] = act_pars.sz[1] = act_pars.sz[2] = act_pars.sz[3] = atoi(optarg);
      break;
    case 's':
      nspt_pars.Sweep = atoi(optarg);
      break;
    case 'b':
      nspt_pars.Beat = atoi(optarg);
      break;
    case 'r':
#ifndef __PARALLEL_OMP__
      seed[0] = atol(optarg);
#else
      srand(atol(optarg));
      for(int thr = 0; thr < NTHR; thr++){
	seed[thr] = rand();
      }      
#endif
      break;
    case 'n':
      PTORD = atoi(optarg);
      break;
    case 'g':
      act_pars.tau_g = atof(optarg);
      break;
    case 'a':
      act_pars.alpha = -atof(optarg);
      break;
    case 't':
      nspt_pars.Init = atoi(optarg);
      break;
    case 'p':
      strcpy(nspt_pars.plaqn,optarg);
      break;
    case 'c':
      strcpy(nspt_pars.confn,optarg);
      break;
    case 'd':
      strcpy(nspt_pars.damon,optarg);
      break;
    case 'h':
      printf( "\n" );
      printf( "-l\t\ttaglia (reticolo t==x==y==z)\n" );
      printf( "-s\t\tnumero iterazioni (SWEEP)\n" );
      printf( "-b\t\titerazioni tra letture damocle (BEAT)\n" );
      printf( "-r\t\tseed generatore random\n" );
      printf( "-g\t\ttau gluoni\n" );
      printf( "-a\t\talpha\n" );
      printf( "-t\t\tstato iniziale: 1)da configurazione; 2)freddo\n" );
      printf( "-p\t\tnome file placchetta\n" );
      printf( "-c\t\tnome configurazione salvata\n" );
      printf( "-d\t\tnome file damocle\n" );
      printf( "-n\t\tordine pt\n" );
      printf( "\n" );
      exit(0);
      break;
    }
  }


  {
    FILE* df;
    if( ( df = fopen(nspt_pars.damon,"r") ) == NULL)
      {
	std::cout << "Errore, impossibile trovare file damocle." << std::endl;
	fclose(df);
	exit(-1);
      }
    fclose(df);
  }



  // Embarassing MPI in a farm environment...
#ifdef __MINIMAL_MPI__
  int rank;
  int rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  // Random generator initialization
#ifdef __MINIMAL_MPI__

  time(NULL);

#ifdef USE_RAND55
  for(int thr = 0; thr < NTHR; thr++){
    Rand[thr].init( seed[thr]+clock() );
  }
#else
  int volume = act_pars.sz[0]*act_pars.sz[1]*act_pars.sz[2]*act_pars.sz[3];
  std::vector<int> seeds(volume);
  srand(seed[0]);
  for (long i = 0; i < volume; ++i)
    seeds[i] = rand();
  Rand = std::vector<ranlxd::Rand>(seeds.begin(), seeds.end());
#endif
    
#else

#ifdef USE_RAND55
  for(int thr = 0; thr < NTHR; thr++){
    Rand[thr].init(seed[thr]);
  }
#else
  int volume = act_pars.sz[0]*act_pars.sz[1]*act_pars.sz[2]*act_pars.sz[3];
  std::vector<int> seeds(volume);
  srand(seed[0]);
  for (long i = 0; i < volume; ++i)
    seeds[i] = rand();
  Rand = std::vector<ranlxd::Rand>(seeds.begin(), seeds.end());
#endif
    
#endif


  // improved action parameters
#if GAUGE_ACTION == WIL
  act_pars.c1 = 0.0;
#elif GAUGE_ACTION == TLSYM
  act_pars.c1 = -1.0/12.0;  
#elif GAUGE_ACTION == IWA
  act_pars.c1 = -0.331;  
#elif GAUGE_ACTION == DBW2
  act_pars.c1 = -1.4088;  
#endif
  act_pars.c0 = 1. - 8. * act_pars.c1;



  // Build strings for file names 
  char *type_g = new char[8];
#ifdef _runSU2
#if GAUGE_ACTION == WIL
  strcpy( type_g, "SU2_wil" );
#elif GAUGE_ACTION == TLSYM
  strcpy( type_g, "SU2_tls" );
#elif GAUGE_ACTION == IWA
  strcpy( type_g, "SU2_iwa" );
#elif GAUGE_ACTION == DBW2
  strcpy( type_g, "SU2_dbw" );
#endif
#else
#if GAUGE_ACTION == WIL
    strcpy( type_g, "SU3_wil" );
#elif GAUGE_ACTION == TLSYM
    strcpy( type_g, "SU3_tls" );
#elif GAUGE_ACTION == IWA
    strcpy( type_g, "SU3_iwa" );
#elif GAUGE_ACTION == DBW2
    strcpy( type_g, "SU3_dbw" );
#endif
#endif

#ifdef __MINIMAL_MPI__
  
#ifdef __WILLOOP_2x2__
  sprintf(nspt_pars.name2x2,"%s2x2_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
#endif
#ifdef __MANY_LOOP__
  sprintf(nspt_pars.name_nxm,"%s_nxm_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
#endif

  sprintf(nspt_pars.plaqn,"%s_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
#ifndef APE_CONFIG
  sprintf(nspt_pars.confn,"%s%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.dat",
	  nspt_pars.confn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  allocORD,rank);
#endif
  sprintf(nspt_pars.logn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.log",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);
  sprintf(nspt_pars.normn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.nor",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD,rank);
  sprintf(nspt_pars.trun,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.r%d.trU",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD,rank);

#else
#ifdef __WILLOOP_2x2__
  sprintf(nspt_pars.name2x2,"%s2x2_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
#endif
#ifdef __MANY_LOOP__
  sprintf(nspt_pars.name_nxm,"%s_nxm_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
#endif

  sprintf(nspt_pars.plaqn,"%s_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.txt",
	  nspt_pars.plaqn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2], act_pars.sz[3], fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
#ifndef APE_CONFIG
  sprintf(nspt_pars.confn,"%s%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d",
	  nspt_pars.confn, type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  allocORD);
#endif
  sprintf(nspt_pars.logn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.log",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
  sprintf(nspt_pars.normn,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.nor",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);
  sprintf(nspt_pars.trun,"%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.trU",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g,
	  PTORD);

#endif

  // Files opening (for writing)
  if(!plaqfile.is_open() ) plaqfile.open(nspt_pars.plaqn, std::ios_base::app);
  plaqfile.precision(30);
  plaqfile << scientific;
  if(!normfile.is_open() ) normfile.open(nspt_pars.normn);
  if(!logfile.is_open()  ) logfile.open(nspt_pars.logn);
  if(!trufile.is_open()  ) trufile.open(nspt_pars.trun);
#ifdef __WILLOOP_2x2__
  if(!loop2x2.is_open() ) loop2x2.open(nspt_pars.name2x2, std::ios_base::app);
  loop2x2.precision(30);
  loop2x2 << scientific;
#endif
#ifdef __MANY_LOOP__
  if(!loop_nxm.is_open() ) loop_nxm.open(nspt_pars.name_nxm, std::ios_base::app);
  loop_nxm.precision(30);
  loop_nxm << scientific;
#endif

  // Initialization data on logfile
#ifdef _runSU2
    logfile << "Gauge group\t\tSU2"       << std::endl;
#else
    logfile << "Gauge group\t\tSU3"       << std::endl;
#endif
#ifdef UPDATE_ONTHEFLY
    logfile << "\tgauge field update \"on the fly\""       << std::endl;
#else
    logfile << "\tgauge field update not \"on the fly\""       << std::endl;
#endif
#ifdef TWISTED_BC
    logfile << "Gauge field boundary conditions\t\tTWISTED"       << std::endl;
#ifdef _TWIST_X
    logfile << "\t\t\t\t\t- X direction is twisted"       << std::endl;
#endif
#ifdef _TWIST_Y
    logfile << "\t\t\t\t\t- Y direction is twisted"       << std::endl;
#endif
#ifdef _TWIST_Z
    logfile << "\t\t\t\t\t- Z direction is twisted"       << std::endl;
#endif
#ifdef _TWIST_T
    logfile << "\t\t\t\t\t- T direction is twisted"       << std::endl;
#endif
#endif
#ifdef PERIODIC_BC
    logfile << "Gauge field boundary conditions\t\tPERIODIC"       << std::endl;
#endif
  logfile << "Taglia\t\t"       << act_pars.sz[0]  << " " << act_pars.sz[1] << " " << act_pars.sz[2] << " " << act_pars.sz[3] << std::endl
	  << "Sweep\t\t"        << nspt_pars.Sweep   << std::endl
	  << "Beat\t\t"         << nspt_pars.Beat    << std::endl
	  << "Tau_g\t\t"        << act_pars.tau_g    << std::endl
	  << "Alpha\t\t"        << act_pars.alpha    << std::endl
	  << "Init status\t"    << nspt_pars.Init    << std::endl
	  << "Ordine\t\t"       << PTORD             << std::endl
	  <<                                            std::endl
#ifdef __PARALLEL_OMP__
	  << "Numero Threads\t" << NTHR              << std::endl
#endif
	  << "Placchetta\t"     << nspt_pars.plaqn   << std::endl
#ifdef __WILLOOP_2x2__
	  << "Altri loops\t"    << nspt_pars.name2x2 << std::endl
#endif
#ifdef __MANY_LOOP__
	  << "Altri loops\t"    << nspt_pars.name_nxm << std::endl
#endif
	  << "Configurazione\t" << nspt_pars.confn   << std::endl
	  << "Damocle\t\t"      << nspt_pars.damon   << std::endl
	  << "Norma\t\t"        << nspt_pars.normn   << std::endl
	  << "Tr(U)\t\t"        << nspt_pars.trun    << std::endl << std::endl;

#ifdef __PARALLEL_OMP__

  cout << "\tThreads partitionning:" << endl;


#if GAUGE_ACTION == WIL

  if ( (act_pars.sz[0]/ntt) < 2 ) { cout << "Errore: troppi threads nella direzione 0 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[1]/ntx) < 2 ) { cout << "Errore: troppi threads nella direzione 1 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[2]/nty) < 2 ) { cout << "Errore: troppi threads nella direzione 2 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[3]/ntz) < 2 ) { cout << "Errore: troppi threads nella direzione 3 rispetto alla taglia di reticolo." << endl; exit(-1); }

#else
	
  if ( (act_pars.sz[0]/ntt) < 3 ) { cout << "Errore: troppi threads nella direzione 0 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[1]/ntx) < 3 ) { cout << "Errore: troppi threads nella direzione 1 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[2]/nty) < 3 ) { cout << "Errore: troppi threads nella direzione 2 rispetto alla taglia di reticolo." << endl; exit(-1); }
  if ( (act_pars.sz[3]/ntz) < 3 ) { cout << "Errore: troppi threads nella direzione 3 rispetto alla taglia di reticolo." << endl; exit(-1); }

#endif
	
  for( int tid = 0; tid < NTHR; tid++){
    thr_pars[tid].xi = new int[dim];
    thr_pars[tid].xf = new int[dim];

    thr_pars[tid].xi[0] = ( act_pars.sz[0] / ntt ) * (  tid / (ntx*nty*ntz)       );
    thr_pars[tid].xi[1] = ( act_pars.sz[1] / ntx ) * ( (tid / (nty*ntz)   ) % ntx );
    thr_pars[tid].xi[2] = ( act_pars.sz[2] / nty ) * ( (tid / (ntz)       ) % nty );
    thr_pars[tid].xi[3] = ( act_pars.sz[3] / ntz ) * (  tid                 % ntz );
    
    thr_pars[tid].xf[0] = thr_pars[tid].xi[0] + act_pars.sz[0]/ntt ;
    thr_pars[tid].xf[1] = thr_pars[tid].xi[1] + act_pars.sz[1]/ntx ;
    thr_pars[tid].xf[2] = thr_pars[tid].xi[2] + act_pars.sz[2]/nty ;
    thr_pars[tid].xf[3] = thr_pars[tid].xi[3] + act_pars.sz[3]/ntz ;

    cout << "tid = " << tid << "\t"
	 << "\tseed = " << seed[tid] 
	 << std::endl
	 << thr_pars[tid].xi[0] << "\t"
	 << thr_pars[tid].xi[1] << "\t"
	 << thr_pars[tid].xi[2] << "\t"
	 << thr_pars[tid].xi[3] << "\n"
	 << thr_pars[tid].xf[0] << "\t"
	 << thr_pars[tid].xf[1] << "\t"
	 << thr_pars[tid].xf[2] << "\t"
	 << thr_pars[tid].xf[3] << "\n\n";
      
  }

#endif


#ifdef __K_MOM_ANALYSIS__
  char *kname = new char[100];
#ifndef __LANDAU_GAUGE_FIXING__
  sprintf(kname,"knorm_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.dat",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD);
#else
  sprintf(kname,"knorm_acc_%s_nf0_%dx%dx%dx%d_a%2.2f_tg%3.3f_o%d.dat",
	  type_g, act_pars.sz[0], act_pars.sz[1], 
	  act_pars.sz[2],act_pars.sz[3],fabs(act_pars.alpha), act_pars.tau_g, 
	  PTORD);
#endif
  FPkn = fopen(kname,"wb");
  delete [] kname;
#endif

  return 0;
}




#ifdef UPDATE_ONTHEFLY
int QuenchedAllocate(ptSUXGluon_fld& Umu){
#else
int QuenchedAllocate(ptSUXGluon_fld& Umu,ptSUXGluon_fld& Fmu){
#endif

    // Compute volume and 1/V normalization
  act_pars.iVol = ( act_pars.sz[0] * act_pars.sz[1] *
		    act_pars.sz[2] * act_pars.sz[3] );
  act_pars.rVol = 1.0/(double)act_pars.iVol;

  // Define sqrt(tau_g) and tau_g/(2Nc)
  act_pars.stau = sqrt(act_pars.tau_g); 
#ifdef _SQRT_BETA_EXPANSION_
  act_pars.tau_g /= -(2*NC);
#elif defined _G_EXPANSION_
  act_pars.tau_g *= -1;
#endif

  // Gluonic field initialization, either cold or hot (from file)
  // according to relevant flag in nspt_pars.Init. 
  // In case of hot start, plaquette value read as well
  // to cross check

  if(!nspt_pars.Init){

    for(int i = 0; i < Umu.Z->TSize; i++){
      for(int mu = 0; mu < dim; mu++){
	Umu.W[i].U[mu].id();
      }
    }

  }
  else{
#ifndef APE_CONFIG    
    if( Umu.load(nspt_pars.confn, w1) ){
      std::cout << "Errore nella lettura della configurazione." << std::endl;
      return 1;
    }
#else
    Umu.load_ape();
#endif
    plaquette_measure(Umu, act_pars);
#ifndef APE_CONFIG
    if( plaquette_check(w, w1) ){
      std::cout << "La misura della placchetta non corrisponde." << std::endl;
      for (int i1=0; i1 <= PTORD; i1++){
	printf("%e %e\t",w[i1].re, w[i1].im );
	printf("%e %e\n",w1[i1].re,w1[i1].im);
      }
      //      return 1;
    }
#endif
    std::cout << "\nMisura della placchetta:" << std::endl;
    for (int i1=0; i1 <= PTORD; i1++){
      printf("%e %e\n",w[i1].re,w[i1].im);
    }

  }

  
  U = Umu.handle();
#ifndef UPDATE_ONTHEFLY
  F = Fmu.handle();
#endif

  // Starting time on logfile
  time_t tempo;
  time(&tempo);
  logfile << "Inizio run\t\t" << ctime(&tempo) << std::endl;

#ifdef __K_MOM_ANALYSIS__
  for(int i1 = 0; i1 <= allocORD;i1++)
    knorm[i1] = new double[Umu.Z->Size];
#endif

  return 0;
}


// update parameters
// damocle.dag file always there to update parameters
// or even to kill the job (that's the reason for the name:
// google Damocle's dag...)
int AggiornaParametri(ptSUXGluon_fld& Umu) {

#ifdef __WILLOOP_2x2__
  WL2x2(Umu);
#endif

#ifdef __MANY_LOOP__
  ComputeLoops(Umu);
#endif

  // Scrive log
  time_t tempo;
  time(&tempo);
  nInt += nspt_pars.Beat;
  logfile << "______________________________________"  << std::endl 
	  << std::endl;
  logfile << "Damocle Update\t\t" << ctime(&tempo)  << std::endl;
  logfile << "# Iterazioni\t\t"   << nInt           << std::endl;

  // Read damocle file and update parameters
  FILE *df;
  df = fopen(nspt_pars.damon,"r");

  get_val(df, "SWEEP",       "%d" ,&(nspt_pars.Sweep) );
  get_val(df, "BEAT",        "%d" ,&(nspt_pars.Beat)  );
  get_val(df, "tau_g",       "%lf",&(act_pars.tau_g)  );
  get_val(df, "alpha",       "%lf",&(act_pars.alpha)  );
  get_val(df, "save_config", "%d" ,&(nspt_pars.Save)  );
  get_val(df, "kill_run",    "%d" ,&(nspt_pars.Kill)  );


  logfile << "Tau_g\t\t\t"        << act_pars.tau_g << std::endl;
  logfile << std::endl;
  logfile.flush();

  act_pars.stau = sqrt(act_pars.tau_g);
#ifdef _SQRT_BETA_EXPANSION_
  act_pars.tau_g /= -(2*NC);
#elif defined _G_EXPANSION_
  act_pars.tau_g *= -1;
#endif

  // In case of checkpoint (save_config = 1)
  if(nspt_pars.Save == 1) {
    plaquette_measure(Umu, act_pars);
      
      sprintf(append, "_%d.dat", ++config_counter);
      strcpy(config_name,nspt_pars.confn);
      strcat(config_name,append);
      Umu.save(config_name, w);
  }

  // If you want to kill the run (kill_run = 1)
  if( nspt_pars.Kill == 1) { 
    nspt_pars.Sweep = 0;
  }

  return 0;
}


int NsptFinalize(ptSUXGluon_fld& Umu, int t){

#ifdef __K_MOM_ANALYSIS__
  fclose(FPkn);
#endif

  time_t tempo; 
  time(&tempo);
  logfile << "__________________________________"  << std::endl
	  << "Termine simulazione\t\t"             << ctime(&tempo) << std::endl
	  << "# Iterazioni\t\t"                    <<  t << std::endl;

  plaquette_measure(Umu, act_pars);

  printf("Misura della placchetta:\n");
  for (int i1=0; i1 <= PTORD; i1++){
    printf("%e\t%e\n", w[i1].re, w[i1].im);
  }


#ifdef __WILLOOP_2x2__
  WL2x2(Umu);
  loop2x2.close();
#endif

#ifdef __MANY_LOOP__
  ComputeLoops(Umu);
  loop_nxm.close();
#endif
  

  Umu.save(nspt_pars.confn,w);
  delete [] w;
  delete [] w1;
  delete [] w2;
  delete [] norm;
  delete [] norm1;

  return 0;
}



int main(int argc, char** argv){

#ifdef PERIODIC_BC
#ifdef TWISTED_BC
    std::cout<<"ERROR: BOTH TWISTED AND PERIODIC BC DEFINED"<<std::endl;
    exit(0);
#endif
#endif
    
#ifndef PERIODIC_BC
#ifndef TWISTED_BC
    std::cout<<"ERROR: NO BC DEFINED"<<std::endl;
    exit(0);
#endif
#endif
    
  if( initialize(argc, argv) ){ exit(1); }


  // Lattice allocation and momenta generation
  latt LL(act_pars.sz);  
  LL.p_init();


  // (pt) gluon field allocated
  ptSUXGluon_fld Umu(&LL);
    
#ifndef UPDATE_ONTHEFLY
// Fmu stores a backup configuration
  ptSUXGluon_fld Fmu(&LL);
#endif
    
#ifdef __K_MOM_ANALYSIS__

#ifdef  __PARALLEL_OMP__
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHR);
#endif

  planGluon = new fftw_plan[2];
    
    planGluon[0] = fftw_plan_many_dft(dim, LL.Sz, dim*(NC*NC*allocORD+1),
                                      (fftw_complex *) Umu.W,
                                      NULL, dim*(NC*NC*allocORD+1), 1,
                                      (fftw_complex *) Umu.W,
                                      NULL, dim*(NC*NC*allocORD+1), 1,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    planGluon[1] = fftw_plan_many_dft(dim, LL.Sz, dim*(NC*NC*allocORD+1),
                                      (fftw_complex *) Umu.W,
                                      NULL, dim*(NC*NC*allocORD+1), 1,
                                      (fftw_complex *) Umu.W, 
                                      NULL, dim*(NC*NC*allocORD+1), 1,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);

#endif

#ifdef __LANDAU_GAUGE_FIXING__

#ifdef  __PARALLEL_OMP__
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHR);
#endif

  ptSUX_fld WFA(&LL);
  Wgauge = &WFA;

  planFA = new fftw_plan[2];
  planFA[0] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD+1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     FFTW_FORWARD, FFTW_ESTIMATE);
  planFA[1] = fftw_plan_many_dft(dim, LL.Sz, NC*NC*PTORD+1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     (fftw_complex *) WFA.W, 
  				     NULL, NC*NC*allocORD+1, 1,
  				     FFTW_BACKWARD, FFTW_ESTIMATE);
#endif

#ifdef UPDATE_ONTHEFLY
  if( QuenchedAllocate(Umu) ) { exit(1); }
#else
  if( QuenchedAllocate(Umu,Fmu) ) { exit(1); }
#endif
    
  int t1 = 0;

#ifdef __TIMING__
  std::string name("time_");
  name += nspt_pars.logn;
  ofstream  of_timing(name.c_str());
  of_timing << "Gauge\tFermion\tGaugeFix\tZeroMom" << endl << endl; 
#endif
    
  while( t1 < nspt_pars.Sweep) {
    t1 += nspt_pars.Beat;
#ifdef UPDATE_ONTHEFLY
    NsptEvolve(Umu);
#else
    NsptEvolve(Umu,Fmu);
#endif

    // Read new parameters from damocle
    if( AggiornaParametri(Umu) ){
      return 0;
    }
  }
    
  NsptFinalize(Umu, t1);
    
  return 0;
}
