#include <limits.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
using namespace std;

#include"FiEstAS.h"

//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  clock_t t0 = clock();

  yINFO(("\n--------------------------------------------------------------------------------"));
  yINFO(("\n Compute densities from the ASCII file <data> \n"));
  yINFO(("\n    Usage  : FiEstAS [options] <data> \n"));
  yINFO(("\n    Options: -d=<D> Consider only first <D> columns (default=0, for all columns)"));
  yINFO(("\n             -kernel=<TopHat(default)|TSC|Epanechnikov> "));
  yINFO(("\n             -m0=<M0> (default=2.0) "));
  yINFO(("\n             -balloon=<true(default)|false> return balloon estimate "));
  yINFO(("\n             -at=<locations> ASCII file with the points where the density "));
  yINFO(("\n                             is to be evaluated (default=<data>) \n"));
  yINFO(("\n    Output: FiEstAS.txt \n"));
  yINFO(("\n Yago Ascasibar (AIP, Summer 2007)"));
  yINFO(("\n--------------------------------------------------------------------------------\n\n"));

  ERROR(argc<2,("Wrong syntax"));

  DIM D = option<int>("-d=",0, argc-1,argv);
  yINFO((" D = %d \n", D));
  
  string kernel = option<string>("-kernel=","TopHat", argc-1,argv);
  yKernel *K = NULL;
  if(kernel=="TopHat") K = new TopHat();
  else if(kernel=="TSC") K = new TSC();
  else if(kernel=="Epanechnikov") K = new Epanechnikov();
  ERROR( K==NULL , ("Wrong kernel '%s'",kernel.c_str()) );
  yINFO((" %s Kernel \n", kernel.c_str()));
  
  double M0 = option<double>("-m0=",2., argc-1,argv);
  ERROR( M0<=0. , ("M0=%g (<=0.)",M0) );
  yINFO((" M0 = %g \n", M0));
  
  string bal = option<string>("-balloon=","true", argc-1,argv);
  bool balloon;
  if( bal=="true") balloon = true;
  else
  {
    ERROR( bal!="false", ("Wrong value '%s' for option -balloon",bal.c_str()) );
    balloon = false;
  }
  if(balloon) yINFO((" balloon = true \n"));
  else yINFO((" balloon = false \n"));
  
  string data(argv[argc-1]);
  yINFO((" data = '%s' \n",data.c_str()));
  
  string locations = option<string>("-at=",string(argv[argc-1]), argc-1,argv);
  yINFO((" locations = '%s' \n",locations.c_str()));

  // ---------------------------------------------------------------

  yTree tree(D);
  tree.read_ASCII(data.c_str());
  tree.branch();
  D = tree.get_D();
  tree.set_hsmooth( K, M0 );

  // ---------------------------------------------------------------

  ifstream points( locations.c_str() );
  if(points.good()==false) // open file
    {
      printf("> ERROR: File '%s' not found\n\n\a",locations.c_str());
      return(-1);
    }
  while( points.peek()=='#' ) points.ignore(INT_MAX,'\n'); // ignore comments

  ofstream file_smooth("FiEstAS.txt");

  clock_t t_0 = clock();
  printf(" Computing densities..."); fflush(stdout);

  string str;
  while( getline(points,str) ) // read data
    {
      istringstream s(str);
      DATA x[D];
      for(DIM d=0; d<D; d++) s >> x[d];
      file_smooth << tree.FiEstAS(x,balloon,true) <<endl;
    }
  printf(" Done! (%g s)\n", (clock()-t_0)/double(CLOCKS_PER_SEC) );

  file_smooth.close();
  delete K;
  
  // ---------------------------------------------------------------

  yINFO(("\n-----------------------------------------"));
  yINFO(("\n Program finshed OK :^) "));
  yINFO(("\n (total time = %g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) ));
  yINFO(("\n ... Paranoy@ Rulz! "));
  yINFO(("\n-----------------------------------------\n\n"));
  return(0);
}
