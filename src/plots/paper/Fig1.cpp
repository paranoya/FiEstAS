#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
using namespace std;

#include"FiEstAS.h"
#include"cpgplot.h"

//----------------------------------------------------------------------
void compute_map( yTree &tree, float* rho,
		  float x0,float x1, float y0,float y1,
		  int Npix, char panel )
//----------------------------------------------------------------------
{
  float dx = (x1-x0)/Npix;
  float dy = (y1-y0)/Npix;
  for(int ix=0; ix<Npix; ix++)
    for(int iy=0; iy<Npix; iy++)
      {
	DATA point[2];
	point[0] =  x0+(ix+.5)*dx;
	point[1] =  y0+(iy+.5)*dy;
	switch( panel )
	  {
	    case 'a':
	    case 'b':
	    case 'c': rho[ ix + iy*Npix ] = tree.smooth(point); break;
	    case 'd': rho[ ix + iy*Npix ] = tree.FiEstAS(point); break;
	  }
      }
}

//----------------------------------------------------------------------
void plot_map(float* data, int size)
//----------------------------------------------------------------------
{
  int size2 = size*size;

  double intensity[size2];
  for(int i=0; i<size2; i++ ) intensity[i] = data[i];
  sort( intensity, intensity+size2 );

  double mass[size2], m=0.;
  for(int i=0; i<size2; i++ ) mass[i] = ( m += intensity[i] );
  for(int i=0; i<size2; i++ ) mass[i] /= m;

  double mean=0., sigma=0., min=yHUGE, max=-yHUGE;
  for(int i=0; i<size2; i++ )
    {
      double x = data[i];

      mean += x;
      sigma += x*x;
      if(x<min) min=x;
      if(x>max) max=x;

      data[i] = mass[ lower_bound( intensity, intensity+size2, x )-intensity ];
    }
  mean /= size*size;
  sigma = sqrt( sigma/size/size - mean*mean );
  printf(" Data: mean=%g sigma=%g min=%g max=%g \n",mean,sigma,min,max);

  float x0,x1, y0,y1;
  cpgqwin(&x0,&x1, &y0,&y1);

  cpgswin(0,size, 0,size);

  static const int n_levels = 7;
  float levels[n_levels] = { 0.0, .05, .25, .5, .75, .95, 1. };
  float tr[6]={-.5,1.,0., -.5,0.,1.};
  //float tr[6]={-.5,0.,1., -.5,1.,0.};

  //cpggray(data,size,size, 1,size,1,size, 1.,0., tr);

  float R[n_levels] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0 };
  float G[n_levels] = { 1.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0 };
  float B[n_levels] = { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 };
  cpgctab(levels,R,G,B,n_levels, 1.0,0.5);
  cpgimag(data,size,size, 1,size,1,size, 0.,1., tr);

  cpgcont(data,size,size, 1,size,1,size, levels+1,n_levels-2, tr);

  cpgswin(x0,x1, y0,y1);
}

//----------------------------------------------------------------------
void plot_2D( yTree &tree, int Npix, char panel )
//----------------------------------------------------------------------
{
  float x0 = tree.get_min(0);
  float x1 = tree.get_max(0);
  float lx = .2*(x1-x0); x0 -=lx; x1+=lx;
  float y0 = tree.get_min(1);
  float y1 = tree.get_max(1);
  float ly = .2*(y1-y0); y0 -=ly; y1+=ly;
  cpgswin(x0,x1, y0,y1);
  printf(" x=[%g,%g] y=[%g,%g], 1/V=%g\n",x0,x1,y0,y1,1./(x1-x0)/(y1-y0));

  float rho[Npix*Npix];
  compute_map( tree, rho, x0,x1, y0,y1, Npix, panel );
  plot_map( rho, Npix );
}

//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  clock_t t0 = clock();

  yINFO(("\n-----------------------------------------------------"));
  yINFO(("\n Plot Figure 'Fig1.eps' \n"));
  yINFO(("\n    Usage  : fig1 [options] <file> <panel> \n"));
  yINFO(("\n    Options: -print : (over)writes PostScript file 'Fig1<panel>.eps' \n"));
  yINFO(("\n Yago Ascasibar (UAM, Fall 2008)"));
  yINFO(("\n-----------------------------------------------------\n\n"));

  ERROR( argc<3, ("Wrong syntax") );

  NUM col_x = 0;
  NUM col_y = 1;
  yINFO((" Plotting columns %d and %d of file '%s' \n", col_x, col_y, argv[argc-2] ));

  char panel = argv[argc-1][0];
  ERROR( panel<'a' || panel>'d' , ("Wrong panel '%s'",argv[argc-1]) );

  // ---------------------------------------------------------------
  // Compute FiEstAS tree

  yTree tree(2);

  yASCII_data x_i( argv[argc-2] );
  NUM N = x_i.get_Nrows();
  for(NUM row=0; row<N; row++) tree.add_xd( 0, x_i.get(row,col_x) );
  for(NUM row=0; row<N; row++) tree.add_xd( 1, x_i.get(row,col_y) );

  tree.branch();

  yKernel *K = NULL;
  double M0 = 2.;
  switch( panel )
  {
    case 'a':
    case 'd': K = new TopHat();       M0 = 2.; break;
    case 'b': K = new Epanechnikov(); M0 = 2.; break;
    case 'c': K = new Epanechnikov(); M0 = 10.; break;
  }
  tree.set_hsmooth( K, M0 );

  // ---------------------------------------------------------------
  // Init plot

  int Npix = 600;
  float x,y;
  char cc;
  string print = option<string>( "-prin", "", argc-2,argv );
  if( print == "t" )
  {
    print = "Fig1#.eps/vcps";
    print[4] = panel;
    cpgopen( print.c_str() );
    cc='q';
    Npix=100;
  }
  else
  {
    cpgopen("/xwin");
    cc=0;
  }

  cpgask(0);
  cpgsci(1);
  cpgslw(1);
  cpgsch(.5);
  cpgpap(0.,2.);

  // ---------------------------------------------------------------
  // Plot 2D color map with contours

  cpgsvp( .1,.9, .5,.9 );

  plot_2D( tree, Npix, panel );
  cpgbox("bc",0.,0,"bc",0.,0);

  // Plot points
  cpgsci(4); for(NUM n=0; n<N; n++) cpgpt1( tree.get_data(0,n), tree.get_data(1,n), -4); cpgsci(1);

  // Draw ring
//   cpgsfs(0);
//   cpgsls(4); cpgcirc( 0.,0., .95); cpgcirc( 0.,0., 1.05); cpgsls(1);
  
  // ---------------------------------------------------------------
  // Adjust meteric

  vector<NUM> dim(2); dim[0]=0; dim[1]=1;
  vector<double> scale(2); scale[0]=1.; scale[1]=1.;
  tree.set_metric(dim,scale);

  // ---------------------------------------------------------------
  // Plot 2D color map with contours

  cpgsvp( .1,.9, .1,.5 );

  plot_2D( tree, Npix, panel );
  cpgbox("bc",0.,0,"bc",0.,0);

  // Plot points
  cpgsci(4); for(NUM n=0; n<N; n++) cpgpt1( tree.get_data(0,n), tree.get_data(1,n), -4); cpgsci(1);

  // Draw ring
//   cpgsfs(0);
//   cpgsls(4); cpgcirc( 0.,0., .95); cpgcirc( 0.,0., 1.05); cpgsls(1);
  
  // ---------------------------------------------------------------
  // Write label

  cpgsch(1.5);
  switch( panel )
    {
      case 'a': cpgsvp( .1,.9, .5,.9 ); cpglab("","","a) Top-Hat kernel, M\\d0\\u=2.0"); break;
      case 'b': cpgsvp( .1,.9, .5,.9 ); cpglab("","","b) Epanechnikov kernel, M\\d0\\u=2.0"); break;
      case 'c': cpgsvp( .1,.9, .5,.9 ); cpglab("","","c) Epanechnikov kernel, M\\d0\\u=10.0"); break;
      case 'd': cpgsvp( .1,.9, .5,.9 ); cpglab("","","d) Top-hat + balloon, M\\d0\\u=2.0"); break;
    }

  // ---------------------------------------------------------------
  // bye

  while(cc!='q' && cc!='Q' && cc!=27) cpgcurs(&x,&y,&cc);
  cpgclos();
  delete K;
  
  yINFO(("\n-----------------------------------------"));
  yINFO(("\n Program finshed OK :^) "));
  yINFO(("\n (total time = %g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) ));
  yINFO(("\n ... Paranoy@ Rulz! "));
  yINFO(("\n-----------------------------------------\n\n"));
  return(0);
}
