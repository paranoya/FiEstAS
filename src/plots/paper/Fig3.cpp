#include<algorithm>
using namespace std;

#include"FiEstAS.h"
#include"plot.h"

static const int NMAX = 6;
double average[2][NMAX][5];
double dispersion[2][NMAX][5];

//----------------------------------------------------------------------
void histogram(const char * base_dir, const char panel, bool Euclidean, const char * file_true)
//----------------------------------------------------------------------
{
  // ---------------------------------------------------------------
  // Read files
  
  char file[100];
  if( Euclidean ) snprintf(file,99,"%s/F%c_metric.txt", base_dir, panel);
  else snprintf(file,99,"%s/F%c.txt", base_dir, panel);
  yASCII_data rho( file );
  snprintf(file,99,"%s/%s", base_dir, file_true);
  yASCII_data rho0( file );
  
  NUM N = rho.get_Nrows();
  ERROR( N!=(NUM)rho0.get_Nrows() , ("Nrows=(%d,%d)",N,rho0.get_Nrows()) );
  // ---------------------------------------------------------------
  // Compute ratio rho / rho_true
  
  vector<double> ratio(N);
  double ave=0., sig=0.;
  for(NUM row=0; row<N; row++)
  {
    double q = log10( rho.get(row,0) / rho0.get(row,0) );
    ratio[row] = q ;
    ave += q;
    sig += q*q;
  }
  ave /= N;
  sig = sqrt( sig/N - ave*ave );
  printf(" <log(q)> = %g +- %g \n", ave, sig );
  
  int index_n = (int)(log10(N)+.5);
  int index_p = (int)(panel-'a');
  average[(int)Euclidean][index_n][index_p] = ave;
  dispersion[(int)Euclidean][index_n][index_p] = sig;
  
  // ---------------------------------------------------------------

  plot_histogram( ratio, pow(N,.5) );

}

//----------------------------------------------------------------------
void plot_panel(float x0, float x1, char panel)
//----------------------------------------------------------------------
{
  // ---------------------------------------------------------------
  // Plot upper panel
  
  cpgsvp( x0,x1, .5,.9 );
  if(panel=='a')
  {
    cpgbox("bcts",0.,0,"bctsvn",0.,0);
    legend();
  }
  else cpgbox("bcts",0.,0,"bctsv",0.,0);
  
  cpgsci(4); cpgsls(4); histogram( "../data/Hernquist/1e2", panel,false, "D6_trueRho.txt" );
  cpgsci(3); cpgsls(1); histogram( "../data/Hernquist/1e3", panel,false, "D6_trueRho.txt" );
  cpgsci(2); cpgsls(3); histogram( "../data/Hernquist/1e4", panel,false, "D6_trueRho.txt" );
  cpgsci(1); cpgsls(2); histogram( "../data/Hernquist/1e5", panel,false, "D6_trueRho.txt" );
  
  cpgsls(4); cpgmove(0.,0.); cpgdraw(0.,10.); cpgsls(1);

  cout<<endl;
  
  // ---------------------------------------------------------------
  // Plot lower panel

  cpgsvp( x0,x1, .1,.5 );
  if(panel=='a')
  {
    cpgbox("bctsn",0.,0,"bctsvn",0.,0);
//     legend();
  }
  else cpgbox("bctsn",0.,0,"bctsv",0.,0);
  
  cpgsci(4); cpgsls(4); histogram( "../data/Hernquist/1e2", panel,true, "D6_trueRho.txt" );
  cpgsci(3); cpgsls(1); histogram( "../data/Hernquist/1e3", panel,true, "D6_trueRho.txt" );
  cpgsci(2); cpgsls(3); histogram( "../data/Hernquist/1e4", panel,true, "D6_trueRho.txt" );
  cpgsci(1); cpgsls(2); histogram( "../data/Hernquist/1e5", panel,true, "D6_trueRho.txt" );
  
  cpgsls(4); cpgmove(0.,0.); cpgdraw(0.,10.); cpgsls(1);

  // ---------------------------------------------------------------
  // Write label

//   cpgsch(1.5);
  switch( panel )
  {
    case 'a': cpgsvp( x0,x1, .1,.88); cpglab("","","a) Top-hat kernel, M\\d0\\u=2.0"); break;
    case 'b': cpgsvp( x0,x1, .1,.88); cpglab("","","b) Epanechnikov kernel, M\\d0\\u=2.0"); break;
    case 'c': cpgsvp( x0,x1, .1,.88); cpglab("","","c) Epanechnikov kernel, M\\d0\\u=10.0"); break;
    case 'd': cpgsvp( x0,x1, .1,.88); cpglab("","","d) Top-hat + balloon, M\\d0\\u=2.0"); break;
  }
}

//----------------------------------------------------------------------
void print_block(int b)
//----------------------------------------------------------------------
{
  printf(" $100$");
  for(int i=0;i<4;i++) printf(" & $%5.2f \\pm %4.2f $", average[b][2][i], dispersion[b][2][i]);
  printf(" \\\\\n$1000$");
  for(int i=0;i<4;i++) printf(" & $%5.2f \\pm %4.2f $", average[b][3][i], dispersion[b][3][i]);
  for(int n=4; n<NMAX; n++)
  {
    printf(" \\\\\n$10^%d$",n);
    for(int i=0;i<4;i++) printf(" & $%5.2f \\pm %4.2f $", average[b][n][i], dispersion[b][n][i]);
  }
  printf(" \\\\ \\hline\n");
}
//----------------------------------------------------------------------
void print_tex()
//----------------------------------------------------------------------
{
  printf("$N$ & Top-hat & Epanechnikov & Epa., $M_0=10$ & Top-hat+balloon \\\\ \\hline\n");
  print_block(0);
  print_block(1);
}

//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  clock_t t0 = clock();

  yINFO(("\n-----------------------------------------------------"));
  yINFO(("\n Plot Figure 'Fig3.eps' \n"));
  yINFO(("\n    Usage  : fig3 [-print] \n"));
  yINFO(("\n    Output : Fig3.eps \n"));
  yINFO(("\n Yago Ascasibar (UAM, Fall 2008)"));
  yINFO(("\n-----------------------------------------------------\n\n"));

  // ---------------------------------------------------------------
  // Init plot

  float x,y;
  char cc;
  string print = option<string>( "-prin", "", argc,argv );
  if( print == "t" )
  {
    cpgopen("Fig3.eps/vcps");
    cc='q';
  }
  else
  {
    cpgopen("/xwin");
    cc=0;
  }

  cpgask(0);
  cpgsci(1);
  cpgslw(1);
  cpgpap(0.,.5);
  cpgsch(.9);
  
  cpgswin( -1.25,1.25, 0.,4.5 );
  
  plot_panel(.1,.3,'a');
  plot_panel(.3,.5,'b');
  plot_panel(.5,.7,'c');
  plot_panel(.7,.9,'d');
  
  cpgsvp( .1,.9, .14,.9 );
  cpgsch(1.5);
  cpglab(" q\\di\\u = log( f\\dFiEstAS\\u / f\\dtrue\\u )", "p(q\\di\\u) dq\\di\\u", "" );

  // ---------------------------------------------------------------
  // bye

  while(cc!='q' && cc!='Q' && cc!=27) cpgcurs(&x,&y,&cc);
  cpgclos();
  
  print_tex();
  
  yINFO(("\n-----------------------------------------"));
  yINFO(("\n Program finshed OK :^) "));
  yINFO(("\n (total time = %g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) ));
  yINFO(("\n ... Paranoy@ Rulz! "));
  yINFO(("\n-----------------------------------------\n\n"));
  return(0);
}
