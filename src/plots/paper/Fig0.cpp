#include"FiEstAS.h"
#include"cpgplot.h"

//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  // ---------------------------------------------------------------
  // Init plot

  float x,y;
  char cc;
  string print = option<string>( "-prin", "", argc,argv );
  if( print == "t" )
  {
    cpgopen("Fig0.eps/vcps");
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
  cpgsch(.5);
  cpgpap(0.,1.);

  cpgswin(-.25,1.25,-1,2);
  cpgsfs(2);

  // ---------------------------------------------------------------
  // Data

  static const NUM N = 7;
  float xx[N] = { .5, .2, .35, .4, .6, .7, 1.1 };
  float yy[N] = { .2, .2, .05, .3, 2., .4, .15 };
  
  // ---------------------------------------------------------------
  // Left panel
  
  cpgsvp( 0,.4, .1,.9 );

  cpgsci(1); for(NUM n=1; n<N; n++) cpgpt1( xx[n], yy[n], -4);

  cpgsci(2); cpgpt1( xx[0], yy[0], -4);
  
  double x_ave   = 0; for(NUM n=0; n<N; n++) x_ave   += xx[n];       x_ave /= N;
  double y_ave   = 0; for(NUM n=0; n<N; n++) y_ave   += yy[n];       y_ave /= N;
  double sigma_x = 0; for(NUM n=0; n<N; n++) sigma_x += xx[n]*xx[n]; sigma_x = sqrt( sigma_x/N - x_ave*x_ave );
  double sigma_y = 0; for(NUM n=0; n<N; n++) sigma_y += yy[n]*yy[n]; sigma_y = sqrt( sigma_y/N - y_ave*y_ave );
  printf(" (%g,%g) +- (%g,%g) : %g \n", x_ave, y_ave, sigma_x, sigma_y, sigma_y/sigma_x );

  cpgrect( xx[0]-sigma_x, xx[0]+sigma_x, yy[0]-sigma_y, yy[0]+sigma_y );

  // ---------------------------------------------------------------
  // Right panel
  
  cpgsvp( .59,.99, .1,.9 );

  cpgsci(1); for(NUM n=1; n<N; n++) cpgpt1( xx[n], yy[n], -4);

  cpgsci(2); cpgpt1( xx[0], yy[0], -4);

  double x0, y0, hx, hy, W;  x0=y0=hx=hy=W= 0;
  for(NUM n=1; n<N; n++)
  {
    double dx = ( xx[n] - xx[0] ) / sigma_x;
    double dy = ( yy[n] - yy[0] ) / sigma_y;
    double ww = exp(-dx*dx/2)/sigma_x * exp(-dy*dy/2)/sigma_y;
    x0 += ww * xx[n];
    y0 += ww * yy[n];
    hx += ww * xx[n]*xx[n];
    hy += ww * yy[n]*yy[n];
    W += ww;
  }
  x0 /= W;
  y0 /= W;
  hx = sqrt( hx/W - x0*x0 );
  hy = sqrt( hy/W - y0*y0 );
  printf(" (%g,%g) +- (%g,%g) : %g \n", x0, y0, hx, hy, hy/hx );

  cpgrect( xx[0]-hx, xx[0]+hx, yy[0]-hy, yy[0]+hy );
  
  // ---------------------------------------------------------------
  // bye

  while(cc!='q' && cc!='Q' && cc!=27) cpgcurs(&x,&y,&cc);
  cpgclos();
}