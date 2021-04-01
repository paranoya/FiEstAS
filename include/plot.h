#include"cpgplot.h"

//----------------------------------------------------------------------
void legend()
//----------------------------------------------------------------------
{
  cpgsci(4); cpgsls(4); cpgmove(-1.1,4.0); cpgdraw(-.75,4.0); cpgtext(-.7,4.0,"N=100");
  cpgsci(3); cpgsls(1); cpgmove(-1.1,3.7); cpgdraw(-.75,3.7); cpgtext(-.7,3.7,"N=1000");
  cpgsci(2); cpgsls(3); cpgmove(-1.1,3.4); cpgdraw(-.75,3.4); cpgtext(-.7,3.4,"N=10\\u4");
  cpgsci(1); cpgsls(2); cpgmove(-1.1,3.1); cpgdraw(-.75,3.1); cpgtext(-.7,3.1,"N=10\\u5");
  
  cpgsci(1); cpgsls(1);
}

//----------------------------------------------------------------------
void plot_histogram(vector<double> x, double N_steps)
//----------------------------------------------------------------------
{
  NUM N_points = x.size();
//   NUM n_step = (NUM)( N_points / N_steps );
  
  sort( x.begin(), x.end() );
  double x_step = ( x[N_points-1] - x[0] ) / N_steps;
  
  double old_xx = x[0];
  NUM old_n = 0 ;
  cpgmove( old_xx, 0 );
  
  NUM n = 0;
  while( n < N_points )
{
    n++;
    double xx = ( x[n] + x[n-1] )/2;
    NUM N_bin = n-old_n;
//     if( N_bin>=n_step || xx-old_xx>=x_step )
    if( xx-old_xx>=x_step )
{
      double yy = N_bin / ( xx - old_xx ) / N_points;
    
//       cpgdraw( old_xx, yy ); cpgdraw( xx, yy );
      cpgdraw( (xx+old_xx)/2, yy );
    
      old_xx = xx;
      old_n = n;
}
}
  old_xx = x[N_points-1];
//   cpgdraw( old_xx, old_yy );
  cpgdraw( old_xx, 0. );
}

