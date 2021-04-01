#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

#include"FiEstAS.h"

static const double M_Hernq = 1.; // Msun
static const double a_Hernq = 1.; // kpc
static const double G = 4.301179e-6; // (km/s)^2 kpc/Msun
static const double GMa = G*M_Hernq/a_Hernq;

//----------------------------------------------------------------------
double Hernquist_3D(DATA *x, double centre)
//----------------------------------------------------------------------
{
  double dx, r=0.;
  for(DIM d=0;d<3;d++){ dx= *(x++)-centre; r+=dx*dx; }
  r=sqrt(r)/a_Hernq;
  dx=(1+r)*a_Hernq;
  return( M_Hernq/2/M_PI/r/dx/dx/dx );
}
//----------------------------------------------------------------------
double Hernquist_6D(DATA *x, double centre)
//----------------------------------------------------------------------
{
  double dx, r2=0., v2=0.;
  for(DIM d=0;d<3;d++){ dx= *(x++)-centre; r2+=dx*dx; }
  for(DIM d=0;d<3;d++){ dx= *(x++);        v2+=dx*dx; }
  double phi = a_Hernq / ( a_Hernq+sqrt(r2) );
  double E = phi - v2/GMa/2;

  if(E<=0.)
    {
      for(DIM d=0;d<6;d++) cerr<<' '<<*(--x);
      cerr<<" : "<<"phi("<<sqrt(r2)<<")="<<phi<<", v2="<<v2/GMa<<"; E="<<E<<endl;
      return(1e-99);
    }

  double f;
  dx = a_Hernq*M_PI;
  f = M_Hernq/4/dx/dx/dx * pow(2*GMa,-1.5) * pow( 1-E, -2.5 );
  f *= 3*asin(sqrt(E)) + sqrt( E*(1-E) )*( 1-2*E )*( 8*E*(E-1)-3 );
  return(f);
}


//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  clock_t t0 = clock();
  srand48(time(NULL));

  yINFO(("\n------------------------------------------------------------------------------"));
  yINFO(("\n Compute densities from the Gadget file <data>, evaluated at <locations>. \n"));
  yINFO(("\n    Usage : FiEstAS_Gadget <snapshot> <part_type> <data> ... [-at<part_type>]\n"));
  yINFO(("\n    Output: D<D>_points.txt D<D>_FiEstAS.txt D<D>_trueRho.txt \n"));
  yINFO(("\n Yago Ascasibar (UAM, November 2007)"));
  yINFO(("\n------------------------------------------------------------------------------\n\n"));

  ERROR(argc<4,("Wrong syntax"));

  ySnapshot::PartType part_type = ySnapshot::get_ptype(argv[2]);
  ySnapshot::PartType at_part_type = part_type;
  if( strncmp("-at",argv[argc-1],3)==0 )
    at_part_type = ySnapshot::get_ptype( argv[--argc]+3 );

  vector<ySnapshot::Data> data;
  DIM D=0;
  for(int i=3; i<argc; i++)
    {
      ySnapshot::Data d = ySnapshot::get_data(argv[i]);
      data.push_back(d);
      D += ySnapshot::get_dim(d);
    }

  yINFO((" <snapshot> = %s \n",argv[1]));
  yINFO((" <part_type> = %s \n",argv[2]));
  yINFO((" <D> = %d \n", D));

  // ---------------------------------------------------------------

  yTree tree(D);
  tree.read_Gadget(argv[1],data,part_type);
  tree.branch();
  
//   TopHat kernel;
  Epanechnikov kernel;
  tree.set_hsmooth(&kernel, 2.);
  bool balloon = false;
  
//   vector<NUM> dim(3);
//   vector<double> scale(3); scale[0] = scale[1] = scale[2] = 1.;
//   dim[0]=0; dim[1]=1; dim[2]=2; tree.set_metric(dim,scale);
//   dim[0]=3; dim[1]=4; dim[2]=5; tree.set_metric(dim,scale);

  // ---------------------------------------------------------------
  char buf[50];
  snprintf(buf,49,"D%d_points.txt",D); ofstream file_points(buf);
  snprintf(buf,49,"D%d_trueRho.txt",D); ofstream file_true(buf);
  snprintf(buf,49,"D%d_FiEstAS.txt",D); ofstream file_rho(buf);
  snprintf(buf,49,"D%d_smooth.txt",D); ofstream file_smooth(buf);

  yGadget snap(argv[1]);

  NUM N = snap.Npart( at_part_type );
  NUM Nread = 0;
  static const NUM Npage = N/200 + 1000;
  float points[D*Npage];
  NUM Ndata = data.size();

  static const double m_part = snap.Mpart(part_type);
  yINFO((" m_part = %g \n", m_part));

  clock_t t_0 = clock();
  printf(" Computing densities...   %%"); fflush(stdout);
  while(Nread<N)
    {
      printf("\b\b\b\b%3d%%",(100*Nread)/N); fflush(stdout);

      NUM Nmin = Npage;
      NUM offset[Ndata], off=0;
      for(NUM i=0; i<Ndata; i++)
	{
	  snap.read( data[i], at_part_type, Nread );
	  NUM n = snap.next( (char *)(points+off), Npage );
	  offset[i]=off;
	  off += n*snap.get_dim(data[i]);
	  if(n<Nmin) Nmin=n;
	}
      Nread += Nmin;

      for(NUM i=0; i<Nmin; i++)
	{
	  DATA x[D];
	  DIM d=0;
	  for(NUM j=0; j<Ndata; j++)
	    {
	      DIM Dj=snap.get_dim(data[j]);
	      for(DIM dj=0; dj<Dj; dj++)
		{
		  x[d] = points[ offset[j]+Dj*i+dj ];
		  file_points<<' '<<x[d];
		  d++;
		}
	    }
	  file_points<<endl;

	  switch(D)
	    {
	    case 3: file_true << Hernquist_3D(x,snap.Lbox()/2) <<endl; break;
	    case 6: file_true << Hernquist_6D(x,snap.Lbox()/2) <<endl; break;
	    }

	  vector<DATA> box; tree.push_box(x,&box);
	  double V=1;
	  for(DIM d=0; d<D; d++) V *= box[2*d+1]-box[2*d];
	  file_rho << m_part/V <<endl;

	  //file_smooth << m_part * tree.smooth(x) <<endl;
	  //printf("%d\n",i);
	  file_smooth << tree.FiEstAS(x,balloon,true) <<endl;
	}

    }
  printf("\b\b\b\b Done! (%g s)\n", (clock()-t_0)/double(CLOCKS_PER_SEC) );

  file_points.close();
  file_true.close();
  file_rho.close();
  file_smooth.close();

  // ---------------------------------------------------------------

  yINFO(("\n-----------------------------------------"));
  yINFO(("\n Program finshed OK :^) "));
  yINFO(("\n (total time = %g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) ));
  yINFO(("\n ... Paranoy@ Rulz! "));
  yINFO(("\n-----------------------------------------\n\n"));

  return(0);
}
