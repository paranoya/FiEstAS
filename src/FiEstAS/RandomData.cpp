#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

//----------------------------------------------------------------------
double Poisson(double* x, int D)
// Generate a random Poisson distribution [0,1] and return density
//----------------------------------------------------------------------
{
  for(int d=0;d<D;d++) x[d] = drand48();
  return(1.);
}

//----------------------------------------------------------------------
double Ring(double* x, int D)
// Generate a hyper-cylinder of radius 1+-0.05 centered at (0,0)
// All other coordinates run from -0.05 to 0.05
// Return density
//----------------------------------------------------------------------
{
  static const double R0 = .95;
  static const double R1 = 1.05;
  static const double R02 = R0*R0;
  static const double dR2 = R1*R1 - R02;
  double r = sqrt( R02 + dR2*drand48() );
  double z = 2.*M_PI*drand48();
  x[0] = r*cos(z);
  x[1] = r*sin(z);
  for(int d=2;d<D;d++) x[d] = .1*drand48()-.05;

  return( 1. / ( M_PI*dR2 * pow(.1,D-2.) ) );
}

//----------------------------------------------------------------------
double Hernquist(double* xx, int D)
// Generate a 3D Hernquist sphere with a=1,
// centered at (0,.75,0) and clipped by [-.5,1].
// TO DO: If D>=6, generate also velocities in units
//        of GM/a=1, clipped by [-1,1].
// Any extra dimensions are set to uniform [0,1].
// Returns rho=0 if particle has been clipped.
//----------------------------------------------------------------------
{
  double mu, r, fi, cosz,sinz;
  double x,y,z;

  mu = sqrt( drand48() );
  r = mu/(1.-mu);

  fi = 2.*M_PI*drand48();
  cosz = 2.*drand48()-1.;
  sinz = sqrt(1.-cosz*cosz);

  x = r*sinz*cos(fi);
  y = r*sinz*sin(fi);
  z = r*cosz;

  xx[0]=x; xx[1]=y; xx[2]=z;
  for(int d=3;d<D;d++) xx[d] = drand48();
  
  if(x<-.5 || x>1. || y<-.5 || y>1. || z<-.5 || z>1.) return(0.);
  else return( 1./2./M_PI/r/pow(1.+r,3.) );
}

//----------------------------------------------------------------------
int main(int argc,char** argv)
//----------------------------------------------------------------------
{
  static const char out1[15] = "RandomData.txt";
  static const char out2[12] = "TrueRho.txt";
  clock_t t0 = clock();

  printf("\n-----------------------------------------------------");
  printf("\n Generates a sample of <N> random points from the ");
  printf("\n distribution <distrib> = Poisson / Ring / Hernquist ");
  printf("\n in <D> dimensions. \n");
  printf("\n    Usage : RandomData <N> <distrib> <D> \n");
  printf("\n    Output: %s   ",out1);
  printf("\n            %s \n",out2);
  printf("\n Yago Ascasibar (AIP, Summer 2007)");
  printf("\n-----------------------------------------------------\n\n");

  if(argc!=4)
    {
      printf(" ERROR: Wrong syntax \n\n\a");
      return(-1);
    }

  // ---------------------------------------------------------------

  int N = atoi(argv[1]);
  int D = atoi(argv[3]);
  if(N<1 || D<1)
    {
      printf("\aERROR: N=%d, D=%d !\n",N,D);
      return(-1);
    }

  double (*dist)(double*,int);

  if(strncmp(argv[2],"Ring",5)==0)
    {
      if(D<2)
	{
	  printf("\aERROR: D must be >=2 for this distribution\n");
	  return(-1);
	}
      else dist = Ring;
    }
  else
    if(strncmp(argv[2],"Hernquist",5)==0)
      {
	if(D<3)
	  {
	    printf("\aERROR: D must be >=3 for this distribution\n");
	    return(-1);
	  }
	else dist = Hernquist;
      }
    else
      if(strncmp(argv[2],"Poisson",8)==0) dist = Poisson;
      else
	{
	  printf("\aERROR: Invalid distribution '%s' \n",argv[2]);
	  return(-1);
	}

  printf(" Generating n=%d points in d=%d dimensions, \n",N,D);
  printf(" according to a %s distribution. \n",argv[2]);

  // ---------------------------------------------------------------

  long int seed=time(NULL);
  srand48(seed);

  FILE *fout1;
  if((fout1=fopen(out1,"w"))==NULL)
    {
      printf("\aERROR: Cannot open '%s'. \n",out1);
      return(-1);
    }
  fprintf(fout1,"#---- RANDOM DATA ------\n");
  fprintf(fout1,"# %s distribution\n",argv[2]);
  fprintf(fout1,"# D=%d, N=%d\n",D,N);
  fprintf(fout1,"# seed=%ld\n",seed);
  fprintf(fout1,"#-----------------------\n");

  FILE *fout2;
  if((fout2=fopen(out2,"w"))==NULL)
    {
      printf("\aERROR: Cannot open '%s'. \n",out2);
      return(-1);
    }
  fprintf(fout2,"#---- TRUE DENSITY -----\n");
  fprintf(fout2,"# %s distribution\n",argv[2]);
  fprintf(fout2,"# D=%d, N=%d\n",D,N);
  fprintf(fout2,"#                       \n");
  fprintf(fout2,"#-----------------------\n");

  for(int n=0; n<N; n++)
    {
      double x[D], rho;
      rho = dist(x,D);
      if(rho>0.)
	{
	  for(int d=0;d<D;d++) fprintf(fout1," %g",x[d]);
	  fprintf(fout1,"\n");
	  fprintf(fout2,"%g\n",rho);
	}
    }
  fclose(fout1);
  fclose(fout2);

  printf("\n-----------------------------------------");
  printf("\n Program finshed OK :^) ");
  printf("\n (total time = %g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) );
  printf("\n ... Paranoy@ Rulz! ");
  printf("\n-----------------------------------------\n\n");
  fflush(stdout);
  return(0);
}
