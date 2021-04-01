#include"FiEstAS.h"

#include<math.h>

//----------------------------------------------------------------------
void yTree::branch()
//----------------------------------------------------------------------
{
  Ntot = x[0].size();

  bool error=false;
  DEBUG((">> Branching starts: Nx ="));

  for(DIM d=0; d<D; d++)
    {
      NUM n = x[d].size();
      DEBUG((" %d", n));
      if(n!=Ntot) error=true;
    }
  DEBUG(("\n>>                   Nextra ="));
  for(DIM d=0; d<D_extra; d++)
    {
      NUM n = extra[d].size();
      DEBUG((" %d", n));
      if(n!=Ntot) error=true;
    }
  DEBUG(("\n"));
  ERROR(error,(" Not all vectors have the same length !!!"));

  // ---------------------------------------------------------------

  init_scratch1( (NUM)( 1+sqrt((double)Ntot) ) );
  node = 0;

  // ---------------------------------------------------------------

  clock_t t0 = clock();
  DEBUG(("\n---------------------------------------\n"));
  if(!silent) yINFO((" yTree: D=%d, N=%d ",D,Ntot); fflush(stdout));

  branch(0,Ntot,min,max);

  if(!silent) yINFO((" (%g s) \n", (clock()-t0)/double(CLOCKS_PER_SEC) ));
  DEBUG(("---------------------------------------\n\n"));

  // ---------------------------------------------------------------

  init_scratch2();
}

//----------------------------------------------------------------------
void yTree::branch(NUM n0, NUM n1, DATA *x0, DATA *x1)
//----------------------------------------------------------------------
{
  ERROR(n0==n1,("IMPOSSIBLE! n0=%d = n1=%d, x[0]=[%g,%g]",n0,n1,x0[0],x1[0]));
  if(n1==n0+1)
    {
      PANIC(
      printf("LEAF: x(%d) = (",n0); for(DIM d=0; d<D; d++) printf(" %10g",x[d][n0]);
      printf(") Nsplits = ("); for(DIM d=0; d<D; d++) printf(" %3d",Nsplits[d]);
      printf(") \n");
      );
    }
  else
    {
      choose_axis(n0,n1,x0,x1);
      spawn_children(n0,n1);
    }
}

//----------------------------------------------------------------------
double lfact(NUM n) // log(n!)
{
  static const NUM N = 20;
  static const double lf[N] ={  0., .693147181, 1.79175947, 3.17805383, 4.78749174, 6.57925121, 8.52516136, 10.6046029, 12.8018275, 15.1044126, 17.5023078, 19.9872145, 22.5521639, 25.1912212, 27.8992714, 30.6718601, 33.5050735, 36.3954452, 39.3398842 };
  if(n<N) return(lf[n]);
  else return( n*log((double)n) -n+1 ); // Stirling's formula
}
//----------------------------------------------------------------------
void yTree::choose_axis(NUM n0, NUM n1, DATA *x_0, DATA *x_1)
//----------------------------------------------------------------------
{
  const NUM N = n1-n0; // Number of particles in the node
  const NUM B = (NUM)( 1+sqrt((double)N) ); // Number of bins

  PANIC( printf("Node %d-%d (N=%d, B=%d) \n",n0,n1, N,B); );

  double Lmin = 1e300; // log(Likelihood)
  for(DIM d=0; d<D; d++)
    {
      NUM b;
      for(b=0; b<B; b++){ Nb[b]=0; xmin[b]=yHUGE; xmax[b]=-yHUGE; } // init histogram

      const DATA x0=x_0[d], x1=x_1[d];
      const double l = x1-x0;
      if(l<=0)
	{
	  PANIC( printf("Skipping dimension %d (l=%g) \n",d,l); );
	  continue;
	}

      for(NUM i=n0; i<n1; i++ ) // add x_i to histogram
	{
	  DATA xi = x[d][i];
	  DEBUG_ERROR(xi<x0 || xi>x1, (" IMPOSSIBLE! x=%g, out of [%g,%g]",xi,x0,x1) );

	  b = (NUM)( B*((xi-x0)/l) );
	  //PANIC( printf(" x_%d[%d]=%g [%g,%g] => bin=%d/%d \n", i,d,xi, x0,x1, b,B); );

	  DEBUG_ERROR(b<0 || b>B, (" IMPOSSIBLE! b = %d",b) );
	  if(b==B)
	    {
	      b--;
	      yWARNING(xi!=x1, ("b=%d (%g!=%g)",b,xi,x1) );
	    }
	  if(xi==x0) b=0;
	  if(xi==x1) b=B-1;

	  Nb[b]++;
	  if(xi>xmax[b]) xmax[b]=xi;
	  if(xi<xmin[b]) xmin[b]=xi;
	}

      // double L = lfact(N) - N*log((double)B); // log(Likelihood)
      double L = 0.; // I don't care about the offset ( it's equal for all d )
      NUM n;
      for(b=0; b<B; b++) if( (n=Nb[b])>0 ) L -= lfact(n);

//       L += Nsplits[d]; // penalize dimensions that have been divided many times
      PANIC( printf("   L[%d]=%g; \n", d,L); );

      if(L < Lmin)
	{
	  n=0; b=B;
	  while(n<=N/2 && b>1) n+=Nb[--b];
	  if(Nb[b]<2*n-N && b<B-1) n-=Nb[b++];
	  NUM bl=b-1, br=b;
	  while(Nb[bl]==0) bl--;
	  while(Nb[br]==0) br++;
	  //PANIC( printf("   Nr=%d/%d (b=[%d-%d]/%d, Nb=%d, l=%g) \n",n,N, bl,br,B, Nb[b], l); );

	  Lmin = L;
	  ds = d;
	  Nr = n;
	  xs = (double)(xmax[bl]+xmin[br])/2.;
	}
    }
  PANIC( printf("   Likelihood[%d]=%g; x_split=%g (Nr=%d/%d) \n", ds,Lmin, xs, Nr,N); );
  if(Lmin == 1e300)
    {
      yWARNING(true,("%d equal points", N));
      ds = D;
      Nr = 1;
      xs = 0.;
      return;
    }

  if(xs==x_1[ds])
    {
      //printf("   x_split(=%g) == x_right(=%g) \n",xs,x_1[ds]);
      //printf("   Nsplits = ("); for(DIM d=0; d<D; d++) printf(" %3d",Nsplits[d]);
      //printf(") \n");
      tilt_left=false;
      ERROR(xs==x_0[ds],(" ERROR: x_split=%g [%g-%g]",xs,x_0[ds],x_1[ds]));
    }
  if(xs==x_0[ds])
    {
      //printf("   x_split(=%g) == x_left(=%g) \n",xs,x_0[ds]);
      //printf("   Nsplits = ("); for(DIM d=0; d<D; d++) printf(" %3d",Nsplits[d]);
      //printf(") \n");
      tilt_left=true;
      ERROR(xs==x_1[ds],("x_split=%g [%g-%g]",xs,x_0[ds],x_1[ds]));
    }
}

//----------------------------------------------------------------------
void yTree::spawn_children(NUM n0, NUM n1)
//----------------------------------------------------------------------
{
  NUM l=n0, r=n1;
  NUM nr=0;
  DATA lmin[D], lmax[D];
  DATA rmin[D], rmax[D];
  for(DIM d=0; d<D; d++){ lmin[d]=yHUGE; lmax[d]=-yHUGE; rmin[d]=yHUGE; rmax[d]=-yHUGE; }

  if(ds==D) // Several equal points
    {
      nr++; r--; // one to the right
      l = n1-1; // rest to the left
    }
  else
    {
      for(NUM i=0; i<(n1-n0); i++)
	{
	  DATA xi=x[ds][l];
	  // printf("   x=%g",xi); fflush(stdout);
	  if( xi<xs || (xi==xs && tilt_left) )
	    {
	      // printf(" => Left (%d) \n",n1-n0-nr);
	      for(DIM d=0; d<D; d++)
		{
		  xi=x[d][l];
		  if(xi<lmin[d]) lmin[d]=xi;
		  if(xi>lmax[d]) lmax[d]=xi;
		}
	      l++;
	    }
	  else
	    {
	      // printf(" => Right (%d) \n",nr);
	      nr++;
	      r--;
	      for(DIM d=0; d<D; d++)
		{
		  xi=x[d][l];
		  if(xi<rmin[d]) rmin[d]=xi;
		  if(xi>rmax[d]) rmax[d]=xi;
		  x[d][l]=x[d][r];
		  x[d][r]=xi;
		}
	      for(DIM d=0; d<D_extra; d++)
		{
		  xi=extra[d][l];
		  extra[d][l]=extra[d][r];
		  extra[d][r]=xi;
		}
	    }
	}
      //if(nr!=Nr) fprintf(stderr,"   Warning: nr=%d (!=%d) \n",nr,Nr);
    }

  d_s[node] = ds;
  x_s[node] = xs;
  N_l[node] = n1-n0-nr;
  if(nr==0 || N_l[node]==0) for(NUM i=n0;i<n1;i++) printf("%g\n",x[ds][i]);
  ERROR(nr==0 || N_l[node]==0,("\n IMPOSSIBLE ! N_left=%d, N_right=%d (ds=%d xs=%g)",N_l[node],nr,ds,xs));

  node++;

  nr=ds;
  Nsplits[nr]++;
  branch(n0,l,lmin,lmax);
  branch(r,n1,rmin,rmax);
  Nsplits[nr]--;
}

//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
