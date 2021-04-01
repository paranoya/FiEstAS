#include"FiEstAS.h"

static const double EPS = .1;

static const double min_s=.01, max_s=100.;

//##############################################################################

//------------------------------------------------------
void yTree::set_hsmooth(yKernel *K, double M_0)
// compute h for each point (public, must be called before FiEstAS)
//------------------------------------------------------
{
  printf(" Set h_smooth...   %%"); fflush(stdout);
  clock_t t0 = clock();
  init_smooth();

  ERROR( M_0<=0., ("M0=%g",M_0) );
  M0 = M_0;
  Kernel = K;
  bias = 1 + pow( 2*K->Kernel(0), D ) / M0;
  printf("\n\n BIAS=%g \n\n",log10(bias) );
  
  for(NUM i=0, di=0; i<Ntot; i++)
    {
      if(i>di)
	{
	  printf("\b\b\b\b%3d%%",(100*i)/Ntot); fflush(stdout);
	  di+=Ntot/100;
	}
      vector<DATA> bb;
      push_box(i,&bb);
      for(DIM d=0;d<D;d++) h[d][i] = bb[2*d+1] - bb[2*d];

      vector<DATA> cells[2*D];
      for(DIM d=0;d<D;d++)
	{
	  cells[2*d  ].push_back( bb[2*d  ] );
	  cells[2*d+1].push_back( bb[2*d+1] );
	}
      vector<NUM> id_nei;
      NUM N_nei = push_cells(cells,&id_nei,false);
      //printf("Nnei_%d=%d \n",i,N_nei);
      if(N_nei<2) printf("AAAA! Nnei=%d \n",N_nei);

      double sum[D]; for(DIM d=0;d<D;d++) sum[d] = 0.;
      double sum2[D]; for(DIM d=0;d<D;d++) sum2[d] = 0.;
      double sumW = 0.;
      DATA hh[D];
      for(DIM d=0;d<D;d++)
	{
	  double s=0., s2 =0.;
	  for(NUM n=0;n<N_nei;n++)
	    {
	      DATA x_nei = x[d][ id_nei[n] ];
	      s  += x_nei;
	      s2 += x_nei*x_nei;
	    }
	  s /= N_nei;
	  hh[d] = sqrt( std::max(0., s2/N_nei-s*s) );
	}
      for(NUM n=0;n<N_nei;n++)
	{
	  DATA xx[D];
	  for(DIM d=0;d<D;d++) xx[d] = x[d][ id_nei[n] ];
    
	  double w = 1.;
	  for(DIM d=0;d<D;d++)
	    {
	      double y = ( xx[d] - x[d][i] ) / hh[d];
	      w *= exp(-y*y/2)/ hh[d];
	    }
	  sumW += w;
    
	  for(DIM d=0;d<D;d++)
	    {
	      sum[d] += w*xx[d];
	      sum2[d] += w*xx[d]*xx[d];
	    }
	}
    
      if(sumW<=0.) sumW=1.;
      for(DIM d=0;d<D;d++)
	{
	  sum[d] /= sumW;
	  h[d][i] = sqrt( std::max(0., sum2[d]/sumW-sum[d]*sum[d]) );
	}

      DATA h0[D]; for(DIM d=0;d<D;d++) h0[d]=h[d][i];
      DATA x0[D]; for(DIM d=0;d<D;d++) x0[d]=x[d][i];
      double M, s, s0=min_s, s1=max_s;
      do
	{
	  s = sqrt(s0*s1);

	  DATA box0[2*D];
	  DATA *ptr = box0;
	  for(DIM d=0;d<D;d++){ *(ptr++)=x0[d]-s*h0[d]; *(ptr++)=x0[d]+s*h0[d]; }
	  //printf("box=");for(DIM d=0;d<D;d++)printf(" %.4e", box0[2*d]); printf("\n");
	  //printf("    ");for(DIM d=0;d<D;d++)printf(" %.4e",box0[2*d+1]); printf("\n");

	  node = 0;
	  Nr = Ntot;
	  id = 0;
	  M_found = 0.;
	  M = compute_FiEstAS_plain( box0, M0*(1.+EPS+EPS) );
	  //printf("\t s=%.3e M_F=%6.3f \n",s,M);

	  if(M>M0) s1=s;
	  else s0=s;
	}
      while( fabs(M-M0)>EPS*M0 && s1-s0>1e-8*s0 );

      for(DIM d=0;d<D;d++) h[d][i] = s*h0[d];

    }
  printf("\b\b\b\bDone! (t=%g s)\n", (clock()-t0)/double(CLOCKS_PER_SEC) );

  DATA *ptr = box; for(DIM d=0;d<D;d++){ *(ptr++)=min[d]; *(ptr++)=max[d]; }
  node = 0;
  Nr = Ntot;
  id = 0;
  set_h(box,box);
}
//------------------------------------------------------
void yTree::set_h(DATA *h_min, DATA *h_max)//, vector<double> &h0)
// compute hmin, hmax for every node (private)
//------------------------------------------------------
{
  PANIC( printf("Visit node %d (N=%d, Nl=%d) ",node,Nr,N_l[node]); );
  //printf("Visit node %d (N=%d, Nl=%d) \n",node,Nr,N_l[node]);
  //printf("  box ="); for(DIM d=0; d<D; d++) printf(" %.4e",box[2*d]); printf("\n");
  //printf("       "); for(DIM d=0; d<D; d++) printf(" %.4e",box[2*d+1]); printf("\n");

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("Particle %d is in! \n",id); );

      for(DIM d=0; d<D; d++)
	{
	  DATA h_d = h[d][id];

	  h[d][id] = h_d;

	  DATA x_d = x[d][id];
	  h_min[d] = x_d - h_d;
	  h_max[d] = x_d + h_d;
	}
    }
  else // ----------------------- Branch node
    {
      DATA lmin[D], rmin[D];
      DATA lmax[D], rmax[D];

      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];

      // if(ds==D) return;

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      // --------- explore left side
      PANIC( printf("%d -> L: ",node); );

      box[2*ds+1] = xs;
      Nr = Nl;
      node++;

      set_h(lmin,lmax);//,h0);

      box[2*ds+1] = old_box1;
      node = old_node;
      Nr = old_Nr;

      // --------- explore right side
      PANIC( printf("%d -> R: ",node); );

      box[2*ds] = xs;
      Nr -= Nl;
      node += Nl;
      id += Nl;

      set_h(rmin,rmax);//,h0);
      
      box[2*ds] = old_box0;
      node = old_node;
      Nr = old_Nr;
      id -= Nl;

      // --------- Compute hmin, hmax
      for(DIM d=0; d<D; d++)
	{
	  if(lmin[d]<rmin[d]) hmin[d][node] = h_min[d] = lmin[d];
	  else                hmin[d][node] = h_min[d] = rmin[d];
	  if(lmax[d]>rmax[d]) hmax[d][node] = h_max[d] = lmax[d];
	  else                hmax[d][node] = h_max[d] = rmax[d];
	}
    }

}

//------------------------------------------------------
double yTree::set_metric(vector<NUM> &dim, vector<double> &scale)
//------------------------------------------------------
{
  
  printf(" Adjust metric \n");
  NUM N_dim = dim.size();
  ERROR( scale.size()!=N_dim, ("N_dim=%d, N_scale=%d",N_dim,(int)scale.size()) );

  printf("  d scale \n");
  for(NUM n=0; n<N_dim; n++) printf("%3d %g \n", dim[n], scale[n]);

  double S = 1.;
  for(NUM n=0; n<N_dim; n++) S *= scale[n];

////////////////////////////
  
  DATA new_h[D][Ntot];
  double max_change = 0.;
  clock_t t0 = clock();
  printf(" Adjusting...   %%"); fflush(stdout);

  for(NUM i=0, di=0; i<Ntot; i++)
    {
      if(i>di)
	{
	  printf("\b\b\b\b%3d%%",(100*i)/Ntot); fflush(stdout);
	  di+=Ntot/100;
	}
  
      double V = 1.;
      for(NUM n=0; n<N_dim; n++) V *= h[ dim[n] ][i];
      V = pow( V/S, 1./N_dim );

      for(DIM d=0; d<D; d++) new_h[d][i] = h[d][i];
      for(NUM n=0; n<N_dim; n++)
	{
	  new_h[ dim[n] ][i] = V * scale[n];
	  double eps = fabs( new_h[ dim[n] ][i] / h[ dim[n] ][i] -1. );
	  if( eps > max_change ) max_change = eps;
	}
    }
  printf("\b\b\b\bDone! (t=%g s) max_change=%g\n", (clock()-t0)/double(CLOCKS_PER_SEC), max_change );

  for(DIM d=0;d<D;d++) for(NUM i=0; i<Ntot; i++) h[d][i] = new_h[d][i];

  DATA *ptr = box; for(DIM d=0;d<D;d++){ *(ptr++)=min[d]; *(ptr++)=max[d]; }
  node = 0;
  Nr = Ntot;
  id = 0;
  set_h(box,box);
/*

  for(NUM i=0, di=0; i<Ntot; i++)
    {
      if(i>di)
	{
	  printf("\b\b\b\b%3d%%",(100*i)/Ntot); fflush(stdout);
	  di+=Ntot/100;
	}

      DATA h0[D]; for(DIM d=0;d<D;d++) h0[d]=h[d][i];
      DATA x0[D]; for(DIM d=0;d<D;d++) x0[d]=x[d][i];
      double M, s, s0=min_s, s1=max_s;
      do
	{
	  s = sqrt(s0*s1);

	  DATA box0[2*D];
	  for(DIM d=0;d<D;d++){ box0[2*d]=x0[d]-s*h0[d]; box0[2*d+1]=x0[d]+s*h0[d]; }
	  //printf("box=");for(DIM d=0;d<D;d++)printf(" %.4e", box0[2*d]); printf("\n");
	  //printf("    ");for(DIM d=0;d<D;d++)printf(" %.4e",box0[2*d+1]);printf("\n");

	  node = 0;
	  Nr = Ntot;
	  id = 0;
	  M_found = 0.;
	  M = compute_FiEstAS_plain( box0, M0*(1.+EPS+EPS) );
	  //printf("\t s=%.3e M_F=%6.3f \n",s,M);

	  if(M>M0) s1=s;
	  else s0=s;
	}
      while( fabs(M-M0)>EPS*M0 && s1-s0>1e-8*s0 );

      for(DIM d=0;d<D;d++) new_h[d][i] = s*h0[d];
    }
  printf("\b\b\b\bDone! (t=%g s) max_change=%g\n", (clock()-t0)/double(CLOCKS_PER_SEC), max_change );

  for(DIM d=0;d<D;d++) for(NUM i=0; i<Ntot; i++) h[d][i] = new_h[d][i];

  ptr = box; for(DIM d=0;d<D;d++){ *(ptr++)=min[d]; *(ptr++)=max[d]; }
  node = 0;
  Nr = Ntot;
  id = 0;
  set_h(box,box);
*/
  return(max_change);
}

//##############################################################################

//------------------------------------------------------
double yTree::smooth(DATA *x0, bool subtract_bias)
// return scatter kernel density estimate (public)
//------------------------------------------------------
{
  ERROR( hmin==NULL, ("hsmooth not set!") );

  DATA *ptr = box; for(DIM d=0;d<D;d++){ *(ptr++)=min[d]; *(ptr++)=max[d]; }
  node = 0;
  Nr = Ntot;
  id = 0;

  double result = compute_smooth(x0);
  // if( fsmooth!=NULL ) result = (*fsmooth)[Ntot]/result;

  result /= Ntot;
  if( subtract_bias ) result /= bias;
  
  return( result );
}
//------------------------------------------------------
double yTree::compute_smooth(DATA *x0)
//------------------------------------------------------
{
  double result=0.;

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("particle %d \n",id); );

      double Wi =1.;
      for(DIM d=0; d<D; d++)
	{
	  //DATA min = box[2*d];
	  //DATA max = box[2*d+1];
	  DATA hi = h[d][id];
	  DATA xi = x[d][id];
	  Wi *= Kernel->Kernel( (x0[d]-xi)/hi ) /hi;
	}

      // if(fsmooth!=NULL) *(fsmooth)[Ntot] += Wi * (*fsmooth)[id];
      result = Wi;

      //printf("point %d: x =",id); for(DIM d=0; d<D; d++) printf(" %6.3f", x[d][id]);
      //printf(", h ="); for(DIM d=0; d<D; d++) printf(" %6.3f", h[d][id]);
      //printf(", M = %.3e \n",Wi);
    }
  else // ----------------------- Branch node
    {
      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];
      PANIC( printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id); );
      //printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id);
      //printf("\t hmin ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmin[d][node]); printf("\n");
      //printf("\t hmax ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmax[d][node]); printf("\n");

      // if(ds==D) return;

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      if( Nl==1 || x0[ds]<hmax[ds][node+1] ) // explore left
	{
	  PANIC( printf("%d -> L: ",node); );
	  //printf("%d -> L: ",node);

	  box[2*ds+1] = xs;
	  Nr = Nl;
	  node++;

	  result += compute_smooth(x0);

	  box[2*ds+1] = old_box1;
	  node = old_node;
	  Nr = old_Nr;
	}

      if( Nr-Nl==1 || x0[ds]>hmin[ds][node+Nl] ) // explore right
	{
	  PANIC( printf("%d -> R: ",node); );
	  //printf("%d -> R: ",node);

	  box[2*ds] = xs;
	  Nr -= Nl;
	  node += Nl;
	  id += Nl;

	  result += compute_smooth(x0);
      
	  box[2*ds] = old_box0;
	  node = old_node;
	  Nr = old_Nr;
	  id -= Nl;
	}
    }

  return(result);
}

//##############################################################################

//------------------------------------------------------
double yTree::get_hsmooth(DATA *x0, DATA *h0)
// returns smoothing box at an arbitrary point x0 (public)
//------------------------------------------------------
{
  ERROR( hmin==NULL, ("hsmooth not set!") );

  for(DIM d=0;d<D;d++) h0[d] = 0.;
  DATA *ptr = box; for(DIM d=0;d<D;d++){ *(ptr++)=min[d]; *(ptr++)=max[d]; }
  node = 0;
  Nr = Ntot;
  id = 0;
  double result = compute_hsmooth(x0,h0);
  if(result>0.) for(DIM d=0;d<D;d++) h0[d] /= result;
  else return(0.);

  return(result);
}
//------------------------------------------------------
double yTree::compute_hsmooth(DATA *x0, DATA *h0)
//------------------------------------------------------
{
  double result=0.;

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("particle %d \n",id); );

      double Wi =1.;
      for(DIM d=0; d<D; d++)
	{
	  DATA hi = h[d][id];
	  DATA xi = x[d][id];
	  Wi *= Kernel->Kernel( (x0[d]-xi)/hi ) /hi;
	}

      for(DIM d=0; d<D; d++) h0[d] += Wi * h[d][id];
      result = Wi;

      //printf("point %d: x =",id); for(DIM d=0; d<D; d++) printf(" %6.3f", x[d][id]);
      //printf(", h ="); for(DIM d=0; d<D; d++) printf(" %6.3f", h[d][id]);
      //printf(", M = %.3e \n",Wi);
    }
  else // ----------------------- Branch node
    {
      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];
      PANIC( printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id); );
      //printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id);
      //printf("\t hmin ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmin[d][node]); printf("\n");
      //printf("\t hmax ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmax[d][node]); printf("\n");

      // if(ds==D) return;

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      if( Nl==1 || x0[ds]<hmax[ds][node+1] ) // explore left
	{
	  PANIC( printf("%d -> L: ",node); );
	  //printf("%d -> L: ",node);

	  box[2*ds+1] = xs;
	  Nr = Nl;
	  node++;

	  result += compute_hsmooth(x0,h0);

	  box[2*ds+1] = old_box1;
	  node = old_node;
	  Nr = old_Nr;
	}

      if( Nr-Nl==1 || x0[ds]>hmin[ds][node+Nl] ) // explore right
	{
	  PANIC( printf("%d -> R: ",node); );
	  //printf("%d -> R: ",node);

	  box[2*ds] = xs;
	  Nr -= Nl;
	  node += Nl;
	  id += Nl;

	  result += compute_hsmooth(x0,h0);
      
	  box[2*ds] = old_box0;
	  node = old_node;
	  Nr = old_Nr;
	  id -= Nl;
	}
    }

  return(result);
}

//##############################################################################

//------------------------------------------------------
double yTree::FiEstAS(DATA *x0, bool balloon, bool subtract_bias)
// Return balloon estimate (public)
//------------------------------------------------------
{
  if( balloon==false ) return( smooth(x0,subtract_bias) );
  
  ERROR( hmin==NULL, ("hsmooth not set!") );

  DATA h0[D];
  get_hsmooth(x0,h0);
  if(h0[0]<=0.) return(0.);

  DATA box0[2*D];
  for(DIM d=0;d<D;d++)
  {
    box0[2*d  ] = x0[d] - h0[d];
    box0[2*d+1] = x0[d] + h0[d];
  }
  node = 0;
  Nr = Ntot;
  id = 0;
  M_found = 0.;
  double rho = compute_FiEstAS( box0, 1e30 );
  for(DIM d=0; d<D; d++) rho /= ( box0[2*d+1] - box0[2*d] );

  if( subtract_bias ) rho /= 1.+1./M0;
  
  return( rho / Ntot );
}

//------------------------------------------------------
double yTree::compute_FiEstAS(DATA *box0, double Mmax)
// returns mass within box?
//------------------------------------------------------
{
  double result=0.;

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("particle %d \n",id); );

      double Mi =1.;
      for(DIM d=0; d<D; d++)
	{
	  DATA hi = h[d][id];
	  DATA xi = x[d][id];
	  DATA b0 = box0[2*d];
	  DATA b1 = box0[2*d+1];
	  Mi *= Kernel->Kernel( (b0-xi)/hi, (b1-xi)/hi ); // /( b1 - b0 );
	}

      M_found += Mi;
      for(DIM d=0; d<D; d++) h_found[d] += Mi * h[d][id];

      // if(fsmooth!=NULL) *(fsmooth)[Ntot] += Wi * (*fsmooth)[id];
      result = Mi;

      //printf("point %d: x =",id); for(DIM d=0; d<D; d++) printf(" %6.3f", x[d][id]);
      //printf(", h ="); for(DIM d=0; d<D; d++) printf(" %6.3f", h[d][id]);
      //printf(", M = %.3e \n",Mi);
    }
  else // ----------------------- Branch node
    {
      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];
      PANIC( printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id); );
      //printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id);
      //printf("\t hmin ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmin[d][node]); printf("\n");
      //printf("\t hmax ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmax[d][node]); printf("\n");

      // if(ds==D) return;

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      if(
	 M_found<Mmax &&
	 ( Nl==1 || box0[2*ds]<hmax[ds][node+1] )
	 ) // explore left
	{
	  PANIC( printf("%d -> L: ",node); );
	  //printf("%d -> L: ",node);

	  box[2*ds+1] = xs;
	  Nr = Nl;
	  node++;

	  result += compute_FiEstAS(box0, Mmax);

	  box[2*ds+1] = old_box1;
	  node = old_node;
	  Nr = old_Nr;
	}
      //else printf("%d -> L out: (%g)>(%g)\n",node, box0[2*ds],hmax[ds][node+1]);

      if(
	 M_found<Mmax &&
	 ( Nr-Nl==1 || box0[2*ds+1]>hmin[ds][node+Nl])
	 ) // explore right
	{
	  PANIC( printf("%d -> R: ",node); );
	  //printf("%d -> R: ",node);

	  box[2*ds] = xs;
	  Nr -= Nl;
	  node += Nl;
	  id += Nl;

	  result += compute_FiEstAS(box0, Mmax);
      
	  box[2*ds] = old_box0;
	  node = old_node;
	  Nr = old_Nr;
	  id -= Nl;
	}
      //else printf("%d -> R out: (%g)<(%g)\n",node, box0[2*ds+1],hmin[ds][node+Nl]);

    }

  return(result);
}
//------------------------------------------------------
double yTree::compute_FiEstAS_plain(DATA *box0, double Mmax)
// returns mass within box?
//------------------------------------------------------
{
  double result=0.;

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("particle %d \n",id); );

      double Mi =1.;
      for(DIM d=0; d<D; d++)
	{
	  DATA min_i = box[2*d];
	  DATA max_i = box[2*d+1];
	  DATA l_i = max_i-min_i;
	  DATA min_box0 = box0[2*d];
	  DATA max_box0 = box0[2*d+1];
	  if(min_box0>max_i){ Mi=0.; break; }
	  if(max_box0<min_i){ Mi=0.; break; }
	  if(l_i>0.)
	    Mi *= ( std::min(max_box0,max_i) - std::max(min_box0,min_i) ) / l_i;
	  //else printf("TODO!\n");
	}

      M_found += Mi;

      // if(fsmooth!=NULL) *(fsmooth)[Ntot] += Wi * (*fsmooth)[id];
      result = Mi;

      //printf("point %d: x =",id); for(DIM d=0; d<D; d++) printf(" %6.3f", x[d][id]);
      //printf(", h ="); for(DIM d=0; d<D; d++) printf(" %6.3f", h[d][id]);
      //printf(", M = %.3e \n",Mi);
    }
  else // ----------------------- Branch node
    {
      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];
      PANIC( printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id); );
      //printf("Visit node %d (N=%d, Nl=%d, id=%d) \n",node,Nr,Nl,id);
      //printf("\t hmin ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmin[d][node]); printf("\n");
      //printf("\t hmax ="); for(DIM d=0; d<D; d++) printf(" %.4e",hmax[d][node]); printf("\n");

      // if(ds==D) return;

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      if( M_found<Mmax && box0[2*ds]<xs ) // explore left
	{
	  PANIC( printf("%d -> L: ",node); );
	  //printf("%d -> L: ",node);

	  box[2*ds+1] = xs;
	  Nr = Nl;
	  node++;

	  result += compute_FiEstAS_plain(box0, Mmax);

	  box[2*ds+1] = old_box1;
	  node = old_node;
	  Nr = old_Nr;
	}
      //else printf("%d -> L out: (%g)>(%g)\n",node, box0[2*ds],hmax[ds][node+1]);

      if( M_found<Mmax && box0[2*ds+1]>xs ) // explore right
	{
	  PANIC( printf("%d -> R: ",node); );
	  //printf("%d -> R: ",node);

	  box[2*ds] = xs;
	  Nr -= Nl;
	  node += Nl;
	  id += Nl;

	  result += compute_FiEstAS_plain(box0, Mmax);
      
	  box[2*ds] = old_box0;
	  node = old_node;
	  Nr = old_Nr;
	  id -= Nl;
	}
      //else printf("%d -> R out: (%g)<(%g)\n",node, box0[2*ds+1],hmin[ds][node+Nl]);

    }

  return(result);
}

//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
