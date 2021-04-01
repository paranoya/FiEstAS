#include"FiEstAS.h"

//------------------------------------------------------
DATA yTree::get_min(DIM d)
//------------------------------------------------------
{
  DEBUG_ERROR(d>=D, ("> Out of bounds in getMin: (d=%d) >= (D=%d)",d,D) );
  return( min[d] );
}

//------------------------------------------------------
DATA yTree::get_max(DIM d)
//------------------------------------------------------
{
  DEBUG_ERROR(d>=D, ("> Out of bounds in getMax: (d=%d) >= (D=%d)",d,D) );
  return( max[d] );
}

//------------------------------------------------------
void yTree::get_Bbox(DATA* b)
//------------------------------------------------------
{
  for(DIM d=0; d<D; d++){ b[2*d]=min[d];  b[2*d+1]=max[d]; }
}

//------------------------------------------------------
NUM yTree::get_n(DATA *x0)
//------------------------------------------------------
{
  node = 0;
  Nl = 0;
  Nr = Ntot;
  if(Ntot>1) compute_n(x0);

  return( Nl );
}
//------------------------------------------------------
void yTree::compute_n(DATA *x0)
//------------------------------------------------------
{
  ds = d_s[node];
  xs = x_s[node];

  PANIC( printf("Visit node %d (N=%d, Nl=%d) \n",node,Nr,N_l[node]); );

  if(ds==D) return;

  if( x0[ds]<=xs ) // always tilt left
    {
      if( (Nr=N_l[node]) >1)
	{
	  node++;
	  compute_n(x0);
	}
    }
  else
    {
      Nl += N_l[node];
      if( (Nr-=N_l[node]) >1)
	{
	  node += N_l[node];
	  compute_n(x0);
	}
    }
}


//------------------------------------------------------
DATA yTree::get_data(DIM d, NUM n)
//------------------------------------------------------
{
  DEBUG_ERROR(d>=D, ("> Out of bounds in get_data: (d=%d) >= (D=%d)",d,D) );
  DEBUG_ERROR(n>=x[d].size(),
	      ("> Out of bounds in x[%d] (n=%d) >= %d", d,n,(int)x[d].size())
	      );
  return( x[d][n] );
}

//------------------------------------------------------
DATA yTree::get_extra(DIM d, NUM n)
//------------------------------------------------------
{
  DEBUG_ERROR(d>=D_extra, ("> Out of bounds in get_extra: (d=%d) >= (D=%d)",d,D_extra) );
  DEBUG_ERROR(n>=extra[d].size(),
	      ("> Out of bounds in extra[%d] (n=%d) >= %d", d,n,(int)extra[d].size())
	      );
  return( extra[d][n] );
}

//------------------------------------------------------
void yTree::push_box(DATA *x0, vector<DATA> *box_out)
//------------------------------------------------------
{
  if(box!=NULL)
    {
      DATA *ptr = box;
      for(DIM d=0;d<D;d++){ *(ptr++)=min[d];  *(ptr++)=max[d]; }
      node = 0;
      Nr = Ntot;
      if(Ntot>1) compute_box(x0);

      for(DIM i=0; i<2*D; i++) box_out->push_back( box[i] );
    }
}
//------------------------------------------------------
void yTree::compute_box(DATA *x0)
//------------------------------------------------------
{
  ds = d_s[node];
  xs = x_s[node];

  PANIC( printf("Visit node %d (N=%d, Nl=%d)\n",node,Nr,N_l[node]); );

  if(ds==D) return;

  if( x0[ds]<=xs ) // always tilt left
    {
      box[2*ds+1] = xs;
      if( (Nr=N_l[node]) >1)
	{
	  node++;
	  compute_box(x0);
	}
    }
  else
    {
      box[2*ds] = xs;
      if( (Nr-=N_l[node]) >1)
	{
	  node += N_l[node];
	  compute_box(x0);
	}
    }
}

//------------------------------------------------------
void yTree::push_box(NUM n, vector<DATA> *box_out)
//------------------------------------------------------
{
  PANIC( printf("push_box called (n=%d) \n",n); );

  if(box!=NULL)
    {
      DATA *ptr = box;
      for(DIM d=0;d<D;d++){ *(ptr++)=min[d];  *(ptr++)=max[d]; }
      node = 0;
      Nr = Ntot;
      if(Ntot>1) compute_box(n+1);
    }

  PANIC(
	printf(" box(n=%d):",n);
	for(DIM d=0; d<D; d++) printf(" [%g,%g]",box[2*d],box[2*d+1]);
	printf("\n");
	);

  for(DIM i=0; i<2*D; i++) box_out->push_back( box[i] );
}

//------------------------------------------------------
void yTree::compute_box(NUM n)
//------------------------------------------------------
{
  ds = d_s[node];
  xs = x_s[node];

  PANIC( printf("Visit node %d (N=%d, Nl=%d) n=%d ",node,Nr,N_l[node],n); );

  if(ds==D) return;

  NUM Nl = N_l[node];
  if( n > Nl )
    {
      PANIC( printf("go right\n"); );

      box[2*ds] = xs;
      n -= Nl;
      Nr-= Nl;
      if( Nr > 1 )
	{
	  node += Nl;
	  compute_box(n);
	}
    }
  else
    {
      PANIC( printf("go left\n"); );

      box[2*ds+1] = xs;
      Nr = Nl;
      if( Nr > 1 )
	{
	  node++;
	  compute_box(n);
	}
    }
}

//------------------------------------------------------
void yTree::get_split(NUM n, DIM &d, DATA &x)
//------------------------------------------------------
{
  ERROR(n>Ntot-1,("n=%d is out of bounds! (N=%d)",n,Ntot));
  // CHECK ERROR CHECK !

  d = d_s[n];
  x = x_s[n];
}

//------------------------------------------------------
NUM yTree::push_cells(vector<DATA> *cells, vector<NUM> *points, bool saveCells)
//------------------------------------------------------
{
  PANIC( printf("get_cells called \n"); );
  PANIC(
	printf(" search box:");
	for(DIM d=0; d<D; d++) printf(" [%g,%g]",cells[2*d][0],cells[2*d+1][0]);
	printf("\n");
	);

  NUM n = 0;

  if(box!=NULL)
    {
      DATA *ptr = box;
      for(DIM d=0;d<D;d++){ *(ptr++)=min[d];  *(ptr++)=max[d]; }
      node = 0;
      Nr = Ntot;
      id = 0;
      if(Ntot>1) n = compute_cells(cells,points,saveCells);
    }

  return(n);
}

//------------------------------------------------------
NUM yTree::compute_cells(vector<DATA> *cells, vector<NUM> *points, bool saveCells)
//------------------------------------------------------
{
  NUM n = 0;

  if(Nr==1) // ----------------------- Leaf node
    {
      PANIC( printf("Particle %d is in! \n",id); );
      //printf("Particle %d is in! \n",id);
      n = 1;

      if(saveCells)
	for(DIM d=0; d<D; d++)
	  {
	    cells[2*d  ].push_back( box[2*d  ] );
	    cells[2*d+1].push_back( box[2*d+1] );
	  }
      //printf("l="); for(DIM d=0; d<D; d++) printf(" %.3e",box[2*d  ]); printf("\n");
      //printf("x="); for(DIM d=0; d<D; d++) printf(" %.3e",  x[d][id]); printf("\n");
      //printf("r="); for(DIM d=0; d<D; d++) printf(" %.3e",box[2*d+1]); printf("\n");

      if(points!=NULL) points->push_back(id);
    }
  else // ----------------------- Branch node
    {

      DIM ds = d_s[node];
      DATA xs = x_s[node];
      NUM Nl = N_l[node];

      PANIC( printf("Visit node %d (N=%d, Nl=%d, ds=%d, xs=%g) \n",node,Nr,Nl,ds,xs); );
      //printf("Visit node %d (N=%d, Nl=%d, ds=%d, xs=%g) \n",node,Nr,Nl,ds,xs);

      if(ds==D)
	{
	  yWARNING(true,("Ignoring %d points at the same place!",Nr));
	  return(0);
	}

      DATA old_box0 =box[2*ds  ];
      DATA old_box1 =box[2*ds+1];
      NUM old_Nr = Nr;
      NUM old_node = node;

      if( cells[2*ds][0] <= xs ) // explore left side
	{
	  PANIC( printf("%d -> L: ",node); );
	  box[2*ds+1] = xs;
	  Nr = Nl;
	  node++;

	  n += compute_cells(cells,points,saveCells);

	  box[2*ds+1] = old_box1;
	  node = old_node;
	  Nr = old_Nr;
	}

      if( cells[2*ds+1][0] >= xs ) // explore right side
	{
	  PANIC( printf("%d -> R: ",node); );
	  box[2*ds] = xs;
	  Nr -= Nl;
	  node += Nl;
	  id += Nl;

	  n += compute_cells(cells,points,saveCells);

	  box[2*ds] = old_box0;
	  node = old_node;
	  Nr = old_Nr;
	  id -= Nl;
	}

    }

  return(n);
}


//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
