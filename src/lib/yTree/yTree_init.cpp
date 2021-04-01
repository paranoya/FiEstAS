#include"FiEstAS.h"

//----------------------------------------------------------------------
yTree::yTree(DIM d, bool quiet):silent(quiet)
//----------------------------------------------------------------------
{
  PANIC( printf("yTree called\n"); );

  null_data();
  null_extra();
  null_scratch1();
  null_scratch2();

  null_smooth();

  init_data(d);
  init_extra(0);
}

//----------------------------------------------------------------------
yTree::~yTree()
//----------------------------------------------------------------------
{
  PANIC( printf("~yTree called\n"); );

  delete_data();
  delete_extra();
  delete_scratch1();
  delete_scratch2();
  delete_smooth();
}

//----------------------------------------------------------------------
void yTree::init_data(DIM d)
//----------------------------------------------------------------------
{
  PANIC( printf("init_data called (D=%d)\n",d); );

  delete_data();
  if(d>0)
    {
      D = d;
      x = new vector<DATA>[D];
      min = new DATA[D];
      max = new DATA[D];
      Nsplits = new NUM[D];
      for(DIM dd=0;dd<D;dd++)
	{
	  min[dd]= yHUGE;
	  max[dd]=-yHUGE;
	  Nsplits[dd]=0;
	}
    }
}

//----------------------------------------------------------------------
void yTree::delete_data()
//----------------------------------------------------------------------
{
  PANIC( printf("delete_data called\n"); );

  if(x      !=NULL) delete[] x;
  if(min    !=NULL) delete[] min;
  if(max    !=NULL) delete[] max;
  if(Nsplits!=NULL) delete[] Nsplits;

  null_data();
}
//----------------------------------------------------------------------
void yTree::null_data()
//----------------------------------------------------------------------
{
  x       = NULL;
  min     = NULL;
  max     = NULL;
  Nsplits = NULL;
  D = 0;
  Ntot = 0;
  tilt_left = true;
}

//----------------------------------------------------------------------
void yTree::init_extra(DIM d)
//----------------------------------------------------------------------
{
  PANIC( printf("init_extra called (d=%d)\n",d); );

  delete_extra();
  if(d>0)
    {
      D_extra = d;
      extra = new vector<DATA>[D_extra];
    }
}

//----------------------------------------------------------------------
void yTree::delete_extra()
//----------------------------------------------------------------------
{
  PANIC( printf("delete_extra called\n"); );

  if(extra!=NULL) delete[] extra;
  null_extra();
}
//----------------------------------------------------------------------
void yTree::null_extra()
//----------------------------------------------------------------------
{
  extra = NULL;
  D_extra = 0;
}

//----------------------------------------------------------------------
void yTree::init_scratch1(NUM B)
//----------------------------------------------------------------------
{
  PANIC( printf("init_scratch1 called (B=%d)\n",B); );

  delete_scratch1();
  DEBUG_ERROR(B<=0,("> IMPOSSIBLE! (B=%d)",B));

  Nb = new NUM[B];
  xmin = new DATA[B];
  xmax = new DATA[B];

  d_s = new DIM[Ntot];
  x_s = new DATA[Ntot];
  N_l = new NUM[Ntot];
}

//----------------------------------------------------------------------
void yTree::delete_scratch1()
//----------------------------------------------------------------------
{
  PANIC( printf("delete_scratch1 called\n"); );

  if(Nb  !=NULL) delete[] Nb;
  if(xmin!=NULL) delete[] xmin;
  if(xmax!=NULL) delete[] xmax;
  if(d_s !=NULL) delete[] d_s;
  if(x_s !=NULL) delete[] x_s;
  if(N_l !=NULL) delete[] N_l;

  null_scratch1();
}
//----------------------------------------------------------------------
void yTree::null_scratch1()
//----------------------------------------------------------------------
{
  Nb   = NULL;
  xmin = NULL;
  xmax = NULL;
  d_s  = NULL;
  x_s  = NULL;
  N_l  = NULL;
}

//----------------------------------------------------------------------
void yTree::init_scratch2()
//----------------------------------------------------------------------
{
  PANIC( printf("init_scratch2 called\n"); );

  delete_scratch2();
  box = new DATA[2*D];
}

//----------------------------------------------------------------------
void yTree::delete_scratch2()
//----------------------------------------------------------------------
{
  PANIC( printf("delete_scratch2 called\n"); );

  if(box!=NULL) delete[] box;
  null_scratch2();
}
//----------------------------------------------------------------------
void yTree::null_scratch2()
//----------------------------------------------------------------------
{
  box = NULL;
}

//----------------------------------------------------------------------
void yTree::init_smooth()
//----------------------------------------------------------------------
{
  PANIC( printf("init_smooth called\n"); );

  delete_smooth();
  /*
  hmin = new vector<DATA>[D](Ntot);
  hmax = new vector<DATA>[D](Ntot);
  h    = new vector<DATA>[D](Ntot);
  */
  hmin = new vector<DATA>[D];
  hmax = new vector<DATA>[D];
  h    = new vector<DATA>[D];
  for(DIM d=0; d<D; d++)
    for(NUM i=0; i<Ntot; i++)
      {
	hmin[d].push_back(0);
	hmax[d].push_back(0);
	h   [d].push_back(0);
      }

  h_found = new double[D];
}

//----------------------------------------------------------------------
void yTree::delete_smooth()
//----------------------------------------------------------------------
{
  PANIC( printf("delete_smooth called\n"); );

  if(hmin!=NULL) delete[] hmin;
  if(hmax!=NULL) delete[] hmax;
  if(h   !=NULL) delete[] h;
  if(h_found!=NULL) delete[] h_found;
  null_smooth();
}
//----------------------------------------------------------------------
void yTree::null_smooth()
//----------------------------------------------------------------------
{
  hmin = NULL;
  hmax = NULL;
  h    = NULL;
  h_found = NULL;
}

//----------------------------------------------------------------------
void yTree::add_xd(DIM d, DATA xd)
//----------------------------------------------------------------------
{
  DEBUG_ERROR(d>=D, ("> Wrong number of dimensions (%d>=%d)",d,D) );

  x[d].push_back(xd);
  if(xd<min[d]) min[d]=xd;
  if(xd>max[d]) max[d]=xd;
}

//----------------------------------------------------------------------
void yTree::add_extra(DIM d, DATA xd)
//----------------------------------------------------------------------
{
  DEBUG_ERROR(d>=D, ("> Wrong number of dimensions (%d>=%d)",d,D_extra) );

  extra[d].push_back(xd);
}

//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
