#include"FiEstAS.h"

#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

//----------------------------------------------------------------------
void yTree::read_Gadget(const char *name, vector<yGadget::Data> data, yGadget::PartType part_type)
//----------------------------------------------------------------------
{
  yGadget snap(name);
  snap.printInfo();

  NUM Ndata = data.size();

  D=0;
  for(NUM i=0; i<Ndata; i++) D+=snap.get_dim(data[i]);
  init_data(D);
  init_extra(0);

  NUM N = snap.Npart( part_type );
  NUM Nread = 0;
  static const NUM Npage = 10000;
  float points[D*Npage];
  while(Nread<N)
    {

      NUM Nmin = Npage;
      NUM offset[Ndata], off=0;
      for(NUM i=0; i<Ndata; i++)
	{
	  snap.read( data[i], part_type, Nread );
	  NUM n = snap.next( (char *)(points+off), Npage );
	  offset[i]=off;
	  off += n*snap.get_dim(data[i]);
	  if(n<Nmin) Nmin=n;
	}
      Nread += Nmin;

      DIM d=0;
      for(NUM j=0; j<Ndata; j++)
	{
	  DIM Dj=snap.get_dim(data[j]);
	  for(DIM dj=0; dj<Dj; dj++)
	    {
	      for(NUM i=0; i<Nmin; i++)
		add_xd( d, points[offset[j]+Dj*i+dj] );
	      d++;
	    }
	}

    }
}

//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
