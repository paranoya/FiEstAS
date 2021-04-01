#include"FiEstAS.h"

#include<limits.h>

#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

//----------------------------------------------------------------------
void yTree::read_ASCII(const char * name)
//----------------------------------------------------------------------
{
  ifstream file(name); // Open file
  ERROR(!file.good(),("File '%s' not found", name));

  while( file.peek()=='#' ) file.ignore(INT_MAX,'\n'); // Ignore comments

  string str;
  getline(file,str);
  istringstream s(str);
  DIM d=0;
  DATA xd;
  while(s.peek()!=EOF){ s>>xd; d++; }
  yINFO((" File '%s' has %d columns \n",name,d));
  ERROR(d==0,(" No data could be read from '%s' \n",name));

  if(D>0 && D<d)
    {
      yINFO((" (only the first %d will be used by yTree)\n",D));
    }
  else
    {
      yWARNING(D>d,(" D(=%d) > d(=%d) => reducing D ",D,d));
      D = d;
    }

  init_data(D);
  init_extra(d-D);

  do
    {
      s.clear();
      s.str(str);
      for(d=0; d<D; d++){ s>>xd; add_xd(d,xd); }
      for(d=0; d<D_extra; d++){ s>>xd; add_extra(d,xd); }
    }
  while(getline(file,str));
  file.close();
}

//----------------------------------------------------------------------
void yTree::copy_data_to(yTree *tree)
//----------------------------------------------------------------------
{
  for(DIM d=0; d<D; d++)
    {
      for(vector<DATA>::iterator i=x[d].begin(); i<x[d].end(); i++)
	tree->add_xd(d,*i);
    }
}

//----------------------------------------------------------------------
void yTree::copy_extra_to(yTree *tree)
//----------------------------------------------------------------------
{
  for(DIM d=0; d<D_extra; d++) copy_extra_d_to(d,tree,d);
}
//----------------------------------------------------------------------
void yTree::copy_extra_d_to(DIM d1, yTree *tree,DIM d2)
//----------------------------------------------------------------------
{
  for(vector<DATA>::iterator i=extra[d1].begin(); i<extra[d1].end(); i++)
    tree->add_extra(d2,*i);
}

//----------------------------------------------------------------------
void yTree::write_data(const char *fname)
//----------------------------------------------------------------------
{
  ofstream file(fname);
  for(NUM i=0; i<Ntot; i++)
    {
      for(DIM d=0; d<D; d++) file <<' '<< x[d][i];
      file << endl;
    }
  file.close();
}

//----------------------------------------------------------------------
void yTree::write_extra(const char *fname)
//----------------------------------------------------------------------
{
  ofstream file(fname);
  for(NUM i=0; i<Ntot; i++)
    {
      for(DIM d=0; d<D_extra; d++) file <<' '<< extra[d][i];
      file << endl;
    }
  file.close();
}

//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
