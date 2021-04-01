#include"FiEstAS.h"

#include<limits.h>

#include<iostream>
#include<fstream>
#include<sstream>
#include<math.h>
using namespace std;

//----------------------------------------------------------------------
int yASCII_data::read(const char *name)
//----------------------------------------------------------------------
{
  if(Nrows>0) reset();

  // --- Open file

  ifstream file(name);
  if(file.good()==false)
    {
      cout << "File '"<<name<<"' not found\n";
      reset();
      return(0);
    }
  file_name = name;

  while( file.peek()=='#' ) // Ignore comments
    {
      Ncomments++;
      file.ignore(INT_MAX,'\n');
    }

  // --- Find number of columns

  string str;
  getline(file,str);
  istringstream s(str);
  double xd;
  //cout<<"str='"<<str<<"'\n";
  while(!(s>>xd).fail()) Ncolumns++; // cout<<" '"<<xd<<"'"; } cout<<endl;
  if(Ncolumns==0)
    {
      cout << "File '"<<name<<"' has 0 data columns\n";
      reset();
      return(0);
    }

  // --- Initialize

  column = new vector<double>[Ncolumns];
  max = new double[Ncolumns];
  min = new double[Ncolumns];
  ave = new double[Ncolumns];
  sig = new double[Ncolumns];
  for(int c=0; c<Ncolumns; c++)
    {
      max[c] = -yHUGE;
      min[c] =  yHUGE;
      ave[c] = 0.;
      sig[c] = 0.;
    }

  // --- Read data

  do
    {
      s.clear();
      s.str(str);
      for(int c=0; c<Ncolumns; c++)
	{
	  s>>xd;
	  column[c].push_back(xd);
	  if(xd>max[c]) max[c]=xd;
	  if(xd<min[c]) min[c]=xd;
	  ave[c] += xd;
	  sig[c] += xd*xd;
	}
      Nrows++;
    }
  while(getline(file,str));
  file.close();

  for(int c=0; c<Ncolumns; c++)
    {
      ave[c] /= Nrows;
      sig[c] = sqrt( sig[c]/Nrows - ave[c]*ave[c] );
    }
  cout << Nrows<<" rows read from '"<<name<<"' ("<<Ncolumns<<" columns)\n";
  return(Nrows);
}

//---------------------------------------------------------------------
