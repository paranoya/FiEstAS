#ifndef FIESTAS_H
#error "#---------------------------------------------------------#"
#error "|                                                         |"
#error "| ERROR: ySnapshot.h should only be included by FiEstAS.h |"
#error "|                                                         |"
#error "#---------------------------------------------------------#"
#else

#include <string.h>

#include <iostream>
#include <fstream>
#include<vector>
using namespace std;

//----------------------------------------------------------------------
class ySnapshot
//----------------------------------------------------------------------
{
 public:

  ySnapshot(){ init(); }
  ySnapshot(const char * name){ open(name); }
  virtual ~ySnapshot(){ close(); }

  const char *filename;
  enum PartType{ no_type=0, all_types, gas, dm, dm1,dm2,dm3,dm4,dm5, stars };
  enum Data{ no_data=0, pos,x,y,z, vel,vx,vy,vz };

  static PartType get_ptype(const char *str)
    {
      if(strcmp("gas",str)==0) return(ySnapshot::gas);
      if(strcmp("dm",str)==0) return(ySnapshot::dm);
      if(strcmp("dm1",str)==0) return(ySnapshot::dm1);
      if(strcmp("dm2",str)==0) return(ySnapshot::dm2);
      if(strcmp("dm3",str)==0) return(ySnapshot::dm3);
      if(strcmp("dm4",str)==0) return(ySnapshot::dm4);
      if(strcmp("dm5",str)==0) return(ySnapshot::dm5);
      if(strcmp("stars",str)==0) return(ySnapshot::stars);
      if(strcmp("all",str)==0) return(ySnapshot::all_types);

      cerr<<" ERROR: part_type '"<<str<<"' not recognized\n\a";
      throw(-1);
    }

  static Data get_data(const char *str)
    {
      if(strcmp("pos",str)==0) return(ySnapshot::pos);
      if(strcmp("x",str)==0) return(ySnapshot::x);
      if(strcmp("y",str)==0) return(ySnapshot::y);
      if(strcmp("z",str)==0) return(ySnapshot::z);

      if(strcmp("vel",str)==0) return(ySnapshot::vel);
      if(strcmp("x",str)==0) return(ySnapshot::x);
      if(strcmp("y",str)==0) return(ySnapshot::y);
      if(strcmp("z",str)==0) return(ySnapshot::z);
      
      cerr<<" ERROR: data '"<<str<<"' not recognized\n\a";
      throw(-1);
    }

  static DIM get_dim(ySnapshot::Data d)
    {
      switch(d)
	{
	case ySnapshot::pos:
	case ySnapshot::vel:
	  return(3);
	case ySnapshot::x:
	case ySnapshot::y:
	case ySnapshot::z:
	case ySnapshot::vx:
	case ySnapshot::vy:
	case ySnapshot::vz:
	  return(1);
	default:
	  break;
	}
      cerr<<" ERROR: data '"<<d<<"' not recognized\n\a";
      throw(-1);
    }

  virtual bool open(const char * name)
    {
      close();
      file.open(name,ios::binary);
      if(!file.is_open())
	{
	  cerr<<"Could not open '"<<name<<"'\n";
	  return(false);
	}
      filename = name;
      PANIC( cerr<<"OPEN: '"<<name<<"'\n"; );
      return(true);
    }

  virtual void printInfo()=0;

  virtual NUM Npart(PartType t)=0;
  virtual double Mpart(PartType t)=0;
  virtual double Lbox()=0;

  virtual void read(Data d, PartType t, NUM Nread=0)=0;
  virtual NUM next(char *dest, NUM n)=0;

 protected:

  ifstream file;
  bool swap_endian;

  virtual void init(){ filename=""; }
  virtual void close(){ if(file.is_open()) file.close(); init(); }

};

//----------------------------------------------------------------------
class yGadget: public ySnapshot
//----------------------------------------------------------------------
{
 public:

  yGadget(){ init(); }
  yGadget(const char *name){ open(name); }
  // ~yGadget(){ delete[] sel; }

  struct
  {
    int      Npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      Nall[6];
    int      flag_cooling;
    int      Nfiles;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  // fills to 256 Bytes
  } hdr;

  bool open(const char *name);

  void printInfo();

  NUM Npart(PartType t){ return( hdr.Npart[type(t)] ); }
  double Mpart(PartType t){ return( hdr.mass[type(t)]*1e10 ); }
  double Lbox(){ return(hdr.BoxSize); }

  void read(Data d, PartType t, NUM Nskip=0);
  NUM next(char *dest, NUM n);

 private:

  NUM type(PartType t)
    {
      switch(t)
	{
	case gas: return(0);
	case dm1: return(1);
	case dm2: return(2);
	case dm3: return(3);
	case dm4: return(4);
	case dm5: return(5);
	case stars: return(4);
	default:
	  cerr << " yGadget: Type "<<t<<" not valid!\n";
	  throw(-1);
	}
    }

  Data data_to_read;
  NUM N_to_read;

  void init(){ ySnapshot::init(); data_to_read=no_data; N_to_read=0; }

};

#endif
//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
