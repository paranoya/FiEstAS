#include"FiEstAS.h"

#include <iomanip>
using namespace std;

//----------------------------------------------------------------------
bool yGadget::open(const char * name)
//----------------------------------------------------------------------
{
  if( ! ySnapshot::open(name) ) return(false);

  unsigned int i;
  file.read( (char*)&i, 4);
  switch(i)
    {
    case 256: swap_endian=false; break;
    case 0x00010000: swap_endian=true; break;
    default:
      cerr<<name<<" is not a Gadget snapshot! \n";
      close();
      return(false);
    }

  file.read( (char*)&hdr, 256);
  if(swap_endian)
    {
      cerr<<"Swap TO DO \n";
      throw(-1);
    }

  return(true);
}

//----------------------------------------------------------------------
void yGadget::printInfo()
//----------------------------------------------------------------------
{
  unsigned int i;

  cout << "\n---------------------------------------\n";
  cout << " yGadget snapshot: '" << filename << "'\n";
  cout << "---------------------------------------\n";

  cout <<" Npart ="; for(i=0;i<5;i++) cout<<setw(10)<<hdr.Npart[i]; cout<<endl;
  cout <<" mass  ="; for(i=0;i<5;i++) cout<<setw(10)<<hdr.mass[i];  cout<<endl;

  cout << " time = " << hdr.time;
  cout << "; z = " << hdr.redshift << endl;

  cout << " cool = " << hdr.flag_cooling;
  cout << "; sfr = " << hdr.flag_sfr;
  cout << "; fb = " << hdr.flag_feedback << endl;

  cout << " Nfiles = " << hdr.Nfiles << endl;
  if(hdr.Nfiles>1)
    {
      cout<<" Nall ="; for(i=0;i<5;i++) cout<<setw(10)<<hdr.Nall[i]; cout<<endl;
    }

  cout << " Lbox = " << hdr.BoxSize << endl;

  cout << " Omega0 = " << hdr.Omega0;
  cout << "; OmegaL = " << hdr.OmegaLambda;
  cout << "; h = " << hdr.HubbleParam << endl;

  cout << "---------------------------------------\n";
}

//----------------------------------------------------------------------
void yGadget::read(Data d, PartType t, NUM Nskip)
//----------------------------------------------------------------------
{
  PANIC( clog << " yGadget: Read data " << d << "type "<< t << endl; );

  streamoff offset=268; // 4+256+4+4
  data_to_read = d;
  N_to_read = Npart(t) - Nskip;

  for(NUM i=0; i<type(t); i++) Nskip += hdr.Npart[i];

  switch(d)
    {
    case pos:
    case x:
    case y:
    case z:
      offset += 12*Nskip;
      break;

    case vel:
    case vx:
    case vy:
    case vz:
      for(unsigned int i=0;i<5;i++) offset+=12*hdr.Npart[i];
      offset += 8;
      offset += 12*Nskip;
      break;
    default:
      cerr << " yGadget: Data "<<d<<" not recognized!\n";
    }

  file.seekg(offset,ios::beg);

  //  clog<<" read( data="<<d<<", type="<<t<<" ) -> offset="<<offset<<endl;
}

//----------------------------------------------------------------------
NUM yGadget::next(char *dest, NUM n)
//----------------------------------------------------------------------
{
  if( n > N_to_read ) n = N_to_read;
  if( n==0 || data_to_read==no_data ) return(0);

  NUM Nread = 0;
  switch(data_to_read)
    {
    case pos:
    case vel:
      Nread = file.readsome(dest,n*12) /12;
      break;
    default:
      cerr << " yGadget: Data "<<data_to_read<<" not implemented!\n";
    }

  N_to_read -= Nread;

  if(swap_endian)
    {
      cerr<<"Swap TO DO \n";
      throw(-1);
    }

  return(Nread);
}
