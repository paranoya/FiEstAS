#ifndef FIESTAS_H
#error "#-----------------------------------------------------#"
#error "|                                                     |"
#error "| ERROR: yTree.h should only be included by FiEstAS.h |"
#error "|                                                     |"
#error "#-----------------------------------------------------#"
#else

#include<cmath>
#include<vector>
using namespace std;

//----------------------------------------------------------------------
class yKernel
//----------------------------------------------------------------------
{
  public:
  yKernel(){};
  virtual ~yKernel(){};
  
  virtual double Kernel(double u){ return(0.); };
  virtual double Int_K(double u){ return(0.); };
  
  double Kernel(double u0, double u1)
  {
    if(u0>=1. || u1<=-1.) return(0.);
    if(u0<-1.) u0=-1.;
    if(u1>1.) u1=1.;
    return( Int_K(u1) - Int_K(u0) );
  }
};

//----------------------------------------------------------------------
class TopHat : public yKernel
//----------------------------------------------------------------------
{
  public:

  double Kernel(double u){ if(fabs(u)<1.) return(.5); else return(0.); }
  double Int_K(double u){ return(.5*u); }
};

//----------------------------------------------------------------------
class TSC : public yKernel
//----------------------------------------------------------------------
{
  public:
  double Kernel(double u){ return( std::max(0., 1.-fabs(u)) ); }
  double Int_K(double u){ if(u<0.) return( (u+1)*(u+1)/2 ); else return( .5 + u*(1-u/2) ); }
};

//----------------------------------------------------------------------
class Epanechnikov : public yKernel
//----------------------------------------------------------------------
{
  public:
  double Kernel(double u){ return( std::max(0., .75*(1.-u*u)) ); }
  double Int_K(double u){ return( u*(3-u*u)/4 ); }
};

//----------------------------------------------------------------------
class yTree
//----------------------------------------------------------------------
{
 public:
  yTree(DIM d=0, bool quiet=false);
  ~yTree();

  void read_ASCII(const char * name);
  void read_Gadget(const char *name, vector<yGadget::Data> data, yGadget::PartType part_type);

  void init_data(DIM d);
  void delete_data();
  void null_data();
  void init_extra(DIM d);
  void delete_extra();
  void null_extra();

  void add(vector<DATA> &xn){ for(DIM d=0;d<D;d++) add_xd(d,xn[d]); }
  void add(DATA *xn){ for(DIM d=0;d<D;d++) add_xd(d,xn[d]); }
  void add_xd(DIM d, DATA xd);
  void add_extra(DIM d, DATA xd);

  NUM get_N(){ return(Ntot); }
  NUM get_D(){ return(D); }
  DATA get_min(DIM d);
  DATA get_max(DIM d);
  void get_Bbox(DATA *b);

  void branch();

  NUM get_n(DATA* x0);
  DATA get_data(DIM d, NUM n);
  DATA get_extra(DIM d, NUM n);

  void copy_data_to(yTree *tree);
  void copy_extra_to(yTree *tree);
  void copy_extra_d_to(DIM d1, yTree *tree,DIM d2);

  void push_box(DATA *x0, vector<DATA> *box_out);
  void push_box(NUM n,    vector<DATA> *box_out);
  NUM push_cells(vector<DATA> *cells, vector<NUM> *points=NULL, bool saveCells=true);

  void get_split(NUM n, DIM &d, DATA &x);

  void write_data(const char *fname);
  void write_extra(const char *fname);

  double FiEstAS(DATA *x0, bool balloon=true, bool subtract_bias=false);
  void set_hsmooth(yKernel *K, double M_0);
  double smooth(DATA *x0, bool subtract_bias=false);
  double get_hsmooth(DATA *x0, DATA *h0);
  double set_metric(vector<NUM> &dim, vector<double> &scale);

 private:

  bool silent;

  // --------------------------- data
  DIM D;
  NUM Ntot;
  vector<DATA> *x;
  DATA *min, *max;
  NUM *Nsplits;     // number of times each dimension has been split
  bool tilt_left;  // what happens when x==xs

  // --------------------------- extra
  DIM D_extra;
  vector<DATA> *extra;

  // --------------------------- scratch1
  NUM *Nb;              // histogram
  DATA *xmin, *xmax; // boundaries
  DIM ds, *d_s;        // dimension to split
  DATA xs, *x_s;       // split point
  NUM Nl,Nr, *N_l;    // particles on each side (tentative or real)

  void init_scratch1(NUM B);
  void delete_scratch1();
  void null_scratch1();

  void branch(NUM n0, NUM n1, DATA *x0, DATA *x1);
  void choose_axis(NUM n0, NUM n1, DATA *x_0, DATA *x_1);
  void spawn_children(NUM n0, NUM n1);

  // --------------------------- scratch2
  NUM node;   // current node
  NUM id;     // current particle id
  DATA *box;  // current boundaries

  void init_scratch2();
  void delete_scratch2();
  void null_scratch2();

  // --------------------------- smooth
  vector<DATA> *hmin, *hmax, *h;
  double M0, bias, M_found, *h_found;
  yKernel *Kernel;
  
  void init_smooth();
  void delete_smooth();
  void null_smooth();

  void set_h(DATA *h_min, DATA *h_max);//, vector<double> &h0);
  double compute_smooth(DATA *x0);
  double compute_hsmooth(DATA *x0, DATA *h0);
  double compute_FiEstAS(DATA *box0, double Mmax);
  double compute_FiEstAS_plain(DATA *box0, double Mmax);

  // ----------------------------------------------

  void compute_n(DATA *x0);
  void compute_box(DATA *x0);
  void compute_box(NUM n);
  NUM compute_cells(vector<DATA> *cells, vector<NUM> *points, bool saveCells);

};

#endif
//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
