#ifndef FIESTAS_H
#error "#------------------------------------------------------#"
#error "|                                                      |"
#error "| ERROR: yASCII.h should only be included by FiEstAS.h |"
#error "|                                                      |"
#error "#------------------------------------------------------#"
#else

//---------------------------------------------------------------------
class yASCII_data
//---------------------------------------------------------------------
{
 public:
  yASCII_data(){ init(); }
  yASCII_data(yASCII_data &f);
  yASCII_data(const char* fname){ init(); read(fname); }

  ~yASCII_data(){ reset(); }

  const char *file_name;
  int read(const char* name);

  int get_Nrows(){ return(Nrows); }
  int get_Ncolumns(){ return(Ncolumns); }
  int get_Ncomments(){ return(Ncomments); }

  double get_min(int col){ return(min[col]); }
  double get_max(int col){ return(max[col]); }
  double get_ave(int col){ return(ave[col]); }
  double get_sig(int col){ return(sig[col]); }

  double get(int row, int col){ return(column[col][row]); }

 private:
  vector<double> *column;
  double *max, *min, *ave, *sig;
  int Nrows, Ncolumns, Ncomments;

  void init()
  {
    file_name="";
    column = NULL;
    max = min = ave = sig = NULL;
    Nrows = Ncolumns = Ncomments = 0;
  }

  void reset()
  {
    if(column!=NULL)
      {
	// for(int c=0; c<Ncolumns; c++) column[c].clear(); // NEEDED?
	delete []column;
      }
    if(max!=NULL) delete []max;
    if(min!=NULL) delete []min;
    if(ave!=NULL) delete []ave;
    if(sig!=NULL) delete []sig;

    init();
  }
};

#endif
//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
