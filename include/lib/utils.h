#ifndef FIESTAS_H
#error "#-----------------------------------------------------#"
#error "|                                                     |"
#error "| ERROR: utils.h should only be included by FiEstAS.h |"
#error "|                                                     |"
#error "#-----------------------------------------------------#"
#else

#include<sstream>
//----------------------------------------------------------------------
template <class T>
T option(const char *opt, T deflt, int argmax,char** argv)
//----------------------------------------------------------------------
{
  T result = deflt;
  int l=strlen(opt);
  int arg=1;
  while( arg<argmax && strncmp(opt,argv[arg],l)!=0 ) arg++;
  if( arg<argmax )
  {
    string rest(argv[arg]+l);
    istringstream s(rest);
    s >> result;
  }
  return(result);
}

#endif
//------------------------------------------------------
//                              ... Paranoy@ Rulz! ;^D
//------------------------------------------------------
