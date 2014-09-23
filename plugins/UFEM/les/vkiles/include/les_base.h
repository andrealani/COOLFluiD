#include <les_interface_v4.h>

#ifndef LES_BASE_H
#define LES_BASE_H

/* Helping things */
#define MSG(formatstr,...) \
{\
  printf("%s:%d:%s ",(char*)__FILE__,__LINE__,(char*)__FUNCTION__);\
  printf(formatstr,__VA_ARGS__);\
  fflush(stdout);\
}

#define ERRMSG(ret_code,formatstr,...) \
{\
  MSG(formatstr,__VA_ARGS__);\
  exit(ret_code);\
}

#define PUTPARAM(paramstr,formatstr,defaultstr,...) \
{\
  char paramval[100];\
  les_put_config_param(paramstr,paramval);\
  if (strcmp(paramval,"NOT_FOUND")==0){\
    MSG("'%s': value '%s', switching to default '%s' with '%s'.\n",paramstr,paramval,defaultstr,formatstr);\
    sprintf(paramval,"%s",defaultstr);\
  }\
  sscanf(paramval,formatstr,__VA_ARGS__);\
}

#define GETPARAM(paramstr,formatstr,...) \
{\
  if (strcmp(param_name,paramstr)==0) sprintf(param_value,formatstr,__VA_ARGS__);\
}


#endif /* LES_BASE_H */
