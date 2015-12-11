// a smart parameter is just like a regular parameter, but it can be
// re-defined at the shell-environment.  Instead of
//
//      #define NUMBER 42
//      ...
//      printf("%g", NUMBER);
//
// do
//      SMART_PARAMETER(NUMBER,42)
//      ...
//      printf("%g", NUMBER());
//
// Notice that the environment only gets queried once, at the first use.
//

#define SMART_PARAMETER_INT(n, v) static int n(void)\
{\
  static int smapa_known_ ## n = 0;\
  static int smapa_value_ ## n = v;\
  if (!smapa_known_ ## n) {\
    int r;\
    char *sv = getenv(#n);\
    int y;\
    if (sv)\
    r = sscanf(sv, "%d", &y);\
    if (sv && r == 1)\
    smapa_value_ ## n = y;\
    smapa_known_ ## n = 1;\
  }\
  return smapa_value_ ## n;\
}

#define SMART_PARAMETER_FLOAT(n,v) static float n(void)\
{\
  static bool smapa_known_ ## n = false;\
  static float smapa_value_ ## n = v;\
  if (!smapa_known_ ## n) {\
    int r;\
    char *sv = getenv(#n);\
    float y;\
    if (sv)\
    r = sscanf(sv, "%f", &y);\
    if (sv && r == 1)\
    smapa_value_ ## n = y;\
    smapa_known_ ## n = true;\
  }\
  return smapa_value_ ## n;\
}

#define SMART_PARAMETER_DOUBLE(n,v) static double n(void)\
{\
  static bool smapa_known_ ## n = false;\
  static double smapa_value_ ## n = v;\
  if (!smapa_known_ ## n) {\
    int r;\
    char *sv = getenv(#n);\
    double y;\
    if (sv)\
    r = sscanf(sv, "%lf", &y);\
    if (sv && r == 1)\
    smapa_value_ ## n = y;\
    smapa_known_ ## n = true;\
  }\
  return smapa_value_ ## n;\
}
