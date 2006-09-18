#ifndef GET_TIME_H
#define GET_TIME_H

double get_time(void);

#ifdef WIN32 
typedef struct _FILETIME {
  unsigned long dwLowDateTime;
  unsigned long dwHighDateTime;
} FILETIME;

void __stdcall GetSystemTimeAsFileTime(FILETIME*);
#endif

#endif //GET_TIME_H
