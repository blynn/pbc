#include <sys/time.h>

#ifdef WIN32
#include "get_time.h"

int __cdecl gettimeofday(struct timeval* p, void* tz /* IGNORED */){
  union {
    long long ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } _now;
  
  GetSystemTimeAsFileTime( &(_now.ft) );
  p->tv_usec=(long)((_now.ns100 / 10LL) % 1000000LL );
  p->tv_sec= (long)((_now.ns100-(116444736000000000LL))/10000000LL);
  return;
}
#else
#include <time.h>
#endif

double get_time(void)
{
    static struct timeval last_tv, tv;
    static int first = 1;
    static double res = 0;

    if (first) {
        gettimeofday(&last_tv, NULL);
        first = 0;
        return 0;
    } else {
        gettimeofday(&tv, NULL);
        res += tv.tv_sec - last_tv.tv_sec;
        res += (tv.tv_usec - last_tv.tv_usec) / 1000000.0;
        last_tv = tv;

        return res;
    }
}


