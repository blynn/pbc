#include <sys/time.h>
#include <time.h>

double pbc_get_time(void) {
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
