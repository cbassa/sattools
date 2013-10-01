#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

int main(int argc,char *argv[])
{
  struct timeval tv;

  for (;;) {
    gettimeofday(&tv,0);
    if (tv.tv_usec>999000)
      break;
  }

  return 0;
}
