/* > sattools_util.c
 *
 * sattools-specific pieces on top of the SGP4 library in extern/dundee_sgp4.
 */

#include "sgdp4h.h"
#include "elset.h"
#include "sattools_util.h"

#include <ctype.h>
#include <stdarg.h>
#include <string.h>

/* The SGP4 library reports errors through fatal_error(), and its own
 * version exits. The sattools tools report the error and carry on (they
 * skip an element set that fails to initialise and continue with the rest
 * of the catalogue), so they keep their own definition here; with this one
 * linked, the library's ferror.o is not pulled in. */
void fatal_error(const char *format, ...)
{
va_list arg_ptr;

    fflush(stdout);

    fprintf(stderr, "\nDundee Satellite Lab fatal run-time error:\n");

    va_start(arg_ptr, format);
    vfprintf(stderr, format, arg_ptr);
    va_end(arg_ptr);

    return;
}

int format_tle_report(const orbit_t *orb, char line1[70], char line2[70])
{
  if (format_tle(orb, line1, line2) == 0)
    return 0;

  if (orb->satno < 0 || orb->satno > ALPHA5_MAX)
    fprintf(stderr, "Object %d cannot be written as a TLE: catalog numbers stop at %ld (alpha-5 Z9999)\n",
            orb->satno, ALPHA5_MAX);
  else if (!isfinite(orb->rev) || !isfinite(orb->ecc) || !isfinite(orb->eqinc) || !isfinite(orb->ascn)
           || !isfinite(orb->argp) || !isfinite(orb->mnan) || !isfinite(orb->bstar))
    fprintf(stderr, "Elements of object %d are not valid numbers, no TLE written\n", orb->satno);
  else
    fprintf(stderr, "Elements of object %d do not fit in a TLE\n", orb->satno);

  return -1;
}

long parse_satno(const char *s)
{
  while (isspace((unsigned char)*s))
    s++;

  // Five characters with a leading letter: alpha-5
  if (strlen(s)==5 && isalpha((unsigned char)s[0]))
    return alpha5_to_number(s);

  return atol(s);
}
