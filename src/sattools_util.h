/* > sattools_util.h
 *
 * sattools-specific pieces on top of the SGP4 library in extern/dundee_sgp4.
 */

#ifndef _SATTOOLS_UTIL_H
#define _SATTOOLS_UTIL_H

#include "sgdp4h.h"

/* Catalogue number from the command line: plain (101234) or alpha-5
 * (A1234). Returns -1 for an alpha-5-looking string that is not valid. */
long parse_satno(const char *s);

/* format_tle(), saying on stderr why the elements cannot be written as a
 * TLE when it fails. Returns 0 or -1. */
int format_tle_report(const orbit_t *orb, char line1[70], char line2[70]);

#endif /* _SATTOOLS_UTIL_H */
