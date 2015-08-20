#ifndef _asttime_h_
#define _asttime_h_
/* --------------------------------------------------------------------
*
*                                asttime.h
*
*    this file contains miscallaneous time functions.
*
*    current :
*               4 may 09  david vallado
*                           misc updates
*    changes :
*              13 feb 08  david vallado
*                           add getmon
*              15 mar 07  david vallado
*                           3rd edition baseline
*              20 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
  ---------------------------------------------------------------------- */

#include "astmath.h"

// object definitions
typedef enum
        {
          eTo, eFrom
        } edirection;



int getmon
    (
      char instr[3]
    );

int getday
    (
      char instr[3]
    );

int dayofweek
    (
      double jd
    );
    
void    jday
        (
          int year, int mon, int day, int hr, int minute, double sec,
          double& jd
        );

void    jdayall
        (
          int year, int mon, int day, int hr, int minute, double sec,
          char whichtype, double& jd
        );

void    days2mdhms
        (
          int year, double days,
          int& mon, int& day, int& hr, int& minute, double& sec
        );

void    invjday
        (
          double jd,
          int& year, int& mon, int& day,
          int& hr, int& minute, double& sec
        );

void    finddays
        (
          int year, int month, int day, int hr, int minute,
          double sec, double& days
        );

double  gstime
        (
          double jdut1
        );

void    lstime
        (
         double lon, double jdut1, double& lst, double& gst
        );

void    hms_sec
        (
          int& hr, int& min, double& sec, edirection direct, double& utsec
        );

void    hms_ut
        (
          int& hr, int& min, double& sec, edirection direct,  double& ut
        );

void    hms_rad
        (
          int& hr, int& min, double& sec, edirection direct, double& hms
        );

void    dms_rad
        (
          int& deg, int& min, double& sec, edirection direct, double& dms
        );

double jd2sse
        (
          double jd
        );

void    convtime
        (
          int year, int mon, int day, int hr, int min, double sec, int timezone,
          double dut1,  int dat,
          double& ut1,  double& tut1,  double& jdut1, double& utc, double& tai,
          double& tt,   double& ttt,   double& jdtt,  double& tcg, double& tdb,
          double& ttdb, double& jdtdb, double& tcb
        );

#endif

