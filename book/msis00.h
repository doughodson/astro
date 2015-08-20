#ifndef _msis00_
#define _msis00_
/*     ----------------------------------------------------------------
*
*                              msis00.h
*
*  this file contains common routines for the nrlmsise-00 atmospheric model.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2007
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com
*
*     *****************************************************************
*
*    current :
*              31 mar 08  david vallado
*                           misc updates
*    changes :
*              15 mar 07  david vallado
*                           3rd edition baseline
*              6 aug 04  david vallado
*                          convert to c++
*              6 sep 03  david vallado
*                          fix low alt test cases (long in lpoly)
*             14 feb 03  david vallado
*                          misc updates
*             28 feb 02  david vallado
*                          finish double conversion
*              1 feb 02  nrl
*                          original baseline
*     http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
*     *****************************************************************       */

      #include "msiscom.h"

void gtd7
     (
       msistype& msis00r,
       lpolytype& lpoly,
       fittype& fit,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong, double stl,
       double f107a, double f107, double ap[8], int mass,
       double d[10], double t[3]
     );

void gtd7d
     (
       msistype& msis00r,
       lpolytype& lpoly,
       fittype& fit,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong,
       double stl, double f107a, double f107, double ap[8], int mass,
       double d[10], double t[3]
     );

void msis00init
     (
       msistype& msis00r
     );

#endif
