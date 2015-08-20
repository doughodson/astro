/*       ----------------------------------------------------------------
*
*                              msis00.cpp
*
*  this file contains common routines for the nrlmsise-00 atmospheric model.
*  note that the indices are increased by one so they will appear the same as
*  original fortran source.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2007
*                            by David Vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              31 mar 08  david vallado
*                           misc updates
*    changes :
*              15 mar 07  david vallado
*                           3rd edition baseline
*              12 nov 04  david vallado
*                           fix for ap array
*               4 oct 04  david vallado
*                           misc updates
*               6 aug 04  david vallado
*                           convert to c++
*               6 sep 03  david vallado
*                           fix low alt test cases (llong in lpoly)
*              14 feb 03  david vallado
*                           misc updates
*              28 feb 02  david vallado
*                           finish double conversion
*               1 feb 02  nrl
*                           original baseline
*     http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
*     *****************************************************************       */

      #include "msis00.h"
      #include "stdafx.h"

/*    ******************* local routines ******************************       */
void gts7
     (
       msistype& msis00r,
       lpolytype& lpoly,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
       double f107, double ap[8], int mass, double d[10], double t[3]
     );

double scalh
       (
         parmbtype& parmb,
         double alt, double xm, double temp
       );

double g0
       (
         double a, double p[151]
       );

double sumex
       (
         double ex
       );
double sg0
       (
         double ex, double p[151], double ap[8]
       );

double globe7
       (
         cswtype& csw,
         lpolytype& lpoly,
         double yrd, double sec, double lat, double llong, double tloc, double f107a,
         double f107, double ap[8], double p[151]
       );

double glob7s
       (
         lpolytype& lpoly,
         cswtype& csw,
         double p[101]
       );



/*  ----------------------------------------------------------------------
*
*     nrlmsise-00
*     -----------
*          neutral atmosphere empirical model from the surface to lower
*          exosphere
*
*          new features:
*            *extensive satellite drag database used in model generation
*            *revised o2 (and o) in lower thermosphere
*            *additional nonlinear solar activity term
*            *"anomalous oxygen" number density, output d[9]
*             at high altitudes (> 500 km), hot atomic oxygen or ionized
*             oxygen can become appreciable for some ranges of subroutine
*             inputs, thereby affecting drag on satellites and debris. we
*             group these species under the term "anomalous oxygen," since
*             their individual variations are not presently separable with
*             the drag data used to define this model component.
*
*          subroutines for special outputs:
*
*          high altitude drag: effective total mass density
*          (subroutine gtd7d, output d[6])
*             for atmospheric drag calculations at altitudes above 500 km,
*             call subroutine gtd7d to compute the "effective total mass
*             density" by including contributions from "anomalous oxygen."
*             see "notes on output variables" below on d[6].
*
*          pressure grid (subroutine ghp7)
*            see subroutine ghp7 to specify outputs at a pressure level
*            rather than at an altitude.
*
*          output in m-3 and kg/m3:   call meters(.true.)
*
*     input variables:
*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366))
*                (year ignored in current model)
*          sec - ut(sec)
*          alt - altitude(km)
*          glat - geodetic latitude(deg)
*          glong - geodetic longitude(deg)
*          stl - local apparent solar time(hrs; see note below)
*          f107a - 81 day average of f10.7 flux (centered on day ddd)
*          f107 - daily f10.7 flux for previous day
*          ap - magnetic index(daily) or when msis00r.csw.sw[9]=-1. :
*             - array containing:
*               [1] daily ap
*               [2] 3 hr ap index for current time
*               [3] 3 hr ap index for 3 hrs before current time
*               [4] 3 hr ap index for 6 hrs before current time
*               [5] 3 hr ap index for 9 hrs before current time
*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
*                      to current time
*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
*                      to current time
*          mass - mass number (only density for selected gas is
*                   calculated.  mass 0 is temperature.  mass 48 for all.
*                   mass 17 is anomalous o only.)
*
*     notes on input variables:
*          ut, local time, and longitude are used independently in the
*          model and are not of equal importance for every situation.
*          for the most physically realistic calculation these three
*          variables should be consistent (stl=sec/3600+glong/15).
*          the equation of time departures from the above formula
*          for apparent local time can be included if available but
*          are of minor importance.
*
*          f107 and f107a values used to generate the model correspond
*          to the 10.7 cm radio flux at the actual distance of the earth
*          from the sun rather than the radio flux at 1 au. the following
*          site provides both classes of values:
*          ftp://ftp.ngdc.noaa.gov/stp/solar_data/solar_radio/flux/
*
*          f107, f107a, and ap effects are neither large nor well
*          established below 80 km and these parameters should be set to
*          150., 150., and 4. respectively.
*
*     output variables:
*          d[1] - he number density(cm-3)
*          d[2] - o number density(cm-3)
*          d[3] - n2 number density(cm-3)
*          d[4] - o2 number density(cm-3)
*          d[5] - ar number density(cm-3]
*          d[6] - total mass density(gm/cm3]
*          d[7] - h number density(cm-3]
*          d[8] - n number density(cm-3]
*          d[9] - anomalous oxygen number density(cm-3]
*          t[1] - exospheric temperature
*          t[2] - temperature at alt
*
*     notes on output variables:
*          to get output in m-3 and kg/m3:   call meters(.true.)
*
*          o, h, and n are set to zero below 72.5 km
*
*          t[1], exospheric temperature, is set to global average for
*          altitudes below 120 km. the 120 km gradient is left at global
*          average value for altitudes below 72 km.
*
*          d[6], total mass density, is not the same for subroutines gtd7
*          and gtd7d
*
*            subroutine gtd7 -- d[6] is the sum of the mass densities of the
*            species labeled by indices 1-5 and 7-8 in output variable d.
*            this includes he, o, n2, o2, ar, h, and n but does not include
*            anomalous oxygen (species index 9).
*
*            subroutine gtd7d -- d[6] is the "effective total mass density
*            for drag" and is the sum of the mass densities of all species
*            in this model, including anomalous oxygen.
*
*     msis00r.csw.switches: the following is for test and special purposes:
*
*          to turn on and off particular variations call tselec(msis00r.csw.sw),
*          where msis00r.csw.sw is a 25 element array containing 0. for off, 1.
*          for on, or 2. for main effects off but cross terms on
*          for the following variations
*                 1 - f10.7 effect on mean  2 - time independent
*                 3 - symmetrical annual    4 - symmetrical semiannual
*                 5 - asymmetrical annual   6 - asymmetrical semiannual
*                 7 - diurnal               8 - semidiurnal
*                 9 - daily ap             10 - all ut/long effects
*                11 - longitudinal         12 - ut and mixed ut/long
*                13 - mixed ap/ut/long     14 - terdiurnal
*                15 - departures from diffusive equilibrium
*                16 - all tinf var         17 - all tlb var
*                18 - all tn1 var           19 - all s var
*                20 - all tn2 var           21 - all nlb var
*                22 - all tn3 var           23 - turbo scale height var
*
*          to get current values of msis00r.csw.sw: call tretrv(msis00r.csw.sw)
*
* --------------------------------------------------------------------------- */

void gtd7
     (
       msistype& msis00r,
       lpolytype& lpoly,
       fittype& fit,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong, double stl,
       double f107a, double f107, double ap[8], int mass,
       double d[10], double t[3]
     )
        {
        double ds[10],ts[3], v1, xlat, xmm, altt, dm28m, dmc, dz28, dmr, tz;
        int i, mss;

// make sure and initilize the data through msis00init at the start

	int mn3 = 5;
	double zn3[6] = {0, 32.5,20.0,15.0,10.0,0.0};
	int mn2 = 4;
	double zn2[5] = {0, 72.5,55.0,45.0,32.5};
	double zmix = 62.5;

        double mssl = -999;
        int sv[26] = {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

// treat as local, but save in between calls
        static double alast = 99999.0;

//------------------------------ begin ---------------------------------
        if (msis00r.csw.isw != 64999)
            tselec(msis00r.csw, sv);

        // ----  test for changed input
        v1 = vtst(msis00r.csw, iyd,sec,glat,glong,stl,f107a,f107,ap,1);

        // ---- latitude variation of gravity (none for msis00r.csw.sw[2]=0)
        xlat = glat;
        if (msis00r.csw.sw[2] == 0)
            xlat = 45.0;
        glatf(xlat, msis00r.parmb.gsurf, msis00r.parmb.re);

        xmm = msis00r.lower.pdm[5][3];

        // ---- thermosphere/mesosphere (above zn2[1])
        if (alt >= zn2[1])
            altt = alt;
          else
            altt = zn2[1];

        mss = mass;
        // ---- only calculate n2 in thermosphere if alt in mixed region
        if (alt < zmix && mass > 0)
            mss = 28;

        // ---- only calculate thermosphere if input parameters changed
        // ----   or altitude above zn2[1] in mesosphere
        if (v1 == 1.0 | alt > zn2[1] | alast > zn2[1] | mss != mssl)
         {
            gts7(msis00r, lpoly, lsqv, iyd, sec, altt, glat, glong, stl,
                 f107a, f107, ap, mss, ds, ts);
            dm28m = msis00r.dmix.dm28;
            // ----   metric adjustment
            if (msis00r.metsel.imr == 1)
                dm28m = msis00r.dmix.dm28*1.0e6;
            mssl = mss;
          }
        t[1] = ts[1];
        t[2] = ts[2];
        if (alt >= zn2[1])
          {
            for (i = 1; i <= 9; i++)
                d[i] = ds[i];
            goto ten;
          }

        // ---- lower mesosphere/upper stratosphere [between zn3[1] and zn2[1]]
        // ----   temperature at nodes and gradients at end nodes
        // ----   inverse temperature a linear function of spherical harmonics
        // ----   only calculate nodes if input changed
        if (v1 == 1.0 | alast >= zn2[1])
          {
            msis00r.meso.tgn2[1] = msis00r.meso.tgn1[2];
            msis00r.meso.tn2[1] = msis00r.meso.tn1[5];
            msis00r.meso.tn2[2] = msis00r.parm.pma[1][1]*msis00r.mavg.pavgm[1] /
                  (1.0-msis00r.csw.sw[20] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[1]));
            msis00r.meso.tn2[3] = msis00r.parm.pma[1][2]*msis00r.mavg.pavgm[2] /
                  (1.0-msis00r.csw.sw[20] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[2]));
            msis00r.meso.tn2[4] = msis00r.parm.pma[1][3]*msis00r.mavg.pavgm[3] /
                  (1.0-msis00r.csw.sw[20] * msis00r.csw.sw[22]*
                   glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[3]));
            msis00r.meso.tgn2[2] = msis00r.mavg.pavgm[9]*msis00r.parm.pma[1][10]*
                  (1.0+msis00r.csw.sw[20] * msis00r.csw.sw[22]*
                   glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[10]))
                    *msis00r.meso.tn2[4] * msis00r.meso.tn2[4] /
                    pow((msis00r.parm.pma[1][3] * msis00r.mavg.pavgm[3]),2);
            msis00r.meso.tn3[1] = msis00r.meso.tn2[4];
          }
        if (alt < zn3[1])
          {
            // ---- lower stratosphere and troposphere [below zn3[1]]
            // ----   temperature at nodes and gradients at end nodes
            // ----   inverse temperature a linear function of spherical harmonics
            // ----   only calculate nodes if input changed
            if (v1 == 1.0 | alast >= zn3[1])
              {
                msis00r.meso.tgn3[1] = msis00r.meso.tgn2[2];
                msis00r.meso.tn3[2]  = msis00r.parm.pma[1][4]*msis00r.mavg.pavgm[4]
                    /(1.0-msis00r.csw.sw[22]*glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[4]));
                msis00r.meso.tn3[3]  = msis00r.parm.pma[1][5]*msis00r.mavg.pavgm[5]
                    /(1.0-msis00r.csw.sw[22]*glob7s(lpoly, msis00r.csw,msis00r.parm.pmadav[5]));
                msis00r.meso.tn3[4]  = msis00r.parm.pma[1][6]*msis00r.mavg.pavgm[6]
                    /(1.0-msis00r.csw.sw[22]*glob7s(lpoly, msis00r.csw,msis00r.parm.pmadav[6]));
                msis00r.meso.tn3[5]  = msis00r.parm.pma[1][7]*msis00r.mavg.pavgm[7]
                    /(1.0-msis00r.csw.sw[22]*glob7s(lpoly, msis00r.csw,msis00r.parm.pmadav[7]));
                msis00r.meso.tgn3[2] = msis00r.parm.pma[1][8]*msis00r.mavg.pavgm[8]*
                    (1.0+msis00r.csw.sw[22]*glob7s(lpoly, msis00r.csw,msis00r.parm.pmadav[8]))
                      *msis00r.meso.tn3[5]*msis00r.meso.tn3[5] /
                      pow((msis00r.parm.pma[1][7]*msis00r.mavg.pavgm[7]),2);
              }
          }

        if (mass == 0)
            goto fifty;
        // ----    linear transition to full mixing below zn2[1]
        dmc = 0;
        if (alt > zmix)
            dmc = 1.0-(zn2[1]-alt)/(zn2[1]-zmix);
        dz28 = ds[3];

        // ----***** n2 density ****
        dmr  = ds[3]/dm28m-1.0;
        d[3] = densm(msis00r.parmb, fit, lsqv, alt, dm28m, xmm, tz, mn3, zn3,
                msis00r.meso.tn3, msis00r.meso.tgn3, mn2, zn2, msis00r.meso.tn2,
                msis00r.meso.tgn2);
        d[3] = d[3]*(1.0+dmr*dmc);

        // ----***** he density ****
        d[1] = 0.0;
        if (mass == 4 | mass == 48)
          {
            dmr  = ds[1]/(dz28*msis00r.lower.pdm[2][1])-1.0;
            d[1] = d[3]*msis00r.lower.pdm[2][1]*(1.0+dmr*dmc);
          }

        // ----**** o density ****
        d[2] = 0.0;
        d[9] = 0.0;

        // ----***** o2 density ****
        d[4] = 0.0;
        if (mass == 32 | mass == 48)
        {
            dmr  = ds[4]/(dz28*msis00r.lower.pdm[2][4])-1.0;
            d[4] = d[3]*msis00r.lower.pdm[2][4]*(1.0+dmr*dmc);
        }

        // ----***** ar density ****
        d[5] = 0.0;
        if (mass == 40 | mass == 48)
          {
            dmr  = ds[5]/(dz28*msis00r.lower.pdm[2][5])-1.0;
            d[5] = d[3]*msis00r.lower.pdm[2][5]*(1.0+dmr*dmc);
          }

        // ----***** hydrogen density ****
        d[7] = 0.0;

        // ----***** atomic nitrogen density ****
        d[8] = 0.0;

        // ---- total mass density

        if (mass == 48)
          {
            d[6] = 1.66e-24*(4.0*d[1]+16.0*d[2]+28.0*d[3]+ 32.0*d[4]+40.0*d[5]+d[7]+14.0*d[8]);
            if (msis00r.metsel.imr == 1)
                d[6] = d[6]/1000.0;
          }
        t[2] = tz;
   ten:
        goto ninety;
   fifty:
        msis00r.gts3c.dd = densm(msis00r.parmb, fit, lsqv,
                   alt, 1.0, 0.0, tz, mn3, zn3, msis00r.meso.tn3, msis00r.meso.tgn3, mn2, zn2,
                   msis00r.meso.tn2, msis00r.meso.tgn2);
        t[2] = tz;
   ninety:
        alast = alt;
        }

/* -----------------------------------------------------------------------
*
*     nrlmsise-00
*     -----------
*          this void provides effective total mass density for
*          output d[6] which includes contributions from "anomalous
*          oxygen" which can affect satellite drag above 500 km.  this
*          void is part of the distribution package for the
*          neutral atmosphere empirical model from the surface to lower
*          exosphere.  see void gtd7 for more extensive comments.
*
*     input variables:
*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366])
*                (year ignored in current model)
*          sec - ut(sec)
*          alt - altitude(km)
*          glat - geodetic latitude(deg)
*          glong - geodetic longitude(deg)
*          stl - local apparent solar time(hrs; see note below)
*          f107a - 81 day average of f10.7 flux (centered on day ddd)
*          f107 - daily f10.7 flux for previous day
*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
*             - array containing:
*               [1] daily ap
*               [2] 3 hr ap index for current time
*               [3] 3 hr ap index for 3 hrs before current time
*               [4] 3 hr ap index for 6 hrs before current time
*               [5] 3 hr ap index for 9 hrs before current time
*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
*                      to current time
*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
*                      to current time
*          mass - mass number (only density for selected gas is
*                   calculated.  mass 0 is temperature.  mass 48 for all.
*                   mass 17 is anomalous o only.)
*
*     notes on input variables:
*          ut, local time, and longitude are used independently in the
*          model and are not of equal importance for every situation.
*          for the most physically realistic calculation these three
*          variables should be consistent (stl = sec/3600+glong/15).
*          the equation of time departures from the above formula
*          for apparent local time can be included if available but
*          are of minor importance.
*
*          f107 and f107a values used to generate the model correspond
*          to the 10.7 cm radio flux at the actual distance of the earth
*          from the sun rather than the radio flux at 1 au.
*
*     output variables:
*          d[1] - he number density(cm-3]
*          d[2] - o number density(cm-3]
*          d[3] - n2 number density(cm-3]
*          d[4] - o2 number density(cm-3]
*          d[5] - ar number density(cm-3]
*          d[6] - total mass density(gm/cm3] [includes anomalous oxygen]
*          d[7] - h number density(cm-3]
*          d[8] - n number density(cm-3]
*          d[9] - anomalous oxygen number density(cm-3]
*          t[1] - exospheric temperature
*          t[2] - temperature at alt
* ----------------------------------------------------------------------      */

void gtd7d
     (
       msistype& msis00r,
       lpolytype& lpoly,
       fittype& fit,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong,
       double stl, double f107a, double f107, double ap[8], int mass,
       double d[10], double t[3]
     )
        {
//------------------------------ begin ---------------------------------
         gtd7( msis00r, lpoly, fit, lsqv,
               iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t);

         // ---- total mass density
         if (mass == 48)
           {
             d[6] = 1.66e-24 * (4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] +
                                40.0*d[5] + d[7] + 14.0*d[8] + 16.0*d[9]);
             if (msis00r.metsel.imr == 1)
                 d[6] = d[6]/1000.0;
           }
        }

/* ----------------------------------------------------------------------
*         find altitude of pressure surface (press) from gtd7
*     input:
*          iyd - year and day as yyddd
*          sec - ut(sec)
*          glat - geodetic latitude(deg)
*          glong - geodetic longitude(deg)
*          stl - local apparent solar time(hrs)
*          f107a - 3 month average of f10.7 flux
*          f107 - daily f10.7 flux for previous day
*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
*             - array containing:
*               [1] daily ap
*               [2] 3 hr ap index for current time
*               [3] 3 hr ap index for 3 hrs before current time
*               [4] 3 hr ap index for 6 hrs before current time
*               [5] 3 hr ap index for 9 hrs before current time
*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
*                      to current time
*               [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
*                      to current time
*          press - pressure level(mb)
*     output:
*          alt - altitude(km)
*          d[1] - he number density(cm-3]
*          d[2] - o number density(cm-3]
*          d[3] - n2 number density(cm-3]
*          d[4] - o2 number density(cm-3]
*          d[5] - ar number density(cm-3]
*          d[6] - total mass density(gm/cm3]
*          d[7] - h number density(cm-3]
*          d[8] - n number density(cm-3]
*          d[9] - hot o number density(cm-3]
*          t[1] - exospheric temperature
*          t[2] - temperature at alt
* ----------------------------------------------------------------------      */

      void ghp7
          (
            msistype& msis00r,
            lpolytype& lpoly,
            fittype& fit,
            lsqvtype& lsqv,
            int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
            double f107, double ap[8],
            double d[10], double t[3], double press
           )
        {
        double pl, zi, cl, cl2, cd, ca, z, l, xm, g, xn, p, sh, diff;
        int iday;

        double bm   = 1.3806e-19;
        double rgas = 831.4;
        double test = 0.00043;
        int   ltest = 12;

//------------------------------ begin ---------------------------------
        pl = log10(press);

        //       initial altitude estimate
        if (pl >= -5.0)
          {
            if (pl > 2.5)   zi = 18.06*(3.00-pl);
            if (pl > 0.75 &&  pl <= 2.5) zi = 14.98*(3.08-pl);
            if (pl > -1.0 &&  pl <= 0.75) zi = 17.8*(2.72-pl);
            if (pl > -2.0 &&  pl <= -1.0) zi = 14.28*(3.64-pl);
            if (pl > -4.0 &&  pl <= -2.0) zi = 12.72*(4.32-pl);
            if (pl <= -4.0) zi = 25.3*(0.11-pl);
            iday = iyd - int(iyd/1000.0)*1000.0;
      printf("not iday %12i \n", iday);
            iday = iyd % 1000;
      printf("not iday %12i \n", iday);
            cl  = glat/90.0;
            cl2 = cl*cl;
            if (iday < 182)
                cd = 1.0-iday/91.25;
            if (iday >= 182)
                cd = iday/91.25 - 3.0;
            ca = 0;
            if (pl > -1.11 && pl <= -.23)
                ca = 1.0;
            if (pl > -0.23)
                ca = (2.79-pl)/(2.79 + 0.23);
            if (pl <= -1.11 && pl > -3.0)
                ca = (-2.93-pl)/(-2.93 + 1.11);
            z = zi - 4.87*cl*cd*ca - 1.64*cl2*ca + 0.31*ca*cl;
          }
        if (pl < -5.0)
            z = 22.0*pow( pl+4.0,2 ) + 110.0;

//       iteration loop
        l = 0;
   ten:
        l = l+1;
        gtd7( msis00r, lpoly, fit, lsqv,
              iyd, sec, z, glat, glong, stl, f107a, f107, ap, 48, d, t);
        xn = d[1] + d[2] + d[3] + d[4] + d[5] + d[7] + d[8];
        p  = bm*xn*t[2];
        if (msis00r.metsel.imr == 1)
            p = p*1.0e-6;
        diff = pl - log10(p);
        if (fabs(diff) < test  |  l == ltest)
            goto twenty;
        xm = d[6]/xn/1.66e-24;
        if (msis00r.metsel.imr == 1)
            xm = xm*1.0e3;
        g  = msis00r.parmb.gsurf/pow( (1.0+z/msis00r.parmb.re),2 );
        sh = rgas*t[2]/(xm*g);

        // ----   new altitude estimate using scale height
        if (l < 6)
            z = z-sh*diff*2.302;
          else
            z = z-sh*diff;
        goto ten;
   twenty:
        if (l == ltest)
            printf("not converging for press %12.2f %12.2f \n", press,diff);
        alt = z;
        }

/*  ----------------------------------------------------------------------
*     thermospheric portion of nrlmsise-00
*     see gtd7 for more extensive comments
*
*          output in m-3 and kg/m3:   call meters(.true.)
*
*     input variables:
*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366])
*                (year ignored in current model)
*          sec - ut(sec)
*          alt - altitude(km) (>72.5 km)
*          glat - geodetic latitude(deg)
*          glong - geodetic longitude(deg)
*          stl - local apparent solar time(hrs; see note below)
*          f107a - 81 day average of f10.7 flux (centered on day ddd)
*          f107 - daily f10.7 flux for previous day
*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
*             - array containing:
*               [1] daily ap
*               [2] 3 hr ap index for current time
*               [3] 3 hr ap index for 3 hrs before current time
*               [4] 3 hr ap index for 6 hrs before current time
*               [5] 3 hr ap index for 9 hrs before current time
*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
*                      to current time
*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
*                      to current time
*          mass - mass number (only density for selected gas is
*                   calculated.  mass 0 is temperature.  mass 48 for all.
*                   mass 17 is anomalous o only.)
*
*     notes on input variables:
*          ut, local time, and longitude are used independently in the
*          model and are not of equal importance for every situation.
*          for the most physically realistic calculation these three
*          variables should be consistent (stl = sec/3600+glong/15).
*          the equation of time departures from the above formula
*          for apparent local time can be included if available but
*          are of minor importance.
*
*          f107 and f107a values used to generate the model correspond
*          to the 10.7 cm radio flux at the actual distance of the earth
*          from the sun rather than the radio flux at 1 au. the following
*          site provides both classes of values:
*          ftp:*ftp.ngdc.noaa.gov/stp/solar_data/solar_radio/flux/
*
*          f107, f107a, and ap effects are neither large nor well
*          established below 80 km and these parameters should be set to
*          150.0, 150.0, and 4.0 respectively.
*
*     output variables:
*          d[1] - he number density(cm-3]
*          d[2] - o number density(cm-3]
*          d[3] - n2 number density(cm-3]
*          d[4] - o2 number density(cm-3]
*          d[5] - ar number density(cm-3]
*          d[6] - total mass density(gm/cm3] [anomalous o not included]
*          d[7] - h number density(cm-3]
*          d[8] - n number density(cm-3]
*          d[9] - anomalous oxygen number density(cm-3]
*          t[1] - exospheric temperature
*          t[2] - temperature at alt
* -----------------------------------------------------------------------     */

void gts7
     (
       msistype& msis00r,
       lpolytype& lpoly,
       lsqvtype& lsqv,
       int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
       double f107, double ap[8], int mass, double d[10], double t[3]
     )
        {
        double v2, yrd, tinf, ps,
               b01    , b04    , b14    , b16    , b28    , b32    , b40   ,
               db16h  , dbh16h , dm14   ,
               g1     , g14    , g16    , g16h   , g28    , g32    , g4    , g40   ,
               hc01   , hc04   , hc14   , hc16   , hc216  , hc32   , hc40,
               hcc01  , hcc14  , hcc16  , hcc232 , hcc32  ,
               rc01   , rc14   , rc16   , rc32   , t2     , tho    , tz    , xmd   , xmm   ,
               z      , zc01   , zc04   , zc14   , zc16   , zc32   , zc40,
               zcc01  , zcc14  , zcc16  , zcc32  ,
               zh01   , zh04   , zh14   , zh16   , zh28   , zh32   , zh40  , zhf,
               zhm01  , zhm04  , zhm14  , zhm16  , zhm28  , zhm32 , zhm40  , zmho  , zsho  , zsht;
        int i;

        double mt[12]  = {0, 48,0,4,16,28,32,40,1,49,14,17};
        int mn1 = 5;
        double zn1[6]  = {0, 120.0,110.0,100.0,90.0,72.5};

	double dgtr    = 1.74533e-2;
	double dr      = 1.72142e-2;
	double alpha[10]= {0, -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0};
	double altl[9] = {0, 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0};

// treat as local, but save in between calls
        static double alast = -999.0;

//------------------------------ begin ---------------------------------
        // ----  test for changed input
        v2 = vtst(msis00r.csw, iyd,sec,glat,glong,stl,f107a,f107,ap,2);

        yrd = iyd ;
        msis00r.gts3c.za  = msis00r.parm.pdl[16][2];
        zn1[1] = msis00r.gts3c.za;
        for (i = 1; i <= 9; i++)
            d[i] = 0.0;

        // ----  tinf variations not important below za or zn1[1]
        if (alt > zn1[1])
          {
            if (v2 == 1.0 | alast <= zn1[1])
                tinf = msis00r.lower.ptm[1]*msis00r.parm.pt[1] *
                    ( 1.0+msis00r.csw.sw[16] *
                       globe7(msis00r.csw, lpoly, yrd,sec,glat,glong,stl,f107a,f107,
                              ap,msis00r.parm.pt) );
          }
          else
            tinf = msis00r.lower.ptm[1]*msis00r.parm.pt[1];

        t[1] = tinf;
        // ----    gradient variations not important below zn1[5]
        if (alt > zn1[5])
          {
            if (v2 == 1 | alast <= zn1[5])
                msis00r.gts3c.g0 = msis00r.lower.ptm[4]*msis00r.parm.ps[1]*(1.0+msis00r.csw.sw[19]*
                     globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,
                            ap,msis00r.parm.ps));
          }
          else
            msis00r.gts3c.g0 = msis00r.lower.ptm[4]*msis00r.parm.ps[1];

        //  calculate these temperatures only if input changed
        if (v2 == 1.0  |  alt < 300.0)
            msis00r.gts3c.tlb = msis00r.lower.ptm[2]*
                  (1.0 + msis00r.csw.sw[17] *
                   globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,
                          ap,msis00r.parm.pddav[4]) ) * msis00r.parm.pd[1][4];

        msis00r.gts3c.s = msis00r.gts3c.g0/(tinf-msis00r.gts3c.tlb);
        // ---- lower thermosphere temp variations not significant for
        // ----  density above 300 km
        if (alt < 300.0)
          {
            if (v2 == 1.0 | alast >= 300.0)
              {
                msis00r.meso.tn1[2] = msis00r.lower.ptm[7]*msis00r.parm.ptl[1][1] /
                    (1.0 - msis00r.csw.sw[18]*glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[1]));
                msis00r.meso.tn1[3] = msis00r.lower.ptm[3]*msis00r.parm.ptl[1][2] /
                    (1.0 - msis00r.csw.sw[18]*glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[2]));
                msis00r.meso.tn1[4] = msis00r.lower.ptm[8]*msis00r.parm.ptl[1][3] /
                    (1.0 - msis00r.csw.sw[18]*glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[3]));
                msis00r.meso.tn1[5] = msis00r.lower.ptm[5]*msis00r.parm.ptl[1][4] /
                    (1.0 - msis00r.csw.sw[18]*msis00r.csw.sw[20]*
                     glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[4]));
                msis00r.meso.tgn1[2] = msis00r.lower.ptm[9]*msis00r.parm.pma[1][9] *
                    (1.0 + msis00r.csw.sw[18]*msis00r.csw.sw[20]*
                     glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[9])) *
                     msis00r.meso.tn1[5]*msis00r.meso.tn1[5] /
                     pow( (msis00r.lower.ptm[5]*msis00r.parm.ptl[1][4]),2 );
              }
          }
          else
          {
            msis00r.meso.tn1[2]  = msis00r.lower.ptm[7]*msis00r.parm.ptl[1][1];
            msis00r.meso.tn1[3]  = msis00r.lower.ptm[3]*msis00r.parm.ptl[1][2];
            msis00r.meso.tn1[4]  = msis00r.lower.ptm[8]*msis00r.parm.ptl[1][3];
            msis00r.meso.tn1[5]  = msis00r.lower.ptm[5]*msis00r.parm.ptl[1][4];
            msis00r.meso.tgn1[2] = msis00r.lower.ptm[9]*msis00r.parm.pma[1][9]*msis00r.meso.tn1[5]*
                                    msis00r.meso.tn1[5] /
                                    pow( (msis00r.lower.ptm[5]*msis00r.parm.ptl[1][4]),2 );
          }

        msis00r.gts3c.z0   = zn1[4];
        msis00r.gts3c.t0   = msis00r.meso.tn1[4];
        msis00r.gts3c.tr12 = 1.0;

        if (mass == 0)
            goto fifty;

        // ---- n2 variation factor at zlb
        g28 = msis00r.csw.sw[21] *
              globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[3]);
        lpoly.day = int(yrd) % 1000;

        // ----  variation of turbopause height
        zhf = msis00r.parm.pdl[25][2]*(1.0+msis00r.csw.sw[5]*msis00r.parm.pdl[25][1]*sin(dgtr*glat)*
              cos(dr*(lpoly.day-msis00r.parm.pt[14])));
        yrd  = iyd;
        t[1] = tinf;
        xmm = msis00r.lower.pdm[5][3];
        z   = alt;

        for(i = 1; i <= 11; i++)
            if (mass == mt[i])   goto fifteen;

        printf("mass %i5 not valid\n", mass);
        goto ninety;
   fifteen:
        if (z <= altl[6] | mass == 28 | mass == 48)
          {
            // ---- **** n2 density ****
            //       diffusive density at zlb
            msis00r.gts3c.db28 = msis00r.lower.pdm[1][3]*exp(g28)*msis00r.parm.pd[1][3];
            //       diffusive density at alt
            d[3] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db28,tinf,msis00r.gts3c.tlb,
                         28.0,alpha[3],t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,
                         mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
            msis00r.gts3c.dd = d[3];
            //       turbopause
            zh28  = msis00r.lower.pdm[3][3]*zhf;
            zhm28 = msis00r.lower.pdm[4][3]*msis00r.parm.pdl[6][2];
            xmd   = 28.0 - xmm;
            //       mixed density at zlb
            b28   = densu(msis00r.parmb, lsqv, zh28,msis00r.gts3c.db28,tinf,msis00r.gts3c.tlb,
                          xmd,alpha[3]-1.0,tz,msis00r.lower.ptm[6],
                          msis00r.gts3c.s,mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);

            if (z <= altl[3] && msis00r.csw.sw[15] != 0)
              {
                //       mixed density at alt
                msis00r.dmix.dm28 = densu(msis00r.parmb, lsqv, z,b28,tinf,msis00r.gts3c.tlb,
                                    xmm,alpha[3],tz,msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                                    zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
                //       net density at alt
                d[3] = dnet(d[3],msis00r.dmix.dm28,zhm28,xmm,28.0);
              }
          }
        switch(i)
          {
          case 1 : goto twenty;
          case 2 : goto fifty;
          case 3 : goto twenty;
          case 4 : goto twentyfive;
          case 5 : goto ninety;
          case 6 : goto thirtyfive;
          case 7 : goto fourty;
          case 8 : goto fourtyfive;
          case 9 : goto twentyfive;
          case 10 : goto fourtyeight;
          case 11 : goto fourtysix;
          }
   twenty:
        // ---- **** he density ****
        // ---- density variation factor at zlb
        g4 = msis00r.csw.sw[21] *
             globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[1]);

        //       diffusive density at zlb
        msis00r.gts3c.db04 = msis00r.lower.pdm[1][1]*exp(g4)*msis00r.parm.pd[1][1];

        //      diffusive density at alt
        d[1] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db04,tinf,msis00r.gts3c.tlb, 4.0,alpha[1],
                     t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        msis00r.gts3c.dd = d[1];
        if (z <= altl[1] && msis00r.csw.sw[15] != 0)
          {
            //      turbopause
            zh04 = msis00r.lower.pdm[3][1];

            //      mixed density at zlb
            b04 = densu(msis00r.parmb, lsqv, zh04,msis00r.gts3c.db04,tinf,msis00r.gts3c.tlb,
                        4.0-xmm,alpha[1]-1.0, t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,
                        mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);

            //      mixed density at alt
            msis00r.dmix.dm04 = densu(msis00r.parmb, lsqv, z,b04,tinf,msis00r.gts3c.tlb,
                                xmm,0.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                                zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
            zhm04 = zhm28;

            //      net density at alt
            d[1] = dnet(d[1],msis00r.dmix.dm04,zhm04,xmm,4.0);

            //      correction to specified mixing ratio at ground
            msis00r.gts3c.rl   = log(b28*msis00r.lower.pdm[2][1]/b04);
            zc04 = msis00r.lower.pdm[5][1]*msis00r.parm.pdl[1][2];
            hc04 = msis00r.lower.pdm[6][1]*msis00r.parm.pdl[2][2];

            //      net density corrected at alt
            d[1] = d[1]*ccor(z,msis00r.gts3c.rl,hc04,zc04);
          }

      if (mass != 48)   goto ninety;

   twentyfive:
        // ----**** o density ****
        // ---- density variation factor at zlb
        g16 = msis00r.csw.sw[21] *
              globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,
              msis00r.parm.pddav[2]);
        //       diffusive density at zlb
        msis00r.gts3c.db16 =  msis00r.lower.pdm[1][2]*exp(g16)*msis00r.parm.pd[1][2];
        // ---- diffusive density at alt
        d[2] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db16,tinf,msis00r.gts3c.tlb,
                     16.0,alpha[2],t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        msis00r.gts3c.dd = d[2];

        if (z <= altl[2] && msis00r.csw.sw[15] != 0)
          {
            //  corrected from msis00r.lower.pdm[3,1) to msis00r.lower.pdm[3][2]  12/2/85
            // ---- turbopause
            zh16 = msis00r.lower.pdm[3][2];
            //       mixed density at zlb
            b16 = densu(msis00r.parmb, lsqv, zh16,msis00r.gts3c.db16,tinf,msis00r.gts3c.tlb,
                        16.0-xmm,alpha[2]-1.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,
                        mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
            //       mixed density at alt
            msis00r.dmix.dm16 = densu(msis00r.parmb, lsqv, z,b16,tinf,msis00r.gts3c.tlb,xmm,0.0,
                                t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,msis00r.meso.tn1,
                                msis00r.meso.tgn1);
            zhm16 = zhm28;
            //       net density at alt
            d[2] = dnet(d[2],msis00r.dmix.dm16,zhm16,xmm,16.0);
            //   3/16/99 change form to match o2 departure from diff equil near 150
            //   km and add dependence on f10.7
            //       rl = dlog(b28*msis00r.lower.pdm[2][2]*fabs(msis00r.parm.pdl[17][2])/b16]
            msis00r.gts3c.rl = msis00r.lower.pdm[2][2]*msis00r.parm.pdl[17][2]*
                              (1.0+msis00r.csw.sw[1]*msis00r.parm.pdl[24][1]*(f107a-150.0));
            hc16 = msis00r.lower.pdm[6][2]*msis00r.parm.pdl[4][2];
            zc16 = msis00r.lower.pdm[5][2]*msis00r.parm.pdl[3][2];
            hc216 = msis00r.lower.pdm[6][2] * msis00r.parm.pdl[5][2];
            d[2] = d[2]*ccor2(z,msis00r.gts3c.rl,hc16,zc16,hc216);

            // ---- chemistry correction
            hcc16 = msis00r.lower.pdm[8][2]*msis00r.parm.pdl[14][2];
            zcc16 = msis00r.lower.pdm[7][2]*msis00r.parm.pdl[13][2];
            rc16  = msis00r.lower.pdm[4][2]*msis00r.parm.pdl[15][2];
            //       net density corrected at alt
            d[2]  = d[2]*ccor(z,rc16,hcc16,zcc16);
          }
        if (mass != 48 && mass != 49) goto ninety;

   thirtyfive:
        // ---- **** o2 density ****
        // ---- density variation factor at zlb
        g32 = msis00r.csw.sw[21] *
              globe7(msis00r.csw, lpoly, yrd,sec,glat,glong,stl,f107a,f107,ap,
                     msis00r.parm.pddav[5]);

        //       diffusive density at zlb
        msis00r.gts3c.db32 = msis00r.lower.pdm[1][4] * exp(g32)*msis00r.parm.pd[1][5];

        // ---- diffusive density at alt
        d[4] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db32,tinf,msis00r.gts3c.tlb,
                     32.0,alpha[4],t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        if (mass == 49)
           msis00r.gts3c.dd = msis00r.gts3c.dd+2.0*d[4];
        else
           msis00r.gts3c.dd = d[4];

        if (msis00r.csw.sw[15] != 0)
          {
            if (z <= altl[4])
              {
                // ---- turbopause
                zh32 = msis00r.lower.pdm[3][4];
                //       mixed density at zlb
                b32 = densu(msis00r.parmb, lsqv, zh32,msis00r.gts3c.db32,tinf,msis00r.gts3c.tlb,
                            32.0-xmm,alpha[4]-1.0, t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,
                            mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
                //       mixed density at alt
                msis00r.dmix.dm32 = densu(msis00r.parmb, lsqv, z,b32,tinf,msis00r.gts3c.tlb,xmm,0.0,t[2],
                                          msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,msis00r.meso.tn1,
                                          msis00r.meso.tgn1);
                zhm32 = zhm28;
                //       net density at alt
                d[4] = dnet(d[4],msis00r.dmix.dm32,zhm32,xmm,32.0);

                // ---- correction to specified mixing ratio at ground
                msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][4]/b32);
                hc32 = msis00r.lower.pdm[6][4]*msis00r.parm.pdl[8][2];
                zc32 = msis00r.lower.pdm[5][4]*msis00r.parm.pdl[7][2];
                d[4] = d[4]*ccor(z,msis00r.gts3c.rl,hc32,zc32);
              }

            //       correction for general departure from diffusive equilibrium above zlb
            hcc32  = msis00r.lower.pdm[8][4]*msis00r.parm.pdl[23][2];
            hcc232 = msis00r.lower.pdm[8][4]*msis00r.parm.pdl[23][1];
            zcc32  = msis00r.lower.pdm[7][4]*msis00r.parm.pdl[22][2];
            rc32   = msis00r.lower.pdm[4][4]*msis00r.parm.pdl[24][2]*(
                      1.0+msis00r.csw.sw[1]*msis00r.parm.pdl[24][1]*(f107a-150.0));

            //       net density corrected at alt
            d[4] = d[4]*ccor2(z,rc32,hcc32,zcc32,hcc232);
          }

        if (mass != 48)   goto ninety;

   fourty:
        // ---- **** ar density ****
        // ---- density variation factor at zlb
        g40 =  msis00r.csw.sw[21]*
               globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[6]);
        //      diffusive density at zlb
        msis00r.gts3c.db40 = msis00r.lower.pdm[1][5]*exp(g40)*msis00r.parm.pd[1][6];

        // ---- diffusive density at alt
        d[5] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db40,tinf,msis00r.gts3c.tlb,
                     40.0,alpha[5],t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        msis00r.gts3c.dd = d[5];

        if (z <= altl[5] && msis00r.csw.sw[15] != 0)
          {
            // ---- turbopause
            zh40 = msis00r.lower.pdm[3][5];
            //       mixed density at zlb
            b40 = densu(msis00r.parmb, lsqv, zh40,msis00r.gts3c.db40,tinf,msis00r.gts3c.tlb,
                        40.0-xmm,alpha[5]-1.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,
                        msis00r.meso.tn1,msis00r.meso.tgn1);
            //       mixed density at alt
            msis00r.dmix.dm40  = densu(msis00r.parmb, lsqv, z,b40,tinf,msis00r.gts3c.tlb,xmm,
                                       0.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,
                                       msis00r.meso.tn1,msis00r.meso.tgn1);
            zhm40 = zhm28;
            //       net density at alt
            d[5] = dnet(d[5],msis00r.dmix.dm40,zhm40,xmm,40.0);

            // ---- correction to specified mixing ratio at ground
            msis00r.gts3c.rl   = log(b28*msis00r.lower.pdm[2][5]/b40);
            hc40 = msis00r.lower.pdm[6][5]*msis00r.parm.pdl[10][2];
            zc40 = msis00r.lower.pdm[5][5]*msis00r.parm.pdl[9][2];

            //       net density corrected at alt
            d[5] = d[5]*ccor(z,msis00r.gts3c.rl,hc40,zc40);
          }

        if (mass != 48)   goto ninety;

   fourtyfive:
        // ----  **** hydrogen density ****
        // ---- density variation factor at zlb
        g1 = msis00r.csw.sw[21]*
             globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[7]);
        //      diffusive density at zlb
        msis00r.gts3c.db01 = msis00r.lower.pdm[1][6]*exp(g1)*msis00r.parm.pd[1][7];

        // ---- diffusive density at alt
        d[7] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db01,tinf,msis00r.gts3c.tlb,
                     1.0,alpha[7],t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        msis00r.gts3c.dd = d[7];

        if (z <= altl[7] && msis00r.csw.sw[15] != 0)
          {
            // ---- turbopause
            zh01 = msis00r.lower.pdm[3][6];
            //       mixed density at zlb
            b01 = densu(msis00r.parmb, lsqv, zh01,msis00r.gts3c.db01,tinf,msis00r.gts3c.tlb,
                        1.0-xmm,alpha[7]-1.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,
                        msis00r.meso.tn1,msis00r.meso.tgn1);
            //       mixed density at alt
            msis00r.dmix.dm01 = densu(msis00r.parmb, lsqv, z,b01,tinf,msis00r.gts3c.tlb,xmm,
                                0.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,
                                mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
            zhm01 = zhm28;

            //       net density at alt
            d[7] = dnet(d[7],msis00r.dmix.dm01,zhm01,xmm,1.0);

            // ---- correction to specified mixing ratio at ground
            msis00r.gts3c.rl   = log(b28*msis00r.lower.pdm[2][6]*fabs(msis00r.parm.pdl[18][2])/b01);
            hc01 = msis00r.lower.pdm[6][6]*msis00r.parm.pdl[12][2];
            zc01 = msis00r.lower.pdm[5][6]*msis00r.parm.pdl[11][2];
            d[7] = d[7]*ccor(z,msis00r.gts3c.rl,hc01,zc01);

            // ---- chemistry correction
            hcc01 = msis00r.lower.pdm[8][6]*msis00r.parm.pdl[20][2];
            zcc01 = msis00r.lower.pdm[7][6]*msis00r.parm.pdl[19][2];
            rc01  = msis00r.lower.pdm[4][6]*msis00r.parm.pdl[21][2];

            //      net density corrected at alt
            d[7] = d[7]*ccor(z,rc01,hcc01,zcc01);
          }

        if (mass != 48)   goto ninety;

   fourtyeight:
        // ----  **** atomic nitrogen density ****
        // ---- density variation factor at zlb
        g14 = msis00r.csw.sw[21] *
              globe7(msis00r.csw, lpoly,yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[8]);

        //       diffusive density at zlb
        msis00r.gts3c.db14 = msis00r.lower.pdm[1][7]*exp(g14)*msis00r.parm.pd[1][8];

        // ---- diffusive density at alt
        d[8] = densu(msis00r.parmb, lsqv, z,msis00r.gts3c.db14,tinf,msis00r.gts3c.tlb,14.0,alpha[8],
                     t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                     zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        msis00r.gts3c.dd = d[8];

        if (z <= altl[8] && msis00r.csw.sw[15] != 0)
          {
            // ---- turbopause
            zh14 = msis00r.lower.pdm[3][7];

            //       mixed density at zlb
            b14 = densu(msis00r.parmb, lsqv, zh14,msis00r.gts3c.db14,tinf,msis00r.gts3c.tlb,
                        14.0-xmm,alpha[8]-1.0,t[2],msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,zn1,
                        msis00r.meso.tn1,msis00r.meso.tgn1);

            //       mixed density at alt
            msis00r.dmix.dm14 = densu(msis00r.parmb, lsqv, z,b14,tinf,msis00r.gts3c.tlb,xmm,0.0,
                                t[2],msis00r.lower.ptm[6],
                                msis00r.gts3c.s,mn1,zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
            zhm14 = zhm28;
            //       net density at alt
            d[8] = dnet(d[8],msis00r.dmix.dm14,zhm14,xmm,14.0);

            // ---- correction to specified mixing ratio at ground
            msis00r.gts3c.rl   = log(b28*msis00r.lower.pdm[2][7]*fabs(msis00r.parm.pdl[3][1])/b14);
            hc14 = msis00r.lower.pdm[6][7]*msis00r.parm.pdl[2][1];
            zc14 = msis00r.lower.pdm[5][7]*msis00r.parm.pdl[1][1];
            d[8] = d[8]*ccor(z,msis00r.gts3c.rl,hc14,zc14);

            // ---- chemistry correction
            hcc14 = msis00r.lower.pdm[8][7]*msis00r.parm.pdl[5][1];
            zcc14 = msis00r.lower.pdm[7][7]*msis00r.parm.pdl[4][1];
            rc14  = msis00r.lower.pdm[4][7]*msis00r.parm.pdl[6][1];

            //       net density corrected at alt
            d[8] = d[8]*ccor(z,rc14,hcc14,zcc14);
          }

        if (mass != 48) goto ninety;

   fourtysix:
        // ----  **** anomalous oxygen density ****
        g16h = msis00r.csw.sw[21] *
               globe7(msis00r.csw, lpoly, yrd,sec,glat,glong,stl,f107a,f107,ap,msis00r.parm.pddav[9]);
        db16h = msis00r.lower.pdm[1][8]*exp(g16h)*msis00r.parm.pd[1][9];
        tho = msis00r.lower.pdm[10][8]*msis00r.parm.pdl[7][1];
        msis00r.gts3c.dd  = densu(msis00r.parmb, lsqv, z,db16h,tho,tho,16.0,alpha[9],t2,
                            msis00r.lower.ptm[6],msis00r.gts3c.s,mn1,
                            zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
        zsht = msis00r.lower.pdm[6][8];
        zmho = msis00r.lower.pdm[5][8];
        zsho = scalh(msis00r.parmb, zmho,16.0,tho);
        d[9] = msis00r.gts3c.dd * exp( -zsht/zsho*( exp(-(z-zmho)/zsht )-1.0) );

        if (mass != 48) goto ninety;

        // ---- total mass density
        d[6] = 1.66e-24 *
              (4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] + 40.0*d[5] + d[7] + 14.0*d[8]);
        msis00r.gts3c.db48 = 1.66e-24*(4.0*msis00r.gts3c.db04+16.0*msis00r.gts3c.db16+
                             28.0*msis00r.gts3c.db28+32.0*msis00r.gts3c.db32+
                             40.0*msis00r.gts3c.db40+msis00r.gts3c.db01+14.0*msis00r.gts3c.db14);

        goto ninety;
        // ---- temperature at altitude
   fifty:
        z = fabs(alt);
//cdav
//        ddum is never used so this call is unnecessary
//        ddum  = densu(msis00r.parmb, lsqv, z,1.0, tinf,msis00r.gts3c.tlb,0.0,0.0,t[2],
//                  msis00r.lower.ptm[6],msis00r.gts3c.msis00r.gts3c.s,mn1,
//                  zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
   ninety:
        // ---- adjust densities from cgs to kgm
        if (msis00r.metsel.imr == 1)
          {
            for (i = 1; i <= 9; i++)
                d[i] = d[i]*1.0e6;

            d[6] = d[6]/1000.0;
          }  

      alast = alt;
      }

/* ----------------------------------------------------------------------
*      calculate scale height (km)
* ----------------------------------------------------------------------      */

double scalh
       (
         parmbtype& parmb,
         double alt, double xm, double temp
       )
        {
        double g;
        double rgas = 831.4;

        g     = parmb.gsurf / pow( (1.0 + alt/parmb.re),2 );
        return rgas*temp /(g*xm);
        }

// cdav these next few functions could be done as inline functions if that speeds
// things up
// ---- 3hr magnetic activity functions
//      eq. a24d
double g0
       (
         double a, double p[151]
       )
          {
          return (a - 4.0 + (p[26]-1.0) *
                 (a - 4.0 + (exp(-fabs(p[25])*(a-4.0)) -1.0) / fabs(p[25]) ));
          }
// ---- eq. a24c
double sumex
       (
         double ex
       )
          {
          return 1.0 + (1.0-pow(ex,19)) / (1.0-ex)*sqrt(ex);
          }
// ---- eq. a24a
double sg0
       (
         double ex, double p[151], double ap[8]
       )
          {
          return ( g0(ap[2],p) +
                  ( g0(ap[3],p)*ex + g0(ap[4],p)*ex*ex + g0(ap[5],p)* pow(ex,3) +
                   ( g0(ap[6],p)*pow(ex,4) + g0(ap[7],p)*pow(ex,12)) * (1.0-pow(ex,8))
                   / (1.0-ex)
                  )
                 ) / sumex(ex);
          }
/* ----------------------------------------------------------------------
*         calculate g(l) function
*         upper thermosphere parameters
* ----------------------------------------------------------------------      */

double globe7
       (
         cswtype& csw,
         lpolytype& lpoly,
         double yrd, double sec, double lat, double llong, double tloc, double f107a,
         double f107, double ap[8], double p[151]
       )
        {
        int i, iyr;
        double c, s, c2, c4, s2, cd14, cd18, cd32, cd39, f1, f2, t71, t72,
               t81, t82, p44, p45, exp1;
        double a, ex, t[16], tinfg;

        double dgtr = 1.74533e-2;
        double dr   = 1.72142e-2;
        double xl   = 1000.0;
        static double tll  = 1000.0;
        static double sw9  = 1.0;
        double dayl = -1.0;
        double p14  = -1000.0;
        double p18  = -1000.0;
        double p32  = -1000.0;
        double hr   = 0.2618;
        double sr   = 7.2722e-5;
        int sv[26] = {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        double nsw = 14;
        double p39 = -1000.0;

//------------------------------ begin ---------------------------------
        if (csw.isw != 64999)
            tselec(csw, sv);

        for (i = 1; i <= 15; i++)
            t[i] = 0.0;

        if (csw.sw[9] > 0)
            sw9 = 1.0;
        if (csw.sw[9] < 0)
            sw9 = -1.0;
        iyr = yrd/1000.0;
        lpoly.day = yrd - iyr*1000.0;
        lpoly.xlong = llong;
        //  eq. a22 (remainder of code)
        if (xl != lat)
          {
            // ----    calculate legendre polynomials
            c  = sin(lat*dgtr);
            s  = cos(lat*dgtr);
            c2 = c*c;
            c4 = c2*c2;
            s2 = s*s;
            lpoly.plg[2][1] = c;
            lpoly.plg[3][1] = 0.5*(3.0*c2 -1.0);
            lpoly.plg[4][1] = 0.5*(5.0*c*c2-3.0*c);
            lpoly.plg[5][1] = (35.0*c4 - 30.0*c2 + 3.0)/8.0;
            lpoly.plg[6][1] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0;
            lpoly.plg[7][1] = (11.0*c*lpoly.plg[6][1] - 5.0*lpoly.plg[5][1])/6.0;
    //        lpoly.plg[8][1] = (13.0*c*lpoly.plg[7][1] - 6.0*lpoly.plg[6][1])/7.0;
            lpoly.plg[2][2] = s;
            lpoly.plg[3][2] = 3.0*c*s;
            lpoly.plg[4][2] = 1.5*(5.0*c2-1.0)*s;
            lpoly.plg[5][2] = 2.5*(7.0*c2*c-3.0*c)*s;
            lpoly.plg[6][2] = 1.875*(21.0*c4 - 14.0*c2 +1.0)*s;
            lpoly.plg[7][2] = (11.0*c*lpoly.plg[6][2]-6.0*lpoly.plg[5][2])/5.0;
    //        lpoly.plg[8][2] = (13.0*c*lpoly.plg[7][2]-7.0*lpoly.plg[6][2])/6.0;
    //        lpoly.plg[9][2] = (15.0*c*lpoly.plg[8][2]-8.0*lpoly.plg[7][2])/7.0;
            lpoly.plg[3][3] = 3.0*s2;
            lpoly.plg[4][3] = 15.0*s2*c;
            lpoly.plg[5][3] = 7.5*(7.0*c2 -1.0)*s2;
            lpoly.plg[6][3] = 3.0*c*lpoly.plg[5][3]-2.0*lpoly.plg[4][3];
            lpoly.plg[7][3] = (11.0*c*lpoly.plg[6][3]-7.0*lpoly.plg[5][3])/4.0;
            lpoly.plg[8][3] = (13.0*c*lpoly.plg[7][3]-8.0*lpoly.plg[6][3])/5.0;
            lpoly.plg[4][4] = 15.0*s2*s;
            lpoly.plg[5][4] = 105.0*s2*s*c;
            lpoly.plg[6][4] = (9.0*c*lpoly.plg[5][4]-7.0*lpoly.plg[4][4])/2.0;
            lpoly.plg[7][4] = (11.0*c*lpoly.plg[6][4]-8.0*lpoly.plg[5][4])/3.0;

            xl = lat;
          }
        if (tll == tloc) goto sixteen;
        if (csw.sw[7] == 0 && csw.sw[8] == 0 && csw.sw[14] == 0) goto sixteen;
        lpoly.stloc  = sin(hr*tloc);
        lpoly.ctloc  = cos(hr*tloc);
        lpoly.s2tloc = sin(2.0*hr*tloc);
        lpoly.c2tloc = cos(2.0*hr*tloc);
        lpoly.s3tloc = sin(3.0*hr*tloc);
        lpoly.c3tloc = cos(3.0*hr*tloc);
        tll = tloc;
   sixteen:
        if (lpoly.day != dayl | p[14] != p14)
            cd14 = cos(dr*(lpoly.day-p[14]));
        if (lpoly.day != dayl | p[18] != p18)
            cd18 = cos(2.0*dr*(lpoly.day-p[18]));

        if (lpoly.day != dayl | p[32] != p32)
            cd32 = cos(dr*(lpoly.day-p[32]));
        if (lpoly.day != dayl | p[39] != p39)
            cd39 = cos(2.0*dr*(lpoly.day-p[39]));

        dayl = lpoly.day;
        p14 = p[14];
        p18 = p[18];
        p32 = p[32];
        p39 = p[39];

        // ----   f10.7 effect
        lpoly.df  = f107  - f107a;
        lpoly.dfa = f107a - 150.0;
        t[1] = p[20]*lpoly.df*(1.0+p[60]*lpoly.dfa) + p[21]*lpoly.df*lpoly.df +
                     p[22]*lpoly.dfa + p[30]*lpoly.dfa*lpoly.dfa;
        f1 = 1.0 + (p[48]*lpoly.dfa + p[20]*lpoly.df+p[21]*lpoly.df*lpoly.df)*csw.swc[1];
        f2 = 1.0 + (p[50]*lpoly.dfa + p[20]*lpoly.df+p[21]*lpoly.df*lpoly.df)*csw.swc[1];

        // ----  time independent
        t[2] = (p[2]*lpoly.plg[3][1] + p[3]*lpoly.plg[5][1]+p[23]*lpoly.plg[7][1])
               +(p[15]*lpoly.plg[3][1])*lpoly.dfa*csw.swc[1] +p[27]*lpoly.plg[2][1];

        // ----  symmetrical annual
        t[3] = (p[19] )*cd32;

        // ----  symmetrical semiannual
        t[4] = (p[16]+p[17]*lpoly.plg[3][1])*cd18;

        // ----  asymmetrical annual
        t[5] =  f1* (p[10]*lpoly.plg[2][1]+p[11]*lpoly.plg[4][1])*cd14;

        // ----   asymmetrical semiannual
        t[6] =    p[38]*lpoly.plg[2][1]*cd39;

        // ----  diurnal
        if (csw.sw[7] != 0)
          {
            t71 = (p[12]*lpoly.plg[3][2])*cd14*csw.swc[5];
            t72 = (p[13]*lpoly.plg[3][2])*cd14*csw.swc[5];
            t[7] = f2* ((p[4]*lpoly.plg[2][2] + p[5]*lpoly.plg[4][2] + p[28]*lpoly.plg[6][2]
                   + t71)*lpoly.ctloc + (p[7]*lpoly.plg[2][2] + p[8]*lpoly.plg[4][2]
                   + p[29]*lpoly.plg[6][2] + t72)*lpoly.stloc);
          }
        // ----  semidiurnal
        if (csw.sw[8] != 0)
          {
            t81  = (p[24]*lpoly.plg[4][3]+p[36]*lpoly.plg[6][3])*cd14*csw.swc[5];
            t82  = (p[34]*lpoly.plg[4][3]+p[37]*lpoly.plg[6][3])*cd14*csw.swc[5];
            t[8] = f2*((p[6]*lpoly.plg[3][3] + p[42]*lpoly.plg[5][3] + t81)*lpoly.c2tloc
                   +(p[9]*lpoly.plg[3][3] + p[43]*lpoly.plg[5][3] + t82)*lpoly.s2tloc);
          }
        // ----  terdiurnal
        if (csw.sw[14] != 0)
          {
            t[14] = f2*
                   ((p[40]*lpoly.plg[4][4] + (p[94]*lpoly.plg[5][4]
                     + p[47]*lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.s3tloc
                     + (p[41]*lpoly.plg[4][4] + (p[95]*lpoly.plg[5][4]
                     + p[49]*lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.c3tloc);
          }

        // ----    magnetic activity based on daily ap
        if (sw9 != -1.0)
          {
            lpoly.apd = (ap[1]-4.0);
            p44 = p[44];
            p45 = p[45];
// cdav at one time i saw an error where p44 was less than zero - unfortunately, i didn't
//      write exactly how i got there. it may have been during debugging when i had other
//      things wrong
            if (p44 <= 0.0)
                p44 = 1.0e-5;
            lpoly.apdf = lpoly.apd+(p45-1.0)*(lpoly.apd+(exp(-p44  *lpoly.apd)-1.0)/p44);
            if (csw.sw[9] == 0)
                goto fourty;
            t[9] = lpoly.apdf*(p[33]+p[46]*lpoly.plg[3][1]+p[35]*lpoly.plg[5][1]+
                (p[101]*lpoly.plg[2][1]+p[102]*lpoly.plg[4][1]+p[103]*lpoly.plg[6][1])*cd14*csw.swc[5]+
                (p[122]*lpoly.plg[2][2]+p[123]*lpoly.plg[4][2]+p[124]*lpoly.plg[6][2])*csw.swc[7]*
                cos(hr*(tloc-p[125])));
            goto fourty;
          }

        if (p[52] == 0)
            goto fourty;
        exp1 = exp(-10800.0*fabs(p[52])/(1.0+p[139]* (45.0-fabs(lat))));
        if (exp1 > .99999)
            exp1 = 0.99999;
        if (p[25] < 1.0e-4)
            p[25] = 1.0e-4;
        lpoly.apt[1] = sg0(exp1, p, ap);
//        apt[2] = sg2[exp1]
//        apt[3] = sg0[exp2]
//        apt[4] = sg2[exp2]

        if (csw.sw[9] == 0) goto fourty;
        t[9] = lpoly.apt[1]*(p[51]+p[97]*lpoly.plg[3][1]+p[55]*lpoly.plg[5][1]+
              (p[126]*lpoly.plg[2][1]+p[127]*lpoly.plg[4][1]+p[128]*lpoly.plg[6][1])*cd14*csw.swc[5]+
              (p[129]*lpoly.plg[2][2]+p[130]*lpoly.plg[4][2]+p[131]*lpoly.plg[6][2])*csw.swc[7]*
               cos(hr*(tloc-p[132])));
  fourty:

        if (csw.sw[10] == 0 | llong <= -1000.0)
            goto fourtynine;

        // ----  longitudinal
        if (csw.sw[11] != 0)
          {
            t[11] = (1.0+p[81]*lpoly.dfa*csw.swc[1])*
                    ((p[65]*lpoly.plg[3][2]+p[66]*lpoly.plg[5][2]+p[67]*lpoly.plg[7][2]
                    +p[104]*lpoly.plg[2][2]+p[105]*lpoly.plg[4][2]+p[106]*lpoly.plg[6][2]
                    +csw.swc[5]*(p[110]*lpoly.plg[2][2]+p[111]*lpoly.plg[4][2]+p[112]*lpoly.plg[6][2])*cd14)*
                    cos(dgtr*llong)
                    +(p[91]*lpoly.plg[3][2]+p[92]*lpoly.plg[5][2]+p[93]*lpoly.plg[7][2]
                    +p[107]*lpoly.plg[2][2]+p[108]*lpoly.plg[4][2]+p[109]*lpoly.plg[6][2]
                    +csw.swc[5]*(p[113]*lpoly.plg[2][2]+p[114]*lpoly.plg[4][2]+p[115]*lpoly.plg[6][2])*cd14)*
                    sin(dgtr*llong));
          }

        // ----  ut and mixed ut,longitude
        if (csw.sw[12] != 0)
          {
            t[12] = (1.0+p[96]*lpoly.plg[2][1])*(1.0+p[82]*lpoly.dfa*csw.swc[1])*
                    (1.0+p[120]*lpoly.plg[2][1]*csw.swc[5]*cd14)*
                    ((p[69]*lpoly.plg[2][1]+p[70]*lpoly.plg[4][1]+p[71]*lpoly.plg[6][1])*
                    cos(sr*(sec-p[72])));
            t[12] = t[12]+csw.swc[11]*
                    (p[77]*lpoly.plg[4][3]+p[78]*lpoly.plg[6][3]+p[79]*lpoly.plg[8][3])*
                    cos(sr*(sec-p[80])+2.0*dgtr*llong)*(1.0+p[138]*lpoly.dfa*
                    csw.swc[1]);
          }

        // ----  ut,longitude magnetic activity
        if (csw.sw[13] == 0) goto fourtyeight;

        if (sw9 != -1.0)
          {
            t[13] =  lpoly.apdf*csw.swc[11]*(1.0+p[121]*lpoly.plg[2][1])*
                    ((p[ 61]*lpoly.plg[3][2]+p[ 62]*lpoly.plg[5][2]+p[ 63]*lpoly.plg[7][2])*
                    cos(dgtr*(llong-p[ 64])))
                    +lpoly.apdf*csw.swc[11]*csw.swc[5]*
                    (p[116]*lpoly.plg[2][2]+p[117]*lpoly.plg[4][2]+p[118]*lpoly.plg[6][2])*
                    cd14*cos(dgtr*(llong-p[119]))
                    + lpoly.apdf*csw.swc[12]*
                   (p[ 84]*lpoly.plg[2][1]+p[ 85]*lpoly.plg[4][1]+p[ 86]*lpoly.plg[6][1])*
                    cos(sr*(sec-p[ 76]));
            goto fourtyeight;
          }

        if (p[52] == 0) goto fourtyeight;
            t[13] = lpoly.apt[1]*csw.swc[11]*(1.0+p[133]*lpoly.plg[2][1])*
                   ((p[53]*lpoly.plg[3][2]+p[99]*lpoly.plg[5][2]+p[68]*lpoly.plg[7][2])*
                   cos(dgtr*(llong-p[98])))
                   +lpoly.apt[1]*csw.swc[11]*csw.swc[5]*
                   (p[134]*lpoly.plg[2][2]+p[135]*lpoly.plg[4][2]+p[136]*lpoly.plg[6][2])*
                    cd14*cos(dgtr*(llong-p[137]))
                   +lpoly.apt[1]*csw.swc[12]*
                   (p[56]*lpoly.plg[2][1]+p[57]*lpoly.plg[4][1]+p[58]*lpoly.plg[6][1])*
                   cos(sr*(sec-p[59]));
   fourtyeight:
//       parms not used: 83, 90,100,140-150
   fourtynine:
        tinfg = p[31];
        for (i = 1; i <= nsw; i++)
            tinfg = tinfg + fabs(csw.sw[i])*t[i];

        return tinfg;

        }
/* ----------------------------------------------------------------------
*      version of globe for lower atmosphere 10/26/99
* ----------------------------------------------------------------------      */

double glob7s
       (
         lpolytype& lpoly,
         cswtype& csw,
         double p[101]
       )
        {
        double cd32, cd18, cd14, cd39, t[15], t71, t72, t81, t82, tt;

        double dr   = 1.72142e-2;
        double dgtr = 1.74533e-2;
        double pset =  2.0;
        double dayl = -1.0;
        double p32 = -1000.0;
        double p18 = -1000.0;
        double p14 = -1000.0;
        double p39 = -1000.0;

        int i;
//------------------------------ begin ---------------------------------
        // ---- confirm parameter set
        if (p[100] == 0)
            p[100] = pset;
        if (p[100] != pset)
          {
            printf("wrong parameter set for glob7s %10.1f %10.1f \n", pset,p[100]);
            return 0;
          }

        for (i = 1; i <= 14; i++)
            t[i] = 0.0;

        if (lpoly.day != dayl | p32 != p[32])
            cd32 = cos(dr*(lpoly.day-p[32]));
        if (lpoly.day != dayl | p18 != p[18])
            cd18 = cos(2.0*dr*(lpoly.day-p[18]));

        if (lpoly.day != dayl | p14 != p[14])
            cd14 = cos(dr*(lpoly.day-p[14]));
        if (lpoly.day != dayl | p39 != p[39])
            cd39 = cos(2.0*dr*(lpoly.day-p[39]));

        dayl = lpoly.day;
        p32  = p[32];
        p18  = p[18];
        p14  = p[14];
        p39  = p[39];

        // ---- f10.7
        t[1] = p[22]*lpoly.dfa;

        // ---- time independent
        t[2] = p[2]*lpoly.plg[3][1]+p[3]*lpoly.plg[5][1]+p[23]*lpoly.plg[7][1]
               +p[27]*lpoly.plg[2][1]+p[15]*lpoly.plg[4][1]+p[60]*lpoly.plg[6][1];

        // ---- symmetrical annual
        t[3] = (p[19]+p[48]*lpoly.plg[3][1]+p[30]*lpoly.plg[5][1])*cd32;

        // ---- symmetrical semiannual
        t[4] = (p[16]+p[17]*lpoly.plg[3][1]+p[31]*lpoly.plg[5][1])*cd18;

        // ---- asymmetrical annual
        t[5] = (p[10]*lpoly.plg[2][1]+p[11]*lpoly.plg[4][1]+p[21]*lpoly.plg[6][1])*cd14;

        // ---- asymmetrical semiannual
        t[6] = (p[38]*lpoly.plg[2][1])*cd39;

        // ----  diurnal
        if (csw.sw[7] != 0)
          {
            t71  = p[12]*lpoly.plg[3][2]*cd14*csw.swc[5];
            t72  = p[13]*lpoly.plg[3][2]*cd14*csw.swc[5];
            t[7] = ((p[4]*lpoly.plg[2][2] + p[5]*lpoly.plg[4][2] + t71)*lpoly.ctloc
                   + (p[7]*lpoly.plg[2][2] + p[8]*lpoly.plg[4][2] + t72)*lpoly.stloc);
          }

        // ----  semidiurnal
        if (csw.sw[8] != 0)
          {
            t81  = (p[24]*lpoly.plg[4][3]+p[36]*lpoly.plg[6][3])*cd14*csw.swc[5];
            t82  = (p[34]*lpoly.plg[4][3]+p[37]*lpoly.plg[6][3])*cd14*csw.swc[5];
            t[8] = ((p[6]*lpoly.plg[3][3] + p[42]*lpoly.plg[5][3] + t81)*lpoly.c2tloc
                   +(p[9]*lpoly.plg[3][3] + p[43]*lpoly.plg[5][3] + t82)*lpoly.s2tloc);
          }

        // ----  terdiurnal
        if (csw.sw[14] != 0)
            t[14] = p[40]*lpoly.plg[4][4]*lpoly.s3tloc +p[41]*lpoly.plg[4][4]*lpoly.c3tloc;

        // ---- magnetic activity
        if (csw.sw[9] != 0)
          {
            if (csw.sw[9] == 1)
                t[9] = lpoly.apdf*(p[33]+p[46]*lpoly.plg[3][1]*csw.swc[2]);
            if (csw.sw[9] == -1)
                t[9] = (p[51]*lpoly.apt[1]+p[97]*lpoly.plg[3][1]*lpoly.apt[1]*csw.swc[2]);
          }
        if (csw.sw[10] == 0 | csw.sw[11] == 0 | lpoly.xlong <= -1000.0)
            goto fourtynine;

        // ----  longitudinal
        t[11] = (1.0+lpoly.plg[2][1]*(p[81]*csw.swc[5]*cos(dr*(lpoly.day-p[82]))
                +p[86]*csw.swc[6]*cos(2.0*dr*(lpoly.day-p[87])))
                +p[84]*csw.swc[3]*cos(dr*(lpoly.day-p[85]))
                +p[88]*csw.swc[4]*cos(2.0*dr*(lpoly.day-p[89])))
                *((p[65]*lpoly.plg[3][2]+p[66]*lpoly.plg[5][2]+p[67]*lpoly.plg[7][2]
                +p[75]*lpoly.plg[2][2]+p[76]*lpoly.plg[4][2]+p[77]*lpoly.plg[6][2])*cos(dgtr*lpoly.xlong)
                +(p[91]*lpoly.plg[3][2]+p[92]*lpoly.plg[5][2]+p[93]*lpoly.plg[7][2]
                +p[78]*lpoly.plg[2][2]+p[79]*lpoly.plg[4][2]+p[80]*lpoly.plg[6][2] )*sin(dgtr*lpoly.xlong));
   fourtynine:
        tt = 0.0;
        for (i = 1; i <= 14; i++)
            tt = tt + fabs(csw.sw[i])*t[i];

        return tt;
        }

/* ----------------------------------------------------------------------
* nrlmsise-00 01-feb-02 dist 17
*
* this function loads all the iniital data for the msis00 model
*
* ----------------------------------------------------------------------      */

void msis00init
     (
       msistype& msis00r
     )
        {

        int i,j;

// cdav i'm not string literate yet in c++ :-)
//        msis00r.datime.isdate = "01-feb-02";
//        msis00r.datime.istime = "15:49:27 ";
//        msis00r.datime.name   = "msise-00 ";

        msis00r.csw.isw = 0;
//
// cdav   set this to output in
        msis00r.metsel.imr = 1;

      // ----   temperature
       static double pt[151] =
        { 0,
        9.86573e-01, 1.62228e-02, 1.55270e-02,-1.04323e-01,-3.75801e-03,
       -1.18538e-03,-1.24043e-01, 4.56820e-03, 8.76018e-03,-1.36235e-01,
       -3.52427e-02, 8.84181e-03,-5.92127e-03,-8.61650e+00, 0.00000e+00,
        1.28492e-02, 0.00000e+00, 1.30096e+02, 1.04567e-02, 1.65686e-03,
       -5.53887e-06, 2.97810e-03, 0.00000e+00, 5.13122e-03, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,-7.27026e-06,
        0.00000e+00, 6.74494e+00, 4.93933e-03, 2.21656e-03, 2.50802e-03,
        0.00000e+00, 0.00000e+00,-2.08841e-02,-1.79873e+00, 1.45103e-03,
        2.81769e-04,-1.44703e-03,-5.16394e-05, 8.47001e-02, 1.70147e-01,
        5.72562e-03, 5.07493e-05, 4.36148e-03, 1.17863e-04, 4.74364e-03,
        6.61278e-03, 4.34292e-05, 1.44373e-03, 2.41470e-05, 2.84426e-03,
        8.56560e-04, 2.04028e-03, 0.00000e+00,-3.15994e+03,-2.46423e-03,
        1.13843e-03, 4.20512e-04, 0.00000e+00,-9.77214e+01, 6.77794e-03,
        5.27499e-03, 1.14936e-03, 0.00000e+00,-6.61311e-03,-1.84255e-02,
       -1.96259e-02, 2.98618e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        6.44574e+02, 8.84668e-04, 5.05066e-04, 0.00000e+00, 4.02881e+03,
       -1.89503e-03, 0.00000e+00, 0.00000e+00, 8.21407e-04, 2.06780e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -1.20410e-02,-3.63963e-03, 9.92070e-05,-1.15284e-04,-6.33059e-05,
       -6.05545e-01, 8.34218e-03,-9.13036e+01, 3.71042e-04, 0.00000e+00,
        4.19000e-04, 2.70928e-03, 3.31507e-03,-4.44508e-03,-4.96334e-03,
       -1.60449e-03, 3.95119e-03, 2.48924e-03, 5.09815e-04, 4.05302e-03,
        2.24076e-03, 0.00000e+00, 6.84256e-03, 4.66354e-04, 0.00000e+00,
       -3.68328e-04, 0.00000e+00, 0.00000e+00,-1.46870e+02, 0.00000e+00,
        0.00000e+00, 1.09501e-03, 4.65156e-04, 5.62583e-04, 3.21596e+00,
        6.43168e-04, 3.14860e-03, 3.40738e-03, 1.78481e-03, 9.62532e-04,
        5.58171e-04, 3.43731e+00,-2.33195e-01, 5.10289e-04, 0.00000e+00,
        0.00000e+00,-9.25347e+04, 0.00000e+00,-1.99639e-03, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        };
        // ----   he density
      static double pd[10][151] = {
        { 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        },
        { 0,
        1.09979e+00,-4.88060e-02,-1.97501e-01,-9.10280e-02,-6.96558e-03,
        2.42136e-02, 3.91333e-01,-7.20068e-03,-3.22718e-02, 1.41508e+00,
        1.68194e-01, 1.85282e-02, 1.09384e-01,-7.24282e+00, 0.00000e+00,
        2.96377e-01,-4.97210e-02, 1.04114e+02,-8.61108e-02,-7.29177e-04,
        1.48998e-06, 1.08629e-03, 0.00000e+00, 0.00000e+00, 8.31090e-02,
        1.12818e-01,-5.75005e-02,-1.29919e-02,-1.78849e-02,-2.86343e-06,
        0.00000e+00,-1.51187e+02,-6.65902e-03, 0.00000e+00,-2.02069e-03,
        0.00000e+00, 0.00000e+00, 4.32264e-02,-2.80444e+01,-3.26789e-03,
        2.47461e-03, 0.00000e+00, 0.00000e+00, 9.82100e-02, 1.22714e-01,
       -3.96450e-02, 0.00000e+00,-2.76489e-03, 0.00000e+00, 1.87723e-03,
       -8.09813e-03, 4.34428e-05,-7.70932e-03, 0.00000e+00,-2.28894e-03,
       -5.69070e-03,-5.22193e-03, 6.00692e-03,-7.80434e+03,-3.48336e-03,
       -6.38362e-03,-1.82190e-03, 0.00000e+00,-7.58976e+01,-2.17875e-02,
       -1.72524e-02,-9.06287e-03, 0.00000e+00, 2.44725e-02, 8.66040e-02,
        1.05712e-01, 3.02543e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -6.01364e+03,-5.64668e-03,-2.54157e-03, 0.00000e+00, 3.15611e+02,
       -5.69158e-03, 0.00000e+00, 0.00000e+00,-4.47216e-03,-4.49523e-03,
        4.64428e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        4.51236e-02, 2.46520e-02, 6.17794e-03, 0.00000e+00, 0.00000e+00,
       -3.62944e-01,-4.80022e-02,-7.57230e+01,-1.99656e-03, 0.00000e+00,
       -5.18780e-03,-1.73990e-02,-9.03485e-03, 7.48465e-03, 1.53267e-02,
        1.06296e-02, 1.18655e-02, 2.55569e-03, 1.69020e-03, 3.51936e-02,
       -1.81242e-02, 0.00000e+00,-1.00529e-01,-5.10574e-03, 0.00000e+00,
        2.10228e-03, 0.00000e+00, 0.00000e+00,-1.73255e+02, 5.07833e-01,
       -2.41408e-01, 8.75414e-03, 2.77527e-03,-8.90353e-05,-5.25148e+00,
       -5.83899e-03,-2.09122e-02,-9.63530e-03, 9.77164e-03, 4.07051e-03,
        2.53555e-04,-5.52875e+00,-3.55993e-01,-2.49231e-03, 0.00000e+00,
        0.00000e+00, 2.86026e+01, 0.00000e+00, 3.42722e-04, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
      // ----   o density
        { 0,
        1.02315e+00,-1.59710e-01,-1.06630e-01,-1.77074e-02,-4.42726e-03,
        3.44803e-02, 4.45613e-02,-3.33751e-02,-5.73598e-02, 3.50360e-01,
        6.33053e-02, 2.16221e-02, 5.42577e-02,-5.74193e+00, 0.00000e+00,
        1.90891e-01,-1.39194e-02, 1.01102e+02, 8.16363e-02, 1.33717e-04,
        6.54403e-06, 3.10295e-03, 0.00000e+00, 0.00000e+00, 5.38205e-02,
        1.23910e-01,-1.39831e-02, 0.00000e+00, 0.00000e+00,-3.95915e-06,
        0.00000e+00,-7.14651e-01,-5.01027e-03, 0.00000e+00,-3.24756e-03,
        0.00000e+00, 0.00000e+00, 4.42173e-02,-1.31598e+01,-3.15626e-03,
        1.24574e-03,-1.47626e-03,-1.55461e-03, 6.40682e-02, 1.34898e-01,
       -2.42415e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.13666e-04,
       -5.40373e-03, 2.61635e-05,-3.33012e-03, 0.00000e+00,-3.08101e-03,
       -2.42679e-03,-3.36086e-03, 0.00000e+00,-1.18979e+03,-5.04738e-02,
       -2.61547e-03,-1.03132e-03, 1.91583e-04,-8.38132e+01,-1.40517e-02,
       -1.14167e-02,-4.08012e-03, 1.73522e-04,-1.39644e-02,-6.64128e-02,
       -6.85152e-02,-1.34414e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        6.07916e+02,-4.12220e-03,-2.20996e-03, 0.00000e+00, 1.70277e+03,
       -4.63015e-03, 0.00000e+00, 0.00000e+00,-2.25360e-03,-2.96204e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        3.92786e-02, 1.31186e-02,-1.78086e-03, 0.00000e+00, 0.00000e+00,
       -3.90083e-01,-2.84741e-02,-7.78400e+01,-1.02601e-03, 0.00000e+00,
       -7.26485e-04,-5.42181e-03,-5.59305e-03, 1.22825e-02, 1.23868e-02,
        6.68835e-03,-1.03303e-02,-9.51903e-03, 2.70021e-04,-2.57084e-02,
       -1.32430e-02, 0.00000e+00,-3.81000e-02,-3.16810e-03, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-9.05762e-04,-2.14590e-03,-1.17824e-03, 3.66732e+00,
       -3.79729e-04,-6.13966e-03,-5.09082e-03,-1.96332e-03,-3.08280e-03,
       -9.75222e-04, 4.03315e+00,-2.52710e-01, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----   n2 density
        { 0,
        1.16112e+00, 0.00000e+00, 0.00000e+00, 3.33725e-02, 0.00000e+00,
        3.48637e-02,-5.44368e-03, 0.00000e+00,-6.73940e-02, 1.74754e-01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
        1.26733e-01, 0.00000e+00, 1.03154e+02, 5.52075e-02, 0.00000e+00,
        0.00000e+00, 8.13525e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-2.50482e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.48894e-03,
        6.16053e-04,-5.79716e-04, 2.95482e-03, 8.47001e-02, 1.70147e-01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----   msis00r.gts3c.tlb
        { 0,
        9.44846e-01, 0.00000e+00, 0.00000e+00,-3.08617e-02, 0.00000e+00,
       -2.44019e-02, 6.48607e-03, 0.00000e+00, 3.08181e-02, 4.59392e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
        2.13260e-02, 0.00000e+00,-3.56958e+02, 0.00000e+00, 1.82278e-04,
        0.00000e+00, 3.07472e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 3.83054e-03, 0.00000e+00, 0.00000e+00,
       -1.93065e-03,-1.45090e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.23493e-03, 1.36736e-03, 8.47001e-02, 1.70147e-01,
        3.71469e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        5.10250e-03, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 3.68756e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----   o2 density
        { 0,
        1.35580e+00, 1.44816e-01, 0.00000e+00, 6.07767e-02, 0.00000e+00,
        2.94777e-02, 7.46900e-02, 0.00000e+00,-9.23822e-02, 8.57342e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 2.38636e+01, 0.00000e+00,
        7.71653e-02, 0.00000e+00, 8.18751e+01, 1.87736e-02, 0.00000e+00,
        0.00000e+00, 1.49667e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-3.67874e+02, 5.48158e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
        1.22631e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        8.17187e-03, 3.71617e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.10826e-03,
       -3.13640e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -7.35742e-02,-5.00266e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 1.94965e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----   ar density
        { 0,
        1.04761e+00, 2.00165e-01, 2.37697e-01, 3.68552e-02, 0.00000e+00,
        3.57202e-02,-2.14075e-01, 0.00000e+00,-1.08018e-01,-3.73981e-01,
        0.00000e+00, 3.10022e-02,-1.16305e-03,-2.07596e+01, 0.00000e+00,
        8.64502e-02, 0.00000e+00, 9.74908e+01, 5.16707e-02, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 3.46193e+02, 1.34297e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-3.48509e-03,
       -1.54689e-04, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
        1.47753e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        1.89320e-02, 3.68181e-05, 1.32570e-02, 0.00000e+00, 0.00000e+00,
        3.59719e-03, 7.44328e-03,-1.00023e-03,-6.50528e+03, 0.00000e+00,
        1.03485e-02,-1.00983e-03,-4.06916e-03,-6.60864e+01,-1.71533e-02,
        1.10605e-02, 1.20300e-02,-5.20034e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -2.62769e+03, 7.13755e-03, 4.17999e-03, 0.00000e+00, 1.25910e+04,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-2.23595e-03, 4.60217e-03,
        5.71794e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -3.18353e-02,-2.35526e-02,-1.36189e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.03522e-02,-6.67837e+01,-1.09724e-03, 0.00000e+00,
       -1.38821e-02, 1.60468e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.51574e-02,
       -5.44470e-04, 0.00000e+00, 7.28224e-02, 6.59413e-02, 0.00000e+00,
       -5.15692e-03, 0.00000e+00, 0.00000e+00,-3.70367e+03, 0.00000e+00,
        0.00000e+00, 1.36131e-02, 5.38153e-03, 0.00000e+00, 4.76285e+00,
       -1.75677e-02, 2.26301e-02, 0.00000e+00, 1.76631e-02, 4.77162e-03,
        0.00000e+00, 5.39354e+00, 0.00000e+00,-7.51710e-03, 0.00000e+00,
        0.00000e+00,-8.82736e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----    h density
        { 0,
        1.26376e+00,-2.14304e-01,-1.49984e-01, 2.30404e-01, 2.98237e-02,
        2.68673e-02, 2.96228e-01, 2.21900e-02,-2.07655e-02, 4.52506e-01,
        1.20105e-01, 3.24420e-02, 4.24816e-02,-9.14313e+00, 0.00000e+00,
        2.47178e-02,-2.88229e-02, 8.12805e+01, 5.10380e-02,-5.80611e-03,
        2.51236e-05,-1.24083e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01,-3.48190e-02, 0.00000e+00, 0.00000e+00, 2.89885e-05,
        0.00000e+00, 1.53595e+02,-1.68604e-02, 0.00000e+00, 1.01015e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.84552e-04,
       -1.22181e-03, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
       -1.04927e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,-5.91313e-03,
       -2.30501e-02, 3.14758e-05, 0.00000e+00, 0.00000e+00, 1.26956e-02,
        8.35489e-03, 3.10513e-04, 0.00000e+00, 3.42119e+03,-2.45017e-03,
       -4.27154e-04, 5.45152e-04, 1.89896e-03, 2.89121e+01,-6.49973e-03,
       -1.93855e-02,-1.48492e-02, 0.00000e+00,-5.10576e-02, 7.87306e-02,
        9.51981e-02,-1.49422e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        2.65503e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 6.37110e-03, 3.24789e-04,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        6.14274e-02, 1.00376e-02,-8.41083e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.27099e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -3.94077e-03,-1.28601e-02,-7.97616e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-6.71465e-03,-1.69799e-03, 1.93772e-03, 3.81140e+00,
       -7.79290e-03,-1.82589e-02,-1.25860e-02,-1.04311e-02,-3.02465e-03,
        2.43063e-03, 3.63237e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----    n density
        { 0,
        7.09557e+01,-3.26740e-01, 0.00000e+00,-5.16829e-01,-1.71664e-03,
        9.09310e-02,-6.71500e-01,-1.47771e-01,-9.27471e-02,-2.30862e-01,
       -1.56410e-01, 1.34455e-02,-1.19717e-01, 2.52151e+00, 0.00000e+00,
       -2.41582e-01, 5.92939e-02, 4.39756e+00, 9.15280e-02, 4.41292e-03,
        0.00000e+00, 8.66807e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 9.74701e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 6.70217e+01,-1.31660e-03, 0.00000e+00,-1.65317e-02,
        0.00000e+00, 0.00000e+00, 8.50247e-02, 2.77428e+01, 4.98658e-03,
        6.15115e-03, 9.50156e-03,-2.12723e-02, 8.47001e-02, 1.70147e-01,
       -2.38645e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.37380e-03,
       -8.41918e-03, 2.80145e-05, 7.12383e-03, 0.00000e+00,-1.66209e-02,
        1.03533e-04,-1.68898e-02, 0.00000e+00, 3.64526e+03, 0.00000e+00,
        6.54077e-03, 3.69130e-04, 9.94419e-04, 8.42803e+01,-1.16124e-02,
       -7.74414e-03,-1.68844e-03, 1.42809e-03,-1.92955e-03, 1.17225e-01,
       -2.41512e-02, 1.50521e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        1.60261e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-3.54403e-04,-1.87270e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        2.76439e-02, 6.43207e-03,-3.54300e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-2.80221e-02, 8.11228e+01,-6.75255e-04, 0.00000e+00,
       -1.05162e-02,-3.48292e-03,-6.97321e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.45546e-03,-1.31970e-02,-3.57751e-03,-1.09021e+00,
       -1.50181e-02,-7.12841e-03,-6.64590e-03,-3.52610e-03,-1.87773e-02,
       -2.22432e-03,-3.93895e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        // ----  hot o density
        { 0,
        6.04050e-02, 1.57034e+00, 2.99387e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-1.51018e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-8.61650e+00, 1.26454e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 5.50878e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 6.23881e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
       -9.45934e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        }
        };
        // ----    msis00r.gts3c.ps param
      static double ps[151] =
        { 0,
        9.56827e-01, 6.20637e-02, 3.18433e-02, 0.00000e+00, 0.00000e+00,
        3.94900e-02, 0.00000e+00, 0.00000e+00,-9.24882e-03,-7.94023e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.74677e-03, 0.00000e+00, 1.54951e-02, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-6.99007e-04, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 1.24362e-02,-5.28756e-03, 8.47001e-02, 1.70147e-01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        };
        // ----    turbo
      static double pdl[3][26] = {
        { 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0
        },
        { 0,
        1.09930e+00, 3.90631e+00, 3.07165e+00, 9.86161e-01, 1.63536e+01,
        4.63830e+00, 1.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 1.28840e+00, 3.10302e-02, 1.18339e-01
        },
        { 0,
        1.00000e+00, 7.00000e-01, 1.15020e+00, 3.44689e+00, 1.28840e+00,
        1.00000e+00, 1.08738e+00, 1.22947e+00, 1.10016e+00, 7.34129e-01,
        1.15241e+00, 2.22784e+00, 7.95046e-01, 4.01612e+00, 4.47749e+00,
        1.23435e+02,-7.60535e-02, 1.68986e-06, 7.44294e-01, 1.03604e+00,
        1.72783e+02, 1.15020e+00, 3.44689e+00,-7.46230e-01, 9.49154e-01
        }
        };
        // ----   lower bouneary
      static double ptm[11] =
        { 0,
        1.04130e+03, 3.86000e+02, 1.95000e+02, 1.66728e+01, 2.13000e+02,
        1.20000e+02, 2.40000e+02, 1.87000e+02,-2.00000e+00, 0.00000e+00
        };
      static double pdm[9][11] = {
        { 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        },
        { 0,
        2.45600e+07, 6.71072e-06, 1.00000e+02, 0.00000e+00, 1.10000e+02,
        1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        8.59400e+10, 1.00000e+00, 1.05000e+02,-8.00000e+00, 1.10000e+02,
        1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        2.81000e+11, 0.00000e+00, 1.05000e+02, 2.80000e+01, 2.89500e+01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        3.30000e+10, 2.68270e-01, 1.05000e+02, 1.00000e+00, 1.10000e+02,
        1.00000e+01, 1.10000e+02,-1.00000e+01, 0.00000e+00, 0.00000e+00
        },
        { 0,
        1.33000e+09, 1.19615e-02, 1.05000e+02, 0.00000e+00, 1.10000e+02,
        1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        1.76100e+05, 1.00000e+00, 9.50000e+01,-8.00000e+00, 1.10000e+02,
        1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        1.00000e+07, 1.00000e+00, 1.05000e+02,-8.00000e+00, 1.10000e+02,
        1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
        },
        { 0,
        1.00000e+06, 1.00000e+00, 1.05000e+02,-8.00000e+00, 5.50000e+02,
        7.60000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 4.00000e+03
        }
        };
        // ----   msis00r.meso.tn1[2]
      static double ptl[5][101] = {
        { 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        },
        { 0,
        1.00858e+00, 4.56011e-02,-2.22972e-02,-5.44388e-02, 5.23136e-04,
       -1.88849e-02, 5.23707e-02,-9.43646e-03, 6.31707e-03,-7.80460e-02,
       -4.88430e-02, 0.00000e+00, 0.00000e+00,-7.60250e+00, 0.00000e+00,
       -1.44635e-02,-1.76843e-02,-1.21517e+02, 2.85647e-02, 0.00000e+00,
        0.00000e+00, 6.31792e-04, 0.00000e+00, 5.77197e-03, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-8.90272e+03, 3.30611e-03, 3.02172e-03, 0.00000e+00,
       -2.13673e-03,-3.20910e-04, 0.00000e+00, 0.00000e+00, 2.76034e-03,
        2.82487e-03,-2.97592e-04,-4.21534e-03, 8.47001e-02, 1.70147e-01,
        8.96456e-03, 0.00000e+00,-1.08596e-02, 0.00000e+00, 0.00000e+00,
        5.57917e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 9.65405e-03, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----   msis00r.meso.tn1[3]
        { 0,
        9.39664e-01, 8.56514e-02,-6.79989e-03, 2.65929e-02,-4.74283e-03,
        1.21855e-02,-2.14905e-02, 6.49651e-03,-2.05477e-02,-4.24952e-02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 1.19148e+01, 0.00000e+00,
        1.18777e-02,-7.28230e-02,-8.15965e+01, 1.73887e-02, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-1.44691e-02, 2.80259e-04, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.16584e+02, 3.18713e-03, 7.37479e-03, 0.00000e+00,
       -2.55018e-03,-3.92806e-03, 0.00000e+00, 0.00000e+00,-2.89757e-03,
       -1.33549e-03, 1.02661e-03, 3.53775e-04, 8.47001e-02, 1.70147e-01,
       -9.17497e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        3.56082e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.00902e-02, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----   msis00r.meso.tn1[4]
        { 0,
        9.85982e-01,-4.55435e-02, 1.21106e-02, 2.04127e-02,-2.40836e-03,
        1.11383e-02,-4.51926e-02, 1.35074e-02,-6.54139e-03, 1.15275e-01,
        1.28247e-01, 0.00000e+00, 0.00000e+00,-5.30705e+00, 0.00000e+00,
       -3.79332e-02,-6.24741e-02, 7.71062e-01, 2.96315e-02, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 6.81051e-03,-4.34767e-03, 8.66784e-02,
        1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 1.07003e+01,-2.76907e-03, 4.32474e-04, 0.00000e+00,
        1.31497e-03,-6.47517e-04, 0.00000e+00,-2.20621e+01,-1.10804e-03,
       -8.09338e-04, 4.18184e-04, 4.29650e-03, 8.47001e-02, 1.70147e-01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -4.04337e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-9.52550e-04,
        8.56253e-04, 4.33114e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.21223e-03,
        2.38694e-04, 9.15245e-04, 1.28385e-03, 8.67668e-04,-5.61425e-06,
        1.04445e+00, 3.41112e+01, 0.00000e+00,-8.40704e-01,-2.39639e+02,
        7.06668e-01,-2.05873e+01,-3.63696e-01, 2.39245e+01, 0.00000e+00,
       -1.06657e-03,-7.67292e-04, 1.54534e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----   msis00r.meso.tn1[5] msis00r.meso.tn2[1]
        { 0,
        1.00320e+00, 3.83501e-02,-2.38983e-03, 2.83950e-03, 4.20956e-03,
        5.86619e-04, 2.19054e-02,-1.00946e-02,-3.50259e-03, 4.17392e-02,
       -8.44404e-03, 0.00000e+00, 0.00000e+00, 4.96949e+00, 0.00000e+00,
       -7.06478e-03,-1.46494e-02, 3.13258e+01,-1.86493e-03, 0.00000e+00,
       -1.67499e-02, 0.00000e+00, 0.00000e+00, 5.12686e-04, 8.66784e-02,
        1.58727e-01,-4.64167e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        4.37353e-03,-1.99069e+02, 0.00000e+00,-5.34884e-03, 0.00000e+00,
        1.62458e-03, 2.93016e-03, 2.67926e-03, 5.90449e+02, 0.00000e+00,
        0.00000e+00,-1.17266e-03,-3.58890e-04, 8.47001e-02, 1.70147e-01,
        0.00000e+00, 0.00000e+00, 1.38673e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.60571e-03,
        6.28078e-04, 5.05469e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-1.57829e-03,
       -4.00855e-04, 5.04077e-05,-1.39001e-03,-2.33406e-03,-4.81197e-04,
        1.46758e+00, 6.20332e+00, 0.00000e+00, 3.66476e-01,-6.19760e+01,
        3.09198e-01,-1.98999e+01, 0.00000e+00,-3.29933e+02, 0.00000e+00,
       -1.10080e-03,-9.39310e-05, 1.39638e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        }
        };
        // ----    msis00r.meso.tn2[2]
      static double pma[11][101] = {
        { 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        },
        { 0,
        9.81637e-01,-1.41317e-03, 3.87323e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-3.58707e-02,
       -8.63658e-03, 0.00000e+00, 0.00000e+00,-2.02226e+00, 0.00000e+00,
       -8.69424e-03,-1.91397e-02, 8.76779e+01, 4.52188e-03, 0.00000e+00,
        2.23760e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-7.07572e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
       -4.11210e-03, 3.50060e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-8.36657e-03, 1.61347e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-1.45130e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.24152e-03,
        6.43365e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.33255e-03,
        2.42657e-03, 1.60666e-03,-1.85728e-03,-1.46874e-03,-4.79163e-06,
        1.22464e+00, 3.53510e+01, 0.00000e+00, 4.49223e-01,-4.77466e+01,
        4.70681e-01, 8.41861e+00,-2.88198e-01, 1.67854e+02, 0.00000e+00,
        7.11493e-04, 6.05601e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    msis00r.meso.tn2[3]
        { 0,
        1.00422e+00,-7.11212e-03, 5.24480e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-5.28914e-02,
       -2.41301e-02, 0.00000e+00, 0.00000e+00,-2.12219e+01,-1.03830e-02,
       -3.28077e-03, 1.65727e-02, 1.68564e+00,-6.68154e-03, 0.00000e+00,
        1.45155e-02, 0.00000e+00, 8.42365e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-4.34645e-03, 0.00000e+00, 0.00000e+00, 2.16780e-02,
        0.00000e+00,-1.38459e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 7.04573e-03,-4.73204e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 1.08767e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-8.08279e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.21769e-04,
       -2.27387e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.26769e-03,
        3.16901e-03, 4.60316e-04,-1.01431e-04, 1.02131e-03, 9.96601e-04,
        1.25707e+00, 2.50114e+01, 0.00000e+00, 4.24472e-01,-2.77655e+01,
        3.44625e-01, 2.75412e+01, 0.00000e+00, 7.94251e+02, 0.00000e+00,
        2.45835e-03, 1.38871e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    msis00r.meso.tn2[4] tn3[1]
        { 0,
        1.01890e+00,-2.46603e-02, 1.00078e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-6.70977e-02,
       -4.02286e-02, 0.00000e+00, 0.00000e+00,-2.29466e+01,-7.47019e-03,
        2.26580e-03, 2.63931e-02, 3.72625e+01,-6.39041e-03, 0.00000e+00,
        9.58383e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.85291e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 1.39717e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 9.19771e-03,-3.69121e+02, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-1.57067e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-7.07265e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.92953e-03,
       -2.77739e-03,-4.40092e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.47280e-03,
        2.95035e-04,-1.81246e-03, 2.81945e-03, 4.27296e-03, 9.78863e-04,
        1.40545e+00,-6.19173e+00, 0.00000e+00, 0.00000e+00,-7.93632e+01,
        4.44643e-01,-4.03085e+02, 0.00000e+00, 1.15603e+01, 0.00000e+00,
        2.25068e-03, 8.48557e-04,-2.98493e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tn3[2]
        { 0,
        9.75801e-01, 3.80680e-02,-3.05198e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.85575e-02,
        5.04057e-02, 0.00000e+00, 0.00000e+00,-1.76046e+02, 1.44594e-02,
       -1.48297e-03,-3.68560e-03, 3.02185e+01,-3.23338e-03, 0.00000e+00,
        1.53569e-02, 0.00000e+00,-1.15558e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 4.89620e-03, 0.00000e+00, 0.00000e+00,-1.00616e-02,
       -8.21324e-03,-1.57757e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 6.63564e-03, 4.58410e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-2.51280e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 9.91215e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-8.73148e-04,
       -1.29648e-03,-7.32026e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-4.68110e-03,
       -4.66003e-03,-1.31567e-03,-7.39390e-04, 6.32499e-04,-4.65588e-04,
       -1.29785e+00,-1.57139e+02, 0.00000e+00, 2.58350e-01,-3.69453e+01,
        4.10672e-01, 9.78196e+00,-1.52064e-01,-3.85084e+03, 0.00000e+00,
       -8.52706e-04,-1.40945e-03,-7.26786e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tn3[3]
        { 0,
        9.60722e-01, 7.03757e-02,-3.00266e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.22671e-02,
        4.10423e-02, 0.00000e+00, 0.00000e+00,-1.63070e+02, 1.06073e-02,
        5.40747e-04, 7.79481e-03, 1.44908e+02, 1.51484e-04, 0.00000e+00,
        1.97547e-02, 0.00000e+00,-1.41844e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 5.77884e-03, 0.00000e+00, 0.00000e+00, 9.74319e-03,
        0.00000e+00,-2.88015e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-4.44902e-03,-2.92760e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 2.34419e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.36685e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-4.65325e-04,
       -5.50628e-04, 3.31465e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.06179e-03,
       -3.08575e-03,-7.93589e-04,-1.08629e-04, 5.95511e-04,-9.05050e-04,
        1.18997e+00, 4.15924e+01, 0.00000e+00,-4.72064e-01,-9.47150e+02,
        3.98723e-01, 1.98304e+01, 0.00000e+00, 3.73219e+03, 0.00000e+00,
       -1.50040e-03,-1.14933e-03,-1.56769e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tn3[4]
        { 0,
        1.03123e+00,-7.05124e-02, 8.71615e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-3.82621e-02,
       -9.80975e-03, 0.00000e+00, 0.00000e+00, 2.89286e+01, 9.57341e-03,
        0.00000e+00, 0.00000e+00, 8.66153e+01, 7.91938e-04, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 4.68917e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 7.86638e-03, 0.00000e+00, 0.00000e+00, 9.90827e-03,
        0.00000e+00, 6.55573e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-4.00200e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 7.07457e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.72268e-03,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.04970e-04,
        1.21560e-03,-8.05579e-06, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.49941e-03,
       -4.57256e-04,-1.59311e-04, 2.96481e-04,-1.77318e-03,-6.37918e-04,
        1.02395e+00, 1.28172e+01, 0.00000e+00, 1.49903e-01,-2.63818e+01,
        0.00000e+00, 4.70628e+01,-2.22139e-01, 4.82292e-02, 0.00000e+00,
       -8.67075e-04,-5.86479e-04, 5.32462e-04, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tn3[5] surface temp tsl
        { 0,
        1.00828e+00,-9.10404e-02,-2.26549e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-2.32420e-02,
       -9.08925e-03, 0.00000e+00, 0.00000e+00, 3.36105e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-1.24957e+01,-5.87939e-03, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.79765e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 2.01237e+03, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-1.75553e-02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.29699e-03,
        1.26659e-03, 2.68402e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.17894e-03,
        1.48746e-03, 1.06478e-04, 1.34743e-04,-2.20939e-03,-6.23523e-04,
        6.36539e-01, 1.13621e+01, 0.00000e+00,-3.93777e-01, 2.38687e+03,
        0.00000e+00, 6.61865e+02,-1.21434e-01, 9.27608e+00, 0.00000e+00,
        1.68478e-04, 1.24892e-03, 1.71345e-03, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tgn3[2] surface grae tslg
        { 0,
        1.57293e+00,-6.78400e-01, 6.47500e-01, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-7.62974e-02,
       -3.60423e-01, 0.00000e+00, 0.00000e+00, 1.28358e+02, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 4.68038e+01, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.67898e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.90994e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 3.15706e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tgn2[1] tgn1[2]
        { 0,
        8.60028e-01, 3.77052e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-1.17570e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 7.77757e-03, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 1.01024e+02, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 6.54251e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,-1.56959e-02,
        1.91001e-02, 3.15971e-02, 1.00982e-02,-6.71565e-03, 2.57693e-03,
        1.38692e+00, 2.82132e-01, 0.00000e+00, 0.00000e+00, 3.81511e+02,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        },
        // ----    tgn3[1] tgn2[2]
        { 0,
        1.06029e+00,-5.25231e-02, 3.73034e-01, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.31072e-02,
       -3.88409e-01, 0.00000e+00, 0.00000e+00,-1.65295e+02,-2.13801e-01,
       -4.38916e-02,-3.22716e-01,-8.82393e+01, 1.18458e-01, 0.00000e+00,
       -4.35863e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00,-1.19782e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 2.62229e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00,-5.37443e+01, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00,-4.55788e-01, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.84009e-02,
        3.96733e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.05494e-02,
        7.39617e-02, 1.92200e-02,-8.46151e-03,-1.34244e-02, 1.96338e-02,
        1.50421e+00, 1.88368e+01, 0.00000e+00, 0.00000e+00,-5.13114e+01,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        5.11923e-02, 3.61225e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
        }
        };
        // ----    semiannual mult sam
      static double sam[101] =
        { 0,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
        };
        // ----   middle atmosphere averages
      static double pavgm[11] =
        { 0,
        2.61000e+02, 2.64000e+02, 2.29000e+02, 2.17000e+02, 2.17000e+02,
        2.23000e+02, 2.86760e+02,-2.93940e+00, 2.50000e+00, 0.00000e+00
        };

      // now assign the record values
      memcpy(msis00r.parm.pt, pt, 151*sizeof(double));
      memcpy(msis00r.parm.ps, ps, 151*sizeof(double));
      memcpy(msis00r.lower.ptm, ptm, 11*sizeof(double));
      memcpy(msis00r.parm.sam, sam, 101*sizeof(double));
      memcpy(msis00r.mavg.pavgm, pavgm, 11*sizeof(double));

      // since c reads in rows first, while the fortran code assumed columns,
      // we do a loop

      for (i=0; i<= 8; i++)
          for (j=0; j<= 10; j++)
              msis00r.lower.pdm[j][i] = pdm[i][j];
      for (i=0; i<= 4; i++)
          for (j=0; j<= 100; j++)
              msis00r.parm.ptl[j][i] = ptl[i][j];
      for (i=0; i<= 9; i++)
          for (j=0; j<= 150; j++)
              msis00r.parm.pd[j][i] = pd[i][j];

      for (i=0; i<= 2; i++)
          for (j=0; j<= 25; j++)
              msis00r.parm.pdl[j][i] = pdl[i][j];
      for (i=0; i<= 10; i++)
          for (j=0; j<= 100; j++)
              msis00r.parm.pma[j][i] = pma[i][j];

// cdav
// these are the extra arrays needed so a whole array can be accessed, but in
// reverse order from the standard [row,col] approach in the fortran.
//    extras to aid function calls for a subarry
      memcpy(msis00r.parm.pddav, pd, 10*151*sizeof(double));
      memcpy(msis00r.parm.ptldav, ptl, 5*101*sizeof(double));
      memcpy(msis00r.parm.pmadav, pma, 11*101*sizeof(double));

     }

