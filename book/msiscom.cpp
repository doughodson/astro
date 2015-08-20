/*       ----------------------------------------------------------------
*
*                              msiscom.cpp
*
*  this file contains common routines between the msis models.
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
*               4 oct 04  david vallado
*                           misc updates
*               6 aug 04  david vallado
*                           convert to c++
*              14 feb 03  david vallado
*                           misc updates
*              28 may 02  david vallado
*                           fix densu per l schwartz correction - tests correct
*              15 apr 02  david vallado
*                           work on gsurf
*              24 feb 02  david vallado
*                           finish double conversion
*               1 feb 02  nrl
*                           original baseline
*     http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
*     *****************************************************************       */

      #include "MSIS_MSIScom.h"

/*    *****************************************************************
*        set up routines that are common to all the msis models
*     *****************************************************************       */

//----------------------------------------------------------------------
//     set switches
//     output in  common/csw/sw(25),isw,swc(25)
//     sw for main terms, swc for cross terms
//
//     to turn on and off particular variations call tselec(sv),
//     where sv is a 25 element array containing 0.0 for off, 1.0
//     for on, or 2.0 for main effects off but cross terms on
//
//  cdav this was deleted since one can simply examine the csw structure
//     to get current values of sw: call tretrv(sw)
//----------------------------------------------------------------------

      void tselec
           (
            cswtype& csw, int sv[26]
           )
        {
//        int sav[26];
        int i;

//------------------------------ begin ---------------------------------
        for (i = 1; i <= 25; i++)
          {
//            sav[i] = sv[i];
            csw.sw[i] = fmod(sv[i],2);
            if (fabs(sv[i]) == 1  |  fabs(sv[i]) == 2)
                csw.swc[i] = 1;
              else
                csw.swc[i] = 0;
          }

        csw.isw = 64999;
        }

//----------------------------------------------------------------------
//     turbopause correction for msis models
//     eq. a12b
//         root mean density
//          dd   - diffusive density
//          dm   - full mixed density
//          zhm  - transition scale length
//          xmm  - full mixed molecular weight
//          xm   - species molecular weight
//          dnet - combined density
//       8/20/80
//----------------------------------------------------------------------

      double dnet
           (
            double dd, double dm, double zhm, double xmm, double xm
           )
        {
        double a, ylog;

//------------------------------ begin ---------------------------------
        a = zhm/(xmm-xm);

        if (dm <= 0.0  |  dd <= 0.0 )
          {
// cdav this loop is exercised with test case # 11 and msis86. i changed the return
// from the first if, and then used else's for the remaining ones.
// essentially, it appears that msis86 does not properly work down to 0km alt.
// the test cases match for 100 km, but no others were given. I expanded the
// original test cases and this cpp version matches 30, 50, and 70 km. the msis86
//  model gets nan's for 0 and 10 km. the answers below 100 km are all suspect.
            printf("dnet log error %11.7g %11.7g %11.7f %11.7f \n",dm,dd,xm,a);
            printf("dd and dm may need to be output vars for this case\n");

            if (dd == 0.0  &&  dm == 0.0)
                return 1.0;    // dd = 1.0
              else
              {
                if (dm == 0.0)
                    return dd;
                  else
//                     if (dd == 0.0)
                         return dm;
              }
          }
          else
          {
            // ---- eq. a12a ----
            ylog = a*log(dm/dd);
            if (ylog < -10.0)
                return dd;
              else
              {
                if (ylog > 10.0)
                    return dm;
                  else
                    return dd*pow(1.0 + exp(ylog),(1.0/a));
              }
          }

        }

//----------------------------------------------------------------------
//     chemistry/dissociation correction for msis models
//        alt - altitude
//        r   - target ratio
//        h1  - transition scale length
//        zh  - altitude of 1/2 r
//     eq. a20a or eq. a21
//----------------------------------------------------------------------

      double ccor
           (
            double alt, double r, double h1, double zh
           )
        {
        double e, ex, ccor1;
//------------------------------ begin ---------------------------------
        e = (alt-zh)/h1;
        if (e > 70.0)
            ccor1 = 0.0;
          else
          {
            if (e < -70.0)
                ccor1 = r;
              else
                {
                ex = exp(e);
                ccor1 = r/(1.0 + ex);
                }
           }
        return exp(ccor1);
        }

      double ccor2
           (
            double alt, double r, double h1, double zh, double h2
           )
        {
        double e1, e2, ccor1;
//------------------------------ begin ---------------------------------
        e1 = (alt-zh)/h1;
        e2 = (alt-zh)/h2;
        if (e1 > 70.0 | e2 > 70.0)
            ccor1 = 0.0;
          else
          {
            if (e1 < -70.0 && e2 < -70.0)
                ccor1 = r;
              else
                ccor1 = r/( 1.0 + 0.5*(exp(e1) + exp(e2)) );
           }
        return exp(ccor1);
        }

//-----------------------------------------------------------------------
//      control option of converting to kg/m3 or lb/ft3.
//-----------------------------------------------------------------------

      void meters
           (
             metseltype& metsel
           )
        {
        if (metsel.imr == 0)
            metsel.imr = 1;
          else
            metsel.imr = 0;
        }

/*     *****************************************************************
*       set up routines that are common to the msis-90 and 00 models
*     **************************************************************** */

//----------------------------------------------------------------------
//       calculate temperature and density profiles for msis models
//       new lower thermo polynomial 10/30/89
//
//----------------------------------------------------------------------

// do as in anline function/////////
      double zeta
           (
             double zz, double zl, double re
           )
           {
           return (zz-zl)*(re+zl)/(re+zz);
           }

      double densu
           (
            parmbtype& parmb, lsqvtype& lsqv,
            double alt, double dlb, double tinf, double tlb, double xm, double alpha, double& tz,
            double zlb, double s2, int mn1, double zn1[6], double tn1[6], double tgn1[3]
           )
        {
        double zz, zl, z, za, zg2, dta, t1, t2, densu1,
               zgdif, zg, z1, yd1, yd2, y, y2out[6], glb, gamma,
               expl, densa, gamm, tt, ta, z2, x;
         static double  yi, xs[6], ys[6];
        int k;
        static int mn;

//        dimension zn1(mn1),tn1(mn1)

        const double rgas = 831.4;

//       statement function
//cdav        zeta(zz,zl) = (zz-zl)*(re+zl)/(re+zz);

//------------------------------ begin ---------------------------------
//       printf(6,*) 'db',alt,dlb,tinf,tlb,xm,alpha,zlb,s2,mn1,zn1,tn1
        densu1 = 1.0;

        // -------- joining altitude of bates and spline
        za = zn1[1];
        if (alt >= za)
            z = alt;
          else
            z = za;
//cdav        z = max(alt,za);

        // -------- geopotential altitude difference from zlb
        zg2 = zeta(z,zlb,parmb.re);

        // -------- bates temperature
        tt = tinf-(tinf-tlb)*exp(-s2*zg2);
        ta = tt;
        tz = tt;
        densu1 = tz;

      // --------------- calculate temperature below za ----------------
      // -------- temperature gradient at za from bates profile
      if (alt < za)
        {
          dta = (tinf-ta)*s2*(pow((parmb.re+zlb)/(parmb.re+za),2));
          tgn1[1] = dta;
          tn1[1] = ta;
          if (alt >= zn1[mn1])
              z = alt;
            else
              z = zn1[mn1];
//cdav          z  = max(alt,zn1[mn1]);

          mn = mn1;
          z1 = zn1[1];
          z2 = zn1[mn];
          t1 = tn1[1];
          t2 = tn1[mn];

          // -------- geopotental difference from z1
          zg = zeta(z,z1,parmb.re);
          zgdif = zeta(z2,z1,parmb.re);

          // -------- set up spline nodes
          for (k = 1; k <=mn; k++)
            {
              xs[k] = zeta(zn1[k],z1,parmb.re)/zgdif;
              ys[k] = 1.0/tn1[k];
            }

          // -------- end node derivatives
          yd1 = -tgn1[1]/(t1*t1)*zgdif;
          yd2 = -tgn1[2]/(t2*t2)*zgdif*pow( ((parmb.re+z2)/(parmb.re+z1)),2 );

          // -------- calculate spline coefficients
          spline(xs,ys,mn,yd1,yd2,y2out);
          x = zg/zgdif;
          splint(xs,ys,y2out,mn,x,y);

          // -------- temperature at altitude
          tz = 1.0/y;
          densu1 = tz;
//cdav
        }

          // ------------ calculate density above za -------------------
          if (xm != 0.0)
            {
              glb  = parmb.gsurf / pow( (1.0 + zlb/parmb.re),2 );
              gamma= xm*glb / (s2*rgas*tinf);
              expl = exp(-s2*gamma*zg2);
              if (expl > 50  |  tt <= 0.0)
                expl = 50.0;

              // -------- density at altitude
              densa = dlb*pow( (tlb/tt),(1.0+alpha+gamma))*expl;
              densu1 = densa;

              if (alt < za)
                {
                  // ---------- calculate density below za -------------
                  glb  = parmb.gsurf/pow( (1.0 + z1/parmb.re), 2 );
                  gamm = xm*glb*zgdif/rgas;

                  // -------- integrate spline temperatures
                  splini (xs,ys,y2out,mn,x,yi);
                  expl = gamm*yi;
                  if (expl > 50.0 | tz <= 0.0)
                      expl = 50.0;

                  // -------- density at altitude
                  densu1 = densu1*pow( (t1/tz),(1.0+alpha) ) *exp(-expl);
                }
            }

         return densu1;
        }

//----------------------------------------------------------------------
//
//
// calculate temperature and density profiles for lower atmos.
//
//----------------------------------------------------------------------

      double densm
            (
              parmbtype& parmb, fittype& fit, lsqvtype& lsqv,
              double alt, double d0, double xm, double& tz, int mn3, double zn3[6],
              double tn3[6], double tgn3[3], int mn2,
              double zn2[5], double tn2[5], double tgn2[3]
            )
        {
        double zz, zl, z, t1, t2, densm1,
               zgdif, zg, xs[11], z1, yd1, yd2, y, y2out[11], glb,
               ys[11], gamm, z2, x, yi, expl;
        int mn, k, kk;

//        dimension zn3(mn3),tn3(mn3)
//        dimension zn2(mn2),tn2(mn2)

        const double rgas = 831.4;

//       statement function
//cdav        zeta(zz,zl) = (zz-zl)*(re+zl)/(re+zz);

//------------------------------ begin ---------------------------------
      densm1 = d0;

      if (alt <= zn2[1])
        {
          // -------- stratosphere/mesosphere temperature
          if (alt >= zn2[mn2])
              z = alt;
            else
              z = zn2[mn2];
          mn = mn2;
          z1 = zn2[1];
          z2 = zn2[mn];
          t1 = tn2[1];
          t2 = tn2[mn];
          zg = zeta(z,z1,parmb.re);
          zgdif = zeta(z2,z1,parmb.re);

          // -------- set up spline nodes
          for (k = 1; k <= mn; k++)
            {
              xs[k] = zeta(zn2[k],z1,parmb.re)/zgdif;
              ys[k] = 1.0/tn2[k];
            }

          yd1 = -tgn2[1]/(t1*t1)*zgdif;
          yd2 = -tgn2[2]/(t2*t2)*zgdif*pow( (parmb.re+z2)/(parmb.re+z1),2 );

          // -------- calculate spline coefficients
          spline(xs,ys,mn,yd1,yd2,y2out);
          x = zg/zgdif;
          splint(xs,ys,y2out,mn,x,y);

          // -------- temperature at altitude
          tz = 1.0/y;

          if (xm != 0.0)
            {
              // ------ calculate stratosphere/mesosphere density ------
              glb = parmb.gsurf/pow( (1.0+z1/parmb.re),2 );
              gamm = xm*glb*zgdif/rgas;

              // -------- integrate temperature profile
              splini(xs,ys,y2out,mn,x,yi);
              expl = gamm*yi;
              if (expl > 50.0)
                  expl = 50.0;

              // -------- density at altitude
              densm1 = densm1*(t1/tz)*exp(-expl);
            }

          if (alt <= zn3[1])
            {
              // -------- troposphere/stratosphere temperature ---------
              z = alt;
              mn = mn3;
              z1 = zn3[1];
              z2 = zn3[mn];
              t1 = tn3[1];
              t2 = tn3[mn];
              zg = zeta(z,z1,parmb.re);
              zgdif = zeta(z2,z1,parmb.re);

              // -------- set up spline nodes
              for (k = 1; k <= mn; k++)
                {
                  xs[k] = zeta(zn3[k],z1,parmb.re)/zgdif;
                  ys[k] = 1.0/tn3[k];
                }
              yd1 = -tgn3[1]/(t1*t1)*zgdif;
              yd2 = -tgn3[2]/(t2*t2)*zgdif*pow( (parmb.re+z2)/(parmb.re+z1),2);

              // -------- calculate spline coefficients
              spline(xs,ys,mn,yd1,yd2,y2out);
              x = zg/zgdif;
              splint(xs,ys,y2out,mn,x,y);

              // -------- temperature at altitude
              tz = 1.0/y;
              if (xm != 0.0)
                {
                  // --- calculate tropospheric/stratosphere density ---
                  glb = parmb.gsurf/pow((1.0+z1/parmb.re),2);
                  gamm = xm*glb*zgdif/rgas;

                  // -------- integrate temperature profile
                  splini(xs,ys,y2out,mn,x,yi);
                  expl = gamm*yi;
                  if (expl > 50.0)
                      expl = 50.0;

                  // -------- density at altitude
                  densm1 = densm1*(t1/tz)*exp(-expl);
                }
            }
        }

      if (xm == 0.0)
          densm1 = tz;

      return densm1;
        }

//----------------------------------------------------------------------
//  calculate 2nd derivatives of cubic spline interp function
//  adapted from numerical recipes by press et al
//  x,y    : arrays of tabulated function in ascending order by x
//  n      : size of arrays x,y
//  yp1,ypn: specified derivatives at x[1] and x[n]; values
//           > =  1d30 signal signal second derivative zero
//  y2     : output array of second derivatives
//----------------------------------------------------------------------

      void spline
           (
            double x[6], double y[6], int n, double yp1, double ypn, double y2[6]
           )
        {
        int nmax, i, k;

        double sig, p, qn, un, u[101];
//        parameter (nmax = 101)
//        dimension x[n],y[n],y2[n],u(nmax)

//------------------------------ begin ---------------------------------
        if (yp1 > 0.99e30)
          {
            y2[1] = 0.0;
            u[1]  = 0.0;
          }
          else
          {
            y2[1] = -0.5;
            u[1]  = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
          }

        for (i = 2; i <= n-1; i++)
          {
            sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
            p   = sig*y2[i-1] + 2.0;
            y2[i] = (sig-1.0)/p;
            u[i]  = (6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])
                   /(x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1]) / p;
          }

        if (ypn > 0.99e30)
          {
            qn = 0.0;
            un = 0.0;
          }
          else
          {
            qn = 0.5;
            un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/ (x[n]-x[n-1]));
          }

        y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
        for (k = n-1; k >=1; k--)
            y2[k] = y2[k]*y2[k+1]+u[k];
        }

//----------------------------------------------------------------------
//  calculate cubic spline interp value
//  adapted from numberical recipes by press et al.
//  xa,ya: arrays of tabulated function in ascending order by x
//  y2a  : array of second derivatives
//  n    : size of arrays xa,ya,y2a
//  x    : abscissa for interpolation
//  y    : output value
//----------------------------------------------------------------------

      void splint
           (
            double xa[6], double ya[6], double y2a[6], int n, double x, double& y
           )
        {
        int k, khi, klo;
        double h, a, b;
//        dimension xa[n],ya[n],y2a[n]

//------------------------------ begin ---------------------------------
        klo = 1;
        khi = n;
        while (khi-klo > 1)
          {
            k = (khi+klo) * 0.5;
            if (xa[k] > x)
                khi = k;
              else
                klo = k;
          }

        h = xa[khi]-xa[klo];
        if (h == 0.0)
            printf("bad xa input to splint\n");
        a = (xa[khi]-x)/h;
        b = (x-xa[klo])/h;
        y = a*ya[klo]+b*ya[khi]+
            ((pow(a,3)-a)*y2a[klo] + (pow(b,3)-b)*y2a[khi])*h*h/6.0;

        }

//----------------------------------------------------------------------
//     integrate cubic spline function from xa[1] to x
//  xa,ya: arrays of tabulated function in ascending order by x
//  y2a  : array of second derivatives
//  n    : size of arrays xa,ya,y2a
//  x    : abscissa endpoint for integration
//  y    : output value
//----------------------------------------------------------------------

      void splini
           (
            double xa[6], double ya[6], double y2a[6], int n, double x, double& yi
           )
        {
        int khi, klo;
        double h, a, b, xx, a2, b2;

//------------------------------ begin ---------------------------------
        yi = 0;
        klo = 1;
        khi = 2;

        while (x > xa[klo] &&  khi <= n)
          {
            xx = x;
            if (khi < n)
              {
                  if (x <= xa[khi])
                      xx = x;
                    else
                      xx = xa[khi];
//                xx = min(x,xa[khi]);
              }

            h  = xa[khi]-xa[klo];
            a  = (xa[khi]-xx)/h;
            b  = (xx-xa[klo])/h;
            a2 = a*a;
            b2 = b*b;
            yi = yi + ((1.0 - a2)*ya[klo]*0.5 + b2*ya[khi]*0.5 +
                 ((-(1.0 + a2*a2)*0.25 + a2*0.5)*y2a[klo] +
                 (b2*b2/4.0 - b2*0.5)*y2a[khi])*h*h/6.0) * h;
            klo = klo + 1;
            khi = khi + 1;
          }

        }

//----------------------------------------------------------------------
//      calculate latitude variable gravity (gv) and effective
//      radius (reff)
//
//dav chg to *8
//----------------------------------------------------------------------

      void glatf
           (
            double lat, double& gv, double& reff
           )
       {
        double c2;

        const double  dgtr = 1.74533e-2;
        double  latl = -999.0;

//------------------------------ begin ---------------------------------
        if(lat != latl)
            c2 = cos(2.0*dgtr*lat);

        latl = lat;
        gv   = 980.616*(1.0-0.0026373*c2);
        reff = 2.0*gv/(3.085462e-6 + 2.27e-9*c2)*1.0e-5;

        }

//----------------------------------------------------------------------
//       test if geophysical variables or switches changed and save
//----------------------------------------------------------------------

      double vtst
           (
            cswtype& csw,
            int iyd, double sec, double glat, double glong, double stl,
            double f107a, double f107, double ap[8], int ic
           )
        {

        int i, iydl[3];
        double secl[3], glatl[3], gll[3], stll[3], fal[3], fl[3], apl[8][3],
               swl[26][3], swcl[26][3], vtst1;
//cdav
i=1;
        for (i=1; i<=3; i++)
          {
            iydl[i] = -999;
            secl[i] = -999.0;
            glatl[i] = -999.0;
            gll[i] = -999.0;
            stll[i] = -999.0;
            fal[i] = -999.0;
            fl[i] = -999.0;
          }
        for (i=1; i<=8; i++)
          {
            apl[i][1] = -999.0;
            apl[i][2] = -999.0;
            apl[i][3] = -999.0;
          }
        for (i=1; i<=26; i++)
          {
            swl[i][1] = -999.0;
            swl[i][2] = -999.0;
            swl[i][3] = -999.0;
            swcl[i][1] = -999.0;
            swcl[i][2] = -999.0;
            swcl[i][3] = -999.0;
          }

//------------------------------ begin ---------------------------------
        vtst1 = 0;
        if(  iyd != iydl[ic])  goto ten;
        if(  sec != secl[ic])  goto ten;
        if( glat != glatl[ic]) goto ten;
        if(glong != gll[ic])   goto ten;
        if(  stl != stll[ic])  goto ten;
        if(f107a != fal[ic])   goto ten;
        if( f107 != fl[ic])    goto ten;

        for (i = 1; i<=7; i++)
            if(ap[i] != apl[i][ic]) goto ten;

        for (i = 1; i<=25; i++)
          {
            if(csw.sw[i] != swl[i][ic]) goto ten;
            if(csw.swc[i] != swcl[i][ic]) goto ten;
          }

        goto twenty;

   ten:

        vtst1 = 1;

        iydl[ic]  = iyd;
        secl[ic]  = sec;
        glatl[ic] = glat;
        gll[ic]   = glong;
        stll[ic]  = stl;
        fal[ic]   = f107a;
        fl[ic]    = f107;

        for (i = 1; i<= 7; i++)
            apl[i][ic]=ap[i];

        for (i = 1; i<=25; i++)
          {
            swl[i][ic] = csw.sw[i];
            swcl[i][ic] = csw.swc[i];
          }

   twenty:

        return vtst1;

        }


