/* ---------------------------------------------------------------------
*
*                              testcpp.cpp
*
*  this program tests the sgp4 propagator. an stk ephemeris file is generated
*  along with the test output. the code for this is left justified for easy
*  location.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2007
*                            by david vallado
*
*     (w) 719-573-2600, email dvallado@agi.com
*     *****************************************************************
*  current :
*             3 sep 08  david vallado
*                        add switch for afspc compatibility and improved operation
*  changes :
*            14 may 08  david vallado
*                        fixes for linux suggested by brian micek
*                        misc fixes noted by the community - manual operation,
*                        formats, char lengths
*            14 aug 06  david vallado
*                        update mfe for verification time steps, constants
*            20 jul 05  david vallado
*                         fixes for paper, corrections from paul crawford
*             7 jul 04  david vallado
*                         fix record file and get working
*            14 may 01  david vallado
*                         2nd edition baseline
*                   80  norad
*                         original baseline
*       ----------------------------------------------------------------      */

#include <cstdio>
#include <cstring>
#include <cmath>

#include <fstream>

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "sgp4io.h"

namespace sgp4 {

	int main()
	{
		char str[2];
		char infilename[15];
		double ro[3];
		double vo[3];
		char typerun, typeinput, opsmode;
		gravconsttype  whichconst;
		int whichcon;
		FILE* infile;
      FILE* outfile;
      FILE* outfilee;

		// ----------------------------  locals  -------------------------------
		double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
		double sec, jd, rad, tsince, startmfe, stopmfe, deltamin;
		double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
		int  year; int mon; int day; int hr; int min;
		char longstr1[130];
		typedef char str3[4];
		str3 monstr[13];
		char outname[64];
		char longstr2[130];
		elsetrec satrec;

		rad = 180.0 / pi;
		// ------------------------  implementation   --------------------------
		std::strcpy(monstr[1], "Jan");
		std::strcpy(monstr[2], "Feb");
		std::strcpy(monstr[3], "Mar");
		std::strcpy(monstr[4], "Apr");
		std::strcpy(monstr[5], "May");
		std::strcpy(monstr[6], "Jun");
		std::strcpy(monstr[7], "Jul");
		std::strcpy(monstr[8], "Aug");
		std::strcpy(monstr[9], "Sep");
		std::strcpy(monstr[10], "Oct");
		std::strcpy(monstr[11], "Nov");
		std::strcpy(monstr[12], "Dec");

		std::printf("%s\n", SGP4Version);

		//opsmode = 'a' best understanding of how afspc code works
		//opsmode = 'i' imporved sgp4 resulting in smoother behavior
      std::printf("input operation mode a, i \n\n");
		opsmode = getchar();
      std::fflush(stdin);

		//typerun = 'c' compare 1 year of full satcat data
		//typerun = 'v' verification run, requires modified elm file with
		//              start, stop, and delta times
		//typerun = 'm' maunual operation- either mfe, epoch, or dayof yr also
      std::printf("input type of run c, v, m \n\n");
		typerun = getchar();
      std::fflush(stdin);

		//typeinput = 'm' input start stop mfe
		//typeinput = 'e' input start stop ymd hms
		//typeinput = 'd' input start stop yr dayofyr
		if ((typerun != 'v') && (typerun != 'c')) {
			printf("input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d \n\n");
			typeinput = getchar();
      } else {
			typeinput = 'e';
      }
      std::printf("input which constants 721 72 84 \n");
      std::scanf("%i", &whichcon);
		if (whichcon == 721) whichconst = wgs72old;
		if (whichcon == 72) whichconst = wgs72;
		if (whichcon == 84) whichconst = wgs84;

		getgravconst(whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2);

		// ---------------- setup files for operation ------------------
		// input 2-line element set file
      std::printf("input elset filename: \n");
      std::scanf("%s", infilename);
      infile = std::fopen(infilename, "r");
		if (infile == NULL) {
         std::printf("Failed to open file: %s\n", infilename);
			return 1;
		}

		if (typerun == 'c')
         outfile = std::fopen("tcppall.out", "w");
		else {
			if (typerun == 'v')
            outfile = std::fopen("tcppver.out", "w");
			else
            outfile = std::fopen("tcpp.out", "w");
		}

		//        dbgfile = fopen("sgp4test.dbg", "w");
		//        fprintf(dbgfile,"this is the debug output\n\n" );

		// ----------------- test simple propagation -------------------
      while (std::feof(infile) == 0)
		{
			do
			{
            std::fgets(longstr1, 130, infile);
            std::strncpy(str, &longstr1[0], 1);
				str[1] = '\0';
         } while ((std::strcmp(str, "#") == 0) && (std::feof(infile) == 0));

         if (std::feof(infile) == 0) {
            std::fgets(longstr2, 130, infile);
				// convert the char string to sgp4 elements
				// includes initialization of sgp4
				twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
					startmfe, stopmfe, deltamin, satrec);
            std::fprintf(outfile, "%ld xx\n", satrec.satnum);
            std::printf(" %ld\n", satrec.satnum);
				// call the propagator to get the initial state vector value
				sgp4(whichconst, satrec, 0.0, ro, vo);

				// generate .e files for stk
				jd = satrec.jdsatepoch;
            std::strncpy(outname, &longstr1[2], 5);
				outname[5] = '.';
				outname[6] = 'e';
				outname[7] = '\0';
				invjday(jd, year, mon, day, hr, min, sec);
            outfilee = std::fopen(outname, "w");
            std::fprintf(outfilee, "stk.v.4.3 \n"); // must use 4.3...
            std::fprintf(outfilee, "\n");
            std::fprintf(outfilee, "BEGIN Ephemeris \n");
            std::fprintf(outfilee, " \n");
            std::fprintf(outfilee, "NumberOfEphemerisPoints		146 \n");
            std::fprintf(outfilee, "ScenarioEpoch	  %3i %3s%5i%3i:%2i:%12.9f \n", day, monstr[mon], year, hr, min, sec);
            std::fprintf(outfilee, "InterpolationMethod		Lagrange \n");
            std::fprintf(outfilee, "InterpolationOrder		5 \n");
            std::fprintf(outfilee, "CentralBody				Earth \n");
            std::fprintf(outfilee, "CoordinateSystem			TEME \n");
            std::fprintf(outfilee, "CoordinateSystemEpoch	%3i %3s%5i%3i:%2i:%12.9f \n", day, monstr[mon], year, hr, min, sec);
            std::fprintf(outfilee, "DistanceUnit			Kilometers \n");
            std::fprintf(outfilee, " \n");
            std::fprintf(outfilee, "EphemerisTimePosVel \n");
            std::fprintf(outfilee, " \n");
            std::fprintf(outfilee, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n", satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

            std::fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n", satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

				tsince = startmfe;
				// check so the first value isn't written twice
            if (std::fabs(tsince) > 1.0e-8)
					tsince = tsince - deltamin;

				// ----------------- loop to perform the propagation ----------------
				while ((tsince < stopmfe) && (satrec.error == 0))
				{
					tsince = tsince + deltamin;

					if (tsince > stopmfe)
						tsince = stopmfe;

					sgp4(whichconst, satrec, tsince, ro, vo);

					if (satrec.error > 0)
                  std::printf("# *** error: t:= %f *** code = %3d\n",
						satrec.t, satrec.error);

					if (satrec.error == 0)
					{
						if ((typerun != 'v') && (typerun != 'c'))
						{
							jd = satrec.jdsatepoch + tsince / 1440.0;
							invjday(jd, year, mon, day, hr, min, sec);

                     std::fprintf(outfile,
								" %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %5i%3i%3i %2i:%2i:%9.6f\n",
								tsince, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2], year, mon, day, hr, min, sec);
							//                            std::fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
							//                                           tsince,ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);
						}
						else
						{
							jd = satrec.jdsatepoch + tsince / 1440.0;
							invjday(jd, year, mon, day, hr, min, sec);

                     std::fprintf(outfilee, " %16.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
								tsince*60.0, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

                     std::fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f",
								tsince, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

							rv2coe(ro, vo, mu, p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper);
                     std::fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f\n",
								a, ecc, incl*rad, node*rad, argp*rad, nu*rad,
								m*rad, year, mon, day, hr, min, sec);
						}
					} // if satrec.error == 0

				} // while propagating the orbit

            std::fprintf(outfilee, " END Ephemeris \n");
            std::fclose(outfilee);

			} // if not eof

		} // while through the input file


		return 0;
	}  // end testcpp

}

int main(int argc, char**)
{
	return sgp4::main();
}

