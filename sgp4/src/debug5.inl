     if (help == 'y')
       //    if (dbgfile != NULL)
       {
        std::printf( "%85s\n",
                         " ------------------after initl  :---------------");
        std::printf( "    inputs : \n");
        std::printf( "%7s%15d%7s%15s%7s%15.9f%7s%15.9f%7s%15.9f\n",
         "satn", satn, "yr", " ", "ecco", ecco, "epoch", epoch, "inclo", inclo);
        std::printf( "    in/out : \n");
        std::printf( "%7s%15.9f\n", "no", no);
        std::printf( "    outputs : \n");
        std::printf(
              "%7s%15c%7s%15.9f%7s%15.9f%7s%15.9f%7s%15.9f%7s%15.9f\n",
               "method", method, "ainv", ainv, "ao", ao, "con41", con41,
               "con42", con42, "cosio", cosio);
        std::printf( "%7s%15.9f\n", "cosio2", cosio2);
        std::printf( "%7s%15.9f%7s%15.9f%7s%15.9f%7s%15.9f%7s%15.9f%7s%15.9f\n",
                         "einx", eccsq, "eccsq", eccsq, "omeosq", omeosq,
                         "posq", posq, "rp", rp, "rteosq", rteosq);
        std::printf( "%7s%15.9f%7s%15.9f\n", "sinio", sinio, "gsto", gsto);
       }

