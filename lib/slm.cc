/* -*- c++ -*- */
/*
 * Copyright 2018 Michel Barbeau, Carleton University.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

 #include "slm.h"
 #include <math.h>

 #ifndef DEBUG
 #define DEBUG 1 // set debug mode
 #endif
 #include "debugmacro.h"

 namespace gr {
  namespace uwspr {
    // Returns frequency drift (in Hz) due to the Doppler effect according to
    // the straight line model, described in paper:
    // Abdel-Mehsen Ahmad, Michel Barbeau, Joaquin Garcia-Alfaro, Jamil Kassem,
    // Evangelos Kranakis and Steven Porretta, "Low Frequency Mobile
    // Communications in Underwater Networks", 2018.
    // The body of this function was originally coded by:
    // Jamil Kassem, Lebanese International University, 2018
    float SLM::slmFrequencyDrift(mode_nonlinear m_nl, float cf, float t)
     // m_nl = trajectory parameters
     // t = time, in seconds
    {
      const float c = 1500.0; // speed of sound in water, in m/s
      // velocity vetor of vehicle of "a" (Eq. 11)
      double a[] = {m_nl.va / sqrt(1 + m_nl.ma * m_nl.ma),
        m_nl.ma * m_nl.va / sqrt(1 + m_nl.ma * m_nl.ma)};
      // initial position of vehicle "a"
      double aBar[] = {m_nl.xa, m_nl.ya};
      // velocity vetor of vehicle of "b" (Eq. 12)
      double b[] = {m_nl.vb / sqrt(1 + m_nl.mb *m_nl.mb),
        m_nl.mb * m_nl.vb / sqrt(1 + m_nl.mb * m_nl.mb)};
      // initial position vehicle "b"
      double bBar[] = { m_nl.xb, m_nl.yb };
      // sign of difference velocity vector projections on position
      // connecting vector
      float Sign=-((b[0]-a[0]+b[1]-a[1])>0)*2-1;
      //
      double normOfDotProduct=(a[0]-b[0])*((b[0]-a[0])*t+bBar[0]-aBar[0]) +
        (a[1]-b[1])*((b[1]-a[1])*t+bBar[1]-aBar[1]);
      // norm of connectiong vector (Eq. 16)
      double normOfUt=sqrt(pow((b[0]-a[0])*t+bBar[0]-aBar[0],2) +
        pow((b[1]-a[1])*t+bBar[1]-aBar[1],2));
      if(normOfUt==0) {
        return 0.0;
      } else {
        return cf*(Sign*normOfDotProduct)/(normOfUt*c);
      }
    }

    // Generators of trajectory parameters for the straight line model
    // The body of this function was originally coded by:
    // Jamil Kassem, Lebanese International University, 2018
    bool SLM::slmGenerator(mode_nonlinear *m_nl)
    {
      // parameters of vehicle "b" are constant
      m_nl->vb = 0;
      m_nl->mb = 0;
      m_nl->xb = 0;
      m_nl->yb = 0;
      // parameters of vehicle "a" are variable
      const double mamin=-1, mamax=1, mastep=1; // slope variables
      const double vamin=0, vamax=10, vastep=1; // speed variables
      const int xmin=-1000, xmax=1000, xstep=100; //x variables
      const int ymin=0, ymax=0, ystep=200; // y variables
      const int nma=(mamax-mamin)/mastep+1;
      const int nva=(vamax-vamin)/vastep+1;
      const int nx=(xmax-xmin)/xstep+1;
      const int ny=(ymax-ymin)/ystep+1;
      const int laststep=ny*nx*nva*nma; // total number of iterations
      static int iy, ix, iva, ima; // variable indices
      if (step==0) {
        iy=0; ix=0; iva=0; ima=0;
      }
      if (step<laststep) {
       if (iy >= ny) {
         iy = 0;
         ix++;
         if (ix >= nx) {
           ix = 0;
           iva++;
             if (iva >= nva) {
               iva = 0;
               ima++;
             }
         }
       }
       m_nl->va = iva * vastep + vamin;
       m_nl->ma = ima * mastep + mamin;
       m_nl->ya = iy * ystep + ymin;
       m_nl->xa = ix * xstep + xmin;
       iy++;
       step++;
       return true;
      } else {
       return false; // reach the end
      }
    }

    void SLM::slmGeneratorInit()
    {
       step = 0;
    }

  } /* namespace uwspr */
} /* namespace gr */
