/* -*- c++ -*- */
/*
 * Copyright 2019 Michel Barbeau, Carleton University.
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
    // Evangelos Kranakis, "Tuning the Demodulation Frequency Based on a Normalized 
    // Trajectory Model for Mobile Unerwater Acoustic Communications", 2019.
    float SLM::slmFrequencyDrift(mode_nonlinear m_nl, float cf, float t)
     // m_nl = trajectory parameters
     // t = time, in seconds
    {
      const float c = 1500.0; // sound speed (m/s)
      // sign of velocity vector
      float Sign = 
         ( 
           (
              (m_nl.V1 * t + m_nl.p1) * m_nl.V1 
              + 
              (m_nl.V2 * t + m_nl.p2) * m_nl.V2 
           )
           > 
           0  
         ) 
         * 
         2 - 1;
      //
      double numerator = 
         abs(
            m_nl.V1 * ( m_nl.V1 * t+m_nl.p1 ) 
            +
            m_nl.V2 * ( m_nl.V2 * t+m_nl.p2 )
         );
      // norm of connectiong vector (Eq. 16)
      double denominator = 
         sqrt( pow(m_nl.V1 * t + m_nl.p1, 2) 
               +
               pow(m_nl.V2 * t + m_nl.p2, 2)
         );
      if(denominator==0) {
        return 0.0;
      } else {
        // return -Sign;
        return -Sign * numerator/denominator * cf/c;
      }
    }

    // Generators of trajectory parameters for the straight line model
    bool SLM::slmGenerator(mode_nonlinear *m_nl)
    {
      // control variables for init velocity
      const double V1_min = -2, V1_max = 2, V1_step = 1;
      const int nV1 = (V1_max - V1_min) / V1_step + 1;
      const double V2_min = -2, V2_max = 2, V2_step = 1;
      const int nV2 = (V2_max - V2_min) / V2_step + 1;
      // control variables for init position on y-axis
      const int p2_min = 50, p2_max = 850, p2_step = 200; 
      const int np2 = (p2_max - p2_min) / p2_step + 1;
      // number of generated instances
      const int last = nV1 * nV2 * np2; 
      // indices into the instances
      static int ip2, iV1, iV2; 
      if (current==0) {
        ip2 = 0; iV1 = 0; iV2 = 0; 
      }
      if (current < last) { // not generated all instances?
        if (ip2 >= np2) {
          ip2 = 0; // reset index into the list positions
          iV1++; // next horizontal velocity
          if (iV1 >= nV1) {
            iV1 = 0; // reset index into the list horiz. velocities
            iV2++;  // next vertical velocity
          }
        }
        // map horizontal velocity index to horizontal velocity (m/s)
        m_nl->V1 = iV1 * V1_step + V1_min;
        // map vertical velocity index to vertical velocity (m/s)
        m_nl->V2 = iV2 * V2_step + V2_min;
        // init coordinate on x-axis is always null
        m_nl->p1 = 0;
        // map y-axis coordinate index to y-axis coordinate (m)
        m_nl->p2 = ip2 * p2_step + p2_min;
        ip2++; // next position on y-axis
        current++; // index of next instance
        return true;
      } else {
        return false; // reach the end
      }
    }

    void SLM::slmGeneratorInit()
    {
       current = 0;
    }

  } /* namespace uwspr */
} /* namespace gr */
