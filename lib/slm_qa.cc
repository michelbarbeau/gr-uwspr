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

// Tester of Straight Line Model
// Build with:
//    make -f Makefile_slm_qa

#include "slm.h"

using namespace gr::uwspr;
#include <iostream>
using namespace std;

int main()
{
    // create a Straight Line Model instance
    SLM * aSLM = new SLM();
    // carrier frequency (Hz)
    int cf = 1500;
    // time (s)
    float t = 0.0;
    // mobility parameters
    mode_nonlinear m_nl;
    // velocity vector
    m_nl.V1 = 1; // horizontal velocity (m/s)
    m_nl.V2 = -2.0; // vertical velocity (m/s)
    // start coordinate 
    m_nl.p1 = 0;// x-axis, always null
    m_nl.p2 = 50;// y-axis (m)
    // calculate frequency drift
    float drift;
    for (int i = 0; i<120; i++)
    {
       drift = aSLM->slmFrequencyDrift(m_nl, (float)cf, (float)i);
       cout << drift << " ";
    }
    cout << endl << "End of test" << endl;
    return 0;
}
