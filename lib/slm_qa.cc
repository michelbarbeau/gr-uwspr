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
    // carrier frequency (in Hz)
    int cf = 1500;
    // time
    float t = 0.0;
    // mobility parameters
    mode_nonlinear m_nl;
    // vehicle "a"
    m_nl.va = 6.0; // velocity
    m_nl.ma = 1.0; // slope
    m_nl.xa = -200; m_nl.ya = 0; // start coordinates
    // vehicle "b"
    m_nl.vb = 0.0; // velocity
    m_nl.mb = 0.0; // slope
    m_nl.xb = 0; m_nl.yb = 0; // start coordinates
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
