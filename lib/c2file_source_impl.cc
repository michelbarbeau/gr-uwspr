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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "c2file_source_impl.h"
#include <boost/math/special_functions/round.hpp>
#include <volk/volk.h>

namespace gr {
namespace uwspr {

c2file_source::sptr
c2file_source::make(const char *filename, bool repeat, float drift_rate)
{
        return gnuradio::get_initial_sptr
                       (new c2file_source_impl(filename, repeat, drift_rate/375.0));
}

/*
 * The private constructor
 *
 * The code is based on work by Steven Franke, K9AN,
 * which in turn was based on earlier work by K1JT.
 *
 * Copyright 2014-2015, Steven Franke, K9AN
 */
c2file_source_impl::c2file_source_impl(const char* filename, bool repeat, float drift_rate)
        : gr::sync_block("c2file_source",
                         gr::io_signature::make(0, 0, 0),
                         gr::io_signature::make(1,1, sizeof(gr_complex))),
        d_repeat(repeat), d_drift_rate(drift_rate)
{
        char *ptr_to_infile; double *freq; int *wspr_type;

        float *buffer;
        double dfreq;
        int i,ntrmin;
        char *c2file[15];
        FILE* fp;

        const int alignment_multiple =
                volk_get_alignment() / sizeof(float);
        set_alignment(std::max(1,alignment_multiple));

        maxpts=65536;
        sample_idx = 0;
        idat=(float *)malloc(sizeof(float)*maxpts);
        qdat=(float *)malloc(sizeof(float)*maxpts);

        buffer=(float *)malloc(sizeof(float)*2*65536);
        memset(buffer,0,sizeof(float)*2*65536);

        //fp = fopen(ptr_to_infile,"rb");
        fp = fopen(filename,"rb");
        if (fp == NULL) {
                fprintf(stderr, "Cannot open data file '%s'\n", ptr_to_infile);
                throw std::runtime_error("can't open file");
        }
        unsigned long nread=fread(c2file,sizeof(char),14,fp);
        nread=fread(&ntrmin,sizeof(int),1,fp);
        nread=fread(&dfreq,sizeof(double),1,fp);
        //*freq=dfreq;
        nread=fread(buffer,sizeof(float),2*45000,fp);
        fclose(fp);

        //*wspr_type=ntrmin;

        for(i=0; i<45000; i++) {
                idat[i]=buffer[2*i];
                qdat[i]=-buffer[2*i+1];
        }

        if( nread != 2*45000 ) {
                throw std::runtime_error("invalid number of samples");
        }
        free(buffer);
}

/*
 * Our virtual destructor.
 */
c2file_source_impl::~c2file_source_impl()
{
        free(idat); free(qdat);
}

int
c2file_source_impl::work(int noutput_items,
                         gr_vector_const_void_star &input_items,
                         gr_vector_void_star &output_items)
{
        gr_complex *out = (gr_complex *) output_items[0];
        int i;
        int nout;
        nout = 0;
        static float drift = 0.0; // simulated frequency drift
        int maxsamples = 45000;
        if (sample_idx >=maxsamples) {
                if(!d_repeat) {
                        return -1; // done!
                } else {
                        sample_idx = 0;
                }
        }
        for(i = 0; i < noutput_items; i++) {
                if (sample_idx >=maxsamples) {
                        break;
                }
                out[i] = gr_complex(idat[sample_idx],qdat[sample_idx])*
                   exp(gr_complex(0, sample_idx * M_PI * drift / 375.0));
                drift += d_drift_rate; // increment frequency drift
                sample_idx++;
                nout++;
        }
        // Tell runtime system how many output items we produced.
        return nout;
}

}   /* namespace uwspr */
} /* namespace gr */
