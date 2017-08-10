/* -*- c++ -*- */
/*
 * Copyright 2017 Michel Barbeau, Carleton University.
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
#include "WSPR_unpacker_impl.h"

namespace gr {
  namespace uwspr {

    WSPR_unpacker::sptr
    WSPR_unpacker::make()
    {
      return gnuradio::get_initial_sptr
        (new WSPR_unpacker_impl());
    }

    /*
     * The private constructor
     */
    WSPR_unpacker_impl::WSPR_unpacker_impl()
      : gr::block("WSPR_unpacker",
      gr::io_signature::make(0,0,0),
      gr::io_signature::make(0,0,0))
    {
      // interface initialization
      in_port = pmt::mp("in");
      out_port = pmt::mp("out");
      message_port_register_in(in_port);
      message_port_register_out(out_port);
      set_msg_handler(in_port, boost::bind(&WSPR_unpacker_impl::unpack, this ,_1));
      // log file initialization
      fall_wspr=fopen("ALL_WSPR.TXT","a");
      if (fall_wspr==NULL) {
        std::cerr << "In WSPR_unpacker_impl: Failed to open file ALL_WSPR.TXT\n";
        exit(EXIT_FAILURE);
      }
      fwsprd=fopen("wspr_spots.txt","w");
      if (fwsprd==NULL) {
        std::cerr << "In WSPR_unpacker_impl: Failed to open file wspr_spots.txt\n";
        exit(EXIT_FAILURE);
      }
      //if ((ftimer=fopen("wspr_timer.out","r"))) {
         //Accumulate timing data
      //   nr=fscanf(ftimer,"%f %f %f %f %f %f %f",
      //      &treadwav,&tcandidates,&tsync0,&tsync1,&tsync2,&tfano,&ttotal);
      //   fclose(ftimer);
      //}
      //ftimer=fopen("wspr_timer.out","w");
      callsign=(char *)malloc(sizeof(char)*13);
      call_loc_pow=(char *)malloc(sizeof(char)*23);
      // read hash table
      hashtab=(char *)malloc(sizeof(char)*32768*13);
      memset(hashtab,0,sizeof(char)*32768*13);
      int nh;
      usehashtable=1;
      if( usehashtable ) {
         char line[80], hcall[12];
         if( (fhash=fopen("hashtable.txt","r+")) ) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
               sscanf(line,"%d %s",&nh,hcall);
               strcpy(hashtab+nh*13,hcall);
            }
         } else {
           fhash=fopen("hashtable.txt","w+");
            if (fhash==NULL) {
              std::cerr << "In WSPR_unpacker_impl: Failed to open file hashtable.txt\n";
              exit(EXIT_FAILURE);
            }
         }
         fclose(fhash);
      }
    }

    // destructor.
    WSPR_unpacker_impl::~WSPR_unpacker_impl()
    {
      fclose(fall_wspr);
      fclose(fwsprd);
      //fclose(ftimer);
      if( usehashtable ) {
         fhash=fopen("hashtable.txt","w");
         if (fhash==NULL) {
           std::cerr << "In WSPR_unpacker_impl: Failed to open file hashtable.txt\n";
           exit(EXIT_FAILURE);
         }
         for (int i=0; i<32768; i++) {
            if( strncmp(hashtab+i*13,"\0",1) != 0 ) {
               fprintf(fhash,"%5d %s\n",i,hashtab+i*13);
            }
         }
         fclose(fhash);
      }
    }

    void WSPR_unpacker_impl::unpack(pmt::pmt_t msg)
    {
      // extract input pdu
      pmt::pmt_t meta(pmt::car(msg));
      pmt::pmt_t vector(pmt::cdr(msg));
      // get payload size in bytes (unpacked, each byte contains one code bit)
      size_t len = pmt::blob_length(vector);
      std::cout << "pdu_length = " << len << std::endl;

      printf("--- WSPR_unpacker_impl is printing unpacked message\n");
      size_t offset(0);
      // convert polymorphic data to array
      // of unsigned chars (8 bits, range 0..255)
      char* d =
         (char*) pmt::uniform_vector_elements(vector, offset);
      //for(size_t i=0; i<len; i+=16){
      //  printf("%04x: ", ((unsigned int)i));
      //  for(size_t j=i; j<std::min(i+16,len); j++){
      //    printf("%02x ",d[j] );
      //  }
      //  std::cout << std::endl;
      //}
      int noprint=0;
      noprint=unpk_(d,hashtab,call_loc_pow,callsign);
      if( !noprint) {
         printf("WSPR_unpacker_impl::unpack: %-s \n", call_loc_pow);
       }
      // prepare output vector
      std::vector<unsigned char> vec(len);
      for (int i=0; i<len; i++) {
         vec[i] = d[i];
      }
      // send the vector
      pmt::pmt_t vecpmt(pmt::make_blob(&vec[0], len));
      pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, vecpmt));
      message_port_pub(out_port, pdu);
    }

    void
    WSPR_unpacker_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      /* <+forecast+> e.g. ninput_items_required[0] = noutput_items */
    }

    int
    WSPR_unpacker_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      //const <+ITYPE+> *in = (const <+ITYPE+> *) input_items[0];
      //<+OTYPE+> *out = (<+OTYPE+> *) output_items[0];

      // Do <+signal processing+>
      // Tell runtime system how many input items we consumed on
      // each input stream.
      consume_each (noutput_items);

      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace uwspr */
} /* namespace gr */
