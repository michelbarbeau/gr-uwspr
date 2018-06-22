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


#ifndef INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_H
#define INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_H

#include <uwspr/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace uwspr {

    /*!
     * \brief Accepts a continous stream of time-domain samples.
     *   The samples are pushed in ring buffer.
     *   When 120~seconds of samples are buffered,
     *   they become the payload of a PDU that it posted on the output port.
     * \ingroup uwspr
     *
     */
    class UWSPR_API sliding_window_stream_to_pdu : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<sliding_window_stream_to_pdu> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of uwspr::sliding_window_stream_to_pdu.
       *
       * To avoid accidental use of raw pointers, uwspr::sliding_window_stream_to_pdu's
       * constructor is in a private implementation
       * class. uwspr::sliding_window_stream_to_pdu::make is the public interface for
       * creating new instances.
       */
      static sptr make(int fs, int fl, int shift, int C);
    };

  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_H */
