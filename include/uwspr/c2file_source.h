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


#ifndef INCLUDED_UWSPR_C2FILE_SOURCE_H
#define INCLUDED_UWSPR_C2FILE_SOURCE_H

#include <uwspr/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace uwspr {

    /*!
     * \brief <+description of block+>
     * \ingroup uwspr
     *
     */
    class UWSPR_API c2file_source : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<c2file_source> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of uwspr::c2file_source.
       *
       * To avoid accidental use of raw pointers, uwspr::c2file_source's
       * constructor is in a private implementation
       * class. uwspr::c2file_source::make is the public interface for
       * creating new instances.
       * \param filename The .c2 file to be opened
       */
      static sptr make(const char *filename, bool repeat = false);
    };

  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_C2FILE_SOURCE_H */
