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

#ifndef INCLUDED_HELPERS_H
#define INCLUDED_HELPERS_H

#include <string.h>
#include <stdint.h>

namespace gr {
   namespace uwspr {

class helpers
{
   public:
      uint32_t nhash( const void * key, size_t length, uint32_t initval);
      void unpack50(  char *dat, int32_t *n1, int32_t *n2 );
      int unpackcall( int32_t ncall, char *call );
      int unpackgrid( int32_t ngrid, char *grid);
      int unpackpfx( int32_t nprefix, char *call);
      int unpk_( char *message, char *hashtab,
                char *call_loc_pow, char *callsign);
};

   }   // namespace uwspr
} // namespace gr

#endif /* INCLUDED_HELPERS_H */
