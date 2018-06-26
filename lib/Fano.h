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

#ifndef INCLUDED_FANO_H
#define INCLUDED_FANO_H

extern unsigned char Partab[];

namespace gr {
  namespace uwspr {
    class Fano
    {
    public:
       // Fano FEC structure
       struct node {
          unsigned long encstate; /* Encoder state of next node */
          long gamma; /* Cumulative metric to this node */
          int metrics[4]; /* Metrics indexed by all possible tx syms */
          int tm[2]; /* Sorted metrics for current hypotheses */
          int i; /* Current branch being tested */
       };
       // metric table
       float bias; //Fano metric bias (used for both Fano and stack algorithms)
       int mettab[2][256];
       // Convolutional FEC (Fano) methods
       int encode(
          unsigned char *symbols, // Output buffer, 2*nbytes*8
          unsigned char *data, // Input buffer, nbytes
          unsigned int nbytes); // Number of bytes in data
       // FANO FEC decoder
       int fano(unsigned int *metric, unsigned int *cycles, unsigned int *maxnp,
          unsigned char *data,unsigned char *symbols, unsigned int nbits,
            int mettab[2][256],int delta,unsigned int maxcycles);
        Fano();
        ~Fano();
    };
  } // namespace uwspr
} // namespace gr
#endif /* INCLUDED_FANO_H */
