/* -*- c++ -*- */

#define UWSPR_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "uwspr_swig_doc.i"

%{
#include "uwspr/c2file_source.h"
#include "uwspr/decoder.h"
#include "uwspr/WSPR_unpacker.h"
%}


%include "uwspr/c2file_source.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, c2file_source);
%include "uwspr/decoder.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, decoder);
%include "uwspr/WSPR_unpacker.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, WSPR_unpacker);
