/* -*- c++ -*- */

#define UWSPR_API

%include "gnuradio.i"			// the common stuff

//load generated python docstrings
%include "uwspr_swig_doc.i"

%{
#include "uwspr/c2file_source.h"
#include "uwspr/WSPR_unpacker.h"
#include "uwspr/sync_and_demodulate.h"
#include "uwspr/sliding_window_stream_to_pdu.h"
#include "uwspr/FDR.h"
%}


%include "uwspr/c2file_source.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, c2file_source);

%include "uwspr/WSPR_unpacker.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, WSPR_unpacker);
%include "uwspr/sync_and_demodulate.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, sync_and_demodulate);
%include "uwspr/sliding_window_stream_to_pdu.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, sliding_window_stream_to_pdu);
%include "uwspr/FDR.h"
GR_SWIG_BLOCK_MAGIC2(uwspr, FDR);
