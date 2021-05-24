#ifndef __register_h__
#define __register_h__

// New measurement headers go here

#include "inline_linear_divergence2.h"
#include "inline_tests.h"
#include "inline_uky_read_obj.h"
#include "inline_cluster_dec.h"
#include "inline_linear_divergence.h"
#include "inline_uky_write_obj.h"
#include "aniso_spectrum_gaugeact2.h"

namespace Chroma {
void register_new(bool &);
}

#endif