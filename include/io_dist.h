#ifndef __IO_DIST_H_
#define __IO_DIST_H_

#include "key_timeslice_colorvec.h"
#include "qdp_disk_map_slice.h"
#include "qdp_map_obj_disk.h"
namespace IHEP {

multi1d<LatticeColorVector> read_dist_vec(std::string filename, int Nt,
                                          int Num_vecs);
}
#endif
