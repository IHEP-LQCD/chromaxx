/*! \file
 * \brief Key for time-sliced color eigenvectors
 */

#ifndef __key_timeslice_colorvec_h__
#define __key_timeslice_colorvec_h__

#include "qdp_map_obj_disk.h"

namespace RH_qcd {
using namespace QDP;

//----------------------------------------------------------------------------
//! Prop operator
/*! \ingroup ferm */
struct KeyTimeSliceColorVec_t {
  int t_slice;  /*!< Source time slice */
  int colorvec; /*!< Colorstd::vector index */
};

//----------------------------------------------------------------------------
//! Diagnostics
StandardOutputStream &operator<<(StandardOutputStream &os,
                                 const KeyTimeSliceColorVec_t &param);

//----------------------------------------------------------------------------
//! KeyTimeSliceColorVec read
void read(BinaryReader &bin, KeyTimeSliceColorVec_t &param);

//! KeyTimeSliceColorVec write
void write(BinaryWriter &bin, const KeyTimeSliceColorVec_t &param);

//! KeyTimeSliceColorVec reader
void read(XMLReader &xml, const std::string &path,
          KeyTimeSliceColorVec_t &param);

//! KeyTimeSliceColorVec writer
void write(XMLWriter &xml, const std::string &path,
           const KeyTimeSliceColorVec_t &param);

//----------------------------------------------------------------------------
//! Diagnostics
StandardOutputStream &operator<<(StandardOutputStream &os,
                                 const KeyTimeSliceColorVec_t &param) {
  os << "KeyTimeSliceColorVec_t:";
  os << " t_slice= " << param.t_slice;
  os << " colorvec= " << param.colorvec;
  os << std::endl;

  return os;
}

//----------------------------------------------------------------------------
// KeyTimeSliceDist read
void read(BinaryReader &bin, KeyTimeSliceColorVec_t &param) {
  read(bin, param.t_slice);
  read(bin, param.colorvec);
}

// KeyTimeSliceDist write
void write(BinaryWriter &bin, const KeyTimeSliceColorVec_t &param) {
  write(bin, param.t_slice);
  write(bin, param.colorvec);
}

//! KeyTimeSliceDist reader
void read(XMLReader &xml, const std::string &path,
          KeyTimeSliceColorVec_t &param) {
  XMLReader paramtop(xml, path);

  read(paramtop, "t_slice", param.t_slice);
  read(paramtop, "colorvec", param.colorvec);
}

// KeyTimeSliceDist writer
void write(XMLWriter &xml, const std::string &path,
           const KeyTimeSliceColorVec_t &param) {
  push(xml, path);

  write(xml, "t_slice", param.t_slice);
  write(xml, "colorvec", param.colorvec);

  pop(xml);
}
}
#endif
