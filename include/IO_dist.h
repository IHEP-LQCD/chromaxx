#include "key_timeslice_colorvec.h"
#include "qdp_disk_map_slice.h"
#include "qdp_map_obj_disk.h"
namespace RH_qcd {

multi1d<LatticeColorVector> read_dist_vec(std::string filename, int Nt,
                                          int Num_vecs) {
  multi1d<LatticeColorVector> vec_array(Num_vecs);
  QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> >
  eigen_source;
  eigen_source.open(filename);
  std::string user_str;
  eigen_source.getUserdata(user_str);

  for (int i = 0; i < Num_vecs; i++) {
    LatticeColorVectorF tmpvec;
    for (int t = 0; t < Nt; t++) {
      KeyTimeSliceColorVec_t key_vec;
      key_vec.t_slice = t;
      key_vec.colorvec = i;
      TimeSliceIO<LatticeColorVectorF> time_slice_io(tmpvec, t);
      eigen_source.get(key_vec, time_slice_io);
    }
    vec_array[i] = tmpvec;
    // QDPIO::cout<<"norm(vec)="<<sqrt(norm2(tmpvec))<<std::endl;
  }

  return vec_array;
}
}
