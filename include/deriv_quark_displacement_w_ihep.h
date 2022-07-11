#ifndef __DERIV_QUARK_DISPLACEMENT_W_H_
#define __DERIV_QUARK_DISPLACEMENT_W_H_

#include "meas/smear/deriv_quark_displacement_w.h"
namespace Chroma {
namespace DerivQuarkDisplacementEnv {
namespace IHEP {
bool registerAll();
}
//! Construct (RhoxNabla_E) source
/*!
 * \ingroup sources
 *
 * Operator is  rho x nabla_E
 * The sink interpolator is
 * \f$\Gamma_f \equiv S_{ijk}\gamma_j D_k\f$
 */
template <typename T>
class MesRhoxNablaEDisplace : public QuarkDisplacement<T> {
public:
  //! Full constructor
  MesRhoxNablaEDisplace(const ParamsDir &p) : params(p) {}

  //! Displace the quark
  void operator()(T &quark, const multi1d<LatticeColorMatrix> &u,
                  enum PlusMinus isign) const;

private:
  ParamsDir params; /*!< source params */
};
}
}
#endif
