/*! \file
 *  \brief Derivative displacements
 */

#include "chromabase.h"

#include "deriv_quark_displacement_w_ihep.h"
#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/displace.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"

namespace Chroma {
namespace DerivQuarkDisplacementEnv {
namespace {
//! Determine sign of plusminus
/*!
 * \ingroup sources
 */
int plusMinus(enum PlusMinus isign) {
  int is = 0;
  switch (isign) {
  case PLUS:
    is = +1;
    break;

  case MINUS:
    is = -1;
    break;

  default:
    QDP_error_exit("illegal isign in plusminus");
  }
  return is;
}

//! Construct (RhoxNabla_E) source
QuarkDisplacement<LatticePropagator> *
mesRhoxNablaEDisplace(XMLReader &xml_in, const std::string &path) {
  return new MesRhoxNablaEDisplace<LatticePropagator>(ParamsDir(xml_in, path));
}
}

// Construct (RhoxNabla_E) source
template <>
void MesRhoxNablaEDisplace<LatticePropagator>::
operator()(LatticePropagator &tmp, const multi1d<LatticeColorMatrix> &u,
           enum PlusMinus isign) const {
  START_CODE();

  LatticePropagator fin = zero;
  int length = plusMinus(isign) * params.deriv_length;

  // \f$\Gamma_f \equiv S_{\alpha jk}\gamma_j D_k\f$
  for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k) {
      Real e = ETensor3d(params.deriv_dir, j, k);
      if (toBool(e != 0.0))
        fin += e * (Gamma(1 << j) * rightNabla(tmp, u, k, length));
    }

  tmp = fin;

  END_CODE();
}
//! Local registration flag
static bool registered = false;
namespace IHEP {
//! Register all the possible deriv mesons
bool registerAll() {
  bool success = true;
  if (!registered) {
    //! Register all the factories

    success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(
        std::string("RHOxNABLA_E-DERIV"), mesRhoxNablaEDisplace);
    registered = true;
  }
  return success;
}
}
}
}
