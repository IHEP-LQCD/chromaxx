// -*- C++ -*-
/*! \file
 *  \brief Fixed gauge boundary conditions at the edge with unit link
 */

#ifndef __fixed_gaugebc_h__
#define __fixed_gaugebc_h__

#include "gaugebc.h"

namespace Chroma {

/*! @ingroup gaugebcs */
namespace FixedGaugeBCEnv {
extern const std::string name;
bool registerAll();
}

/*! @ingroup gaugebcs */
struct FixedGaugeBCParams {
  FixedGaugeBCParams();
  FixedGaugeBCParams(XMLReader &xml, const std::string &path);
  multi1d<Complex> boundary;
  int maxlink;
};

/*! @ingroup gaugebcs */
void read(XMLReader &xml, const std::string &path, FixedGaugeBCParams &p);

/*! @ingroup gaugebcs */
void write(XMLWriter &xml, const std::string &path,
           const FixedGaugeBCParams &p);

//! Concrete class for gauge actions with fixed boundary conditions
/*! @ingroup gaugebcs
 *
 *  the links on the edge of the lattice is set fixed
 */
template <typename P, typename Q> class FixedGaugeBC : public GaugeBC<P, Q> {
public:
  FixedGaugeBC(const multi1d<Complex> &boundary_) : boundary(boundary_) {}

  //! From param struct
  FixedGaugeBC(const FixedGaugeBCParams &p)
      : boundary(p.boundary), maxlink(p.maxlink) {}

  //! Destructor is automatic
  ~FixedGaugeBC() {}

  //! Modify U fields in place
  void modify(Q &u) const {
    Q BndField;
    BndField.resize(Nd);

    LatticeBoolean mask_m = false;
    for (int m = 0; m < Nd; ++m) {
      LatticeInteger coord_m = Layout::latticeCoordinate(m);

      mask_m |= (coord_m == 0);
      mask_m |= (coord_m == (QDP::Layout::lattSize()[m] - 1));
      switch (maxlink) {
      case 1:
        break;

      case 2:
        mask_m |= (coord_m == 1);
        mask_m |= (coord_m == (QDP::Layout::lattSize()[m] - 2));
        break;

      default:
        QDP_error_exit("FixedGaugeBC: unsupported maxlink = %d\n", maxlink);
      }
    }

    for (int m = 0; m < Nd; ++m) {
      BndField[m] = boundary[m]; // set boundary to a diagnoal matrix

      copymask(u[m], mask_m, BndField[m]);
    }
  }

  //! Zero the U fields in place on the masked links
  void zero(P &u) const {
    LatticeColorMatrix BndField = QDP::zero;

    for (int m = 0; m < Nd; ++m) {
      LatticeInteger coord_m = Layout::latticeCoordinate(m);
      LatticeBoolean mask_m = false;

      mask_m = (coord_m == 0);
      mask_m |= (coord_m == (QDP::Layout::lattSize()[m] - 1));
      switch (maxlink) {
      case 1:
        break;

      case 2:
        mask_m |= (coord_m == 1);
        mask_m |= (coord_m == (QDP::Layout::lattSize()[m] - 2));
        break;

      default:
        QDP_error_exit("FixedGaugeBC: unsupported maxlink = %d\n", maxlink);
      }

      copymask(u[m], mask_m, BndField);
    }
  }

  //! Says if there are non-trivial BC links
  /*!
   * This is a simple implementation - always do the work.
   * Could be improved by checking boundary
   */
  bool nontrivialP() const { return true; }

private:
  // Hide empty constructor
  FixedGaugeBC() {}

  //! Hide assignment
  void operator=(const FixedGaugeBC<P, Q> &a) {}

private:
  multi1d<Complex> boundary;
  int maxlink;
};
}

#endif
