/*! \file
 *  \brief Fixed gauge boundary conditions at the edge with unit link
 */

#include "fixed_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma {

namespace FixedGaugeBCEnv {
//! Calllback function to register with the factory
GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>> *
createGaugeBC(XMLReader &xml, const std::string &path) {
  QDPIO::cout << "Factory Callback: Creating FixedGaugeBC " << std::endl;
  return new FixedGaugeBC<multi1d<LatticeColorMatrix>,
                          multi1d<LatticeColorMatrix>>(
      FixedGaugeBCParams(xml, path));
}

const std::string name = "FIXED_GAUGEBC";

//! Local registration flag
static bool registered = false;

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &=
        TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
    registered = true;
  }
  return success;
}
} // namespace FixedGaugeBCEnv

FixedGaugeBCParams::FixedGaugeBCParams(XMLReader &xml,
                                       const std::string &path) {

  XMLReader paramtop(xml, path);

  try {
    read(paramtop, "./boundary", boundary);
    read(paramtop, "./maxlink", maxlink);
  } catch (const std::string &e) {
    QDPIO::cerr << "Error reading XML: " << e << std::endl;
    QDP_abort(1);
  }

  QDPIO::cout << "Creating FixedGaugeBCParams with boundary: ";
  QDPIO::cout << "Boundary.size = " << boundary.size() << std::endl;
  for (int i = 0; i < boundary.size(); i++) {
    QDPIO::cout << boundary[i] << std::endl;
  }
  QDPIO::cout << std::endl;
}

void read(XMLReader &xml, const std::string &path, FixedGaugeBCParams &p) {
  FixedGaugeBCParams tmp(xml, path);
  p = tmp;
}

void write(XMLWriter &xml, const std::string &path,
           const FixedGaugeBCParams &p) {
  push(xml, path);
  write(xml, "boundary", p.boundary);
  write(xml, "maxlink", p.maxlink);
  pop(xml);
}

} // End namespace Chroma
