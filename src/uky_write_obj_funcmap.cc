/*! \file
 *  \brief Write object function std::map
 */

#include "uky_write_obj_funcmap.h"
#include "kyuqprop_io.h"
#include "meas/inline/io/named_objmap.h"
#include "named_obj.h"

namespace Chroma {

//! IO function std::map environment
/*! \ingroup inlineio */
namespace UKYWriteObjCallMapEnv {
// Anonymous namespace
namespace {
//! Write a propagator
void UKYWriteLatProp(const std::string &buffer_id, const std::string &file) {
  LatticePropagator obj =
      TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
  writeKYUQprop2(obj, file);
}

//! Local registration flag
bool registered = false;

} // end namespace

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= TheUKYWriteObjFuncMap::Instance().registerFunction(
        std::string("LatticePropagator"), UKYWriteLatProp);
    registered = true;
  }
  return success;
}
} // namespace UKYWriteObjCallMapEnv
} // namespace Chroma
