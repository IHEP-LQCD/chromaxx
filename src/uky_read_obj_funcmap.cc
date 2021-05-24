/*! \file
 *  \brief Read object function std::map
 */

#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"
#include "uky_read_obj_funcmap.h"
#include "kyuqprop_io.h"

namespace Chroma {

//! IO function std::map environment
/*! \ingroup inlineio */
namespace UKYReadObjCallMapEnv {
// Anonymous namespace
namespace {
//! Read a propagator
void UKYReadLatProp(const std::string &buffer_id, const std::string &file) {
  LatticePropagator obj;
  XMLReader record_xml;

  /////////////////////////////////////////////
  XMLBufferWriter xml_buf;
  push(xml_buf, "Propagator");
  write(xml_buf, "Kappa", 0.2);

  push(xml_buf, "descendant");
  write(xml_buf, "j_decay", 3);
  pop(xml_buf);

  push(xml_buf, "ForwardProp");
  write(xml_buf, "version", 10);
  push(xml_buf, "InvertParam");
  write(xml_buf, "invType", "CG");
  pop(xml_buf);
  write(xml_buf, "boundary", "1 1 1 -1");
  write(xml_buf, "FermAct", "CLOVER");
  push(xml_buf, "FermionAction");
  pop(xml_buf);
  pop(xml_buf);

  push(xml_buf, "PropSource");
  write(xml_buf, "version", 6);
  pop(xml_buf);

  pop(xml_buf);

  try {
// Temporary XLC failure workaround
#if 0
	  record_xml.open(xml_buf);
#else
    const std::string bufcontent = xml_buf.str() + "\n";
    std::istringstream is(bufcontent);
    record_xml.open(is);
#endif
  }
  catch (const std::string &e) {
    QDP_error_exit("Error in readukyqprop: %s", e.c_str());
  }
  /////////////////////////////////////////////

  readKYUQprop2(obj, file);

  XMLBufferWriter file_xml;
  push(file_xml, "UKY");
  pop(file_xml);

  TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
  TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
  TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
  TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
}

//! Local registration flag
bool registered = false;

} // end anonymous namespace

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= TheUKYReadObjFuncMap::Instance().registerFunction(
        std::string("LatticePropagator"), UKYReadLatProp);
    registered = true;
  }
  return success;
}
}
}
