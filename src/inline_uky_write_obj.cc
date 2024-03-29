/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "inline_uky_write_obj.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "uky_write_obj_funcmap.h"

namespace Chroma {
namespace InlineUKYWriteNamedObjEnv {
namespace {
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineMeas(Params(xml_in, path));
}

//! Local registration flag
bool registered = false;

const std::string name = "UKY_WRITE_NAMED_OBJECT";
} // namespace

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    // Datatype writer
    success &= UKYWriteObjCallMapEnv::registerAll();

    // Inline measurement
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);

    registered = true;
  }
  return success;
}

//! Object buffer
void write(XMLWriter &xml, const std::string &path,
           const Params::NamedObject_t &input) {
  push(xml, path);

  write(xml, "object_id", input.object_id);
  write(xml, "object_type", input.object_type);

  pop(xml);
}

//! File output
void write(XMLWriter &xml, const std::string &path,
           const Params::File_t &input) {
  push(xml, path);

  write(xml, "file_name", input.file_name);

  pop(xml);
}

//! Object buffer
void read(XMLReader &xml, const std::string &path,
          Params::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "object_id", input.object_id);
  read(inputtop, "object_type", input.object_type);
}

//! File output
void read(XMLReader &xml, const std::string &path, Params::File_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "file_name", input.file_name);
}

// Param stuff
Params::Params() { frequency = 0; }

Params::Params(XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Parameters for source construction
    read(paramtop, "NamedObject", named_obj);

    // Read in the destination
    read(paramtop, "File", file);
  } catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

void Params::writeXML(XMLWriter &xml_out, const std::string &path) {
  push(xml_out, path);

  // Parameters for source construction
  write(xml_out, "NamedObject", named_obj);

  // Write out the destination
  write(xml_out, "File", file);

  pop(xml_out);
}

void InlineMeas::operator()(unsigned long update_no, XMLWriter &xml_out) {
  START_CODE();

  push(xml_out, "uky_write_named_obj");
  write(xml_out, "update_no", update_no);

  QDPIO::cout << name << ": object writer" << std::endl;
  StopWatch swatch;

  // Write the object
  // ONLY UKY output format is supported in this task
  // Other tasks could support other disk formats
  QDPIO::cout << "Attempt to write object name = " << params.named_obj.object_id
              << std::endl;
  write(xml_out, "object_id", params.named_obj.object_id);
  try {
    swatch.reset();

    // Write the object
    swatch.start();
    UKYWriteObjCallMapEnv::TheUKYWriteObjFuncMap::Instance().callFunction(
        params.named_obj.object_type, params.named_obj.object_id,
        params.file.file_name);
    swatch.stop();

    QDPIO::cout << "Object successfully written: time= "
                << swatch.getTimeInSeconds() << " secs" << std::endl;
  } catch (std::bad_cast) {
    QDPIO::cerr << name << ": cast error" << std::endl;
    QDP_abort(1);
  } catch (const std::string &e) {
    QDPIO::cerr << name << ": error message: " << e << std::endl;
    QDP_abort(1);
  }

  QDPIO::cout << name << ": ran successfully" << std::endl;

  pop(xml_out); // write_named_obj

  END_CODE();
}
} // namespace InlineUKYWriteNamedObjEnv
} // namespace Chroma
