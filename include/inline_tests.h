// -*- C++ -*-

#ifndef __inline_tests_h__
#define __inline_tests_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma {
namespace InlineTestsEnv {
extern const std::string name;
bool registerAll();
}
struct InlinetestsParams {
  InlinetestsParams();
  InlinetestsParams(XMLReader &xml_in, const std::string &path);
  void writeXML(XMLWriter &xml_out, const std::string &path);

  unsigned long frequency;

  struct NamedObject_t {
    std::string gauge_id; /*!< Input gauge field */

  } named_obj;
  struct Param_t {
    int length; /*!< Input param */
  } param;

  std::string xml_file; // Alternate XML file pattern
  std::string iog_file;
  int cfg_serial;
};

class Inlinetests : public AbsInlineMeasurement {
public:
  ~Inlinetests() {}
  Inlinetests(const InlinetestsParams &p) : params(p) {}
  Inlinetests(const Inlinetests &p) : params(p.params) {}

  unsigned long getFrequency(void) const { return params.frequency; }

  //! Do the measurement
  void operator()(const unsigned long update_no, XMLWriter &xml_out);

protected:
  //! Do the measurement
  void func(const unsigned long update_no, XMLWriter &xml_out);

private:
  InlinetestsParams params;
};
}
#endif