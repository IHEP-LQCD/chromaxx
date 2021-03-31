// -*- C++ -*-

#ifndef __inline_kyu_prop_h__
#define __inline_kyu_prop_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/kyuqprop_io.h"
namespace Chroma
{
  namespace InlineKyuPropEnv
  {
    extern const std::string name;
    bool registerAll();
  } 
  struct InlineKyuPropParams
  { 
    InlineKyuPropParams();
    InlineKyuPropParams(XMLReader& xml_in, const std::string& path);

    unsigned long     frequency;

struct Param_t
{
    param
} param;
    struct NamedObject_t
    {
      std::string gauge_id; /*!< Input gauge field */
    } named_obj;


    std::string xml_file;  // Alternate XML file pattern
    std::string iog_file;
    int cfg_serial;
  };

  class InlineKyuProp : public AbsInlineMeasurement
  {
  public:
    ~InlineKyuProp() {}
    InlineKyuProp(const InlineKyuPropParams& p) : params(p) {}
    InlineKyuProp(const InlineKyuProp& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
            XMLWriter& xml_out);

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
          XMLWriter& xml_out);

  private:
    InlineKyuPropParams params;
  };
}
#endif