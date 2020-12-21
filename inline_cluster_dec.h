// -*- C++ -*-

#ifndef __inline_clusterdec_h__
#define __inline_clusterdec_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace InlineClusterDecEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineClusterDecParams
  {
    InlineClusterDecParams();
    InlineClusterDecParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct NamedObject_t
    {
      std::string gauge_id; /*!< Input gauge field */
    } named_obj;


    std::string xml_file;  // Alternate XML file pattern
    std::string iog_file;
    int cfg_serial;
    int j_decay;
    int bl_level;
  };

  //! Inline propagator calculation
  /*! \ingroup inlinehadron */
  class InlineClusterDec : public AbsInlineMeasurement
  {
  public:
    ~InlineClusterDec() {}
    InlineClusterDec(const InlineClusterDecParams& p) : params(p) {}
    InlineClusterDec(const InlineClusterDec& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out);

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out);

  private:
    InlineClusterDecParams params;
  };

}

#endif