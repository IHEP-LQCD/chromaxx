// -*- C++ -*-
// $Id: inline_propagator_w.h,v 3.5 2007-08-23 19:02:45 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_lineardivergenc2_h__
#define __inline_lineardivergenc2_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace InlineLinearDivergence2Env
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineLinearDivergence2Params
  {
    InlineLinearDivergence2Params();
    InlineLinearDivergence2Params(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    struct NamedObject_t
    {
      std::string gauge_id; /*!< Input gauge field */

      struct Props_t
      {
        std::string first_id;
        std::string second_id;
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;


    std::string xml_file;  // Alternate XML file pattern
    std::string iog_file;
    int cfg_serial;
  };

  //! Inline propagator calculation
  /*! \ingroup inlinehadron */
  class InlineLinearDivergence2 : public AbsInlineMeasurement
  {
  public:
    ~InlineLinearDivergence2() {}
    InlineLinearDivergence2(const InlineLinearDivergence2Params& p) : params(p) {}
    InlineLinearDivergence2(const InlineLinearDivergence2& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out);

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out);

  private:
    InlineLinearDivergence2Params params;
  };

}

#endif
