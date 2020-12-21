// -*- C++ -*-
// $Id: inline_propagator_w.h,v 3.5 2007-08-23 19:02:45 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_lineardivergence_h__
#define __inline_lineardivergence_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma
{
  /*! \ingroup inlinehadron */
  namespace InlineLinearDivergenceEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineLinearDivergenceParams
  {
    InlineLinearDivergenceParams();
    InlineLinearDivergenceParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    multi1d<int> grid;
    multi1d<int> link_lengths;
    int nsets;

    ChromaProp_t      param;

    GroupXML_t inv_param_h;
    GroupXML_t inv_param_l;
    int HL_ratio;

    struct NamedObject_t
    {
      std::string     gauge_id;
    } named_obj;
    bool use_ckpoint;

    std::string xml_file;  // Alternate XML file pattern
    std::string iog_file;
    int cfg_serial;
  };

  //! Inline propagator calculation
  /*! \ingroup inlinehadron */
  class InlineLinearDivergence : public AbsInlineMeasurement
  {
  public:
    ~InlineLinearDivergence() {}
    InlineLinearDivergence(const InlineLinearDivergenceParams& p) : params(p) {}
    InlineLinearDivergence(const InlineLinearDivergence& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out);

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out);

  private:
    InlineLinearDivergenceParams params;
  };

}

#endif
