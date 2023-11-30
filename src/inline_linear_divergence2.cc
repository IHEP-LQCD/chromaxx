// $Id: inline_propagator_w.cc,v 3.11 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "inline_linear_divergence2.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "fermact.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"

#include "meas/inline/io/named_objmap.h"

#include "meas/smear/no_quark_displacement.h"

#include "io_general_class.h"
#include "util/ferm/transf.h"
#include <time.h>

namespace Chroma {
namespace InlineLinearDivergence2Env {
namespace {
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineLinearDivergence2(
      InlineLinearDivergence2Params(xml_in, path));
}

//! Local registration flag
bool registered = false;
} // namespace

const std::string name = "LinearDivergence2";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    registered = true;
  }
  return success;
}
} // namespace InlineLinearDivergence2Env
//! Propagator input
void read(XMLReader &xml, const std::string &path,
          InlineLinearDivergence2Params::NamedObject_t::Props_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "first_id", input.first_id);
  read(inputtop, "second_id", input.second_id);
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const InlineLinearDivergence2Params::NamedObject_t::Props_t &input) {
  push(xml, path);

  write(xml, "first_id", input.first_id);
  write(xml, "second_id", input.second_id);

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path,
          InlineLinearDivergence2Params::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "gauge_id", input.gauge_id);
  read(inputtop, "sink_pairs", input.sink_pairs);
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const InlineLinearDivergence2Params::NamedObject_t &input) {
  push(xml, path);

  write(xml, "gauge_id", input.gauge_id);
  write(xml, "sink_pairs", input.sink_pairs);

  pop(xml);
}

// Param stuff
InlineLinearDivergence2Params::InlineLinearDivergence2Params() {
  frequency = 0;
}

InlineLinearDivergence2Params::InlineLinearDivergence2Params(
    XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    read(paramtop, "cfg_serial", cfg_serial);

    // Read in the output propagator/source configuration info
    read(paramtop, "NamedObject", named_obj);

    read(paramtop, "iog_file", iog_file);
  } catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

void InlineLinearDivergence2Params::writeXML(XMLWriter &xml_out,
                                             const std::string &path) {
  push(xml_out, path);

  write(xml_out, "NamedObject", named_obj);

  pop(xml_out);
}

// Anonymous namespace
namespace {
//! Useful structure holding sink props
struct SinkPropContainer_t {
  ForwardProp_t prop_header;
  std::string sink_prop_1_id;
  Real Mass;

  multi1d<int> bc;

  std::string source_type;
  std::string source_disp_type;
  std::string sink_type;
  std::string sink_disp_type;
};

//! Useful structure holding sink props
struct AllSinkProps_t {
  SinkPropContainer_t sink_prop_1;
  SinkPropContainer_t sink_prop_2;
};

//! Read a sink prop
void readSinkProp(SinkPropContainer_t &s, const std::string &id) {
  try {
    // Try a cast to see if it succeeds
    const LatticePropagator &foo =
        TheNamedObjMap::Instance().getData<LatticePropagator>(id);

    // Snarf the data into a copy
    s.sink_prop_1_id = id;

    // Snarf the prop info. This is will throw if the prop_id is not there
    XMLReader prop_file_xml, prop_record_xml;
    TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
    TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);

    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    {
      std::string xpath;
      read(prop_record_xml, "/SinkSmear", s.prop_header);

      read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType",
           s.source_type);
      xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
      if (prop_record_xml.count(xpath) != 0)
        read(prop_record_xml, xpath, s.source_disp_type);
      else
        s.source_disp_type = NoQuarkDisplacementEnv::getName();

      read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", s.sink_type);
      xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
      if (prop_record_xml.count(xpath) != 0)
        read(prop_record_xml, xpath, s.sink_disp_type);
      else
        s.sink_disp_type = NoQuarkDisplacementEnv::getName();
    }
  } catch (std::bad_cast) {
    QDPIO::cerr << InlineLinearDivergence2Env::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  } catch (const std::string &e) {
    QDPIO::cerr << InlineLinearDivergence2Env::name << ": error message: " << e
                << std::endl;
    QDP_abort(1);
  }

  // Derived from input prop
  // Hunt around to find the mass
  // NOTE: this may be problematic in the future if actions are used with no
  // clear def. of a Mass
  QDPIO::cout << "Try action and mass" << std::endl;
  s.Mass = getMass(s.prop_header.prop_header.fermact);

  // Only baryons care about boundaries
  // Try to find them. If not present, assume dirichlet.
  // This turns off any attempt to time reverse which is the
  // only thing that the BC are affecting.
  s.bc.resize(Nd);
  s.bc = 0;

  try {
    s.bc = getFermActBoundary(s.prop_header.prop_header.fermact);
  } catch (const std::string &e) {
    QDPIO::cerr << InlineLinearDivergence2Env::name
                << ": caught exception. No BC found in these headers. Will "
                   "assume dirichlet: "
                << e << std::endl;
  }

  QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id
              << std::endl;
  QDPIO::cout << "Mass = " << s.Mass << std::endl;
}

//! Read all sinks
void readAllSinks(
    AllSinkProps_t &s,
    InlineLinearDivergence2Params::NamedObject_t::Props_t sink_pair) {
  QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.first_id
              << std::endl;
  readSinkProp(s.sink_prop_1, sink_pair.first_id);
  QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

  QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.second_id
              << std::endl;
  readSinkProp(s.sink_prop_2, sink_pair.second_id);
  QDPIO::cout << "Forward propagator successfully parsed" << std::endl;
}

} // namespace

// Function call
void InlineLinearDivergence2::operator()(unsigned long update_no,
                                         XMLWriter &xml_out) {
  // If xml file not empty, then use alternate
  if (params.xml_file != "") {
    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);
    write(xml_out, "xml_file", xml_file);
    pop(xml_out);

    XMLFileWriter xml(xml_file);
    func(update_no, xml);
  } else {
    func(update_no, xml_out);
  }
}

// Real work done here
void InlineLinearDivergence2::func(unsigned long update_no,
                                   XMLWriter &xml_out) {
  START_CODE();

  QDPIO::cout << InlineLinearDivergence2Env::name
              << ": linear divergence calculation" << std::endl;

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  // Test and grab a reference to the gauge field
  XMLBufferWriter gauge_xml;
  try {
    TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(
        params.named_obj.gauge_id);
    TheNamedObjMap::Instance()
        .get(params.named_obj.gauge_id)
        .getRecordXML(gauge_xml);
  } catch (std::bad_cast) {
    QDPIO::cerr << InlineLinearDivergence2Env::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  } catch (const std::string &e) {
    QDPIO::cerr << InlineLinearDivergence2Env::name
                << ": std::map call failed: " << e << std::endl;
    QDP_abort(1);
  }
  const multi1d<LatticeColorMatrix> &u =
      TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(
          params.named_obj.gauge_id);

  push(xml_out, "propagator");
  write(xml_out, "update_no", update_no);

  proginfo(xml_out); // Print out basic program info

  // Write out the input
  params.writeXML(xml_out, "Input");

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);

  // These pesky variables are needed in the quarkprop call - only chiral dudes
  // need this stuff, but it must be there for the bleeping virtual function
  // to live at the base class

  const InlineLinearDivergence2Params::NamedObject_t::Props_t named_obj =
      params.named_obj.sink_pairs[0];
  AllSinkProps_t all_sinks;
  readAllSinks(all_sinks, named_obj);
  multi1d<int> t_srce =
      all_sinks.sink_prop_1.prop_header.source_header.getTSrce();
  int j_decay = all_sinks.sink_prop_1.prop_header.source_header.j_decay;
  int t0 = all_sinks.sink_prop_1.prop_header.source_header.t_source;
  int it0 = t0;

  SftMom phases_nomom(0, true, j_decay);
  Set timeslice = phases_nomom.getSet();

  int nt = Layout::lattSize()[3];
  int nx = Layout::lattSize()[2];

  int op_size = 4 * nt / 2 + 2;
  multi1d<LatticeComplex> quark_loop(4 * op_size);

  LatticeComplex c;
  LatticeReal rnd1, theta;
  random(rnd1);
  Real twopiN = Chroma::twopi / 3; // Z_4 grid for loop;
  theta = twopiN * floor(3 * rnd1);
  c = cmplx(cos(theta), sin(theta));

  srand(time(0));

  general_data_base result;
  result.add_dimension(dim_conf, 1, &params.cfg_serial);
  result.add_dimension(dim_displacement, op_size);
  result.dim[1].indices[0] = 0;
  {
    int count = 1;
    for (int idr = 0; idr < op_size; idr++) {
      result.dim[1].indices[count++] = idr;
    }
  }
  result.add_dimension(dim_operator, 4);
  result.add_dimension(dim_t, nt);
  result.add_dimension(dim_complex, 2);
  if (Layout::primaryNode())
    result.initialize();

  QDPIO::cout << InlineLinearDivergence2Env::name << ": initialized "
              << std::endl;

  for (int is = 0; is < 4 * op_size; is++)
    quark_loop[is] = zero;

  // References for use later
  const LatticePropagator &sink_prop_1 =
      TheNamedObjMap::Instance().getData<LatticePropagator>(
          all_sinks.sink_prop_1.sink_prop_1_id);
  const LatticePropagator &sink_prop_2 =
      TheNamedObjMap::Instance().getData<LatticePropagator>(
          all_sinks.sink_prop_2.sink_prop_1_id);

  // handle the loop contraction.
  LatticePropagator F, F_mu, F_plus, F_minus, F_tmp;
  LatticeComplex tmp, tmp2, lcone;
  lcone = cmplx(Real(1.0), Real(0.0));
  // Construct time direction links
  LatticeColorMatrix latticeones, tmpu1, tmpu2;
  latticeones = cmplx(Real(1.0), Real(0.0));

  F = sink_prop_1;
  QDPIO::cout << "calculate type I operator" << std::endl;
  //   type I
  int ig = 1;
  int icount;
  for (int imu = 0; imu < 4; imu++) {
    F_mu = sink_prop_1;
    icount = 0;
    for (int iop = 0; iop < nt / 2; iop++) {
      quark_loop[imu + icount * 4] = LatticeComplex(trace(F_mu * Gamma(ig)));
      //          F_tmp = F_mu * shift(adj(u[imu]),FORWARD,imu);
      F_tmp = adj(u[imu]) * shift(F_mu, FORWARD, imu);
      F_mu = F_tmp;
      icount++;
    }
    ig = ig * 2;
  }

  ig = 1;
  for (int imu = 0; imu < 4; imu++) {
    F_mu = sink_prop_1;
    icount = nt / 2;
    for (int iop = 0; iop < nt / 2; iop++) {
      quark_loop[imu + icount * 4] = LatticeComplex(trace(F_mu * Gamma(ig)));
      //        F_tmp = F_mu * shift(u[imu], FORWARD, imu);
      F_tmp = u[imu] * shift(F_mu, FORWARD, imu);
      F_mu = F_tmp;
      icount++;
    }
    ig = ig * 2;
  }

  QDPIO::cout << "calculate type II operator" << std::endl;
  // type II
  ig = 1;
  for (int imu = 0; imu < 4; imu++) {
    tmpu2 = where((((Layout::latticeCoordinate(3)) - it0) == 0), adj(u[imu]),
                  u[imu]);
    F_mu = sink_prop_1;
    icount = 2 * nt / 2;
    for (int iop = 0; iop < nx / 2; iop++) {
      if (icount == 2 * nt / 2) {
        quark_loop[imu + icount * 4] =
            LatticeComplex(trace(F_mu * Gamma(ig) * adj(F_mu) * Gamma(ig)));
      } else {
        quark_loop[imu + icount * 4] =
            LatticeComplex(trace(F_mu * Gamma(ig) * adj(F) * Gamma(ig)));
      }
      icount++;
      F_tmp = tmpu2 * shift(F_mu, FORWARD, imu);
      F_mu = F_tmp;
    }
    ig = ig * 2;
  }

  // calculate meson spectrum
  QDPIO::cout << "calculate meson spectrum" << std::endl;
  icount = 3 * nt / 2;
  ig = 1;
  for (int imu = 0; imu < 4; imu++) {
    quark_loop[imu + icount * 4] = LatticeComplex(
        trace(Gamma(15) * adj(F) * Gamma(15) * Gamma(ig) * F * Gamma(ig)));
    ig = ig * 2;
  }
  icount++;

  ig = 1;
  for (int imu = 0; imu < 4; imu++) {
    F_minus = shift(adj(u[imu]) * F, BACKWARD, imu);
    F_plus = u[imu] * shift(F, FORWARD, imu);
    quark_loop[imu - 1 + icount * 4] =
        0.5 * LatticeComplex(trace(Gamma(15) * adj(F) * Gamma(15) *
                                   (F_plus - F_minus) * Gamma(ig)));
    ig = ig * 2;
  }
  icount++;

  QDPIO::cout << InlineLinearDivergence2Env::name << "quark loop done"
              << std::endl;

  for (int is = 0; is < op_size; is++) {
    for (int ig = 0; ig < 4; ig++) {
      multi1d<DComplex> hsum = sumMulti(quark_loop[is * 4 + ig], timeslice);
      if (Layout::primaryNode())
        for (int t = 0; t < nt;
             t++) // addtional minus sign for the disconnect quark loop
        {
          result.data[(is * 4 + ig) * 2 * nt + 2 * t] -=
              toDouble(hsum[t].elem().elem().elem().real());
          result.data[(is * 4 + ig) * 2 * nt + 2 * t + 1] -=
              toDouble(hsum[t].elem().elem().elem().imag());
        }
    }
  }
  if (Layout::primaryNode()) {
    sprintf(result.name, params.iog_file.c_str());
    result.save();
  }

  snoop.stop();
  QDPIO::cout << InlineLinearDivergence2Env::name
              << ": total time = " << snoop.getTimeInSeconds() << " secs"
              << std::endl;

  QDPIO::cout << InlineLinearDivergence2Env::name << ": ran successfully"
              << std::endl;

  END_CODE();
}

} // namespace Chroma
