// $Id: inline_propagator_w.cc,v 3.11 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "inline_linear_divergence.h"
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

#include "io_general_class.h"
#include "util/ferm/transf.h"
#include <time.h>

namespace Chroma {
namespace InlineLinearDivergenceEnv {
namespace {
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineLinearDivergence(InlineLinearDivergenceParams(xml_in, path));
}

//! Local registration flag
bool registered = false;
} // namespace

const std::string name = "LinearDivergence";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= WilsonTypeFermActsEnv::registerAll();
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    registered = true;
  }
  return success;
}
} // namespace InlineLinearDivergenceEnv

//! Propagator input
void read(XMLReader &xml, const std::string &path,
          InlineLinearDivergenceParams::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "gauge_id", input.gauge_id);
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const InlineLinearDivergenceParams::NamedObject_t &input) {
  push(xml, path);

  write(xml, "gauge_id", input.gauge_id);

  pop(xml);
}

// Param stuff
InlineLinearDivergenceParams::InlineLinearDivergenceParams() { frequency = 0; }

InlineLinearDivergenceParams::InlineLinearDivergenceParams(
    XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Parameters for source construction
    read(paramtop, "Grid", grid);
    if (grid.size() != 3) {
      QDPIO::cerr << InlineLinearDivergenceEnv::name
                  << ": wrong size of grid array: expected length=" << 3
                  << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "link_lengths", link_lengths);
    if (link_lengths.size() != 3) {
      QDPIO::cerr << InlineLinearDivergenceEnv::name
                  << ": wrong size of link_length array: expected length=" << 3
                  << std::endl;
      QDP_abort(1);
    }

    if (paramtop.count("use_ckpoint") != 0) {
      read(paramtop, "use_ckpoint", use_ckpoint);
    } else {
      use_ckpoint = true;
    }
    read(paramtop, "nsets", nsets);
    read(paramtop, "HL_ratio", HL_ratio);
    if (HL_ratio <= 0) {
      QDPIO::cerr << InlineLinearDivergenceEnv::name
                  << ": set HL_ratio<=0 is pointless, should be larger than 0."
                  << std::endl;
      QDP_abort(1);
    }
    if (nsets % HL_ratio != 0) {
      QDPIO::cerr << InlineLinearDivergenceEnv::name << ": HL_ratio ("
                  << HL_ratio << ") should be a factor of nsets (" << nsets
                  << ")." << std::endl;
      QDP_abort(1);
    }
    read(paramtop, "Param", param);
    inv_param_h = readXMLGroup(paramtop, "InvParamH", "invType");
    inv_param_l = readXMLGroup(paramtop, "InvParamL", "invType");
    read(paramtop, "cfg_serial", cfg_serial);

    // Read in the output propagator/source configuration info
    read(paramtop, "NamedObject", named_obj);

    read(paramtop, "iog_file", iog_file);
    // Possible alternate XML file pattern
    if (paramtop.count("xml_file") != 0) {
      read(paramtop, "xml_file", xml_file);
    }
  } catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

void InlineLinearDivergenceParams::writeXML(XMLWriter &xml_out,
                                            const std::string &path) {
  push(xml_out, path);

  write(xml_out, "Param", param);
  write(xml_out, "NamedObject", named_obj);

  pop(xml_out);
}

int get_ckpoint_info2(char *name) {
  int info = -1;
  if (Layout::primaryNode()) {
    FILE *pfile = fopen(name, "r");
    if (pfile != NULL) {
      fclose(pfile);
      general_data_base tmp(name);
      tmp.load_type();
      info = tmp.dim[1].indices[0];
    }
  }
  QDPInternal::broadcast(info);
  return info;
}
// Function call
void InlineLinearDivergence::operator()(unsigned long update_no,
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
void InlineLinearDivergence::func(unsigned long update_no, XMLWriter &xml_out) {
  START_CODE();

  QDPIO::cout << InlineLinearDivergenceEnv::name << ": loop calculation"
              << std::endl;

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
    QDPIO::cerr << InlineLinearDivergenceEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  } catch (const std::string &e) {
    QDPIO::cerr << InlineLinearDivergenceEnv::name
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
  int t0;
  int j_decay;
  SftMom phases(0, true, Nd - 1);

  int op_size = 3 * Layout::lattSize()[2] + 2;
  multi1d<LatticeComplex> quark_loop(16 * op_size);

  LatticeComplex c;
  LatticeReal rnd1, theta;
  random(rnd1);
  Real twopiN = Chroma::twopi / 3; // Z_4 grid for loop;
  theta = twopiN * floor(3 * rnd1);
  c = cmplx(cos(theta), sin(theta));

  srand(time(0));
  multi1d<int> offset(4);
  for (int idr = 0; idr < 3; idr++)
    offset[idr] =
        (params.grid[idr] > 0) ? params.grid[idr] : (Layout::lattSize()[idr]);
  offset[3] = 128;

  int nt = Layout::lattSize()[3];
  int nx = Layout::lattSize()[2];

  general_data_base result;
  result.add_dimension(dim_conf, 1, &params.cfg_serial);
  result.add_dimension(dim_displacement, op_size);
  for (int is = 0; is < op_size; is++)
    result.dim[1].indices[is] = is;
  result.add_dimension(dim_operator, 16);
  result.add_dimension(dim_t, nt);
  result.add_dimension(dim_complex, 2);
  if (Layout::primaryNode())
    result.initialize();

  QDPIO::cout << InlineLinearDivergenceEnv::name << ": initialized "
              << std::endl;

  bool mg_setup_flag = false;

  // create the sources;
  multi1d<int> pos(3);
  int flag = 0;
  for (int i = 0; i < 3; i++) {
    pos[i] = rand() % params.grid[i];
  }
  result.clear();

  for (int is = 0; is < 16 * op_size; is++)
    quark_loop[is] = zero;

  //         for(int it0=0;it0<2;it0++)
  int it0 = 0;
  int ieo = 0;
  LatticeComplex lcone;
  LatticeComplex c1;
  lcone = cmplx(Real(1.0), Real(0.0));
  //  c1 = where(
  //        (((Layout::latticeCoordinate(3)) - it0) % offset[3] == 0),
  //        c, LatticeComplex(zero));
  c1 = where((((Layout::latticeCoordinate(3)) - it0) % offset[3] == 0), lcone,
             LatticeComplex(zero));

  LatticePropagator quark_prop_source;

  for (int color_source = 0; color_source < Nc; ++color_source)
    for (int spin_source = 0; spin_source < Ns; ++spin_source) {
      LatticeColorVector colorvec = zero;
      LatticeFermion chi = zero;
      pokeSpin(chi, pokeColor(colorvec, c1, color_source), spin_source);
      FermToProp(chi, quark_prop_source, color_source, spin_source);
    }

  // original propagator generation starts here

  // Sanity check - write out the norm2 of the source in the Nd-1 direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases

    multi1d<DComplex> source_corr = sumMulti(
        trace(quark_prop_source * adj(quark_prop_source)), phases.getSet());
    if (Layout::primaryNode()) {
      for (int it = 0; it < Layout::lattSize()[3]; it++)
#if defined(QDP_IS_QDPJIT)
        printf("srcNROM:%4d%13.5f\n", it,
               source_corr[it].elem().elem().elem().real().elem());
#else
        printf("srcNROM:%4d%13.5f\n", it,
               source_corr[it].elem().elem().elem().real());
#endif
      fflush(stdout);
    }
  }

  LatticePropagator quark_propagator, quark_propagator_H;
  int ncg_had = 0;

  //
  // Initialize fermion action
  //
  std::istringstream xml_s(params.param.fermact.xml);
  XMLReader fermacttop(xml_s);
  QDPIO::cout << "FermAct = " << params.param.fermact.id << std::endl;

  //
  // Try the factories
  //
  bool success = false;

  if (!success) {
    try {
      StopWatch swatch;
      swatch.reset();
      QDPIO::cout << "Try the various factories" << std::endl;

      // Typedefs to save typing
      typedef LatticeFermion T;
      typedef multi1d<LatticeColorMatrix> P;
      typedef multi1d<LatticeColorMatrix> Q;

      // Generic Wilson-Type stuff
      Handle<FermionAction<T, P, Q>> S_f(
          TheFermionActionFactory::Instance().createObject(
              params.param.fermact.id, fermacttop, params.param.fermact.path));

      Handle<FermState<T, P, Q>> state(S_f->createState(u));

      QDPIO::cout << "Suitable factory found: compute the quark prop"
                  << std::endl;
      swatch.start();

      GroupXML_t invParamR =
          (mg_setup_flag == false) ? params.param.invParam : params.inv_param_l;
      //           	if(params.HL_ratio==1)invParamR=params.inv_param_h;

      QDPIO::cout << "Calling quarkProp" << std::endl;
      S_f->quarkProp(quark_propagator, xml_out, quark_prop_source, t0, j_decay,
                     state, invParamR, params.param.quarkSpinType,
                     params.param.obsvP, ncg_had);
      swatch.stop();
      QDPIO::cout << "Propagator computed: time= " << swatch.getTimeInSeconds()
                  << " secs" << std::endl;

      success = true;
    } catch (const std::string &e) {
      QDPIO::cout << InlineLinearDivergenceEnv::name
                  << ": caught exception around quarkprop: " << e << std::endl;
    }
  }

  if (!success) {
    QDPIO::cerr << "Error: no fermact found" << std::endl;
    QDP_abort(1);
  }

  // handle the loop contraction.
  LatticePropagator F, F_mu, F_plus, F_minus, F_tmp;
  F = quark_propagator;
  LatticeComplex tmp, tmp2;
  // Construct time direction links
  LatticeColorMatrix tmpu1, tmpu2, tmpu3;
  LatticeColorMatrix tmput1;
  LatticeComplex c2 = c1;
  for (int it = 0; it < nt; it++) {
    c2 += shift(c1, FORWARD, 3);
  }
  LatticeComplex c2_dagger = adj(c2);

  QDPIO::cout << "calculate type I operator" << std::endl;
  //   type I
  int ig = 1;
  int icount;
  for (int imu = 0; imu < 4; imu++) {
    F_mu = quark_propagator * c2_dagger;
    icount = 0;
    for (int iop = 0; iop < nx; iop++) {
      quark_loop[imu + icount * 16] = LatticeComplex(trace(F_mu * Gamma(ig)));
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
    tmpu1 = where((((Layout::latticeCoordinate(3)) - it0) % offset[3] == 0),
                  u[imu], adj(u[imu]));
    F_mu = quark_propagator;
    icount = nx;
    for (int iop = 0; iop < nx; iop++) {
      if (icount == nx) {
        quark_loop[ig + icount * 16] = LatticeComplex(
            trace(F_mu * c2_dagger * Gamma(ig) * adj(F_mu) * Gamma(ig)));
      } else {
        quark_loop[ig + icount * 16] =
            LatticeComplex(trace(F_mu * Gamma(ig) * adj(F) * Gamma(ig)));
      }
      icount++;
      F_tmp = tmpu1 * shift(F_mu, FORWARD, imu);
      F_mu = F_tmp;
    }
    ig = ig * 2;
  }

  ig = 1;
  for (int imu = 0; imu < 4; imu++) {
    tmpu2 = where((((Layout::latticeCoordinate(3)) - it0) % offset[3] == 0),
                  adj(u[imu]), u[imu]);
    F_mu = quark_propagator;
    icount = 2 * nx;
    for (int iop = 0; iop < nx; iop++) {
      if (icount == 2 * nx) {
        quark_loop[ig + icount * 16] = LatticeComplex(
            trace(F_mu * c2_dagger * Gamma(ig) * adj(F_mu) * Gamma(ig)));
      } else {
        quark_loop[ig + icount * 16] =
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
  for (int ig = 0; ig < 16; ig++)
    quark_loop[ig + icount * 16] = LatticeComplex(
        trace(Gamma(15) * adj(F) * Gamma(15) * Gamma(ig) * F * Gamma(ig)));
  icount++;

  int gamma_index = 1;
  for (int mu = 0; mu < 4; mu++) {
    F_minus = shift(adj(u[mu]) * F, BACKWARD, mu);
    F_plus = u[mu] * shift(F, FORWARD, mu);
    quark_loop[mu - 1 + icount * 16] =
        0.5 * LatticeComplex(trace(Gamma(15) * adj(F) * Gamma(15) *
                                   (F_plus - F_minus) * Gamma(gamma_index)));
    gamma_index = gamma_index * 2;
  }
  icount++;

  QDPIO::cout << InlineLinearDivergenceEnv::name << "quark loop done"
              << std::endl;

  for (int is = 0; is < op_size; is++) {
    for (int ig = 0; ig < 16; ig++) {
      multi1d<DComplex> hsum =
          sumMulti(quark_loop[is * 16 + ig], phases.getSet());
      if (Layout::primaryNode())
        for (int t = 0; t < nt;
             t++) // addtional minus sign for the disconnect quark loop
        {
          result.data[(is * 16 + ig) * 2 * nt + 2 * t] -=
              toDouble(hsum[t].elem().elem().elem().real());
          result.data[(is * 16 + ig) * 2 * nt + 2 * t + 1] -=
              toDouble(hsum[t].elem().elem().elem().imag());
        }
    }
  }

  // merge the data
  if (Layout::primaryNode()) {
    sprintf(result.name, params.iog_file.c_str());
    result.save();
  }

  snoop.stop();
  QDPIO::cout << InlineLinearDivergenceEnv::name
              << ": total time = " << snoop.getTimeInSeconds() << " secs"
              << std::endl;

  QDPIO::cout << InlineLinearDivergenceEnv::name << ": ran successfully"
              << std::endl;

  END_CODE();
}

} // namespace Chroma
