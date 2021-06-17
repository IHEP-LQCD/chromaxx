/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */

#include "inline_building_blocks_w_ihep.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "BuildingBlocks_w_ihep.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"

namespace Chroma {
namespace InlineBuildingBlocksIHEPEnv {
namespace {
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineBuildingBlocksIHEP(InlineBuildingBlocksIHEPParams(xml_in, path));
}

//! Local registration flag
bool registered = false;
}

const std::string name = "BUILDING_BLOCKS_IHEP";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= CreateFermStateEnv::registerAll();
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    registered = true;
  }
  return success;
}
}

//! Param input
void read(XMLReader &xml, const std::string &path,
          InlineBuildingBlocksIHEPParams::Param_t &input) {
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);
  input.use_sink_offset = false;
  input.canonical = false;
  input.time_reverse = false;
  input.translate = false;

  input.cfs = CreateFermStateEnv::nullXMLGroup();

  switch (version) {
  case 1:
    break;

  case 2:
    read(paramtop, "use_sink_offset", input.use_sink_offset);
    break;

  case 3:
    read(paramtop, "use_sink_offset", input.use_sink_offset);
    read(paramtop, "canonical", input.canonical);
    break;

  case 4:
    read(paramtop, "use_sink_offset", input.use_sink_offset);
    read(paramtop, "canonical", input.canonical);

    if (paramtop.count("FermState") != 0)
      input.cfs = readXMLGroup(paramtop, "FermState", "Name");
    break;

  case 5:
    read(paramtop, "use_sink_offset", input.use_sink_offset);
    read(paramtop, "canonical", input.canonical);
    read(paramtop, "time_reverse", input.time_reverse);
    read(paramtop, "translate", input.translate);

    if (paramtop.count("FermState") != 0)
      input.cfs = readXMLGroup(paramtop, "FermState", "Name");
    break;

  default:
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name << ": input parameter version "
                << version << " unsupported." << std::endl;
    QDP_abort(1);
  }

  read(paramtop, "links_max", input.links_max);
  read(paramtop, "mom2_max", input.mom2_max);
}

//! Param write
void write(XMLWriter &xml, const std::string &path,
           const InlineBuildingBlocksIHEPParams::Param_t &input) {
  push(xml, path);

  int version = 5;
  write(xml, "version", version);
  write(xml, "links_max", input.links_max);
  write(xml, "mom2_max", input.mom2_max);
  write(xml, "canonical", input.canonical);
  write(xml, "time_reverse", input.time_reverse);
  write(xml, "translate", input.translate);
  xml << input.cfs.xml;

  pop(xml);
}

//! Propagator input
void read(XMLReader &xml, const std::string &path,
          InlineBuildingBlocksIHEPParams::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "BkwdPropId", input.BkwdPropId);
  read(inputtop, "BkwdPropG5Format", input.BkwdPropG5Format);
  read(inputtop, "GammaInsertion", input.GammaInsertion);
  read(inputtop, "Flavor", input.Flavor);
}

//! Propagator output
void write(XMLWriter &xml, const std::string &path,
           const InlineBuildingBlocksIHEPParams::NamedObject_t &input) {
  push(xml, path);

  write(xml, "BkwdPropId", input.BkwdPropId);
  write(xml, "BkwdPropG5Format", input.BkwdPropG5Format);
  write(xml, "GammaInsertion", input.GammaInsertion);
  write(xml, "Flavor", input.Flavor);

  pop(xml);
}

//! BB parameters
void read(XMLReader &xml, const std::string &path,
          InlineBuildingBlocksIHEPParams::BB_out_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "GaugeId", input.GaugeId);
  read(inputtop, "FrwdPropId", input.FrwdPropId);
  read(inputtop, "BkwdProps", input.BkwdProps);
}

//! BB parameters
void write(XMLWriter &xml, const std::string &path,
           const InlineBuildingBlocksIHEPParams::BB_out_t &input) {
  push(xml, path);

  write(xml, "GaugeId", input.GaugeId);
  write(xml, "FrwdPropId", input.FrwdPropId);
  write(xml, "BkwdProps", input.BkwdProps);

  pop(xml);
}

// Param stuff
InlineBuildingBlocksIHEPParams::InlineBuildingBlocksIHEPParams() { frequency = 0; }

InlineBuildingBlocksIHEPParams::InlineBuildingBlocksIHEPParams(
    XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Read program parameters
    read(paramtop, "Param", param);

    // Read in the building block setup
    read(paramtop, "BuildingBlocks", bb);

    // Possible alternate XML file pattern
    if (paramtop.count("xml_file") != 0) {
      read(paramtop, "xml_file", xml_file);
    }
  }
  catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

void InlineBuildingBlocksIHEPParams::write(XMLWriter &xml_out,
                                       const std::string &path) {
  push(xml_out, path);

  Chroma::write(xml_out, "Param", param);
  Chroma::write(xml_out, "BuildingBlocks", bb);
  QDP::write(xml_out, "xml_file", xml_file);

  pop(xml_out);
}

//###################################################################################//
// Accept All Link Patterns //
//###################################################################################//

void AllLinkPatternsIHEP(bool &DoThisPattern, bool &DoFurtherPatterns,
                     multi1d<unsigned short int> &LinkPattern) {
  DoThisPattern = true;
  DoFurtherPatterns = true;

  return;
}

// Function call
void InlineBuildingBlocksIHEP::operator()(unsigned long update_no,
                                      XMLWriter &xml_out) {
  // If xml file not empty, then use alternate
  if (params.xml_file != "") {
    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

    push(xml_out, "ExampleBuildingBlocks");
    write(xml_out, "update_no", update_no);
    write(xml_out, "xml_file", xml_file);
    pop(xml_out);

    XMLFileWriter xml(xml_file);
    func(update_no, xml);
  } else {
    func(update_no, xml_out);
  }
}

// Function call
void InlineBuildingBlocksIHEP::func(unsigned long update_no, XMLWriter &XmlOut) {
  START_CODE();

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  push(XmlOut, "BuildingBlocks");
  write(XmlOut, "update_no", update_no);

  QDPIO::cout << " BuildingBlocks" << std::endl;
  QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
  for (int i = 1; i < Nd; ++i) {
    QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
  }
  QDPIO::cout << std::endl;

  //#################################################################################//
  // XML output
  //#################################################################################//

  proginfo(XmlOut); // Print out basic program info

  push(XmlOut, "Output_version");
  write(XmlOut, "out_version", 2);
  pop(XmlOut);

  //###############################################################################//
  // Read Gauge Field //
  //###############################################################################//

  // Grab the gauge field
  multi1d<LatticeColorMatrix> U;
  XMLBufferWriter gauge_xml;

  try {
    U = TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(
        params.bb.GaugeId);
    TheNamedObjMap::Instance().get(params.bb.GaugeId).getRecordXML(gauge_xml);

    // Set the construct state and modify the fields
    {
      QDPIO::cout << "cfs=XX" << params.param.cfs.xml << "XX" << std::endl;
      std::istringstream xml_s(params.param.cfs.xml);
      XMLReader fermtop(xml_s);

      Handle<CreateFermState<LatticeFermion, multi1d<LatticeColorMatrix>,
                             multi1d<LatticeColorMatrix> > >
      cfs(TheCreateFermStateFactory::Instance().createObject(
          params.param.cfs.id, fermtop, params.param.cfs.path));

      Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>,
                       multi1d<LatticeColorMatrix> > > state((*cfs)(U));

      // Pull the u fields back out from the state since they might have been
      // munged with fermBC's
      U = state->getLinks();
    }
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                << ": std::map call failed: " << e << std::endl;
    QDP_abort(1);
  }
  catch (...) {
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                << ": caught generic exception " << std::endl;
    QDP_abort(1);
  }

  // Write out the input
  params.write(XmlOut, "Input");

  // Write out the config info
  write(XmlOut, "Config_info", gauge_xml);

  // check that the gauge field seems normal
  Double ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace;
  MesPlq(U, ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace);

  push(XmlOut, "Observables");
  write(XmlOut, "ave_plaq", ave_plaq);
  write(XmlOut, "ave_spacelike_plaq", ave_spacelike_plaq);
  write(XmlOut, "ave_timelike_plaq", ave_timelike_plaq);
  write(XmlOut, "ave_link_trace", ave_link_trace);
  pop(XmlOut);

  //#################################################################################//
  // Read Forward Propagator //
  //#################################################################################//

  SftMom phases_nomom(0, true,
                      Nd - 1); // used to check props. Fix to Nd-1 direction.

  LatticePropagator F;
  ChromaProp_t prop_header;
  PropSourceConst_t source_header;
  QDPIO::cout << "Attempt to parse forward propagator" << std::endl;

  try {
    // Snarf a copy
    F = TheNamedObjMap::Instance().getData<LatticePropagator>(
        params.bb.FrwdPropId);

    // Snarf the frwd prop info. This is will throw if the frwd prop id is not
    // there
    XMLReader FrwdPropXML, FrwdPropRecordXML;
    TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getFileXML(
        FrwdPropXML);
    TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getRecordXML(
        FrwdPropRecordXML);

    // Try to invert this record XML into a ChromaProp struct
    {
      read(FrwdPropRecordXML, "/Propagator/ForwardProp", prop_header);
      read(FrwdPropRecordXML, "/Propagator/PropSource", source_header);
    }

    // Sanity check - write out the norm2 of the forward prop in the j_decay
    // direction
    // Use this for any possible verification
    {
      multi1d<Double> FrwdPropCheck =
          sumMulti(localNorm2(F), phases_nomom.getSet());

      // Write out the forward propagator header
      push(XmlOut, "ForwardProp");
      write(XmlOut, "FrwdPropId", params.bb.FrwdPropId);
      write(XmlOut, "FrwdPropXML", FrwdPropXML);
      write(XmlOut, "FrwdPropRecordXML", FrwdPropRecordXML);
      write(XmlOut, "FrwdPropCheck", FrwdPropCheck);
      pop(XmlOut);
    }
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                << ": forward prop: error message: " << e << std::endl;
    QDP_abort(1);
  }

  QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

  // Derived from input prop
  int j_decay = source_header.j_decay;
  multi1d<int> t_srce = source_header.getTSrce();

  //#################################################################################//
  // Read Backward (or Sequential) Propagators //
  //#################################################################################//

  StopWatch swatch;
  const int NF = params.bb.BkwdProps.size();

  push(XmlOut, "SequentialSource");

  for (int loop = 0; loop < NF; ++loop) {
    push(XmlOut, "elem");
    write(XmlOut, "loop_ctr", loop);

    QDPIO::cout << "Loop = " << loop << std::endl;

    multi1d<LatticePropagator> B(1);
    SeqSource_t seqsource_header;
    QDPIO::cout << "Attempt to parse backward propagator" << std::endl;
    try {
      // Extract the backward prop
      B[0] = TheNamedObjMap::Instance().getData<LatticePropagator>(
          params.bb.BkwdProps[loop].BkwdPropId);

      // Snarf the bkwd prop info. This is will throw if the bkwd prop id is not
      // there
      XMLReader BkwdPropXML, BkwdPropRecordXML;
      TheNamedObjMap::Instance()
          .get(params.bb.BkwdProps[loop].BkwdPropId)
          .getFileXML(BkwdPropXML);
      TheNamedObjMap::Instance()
          .get(params.bb.BkwdProps[loop].BkwdPropId)
          .getRecordXML(BkwdPropRecordXML);

      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      // NEED SECURITY HERE - need a way to cross check props. Use the ID.
      {
        read(BkwdPropRecordXML, "/SequentialProp/SeqSource", seqsource_header);
      }

      // Sanity check - write out the norm2 of the forward prop in the j_decay
      // direction
      // Use this for any possible verification
      {
        multi1d<Double> BkwdPropCheck =
            sumMulti(localNorm2(B[0]), phases_nomom.getSet());

        // Write out the forward propagator header
        push(XmlOut, "BackwardProp");
        write(XmlOut, "BkwdPropId", params.bb.BkwdProps[loop].BkwdPropId);
        write(XmlOut, "BkwdPropG5Format",
              params.bb.BkwdProps[loop].BkwdPropG5Format);
        write(XmlOut, "SequentialSourceType", seqsource_header.seqsrc.id);
        write(XmlOut, "BkwdPropXML", BkwdPropXML);
        write(XmlOut, "BkwdPropRecordXML", BkwdPropRecordXML);
        write(XmlOut, "BkwdPropCheck", BkwdPropCheck);
        pop(XmlOut);
      }
    }
    catch (std::bad_cast) {
      QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                  << ": forward prop: caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string &e) {
      QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name
                  << ": forward prop: error message: " << e << std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Backward propagator successfully parse" << std::endl;

    //#################################################################################//
    // Additional Gamma Matrix Insertions //
    //#################################################################################//

    multi1d<int> GammaInsertions(1);
    GammaInsertions[0] = params.bb.BkwdProps[loop].GammaInsertion;

    if (GammaInsertions[0] < 0 || GammaInsertions[0] >= Ns * Ns) {
      QDPIO::cerr << "InlineBuildingBlocks: Gamma insertion out of bounds: "
                  << GammaInsertions[0] << std::endl;
      QDP_abort(1);
    }

    //#################################################################################//
    // Some output //
    //#################################################################################//

    QDPIO::cout << "Seqsource name  = " << seqsource_header.seqsrc.id
                << std::endl;
    QDPIO::cout << "Gamma insertion = " << params.bb.BkwdProps[loop]
                                               .GammaInsertion << std::endl;
    QDPIO::cout << "Flavor          = " << params.bb.BkwdProps[loop].Flavor
                << std::endl;

    write(XmlOut, "seq_src", seqsource_header.seqsrc.id);
    write(XmlOut, "gamma_insertion", GammaInsertions[0]);
    write(XmlOut, "flavor", params.bb.BkwdProps[loop].Flavor);
    write(XmlOut, "t_source", source_header.t_source);
    write(XmlOut, "t_sink", seqsource_header.t_sink);
    write(XmlOut, "sink_mom", seqsource_header.sink_mom);

    //#################################################################################//
    // Convert Backward Propagator Format //
    //#################################################################################//

    if (params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B") {
      LatticePropagator Bu = B[0];
      B[0] = Gamma(15) * Bu;
    } else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "B_G5") {
      LatticePropagator Bu = B[0];
      B[0] = Bu * Gamma(15);
    } else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B_G5") {
      LatticePropagator Bu = B[0];
      B[0] = Gamma(15) * Bu * Gamma(15);
    }

    //#################################################################################//
    // Set Momenta //
    //#################################################################################//

    multi1d<int> SnkMom(Nd - 1);
    if (params.param.use_sink_offset)
      SnkMom = seqsource_header.sink_mom;
    else
      SnkMom = 0;

    SftMom Phases(params.param.mom2_max, t_srce, SnkMom, false, j_decay);
    SftMom PhasesCanonical(params.param.mom2_max, t_srce, SnkMom,
                           params.param.canonical, j_decay);

    //#################################################################################//
    // Construct File Names //
    //#################################################################################//

    int NumO = PhasesCanonical.numMom();

    //#################################################################################//
    // Flavor-ology //
    //#################################################################################//

    multi1d<int> Flavors(1);

    //
    // Dru puts a Flavor into the BB
    // Telephone book std::map of flavor name to a number
    //
    if (params.bb.BkwdProps[loop].Flavor == "U")
      Flavors[0] = 0;
    else if (params.bb.BkwdProps[loop].Flavor == "D")
      Flavors[0] = 1;
    else if (params.bb.BkwdProps[loop].Flavor == "S")
      Flavors[0] = 2;
    else if (params.bb.BkwdProps[loop].Flavor == "C")
      Flavors[0] = 3;
    else if (params.bb.BkwdProps[loop].Flavor == "B")
      Flavors[0] = 4;
    else if (params.bb.BkwdProps[loop].Flavor == "T")
      Flavors[0] = 5;
    else {
      QDPIO::cerr << InlineBuildingBlocksIHEPEnv::name << ": invalid flavor tag = "
                  << params.bb.BkwdProps[loop].Flavor
                  << ", should be one of U,D,S,C,T,B" << std::endl;
      QDP_abort(1);
    }

    //#################################################################################//
    // Loop-back tests //
    //#################################################################################//
    {
      //! Test a meson sequential source.
      /*!
       *  For the case of a meson, we have evaluated as the sequential source
       *
       *  H(y, 0; tx, p) = \sum exp{ip.x} U(y,x) \gamma_5\Gamma_f^\dag\gamma_5
       *D(x,0)
       *
       *  H^\dag(y, 0; tx, p) = \sum_x exp{-ip.x} \gamma_5 D(0,x) \Gamma_f
       *U(x,y) \gamma_5
       *
       *  Thus we can see that
       *
       *  Tr[ \gamma_5 H^\dag(0,0; tx, p)\gamma_5 \Gamma_i] =
       *                         \sum_x exp{-ip.x} Tr[ D(0,x)\Gamma_f U(x,0)
       *\Gamma_i ]
       *
       *  which is the desired meson correlator at momentum p and timslice tx
       */

      // assumes any Gamma5 matrices have already been absorbed
      push(XmlOut, "LoopBackTest");
      write(XmlOut, "seq_src", seqsource_header.seqsrc.id);
      write(XmlOut, "gamma_insertion", GammaInsertions[0]);
      write(XmlOut, "flavor", params.bb.BkwdProps[loop].Flavor);
      write(XmlOut, "t_srce", t_srce);
      write(XmlOut, "t_source", source_header.t_source);
      write(XmlOut, "t_sink", seqsource_header.t_sink);
      write(XmlOut, "sink_mom", seqsource_header.sink_mom);
      LatticeComplex tr = trace(adj(B[0]) * Gamma(GammaInsertions[0]));
      Complex seq_src_value = peekSite(tr, t_srce);
      write(XmlOut, "source_value", seq_src_value);
      pop(XmlOut);
    }

    //#################################################################################//
    // Diagnostic tests //
    //#################################################################################//
    {
      // assumes any Gamma5 matrices have already been absorbed
      int GammaInsertion = GammaInsertions[0];

      push(XmlOut, "DiagnosticTest");
      write(XmlOut, "GammaInsertion", GammaInsertion);

      {
        LatticePropagator GFG = Gamma(0) * F * Gamma(GammaInsertion);
        LatticeComplex tr = localInnerProduct(B[0], GFG);
        multi1d<DComplex> pr = sumMulti(tr, Phases.getSet());
        write(XmlOut, "formFactor_G0", pr);
      }
      {
        LatticePropagator GFG = Gamma(8) * F * Gamma(GammaInsertion);
        LatticeComplex tr = localInnerProduct(B[0], GFG);
        multi1d<DComplex> pr = sumMulti(tr, Phases.getSet());
        write(XmlOut, "formFactor_G8", pr);
      }
      {
        LatticePropagator GFG = Gamma(Ns * Ns - 1) * F * Gamma(GammaInsertion);
        LatticeComplex tr = localInnerProduct(B[0], GFG);
        multi1d<DComplex> pr = sumMulti(tr, Phases.getSet());
        write(XmlOut, "formFactor_G15", pr);
      }

      pop(XmlOut);
    }

    //#################################################################################//
    // Construct Building Blocks //
    //#################################################################################//

    swatch.reset();
    QDPIO::cout << "calculating building blocks" << std::endl;

    const signed short int T1 = 0;
    const signed short int T2 = QDP::Layout::lattSize()[j_decay] - 1;
    const signed short int DecayDir = j_decay;
    const signed short int Tsrc = source_header.t_source;
    const signed short int Tsnk = seqsource_header.t_sink;

    swatch.start();
    BuildingBlocks(B, F, U, GammaInsertions, Flavors, params.param.links_max,
                   AllLinkPatternsIHEP, Phases, PhasesCanonical, XmlOut, T1, T2,
                   Tsrc, Tsnk, seqsource_header.seqsrc.id,
                   seqsource_header.sink_mom, DecayDir,
                   params.param.time_reverse, params.param.translate);
    swatch.stop();

    QDPIO::cout << "finished calculating building blocks for loop = " << loop
                << "  time= " << swatch.getTimeInSeconds() << " secs"
                << std::endl;

    pop(XmlOut); // elem
  }              // end loop over sequential sources

  pop(XmlOut); // SequentialSource

  pop(XmlOut); // ExampleBuildingBlocks

  snoop.stop();
  QDPIO::cout << InlineBuildingBlocksIHEPEnv::name
              << ": total time = " << snoop.getTimeInSeconds() << " secs"
              << std::endl;

  QDPIO::cout << InlineBuildingBlocksIHEPEnv::name << ": ran successfully"
              << std::endl;

  END_CODE();
}
}
