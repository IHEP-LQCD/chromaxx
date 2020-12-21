// $Id: inline_propagator_w.cc,v 3.11 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "inline_cluster_dec.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "util/gauge/shift2.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include <time.h>
#include "io_general_class.h"

namespace Chroma
{
  namespace InlineClusterDecEnv
  {
    namespace
    {
      AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                              const std::string &path)
      {
        return new InlineClusterDec(InlineClusterDecParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    } // namespace

    const std::string name = "ClusterDec";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
        success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }
  } // namespace InlineClusterDecEnv

  //! gauge input
  void read(XMLReader &xml, const std::string &path, InlineClusterDecParams::NamedObject_t &input)
  {
    XMLReader inputtop(xml, path);
    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! gauge output
  void write(XMLWriter &xml, const std::string &path, const InlineClusterDecParams::NamedObject_t &input)
  {
    push(xml, path);
    write(xml, "gauge_id", input.gauge_id);
    pop(xml);
  }

  // Param stuff
  InlineClusterDecParams::InlineClusterDecParams() { frequency = 0; }

  InlineClusterDecParams::InlineClusterDecParams(XMLReader &xml_in, const std::string &path)
  {
    try
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
        read(paramtop, "Frequency", frequency);
      else
        frequency = 1;

      read(paramtop, "cfg_serial", cfg_serial);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "bl_level", bl_level);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);

      read(paramtop, "iog_file", iog_file);
    }
    catch (const std::string &e)
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }

  void InlineClusterDecParams::writeXML(XMLWriter &xml_out, const std::string &path)
  {
    push(xml_out, path);
    write(xml_out, "NamedObject", named_obj);
    pop(xml_out);
  }

  multi2d<Double> shiftop(LatticeDouble op, int nrcut, int half_length)
  {
    LatticeDouble shift_tmp, shift_f, shift_b, shift_t, shift_op;
    multi2d<Double> glue(nrcut, half_length);
    shift_op = 0;
    shift_f = op;
    shift_b = op;
    shift_t = op;

    for (int it = 0; it < half_length; it++)
    {
      for (int ircut = 0; ircut < nrcut; ircut++)
      {
        for (int i = 0; i < 3; i++)
        {
          shift_tmp = shift_f;
          shift_f = shift(shift_tmp, FORWARD, i);
          shift_tmp = shift_b;
          shift_b = shift(shift_tmp, BACKWARD, i);
          shift_op = shift_op + shift_f + shift_b;
        }
        glue[ircut][it] = sum(shift_op * op);
      }
      shift_tmp = shift_t;
      shift_t = shift(shift_tmp, FORWARD, 3);
      shift_f = shift_t;
      shift_b = shift_t;
    }
    return glue;
  } // namespace Chroma

  // Function call
  void InlineClusterDec::operator()(unsigned long update_no,
                                    XMLWriter &xml_out)
  {
    func(update_no, xml_out);
  }

  void InlineClusterDec::func(unsigned long update_no,
                              XMLWriter &xml_out)
  {
    START_CODE();

    QDPIO::cout << InlineClusterDecEnv::name << ":glueball cluster decomposition calculation" << std::endl;

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineClusterDecEnv::name << ": caught dynamic cast error"
                  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string &e)
    {
      QDPIO::cerr << InlineClusterDecEnv::name << ": std::map call failed: " << e
                  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix> &u =
        TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);

    // Print out basic program info
    proginfo(xml_out);

    // Write out the input
    params.writeXML(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //copy from gluecor.cc
    // Length of lattice in decay direction
    int nx = Layout::lattSize()[2];
    int nt = Layout::lattSize()[3];
    int j_decay = params.j_decay;
    int bl_level = params.bl_level;
    SftMom phases(0, true, j_decay);
    int length = phases.numSubsets();

    int half_length = length / 2 + 1;
    int nrcut = nx / 2;

    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    multi1d<LatticeColorMatrix> cplaq_tmp(6);
    LatticeDouble rcut;
    LatticeDouble one;
    LatticeDouble op0, op1, op2, op3;
    one = 1.0;

    multi2d<Double> glue0(nrcut, half_length); //scalar
    multi2d<Double> glue1(nrcut, half_length);
    multi2d<Double> glue2(nrcut, half_length);
    multi2d<Double> glue3(nrcut, half_length);
    Double vac0;
    Double dummy;
    multi1d<multi1d<DComplex>> plaq(6);
    int mu;
    int nu;
    int plane;
    int t;
    int t0;
    int t1;
    int rcut_step = (3 * pow(nx / 2, 2)) / nrcut;

    QDPIO::cout << InlineClusterDecEnv::name << ": nrcut = " << nrcut << "; rcut_step = " << rcut_step << std::endl;
    if (Nd != 4)
      QDP_error_exit("Nd for glueball construction has to be 4 but: ", Nd);

    for (int i = 0; i < 6; ++i)
    {
      plaq[i].resize(length);
      plaq[i] = 0;
    }

    /* Initialize glueball operators to zero */
    /* Construct glueball operators */
    //    for (int ircut = 0; ircut < nrcut; ++ircut)
    //    {
    /*      
rcut = where((Layout::latticeCoordinate(0) - nx / 2) * (Layout::latticeCoordinate(0) - nx / 2) +
                           (Layout::latticeCoordinate(1) - nx / 2) * (Layout::latticeCoordinate(1) - nx / 2) +
                           (Layout::latticeCoordinate(2) - nx / 2) * (Layout::latticeCoordinate(2) - nx / 2) <=
                       ircut * rcut_step,
                   one, LatticeDouble(zero));
*/
    op0 = 0;
    op1 = 0;
    op2 = 0;
    op3 = 0;
    plane = 0;
    for (mu = 0; mu < Nd; ++mu)
      for (nu = mu + 1; nu < Nd; ++nu)
      {
        /* Construct (block) plaquettes */
        /* tmp_1(x) = u(x+mu*2^bl_level,nu) */
        //tmp_1 = shift2(u[nu], FORWARD, mu, bl_level);
        tmp_1 = shift(u[nu], FORWARD, mu);

        /* tmp_2(x) = u(x+nu*2^bl_level,mu) */
        //tmp_2 = shift2(u[mu], FORWARD, nu, bl_level);
        tmp_2 = shift(u[mu], FORWARD, nu);

        /* tmp_3 = tmp_1 * tmp_2_dagger */
        tmp_3 = tmp_1 * adj(tmp_2);

        /* tmp_1 = tmp_3 * u_dagger(x,nu) */
        tmp_1 = tmp_3 * adj(u[nu]);

        /* cplaq_tmp = Tr(u(x,mu) * tmp_1) */
        cplaq_tmp[plane] = (u[mu] * tmp_1 - adj(u[mu] * tmp_1)) / 2.0;

        plane = plane + 1;
      } //mu nu
    op0 = real(trace(cplaq_tmp[0] * cplaq_tmp[0] + cplaq_tmp[1] * cplaq_tmp[1] + cplaq_tmp[3] * cplaq_tmp[3]));
    op1 = real(trace(cplaq_tmp[2] * cplaq_tmp[2] + cplaq_tmp[4] * cplaq_tmp[4] + cplaq_tmp[5] * cplaq_tmp[5]));
    op2 = real(trace(cplaq_tmp[3] * cplaq_tmp[2] - cplaq_tmp[1] * cplaq_tmp[4] + cplaq_tmp[0] * cplaq_tmp[5]));
    //    } //ircut

    QDPIO::cout << InlineClusterDecEnv::name << ":op's calculation is done!  " << std::endl;

    std::string op_out_filename = "./op_B2_E2_BE_nx_" + std::to_string(nx) + "_nt_" + std::to_string(nt) +
                                  "_conf_" + std::to_string(params.cfg_serial) + ".dat";
    QDPIO::cout << InlineClusterDecEnv::name << ":save op's filename is  " << op_out_filename << std::endl;
    BinaryFileWriter op_out(op_out_filename);
    write(op_out, op0);
    write(op_out, op1);
    write(op_out, op2);
    op_out.close();

    //    /* Initialize glueball correlations to zero */
    //    vac0 = 0;
    //    glue0 = 0;
    //    glue1 = 0;
    //    glue2 = 0;
    //    glue3 = 0;
    //
    //    glue0 = shiftop(op0,nrcut,half_length);
    //    glue1 = shiftop(op1,nrcut,half_length);
    //    glue2 = shiftop(op2,nrcut,half_length);
    //    /* Normalize */
    //    vac0 /= Double(QDP::Layout::vol());
    //    dummy = Double(1) / Double(QDP::Layout::vol());
    //    for (int ir = 0; ir < nrcut; ir++)
    //      for (t = 0; t < (half_length); ++t)
    //      {
    //        glue0[ir][t] = glue0[ir][t] * dummy;
    //        glue1[ir][t] = glue1[ir][t] * dummy;
    //        glue2[ir][t] = glue2[ir][t] * dummy;
    //      }
    //
    //    general_data_base result;
    //    result.add_dimension(dim_conf, 1, &params.cfg_serial);
    //    result.add_dimension(dim_displacement, nrcut);
    //    result.add_dimension(dim_operator, 4);
    //    result.add_dimension(dim_t, half_length);
    //
    //    if (Layout::primaryNode())
    //    {
    //      result.initialize();
    //      for (int ir = 0; ir < nrcut; ir++)
    //        for (int t = 0; t < half_length; t++)
    //        {
    //          result.data[(ir * 4 + 0) * half_length + t] = glue0[ir][t].elem().elem().elem().elem();
    //          result.data[(ir * 4 + 1) * half_length + t] = glue1[ir][t].elem().elem().elem().elem();
    //          result.data[(ir * 4 + 2) * half_length + t] = glue2[ir][t].elem().elem().elem().elem();
    //          result.data[(ir * 4 + 3) * half_length + t] = glue3[ir][t].elem().elem().elem().elem();
    //        }
    //
    //      sprintf(result.name, params.iog_file.c_str());
    //      result.save();
    //    }
    //
    snoop.stop();
    QDPIO::cout << InlineClusterDecEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    END_CODE();
  } //func

} // namespace Chroma