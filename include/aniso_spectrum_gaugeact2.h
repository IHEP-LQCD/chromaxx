// -*- C++ -*-
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#ifndef __aniso_spectrum_gaugeact2_h__
#define __aniso_spectrum_gaugeact2_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "meas/glue/mesplq.h"
//#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"

namespace Chroma {

/*! @ingroup gaugeacts */
namespace AnisoSpectrumGaugeActEnv2 {
extern const std::string name;
bool registerAll();
}

//! Parameter structure
/*! @ingroup gaugeacts */
struct AnisoSpectrumGaugeAct2Params {
  // Base Constructor
  AnisoSpectrumGaugeAct2Params() {};

  // Read params from some root path
  AnisoSpectrumGaugeAct2Params(XMLReader &xml_in, const std::string &path);

  Real beta;
  Real u_s;
  Real u_t;
  Real omega;
  Real plaq_c_s;
  Real plaq_c_t;
  Real rect_c_s;
  Real rect_c_t_2;
  AnisoParam_t aniso;
};

/*! @ingroup gaugeacts */
void read(XMLReader &xml, const std::string &path,
          AnisoSpectrumGaugeAct2Params &param);

/*! @ingroup gaugeacts */
void write(XMLWriter &xml, const std::string &path,
           const AnisoSpectrumGaugeAct2Params &param);

//! MP gauge action
/*! \ingroup gaugeacts
   *
   *  Anisotropic gaugeact useful for spectrum from hep-lat/9911003
   *
   *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
   *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
   */

class AnisoSpectrumGaugeAct2 : public LinearGaugeAction {
public:
  // Typedefs to save typing
  typedef multi1d<LatticeColorMatrix> P;
  typedef multi1d<LatticeColorMatrix> Q;

  //! Read beta from a param struct
  AnisoSpectrumGaugeAct2(Handle<CreateGaugeState<P, Q> > cgs_,
                         const AnisoSpectrumGaugeAct2Params &p)
      : param(p) {
    init(cgs_);
  }

  //! Is anisotropy used?
  bool anisoP() const { return param.aniso.anisoP; }

  //! Anisotropy factor
  const Real anisoFactor() const { return param.aniso.xi_0; }

  //! Anisotropic direction
  int tDir() const { return param.aniso.t_dir; }

  //! Return the set on which the gauge action is defined
  /*! Defined on the even-off (red/black) set */
  const Set &getSet() const { return rb; }

  //! update coeff of plaq and rect plaq
  void update_coeff(multi1d<LatticeColorMatrix> &u) {
    Double w_plaq, s_plaq, t_plaq, link;
    multi2d<Double> plane_plaq;
    MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
    Real u_s = pow(s_plaq,0.25);
    Real u_s_2 = u_s * u_s;
    Real u_s_4 = u_s_2 * u_s_2;
    Real u_s_6 = u_s_4 * u_s_2;
    Real u_t = pow(t_plaq / u_s_2 ,0.5);
    Real u_t_2 = u_t * u_t;
    QDPIO::cout << "s_plaq t_plaq u_s u_t are: " << s_plaq<< t_plaq<<u_s<< u_t << std::endl;
    param.plaq_c_s = param.beta * Real(5) / (Real(3) * u_s_4);
    if (param.aniso.anisoP)
      param.plaq_c_s /= param.aniso.xi_0;
    param.plaq_c_t = param.beta * Real(4) / (Real(3) * u_s_2 * u_t_2);
    if (param.aniso.anisoP)
      param.plaq_c_t *= param.aniso.xi_0;
    // Coefficients for the rectangle
    param.rect_c_s = -param.beta / (Real(12) * u_s_6);
    if (param.aniso.anisoP)
      param.rect_c_s /= param.aniso.xi_0;
    // Loops that are short in the time direction
    param.rect_c_t_2 =
        -param.beta / (Real(12) * u_s_4 * u_t_2) * param.aniso.xi_0;
    if (param.aniso.anisoP)
      param.rect_c_t_2 *= param.aniso.xi_0;

    QDPIO::cout << "plaq_c_s plaq_c_t rect_c_s rect_c_t_2 are: " << param.plaq_c_s<<param.plaq_c_t<<param.rect_c_s<<param.rect_c_t_2 <<std::endl;
  }

  //! Compute staple
  /*! Default version. Derived class should override this if needed. */
  /*

  */
  void staple(LatticeColorMatrix &result,
              const Handle<GaugeState<P, Q> > &state, int mu, int cb) const {
    START_CODE();
    const multi1d<LatticeColorMatrix> &u = state->getLinks();

    result = zero;
    Complex coeff_p;
    Complex coeff_r;
    bool ifspace;

    for (int nu = 0; nu < Nd; ++nu) {
      if (nu == mu)
        continue;
      ifspace = (nu != param.aniso.t_dir && mu != param.aniso.t_dir);
      if (!ifspace) {
        coeff_p = param.plaq_c_t;
        coeff_r = param.rect_c_t_2;
      } else {
        coeff_p = param.plaq_c_s;
        coeff_r = param.rect_c_s;
      }
      QDPIO::cout << "coeff_p and coeff_r are: "<< coeff_p << coeff_r << std::endl;

      LatticeColorMatrix tmp1 = shift(u[nu], FORWARD, mu);
      LatticeColorMatrix tmp2 = shift(u[mu], FORWARD, mu);
      LatticeColorMatrix tmp3 = shift(u[mu], FORWARD, nu);
      // 1mu 1nu up staple
      result[rb[cb]] += coeff_p * tmp1 * adj(tmp3) * adj(u[nu]);

      // 1mu 1nu down staple
      result[rb[cb]] += coeff_p * adj(shift(tmp1, BACKWARD, nu)) *
                        adj(shift(u[mu], BACKWARD, nu)) *
                        shift(u[nu], BACKWARD, nu);

      if (nu == param.aniso.t_dir || ifspace) {
        // 2 mu * 1nu up staple
        result[rb[cb]] += coeff_r * tmp2 * shift(tmp1, FORWARD, mu) *
                          adj(shift(tmp2, FORWARD, nu)) * adj(tmp3) *
                          adj(u[nu]);

        // 2 mu * 1nu down staple
        LatticeColorMatrix tmp4 = shift(tmp1, FORWARD, mu);

        result[rb[cb]] += coeff_r * tmp3 * adj(shift(tmp4, BACKWARD, nu)) *
                          adj(shift(tmp2, BACKWARD, nu)) *
                          adj(shift(u[mu], BACKWARD, nu)) *
                          shift(u[nu], BACKWARD, nu);
      }

      if (mu == param.aniso.t_dir || ifspace) {
        LatticeColorMatrix tmp5 = shift(u[nu], FORWARD, nu);
        // 2 nu * 1mu up staple
        result[rb[cb]] += coeff_r * tmp1 * shift(tmp1, FORWARD, nu) *
                          adj(shift(tmp3, FORWARD, nu)) * adj(tmp5) *
                          adj(u[nu]);

        LatticeColorMatrix tmp6 = shift(tmp1, BACKWARD, nu);
        LatticeColorMatrix tmp7 = shift(u[mu], BACKWARD, nu);
        LatticeColorMatrix tmp8 = shift(u[nu], BACKWARD, nu);
        // 2 nu * 1mu down staple
        result[rb[cb]] += coeff_r * adj(tmp6) * adj(shift(tmp6, BACKWARD, nu)) *
                          adj(shift(tmp7, BACKWARD, nu)) *
                          shift(tmp8, BACKWARD, nu) *
                          shift(u[nu], BACKWARD, nu);
      }
    }

    // NOTE: a heatbath code should be responsible for resetting links on
    // a boundary. The staple is not really the correct place.

    END_CODE();
  }

  //! Compute dS/dU
  void deriv(multi1d<LatticeColorMatrix> &result,
             const Handle<GaugeState<P, Q> > &state) const {
    plaq->deriv(result, state);
    // plaq_plus_two_plaq->deriv(result,state);

    multi1d<LatticeColorMatrix> tmp;
    rect->deriv(tmp, state);
    result += tmp;
  }

  //! Compute the actions
  Double S(const Handle<GaugeState<P, Q> > &state) const {
    return plaq->S(state) + rect->S(state);
    // return plaq_plus_two_plaq->S(state) + rect->S(state);
  }

  //! Produce a gauge create state object
  const CreateGaugeState<P, Q> &getCreateState() const {
    return plaq->getCreateState();
  }

  //! Destructor is automatic
  ~AnisoSpectrumGaugeAct2() {}

  // Accessors -- non mutable members.
  const Real getBeta(void) const { return param.beta; }

  const Real getUS(void) const { return param.u_s; }

  const Real getUT(void) const { return param.u_t; }

  const Real getOmega(void) const { return param.omega; }

protected:
  //! Private initializer
  void init(Handle<CreateGaugeState<P, Q> > cgs);

  //! Hide assignment
  void operator=(const AnisoSpectrumGaugeAct2 &a) {}

private:
  AnisoSpectrumGaugeAct2Params param; /*!< The couplings and anisotropy*/
  //    Handle<PlaqGaugeAct>           plaq;     /*!< Hold a plaquette gaugeact
  // */
  Handle<RectGaugeAct> rect; /*!< Hold a rectangle gaugeact */
  Handle<PlaqGaugeAct>
  plaq; /*!< Hold spatial plaquettes separated in time type gaugeact */
};
}

#endif
