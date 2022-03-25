#include "register_new.h"
namespace Chroma {
void register_new(bool &foo) {
  // New measurement registrations go here

  foo &= InlineLinearDivergence2Env::registerAll();
  foo &= InlineTestsEnv::registerAll();
  foo &= InlineUKYReadNamedObjEnv::registerAll();
  foo &= InlineClusterDecEnv::registerAll();
  foo &= InlineLinearDivergenceEnv::registerAll();
  foo &= InlineUKYWriteNamedObjEnv::registerAll();
  foo &= AnisoSpectrumGaugeActEnv2::registerAll();
  foo &= InlineBuildingBlocksIHEPEnv::registerAll();
  foo &= InlinePropAndMatElemDistillationIHEPEnv::registerAll();
  foo &= InlineMesonMatElemColorVecIHEPEnv::registerAll();
  foo &= FixedGaugeBCEnv::registerAll();
}
}
