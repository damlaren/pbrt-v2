#include "wet.h"
#include "core/montecarlo.h"

#include "core/paramset.h"
#include "materials/skin.h"

#include <iostream>

// Monte-Carlo integration of a BRDF over a hemisphere, as in DJ2006.
Spectrum Wet::integrate_BRDF(BSDF *bsdf, const Vector& wi,
			     int sqrtSamples, BxDFType brdfType) {
  int nSamples = sqrtSamples * sqrtSamples;
  
  Spectrum rho_dr = Spectrum(0.0f);
  for (int i = 0; i < nSamples; i++) {
    Vector w0;
    float pdf;
    Spectrum f = bsdf->Sample_f(-wi, &w0, BSDFSample(rng), &pdf, brdfType);
    f *= AbsCosTheta(wi);
    if (pdf != 0) {
      f /= pdf;
    }
    rho_dr += f;
  }
  rho_dr /= nSamples; //duh!
  
  // Take final result and clamp to range [0, 1]
  rho_dr = rho_dr.Clamp(0.0f, 1.0f);

  if (!rho_dr.IsBlack()) {
    // rho_dr.Print(stderr); std::cout << std::endl;
  }

  // HACK! Fetch oiliness from BSDF.
  for (int i = 0; i < bsdf->NumComponents(); i++) {
    BxDF* bxdf = bsdf->getComponent(i);
    SkinMicrofacet *smf = dynamic_cast<SkinMicrofacet*>(bxdf);
    if (smf) {
      oiliness = smf->oiliness;
    }
  }
  
  return rho_dr;
}
