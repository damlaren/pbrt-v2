#include "wet.h"
#include "core/montecarlo.h"

#include <iostream>

// Monte-Carlo integration of a BRDF over a hemisphere, as in DJ2006.
Spectrum Wet::integrate_BRDF(BSDF *bsdf, const Vector& wi,
			     int sqrtSamples, BxDFType brdfType) {
  int nSamples = sqrtSamples * sqrtSamples;
  
  /*
  float *s1 = ALLOCA(float, 2 * nSamples);
  StratifiedSample2D(s1, sqrtSamples, sqrtSamples, rng);
  float *s2 = ALLOCA(float, 2 * nSamples);
  StratifiedSample2D(s2, sqrtSamples, sqrtSamples, rng);

  // Do the MC sampling of BRDFs. wrong way?
  Spectrum rho_dr = Spectrum(0.0f);
  for (int i = 0; i < nSamples; i++) {
    Vector wo = UniformSampleHemisphere(s1[i], s2[i]);
    rho_dr += bsdf->f(wo, wi, brdfType) * CosTheta(wi);
  }

  rho_dr /= UniformHemispherePdf() * nSamples;
  */

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
  return rho_dr;
}
