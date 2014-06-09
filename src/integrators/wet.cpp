#include "wet.h"
#include "core/montecarlo.h"

// x: Surface point-- but unused?
// wo: Outgoing light dir
// wi: Incoming light dir
// F: Fresnel term
float Wet::frTS(const Point& x, const Vector& wo, const Vector& wi, const Normal& N, const float F) {
  //  return OILINESS * F(x, wo, wi) * D(x, wo, wi, ROUGHNESS) * G(x, wo, wi, n) / (4.0f * Dot(wo,n));
  Vector L = -wi; // vector to light source (i.e. towards source of incoming light)
  const Vector& E = wo; // vector to camera/eye (i.e. the outgoing light ray)
  Vector H = Normalize(L + E);
  return (OILINESS * F * D(N, H, ROUGHNESS) * G(L, H, E, N) / (4.0f * Dot(wo, N)));
}

Spectrum Wet::F(const Point& x, const Vector& wo, const Vector& wi) {
  //TODO? not sure about this one. It's the Fresnel term but I'm not sure how to compute.
  //  I'll try recycling Fdr...
  return 1.0f;
}

// Beckmann distribution: http://en.wikipedia.org/wiki/Cook%E2%80%93Torrance#Cook.E2.80.93Torrance_model
float Wet::D(const Normal& N, const Vector& H, const float m) {
  float alpha = acosf(Dot(N, H));
  float t = tanf(alpha);
  float c = cosf(alpha);
  float e = t/m;
  float csquared = c*c;
  return exp(-e * e) / (m*m * csquared*csquared);
}

// Geometry term: http://en.wikipedia.org/wiki/Cook%E2%80%93Torrance#Cook.E2.80.93Torrance_model
float Wet::G(const Vector& L, const Vector& H, const Vector& E, const Normal& N) {
  float g = min(1.0f, 2.0f * Dot(H, N) * Dot(E, N) / Dot(E, H));
  return min(g, 2.0f * Dot(H, N) * Dot(L, N) / Dot(E, H));
}

// Monte-Carlo integration of frTS over a hemisphere
float Wet::integrate_frTS(const Point& x, const Vector& wo, const Normal& N, const float F) {
  RNG rng;
  const int SQRT_SAMPLES = 10; // TODO: how many
  const int N_SAMPLES = SQRT_SAMPLES * SQRT_SAMPLES;
  float *s1 = ALLOCA(float, 2 * N_SAMPLES);
  StratifiedSample2D(s1, SQRT_SAMPLES, SQRT_SAMPLES, rng);
  float *s2 = ALLOCA(float, 2 * N_SAMPLES);
  StratifiedSample2D(s2, SQRT_SAMPLES, SQRT_SAMPLES, rng);

  float rho_dr = 0.0f;
  for (int i = 0; i < N_SAMPLES; i++) {
    Vector wi = UniformSampleHemisphere(s1[i], s2[i]);
    //TODO rotate to account for the normal...

    rho_dr += frTS(x, wo, wi, N, F);
  }
  return rho_dr / (UniformHemispherePdf() * N_SAMPLES);
}
