/*
  An attempt to simulate the appearance of oil and the roughness of human skin, based on:
  * http://graphics.ucsd.edu/~henrik/papers/skin_bssrdf/skin_bssrdf.pdf [DJ2006]
  * http://graphics.ucsd.edu/papers/layered/layered.pdf [DJ2005]
  * "A Reflectance Model for Computer Graphics" [CT82]
 */

#include "integrator.h"
#include "pbrt.h"

class Wet {
 public:

  static const float OILINESS = 1.0;   // rho_s, in [0, 1]
  static const float ROUGHNESS = 0.35; // skin roughness

  // TODO: floats or Spectrums?
  // DJ2006 eq. 1, except that dot(wi, n) on bottom is left out.
  float frTS(const Point& x, const Vector& wo, const Vector& wi, const Normal& N, const float F);

  // Fresnel reflectance [TODO?]
  Spectrum F(const Point& x, const Vector& wo, const Vector& wi);

  // Beckmann microfacet Distribution [CT82]
  float D(const Normal& N, const Vector& H, const float m);

  // Geometry term [CT82]
  float G(const Vector& L, const Vector& H, const Vector& E, const Normal& N);

  // Integrate over a hemisphere as in [DJ2006], eq. 2
  float integrate_frTS(const Point& x, const Vector& wo, const Normal& N, const float F);
};
