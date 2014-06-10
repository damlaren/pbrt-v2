/*
 * A material for the tastiest part of the human body.
 * Based on "A Spectral BSSRDF for Shading Human Skin" [Donner & Jensen 2006].
 *
 * The material is adapted from the Subsurface material, but has a different
 * BRDF and parameterization.
 */


#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKIN_H
#define PBRT_MATERIALS_SKIN_H

// materials/subsurface.h*
#include "pbrt.h"
#include "material.h"
#include "core/reflection.h"

// The BRDF defined in DF2006 is a Torrance-Sparrowe micofacet
// BRDF modulated by an additional "oiliness" parameter.
class SkinMicrofacet : public Microfacet {
 public:
  SkinMicrofacet(const Spectrum &reflectance, Fresnel *f, MicrofacetDistribution *d, float oil);
  Spectrum f(const Vector &wo, const Vector &wi) const;
  Spectrum Sample_f(const Vector &wo, Vector *wi,
		    float u1, float u2, float *pdf) const;
  float Pdf(const Vector &wo, const Vector &wi) const;

  float oiliness; // Unctuousness. Must be between 0 and 1.
};

/*
class Beckmann : public MicrofacetDistribution {
public:
    Beckmann(float roughness) {
      m = roughness;
    }
    // Beckmann Public Methods
    float D(const Vector &wh) const;
    virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
    virtual float Pdf(const Vector &wi, const Vector &wo) const;
private:
    float m;
};
*/

// SubsurfaceMaterial Declarations
class SkinMaterial : public Material {
public:
    // SkinMaterial Public Methods
    SkinMaterial(float sc, Reference<Texture<Spectrum> > kr,
		 Reference<Texture<Spectrum> > sa,
		 Reference<Texture<Spectrum> > sps,
		 Reference<Texture<float> > e,
		 Reference<Texture<float> > bump,
		 float oil) {
        scale = sc;
        Kr = kr;
        sigma_a = sa;
        sigma_prime_s = sps;
        eta = e;
        bumpMap = bump;
	oiliness = oil;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom,
                      const DifferentialGeometry &dgShading,
                      MemoryArena &arena) const;
private:
    // SkinMaterial Private Data
    float scale;
    float oiliness;
    Reference<Texture<Spectrum> > Kr, sigma_a, sigma_prime_s;
    Reference<Texture<float> > eta, bumpMap;
};


SkinMaterial *CreateSkinMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_SKIN_H
