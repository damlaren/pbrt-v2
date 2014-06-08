/*
  An implementation of subsurface scattering as described in:
  "Photon Beam Diffusion: A Hybrid Monte Carlo Method for Subsurface Scattering"
  by Habel, Christensen, and Jarosz, 2013:
  http://graphics.pixar.com/library/PhotonBeamDiffusion/paper.pdf

  TODO: right now this is just a straight copy of the DipoleSubsurfaceIntegrator
 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
#define PBRT_INTEGRATORS_DIPOLESUBSURFACE_H

// integrators/pbdsubsurface.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "renderers/surfacepoints.h"
struct SubsurfaceOctreeNode;

// PBDSubsurfaceIntegrator Helper Declarations
struct IrradiancePoint {
    IrradiancePoint() { }
    IrradiancePoint(const SurfacePoint &sp, const Spectrum &ee)
        : p(sp.p), n(sp.n), E(ee), area(sp.area),
          rayEpsilon(sp.rayEpsilon) { }
    Point p;
    Normal n;
    Spectrum E;
    float area, rayEpsilon;
};



// PBDSubsurfaceIntegrator Declarations
class PBDSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // PBDSubsurfaceIntegrator Public Methods
    PBDSubsurfaceIntegrator(int mdepth, float merror, float mindist,
			    const string &fn) {
        maxSpecularDepth = mdepth;
        maxError = merror;
        minSampleDist = mindist;
        filename = fn;
        octree = NULL;
    }
    ~PBDSubsurfaceIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);
private:
    // PBDSubsurfaceIntegrator Private Data
    int maxSpecularDepth;
    float maxError, minSampleDist;
    string filename;
    vector<IrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    SubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


PBDSubsurfaceIntegrator *CreatePBDSubsurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
