/**
 * @class SimpleFinder
 * @brief V0 particle (Lambda) reconstruction algorithm
 * @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 *
 * Simplified version of KFParticleFinder class. At the current moment it is developed to reconstruct 2- and 3-body-decays. \n
 * SimpleFinder is based on KFParticle package, it uses mathematical apparatus implemented in KFParticle.
 * The basic idea and reconstruction algorithm also reproduce the KFParticle, but SimpleFinder is free of overloading of too complicated procedure as in KFParticleFinder. \n
 * Also the advantage of SimpleFinder is that reconstruction procedure is clear and under full control of user, almost in the "hand mode".
 * It gives a possibility of detailed analysis of V0 reconstruction, in particular of decay parameters distributions in order to optomize selection criterias (cuts).\n
 */

#ifndef KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_

#include <vector>

#include <KFParticleSIMD.h>
#include <KFVertex.h>
#include <TVector3.h>
#include <map>

#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "Constants.hpp"
#include "Decay.hpp"
#include "NonLinearCutBase.hpp"

class SimpleFinder {

  typedef std::array<float, 8> Param_t;
  typedef std::vector<Param_t> Parameters_t;

 public:
  SimpleFinder() = default;
  SimpleFinder(const SimpleFinder&) = default;
  SimpleFinder(SimpleFinder&&) = default;
  SimpleFinder& operator=(SimpleFinder&&) = default;
  SimpleFinder& operator=(const SimpleFinder&) = default;
  ~SimpleFinder() = default;

  void Init(std::vector<KFParticle>&& tracks, const KFVertex& pv);///< Initialize SimpleFinder object with PV and set of tracks of the current event
  void Init(const InputContainer& input);

  void InitIndexesMap();

  void FindParticles() {
    for (const auto& decay : decays_) {
      ReconstructDecay(decay);
    }
  }

  void ReconstructDecay(const Decay& decay);
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void AddDecay(const Decay& decay) { decays_.emplace_back(decay); }
  const std::vector<OutputContainer>& GetCandidates() const { return output_; }

 private:
  std::vector<KFParticle> tracks_;    ///< input information: vector of tracks
  KFVertex prim_vx_;                  ///< input information: primiry vertex
  std::array<float, 3> sec_vx_;       ///< input information: primiry vertex
  std::vector<Decay> decays_{};       ///< input information: list of decays to reconstruct
  NonLinearCutBase* ml_cuts_{nullptr};///< input information: non-linear cuts class (optional)

  std::map<Pdg_t, std::vector<int>> indexes_{};///< map of indexes for a given particle specie

  Parameters_t params_{};   ///< vector of daughter parameters at current SV estimation
  SelectionValues values_{};///< struct with mother and daughters properties used to apply cuts

  std::vector<OutputContainer> output_{};///< output information: vector of candidates

  /**
  * Find indexes of good daughters
  * @param cuts daughter particle cuts container
  * @return vector of indexes
  */
  std::vector<int> GetIndexes(const Daughter& cuts);

  bool IsGoodDaughter(const KFParticle& track, const Daughter& cuts);
  bool IsGoodPair(const KFParticle& track1, const KFParticle& track2, const Decay& decay);
  bool IsGoodThree(const KFParticle& track1, const KFParticle& track2, const KFParticle& track3, const Decay& decay);
  bool IsGoodMother(const KFParticleSIMD& mother, const Mother& cuts, int id_mother);
  bool IsMotherFromPV(const KFParticleSIMD& mother, const Mother& cuts, int id_mother);
  bool IsGoodDecayLength(const KFParticleSIMD& mother, const Mother& cuts);
  bool IsGoodCos(const KFParticleSIMD& mother, const Parameters_t& daughter_pars, const Decay& decay);

  KFParticleSIMD ConstructMother(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs);
  void CalculateSecondaryVertex();

  void CalculateParamsInPCA(const KFParticle& track1, int pid1, const KFParticle& track2, int pid2);
  void CalculateParamsInSV(const KFParticle& track, int pid, int id);
  float CalculateChiToPrimaryVertex(const KFParticle& track, Pdg_t pid) const;
  float CalculateCosTopo(const KFParticleSIMD& mother) const;
  static float CalculateDistanceBetweenParticles(const Parameters_t& parameters);
  float CalculateDistanceToSV(int id) const;

  void FillDaughtersInfo(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs);
  void SaveParticle(KFParticleSIMD& particle_simd, const Decay& decay);

  bool ApplyNonLinearCut() const {
    return ml_cuts_ == nullptr || ml_cuts_->ApplyCut(values_);
  }

  KFParticle GetTrack(int i) {
    return tracks_.at(i);
  }

  static KFPTrack ToKFPTrack(const KFParticle& particle);
  static void SetKFParticleEnergy(KFParticle& particle, int pdg) ;
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDERNEW_HPP_
