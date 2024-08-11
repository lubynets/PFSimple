#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_

#include "Decay.hpp"
#include "OutputContainer.hpp"

#include "AnalysisTree/Cuts.hpp"
#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/EventHeader.hpp"
#include "AnalysisTree/Task.hpp"

class PFSimpleTask;

class ConverterOut : public AnalysisTree::Task {
 public:
  explicit ConverterOut() = default;
  ~ConverterOut() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override {}

  void SetPFSimpleTask(PFSimpleTask* pfsimple_task) { pfsimple_task_ = pfsimple_task; }

  void CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const;
  void SetDecay(const Decay& decay) { decay_ = decay; }
  void SetOutputCuts(AnalysisTree::Cuts* output_cuts) { output_cuts_ = output_cuts; }

  void SetSimEventHeaderName(const std::string& name) { sim_events_name_ = name; }
  void SetRecTracksName(const std::string& name) { rec_tracks_name_ = name; }
  void SetSimTracksName(const std::string& name) { mc_particles_name_ = name; }
  void SetIsWriteDetailedBG(bool is=true) { is_detailed_bg_ = is; }
  void SetIsSaveTopologicalVariables(bool is=true) { is_save_topo_vars_ = is; }
  void SetIsSaveErrorValues(bool is=true) { is_save_error_values_ = is; }

 protected:
  void InitIndexes();
  void MatchWithMc(AnalysisTree::Particle& particle);
  int GetMothersSimId(AnalysisTree::Particle& lambdarec);
  int DetermineGeneration(int mother_sim_id);

  // first returned value is:
  // 1 - reco daughter is unmatched to mc
  // 2 - reco daughter is matched, but primary
  // 3 - reco daughter is secondary, produced not in decay from mother with not expected pdg
  // 4 - reco daughter is secondary, produced not in decay from mother with expected pdg
  // 5 - reco daughter is secondary, produced in decay from mother with not expected pdg
  // 6 - reco daughter is secondary, produced in decay from mother with expected pdg
  // second returned value is mother's sim id (if the first value is >2) or -999 (if mother doesn't exist)
  std::pair<int, int> DetermineDaughtersMCStatus(int daughter_rec_id, Pdg_t mother_expected_pdg) const;

  // 1 - daughters have the same mother
  // 2 - daughters have different mothers
  // 0 - at least one daughter does not have mother (e.g. primary)
  static int DetermineMotherMCStatus(int mid1, int mid2);

  int DetermineBGType(AnalysisTree::Particle& particle);

  // output branches
  AnalysisTree::EventHeader* events_{nullptr};
  AnalysisTree::Particles* lambda_reco_{nullptr};
  AnalysisTree::Particles* lambda_sim_{nullptr};
  AnalysisTree::Matching* lambda_reco2sim_{nullptr};

  // input branches
  std::string mc_particles_name_{"SimParticles"};
  std::string rec_tracks_name_{"RecParticles"};
  std::string sim_events_name_{"SimEventHeader"};

  AnalysisTree::Particles* mc_particles_{nullptr};
  AnalysisTree::Matching* rec_to_mc_{nullptr};
  AnalysisTree::EventHeader* sim_events_{nullptr};
  AnalysisTree::Cuts* output_cuts_{nullptr};
  Decay decay_{};

  std::vector<OutputContainer> candidates_;

  PFSimpleTask* pfsimple_task_{nullptr};

  // field ids of simulated events
  int b_field_id_{AnalysisTree::UndefValueInt};

  // field ids of input simulated mother
  int mother_id_field_id_{AnalysisTree::UndefValueInt};

  int x_decay_field_id_{AnalysisTree::UndefValueInt};
  int y_decay_field_id_{AnalysisTree::UndefValueInt};
  int z_decay_field_id_{AnalysisTree::UndefValueInt};
  int x_pca_field_id_{AnalysisTree::UndefValueInt};
  int y_pca_field_id_{AnalysisTree::UndefValueInt};
  int z_pca_field_id_{AnalysisTree::UndefValueInt};

  int x_decay_error_field_id_{AnalysisTree::UndefValueInt};
  int y_decay_error_field_id_{AnalysisTree::UndefValueInt};
  int z_decay_error_field_id_{AnalysisTree::UndefValueInt};
  int x_pca_error_field_id_{AnalysisTree::UndefValueInt};
  int y_pca_error_field_id_{AnalysisTree::UndefValueInt};
  int z_pca_error_field_id_{AnalysisTree::UndefValueInt};

  int daughter_id_field_id_{AnalysisTree::UndefValueInt};
  int generation_field_id_{AnalysisTree::UndefValueInt};
  int g4process_field_id_{AnalysisTree::UndefValueInt};
  int g4process_field_id_w_{-1};

  int pt_err_field_id_{AnalysisTree::UndefValueInt};
  int phi_err_field_id_{AnalysisTree::UndefValueInt};
  int eta_err_field_id_{AnalysisTree::UndefValueInt};
  int invmass_err_field_id_{AnalysisTree::UndefValueInt};

  int invmass_field_id_{AnalysisTree::UndefValueInt};
  int chi2prim_field_id_{AnalysisTree::UndefValueInt};
  int distance_field_id_{AnalysisTree::UndefValueInt};
  int cosine_field_id_{AnalysisTree::UndefValueInt};

  int chi2geo_field_id_{AnalysisTree::UndefValueInt};
  int l_field_id_{AnalysisTree::UndefValueInt};
  int l_over_dl_field_id_{AnalysisTree::UndefValueInt};
  int chi2_topo_field_id_{AnalysisTree::UndefValueInt};
  int cosine_topo_field_id_{AnalysisTree::UndefValueInt};

  int chi2geo_sm_field_id_{AnalysisTree::UndefValueInt};
  int chi2topo_sm_field_id_{AnalysisTree::UndefValueInt};
  int cosine_topo_sm_field_id_{AnalysisTree::UndefValueInt};

  bool is_detailed_bg_{false};
  bool is_save_topo_vars_{true};
  bool is_save_error_values_{false};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_CONVERTEROUT_H_
