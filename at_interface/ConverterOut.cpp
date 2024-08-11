#include "ConverterOut.hpp"

#include "PFSimpleTask.hpp"

#include "AnalysisTree/DataHeader.hpp"
#include "AnalysisTree/Matching.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <algorithm>

void ConverterOut::CopyParticle(const OutputContainer& kf_particle, AnalysisTree::Particle& particle) const {

  particle.SetMomentum(kf_particle.GetPx(), kf_particle.GetPy(), kf_particle.GetPz());
  particle.SetField(kf_particle.GetMass(), invmass_field_id_);
  particle.SetPid(kf_particle.GetPdg());

  particle.SetField(kf_particle.GetXDecay(), x_decay_field_id_);
  particle.SetField(kf_particle.GetYDecay(), y_decay_field_id_);
  particle.SetField(kf_particle.GetZDecay(), z_decay_field_id_);
  particle.SetField(kf_particle.GetXPCA(), x_pca_field_id_);
  particle.SetField(kf_particle.GetYPCA(), y_pca_field_id_);
  particle.SetField(kf_particle.GetZPCA(), z_pca_field_id_);
  if(is_save_error_values_) {
    particle.SetField(kf_particle.GetXDecayError(), x_decay_error_field_id_);
    particle.SetField(kf_particle.GetYDecayError(), y_decay_error_field_id_);
    particle.SetField(kf_particle.GetZDecayError(), z_decay_error_field_id_);
    particle.SetField(kf_particle.GetXPCAError(), x_pca_error_field_id_);
    particle.SetField(kf_particle.GetYPCAError(), y_pca_error_field_id_);
    particle.SetField(kf_particle.GetZPCAError(), z_pca_error_field_id_);
    particle.SetField(kf_particle.GetPtError(), pt_err_field_id_);
    particle.SetField(kf_particle.GetPhiError(), phi_err_field_id_);
    particle.SetField(kf_particle.GetEtaError(), eta_err_field_id_);
    particle.SetField(kf_particle.GetMassError(), invmass_err_field_id_);
  }

  for (int i = 0; i < decay_.GetNDaughters(); ++i) {
    if(is_save_topo_vars_) {
      particle.SetField(kf_particle.GetChi2Prim(i), chi2prim_field_id_ + i);
      particle.SetField(kf_particle.GetCos(i), cosine_field_id_ + i);
    }
    particle.SetField(kf_particle.GetDaughterIds().at(i), daughter_id_field_id_ + i);
  }

  if(is_save_topo_vars_) particle.SetField(kf_particle.GetDistance(), distance_field_id_);

  if (decay_.GetNDaughters() == 3 && is_save_topo_vars_) {
    particle.SetField(kf_particle.GetDistanceToSV(), distance_field_id_ + 1);
    for (int i = 0; i < decay_.GetNDaughters(); ++i) {
      particle.SetField(kf_particle.GetChi2Geo(i + 1), chi2geo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetChi2Topo(i + 1), chi2topo_sm_field_id_ + i);
      particle.SetField(kf_particle.GetCosineTopo(i + 1), cosine_topo_sm_field_id_ + i);
    }
  }

  if(is_save_topo_vars_) {
    particle.SetField(kf_particle.GetChi2Geo(0), chi2geo_field_id_);
    particle.SetField(kf_particle.GetL(), l_field_id_);
    particle.SetField(kf_particle.GetLdL(), l_over_dl_field_id_);
    particle.SetField(kf_particle.GetChi2Topo(0), chi2_topo_field_id_);
    particle.SetField(kf_particle.GetCosineTopo(0), cosine_topo_field_id_);
  }
}

void ConverterOut::Exec() {

  candidates_ = pfsimple_task_->GetSimpleFinder()->GetCandidates();

  lambda_reco_->ClearChannels();
  lambda_sim_->ClearChannels();
  lambda_reco2sim_->Clear();

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  events_->SetField(float(sim_events_->GetField<float>(b_field_id_)), 0);//TODO

  const auto& br_conf = out_config->GetBranchConfig(lambda_reco_->GetId());

  for (const auto& candidate : candidates_) {

    AnalysisTree::Particle particle(lambda_reco_->GetNumberOfChannels(), br_conf);
    CopyParticle(candidate, particle);
    if (mc_particles_) {
      MatchWithMc(particle);
    }

    bool is_write = true;
    if (output_cuts_) {
      is_write = output_cuts_->Apply(particle);
    }

    if (is_write) {
      auto& lambdarec = lambda_reco_->AddChannel(br_conf);
      lambdarec = particle;
    }
  }
  delete pfsimple_task_->GetSimpleFinder();
}

void ConverterOut::Init() {

  this->SetInputBranchNames({sim_events_name_, rec_tracks_name_, mc_particles_name_});

  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();

  sim_events_ = ANALYSISTREE_UTILS_GET<AnalysisTree::EventHeader*>(chain->GetPointerToBranch(sim_events_name_));
  mc_particles_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch(mc_particles_name_));
  rec_to_mc_ = chain->GetMatchPointers().find(config_->GetMatchName(rec_tracks_name_, mc_particles_name_))->second;

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  std::string out_branch_event = "Events";
  std::string out_branch = std::string("Candidates");
  std::string out_branch_sim = std::string("Simulated");
  std::string out_branch_reco2sim = out_branch + "2" + out_branch_sim;

  AnalysisTree::BranchConfig EventBranch(out_branch_event, AnalysisTree::DetType::kEventHeader);
  EventBranch.AddField<float>("b", "Impact parameter, fm");

  AnalysisTree::BranchConfig out_particles(out_branch, AnalysisTree::DetType::kParticle);
  out_particles.AddField<float>("mass_inv", "Invariant mass of the candidate, GeV/c^2");

  out_particles.AddField<float>("x_decay", "X coordinate of the decay point of the candidate, cm");
  out_particles.AddField<float>("y_decay", "Y coordinate of the decay point of the candidate, cm");
  out_particles.AddField<float>("z_decay", "Z coordinate of the decay point of the candidate, cm");
  out_particles.AddField<float>("x_pca", "X coordinate of the point of candidate's closest approach to the primary vertex, cm");
  out_particles.AddField<float>("y_pca", "Y coordinate of the point of candidate's closest approach to the primary vertex, cm");
  out_particles.AddField<float>("z_pca", "Z coordinate of the point of candidate's closest approach to the primary vertex, cm");
  if(is_save_error_values_) {
    out_particles.AddField<float>("x_decay_err", "Estimate of x_decay error, cm");
    out_particles.AddField<float>("y_decay_err", "Estimate of y_decay error, cm");
    out_particles.AddField<float>("z_decay_err", "Estimate of z_decay error, cm");
    out_particles.AddField<float>("x_pca_err", "Estimate of x_pca error, cm");
    out_particles.AddField<float>("y_pca_err", "Estimate of y_pca error, cm");
    out_particles.AddField<float>("z_pca_err", "Estimate of z_pca error, cm");
    out_particles.AddField<float>("pT_err", "Estimate of transverse momentum error, GeV/c");
    out_particles.AddField<float>("phi_err", "Estimate of azimuthal angle error, rad");
    out_particles.AddField<float>("eta_err", "Estimate of pseudorapidity error");
    out_particles.AddField<float>("mass_inv_err", "Estimate of invariant mass error, GeV/c^2");
  }

  if (decay_.GetNDaughters() == 3) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id", "daughter3_id"}, "Daughter's id according to enumeration in the input file with reconstructed tracks");
    if(is_save_topo_vars_) {
      out_particles.AddFields<float>({"chi2_daughter1", "chi2_daughter2", "chi2_daughter3"}, "Chi2 between daughter track and primary vertex");
      out_particles.AddField<float>("distance", "Distance of closest approach between first two daughters, cm");
      out_particles.AddField<float>("distance_sv", "Distance of closest approach between third daughter and\n"
                                                   "temporary secondary vertex constructed from first two daughters, cm");
      out_particles.AddFields<float>({"cosine_daughter1", "cosine_daughter2", "cosine_daughter3"}, "Cosine of the angle between mother and daughter momenta");
      out_particles.AddFields<float>({"chi2_geo_sm1", "chi2_geo_sm2", "chi2_geo_sm3"}, "Chi2 between two daughter tracks:\n"
                                                                                       "sm1 - {first and second}, sm2 - {first and third}, sm3 - {second and third}");
      out_particles.AddFields<float>({"chi2_topo_sm1", "chi2_topo_sm2", "chi2_topo_sm3"}, "Chi2 between intermediate mother track, constructed from 2 daughters, and primary vertex:\n"
                                                                                          "sm1 - {first and second}, sm2 - {first and third}, sm3 - {second and third}");
      out_particles.AddFields<float>({"cosine_topo_sm1", "cosine_topo_sm2", "cosine_topo_sm3"}, "Cosine of the pointing angle (intermediate mother track, constructed from 2 daughters, and shortest line connecting primary and decay vertices):\n"
                                                                                                "sm1 - {first and second}, sm2 - {first and third}, sm3 - {second and third}");
    }
  } else if (decay_.GetNDaughters() == 2) {
    out_particles.AddFields<int>({"daughter1_id", "daughter2_id"}, "Daughter's id according to enumeration in the input file with reconstructed tracks");
    if(is_save_topo_vars_) {
      out_particles.AddFields<float>({"chi2_daughter1", "chi2_daughter2"}, "Chi2 between daughter track and primary vertex");
      out_particles.AddField<float>("distance", "Distance of closest approach between two daughters, cm");
      out_particles.AddFields<float>({"cosine_daughter1", "cosine_daughter2"}, "Cosine of the angle between mother and daughter momenta");
    }
  }

  if(is_save_topo_vars_) {
    out_particles.AddField<float>("chi2_geo", "Chi2 between two daughter tracks");
    out_particles.AddField<float>("l", "Signed distance along mother's trajectory from the point of closaest approach\nto the primary vertex to the decay vertex, cm");
    out_particles.AddField<float>("l_over_dl", "Value of l divided by its error estimate");
    out_particles.AddField<float>("chi2_topo", "Chi2 between mother track and primary vertex");
    out_particles.AddField<float>("cosine_topo", "Cosine of the pointing angle (mother track and shortest line connecting primary and decay vertices)");
  }

  AnalysisTree::BranchConfig LambdaSimBranch(out_branch_sim, AnalysisTree::DetType::kParticle);

  if (mc_particles_) {
    out_particles.AddField<int>("type", "MC-status of the candidate;\n"
                                        "Non-positive for BG, 1 for primary particles,\n"
                                        "2 for secondary, 3 for tertiary, etc. For details see\n"
                                        "https://github.com/HeavyIonAnalysis/PFSimple/blob/master/docs/type_field_description.pdf");
    LambdaSimBranch.AddField<int>("geant_process_id", "Id of the process of candidate particle production, see https://root.cern/doc/v624/TMCProcess_8h.html");
  }

  man->AddBranch(events_, EventBranch);
  man->AddBranch(lambda_reco_, out_particles);
  man->AddBranch(lambda_sim_, LambdaSimBranch);
  man->AddMatching(out_branch, out_branch_sim, lambda_reco2sim_);

  if (output_cuts_)
    output_cuts_->Init(*out_config);

  events_->Init(EventBranch);
  InitIndexes();
}

int ConverterOut::GetMothersSimId(AnalysisTree::Particle& lambdarec) {
  std::vector<int> daughter_sim_id;
  for (int i = 0; i < decay_.GetNDaughters(); i++)
    daughter_sim_id.push_back(rec_to_mc_->GetMatch(lambdarec.GetField<int>(daughter_id_field_id_ + i)));

  if (*std::min_element(daughter_sim_id.begin(), daughter_sim_id.end()) < 0)// at least one daughter has no matching with mc
    return -1;

  std::vector<int> mother_sim_id;
  for (int i = 0; i < decay_.GetNDaughters(); i++)
    mother_sim_id.push_back(mc_particles_->GetChannel(daughter_sim_id.at(i)).GetField<int>(mother_id_field_id_));

  if (*std::min_element(mother_sim_id.begin(), mother_sim_id.end()) != *std::max_element(mother_sim_id.begin(), mother_sim_id.end()))// daughters belong to not the same mother
    return -1;

  if (mother_sim_id.at(0) < 0)// mother has negative id
    return -1;

  if (mc_particles_->GetChannel(mother_sim_id.at(0)).GetPid() != lambdarec.GetPid())// mother has not PDG which was supposed
    return -1;

  return mother_sim_id.at(0);
}

int ConverterOut::DetermineGeneration(int mother_sim_id) {
  int generation = 0;
  int older_id = mother_sim_id;
  while (older_id >= 0) {
    const auto& simtrackolder = mc_particles_->GetChannel(older_id);
    older_id = simtrackolder.GetField<int>(mother_id_field_id_);
    generation++;
  }

  return generation;
}

void ConverterOut::MatchWithMc(AnalysisTree::Particle& lambdarec) {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  int mother_id = GetMothersSimId(lambdarec);
  int generation = DetermineGeneration(mother_id);
  if(is_detailed_bg_ && generation==0) generation = DetermineBGType(lambdarec);
  lambdarec.SetField(generation, generation_field_id_);

  if (generation < 1) return;

  const AnalysisTree::Particle& simtrackmother = mc_particles_->GetChannel(mother_id);

  auto& lambdasim = lambda_sim_->AddChannel(out_config->GetBranchConfig(lambda_sim_->GetId()));

  lambdasim.SetMomentum(simtrackmother.GetPx(), simtrackmother.GetPy(), simtrackmother.GetPz());
  lambdasim.SetPid(simtrackmother.GetPid());
  lambdasim.SetField(simtrackmother.GetField<int>(g4process_field_id_), g4process_field_id_w_);
  lambda_reco2sim_->AddMatch(lambdarec.GetId(), lambdasim.GetId());
}

void ConverterOut::InitIndexes() {

  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();

  const auto& out_branch_reco = out_config->GetBranchConfig(lambda_reco_->GetId());
  const auto& out_branch_sim = out_config->GetBranchConfig(lambda_sim_->GetId());

  auto branch_conf_sim_event = config_->GetBranchConfig(sim_events_name_);
  b_field_id_ = branch_conf_sim_event.GetFieldId("b");

  x_decay_field_id_ = out_branch_reco.GetFieldId("x_decay");
  y_decay_field_id_ = out_branch_reco.GetFieldId("y_decay");
  z_decay_field_id_ = out_branch_reco.GetFieldId("z_decay");
  x_pca_field_id_ = out_branch_reco.GetFieldId("x_pca");
  y_pca_field_id_ = out_branch_reco.GetFieldId("y_pca");
  z_pca_field_id_ = out_branch_reco.GetFieldId("z_pca");
  x_decay_error_field_id_ = out_branch_reco.GetFieldId("x_decay_err");
  y_decay_error_field_id_ = out_branch_reco.GetFieldId("y_decay_err");
  z_decay_error_field_id_ = out_branch_reco.GetFieldId("z_decay_err");
  x_pca_error_field_id_ = out_branch_reco.GetFieldId("x_pca_err");
  y_pca_error_field_id_ = out_branch_reco.GetFieldId("y_pca_err");
  z_pca_error_field_id_ = out_branch_reco.GetFieldId("z_pca_err");
  daughter_id_field_id_ = out_branch_reco.GetFieldId("daughter1_id");

  pt_err_field_id_ = out_branch_reco.GetFieldId("pT_err");
  phi_err_field_id_ = out_branch_reco.GetFieldId("phi_err");
  eta_err_field_id_ = out_branch_reco.GetFieldId("eta_err");
  invmass_err_field_id_ = out_branch_reco.GetFieldId("mass_inv_err");

  if (mc_particles_) {
    auto branch_conf_sim = config_->GetBranchConfig(mc_particles_name_);
    mother_id_field_id_ = branch_conf_sim.GetFieldId("mother_id");
    g4process_field_id_ = branch_conf_sim.GetFieldId("geant_process_id");
    generation_field_id_ = out_branch_reco.GetFieldId("type");
    g4process_field_id_w_ = out_branch_sim.GetFieldId("geant_process_id");
  }

  invmass_field_id_ = out_branch_reco.GetFieldId("mass_inv");
  chi2prim_field_id_ = out_branch_reco.GetFieldId("chi2_daughter1");
  distance_field_id_ = out_branch_reco.GetFieldId("distance");
  cosine_field_id_ = out_branch_reco.GetFieldId("cosine_daughter1");

  chi2geo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_geo_sm1");
  chi2topo_sm_field_id_ = out_branch_reco.GetFieldId("chi2_topo_sm1");
  cosine_topo_sm_field_id_ = out_branch_reco.GetFieldId("cosine_topo_sm1");

  chi2geo_field_id_ = out_branch_reco.GetFieldId("chi2_geo");
  l_field_id_ = out_branch_reco.GetFieldId("l");
  l_over_dl_field_id_ = out_branch_reco.GetFieldId("l_over_dl");
  chi2_topo_field_id_ = out_branch_reco.GetFieldId("chi2_topo");
  cosine_topo_field_id_ = out_branch_reco.GetFieldId("cosine_topo");
}

std::pair<int, int> ConverterOut::DetermineDaughtersMCStatus(const int daughter_rec_id, const Pdg_t mother_expected_pdg) const {
  const int daughter_sim_id = rec_to_mc_->GetMatch(daughter_rec_id);
  if(daughter_sim_id < 0) return std::make_pair(1, -999); // no match to MC

  auto& daughter_sim = mc_particles_->GetChannel(daughter_sim_id);
  const int mother_sim_id = daughter_sim.GetField<int>(mother_id_field_id_);
  if(mother_sim_id < 0) return std::make_pair(2, -999); // daughter is primary

  const int geant_process = daughter_sim.GetField<int>(g4process_field_id_);
  auto& mother_sim = mc_particles_->GetChannel(mother_sim_id);
  auto mother_pdg = mother_sim.GetPid();

  int daughter_status{-999};

  if(geant_process != 4 && mother_pdg != mother_expected_pdg) daughter_status = 3; // daughter not from decay, mother's pdg unexpected
  if(geant_process != 4 && mother_pdg == mother_expected_pdg) daughter_status = 4; // daughter not from decay, mother's pdg expected
  if(geant_process == 4 && mother_pdg != mother_expected_pdg) daughter_status = 5; // daughter from decay, mother's pdg unexpected
  if(geant_process == 4 && mother_pdg == mother_expected_pdg) daughter_status = 6; // daughter from decay, mother's pdg expected

  return std::make_pair(daughter_status, mother_sim_id);
}

int ConverterOut::DetermineMotherMCStatus(const int mid1, const int mid2) {
  if(mid1 == -999 || mid2 == -999) return 0;
  if(mid1 == mid2) return 1;
  else             return 2;
}

int ConverterOut::DetermineBGType(AnalysisTree::Particle& particle) {
  std::vector<std::pair<int, int>> daughters_statuses;
  for (int i = 0; i < decay_.GetNDaughters(); i++) {
    auto daughter_rec_id = particle.GetField<int>(daughter_id_field_id_ + i);
    daughters_statuses.emplace_back(DetermineDaughtersMCStatus(daughter_rec_id, particle.GetPid()));
  }

  int result{0};
  int decimal{1};
  for(auto& ds : daughters_statuses) {
    result += decimal * ds.first;
    decimal *= 10;
  }

  std::vector<int> common_mother_statuses;
  common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(0).second, daughters_statuses.at(1).second));

  if(daughters_statuses.size() == 3) {
    common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(1).second, daughters_statuses.at(2).second));
    common_mother_statuses.emplace_back(DetermineMotherMCStatus(daughters_statuses.at(0).second, daughters_statuses.at(2).second));
  }

  decimal = 1000;
  for(auto& cms : common_mother_statuses) {
    result += decimal * cms;
    decimal *= 10;
  }

  return -result;
}
