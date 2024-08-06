/**
 ** @class OutputContainer
 ** @brief Container with output information about reconstructed particles and geometrical decay parameters (quantities to be cut in order to select particles)
 ** @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 **
 ** Each particle candidate is characterized with set of geometrical decay parameters. Depending on the
 ** value of each parameter the candidate is saved or rejected.
 ** In order to save the reconstructed particle, the KFParticle object is used. It contains all
 ** information about the particle (mass, momentum etc), and access to this information is
 ** possible via KFParticle methods.
 **/

#ifndef OutputContainer_H
#define OutputContainer_H

#include "KFParticle.h"

#include "Constants.hpp"

class OutputContainer {
 public:
  OutputContainer() = default;
  explicit OutputContainer(const KFParticle& particle,
                           std::array<float, 3> decay_point,
                           std::array<float, 3> decay_point_err);

  virtual ~OutputContainer() = default;

  [[nodiscard]] const std::vector<int>& GetDaughterIds() const { return daughter_ids_; }
  [[nodiscard]] float GetPx() const { return px_; }
  [[nodiscard]] float GetPy() const { return py_; }
  [[nodiscard]] float GetPz() const { return pz_; }
  [[nodiscard]] float GetMass() const { return mass_; }
  [[nodiscard]] float GetPtError() const { return pt_error_; }
  [[nodiscard]] float GetPhiError() const { return phi_error_; }
  [[nodiscard]] float GetEtaError() const { return eta_error_; }
  [[nodiscard]] float GetMassError() const { return mass_error_; }
  [[nodiscard]] Pdg_t GetPdg() const { return pdg_; }

  [[nodiscard]] float GetChi2Prim(int i) const { return values_.chi2_prim[i]; }
  [[nodiscard]] float GetCos(int i) const { return values_.cos[i]; }
  [[nodiscard]] float GetChi2Geo(int i) const { return values_.chi2_geo[i]; }
  [[nodiscard]] float GetChi2Topo(int i) const { return values_.chi2_topo[i]; }
  [[nodiscard]] float GetDistance() const { return values_.distance; }
  [[nodiscard]] float GetDistanceToSV() const { return values_.distance_sv; }
  [[nodiscard]] float GetL() const { return values_.l; }
  [[nodiscard]] float GetLdL() const { return values_.l_over_dl; }
  [[nodiscard]] float GetDistanceToPVLine() const { return values_.distance_pv; }
  [[nodiscard]] float GetCosineTopo(int i) const { return values_.cos_topo[i]; }

  [[nodiscard]] float GetXDecay() const { return x_decay_; }
  [[nodiscard]] float GetYDecay() const { return y_decay_; }
  [[nodiscard]] float GetZDecay() const { return z_decay_; }
  [[nodiscard]] float GetXDecayError() const { return x_decay_error_; }
  [[nodiscard]] float GetYDecayError() const { return y_decay_error_; }
  [[nodiscard]] float GetZDecayError() const { return z_decay_error_; }

  [[nodiscard]] float GetXPCA() const { return x_pca_; }
  [[nodiscard]] float GetYPCA() const { return y_pca_; }
  [[nodiscard]] float GetZPCA() const { return z_pca_; }
  [[nodiscard]] float GetXPCAError() const { return x_pca_error_; }
  [[nodiscard]] float GetYPCAError() const { return y_pca_error_; }
  [[nodiscard]] float GetZPCAError() const { return z_pca_error_; }

  [[nodiscard]] int GetId() const { return id_; }
  [[nodiscard]] bool IsFromPV() const { return values_.is_from_PV; }

  void SetId(int id) { id_ = id; }

  void SetSelectionValues(const SelectionValues& v) { values_ = v; }

 protected:
  int id_{-1};
  std::vector<int> daughter_ids_{};

  float px_{-1.};
  float py_{-1.};
  float pz_{-1.};
  float mass_{-1.};
  Pdg_t pdg_{-1};

  float pt_error_{-1.};
  float phi_error_{-1.};
  float eta_error_{-1.};
  float mass_error_{-1.};

  float x_decay_{-1.};
  float y_decay_{-1.};
  float z_decay_{-1.};
  float x_decay_error_{-1.};
  float y_decay_error_{-1.};
  float z_decay_error_{-1.};

  float x_pca_{-1.};
  float y_pca_{-1.};
  float z_pca_{-1.};
  float x_pca_error_{-1.};
  float y_pca_error_{-1.};
  float z_pca_error_{-1.};

  //  int n_hits_{-1};

  SelectionValues values_{};
};

#endif// OutputContainer_H
