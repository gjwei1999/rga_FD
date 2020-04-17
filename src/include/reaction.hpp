/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;

  std::unique_ptr<TLorentzVector> _x_mu;

  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _boosted_gamma;
  std::unique_ptr<TLorentzVector> _boosted_prot;
  std::unique_ptr<TLorentzVector> _boosted_pip;
  std::unique_ptr<TLorentzVector> _boosted_pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;

  float _weight = NAN;

  bool _mc = false;

  bool _is_boosted = false;

  bool _hasE = false;
  bool _hasP = false;

  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  short _numPart = 0;
  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  short _sector = -1;

  float _MM = NAN;
  float _MM2 = NAN;

  float _nu;
  float _theta_e;
  float _epsilont;
  float _flux;
  
  float _W = NAN;
  float _Q2 = NAN;

  float _W_e_Prot = NAN;
  float _Q2_e_Prot = NAN;

  float _inv_Ppip = NAN;
  float _inv_Ppim = NAN;
  float _inv_pip_pim = NAN;
  float _W_P2pi = NAN;

  float _phi_gamma = NAN;
  float _phi_prot = NAN;
  float _phi_pip = NAN;
  float _phi_pim = NAN;

  float _alpha_ppip_pipim = NAN;
  float _alpha_pippim_pipf = NAN;
  float _alpha_ppim_pipip = NAN;

  void SetElec();

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();

  inline bool mc() { return _mc; }
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  // float q_3_();
  // TLorentzVector p_mu_prime_cm();
  // TLorentzVector pip_mu_prime_cm();
  // TLorentzVector pim_mu_prime_cm();
  // float theta_();

  void boost();
  void CalcMissMass();
  float MM();
  float MM2();

  float weight();
//  inline float weight() {
//    return  //_data->mc_weight();  //
//        1.0;
//  }
  void epsilont();
  float flux();

  float inv_Ppip();
  float inv_Ppim();
  float inv_pip_pim();
  float w_P2pi_rec();

  void W_2pi_P();
  void invMassPpip();
  void invMassPpim();
  void invMasspippim();

  float prot_theta_lab();
  float pip_theta_lab();
  float pim_theta_lab();

//
//
//modified by Jiawei
//

  float theta_Ppip_lab();
  float theta_Ppim_lab();
  float theta_pippim_lab();

  float prot_phi_lab();
  float pip_phi_lab();
  float pim_phi_lab();
    
  float prot_p();
  float pip_p();
  float pim_p();

////
////
////

  float prot_theta();
  float pip_theta();
  float pim_theta();

  void AlphaCalc();

  float gamma_Phi();
  float prot_Phi();
  float pip_Phi();
  float pim_Phi();
  float alpha_ppip_pipim();
  float alpha_pippim_pipf();
  float alpha_ppim_pipip();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  float theta_x_mu();
  float E_x_mu();
  float P_x_mu();
  float theta_beam();

  inline short sec() { return _data->dc_sec(0); }
  inline int det() { return abs(_data->status(0) / 1000); }

  inline bool TwoPion() {
    return ((_numPip == 1 && _numPim == 1) && (_hasE && _hasP && _hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool ProtonPim() {
    return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && !_hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool SinglePip() { return ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim && !_hasOther)); }
  inline bool SingleP() {
    return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim && !_hasNeutron && !_hasOther));
  }

  inline bool NeutronPip() {
    bool _channel = true;
    _channel &= ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim && _hasNeutron)) ||
                (Reaction::SinglePip() && Reaction::MM() >= 0.85 && Reaction::MM() <= 1.0);
    return _channel;
  }

  const TLorentzVector &e_mu() { return *_beam; }
  const TLorentzVector &e_mu_prime() { return *_elec; }
  const TLorentzVector &gamma() { return *_gamma; }
};

class MCReaction : public Reaction {
 private:
  float _weight_mc = NAN;
  float _W_mc = NAN;
  float _Q2_mc = NAN;

  float _MCinv_Ppip = NAN;
  float _MCinv_Ppim = NAN;
  float _MCinv_pip_pim = NAN;

  std::unique_ptr<TLorentzVector> _elec_mc;
  std::unique_ptr<TLorentzVector> _gamma_mc;
  std::unique_ptr<TLorentzVector> _prot_mc;
  std::unique_ptr<TLorentzVector> _pip_mc;
  std::unique_ptr<TLorentzVector> _pim_mc;
  std::unique_ptr<TLorentzVector> _boosted_gamma_mc;
  std::unique_ptr<TLorentzVector> _boosted_prot_mc;
  std::unique_ptr<TLorentzVector> _boosted_pip_mc;
  std::unique_ptr<TLorentzVector> _boosted_pim_mc;

  bool _is_boosted_mc = false;

  float _MM_mc = NAN;
  float _MM2_mc = NAN;

  float _alpha_ppip_pipim_mc = NAN;
  float _alpha_pippim_pipf_mc = NAN;
  float _alpha_ppim_pipip_mc = NAN;

  float _alpha_ppip_pipim_thrown_mc = NAN;
  float _alpha_pippim_pipf_thrown_mc = NAN;
  float _alpha_ppim_pipip_thrown_mc = NAN;

 public:
  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  void SetMCElec();
  inline float weight() { return _data->mc_weight(); }
  inline float W_mc() { return _W_mc; }
  inline float Q2_mc() { return _Q2_mc; }
  void CalcMissMass_mc();
  float MM_mc();
  float MM2_mc();
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);

  void boost_mc();

  float MCinv_Ppip();
  float MCinv_Ppim();
  float MCinv_pip_pim();

  float MCprot_theta();
  float MCpip_theta();
  float MCpim_theta();

  float MCprot_theta_lab();
  float MCpip_theta_lab();
  float MCpim_theta_lab();

  void MCinvMassPpip();
  void MCinvMassPpim();
  void MCinvMasspippim();

  void MCAlphaCalc();

  float MCgamma_Phi();
  float MCprot_Phi();
  float MCpip_Phi();
  float MCpim_Phi();

  float MCalpha_ppip_pipim();
  float MCalpha_pippim_pipf();
  float MCalpha_ppim_pipip();

  float MCalpha_ppip_pipim_thrown();
  float MCalpha_pippim_pipf_thrown();
  float MCalpha_ppim_pipip_thrown();

  float MCprot_theta_thrown();
  float MCpip_theta_thrown();
  float MCpim_theta_thrown();

  float MCgamma_Phi_thrown();
  float MCprot_Phi_thrown();
  float MCpip_Phi_thrown();
  float MCpim_Phi_thrown();
};

#endif
