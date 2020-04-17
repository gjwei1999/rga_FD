/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;  // atof(getenv("CLAS12_E"));

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();
  _prot = std::make_unique<TLorentzVector>();
  _x_mu = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

  _weight = // _data->mc_weight();  //
      1.0;
}

Reaction::~Reaction() {}

void Reaction::epsilont() {
    _nu = _beam_energy - sqrt(_data->px(0)*_data->px(0) + _data->py(0)*_data->py(0) + _data->pz(0)*_data->pz(0) + MASS_E*MASS_E);
    _theta_e = _elec->Theta();
    _epsilont = 1/(1 + 2*(1+_nu*_nu/_Q2)*TMath::Tan(_theta_e/2)*TMath::Tan(_theta_e/2));
}

float Reaction::flux() {
    _flux = 1/4/137/PI/_beam_energy/_beam_energy/MASS_P/MASS_P * _W*(_W*_W-MASS_P*MASS_P);
    _flux = _flux/_Q2 /(1-_epsilont);
    return _flux;
}

void Reaction::SetElec() {
  _numPart++;
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;
  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}
void Reaction::SetProton(int i) {
  _numPart++;
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}

void Reaction::SetPip(int i) {
  _numPart++;
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPart++;
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numPart++;
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->pid(i) == NEUTRON && abs(_data->chi2pid(i)) < 0.5)
    SetNeutron(i);
  else {
    _numPart++;
    _numOther++;
    _hasOther = true;
    _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();

  *mm += (*_gamma + *_target);
  *_x_mu += (*_gamma + *_target);

  if (SinglePip() || NeutronPip()) {
    *mm -= *_pip;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (TwoPion()) {
    *mm -= *_prot;
    *mm -= *_pip;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
    *_x_mu -= *_prot;
    *_x_mu -= *_pip;
    *_x_mu -= *_pim;
  } else if (ProtonPim()) {
    *mm -= *_prot;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (SingleP()) {
    *mm -= *_prot;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

float Reaction::P_x_mu() {
  if (TwoPion()) {
    if (_x_mu->P() != _x_mu->P()) CalcMissMass();
    return _x_mu->P();
  } else
    return NAN;
}
float Reaction::E_x_mu() {
  if (TwoPion()) {
    if (_x_mu->E() != _x_mu->E()) CalcMissMass();
    return _x_mu->E();
  } else
    return NAN;
}

float Reaction::theta_x_mu() {
  if (TwoPion()) {
    if (_x_mu->Theta() != _x_mu->Theta()) CalcMissMass();

    return (_x_mu->Theta() * 180.0 / PI);
  } else

    return NAN;
}
float Reaction::theta_beam() {
  if (TwoPion()) {
    return (_beam->Theta() * 180.0 / PI);
  } else

    return NAN;
}

float Reaction::prot_theta_lab() {
  if (TwoPion())
    return _prot->Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pip_theta_lab() {
  if (TwoPion())
    return _pip->Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pim_theta_lab() {
  if (TwoPion())
    return _pim->Theta() * 180.0 / PI;
  else
    return NAN;
}
//
//
//modified by Jiawei
//
    float Reaction::theta_Ppip_lab() {
        if (TwoPion()){
            double angle;
	    angle = _prot->Angle(_pip->Vect()) * 180.0 / PI;
            return angle;
        }
        else
            return NAN;
        
    }
    float Reaction::theta_Ppim_lab() {
        if (TwoPion()){
	    double angle;
            angle = _prot->Angle(_pim->Vect()) * 180.0 / PI;
            return angle;
	}
        else
            return NAN;
    }
    float Reaction::theta_pippim_lab(){
        if (TwoPion()){
	    double angle;
            angle = _pip->Angle(_pim->Vect()) * 180.0 / PI;
            return angle;
	}
        else
            return NAN;
    }

  float Reaction::prot_phi_lab() {
    if (TwoPion())
      return _prot->Phi() * 180.0 / PI;
    else
      return NAN;
}
  float Reaction::pip_phi_lab() {
    if (TwoPion())
      return _pip->Phi() * 180.0 / PI;
    else
      return NAN;
}
  float Reaction::pim_phi_lab() {
    if (TwoPion())
      return _pim->Phi() * 180.0 / PI;
    else
      return NAN;
}

  float Reaction::prot_p(){
    if(TwoPion())
        return _prot->P();
    else
        return NAN;
}
  float Reaction::pip_p(){
    if(TwoPion())
        return _pip->P();
    else
        return NAN;
}
  float Reaction::pim_p(){
    if(TwoPion())
        return _pim->P();
    else
        return NAN;
}
////
////
////



void Reaction::invMassPpim() {
  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_prot;
  *inv_Ppim += *_pim;
  if (TwoPion()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_pip;
  *inv_pip_pim += *_pim;
  if (TwoPion()) _inv_pip_pim = inv_pip_pim->M();
}
void Reaction::invMassPpip() {
  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_prot;
  *inv_Ppip += *_pip;
  if (TwoPion()) _inv_Ppip = inv_Ppip->M();
}

void Reaction::W_2pi_P() {
  auto W_P2pi = std::make_unique<TLorentzVector>();
  *W_P2pi += *_prot;
  *W_P2pi += *_pip;
  *W_P2pi += *_pim;

  if (TwoPion()) _W_P2pi = W_P2pi->M();
}

// float P_P_prime_calc_e(TLorentzVector e_mu_prime) {
//   TLorentzVector e_mu(0.0, 0.0, 2.2, 2.2);
//   TVector3 p_mu_3(0, 0, 0);
//   TLorentzVector p_mu;
//   p_mu.SetVectM(p_mu_3, MASS_P);
//   float E_prime = (e_mu.E() / (1 + ((e_mu.E() / MASS_P) * (1 - cos(e_mu_prime.Theta())))));
//
//   float E_p_prime = e_mu.E() + p_mu.E() - E_prime;
//   return sqrt(E_p_prime * E_p_prime -
//               MASS_P *
//                   MASS_P return sqrt(
//                       (e_mu.E() - e_mu_prime.E() + MASS_P) * (e_mu.E() - e_mu_prime.E() + MASS_P) +
//                       (e_mu.E() - e_mu_prime.E()) *
//                           (e_mu.E() - e_mu_prime.E())(e_mu.E() * e_mu.E() + e_mu_prime.E() * e_mu_prime.E()) +
//                       2 * cos(e_mu_prime.Theta()))(e_mu.E() * e_mu_prime.E() - (MASS_P * (e_mu.E() -
//                       e_mu_prime.E()))));
// }
float Reaction::inv_Ppip() {
  if (_inv_Ppip != _inv_Ppip) invMassPpip();
  return _inv_Ppip;
}
float Reaction::inv_Ppim() {
  if (_inv_Ppim != _inv_Ppim) invMassPpim();
  return _inv_Ppim;
}
float Reaction::inv_pip_pim() {
  if (_inv_pip_pim != _inv_pip_pim) invMasspippim();
  return _inv_pip_pim;
}
float Reaction::w_P2pi_rec() {
  if (_W_P2pi != _W_P2pi) W_2pi_P();
  return _W_P2pi;
}

void Reaction::boost() {
  _is_boosted = true;
  _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
  _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
  _boosted_pim = std::make_unique<TLorentzVector>(*_pim);
  TRotation rot;
  float_t beta_1 = ((sqrt(_gamma->E() * _gamma->E() + _Q2)) / (_gamma->E() + MASS_P));
  TVector3 uz = _gamma->Vect().Unit();                          // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
  ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot
  _boosted_prot->Transform(rot);
  _boosted_prot->Boost(0, 0, -beta_1);
  _boosted_pip->Transform(rot);
  _boosted_pip->Boost(0, 0, -beta_1);
  _boosted_pim->Transform(rot);
  _boosted_pim->Boost(0, 0, -beta_1);  // -beta ko value (0.5 to -0.5 huda samma
                                       // value aauchha nattra aaudyna)
}
 float Reaction::weight() { return _weight; }

float Reaction::prot_theta() {
  if (!_is_boosted) boost();
  if (TwoPion())
    return _boosted_prot->Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pip_theta() {
  if (!_is_boosted) boost();
  if (TwoPion())
    return _boosted_pip->Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pim_theta() {
  if (!_is_boosted) boost();
  if (TwoPion())
    return _boosted_pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::gamma_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion()) {
    if (_gamma->Phi() > 0)
      return _gamma->Phi() * 180 / PI;
    else if (_gamma->Phi() < 0)
      return (_gamma->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::prot_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion()) {
    if (_boosted_prot->Phi() > 0)
      return _boosted_prot->Phi() * 180 / PI;
    else if (_boosted_prot->Phi() < 0)
      return (_boosted_prot->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pip_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion()) {
    if (_boosted_pip->Phi() > 0)
      return _boosted_pip->Phi() * 180 / PI;
    else if (_boosted_pip->Phi() < 0)
      return (_boosted_pip->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pim_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion()) {
    if (_boosted_pim->Phi() > 0)
      return _boosted_pim->Phi() * 180 / PI;
    else if (_boosted_pim->Phi() < 0)
      return (_boosted_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}

void Reaction::AlphaCalc() {
  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm;
  float alpha_PIpPIm_pipf;
  float alpha_PPIm_piPIp;

  if (!_is_boosted) boost();

  // 1 this one is used for α[π−]
  a_gamma = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * V3_anti_z), 2)));  // V3_anti_z(0,0,-1);
  b_gamma = -(_boosted_pim->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pim->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_pim->Vect().Unit();

  alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pim->Vect() < 0) alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

  //α[pπ+][p'π−]
  /// 2
  a_gamma = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_prot->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_prot->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_prot->Vect().Unit();

  alpha_PIpPIm_pipf = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_prot->Vect() < 0) alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf;
  //α[pp'][π+π−]

  /// 3
  a_gamma = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_pip->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pip->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()), 2)));
  b_beta = -(_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pim->Vect().Unit() + b_beta * _boosted_pip->Vect().Unit();

  alpha_PPIm_piPIp = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pip->Vect() < 0) alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp;

  _alpha_ppip_pipim = alpha_PPIp_piPIm;
  _alpha_pippim_pipf = alpha_PIpPIm_pipf;
  _alpha_ppim_pipip = alpha_PPIm_piPIp;
}

float Reaction::alpha_ppip_pipim() {  // pipim bhaneko proton initial  pim ho?
  if (_alpha_ppip_pipim != _alpha_ppip_pipim) AlphaCalc();
  if (TwoPion())
    return _alpha_ppip_pipim;
  else
    return NAN;
}
float Reaction::alpha_pippim_pipf() {  // alpha P (proton initial proton final)
  if (_alpha_pippim_pipf != _alpha_pippim_pipf) AlphaCalc();
  if (TwoPion())
    return _alpha_pippim_pipf;
  else
    return NAN;
}
float Reaction::alpha_ppim_pipip() {  // alpha pip (proton initial pip)
  if (_alpha_ppim_pipip != _alpha_ppim_pipip) AlphaCalc();
  if (TwoPion())
    return _alpha_ppim_pipip;
  else
    return NAN;
}

MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_enrgy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_enrgy;
  _weight_mc = _data->mc_weight();
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _gamma_mc = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  _elec_mc = std::make_unique<TLorentzVector>();
  this->SetElec();
  this->SetMCElec();
  _prot_mc = std::make_unique<TLorentzVector>();
  _pip_mc = std::make_unique<TLorentzVector>();
  _pim_mc = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

void MCReaction::SetMCElec() {
  _hasE = true;
  _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);

  *_gamma_mc += *_beam - *_elec_mc;

  // Can calculate W and Q2 here
  _W_mc = physics::W_calc(*_beam, *_elec_mc);
  _Q2_mc = physics::Q2_calc(*_beam, *_elec_mc);
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }

void MCReaction::CalcMissMass_mc() {
  auto mm_mc = std::make_unique<TLorentzVector>();

  *mm_mc += (*_gamma_mc + *_target);
  *mm_mc -= *_prot_mc;
  *mm_mc -= *_pip_mc;
  *mm_mc -= *_pim_mc;
  _MM_mc = mm_mc->M();
  _MM2_mc = mm_mc->M2();
}
float MCReaction::MM_mc() {
  if (_MM_mc != _MM_mc) CalcMissMass_mc();
  return _MM_mc;
}
float MCReaction::MM2_mc() {
  if (_MM2_mc != _MM2) CalcMissMass_mc();
  return _MM2_mc;
}
float MCReaction::MCprot_theta_lab() { return _prot_mc->Theta() * 180.0 / PI; }
float MCReaction::MCpip_theta_lab() { return _pip_mc->Theta() * 180.0 / PI; }
float MCReaction::MCpim_theta_lab() { return _pim_mc->Theta() * 180.0 / PI; }

void MCReaction::boost_mc() {
  _is_boosted_mc = true;
  _boosted_prot_mc = std::make_unique<TLorentzVector>(*_prot_mc);
  _boosted_pip_mc = std::make_unique<TLorentzVector>(*_pip_mc);
  _boosted_pim_mc = std::make_unique<TLorentzVector>(*_pim_mc);
  TRotation rot_mc;
  float_t beta_1_mc = ((sqrt(_gamma_mc->E() * _gamma_mc->E() + _Q2_mc)) / (_gamma_mc->E() + MASS_P));
  TVector3 uz_mc = _gamma_mc->Vect().Unit();                          // uit vector along virtual photon
  TVector3 ux_mc = ((_beam->Vect()).Cross(_elec_mc->Vect())).Unit();  // unit vector along e cross e'
  ux_mc.Rotate(3. * PI / 2, uz_mc);        // rotating ux by 3pi/2 with uz as axis of roration
  rot_mc.SetZAxis(uz_mc, ux_mc).Invert();  // setting TRotation rot
  _boosted_prot_mc->Transform(rot_mc);
  _boosted_prot_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pip_mc->Transform(rot_mc);
  _boosted_pip_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pim_mc->Transform(rot_mc);
  _boosted_pim_mc->Boost(0, 0, -beta_1_mc);  // -beta ko value (0.5 to -0.5 huda samma
                                             // value aauchha nattra aaudyna)
}

float MCReaction::MCgamma_Phi() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion()) {
    if (_gamma_mc->Phi() > 0)
      return _gamma_mc->Phi() * 180 / PI;
    else if (_gamma_mc->Phi() < 0)
      return (_gamma_mc->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float MCReaction::MCprot_Phi() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion()) {
    if (_boosted_prot_mc->Phi() > 0)
      return _boosted_prot_mc->Phi() * 180 / PI;
    else if (_boosted_prot_mc->Phi() < 0)
      return (_boosted_prot_mc->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float MCReaction::MCpip_Phi() {
  if (TwoPion()) {
    if (_boosted_pip_mc->Phi() > 0)
      return _boosted_pip_mc->Phi() * 180 / PI;
    else if (_boosted_pip_mc->Phi() < 0)
      return (_boosted_pip_mc->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float MCReaction::MCpim_Phi() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion()) {
    if (_boosted_pim_mc->Phi() > 0)
      return _boosted_pim_mc->Phi() * 180 / PI;
    else if (_boosted_pim_mc->Phi() < 0)
      return (_boosted_pim_mc->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float MCReaction::MCprot_theta() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion())
    return _boosted_prot_mc->Theta() * 180.0 / PI;
  else
    return NAN;
}
float MCReaction::MCpip_theta() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion())
    return _boosted_pip_mc->Theta() * 180.0 / PI;
  else
    return NAN;
}
float MCReaction::MCpim_theta() {
  if (!_is_boosted_mc) boost_mc();
  if (TwoPion())
    return _boosted_pim_mc->Theta() * 180.0 / PI;
  else
    return NAN;
}

void MCReaction::MCinvMassPpip() {
  //  if (!_is_boosted_mc) boost_mc();
  auto MCinv_Ppip = std::make_unique<TLorentzVector>();
  *MCinv_Ppip += *_boosted_prot_mc;
  *MCinv_Ppip += *_boosted_pip_mc;
  _MCinv_Ppip = MCinv_Ppip->M();
}

void MCReaction::MCinvMassPpim() {
  //  if (!_is_boosted_mc) boost_mc();
  auto MCinv_Ppim = std::make_unique<TLorentzVector>();
  *MCinv_Ppim += *_boosted_prot_mc;
  *MCinv_Ppim += *_boosted_pim_mc;
  _MCinv_Ppim = MCinv_Ppim->M();
}
void MCReaction::MCinvMasspippim() {
  // if (!_is_boosted_mc) boost_mc();
  auto MCinv_pip_pim = std::make_unique<TLorentzVector>();
  *MCinv_pip_pim += *_boosted_pip_mc;
  *MCinv_pip_pim += *_boosted_pim_mc;
  _MCinv_pip_pim = MCinv_pip_pim->M();
}
float MCReaction::MCinv_Ppip() {
  if (_MCinv_Ppip != _MCinv_Ppip) MCinvMassPpip();
  // if (TwoPion())
  return _MCinv_Ppip;
  // else
  //   return NAN;
}
float MCReaction::MCinv_Ppim() {
  if (_MCinv_Ppim != _MCinv_Ppim) MCinvMassPpim();
  //  if (TwoPion())
  return _MCinv_Ppim;
  // else
  //  return NAN;
}
float MCReaction::MCinv_pip_pim() {
  if (_MCinv_pip_pim != _MCinv_pip_pim) MCinvMasspippim();
  //  if (TwoPion())
  return _MCinv_pip_pim;
  // else
  // return NAN;
}
void MCReaction::MCAlphaCalc() {
  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm_mc;
  float alpha_pippim_pipf_mc;
  float alpha_PPIm_piPIp_mc;

  if (!_is_boosted_mc) boost_mc();
  // 1 this one is used for
  a_gamma = sqrt(1. / (1 - pow((_boosted_pim_mc->Vect().Unit() * V3_anti_z), 2)));  // V3_anti_z(0,0,-1);
  b_gamma = -(_boosted_pim_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pim_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pim_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_pim_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip_mc->Vect().Unit() + b_beta * _boosted_pim_mc->Vect().Unit();

  alpha_PPIp_piPIm_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pim_mc->Vect() < 0) alpha_PPIp_piPIm_mc = 360. - alpha_PPIp_piPIm_mc;

  //α[pπ+][p'π−]
  /// 2
  a_gamma = sqrt(1. / (1 - pow((_boosted_prot_mc->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_prot_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_prot_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_prot_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_prot_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip_mc->Vect().Unit() + b_beta * _boosted_prot_mc->Vect().Unit();

  alpha_pippim_pipf_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_prot_mc->Vect() < 0) alpha_pippim_pipf_mc = 360. - alpha_pippim_pipf_mc;
  //α[pp'][π+π−]

  /// 3
  a_gamma = sqrt(1. / (1 - pow((_boosted_pip_mc->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_pip_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pip_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pip_mc->Vect().Unit() * _boosted_pim_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_pip_mc->Vect().Unit() * _boosted_pim_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pim_mc->Vect().Unit() + b_beta * _boosted_pip_mc->Vect().Unit();

  alpha_PPIm_piPIp_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pip_mc->Vect() < 0) alpha_PPIm_piPIp_mc = 360. - alpha_PPIm_piPIp_mc;
  //α[pπ−][p'π+]

  _alpha_ppip_pipim_mc = alpha_PPIp_piPIm_mc;
  _alpha_pippim_pipf_mc = alpha_pippim_pipf_mc;
  _alpha_ppim_pipip_mc = alpha_PPIm_piPIp_mc;
}

float MCReaction::MCalpha_ppip_pipim() {  // pipim bhaneko proton initial  pim ho?
  if (_alpha_ppip_pipim_mc != _alpha_ppip_pipim_mc) MCAlphaCalc();
  if (TwoPion())
    return _alpha_ppip_pipim_mc;
  else
    return NAN;
}
float MCReaction::MCalpha_pippim_pipf() {  // alpha P (proton initial proton final)
  if (_alpha_pippim_pipf_mc != _alpha_pippim_pipf_mc) MCAlphaCalc();
  if (TwoPion())
    return _alpha_pippim_pipf_mc;
  else
    return NAN;
}
float MCReaction::MCalpha_ppim_pipip() {  // alpha pip (proton initial pip)
  if (_alpha_ppim_pipip_mc != _alpha_ppim_pipip_mc) MCAlphaCalc();
  if (TwoPion())
    return _alpha_ppim_pipip_mc;
  else
    return NAN;
}

float MCReaction::MCalpha_ppip_pipim_thrown() {  // pipim bhaneko proton initial  pim ho?
  if (_alpha_ppip_pipim_mc != _alpha_ppip_pipim_mc) MCAlphaCalc();

  return _alpha_ppip_pipim_mc;
}
float MCReaction::MCalpha_pippim_pipf_thrown() {  // alpha P (proton initial proton final)
  if (_alpha_pippim_pipf_mc != _alpha_pippim_pipf_mc) MCAlphaCalc();

  return _alpha_pippim_pipf_mc;
}
float MCReaction::MCalpha_ppim_pipip_thrown() {  // alpha pip (proton initial pip)
  if (_alpha_ppim_pipip_mc != _alpha_ppim_pipip_mc) MCAlphaCalc();

  return _alpha_ppim_pipip_mc;
}
float MCReaction::MCgamma_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();
  if (_gamma_mc->Phi() > 0)
    return _gamma_mc->Phi() * 180 / PI;
  else if (_gamma_mc->Phi() < 0)
    return (_gamma_mc->Phi() + 2 * PI) * 180 / PI;
  else
    return NAN;
}
float MCReaction::MCprot_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();

  if (_boosted_prot_mc->Phi() > 0)
    return _boosted_prot_mc->Phi() * 180 / PI;
  else if (_boosted_prot_mc->Phi() < 0)
    return (_boosted_prot_mc->Phi() + 2 * PI) * 180 / PI;
  else
    return NAN;
}
float MCReaction::MCpip_Phi_thrown() {
  if (_boosted_pip_mc->Phi() > 0)
    return _boosted_pip_mc->Phi() * 180 / PI;
  else if (_boosted_pip_mc->Phi() < 0)
    return (_boosted_pip_mc->Phi() + 2 * PI) * 180 / PI;
  else
    return NAN;
}
float MCReaction::MCpim_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();

  if (_boosted_pim_mc->Phi() > 0)
    return _boosted_pim_mc->Phi() * 180 / PI;
  else if (_boosted_pim_mc->Phi() < 0)
    return (_boosted_pim_mc->Phi() + 2 * PI) * 180 / PI;
  else
    return NAN;
}
float MCReaction::MCprot_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_prot_mc->Theta() * 180.0 / PI;
}
float MCReaction::MCpip_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_pip_mc->Theta() * 180.0 / PI;
}
float MCReaction::MCpim_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_pim_mc->Theta() * 180.0 / PI;
}
