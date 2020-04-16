/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include <mutex>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"

using namespace std;

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram {
 protected:
  std::shared_ptr<TFile> RootOutputFile;
  std::shared_ptr<TCanvas> def;

  int bins = 500;
  double p_min = 0.0;
  double p_max = 10.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 12.0;
  double w_max = 5.5;

  double zero = 0.0;

  static const short particle_num = 4;  // 0-e 1-Pi 2-P 3-K
  std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
  static const short charge_num = 2;  // 0-pos 1-neg
  std::string charge_name[charge_num] = {"positive", "negative"};
  static const short with_id_num = 3;  // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};
  static const short num_sectors = 6;
  std::string sec_name[num_sectors] = {"1", "2", "3", "4", "5", "6"};

  static const short CUTS = 2;
  enum cuts { before_cut, after_cut };
  std::mutex mutex;

  static const short w_range_num = 1;
  std::string w_range_name[w_range_num] = {
      "all_W_range" /*" W<2.5 ",
                       " 2.5<W<3.5 "*/
  };                //{" W<2.0 ", " 2.0<W<2.5 ", " 2.5<W<3.0 ", " 3.0<W<3.5 "};
  static const short q2_range_num = 10;
  std::string q2_range_name[q2_range_num] = {" Q2<1.0 ",     " 1.0<Q2<2.0 ", " 2.0<Q2<3.0 ", " 3.0<Q2<4.0 ",
                                             " 4.0<Q2<5.0 ", " 5.0<Q2<6.0 ", " 6.0<Q2<7.0 ", " 7.0<Q2<8.0 ",
                                             " 8.0<Q2<9.0 ", " Q2>9.0 "};  //{" Q2<2.0 ", " 2.0<Q2<5.0 ",
                                                                           //" 5.0<Q2<10.0 "};

  static const short inv_Ppip_range_num = 2;
  std::string inv_Ppip_range_name[inv_Ppip_range_num] = {" M[Ppip]<1.5 GeV ", " M[Ppip]>1.5 GeV "};

  static const short inv_pip_pim_range_num = 2;
  std::string inv_pip_pim_range_name[inv_pip_pim_range_num] = {" M[pi+pi-]<0.9 GeV ", " M[pi+pi-]>0.9 GeV "};

  static const short theta_pim_range_num = 2;
  std::string theta_pim_range_name[theta_pim_range_num] = {" theta[pi-]<90 deg ", " theta[pi-]>90 deg "};

  static const short phi_pim_range_num = 2;
  std::string phi_pim_range_name[phi_pim_range_num] = {" phi[pi-]<90 deg ", " phi[pi-]>90 deg "};
  static const short alpha_pim_range_num = 2;
  std::string alpha_pim_range_name[alpha_pim_range_num] = {" alpha[pi-]<180 deg ", " alpha[pi-]>180 deg "};
  // Kinematics

  THnSparse* sevenDHist;
  THnSparse* h8_1[num_sectors];
  static const short NUM_CUT = 2;

  TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500, 0, 0.5);

  TH1D_ptr E_x_mu = std::make_shared<TH1D>("Energy (x_mu = e(p,p')e')", "Energy (x_mu = e(p,p')e')", 500, -1.0, 1.0);
  TH1D_ptr P_x_mu = std::make_shared<TH1D>("mom (x_mu = e(p,p')e')", "mom (x_mu = e(p,p')e')", 500, -0.50, 1.0);
  TH1D_ptr diff_E_P_x_mu =
      std::make_shared<TH1D>("Energy - mom (x_mu = e(p,p')e')", "Energy - mom (x_mu = e(p,p')e')", 500, -0.50, 1.0);
  TH2D_ptr mom_vs_E_x_mu =
      std::make_shared<TH2D>("Mom_vs_Energy", " Mom_vs_Energy_component", 500, -1.0, 2.0, 500, 0.0, 2.0);
  TH1D_ptr diff_theta_in_x_mu = std::make_shared<TH1D>("#Delata#Theta x_mu and initial electron",
                                                       "#Delata#Theta x_mu and initial electron", 500, 0, 180);
  TH2D_ptr Dthtea_vs_Dphi = std::make_shared<TH2D>("#Delata#Theta vs #Delata#Phi", "#Delata#Theta vs #Delata#Phi", 500,
                                                   -180, 180, 500, -180, 180);

  TH1D_ptr momentum;
  TH1D_ptr W_hist;
  TH1D_ptr Q2_hist;
  TH2D_ptr W_vs_q2;
  TH1D_ptr W_P2pi_hist;

  TH1D_ptr W_hist_ftof;
  TH2D_ptr W_vs_q2_ftof;
  TH1D_ptr W_hist_ctof;
  TH2D_ptr W_vs_q2_ctof;

  TH1D_ptr W_thrown;
  TH2D_ptr W_vs_Q2_thrown;

  TH1D_ptr vz_position[CUTS];
  TH2D_ptr pcal_sec[CUTS];
  TH2D_ptr dcr1_sec[CUTS];
  TH2D_ptr dcr2_sec[CUTS];
  TH2D_ptr dcr3_sec[CUTS];

  TH1D_ptr theta_prot;
  TH1D_ptr theta_pip;
  TH1D_ptr theta_pim;

  TH2D_ptr Theta_vs_mom_x_mu;

  TH1D_ptr Phi_gamma;
  TH1D_ptr Phi_prot;
  TH1D_ptr Phi_pip;
  TH1D_ptr Phi_pim;

  TH1D_ptr alpha_pim;
  TH1D_ptr alpha_pip;
  TH1D_ptr alpha_prot;

  TH1D_ptr theta_prot_mc;
  TH1D_ptr theta_pip_mc;
  TH1D_ptr theta_pim_mc;

  TH1D_ptr Phi_gamma_mc;
  TH1D_ptr Phi_prot_mc;
  TH1D_ptr Phi_pip_mc;
  TH1D_ptr Phi_pim_mc;

  TH1D_ptr alpha_pim_mc;
  TH1D_ptr alpha_pip_mc;
  TH1D_ptr alpha_prot_mc;

  TH1D_ptr theta_prot_thrown;
  TH1D_ptr theta_pip_thrown;
  TH1D_ptr theta_pim_thrown;

  TH1D_ptr Phi_gamma_thrown;
  TH1D_ptr Phi_prot_thrown;
  TH1D_ptr Phi_pip_thrown;
  TH1D_ptr Phi_pim_thrown;

  TH1D_ptr alpha_pim_thrown;
  TH1D_ptr alpha_pip_thrown;
  TH1D_ptr alpha_prot_thrown;

  TH2D_ptr W_vs_q2_sec[num_sectors];
  TH1D_ptr W_sec[num_sectors];

  TH1D_ptr W_det[3];
  TH2D_ptr WQ2_det[3];

  TH1D_ptr W_hist_singlePi;
  TH1D_ptr Q2_hist_singlePi;
  TH2D_ptr W_vs_q2_singlePi;
  TH2D_ptr W_vs_q2_singlePi_sec[num_sectors];
  TH1D_ptr W_singlePi_sec[num_sectors];
  TH2D_ptr W_vs_MM_singlePi[num_sectors];

  TH2D_ptr W_vs_q2_Npip_sec[num_sectors];
  TH1D_ptr W_Npip_sec[num_sectors];
  TH1D_ptr MM_Npip_sec[num_sectors];

  TH1D_ptr MM_neutron;
  TH1D_ptr MM_neutron_sec[num_sectors];

  TH2D_ptr W_vs_q2_twoPi_sec[num_sectors];
  TH1D_ptr W_twoPi_sec[num_sectors];
  TH2D_ptr W_vs_MM_twoPi[num_sectors];
  TH1D_ptr MM_twoPi;
  TH1D_ptr MM_twoPi_sec[num_sectors];
  TH1D_ptr MM2_twoPi;
  TH1D_ptr MM2_twoPi_sec[num_sectors];
  TH1D_ptr W_hist_twoPi;
  TH1D_ptr Q2_hist_twoPi;
  TH2D_ptr W_vs_q2_twoPi;
  
  TH1D_ptr W_q2_15_twoPi;
  TH1D_ptr W_q2_2_twoPi;
  TH1D_ptr W_q2_23_twoPi;
  TH1D_ptr W_q2_34_twoPi; 
  
  TH1D_ptr W_q2_15_photon_twoPi;
  TH1D_ptr W_q2_2_photon_twoPi;
  TH1D_ptr W_q2_23_photon_twoPi;
  TH1D_ptr W_q2_34_photon_twoPi; 
  //
  //
  //
  //modified by Jiawei
  TH1D_ptr invmass_Ppip_twoPi;
  TH1D_ptr invmass_Ppim_twoPi;
  TH1D_ptr invmass_pippim_twoPi;

  TH2D_ptr theta_Ppip_vs_invmass_twoPi;
  TH2D_ptr theta_Ppim_vs_invmass_twoPi;
  TH2D_ptr theta_pippim_vs_invmass_twoPi;
   
  TH2D_ptr p_vs_phi_Prot_twoPi;
  TH2D_ptr p_vs_phi_pip_twoPi;
  TH2D_ptr p_vs_phi_pim_twoPi;
  
  TH2D_ptr p_vs_theta_Prot_twoPi;
  TH2D_ptr p_vs_theta_pip_twoPi;
  TH2D_ptr p_vs_theta_pim_twoPi;

  TH2D_ptr theta_vs_phi_Prot_twoPi;
  TH2D_ptr theta_vs_phi_pip_twoPi;
  TH2D_ptr theta_vs_phi_pim_twoPi;
  //
  //
  //

  TH1D_ptr inv_mass_P_pip[w_range_num][q2_range_num];
  TH1D_ptr inv_mass_P_pim[w_range_num][q2_range_num];
  TH1D_ptr inv_mass_pip_pim[w_range_num][q2_range_num];

  TH2D_ptr theta_P_vs_mass_pip_pim[w_range_num][q2_range_num];
  TH2D_ptr theta_pim_vs_mass_Ppip[w_range_num][q2_range_num];
  TH2D_ptr theta_pip_vs_mass_Ppim[w_range_num][q2_range_num];
  TH2D_ptr theta_P_lab_vs_mass_pip_pim[w_range_num][q2_range_num];
  TH2D_ptr theta_pim_lab_vs_mass_Ppip[w_range_num][q2_range_num];
  TH2D_ptr theta_pip_lab_vs_mass_Ppim[w_range_num][q2_range_num];

  TH2D_ptr W_vs_q2_twoPi_sec_thrown[num_sectors];
  TH1D_ptr W_twoPi_sec_thrown[num_sectors];
  TH2D_ptr W_vs_MM_twoPi_thrown[num_sectors];
  TH1D_ptr MM_twoPi_sec_thrown[num_sectors];
  TH1D_ptr MM2_twoPi_sec_thrown[num_sectors];

  TH1D_ptr inv_mass_P_pip_thrown[w_range_num][q2_range_num];
  TH1D_ptr inv_mass_P_pim_thrown[w_range_num][q2_range_num];
  TH1D_ptr inv_mass_pip_pim_thrown[w_range_num][q2_range_num];

  TH2D_ptr theta_P_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];
  TH2D_ptr theta_pim_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
  TH2D_ptr theta_pip_vs_mass_Ppim_thrown[w_range_num][q2_range_num];
  TH2D_ptr theta_P_lab_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];
  TH2D_ptr theta_pim_lab_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
  TH2D_ptr theta_pip_lab_vs_mass_Ppim_thrown[w_range_num][q2_range_num];

  // EC Sampling Fraction
  TH2D_ptr EC_sampling_fraction[CUTS];
  // EC Sampling Fraction

  // Mom vs Beta
  TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
  // Mom vs Beta

  // Delta T
  TH2D_ptr delta_t_hist[particle_num][charge_num][with_id_num][2];
  // Delta T

 public:
  Histogram(const std::string& output_file);
  ~Histogram();
  void Fill_histSevenD(const std::shared_ptr<Reaction>& _e);
  void writeHists7D();
  // W and Q^2
  void makeHists_sector();
  void Fill_WvsQ2(const std::shared_ptr<Reaction>& _e);
  void Fill_WvsQ2(const std::shared_ptr<MCReaction>& _e);
  void Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction>& _e);
  void Fill_WvsQ2_Npip(const std::shared_ptr<Reaction>& _e);
  void Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction>& _e);
  // void Fill_WvsQ2_twoPi(const std::shared_ptr<MCReaction>& _e);
  void Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<MCReaction>& _e);
  void Write_WvsQ2();

  void Fill_W_vs_Q2_all_sec();
  void Fill_W_vs_Q2_thrown();
  void Fill_inv_mass_hist();
  // P and E
  // ecectron cuts
  void makeHists_electron_cuts();
  void FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d);
  void FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d);

  void Write_Electron_cuts();

  // sampling Fraction
  void Fill_SF(const std::shared_ptr<Branches12>& _d);
  void Write_SF();
  void Fill_x_mu(const std::shared_ptr<Reaction>& _e);

  void makeHists_MomVsBeta();
  void Fill_MomVsBeta(const std::shared_ptr<Branches12>& data, int part);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat_pi(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt, int part);
  void Fill_deltat_prot(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt, int part);

  void Write_deltat();

  // EC Sampling Fraction
  //  void Fill_EC(double etot, double momentum);
  // void Write_EC();

  //
  void Write();
};

#endif
