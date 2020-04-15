/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "branches.hpp"

Branches12::Branches12(const std::shared_ptr<TChain> &tree) {
  _tree = tree;
  Branches12::init();
}

Branches12::Branches12(const std::shared_ptr<TChain> &tree, bool mc) {
  _tree = tree;
  _is_mc = mc;
  Branches12::init();
  if (_is_mc) Branches12::initMC();
}

bool Branches12::mc() { return _is_mc; }

void Branches12::mc_branches() {
  _is_mc = true;
  Branches12::initMC();
}

void Branches12::init() {
  _pid = 0;
  _p = 0;
  _p2 = 0;
  _px = 0;
  _py = 0;
  _pz = 0;
  _vx = 0;
  _vy = 0;
  _vz = 0;
  _charge = 0;
  _beta = 0;
  _chi2pid = 0;
  _status = 0;

  _dc_sec = 0;
  _dc_r1_x = 0;
  _dc_r1_y = 0;
  _dc_r1_z = 0;
  _dc_r2_x = 0;
  _dc_r2_y = 0;
  _dc_r2_z = 0;
  _dc_r3_x = 0;
  _dc_r3_y = 0;
  _dc_r3_z = 0;
  _cvt_x = 0;
  _cvt_y = 0;
  _cvt_z = 0;
  _fmt_x = 0;
  _fmt_y = 0;
  _fmt_z = 0;
  _ec_tot_energy = 0;
  _ec_pcal_energy = 0;
  _ec_pcal_sec = 0;
  _ec_pcal_time = 0;
  _ec_pcal_path = 0;
  _ec_pcal_x = 0;
  _ec_pcal_y = 0;
  _ec_pcal_z = 0;
  _ec_pcal_hx = 0;
  _ec_pcal_hy = 0;
  _ec_pcal_hz = 0;
  _ec_pcal_lu = 0;
  _ec_pcal_lv = 0;
  _ec_pcal_lw = 0;
  _ec_pcal_du = 0;
  _ec_pcal_dv = 0;
  _ec_pcal_dw = 0;
  _ec_pcal_m2u = 0;
  _ec_pcal_m2v = 0;
  _ec_pcal_m2w = 0;
  _ec_pcal_m3u = 0;
  _ec_pcal_m3v = 0;
  _ec_pcal_m3w = 0;
  _ec_ecin_energy = 0;
  _ec_ecin_sec = 0;
  _ec_ecin_time = 0;
  _ec_ecin_path = 0;
  _ec_ecin_x = 0;
  _ec_ecin_y = 0;
  _ec_ecin_z = 0;
  _ec_ecin_hx = 0;
  _ec_ecin_hy = 0;
  _ec_ecin_hz = 0;
  _ec_ecin_lu = 0;
  _ec_ecin_lv = 0;
  _ec_ecin_lw = 0;
  _ec_ecin_du = 0;
  _ec_ecin_dv = 0;
  _ec_ecin_dw = 0;
  _ec_ecin_m2u = 0;
  _ec_ecin_m2v = 0;
  _ec_ecin_m2w = 0;
  _ec_ecin_m3u = 0;
  _ec_ecin_m3v = 0;
  _ec_ecin_m3w = 0;
  _ec_ecout_energy = 0;
  _ec_ecout_sec = 0;
  _ec_ecout_time = 0;
  _ec_ecout_path = 0;
  _ec_ecout_x = 0;
  _ec_ecout_y = 0;
  _ec_ecout_z = 0;
  _ec_ecout_hx = 0;
  _ec_ecout_hy = 0;
  _ec_ecout_hz = 0;
  _ec_ecout_lu = 0;
  _ec_ecout_lv = 0;
  _ec_ecout_lw = 0;
  _ec_ecout_du = 0;
  _ec_ecout_dv = 0;
  _ec_ecout_dw = 0;
  _ec_ecout_m2u = 0;
  _ec_ecout_m2v = 0;
  _ec_ecout_m2w = 0;
  _ec_ecout_m3u = 0;
  _ec_ecout_m3v = 0;
  _ec_ecout_m3w = 0;
  _cc_nphe_tot = 0;
  _cc_ltcc_sec = 0;
  _cc_ltcc_nphe = 0;
  _cc_ltcc_time = 0;
  _cc_ltcc_path = 0;
  _cc_ltcc_theta = 0;
  _cc_ltcc_phi = 0;
  _cc_ltcc_x = 0;
  _cc_ltcc_y = 0;
  _cc_ltcc_z = 0;
  _cc_htcc_sec = 0;
  _cc_htcc_nphe = 0;
  _cc_htcc_time = 0;
  _cc_htcc_path = 0;
  _cc_htcc_theta = 0;
  _cc_htcc_phi = 0;
  _cc_htcc_x = 0;
  _cc_htcc_y = 0;
  _cc_htcc_z = 0;
  _cc_rich_sec = 0;
  _cc_rich_nphe = 0;
  _cc_rich_time = 0;
  _cc_rich_path = 0;
  _cc_rich_theta = 0;
  _cc_rich_phi = 0;
  _cc_rich_x = 0;
  _cc_rich_y = 0;
  _cc_rich_z = 0;
  _sc_ftof_1a_sec = 0;
  _sc_ftof_1a_time = 0;
  _sc_ftof_1a_path = 0;
  _sc_ftof_1a_energy = 0;
  _sc_ftof_1a_component = 0;
  _sc_ftof_1a_x = 0;
  _sc_ftof_1a_y = 0;
  _sc_ftof_1a_z = 0;
  _sc_ftof_1a_hx = 0;
  _sc_ftof_1a_hy = 0;
  _sc_ftof_1a_hz = 0;
  _sc_ftof_1b_sec = 0;
  _sc_ftof_1b_time = 0;
  _sc_ftof_1b_path = 0;
  _sc_ftof_1b_energy = 0;
  _sc_ftof_1b_component = 0;
  _sc_ftof_1b_x = 0;
  _sc_ftof_1b_y = 0;
  _sc_ftof_1b_z = 0;
  _sc_ftof_1b_hx = 0;
  _sc_ftof_1b_hy = 0;
  _sc_ftof_1b_hz = 0;
  _sc_ftof_2_sec = 0;
  _sc_ftof_2_time = 0;
  _sc_ftof_2_path = 0;
  _sc_ftof_2_energy = 0;
  _sc_ftof_2_component = 0;
  _sc_ftof_2_x = 0;
  _sc_ftof_2_y = 0;
  _sc_ftof_2_z = 0;
  _sc_ftof_2_hx = 0;
  _sc_ftof_2_hy = 0;
  _sc_ftof_2_hz = 0;
  _sc_ctof_time = 0;
  _sc_ctof_path = 0;
  _sc_ctof_energy = 0;
  _sc_ctof_component = 0;
  _sc_ctof_x = 0;
  _sc_ctof_y = 0;
  _sc_ctof_z = 0;
  _sc_ctof_hx = 0;
  _sc_ctof_hy = 0;
  _sc_ctof_hz = 0;
  _sc_cnd_time = 0;
  _sc_cnd_path = 0;
  _sc_cnd_energy = 0;
  _sc_cnd_component = 0;
  _sc_cnd_x = 0;
  _sc_cnd_y = 0;
  _sc_cnd_z = 0;
  _sc_cnd_hx = 0;
  _sc_cnd_hy = 0;
  _sc_cnd_hz = 0;
  _ft_cal_energy = 0;
  _ft_cal_time = 0;
  _ft_cal_path = 0;
  _ft_cal_x = 0;
  _ft_cal_y = 0;
  _ft_cal_z = 0;
  _ft_cal_dx = 0;
  _ft_cal_dy = 0;
  _ft_cal_radius = 0;
  _ft_hodo_energy = 0;
  _ft_hodo_time = 0;
  _ft_hodo_path = 0;
  _ft_hodo_x = 0;
  _ft_hodo_y = 0;
  _ft_hodo_z = 0;
  _ft_hodo_dx = 0;
  _ft_hodo_dy = 0;
  _ft_hodo_radius = 0;

  _tree->SetBranchStatus("*", 0);
  _tree->SetBranchAddress("pid", &_pid);
  _tree->SetBranchAddress("p", &_p);
  _tree->SetBranchAddress("p2", &_p2);
  _tree->SetBranchAddress("px", &_px);
  _tree->SetBranchAddress("py", &_py);
  _tree->SetBranchAddress("pz", &_pz);
  _tree->SetBranchAddress("vx", &_vx);
  _tree->SetBranchAddress("vy", &_vy);
  _tree->SetBranchAddress("vz", &_vz);
  _tree->SetBranchAddress("charge", &_charge);
  _tree->SetBranchAddress("beta", &_beta);
  _tree->SetBranchAddress("chi2pid", &_chi2pid);
  _tree->SetBranchAddress("status", &_status);

  _tree->SetBranchAddress("dc_sec", &_dc_sec);
  _tree->SetBranchAddress("dc_r1_x", &_dc_r1_x);
  _tree->SetBranchAddress("dc_r1_y", &_dc_r1_y);
  _tree->SetBranchAddress("dc_r1_z", &_dc_r1_z);
  _tree->SetBranchAddress("dc_r2_x", &_dc_r2_x);
  _tree->SetBranchAddress("dc_r2_y", &_dc_r2_y);
  _tree->SetBranchAddress("dc_r2_z", &_dc_r2_z);
  _tree->SetBranchAddress("dc_r3_x", &_dc_r3_x);
  _tree->SetBranchAddress("dc_r3_y", &_dc_r3_y);
  _tree->SetBranchAddress("dc_r3_z", &_dc_r3_z);
  _tree->SetBranchAddress("cvt_x", &_cvt_x);
  _tree->SetBranchAddress("cvt_y", &_cvt_y);
  _tree->SetBranchAddress("cvt_z", &_cvt_z);
  _tree->SetBranchAddress("fmt_x", &_fmt_x);
  _tree->SetBranchAddress("fmt_y", &_fmt_y);
  _tree->SetBranchAddress("fmt_z", &_fmt_z);
  _tree->SetBranchAddress("ec_tot_energy", &_ec_tot_energy);
  _tree->SetBranchAddress("ec_pcal_energy", &_ec_pcal_energy);
  _tree->SetBranchAddress("ec_pcal_sec", &_ec_pcal_sec);
  _tree->SetBranchAddress("ec_pcal_time", &_ec_pcal_time);
  _tree->SetBranchAddress("ec_pcal_path", &_ec_pcal_path);
  _tree->SetBranchAddress("ec_pcal_x", &_ec_pcal_x);
  _tree->SetBranchAddress("ec_pcal_y", &_ec_pcal_y);
  _tree->SetBranchAddress("ec_pcal_z", &_ec_pcal_z);
  _tree->SetBranchAddress("ec_pcal_hx", &_ec_pcal_hx);
  _tree->SetBranchAddress("ec_pcal_hy", &_ec_pcal_hy);
  _tree->SetBranchAddress("ec_pcal_hz", &_ec_pcal_hz);
  _tree->SetBranchAddress("ec_pcal_lu", &_ec_pcal_lu);
  _tree->SetBranchAddress("ec_pcal_lv", &_ec_pcal_lv);
  _tree->SetBranchAddress("ec_pcal_lw", &_ec_pcal_lw);
  _tree->SetBranchAddress("ec_pcal_du", &_ec_pcal_du);
  _tree->SetBranchAddress("ec_pcal_dv", &_ec_pcal_dv);
  _tree->SetBranchAddress("ec_pcal_dw", &_ec_pcal_dw);
  _tree->SetBranchAddress("ec_pcal_m2u", &_ec_pcal_m2u);
  _tree->SetBranchAddress("ec_pcal_m2v", &_ec_pcal_m2v);
  _tree->SetBranchAddress("ec_pcal_m2w", &_ec_pcal_m2w);
  _tree->SetBranchAddress("ec_pcal_m3u", &_ec_pcal_m3u);
  _tree->SetBranchAddress("ec_pcal_m3v", &_ec_pcal_m3v);
  _tree->SetBranchAddress("ec_pcal_m3w", &_ec_pcal_m3w);
  _tree->SetBranchAddress("ec_ecin_energy", &_ec_ecin_energy);
  _tree->SetBranchAddress("ec_ecin_sec", &_ec_ecin_sec);
  _tree->SetBranchAddress("ec_ecin_time", &_ec_ecin_time);
  _tree->SetBranchAddress("ec_ecin_path", &_ec_ecin_path);
  _tree->SetBranchAddress("ec_ecin_x", &_ec_ecin_x);
  _tree->SetBranchAddress("ec_ecin_y", &_ec_ecin_y);
  _tree->SetBranchAddress("ec_ecin_z", &_ec_ecin_z);
  _tree->SetBranchAddress("ec_ecin_hx", &_ec_ecin_hx);
  _tree->SetBranchAddress("ec_ecin_hy", &_ec_ecin_hy);
  _tree->SetBranchAddress("ec_ecin_hz", &_ec_ecin_hz);
  _tree->SetBranchAddress("ec_ecin_lu", &_ec_ecin_lu);
  _tree->SetBranchAddress("ec_ecin_lv", &_ec_ecin_lv);
  _tree->SetBranchAddress("ec_ecin_lw", &_ec_ecin_lw);
  _tree->SetBranchAddress("ec_ecin_du", &_ec_ecin_du);
  _tree->SetBranchAddress("ec_ecin_dv", &_ec_ecin_dv);
  _tree->SetBranchAddress("ec_ecin_dw", &_ec_ecin_dw);
  _tree->SetBranchAddress("ec_ecin_m2u", &_ec_ecin_m2u);
  _tree->SetBranchAddress("ec_ecin_m2v", &_ec_ecin_m2v);
  _tree->SetBranchAddress("ec_ecin_m2w", &_ec_ecin_m2w);
  _tree->SetBranchAddress("ec_ecin_m3u", &_ec_ecin_m3u);
  _tree->SetBranchAddress("ec_ecin_m3v", &_ec_ecin_m3v);
  _tree->SetBranchAddress("ec_ecin_m3w", &_ec_ecin_m3w);
  _tree->SetBranchAddress("ec_ecout_energy", &_ec_ecout_energy);
  _tree->SetBranchAddress("ec_ecout_sec", &_ec_ecout_sec);
  _tree->SetBranchAddress("ec_ecout_time", &_ec_ecout_time);
  _tree->SetBranchAddress("ec_ecout_path", &_ec_ecout_path);
  _tree->SetBranchAddress("ec_ecout_x", &_ec_ecout_x);
  _tree->SetBranchAddress("ec_ecout_y", &_ec_ecout_y);
  _tree->SetBranchAddress("ec_ecout_z", &_ec_ecout_z);
  _tree->SetBranchAddress("ec_ecout_hx", &_ec_ecout_hx);
  _tree->SetBranchAddress("ec_ecout_hy", &_ec_ecout_hy);
  _tree->SetBranchAddress("ec_ecout_hz", &_ec_ecout_hz);
  _tree->SetBranchAddress("ec_ecout_lu", &_ec_ecout_lu);
  _tree->SetBranchAddress("ec_ecout_lv", &_ec_ecout_lv);
  _tree->SetBranchAddress("ec_ecout_lw", &_ec_ecout_lw);
  _tree->SetBranchAddress("ec_ecout_du", &_ec_ecout_du);
  _tree->SetBranchAddress("ec_ecout_dv", &_ec_ecout_dv);
  _tree->SetBranchAddress("ec_ecout_dw", &_ec_ecout_dw);
  _tree->SetBranchAddress("ec_ecout_m2u", &_ec_ecout_m2u);
  _tree->SetBranchAddress("ec_ecout_m2v", &_ec_ecout_m2v);
  _tree->SetBranchAddress("ec_ecout_m2w", &_ec_ecout_m2w);
  _tree->SetBranchAddress("ec_ecout_m3u", &_ec_ecout_m3u);
  _tree->SetBranchAddress("ec_ecout_m3v", &_ec_ecout_m3v);
  _tree->SetBranchAddress("ec_ecout_m3w", &_ec_ecout_m3w);
  _tree->SetBranchAddress("cc_nphe_tot", &_cc_nphe_tot);
  _tree->SetBranchAddress("cc_ltcc_sec", &_cc_ltcc_sec);
  _tree->SetBranchAddress("cc_ltcc_nphe", &_cc_ltcc_nphe);
  _tree->SetBranchAddress("cc_ltcc_time", &_cc_ltcc_time);
  _tree->SetBranchAddress("cc_ltcc_path", &_cc_ltcc_path);
  _tree->SetBranchAddress("cc_ltcc_theta", &_cc_ltcc_theta);
  _tree->SetBranchAddress("cc_ltcc_phi", &_cc_ltcc_phi);
  _tree->SetBranchAddress("cc_ltcc_x", &_cc_ltcc_x);
  _tree->SetBranchAddress("cc_ltcc_y", &_cc_ltcc_y);
  _tree->SetBranchAddress("cc_ltcc_z", &_cc_ltcc_z);
  _tree->SetBranchAddress("cc_htcc_sec", &_cc_htcc_sec);
  _tree->SetBranchAddress("cc_htcc_nphe", &_cc_htcc_nphe);
  _tree->SetBranchAddress("cc_htcc_time", &_cc_htcc_time);
  _tree->SetBranchAddress("cc_htcc_path", &_cc_htcc_path);
  _tree->SetBranchAddress("cc_htcc_theta", &_cc_htcc_theta);
  _tree->SetBranchAddress("cc_htcc_phi", &_cc_htcc_phi);
  _tree->SetBranchAddress("cc_htcc_x", &_cc_htcc_x);
  _tree->SetBranchAddress("cc_htcc_y", &_cc_htcc_y);
  _tree->SetBranchAddress("cc_htcc_z", &_cc_htcc_z);
  _tree->SetBranchAddress("cc_rich_sec", &_cc_rich_sec);
  _tree->SetBranchAddress("cc_rich_nphe", &_cc_rich_nphe);
  _tree->SetBranchAddress("cc_rich_time", &_cc_rich_time);
  _tree->SetBranchAddress("cc_rich_path", &_cc_rich_path);
  _tree->SetBranchAddress("cc_rich_theta", &_cc_rich_theta);
  _tree->SetBranchAddress("cc_rich_phi", &_cc_rich_phi);
  _tree->SetBranchAddress("cc_rich_x", &_cc_rich_x);
  _tree->SetBranchAddress("cc_rich_y", &_cc_rich_y);
  _tree->SetBranchAddress("cc_rich_z", &_cc_rich_z);
  _tree->SetBranchAddress("sc_ftof_1a_sec", &_sc_ftof_1a_sec);
  _tree->SetBranchAddress("sc_ftof_1a_time", &_sc_ftof_1a_time);
  _tree->SetBranchAddress("sc_ftof_1a_path", &_sc_ftof_1a_path);
  _tree->SetBranchAddress("sc_ftof_1a_energy", &_sc_ftof_1a_energy);
  _tree->SetBranchAddress("sc_ftof_1a_component", &_sc_ftof_1a_component);
  _tree->SetBranchAddress("sc_ftof_1a_x", &_sc_ftof_1a_x);
  _tree->SetBranchAddress("sc_ftof_1a_y", &_sc_ftof_1a_y);
  _tree->SetBranchAddress("sc_ftof_1a_z", &_sc_ftof_1a_z);
  _tree->SetBranchAddress("sc_ftof_1a_hx", &_sc_ftof_1a_hx);
  _tree->SetBranchAddress("sc_ftof_1a_hy", &_sc_ftof_1a_hy);
  _tree->SetBranchAddress("sc_ftof_1a_hz", &_sc_ftof_1a_hz);
  _tree->SetBranchAddress("sc_ftof_1b_sec", &_sc_ftof_1b_sec);
  _tree->SetBranchAddress("sc_ftof_1b_time", &_sc_ftof_1b_time);
  _tree->SetBranchAddress("sc_ftof_1b_path", &_sc_ftof_1b_path);
  _tree->SetBranchAddress("sc_ftof_1b_energy", &_sc_ftof_1b_energy);
  _tree->SetBranchAddress("sc_ftof_1b_component", &_sc_ftof_1b_component);
  _tree->SetBranchAddress("sc_ftof_1b_x", &_sc_ftof_1b_x);
  _tree->SetBranchAddress("sc_ftof_1b_y", &_sc_ftof_1b_y);
  _tree->SetBranchAddress("sc_ftof_1b_z", &_sc_ftof_1b_z);
  _tree->SetBranchAddress("sc_ftof_1b_hx", &_sc_ftof_1b_hx);
  _tree->SetBranchAddress("sc_ftof_1b_hy", &_sc_ftof_1b_hy);
  _tree->SetBranchAddress("sc_ftof_1b_hz", &_sc_ftof_1b_hz);
  _tree->SetBranchAddress("sc_ftof_2_sec", &_sc_ftof_2_sec);
  _tree->SetBranchAddress("sc_ftof_2_time", &_sc_ftof_2_time);
  _tree->SetBranchAddress("sc_ftof_2_path", &_sc_ftof_2_path);
  _tree->SetBranchAddress("sc_ftof_2_energy", &_sc_ftof_2_energy);
  _tree->SetBranchAddress("sc_ftof_2_component", &_sc_ftof_2_component);
  _tree->SetBranchAddress("sc_ftof_2_x", &_sc_ftof_2_x);
  _tree->SetBranchAddress("sc_ftof_2_y", &_sc_ftof_2_y);
  _tree->SetBranchAddress("sc_ftof_2_z", &_sc_ftof_2_z);
  _tree->SetBranchAddress("sc_ftof_2_hx", &_sc_ftof_2_hx);
  _tree->SetBranchAddress("sc_ftof_2_hy", &_sc_ftof_2_hy);
  _tree->SetBranchAddress("sc_ftof_2_hz", &_sc_ftof_2_hz);
  _tree->SetBranchAddress("sc_ctof_time", &_sc_ctof_time);
  _tree->SetBranchAddress("sc_ctof_path", &_sc_ctof_path);
  _tree->SetBranchAddress("sc_ctof_energy", &_sc_ctof_energy);
  _tree->SetBranchAddress("sc_ctof_component", &_sc_ctof_component);
  _tree->SetBranchAddress("sc_ctof_x", &_sc_ctof_x);
  _tree->SetBranchAddress("sc_ctof_y", &_sc_ctof_y);
  _tree->SetBranchAddress("sc_ctof_z", &_sc_ctof_z);
  _tree->SetBranchAddress("sc_ctof_hx", &_sc_ctof_hx);
  _tree->SetBranchAddress("sc_ctof_hy", &_sc_ctof_hy);
  _tree->SetBranchAddress("sc_ctof_hz", &_sc_ctof_hz);
  _tree->SetBranchAddress("sc_cnd_time", &_sc_cnd_time);
  _tree->SetBranchAddress("sc_cnd_path", &_sc_cnd_path);
  _tree->SetBranchAddress("sc_cnd_energy", &_sc_cnd_energy);
  _tree->SetBranchAddress("sc_cnd_component", &_sc_cnd_component);
  _tree->SetBranchAddress("sc_cnd_x", &_sc_cnd_x);
  _tree->SetBranchAddress("sc_cnd_y", &_sc_cnd_y);
  _tree->SetBranchAddress("sc_cnd_z", &_sc_cnd_z);
  _tree->SetBranchAddress("sc_cnd_hx", &_sc_cnd_hx);
  _tree->SetBranchAddress("sc_cnd_hy", &_sc_cnd_hy);
  _tree->SetBranchAddress("sc_cnd_hz", &_sc_cnd_hz);
  _tree->SetBranchAddress("ft_cal_energy", &_ft_cal_energy);
  _tree->SetBranchAddress("ft_cal_time", &_ft_cal_time);
  _tree->SetBranchAddress("ft_cal_path", &_ft_cal_path);
  _tree->SetBranchAddress("ft_cal_x", &_ft_cal_x);
  _tree->SetBranchAddress("ft_cal_y", &_ft_cal_y);
  _tree->SetBranchAddress("ft_cal_z", &_ft_cal_z);
  _tree->SetBranchAddress("ft_cal_dx", &_ft_cal_dx);
  _tree->SetBranchAddress("ft_cal_dy", &_ft_cal_dy);
  _tree->SetBranchAddress("ft_cal_radius", &_ft_cal_radius);
  _tree->SetBranchAddress("ft_hodo_energy", &_ft_hodo_energy);
  _tree->SetBranchAddress("ft_hodo_time", &_ft_hodo_time);
  _tree->SetBranchAddress("ft_hodo_path", &_ft_hodo_path);
  _tree->SetBranchAddress("ft_hodo_x", &_ft_hodo_x);
  _tree->SetBranchAddress("ft_hodo_y", &_ft_hodo_y);
  _tree->SetBranchAddress("ft_hodo_z", &_ft_hodo_z);
  _tree->SetBranchAddress("ft_hodo_dx", &_ft_hodo_dx);
  _tree->SetBranchAddress("ft_hodo_dy", &_ft_hodo_dy);
  _tree->SetBranchAddress("ft_hodo_radius", &_ft_hodo_radius);
}

void Branches12::initMC() {
  _mc_run = 0;
  _mc_event = 0;
  _mc_type = 0;
  _mc_helicity = 0;
  _mc_weight = 0;
  _mc_npart = 0;
  _mc_pid = 0;
  _mc_helicity = 0;
  _mc_px = 0;
  _mc_py = 0;
  _mc_pz = 0;
  _mc_vx = 0;
  _mc_vy = 0;
  _mc_vz = 0;
  _mc_vt = 0;

  _tree->SetBranchAddress("mc_helicity", &_mc_helicity);
  _tree->SetBranchAddress("mc_weight", &_mc_weight);
  _tree->SetBranchAddress("mc_npart", &_mc_npart);
  _tree->SetBranchAddress("mc_pid", &_mc_pid);
  _tree->SetBranchAddress("mc_px", &_mc_px);
  _tree->SetBranchAddress("mc_py", &_mc_py);
  _tree->SetBranchAddress("mc_pz", &_mc_pz);
  _tree->SetBranchAddress("mc_vx", &_mc_vx);
  _tree->SetBranchAddress("mc_vy", &_mc_vy);
  _tree->SetBranchAddress("mc_vz", &_mc_vz);
  _tree->SetBranchAddress("mc_vt", &_mc_vt);
}

int Branches12::gpart() { return _pid->size(); }
int Branches12::pid(int i) {
  if (i >= _pid->size())
    return -9999;
  else
    return _pid->at(i);
}
float Branches12::p(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _p->at(i);
}
float Branches12::p2(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _p2->at(i);
}
float Branches12::px(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _px->at(i);
}
float Branches12::py(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _py->at(i);
}
float Branches12::pz(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _pz->at(i);
}
float Branches12::vx(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _vx->at(i);
}
float Branches12::vy(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _vy->at(i);
}
float Branches12::vz(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _vz->at(i);
}
int Branches12::charge(int i) {
  if (i >= _pid->size())
    return -9999;
  else
    return _charge->at(i);
}
float Branches12::beta(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _beta->at(i);
}

float Branches12::chi2pid(int i) {
  if (i >= _pid->size())
    return NAN;
  else
    return _chi2pid->at(i);
}

int Branches12::status(int i) {
  if (i >= _pid->size())
    return -9999;
  else
    return _status->at(i);
}

int Branches12::dc_sec(int i) {
  if (i >= _dc_sec->size())
    return -9999;
  else
    return _dc_sec->at(i);
}
float Branches12::dc_r1_x(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r1_x->at(i);
}
float Branches12::dc_r1_y(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r1_y->at(i);
}
float Branches12::dc_r1_z(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r1_z->at(i);
}
float Branches12::dc_r2_x(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r2_x->at(i);
}
float Branches12::dc_r2_y(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r2_y->at(i);
}
float Branches12::dc_r2_z(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r2_z->at(i);
}
float Branches12::dc_r3_x(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r3_x->at(i);
}
float Branches12::dc_r3_y(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r3_y->at(i);
}
float Branches12::dc_r3_z(int i) {
  if (i >= _dc_sec->size())
    return NAN;
  else
    return _dc_r3_z->at(i);
}

int Branches12::sc_ftof_1a_sec(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_1a_sec->at(i);
}
float Branches12::sc_ftof_1a_time(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_time->at(i);
}
float Branches12::sc_ftof_1a_path(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_path->at(i);
}
float Branches12::sc_ftof_1a_energy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_energy->at(i);
}
int Branches12::sc_ftof_1a_component(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_1a_component->at(i);
}
float Branches12::sc_ftof_1a_x(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_x->at(i);
}
float Branches12::sc_ftof_1a_y(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_y->at(i);
}
float Branches12::sc_ftof_1a_z(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_z->at(i);
}
float Branches12::sc_ftof_1a_hx(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_hx->at(i);
}
float Branches12::sc_ftof_1a_hy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_hy->at(i);
}
float Branches12::sc_ftof_1a_hz(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1a_hz->at(i);
}
int Branches12::sc_ftof_1b_sec(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_1b_sec->at(i);
}
float Branches12::sc_ftof_1b_time(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_time->at(i);
}
float Branches12::sc_ftof_1b_path(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_path->at(i);
}
float Branches12::sc_ftof_1b_energy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_energy->at(i);
}
int Branches12::sc_ftof_1b_component(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_1b_component->at(i);
}
float Branches12::sc_ftof_1b_x(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_x->at(i);
}
float Branches12::sc_ftof_1b_y(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_y->at(i);
}
float Branches12::sc_ftof_1b_z(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_z->at(i);
}
float Branches12::sc_ftof_1b_hx(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_hx->at(i);
}
float Branches12::sc_ftof_1b_hy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_hy->at(i);
}
float Branches12::sc_ftof_1b_hz(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_1b_hz->at(i);
}
int Branches12::sc_ftof_2_sec(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_2_sec->at(i);
}
float Branches12::sc_ftof_2_time(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_time->at(i);
}
float Branches12::sc_ftof_2_path(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_path->at(i);
}
float Branches12::sc_ftof_2_energy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_energy->at(i);
}
int Branches12::sc_ftof_2_component(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ftof_2_component->at(i);
}
float Branches12::sc_ftof_2_x(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_x->at(i);
}
float Branches12::sc_ftof_2_y(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_y->at(i);
}
float Branches12::sc_ftof_2_z(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_z->at(i);
}
float Branches12::sc_ftof_2_hx(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_hx->at(i);
}
float Branches12::sc_ftof_2_hy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_hy->at(i);
}
float Branches12::sc_ftof_2_hz(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ftof_2_hz->at(i);
}
float Branches12::sc_ctof_time(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ctof_time->at(i);
}
float Branches12::sc_ctof_path(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ctof_path->at(i);
}
float Branches12::sc_ctof_energy(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return NAN;
  else
    return _sc_ctof_energy->at(i);
}
int Branches12::sc_ctof_component(int i) {
  if (i >= _sc_ftof_1b_sec->size())
    return -9999;
  else
    return _sc_ctof_component->at(i);
}

int Branches12::mc_run() { return _mc_run; }
int Branches12::mc_event() { return _mc_event; }
int Branches12::mc_type() { return _mc_type; }
int Branches12::mc_helicity() { return _mc_helicity; }
float Branches12::mc_weight() { return _mc_weight; }
int Branches12::mc_npart() { return _mc_npart; }

int Branches12::mc_pid(int i) {
  if (i >= _mc_npart)
    return -9999;
  else
    return _mc_pid->at(i);
}
float Branches12::mc_px(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_px->at(i);
}
float Branches12::mc_py(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_py->at(i);
}
float Branches12::mc_pz(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_pz->at(i);
}
float Branches12::mc_vx(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_vx->at(i);
}
float Branches12::mc_vy(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_vy->at(i);
}
float Branches12::mc_vz(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_vz->at(i);
}
float Branches12::mc_vt(int i) {
  if (i >= _mc_npart)
    return NAN;
  else
    return _mc_vt->at(i);
}

float Branches12::ec_tot_energy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_tot_energy->at(i);
}
float Branches12::ec_pcal_energy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_energy->at(i);
}
int Branches12::ec_pcal_sec(int i) {
  if (i >= _ec_ecin_sec->size())
    return -9999;
  else
    return _ec_pcal_sec->at(i);
};
float Branches12::ec_pcal_time(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_time->at(i);
}
float Branches12::ec_pcal_path(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_path->at(i);
}
float Branches12::ec_pcal_x(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_x->at(i);
}
float Branches12::ec_pcal_y(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_y->at(i);
}
float Branches12::ec_pcal_z(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_z->at(i);
}
float Branches12::ec_pcal_hx(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_hx->at(i);
}
float Branches12::ec_pcal_hy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_hy->at(i);
}
float Branches12::ec_pcal_hz(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_hz->at(i);
}
float Branches12::ec_pcal_lu(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_lu->at(i);
}
float Branches12::ec_pcal_lv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_lv->at(i);
}
float Branches12::ec_pcal_lw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_lw->at(i);
}
float Branches12::ec_pcal_du(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_du->at(i);
}
float Branches12::ec_pcal_dv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_dv->at(i);
}
float Branches12::ec_pcal_dw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_dw->at(i);
}
float Branches12::ec_pcal_m2u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m2u->at(i);
}
float Branches12::ec_pcal_m2v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m2v->at(i);
}
float Branches12::ec_pcal_m2w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m2w->at(i);
}
float Branches12::ec_pcal_m3u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m3u->at(i);
}
float Branches12::ec_pcal_m3v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m3v->at(i);
}
float Branches12::ec_pcal_m3w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_pcal_m3w->at(i);
}

float Branches12::ec_ecin_energy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_energy->at(i);
}
int Branches12::ec_ecin_sec(int i) {
  if (i >= _ec_ecin_sec->size())
    return -9999;
  else
    return _ec_ecin_sec->at(i);
}
float Branches12::ec_ecin_time(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_time->at(i);
}
float Branches12::ec_ecin_path(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_path->at(i);
}
float Branches12::ec_ecin_x(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_x->at(i);
}
float Branches12::ec_ecin_y(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_y->at(i);
}
float Branches12::ec_ecin_z(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_z->at(i);
}
float Branches12::ec_ecin_hx(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_hx->at(i);
}
float Branches12::ec_ecin_hy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_hy->at(i);
}
float Branches12::ec_ecin_hz(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_hz->at(i);
}
float Branches12::ec_ecin_lu(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_lu->at(i);
}
float Branches12::ec_ecin_lv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_lv->at(i);
}
float Branches12::ec_ecin_lw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_lw->at(i);
}
float Branches12::ec_ecin_du(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_du->at(i);
}
float Branches12::ec_ecin_dv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_dv->at(i);
}
float Branches12::ec_ecin_dw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_dw->at(i);
}
float Branches12::ec_ecin_m2u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m2u->at(i);
}
float Branches12::ec_ecin_m2v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m2v->at(i);
}
float Branches12::ec_ecin_m2w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m2w->at(i);
}
float Branches12::ec_ecin_m3u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m3u->at(i);
}
float Branches12::ec_ecin_m3v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m3v->at(i);
}
float Branches12::ec_ecin_m3w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecin_m3w->at(i);
}
float Branches12::ec_ecout_energy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_energy->at(i);
}
int Branches12::ec_ecout_sec(int i) {
  if (i >= _ec_ecin_sec->size())
    return -9999;
  else
    return _ec_ecout_sec->at(i);
}
float Branches12::ec_ecout_time(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_time->at(i);
}
float Branches12::ec_ecout_path(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_path->at(i);
}
float Branches12::ec_ecout_x(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_x->at(i);
}
float Branches12::ec_ecout_y(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_y->at(i);
}
float Branches12::ec_ecout_z(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_z->at(i);
}
float Branches12::ec_ecout_hx(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_hx->at(i);
}
float Branches12::ec_ecout_hy(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_hy->at(i);
}
float Branches12::ec_ecout_hz(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_hz->at(i);
}
float Branches12::ec_ecout_lu(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_lu->at(i);
}
float Branches12::ec_ecout_lv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_lv->at(i);
}
float Branches12::ec_ecout_lw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_lw->at(i);
}
float Branches12::ec_ecout_du(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_du->at(i);
}
float Branches12::ec_ecout_dv(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_dv->at(i);
}
float Branches12::ec_ecout_dw(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_dw->at(i);
}
float Branches12::ec_ecout_m2u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m2u->at(i);
}
float Branches12::ec_ecout_m2v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m2v->at(i);
}
float Branches12::ec_ecout_m2w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m2w->at(i);
}
float Branches12::ec_ecout_m3u(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m3u->at(i);
}
float Branches12::ec_ecout_m3v(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m3v->at(i);
}
float Branches12::ec_ecout_m3w(int i) {
  if (i >= _ec_ecin_sec->size())
    return NAN;
  else
    return _ec_ecout_m3w->at(i);
}

float Branches12::cvt_x(int i) { return _cvt_x->at(i); }
float Branches12::cvt_y(int i) { return _cvt_y->at(i); }
float Branches12::cvt_z(int i) { return _cvt_z->at(i); }
float Branches12::fmt_x(int i) { return _fmt_x->at(i); }
float Branches12::fmt_y(int i) { return _fmt_y->at(i); }
float Branches12::fmt_z(int i) { return _fmt_z->at(i); }

float Branches12::cc_nphe_tot(int i) { return _cc_nphe_tot->at(i); }
int Branches12::cc_ltcc_sec(int i) { return _cc_ltcc_sec->at(i); }
float Branches12::cc_ltcc_nphe(int i) { return _cc_ltcc_nphe->at(i); }
float Branches12::cc_ltcc_time(int i) { return _cc_ltcc_time->at(i); }
float Branches12::cc_ltcc_path(int i) { return _cc_ltcc_path->at(i); }
float Branches12::cc_ltcc_theta(int i) { return _cc_ltcc_theta->at(i); }
float Branches12::cc_ltcc_phi(int i) { return _cc_ltcc_phi->at(i); }
float Branches12::cc_ltcc_x(int i) { return _cc_ltcc_x->at(i); }
float Branches12::cc_ltcc_y(int i) { return _cc_ltcc_y->at(i); }
float Branches12::cc_ltcc_z(int i) { return _cc_ltcc_z->at(i); }
int Branches12::cc_htcc_sec(int i) { return _cc_htcc_sec->at(i); }
float Branches12::cc_htcc_nphe(int i) { return _cc_htcc_nphe->at(i); }
float Branches12::cc_htcc_time(int i) { return _cc_htcc_time->at(i); }
float Branches12::cc_htcc_path(int i) { return _cc_htcc_path->at(i); }
float Branches12::cc_htcc_theta(int i) { return _cc_htcc_theta->at(i); }
float Branches12::cc_htcc_phi(int i) { return _cc_htcc_phi->at(i); }
float Branches12::cc_htcc_x(int i) { return _cc_htcc_x->at(i); }
float Branches12::cc_htcc_y(int i) { return _cc_htcc_y->at(i); }
float Branches12::cc_htcc_z(int i) { return _cc_htcc_z->at(i); }
int Branches12::cc_rich_sec(int i) { return _cc_rich_sec->at(i); }
float Branches12::cc_rich_nphe(int i) { return _cc_rich_nphe->at(i); }
float Branches12::cc_rich_time(int i) { return _cc_rich_time->at(i); }
float Branches12::cc_rich_path(int i) { return _cc_rich_path->at(i); }
float Branches12::cc_rich_theta(int i) { return _cc_rich_theta->at(i); }
float Branches12::cc_rich_phi(int i) { return _cc_rich_phi->at(i); }
float Branches12::cc_rich_x(int i) { return _cc_rich_x->at(i); }
float Branches12::cc_rich_y(int i) { return _cc_rich_y->at(i); }
float Branches12::cc_rich_z(int i) { return _cc_rich_z->at(i); }

float Branches12::sc_cnd_time(int i) { return _sc_cnd_time->at(i); }
float Branches12::sc_cnd_path(int i) { return _sc_cnd_path->at(i); }
float Branches12::sc_cnd_energy(int i) { return _sc_cnd_energy->at(i); }
int Branches12::sc_cnd_component(int i) { return _sc_cnd_component->at(i); }
float Branches12::sc_cnd_x(int i) { return _sc_cnd_x->at(i); }
float Branches12::sc_cnd_y(int i) { return _sc_cnd_y->at(i); }
float Branches12::sc_cnd_z(int i) { return _sc_cnd_z->at(i); }
float Branches12::sc_cnd_hx(int i) { return _sc_cnd_hx->at(i); }
float Branches12::sc_cnd_hy(int i) { return _sc_cnd_hy->at(i); }
float Branches12::sc_cnd_hz(int i) { return _sc_cnd_hz->at(i); }
float Branches12::ft_cal_energy(int i) { return _ft_cal_energy->at(i); }
float Branches12::ft_cal_time(int i) { return _ft_cal_time->at(i); }
float Branches12::ft_cal_path(int i) { return _ft_cal_path->at(i); }
float Branches12::ft_cal_x(int i) { return _ft_cal_x->at(i); }
float Branches12::ft_cal_y(int i) { return _ft_cal_y->at(i); }
float Branches12::ft_cal_z(int i) { return _ft_cal_z->at(i); }
float Branches12::ft_cal_dx(int i) { return _ft_cal_dx->at(i); }
float Branches12::ft_cal_dy(int i) { return _ft_cal_dy->at(i); }
float Branches12::ft_cal_radius(int i) { return _ft_cal_radius->at(i); }
float Branches12::ft_hodo_energy(int i) { return _ft_hodo_energy->at(i); }
float Branches12::ft_hodo_time(int i) { return _ft_hodo_time->at(i); }
float Branches12::ft_hodo_path(int i) { return _ft_hodo_path->at(i); }
float Branches12::ft_hodo_x(int i) { return _ft_hodo_x->at(i); }
float Branches12::ft_hodo_y(int i) { return _ft_hodo_y->at(i); }
float Branches12::ft_hodo_z(int i) { return _ft_hodo_z->at(i); }
float Branches12::ft_hodo_dx(int i) { return _ft_hodo_dx->at(i); }
float Branches12::ft_hodo_dy(int i) { return _ft_hodo_dy->at(i); }
float Branches12::ft_hodo_radius(int i) { return _ft_hodo_radius->at(i); }
