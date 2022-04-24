#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "JT21" model

enum state {
    STATE_m,
    STATE_j,
    STATE_mL,
    STATE_hL,
    STATE_Xr1,
    STATE_Xr2,
    STATE_x_Ks,
    STATE_q,
    STATE_r,
    STATE_d,
    STATE_f,
    STATE_f_Ca_B,
    STATE_xf,
    STATE_r_RyR,
    STATE_cn,
    STATE_cc,
    STATE_cd,
    STATE_csl,
    STATE_cs,
    STATE_bc,
    STATE_bd,
    STATE_bs,
    STATE_bsl,
    STATE_V_m,
    STATE_Na_i,
    NUM_STATES,
};

enum parameter {
    PARAM_g_Na,
    PARAM_lambda_Na,
    PARAM_g_NaL,
    PARAM_thL,
    PARAM_KmKo,
    PARAM_KmNaip,
    PARAM_Q10KmNai,
    PARAM_Q10NaK,
    PARAM_g_NaK,
    PARAM_epi,
    PARAM_g_Ks,
    PARAM_pNaK,
    PARAM_g_KATP,
    PARAM_g_Kr,
    PARAM_lambda_K,
    PARAM_g_to,
    PARAM_g_K1,
    PARAM_g_bCl,
    PARAM_Q10CaL,
    PARAM_g_CaL,
    PARAM_Kdact,
    PARAM_KmCai,
    PARAM_KmCao,
    PARAM_KmNai,
    PARAM_KmNao,
    PARAM_Q10NCX,
    PARAM_g_NaCa,
    PARAM_ksat,
    PARAM_nu,
    PARAM_KmPCa,
    PARAM_Q10SLCaP,
    PARAM_g_pCa,
    PARAM_g_bCa,
    PARAM_E_f,
    PARAM_g_f,
    PARAM_Na_sl,
    PARAM_Nao,
    PARAM_K_i,
    PARAM_Ko,
    PARAM_ce,
    PARAM_K_RyR,
    PARAM_alpha_RyR,
    PARAM_beta_RyR,
    PARAM_eta_RyR,
    PARAM_gamma_RyR,
    PARAM_lambda_RyR,
    PARAM_Vc,
    PARAM_Vd,
    PARAM_Vn,
    PARAM_Vs,
    PARAM_Vsl,
    PARAM_J_SERCA_bar,
    PARAM_K_c,
    PARAM_K_n,
    PARAM_B_tot_c,
    PARAM_B_tot_d,
    PARAM_B_tot_s,
    PARAM_B_tot_sl,
    PARAM_k_off_c,
    PARAM_k_off_d,
    PARAM_k_off_s,
    PARAM_k_off_sl,
    PARAM_k_on_c,
    PARAM_k_on_d,
    PARAM_k_on_s,
    PARAM_k_on_sl,
    PARAM_lambda_B,
    PARAM_lambda_B_c,
    PARAM_alpha_d_c,
    PARAM_alpha_n_s,
    PARAM_alpha_sl_c,
    PARAM_lambda_c_d,
    PARAM_lambda_c_i,
    PARAM_lambda_diff,
    PARAM_Cli,
    PARAM_Clo,
    PARAM_Cm,
    PARAM_Frdy,
    PARAM_R,
    PARAM_Temp,
    PARAM_chi,
    PARAM_lambda_c_e,
    PARAM_stim_amplitude,
    PARAM_stim_duration,
    PARAM_stim_period,
    PARAM_stim_start,
    NUM_PARAMS,
};

enum monitored {
    NUM_MONITORED,
};

// Init state values
void init_state_values(double *states)
{
    states[STATE_m] = 0.0047593841564;
    states[STATE_j] = 0.59186803670583;
    states[STATE_mL] = 0.00056156997976;
    states[STATE_hL] = 0.19825212801003;
    states[STATE_Xr1] = 0.02598687936732;
    states[STATE_Xr2] = 0.46218104728376;
    states[STATE_x_Ks] = 0.00418046147173;
    states[STATE_q] = 0.89186079496414;
    states[STATE_r] = 0.00415791788969;
    states[STATE_d] = 3.47393157e-06;
    states[STATE_f] = 0.99328007454758;
    states[STATE_f_Ca_B] = 0.05601726193501;
    states[STATE_xf] = 0.08014973245989;
    states[STATE_r_RyR] = 0.99999973974165;
    states[STATE_cn] = 0.69858720814252;
    states[STATE_cc] = 0.00011382799663;
    states[STATE_cd] = 0.00020977171788;
    states[STATE_csl] = 0.00011798425266;
    states[STATE_cs] = 0.69792301896013;
    states[STATE_bc] = 0.00839395214247;
    states[STATE_bd] = 0.05570692056887;
    states[STATE_bs] = 31.06660220440676;
    states[STATE_bsl] = 0.10584192890313;
    states[STATE_V_m] = -80.42101165870085;
    states[STATE_Na_i] = 8.07942377455463;
}

// Default parameter values
void init_parameters_values(double *parameters)
{
    parameters[PARAM_g_Na] = 0.36;
    parameters[PARAM_lambda_Na] = 1.0;
    parameters[PARAM_g_NaL] = 0.03;
    parameters[PARAM_thL] = 200.0;
    parameters[PARAM_KmKo] = 1.5;
    parameters[PARAM_KmNaip] = 11.0;
    parameters[PARAM_Q10KmNai] = 1.4;
    parameters[PARAM_Q10NaK] = 1.6;
    parameters[PARAM_g_NaK] = 2.3939999999999997;
    parameters[PARAM_epi] = 0.0;
    parameters[PARAM_g_Ks] = 0.0127;
    parameters[PARAM_pNaK] = 0.018;
    parameters[PARAM_g_KATP] = 0.01;
    parameters[PARAM_g_Kr] = 0.029399999999999996;
    parameters[PARAM_lambda_K] = 1.0;
    parameters[PARAM_g_to] = 0.105;
    parameters[PARAM_g_K1] = 0.27;
    parameters[PARAM_g_bCl] = 0.003;
    parameters[PARAM_Q10CaL] = 1.8;
    parameters[PARAM_g_CaL] = 0.8316000000000001;
    parameters[PARAM_Kdact] = 0.00015;
    parameters[PARAM_KmCai] = 0.0036;
    parameters[PARAM_KmCao] = 1.3;
    parameters[PARAM_KmNai] = 12.3;
    parameters[PARAM_KmNao] = 87.5;
    parameters[PARAM_Q10NCX] = 1.6;
    parameters[PARAM_g_NaCa] = 14.1;
    parameters[PARAM_ksat] = 0.3;
    parameters[PARAM_nu] = 0.3;
    parameters[PARAM_KmPCa] = 0.0005;
    parameters[PARAM_Q10SLCaP] = 2.35;
    parameters[PARAM_g_pCa] = 0.12;
    parameters[PARAM_g_bCa] = 0.0021;
    parameters[PARAM_E_f] = -17.0;
    parameters[PARAM_g_f] = 0.01;
    parameters[PARAM_Na_sl] = 8.0;
    parameters[PARAM_Nao] = 140.0;
    parameters[PARAM_K_i] = 120.0;
    parameters[PARAM_Ko] = 5.0;
    parameters[PARAM_ce] = 0.42;
    parameters[PARAM_K_RyR] = 0.015;
    parameters[PARAM_alpha_RyR] = 0.02075;
    parameters[PARAM_beta_RyR] = 0.042;
    parameters[PARAM_eta_RyR] = 1e-05;
    parameters[PARAM_gamma_RyR] = 0.001;
    parameters[PARAM_lambda_RyR] = 0.63;
    parameters[PARAM_Vc] = 0.917;
    parameters[PARAM_Vd] = 0.001;
    parameters[PARAM_Vn] = 0.05;
    parameters[PARAM_Vs] = 0.004;
    parameters[PARAM_Vsl] = 0.028;
    parameters[PARAM_J_SERCA_bar] = 0.00016;
    parameters[PARAM_K_c] = 0.00025;
    parameters[PARAM_K_n] = 1.7;
    parameters[PARAM_B_tot_c] = 0.063;
    parameters[PARAM_B_tot_d] = 2.7;
    parameters[PARAM_B_tot_s] = 60.0;
    parameters[PARAM_B_tot_sl] = 1.45;
    parameters[PARAM_k_off_c] = 0.03;
    parameters[PARAM_k_off_d] = 1.0;
    parameters[PARAM_k_off_s] = 65.0;
    parameters[PARAM_k_off_sl] = 0.15;
    parameters[PARAM_k_on_c] = 40.0;
    parameters[PARAM_k_on_d] = 100.0;
    parameters[PARAM_k_on_s] = 100.0;
    parameters[PARAM_k_on_sl] = 100.0;
    parameters[PARAM_lambda_B] = 1.0;
    parameters[PARAM_lambda_B_c] = 1.0;
    parameters[PARAM_alpha_d_c] = 0.0027;
    parameters[PARAM_alpha_n_s] = 0.0093;
    parameters[PARAM_alpha_sl_c] = 0.3;
    parameters[PARAM_lambda_c_d] = 1.0;
    parameters[PARAM_lambda_c_i] = 1.0;
    parameters[PARAM_lambda_diff] = 1.0;
    parameters[PARAM_Cli] = 15.0;
    parameters[PARAM_Clo] = 150.0;
    parameters[PARAM_Cm] = 0.01;
    parameters[PARAM_Frdy] = 96.485;
    parameters[PARAM_R] = 8.314;
    parameters[PARAM_Temp] = 310.0;
    parameters[PARAM_chi] = 0.9;
    parameters[PARAM_lambda_c_e] = 1.0;
    parameters[PARAM_stim_amplitude] = 5.0;
    parameters[PARAM_stim_duration] = 20.0;
    parameters[PARAM_stim_period] = 10000.0;
    parameters[PARAM_stim_start] = 0.0;
}

// State index
int state_index(const char name[])
{
    if (strcmp(name, "m") == 0) {
        return STATE_m;
    } else if (strcmp(name, "j") == 0) {
        return STATE_j;
    } else if (strcmp(name, "mL") == 0) {
        return STATE_mL;
    } else if (strcmp(name, "hL") == 0) {
        return STATE_hL;
    } else if (strcmp(name, "Xr1") == 0) {
        return STATE_Xr1;
    } else if (strcmp(name, "Xr2") == 0) {
        return STATE_Xr2;
    } else if (strcmp(name, "x_Ks") == 0) {
        return STATE_x_Ks;
    } else if (strcmp(name, "q") == 0) {
        return STATE_q;
    } else if (strcmp(name, "r") == 0) {
        return STATE_r;
    } else if (strcmp(name, "d") == 0) {
        return STATE_d;
    } else if (strcmp(name, "f") == 0) {
        return STATE_f;
    } else if (strcmp(name, "f_Ca_B") == 0) {
        return STATE_f_Ca_B;
    } else if (strcmp(name, "xf") == 0) {
        return STATE_xf;
    } else if (strcmp(name, "r_RyR") == 0) {
        return STATE_r_RyR;
    } else if (strcmp(name, "cn") == 0) {
        return STATE_cn;
    } else if (strcmp(name, "cc") == 0) {
        return STATE_cc;
    } else if (strcmp(name, "cd") == 0) {
        return STATE_cd;
    } else if (strcmp(name, "csl") == 0) {
        return STATE_csl;
    } else if (strcmp(name, "cs") == 0) {
        return STATE_cs;
    } else if (strcmp(name, "bc") == 0) {
        return STATE_bc;
    } else if (strcmp(name, "bd") == 0) {
        return STATE_bd;
    } else if (strcmp(name, "bs") == 0) {
        return STATE_bs;
    } else if (strcmp(name, "bsl") == 0) {
        return STATE_bsl;
    } else if (strcmp(name, "V_m") == 0) {
        return STATE_V_m;
    } else if (strcmp(name, "Na_i") == 0) {
        return STATE_Na_i;
    }
    return -1;
}

// Parameter index
int parameter_index(const char name[])
{
    if (strcmp(name, "g_Na") == 0) {
        return PARAM_g_Na;
    } else if (strcmp(name, "lambda_Na") == 0) {
        return PARAM_lambda_Na;
    } else if (strcmp(name, "g_NaL") == 0) {
        return PARAM_g_NaL;
    } else if (strcmp(name, "thL") == 0) {
        return PARAM_thL;
    } else if (strcmp(name, "KmKo") == 0) {
        return PARAM_KmKo;
    } else if (strcmp(name, "KmNaip") == 0) {
        return PARAM_KmNaip;
    } else if (strcmp(name, "Q10KmNai") == 0) {
        return PARAM_Q10KmNai;
    } else if (strcmp(name, "Q10NaK") == 0) {
        return PARAM_Q10NaK;
    } else if (strcmp(name, "g_NaK") == 0) {
        return PARAM_g_NaK;
    } else if (strcmp(name, "epi") == 0) {
        return PARAM_epi;
    } else if (strcmp(name, "g_Ks") == 0) {
        return PARAM_g_Ks;
    } else if (strcmp(name, "pNaK") == 0) {
        return PARAM_pNaK;
    } else if (strcmp(name, "g_KATP") == 0) {
        return PARAM_g_KATP;
    } else if (strcmp(name, "g_Kr") == 0) {
        return PARAM_g_Kr;
    } else if (strcmp(name, "lambda_K") == 0) {
        return PARAM_lambda_K;
    } else if (strcmp(name, "g_to") == 0) {
        return PARAM_g_to;
    } else if (strcmp(name, "g_K1") == 0) {
        return PARAM_g_K1;
    } else if (strcmp(name, "g_bCl") == 0) {
        return PARAM_g_bCl;
    } else if (strcmp(name, "Q10CaL") == 0) {
        return PARAM_Q10CaL;
    } else if (strcmp(name, "g_CaL") == 0) {
        return PARAM_g_CaL;
    } else if (strcmp(name, "Kdact") == 0) {
        return PARAM_Kdact;
    } else if (strcmp(name, "KmCai") == 0) {
        return PARAM_KmCai;
    } else if (strcmp(name, "KmCao") == 0) {
        return PARAM_KmCao;
    } else if (strcmp(name, "KmNai") == 0) {
        return PARAM_KmNai;
    } else if (strcmp(name, "KmNao") == 0) {
        return PARAM_KmNao;
    } else if (strcmp(name, "Q10NCX") == 0) {
        return PARAM_Q10NCX;
    } else if (strcmp(name, "g_NaCa") == 0) {
        return PARAM_g_NaCa;
    } else if (strcmp(name, "ksat") == 0) {
        return PARAM_ksat;
    } else if (strcmp(name, "nu") == 0) {
        return PARAM_nu;
    } else if (strcmp(name, "KmPCa") == 0) {
        return PARAM_KmPCa;
    } else if (strcmp(name, "Q10SLCaP") == 0) {
        return PARAM_Q10SLCaP;
    } else if (strcmp(name, "g_pCa") == 0) {
        return PARAM_g_pCa;
    } else if (strcmp(name, "g_bCa") == 0) {
        return PARAM_g_bCa;
    } else if (strcmp(name, "E_f") == 0) {
        return PARAM_E_f;
    } else if (strcmp(name, "g_f") == 0) {
        return PARAM_g_f;
    } else if (strcmp(name, "Na_sl") == 0) {
        return PARAM_Na_sl;
    } else if (strcmp(name, "Nao") == 0) {
        return PARAM_Nao;
    } else if (strcmp(name, "K_i") == 0) {
        return PARAM_K_i;
    } else if (strcmp(name, "Ko") == 0) {
        return PARAM_Ko;
    } else if (strcmp(name, "ce") == 0) {
        return PARAM_ce;
    } else if (strcmp(name, "K_RyR") == 0) {
        return PARAM_K_RyR;
    } else if (strcmp(name, "alpha_RyR") == 0) {
        return PARAM_alpha_RyR;
    } else if (strcmp(name, "beta_RyR") == 0) {
        return PARAM_beta_RyR;
    } else if (strcmp(name, "eta_RyR") == 0) {
        return PARAM_eta_RyR;
    } else if (strcmp(name, "gamma_RyR") == 0) {
        return PARAM_gamma_RyR;
    } else if (strcmp(name, "lambda_RyR") == 0) {
        return PARAM_lambda_RyR;
    } else if (strcmp(name, "Vc") == 0) {
        return PARAM_Vc;
    } else if (strcmp(name, "Vd") == 0) {
        return PARAM_Vd;
    } else if (strcmp(name, "Vn") == 0) {
        return PARAM_Vn;
    } else if (strcmp(name, "Vs") == 0) {
        return PARAM_Vs;
    } else if (strcmp(name, "Vsl") == 0) {
        return PARAM_Vsl;
    } else if (strcmp(name, "J_SERCA_bar") == 0) {
        return PARAM_J_SERCA_bar;
    } else if (strcmp(name, "K_c") == 0) {
        return PARAM_K_c;
    } else if (strcmp(name, "K_n") == 0) {
        return PARAM_K_n;
    } else if (strcmp(name, "B_tot_c") == 0) {
        return PARAM_B_tot_c;
    } else if (strcmp(name, "B_tot_d") == 0) {
        return PARAM_B_tot_d;
    } else if (strcmp(name, "B_tot_s") == 0) {
        return PARAM_B_tot_s;
    } else if (strcmp(name, "B_tot_sl") == 0) {
        return PARAM_B_tot_sl;
    } else if (strcmp(name, "k_off_c") == 0) {
        return PARAM_k_off_c;
    } else if (strcmp(name, "k_off_d") == 0) {
        return PARAM_k_off_d;
    } else if (strcmp(name, "k_off_s") == 0) {
        return PARAM_k_off_s;
    } else if (strcmp(name, "k_off_sl") == 0) {
        return PARAM_k_off_sl;
    } else if (strcmp(name, "k_on_c") == 0) {
        return PARAM_k_on_c;
    } else if (strcmp(name, "k_on_d") == 0) {
        return PARAM_k_on_d;
    } else if (strcmp(name, "k_on_s") == 0) {
        return PARAM_k_on_s;
    } else if (strcmp(name, "k_on_sl") == 0) {
        return PARAM_k_on_sl;
    } else if (strcmp(name, "lambda_B") == 0) {
        return PARAM_lambda_B;
    } else if (strcmp(name, "lambda_B_c") == 0) {
        return PARAM_lambda_B_c;
    } else if (strcmp(name, "alpha_d_c") == 0) {
        return PARAM_alpha_d_c;
    } else if (strcmp(name, "alpha_n_s") == 0) {
        return PARAM_alpha_n_s;
    } else if (strcmp(name, "alpha_sl_c") == 0) {
        return PARAM_alpha_sl_c;
    } else if (strcmp(name, "lambda_c_d") == 0) {
        return PARAM_lambda_c_d;
    } else if (strcmp(name, "lambda_c_i") == 0) {
        return PARAM_lambda_c_i;
    } else if (strcmp(name, "lambda_diff") == 0) {
        return PARAM_lambda_diff;
    } else if (strcmp(name, "Cli") == 0) {
        return PARAM_Cli;
    } else if (strcmp(name, "Clo") == 0) {
        return PARAM_Clo;
    } else if (strcmp(name, "Cm") == 0) {
        return PARAM_Cm;
    } else if (strcmp(name, "Frdy") == 0) {
        return PARAM_Frdy;
    } else if (strcmp(name, "R") == 0) {
        return PARAM_R;
    } else if (strcmp(name, "Temp") == 0) {
        return PARAM_Temp;
    } else if (strcmp(name, "chi") == 0) {
        return PARAM_chi;
    } else if (strcmp(name, "lambda_c_e") == 0) {
        return PARAM_lambda_c_e;
    } else if (strcmp(name, "stim_amplitude") == 0) {
        return PARAM_stim_amplitude;
    } else if (strcmp(name, "stim_duration") == 0) {
        return PARAM_stim_duration;
    } else if (strcmp(name, "stim_period") == 0) {
        return PARAM_stim_period;
    } else if (strcmp(name, "stim_start") == 0) {
        return PARAM_stim_start;
    }
    return -1;
}

// Compute a forward step using the explicit Euler algorithm to the
// JT21 ODE
void step_FE_singlecell(double *__restrict states, const double t, const double dt,
                        const double *__restrict parameters)
{

    // Assign states
    const double m = states[STATE_m];
    const double j = states[STATE_j];
    const double mL = states[STATE_mL];
    const double hL = states[STATE_hL];
    const double Xr1 = states[STATE_Xr1];
    const double Xr2 = states[STATE_Xr2];
    const double x_Ks = states[STATE_x_Ks];
    const double q = states[STATE_q];
    const double r = states[STATE_r];
    const double d = states[STATE_d];
    const double f = states[STATE_f];
    const double f_Ca_B = states[STATE_f_Ca_B];
    const double xf = states[STATE_xf];
    const double r_RyR = states[STATE_r_RyR];
    const double cn = states[STATE_cn];
    const double cc = states[STATE_cc];
    const double cd = states[STATE_cd];
    const double csl = states[STATE_csl];
    const double cs = states[STATE_cs];
    const double bc = states[STATE_bc];
    const double bd = states[STATE_bd];
    const double bs = states[STATE_bs];
    const double bsl = states[STATE_bsl];
    const double V_m = states[STATE_V_m];
    const double Na_i = states[STATE_Na_i];

    // Assign parameters
    const double g_Na = parameters[PARAM_g_Na];
    const double lambda_Na = parameters[PARAM_lambda_Na];
    const double g_NaL = parameters[PARAM_g_NaL];
    const double thL = parameters[PARAM_thL];
    const double KmKo = parameters[PARAM_KmKo];
    const double KmNaip = parameters[PARAM_KmNaip];
    const double g_NaK = parameters[PARAM_g_NaK];
    const double epi = parameters[PARAM_epi];
    const double g_Ks = parameters[PARAM_g_Ks];
    const double pNaK = parameters[PARAM_pNaK];
    const double g_KATP = parameters[PARAM_g_KATP];
    const double g_Kr = parameters[PARAM_g_Kr];
    const double lambda_K = parameters[PARAM_lambda_K];
    const double g_to = parameters[PARAM_g_to];
    const double g_K1 = parameters[PARAM_g_K1];
    const double g_bCl = parameters[PARAM_g_bCl];
    const double Q10CaL = parameters[PARAM_Q10CaL];
    const double g_CaL = parameters[PARAM_g_CaL];
    const double Kdact = parameters[PARAM_Kdact];
    const double KmCai = parameters[PARAM_KmCai];
    const double KmCao = parameters[PARAM_KmCao];
    const double KmNai = parameters[PARAM_KmNai];
    const double KmNao = parameters[PARAM_KmNao];
    const double Q10NCX = parameters[PARAM_Q10NCX];
    const double g_NaCa = parameters[PARAM_g_NaCa];
    const double ksat = parameters[PARAM_ksat];
    const double nu = parameters[PARAM_nu];
    const double KmPCa = parameters[PARAM_KmPCa];
    const double Q10SLCaP = parameters[PARAM_Q10SLCaP];
    const double g_pCa = parameters[PARAM_g_pCa];
    const double g_bCa = parameters[PARAM_g_bCa];
    const double E_f = parameters[PARAM_E_f];
    const double g_f = parameters[PARAM_g_f];
    const double Nao = parameters[PARAM_Nao];
    const double K_i = parameters[PARAM_K_i];
    const double Ko = parameters[PARAM_Ko];
    const double ce = parameters[PARAM_ce];
    const double K_RyR = parameters[PARAM_K_RyR];
    const double alpha_RyR = parameters[PARAM_alpha_RyR];
    const double beta_RyR = parameters[PARAM_beta_RyR];
    const double eta_RyR = parameters[PARAM_eta_RyR];
    const double gamma_RyR = parameters[PARAM_gamma_RyR];
    const double lambda_RyR = parameters[PARAM_lambda_RyR];
    const double Vc = parameters[PARAM_Vc];
    const double Vd = parameters[PARAM_Vd];
    const double Vn = parameters[PARAM_Vn];
    const double Vs = parameters[PARAM_Vs];
    const double Vsl = parameters[PARAM_Vsl];
    const double J_SERCA_bar = parameters[PARAM_J_SERCA_bar];
    const double K_c = parameters[PARAM_K_c];
    const double K_n = parameters[PARAM_K_n];
    const double B_tot_c = parameters[PARAM_B_tot_c];
    const double B_tot_d = parameters[PARAM_B_tot_d];
    const double B_tot_s = parameters[PARAM_B_tot_s];
    const double B_tot_sl = parameters[PARAM_B_tot_sl];
    const double k_off_c = parameters[PARAM_k_off_c];
    const double k_off_d = parameters[PARAM_k_off_d];
    const double k_off_s = parameters[PARAM_k_off_s];
    const double k_off_sl = parameters[PARAM_k_off_sl];
    const double k_on_c = parameters[PARAM_k_on_c];
    const double k_on_d = parameters[PARAM_k_on_d];
    const double k_on_s = parameters[PARAM_k_on_s];
    const double k_on_sl = parameters[PARAM_k_on_sl];
    const double lambda_B = parameters[PARAM_lambda_B];
    const double lambda_B_c = parameters[PARAM_lambda_B_c];
    const double alpha_d_c = parameters[PARAM_alpha_d_c];
    const double alpha_n_s = parameters[PARAM_alpha_n_s];
    const double alpha_sl_c = parameters[PARAM_alpha_sl_c];
    const double lambda_c_d = parameters[PARAM_lambda_c_d];
    const double lambda_c_i = parameters[PARAM_lambda_c_i];
    const double lambda_diff = parameters[PARAM_lambda_diff];
    const double Cli = parameters[PARAM_Cli];
    const double Clo = parameters[PARAM_Clo];
    const double Cm = parameters[PARAM_Cm];
    const double Frdy = parameters[PARAM_Frdy];
    const double R = parameters[PARAM_R];
    const double Temp = parameters[PARAM_Temp];
    const double chi = parameters[PARAM_chi];
    const double lambda_c_e = parameters[PARAM_lambda_c_e];
    const double stim_amplitude = parameters[PARAM_stim_amplitude];
    const double stim_duration = parameters[PARAM_stim_duration];
    const double stim_period = parameters[PARAM_stim_period];
    const double stim_start = parameters[PARAM_stim_start];

    // Expressions for the Reversal potentials component
    const double FoRT = Frdy / (R * Temp);
    const double ena = log(Nao / Na_i) / FoRT;
    const double ek = log(Ko / K_i) / FoRT;
    const double eca_sl = log(ce / csl) / (2. * FoRT);
    const double ecl = log(Cli / Clo) / FoRT;
    const double Qpow = -31. + Temp / 10.;

    // Expressions for the I_Na component
    const double mss = 1.0
                       / ((1. + 0.00177610354573438 * exp(-V_m / 9.))
                          * (1. + 0.00177610354573438 * exp(-V_m / 9.)));
    const double taum =
            0.06 * exp(-((-0.0980392156862745 + V_m / 51.) * (-0.0980392156862745 + V_m / 51.)))
            + 0.13 * exp(-((2.875 + V_m / 16.) * (2.875 + V_m / 16.)));
    const double aj =
            (V_m >= -38. ? 0.
                         : (38. + V_m) * (-25000.0 * exp(0.2 * V_m) - 7.0e-6 * exp(-0.04 * V_m))
                                   / (1. + 19623624323.6513 * exp(0.3 * V_m)));
    const double bj = (V_m >= -38.218 ? 0.6 * exp(0.09 * V_m) / (1. + exp(-40. - V_m))
                                      : 0.02 * exp(-0.01 * V_m)
                                                / (1. + 0.00369786371648293 * exp(-0.14 * V_m)));
    const double tauj = 1.0 / (aj + bj);
    const double jss =
            1.0
            / ((1. + 29310.8866998062 * exp(V_m / 7.)) * (1. + 29310.8866998062 * exp(V_m / 7.)));
    const double I_Na = g_Na * lambda_Na * (m * m * m) * (-ena + V_m) * j;
    const double dm_dt = (-m + mss) / taum;
    states[STATE_m] = dt * dm_dt + m;
    const double dj_dt = (-j + jss) / tauj;
    states[STATE_j] = dt * dj_dt + j;

    // Expressions for the I_NaL component
    const double mLss = 1.0 / (1. + 0.000184105793667579 * exp(-V_m / 5.));
    const double tm =
            1.0 / (9.58097877207697 * exp(V_m / 35.) + 2.29642676604141e-5 * exp(-V_m / 6.));
    const double tmL = tm;
    const double hLss = 1.0 / (1. + 124658.506952 * exp(2. * V_m / 15.));
    const double GNaL = (epi == 1. ? g_NaL : 0.6 * g_NaL);
    const double I_NaL = lambda_Na * (-ena + V_m) * GNaL * hL * mL;
    const double dmL_dt = (-mL + mLss) / tmL;
    states[STATE_mL] = dt * dmL_dt + mL;
    const double dhL_dt = (-hL + hLss) / thL;
    states[STATE_hL] = dt * dhL_dt + hL;

    // Expressions for the I_NaK component
    const double sigma = -0.142857142857143 + exp(Nao / 67.) / 7.;
    const double fNaK =
            1.0 / (1. + 0.12 * exp(-0.1 * FoRT * V_m) + 0.037 * exp(-FoRT * V_m) * sigma);
    const double I_NaK = Ko * g_NaK * fNaK / ((1. + pow(KmNaip, 4.) / pow(Na_i, 4.)) * (KmKo + Ko));

    // Expressions for the I_Kr component
    const double Xr1_inf = 1.0 / (1.0 + 0.0146327985189294 * exp(-0.204081632653061 * V_m));
    const double alpha_Xr1 = 450.0 / (1.0 + 0.0111089965382 * exp(-0.1 * V_m));
    const double beta_Xr1 = 6.0 / (1.0 + 13.5813245226 * exp(0.0869565217391 * V_m));
    const double tau_Xr1 = 1.0 * alpha_Xr1 * beta_Xr1;
    const double Xr2_infinity = 1.0 / (1.0 + exp(44. / 25. + V_m / 50.));
    const double alpha_Xr2 = 3.0 / (1.0 + exp(-3. - V_m / 20.));
    const double beta_Xr2 = 1.12 / (1.0 + exp(-3. + V_m / 20.));
    const double tau_Xr2 = 1.0 * alpha_Xr2 * beta_Xr2;
    const double I_Kr = 0.430331482912 * g_Kr * sqrt(Ko) * (-ek + V_m) * Xr1 * Xr2;
    const double dXr1_dt = (-Xr1 + Xr1_inf) / tau_Xr1;
    states[STATE_Xr1] = dt * dXr1_dt + Xr1;
    const double dXr2_dt = (-Xr2 + Xr2_infinity) / tau_Xr2;
    states[STATE_Xr2] = dt * dXr2_dt + Xr2;

    // Expressions for the I_Ks component
    const double eks = log((Ko + Nao * pNaK) / (K_i + pNaK * Na_i)) / FoRT;
    const double xsss = 1.0 / (1. + 0.76228973079 * exp(-V_m / 14.));
    const double tauxs = 990. / (1. + 0.842460441617 * exp(-V_m / 14.));
    const double I_Ks = g_Ks * (x_Ks * x_Ks) * (-eks + V_m);
    const double dx_Ks_dt = (-x_Ks + xsss) / tauxs;
    states[STATE_x_Ks] = dt * dx_Ks_dt + x_Ks;

    // Expressions for the i_to component
    const double q_inf = 1.0 / (1.0 + 58.9637634804 * exp(0.0769230769231 * V_m));
    const double tau_q =
            6. + 39. / (0.0168716780457 * exp(-0.08 * V_m) + 6.46648051673 * exp(0.1 * V_m));
    const double r_inf = 1.0 / (1.0 + 3.28489055021 * exp(-0.0533333333333 * V_m));
    const double tau_r =
            2.75 + 14.4 / (0.0207698622486 * exp(-0.12 * V_m) + 15.7194688773 * exp(0.09 * V_m));
    const double I_to = g_to * (-ek + V_m) * q * r;
    const double dq_dt = (-q + q_inf) / tau_q;
    states[STATE_q] = dt * dq_dt + q;
    const double dr_dt = (-r + r_inf) / tau_r;
    states[STATE_r] = dt * dr_dt + r;

    // Expressions for the I_K1 component
    const double aK1 = 1.0 / (1. + 7.50455791508e-6 * exp(0.2 * V_m - 0.2 * ek));
    const double bK1 = (0.745912348821 * exp(0.08 * V_m - 0.08 * ek)
                        + 3.32464030033e-16 * exp(0.06 * V_m - 0.06 * ek))
                       / (1. + 0.0820849986239 * exp(0.5 * ek - 0.5 * V_m));
    const double K1ss = aK1 / (aK1 + bK1);
    const double I_K1 = 0.430331482912 * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m) * K1ss;

    // Expressions for the I_bCl component
    const double I_bCl = g_bCl * (-ecl + V_m);

    // Expressions for the I_Ca component
    const double fss = 1.0 / (1. + 48.8565712749872 * exp(V_m / 9.))
                       + 0.6 / (1. + 12.1824939607035 * exp(-V_m / 20.));
    const double dss = 1.0 / (1. + exp(-5. / 6. - V_m / 6.));
    const double taud = (fabs(5. + V_m) < 0.02 ? 2.38095238095238
                                               : (1. - 0.434598208507078 * exp(-V_m / 6.)) * dss
                                                         / (0.175 + 0.035 * V_m));
    const double tauf = 1.0 / (0.02 + 0.02 * exp(-((0.493 + 0.034 * V_m) * (0.493 + 0.034 * V_m))));
    const double ibarca_j = 4. * Frdy * g_CaL * (-0.34 * ce + 0.34 * cd * exp(2. * FoRT * V_m))
                            * FoRT * V_m / (-1. + exp(2. * FoRT * V_m));
    const double I_CaL = lambda_c_d * pow(Q10CaL, Qpow) * (1. - f_Ca_B) * d * f * ibarca_j;
    const double dd_dt = (-d + dss) / taud;
    states[STATE_d] = dt * dd_dt + d;
    const double df_dt = (-f + fss) / tauf;
    states[STATE_f] = dt * df_dt + f;
    const double df_Ca_B_dt = -0.012 * f_Ca_B + (1.7 - 1.7 * f_Ca_B) * cd;
    states[STATE_f_Ca_B] = dt * df_Ca_B_dt + f_Ca_B;

    // Expressions for the I_NCX component
    const double Ka_sl = 1.0 / (1. + (Kdact * Kdact) / (csl * csl));
    const double s1_sl = ce * (Na_i * Na_i * Na_i) * exp(nu * FoRT * V_m);
    const double s2_sl = (Nao * Nao * Nao) * csl * exp((-1. + nu) * FoRT * V_m);
    const double s3_sl =
            KmCao * (Na_i * Na_i * Na_i) + ce * (Na_i * Na_i * Na_i) + (Nao * Nao * Nao) * csl
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_i * Na_i * Na_i) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + csl / KmCai) * csl;
    const double I_NaCa = g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl
                          / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);

    // Expressions for the I_PCa component
    const double I_pCa = g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl)
                         / ((KmPCa * KmPCa) + (csl * csl));

    // Expressions for the I_CaBK component
    const double I_bCa = g_bCa * lambda_c_e * (-eca_sl + V_m);

    // Expressions for the I_f component
    const double xf_inf = 1.0 / (1. + 5956538.01318461 * exp(V_m / 5.));
    const double tau_xf = 1900. / (1. + 4.48168907033806 * exp(V_m / 10.));
    const double I_f = g_f * (-E_f + V_m) * xf;
    const double dxf_dt = (-xf + xf_inf) / tau_xf;
    states[STATE_xf] = dt * dxf_dt + xf;

    // Expressions for the I_KATP component
    const double I_KATP =
            0.60295079490657 * g_KATP * pow(Ko, 0.3) * (-ek + V_m) / (40. + 0.0875 * V_m);

    // Expressions for the Ca Fluxes component
    const double J_CaL = -Cm * chi * I_CaL / (2. * Frdy);
    const double J_pCa = -Cm * chi * I_pCa / (2. * Frdy);
    const double J_bCa = -Cm * chi * I_bCa / (2. * Frdy);
    const double J_NaCa = Cm * chi * I_NaCa / Frdy;
    const double J_e_sl = J_NaCa + J_bCa + J_pCa;
    const double Q10SERCA = 2.60000000000000;
    const double J_SERCA = J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                           * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n))
                           / (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n));
    const double J_n_s = alpha_n_s * lambda_c_i * lambda_diff * (-cs + cn);
    const double J_sl_c = alpha_sl_c * lambda_c_i * lambda_diff * (-cc + csl);
    const double J_d_c = alpha_d_c * lambda_c_d * lambda_diff * (-cc + cd);

    // Expressions for the RyRs component
    const double p = 1.0 / (1. + (K_RyR * K_RyR * K_RyR) / (cd * cd * cd));
    const double J_RyR_active = alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p * r_RyR;
    const double J_leak = alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i * (-csl + cs);
    const double J_RyR = J_RyR_active + J_leak;
    const double dr_RyR_dt =
            eta_RyR * (1. - r_RyR) / p - J_RyR_active / (beta_RyR * lambda_RyR * lambda_c_i);
    states[STATE_r_RyR] = dt * dr_RyR_dt + r_RyR;

    // Expressions for the Ca Buffers component
    const double J_c_b =
            Vc * (-k_off_c * bc + k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c) * cc);
    const double J_d_b =
            Vd * (-k_off_d * bd + k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c) * cd);
    const double J_s_b = Vs * (-k_off_s * bs + k_on_s * (-bs + B_tot_s * lambda_B) * cs);
    const double J_sl_b =
            Vsl * (-k_off_sl * bsl + k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c) * csl);

    // Expressions for the Ca Concentrations component
    const double dcn_dt = (1.0 * J_SERCA - 1.0 * J_n_s) / Vn;
    states[STATE_cn] = dt * dcn_dt + cn;
    const double dcc_dt = (1.0 * J_d_c + 1.0 * J_sl_c - 1.0 * J_SERCA - 1.0 * J_c_b) / Vc;
    states[STATE_cc] = dt * dcc_dt + cc;
    const double dcd_dt = (1.0 * J_CaL - 1.0 * J_d_b - 1.0 * J_d_c) / Vd;
    states[STATE_cd] = dt * dcd_dt + cd;
    const double dcsl_dt = (1.0 * J_RyR + 1.0 * J_e_sl - 1.0 * J_sl_b - 1.0 * J_sl_c) / Vsl;
    states[STATE_csl] = dt * dcsl_dt + csl;
    const double dcs_dt = (1.0 * J_n_s - 1.0 * J_RyR - 1.0 * J_s_b) / Vs;
    states[STATE_cs] = dt * dcs_dt + cs;

    // Expressions for the Ca Buffer Concentrations component
    const double dbc_dt = 1.0 * J_c_b / Vc;
    states[STATE_bc] = dt * dbc_dt + bc;
    const double dbd_dt = 1.0 * J_d_b / Vd;
    states[STATE_bd] = dt * dbd_dt + bd;
    const double dbs_dt = 1.0 * J_s_b / Vs;
    states[STATE_bs] = dt * dbs_dt + bs;
    const double dbsl_dt = 1.0 * J_sl_b / Vsl;
    states[STATE_bsl] = dt * dbsl_dt + bsl;

    // Expressions for the Membrane potential component
    const double i_Stim =
            (V_m < -40. ? 1. : 0.)
            * (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                               && t - stim_period * floor(t / stim_period) >= stim_start
                       ? -stim_amplitude
                       : 0.);
    const double I_tot = I_CaL + I_K1 + I_KATP + I_Kr + I_Ks + I_Na + I_NaCa + I_NaK + I_NaL + I_bCa
                         + I_bCl + I_f + I_pCa + I_to;
    const double dV_m_dt = -I_tot - i_Stim;
    states[STATE_V_m] = dt * dV_m_dt + V_m;

    // Expressions for the Sodium concentration component
    const double I_Na_tot = 3. * I_NaCa + 3. * I_NaK + 0.3293 * I_f + I_Na + I_NaL;
    const double J_Na = -Cm * chi * I_Na_tot / Frdy;
    const double dNa_i_dt = J_Na;
    states[STATE_Na_i] = dt * dNa_i_dt + Na_i;
}

// Compute a forward step using the rush larsen algorithm to the JT21 ODE
void step_GRL1_singlecell(double *__restrict states, const double t, const double dt,
                          const double *__restrict parameters)
{

    // Assign states
    const double m = states[STATE_m];
    const double j = states[STATE_j];
    const double mL = states[STATE_mL];
    const double hL = states[STATE_hL];
    const double Xr1 = states[STATE_Xr1];
    const double Xr2 = states[STATE_Xr2];
    const double x_Ks = states[STATE_x_Ks];
    const double q = states[STATE_q];
    const double r = states[STATE_r];
    const double d = states[STATE_d];
    const double f = states[STATE_f];
    const double f_Ca_B = states[STATE_f_Ca_B];
    const double xf = states[STATE_xf];
    const double r_RyR = states[STATE_r_RyR];
    const double cn = states[STATE_cn];
    const double cc = states[STATE_cc];
    const double cd = states[STATE_cd];
    const double csl = states[STATE_csl];
    const double cs = states[STATE_cs];
    const double bc = states[STATE_bc];
    const double bd = states[STATE_bd];
    const double bs = states[STATE_bs];
    const double bsl = states[STATE_bsl];
    const double V_m = states[STATE_V_m];
    const double Na_i = states[STATE_Na_i];

    // Assign parameters
    const double g_Na = parameters[PARAM_g_Na];
    const double lambda_Na = parameters[PARAM_lambda_Na];
    const double g_NaL = parameters[PARAM_g_NaL];
    const double thL = parameters[PARAM_thL];
    const double KmKo = parameters[PARAM_KmKo];
    const double KmNaip = parameters[PARAM_KmNaip];
    const double g_NaK = parameters[PARAM_g_NaK];
    const double epi = parameters[PARAM_epi];
    const double g_Ks = parameters[PARAM_g_Ks];
    const double pNaK = parameters[PARAM_pNaK];
    const double g_KATP = parameters[PARAM_g_KATP];
    const double g_Kr = parameters[PARAM_g_Kr];
    const double lambda_K = parameters[PARAM_lambda_K];
    const double g_to = parameters[PARAM_g_to];
    const double g_K1 = parameters[PARAM_g_K1];
    const double g_bCl = parameters[PARAM_g_bCl];
    const double Q10CaL = parameters[PARAM_Q10CaL];
    const double g_CaL = parameters[PARAM_g_CaL];
    const double Kdact = parameters[PARAM_Kdact];
    const double KmCai = parameters[PARAM_KmCai];
    const double KmCao = parameters[PARAM_KmCao];
    const double KmNai = parameters[PARAM_KmNai];
    const double KmNao = parameters[PARAM_KmNao];
    const double Q10NCX = parameters[PARAM_Q10NCX];
    const double g_NaCa = parameters[PARAM_g_NaCa];
    const double ksat = parameters[PARAM_ksat];
    const double nu = parameters[PARAM_nu];
    const double KmPCa = parameters[PARAM_KmPCa];
    const double Q10SLCaP = parameters[PARAM_Q10SLCaP];
    const double g_pCa = parameters[PARAM_g_pCa];
    const double g_bCa = parameters[PARAM_g_bCa];
    const double E_f = parameters[PARAM_E_f];
    const double g_f = parameters[PARAM_g_f];
    const double Nao = parameters[PARAM_Nao];
    const double K_i = parameters[PARAM_K_i];
    const double Ko = parameters[PARAM_Ko];
    const double ce = parameters[PARAM_ce];
    const double K_RyR = parameters[PARAM_K_RyR];
    const double alpha_RyR = parameters[PARAM_alpha_RyR];
    const double beta_RyR = parameters[PARAM_beta_RyR];
    const double eta_RyR = parameters[PARAM_eta_RyR];
    const double gamma_RyR = parameters[PARAM_gamma_RyR];
    const double lambda_RyR = parameters[PARAM_lambda_RyR];
    const double Vc = parameters[PARAM_Vc];
    const double Vd = parameters[PARAM_Vd];
    const double Vn = parameters[PARAM_Vn];
    const double Vs = parameters[PARAM_Vs];
    const double Vsl = parameters[PARAM_Vsl];
    const double J_SERCA_bar = parameters[PARAM_J_SERCA_bar];
    const double K_c = parameters[PARAM_K_c];
    const double K_n = parameters[PARAM_K_n];
    const double B_tot_c = parameters[PARAM_B_tot_c];
    const double B_tot_d = parameters[PARAM_B_tot_d];
    const double B_tot_s = parameters[PARAM_B_tot_s];
    const double B_tot_sl = parameters[PARAM_B_tot_sl];
    const double k_off_c = parameters[PARAM_k_off_c];
    const double k_off_d = parameters[PARAM_k_off_d];
    const double k_off_s = parameters[PARAM_k_off_s];
    const double k_off_sl = parameters[PARAM_k_off_sl];
    const double k_on_c = parameters[PARAM_k_on_c];
    const double k_on_d = parameters[PARAM_k_on_d];
    const double k_on_s = parameters[PARAM_k_on_s];
    const double k_on_sl = parameters[PARAM_k_on_sl];
    const double lambda_B = parameters[PARAM_lambda_B];
    const double lambda_B_c = parameters[PARAM_lambda_B_c];
    const double alpha_d_c = parameters[PARAM_alpha_d_c];
    const double alpha_n_s = parameters[PARAM_alpha_n_s];
    const double alpha_sl_c = parameters[PARAM_alpha_sl_c];
    const double lambda_c_d = parameters[PARAM_lambda_c_d];
    const double lambda_c_i = parameters[PARAM_lambda_c_i];
    const double lambda_diff = parameters[PARAM_lambda_diff];
    const double Cli = parameters[PARAM_Cli];
    const double Clo = parameters[PARAM_Clo];
    const double Cm = parameters[PARAM_Cm];
    const double Frdy = parameters[PARAM_Frdy];
    const double R = parameters[PARAM_R];
    const double Temp = parameters[PARAM_Temp];
    const double chi = parameters[PARAM_chi];
    const double lambda_c_e = parameters[PARAM_lambda_c_e];
    const double stim_amplitude = parameters[PARAM_stim_amplitude];
    const double stim_duration = parameters[PARAM_stim_duration];
    const double stim_period = parameters[PARAM_stim_period];
    const double stim_start = parameters[PARAM_stim_start];

    // Expressions for the Reversal potentials component
    const double FoRT = Frdy / (R * Temp);
    const double ena = log(Nao / Na_i) / FoRT;
    const double ek = log(Ko / K_i) / FoRT;
    const double eca_sl = log(ce / csl) / (2. * FoRT);
    const double ecl = log(Cli / Clo) / FoRT;
    const double Qpow = -31. + Temp / 10.;

    // Expressions for the I_Na component
    const double mss = 1.0
                       / ((1. + 0.00177610354573438 * exp(-V_m / 9.))
                          * (1. + 0.00177610354573438 * exp(-V_m / 9.)));
    const double taum =
            0.06 * exp(-((-0.0980392156862745 + V_m / 51.) * (-0.0980392156862745 + V_m / 51.)))
            + 0.13 * exp(-((2.875 + V_m / 16.) * (2.875 + V_m / 16.)));
    const double aj =
            (V_m >= -38. ? 0.
                         : (38. + V_m) * (-25000.0 * exp(0.2 * V_m) - 7.0e-6 * exp(-0.04 * V_m))
                                   / (1. + 19623624323.6513 * exp(0.3 * V_m)));
    const double bj = (V_m >= -38.218 ? 0.6 * exp(0.09 * V_m) / (1. + exp(-40. - V_m))
                                      : 0.02 * exp(-0.01 * V_m)
                                                / (1. + 0.00369786371648293 * exp(-0.14 * V_m)));
    const double tauj = 1.0 / (aj + bj);
    const double jss =
            1.0
            / ((1. + 29310.8866998062 * exp(V_m / 7.)) * (1. + 29310.8866998062 * exp(V_m / 7.)));
    const double I_Na = g_Na * lambda_Na * (m * m * m) * (-ena + V_m) * j;
    const double dm_dt = (-m + mss) / taum;
    const double dm_dt_linearized = -1. / taum;
    states[STATE_m] = (fabs(dm_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized
                               : dt * dm_dt)
                      + m;
    const double dj_dt = (-j + jss) / tauj;
    const double dj_dt_linearized = -1. / tauj;
    states[STATE_j] = (fabs(dj_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized
                               : dt * dj_dt)
                      + j;

    // Expressions for the I_NaL component
    const double mLss = 1.0 / (1. + 0.000184105793667579 * exp(-V_m / 5.));
    const double tm =
            1.0 / (9.58097877207697 * exp(V_m / 35.) + 2.29642676604141e-5 * exp(-V_m / 6.));
    const double tmL = tm;
    const double hLss = 1.0 / (1. + 124658.506952 * exp(2. * V_m / 15.));
    const double GNaL = (epi == 1. ? g_NaL : 0.6 * g_NaL);
    const double I_NaL = lambda_Na * (-ena + V_m) * GNaL * hL * mL;
    const double dmL_dt = (-mL + mLss) / tmL;
    const double dmL_dt_linearized = -1. / tmL;
    states[STATE_mL] = (fabs(dmL_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dmL_dt_linearized)) * dmL_dt / dmL_dt_linearized
                                : dt * dmL_dt)
                       + mL;
    const double dhL_dt = (-hL + hLss) / thL;
    const double dhL_dt_linearized = -1. / thL;
    states[STATE_hL] = (fabs(dhL_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dhL_dt_linearized)) * dhL_dt / dhL_dt_linearized
                                : dt * dhL_dt)
                       + hL;

    // Expressions for the I_NaK component
    const double sigma = -0.142857142857143 + exp(Nao / 67.) / 7.;
    const double fNaK =
            1.0 / (1. + 0.12 * exp(-0.1 * FoRT * V_m) + 0.037 * exp(-FoRT * V_m) * sigma);
    const double I_NaK = Ko * g_NaK * fNaK / ((1. + pow(KmNaip, 4.) / pow(Na_i, 4.)) * (KmKo + Ko));

    // Expressions for the I_Kr component
    const double Xr1_inf = 1.0 / (1.0 + 0.0146327985189294 * exp(-0.204081632653061 * V_m));
    const double alpha_Xr1 = 450.0 / (1.0 + 0.0111089965382 * exp(-0.1 * V_m));
    const double beta_Xr1 = 6.0 / (1.0 + 13.5813245226 * exp(0.0869565217391 * V_m));
    const double tau_Xr1 = 1.0 * alpha_Xr1 * beta_Xr1;
    const double Xr2_infinity = 1.0 / (1.0 + exp(44. / 25. + V_m / 50.));
    const double alpha_Xr2 = 3.0 / (1.0 + exp(-3. - V_m / 20.));
    const double beta_Xr2 = 1.12 / (1.0 + exp(-3. + V_m / 20.));
    const double tau_Xr2 = 1.0 * alpha_Xr2 * beta_Xr2;
    const double I_Kr = 0.430331482912 * g_Kr * sqrt(Ko) * (-ek + V_m) * Xr1 * Xr2;
    const double dXr1_dt = (-Xr1 + Xr1_inf) / tau_Xr1;
    const double dXr1_dt_linearized = -1. / tau_Xr1;
    states[STATE_Xr1] = (fabs(dXr1_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dXr1_dt_linearized))
                                                                     * dXr1_dt / dXr1_dt_linearized
                                                           : dt * dXr1_dt)
                        + Xr1;
    const double dXr2_dt = (-Xr2 + Xr2_infinity) / tau_Xr2;
    const double dXr2_dt_linearized = -1. / tau_Xr2;
    states[STATE_Xr2] = (fabs(dXr2_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dXr2_dt_linearized))
                                                                     * dXr2_dt / dXr2_dt_linearized
                                                           : dt * dXr2_dt)
                        + Xr2;

    // Expressions for the I_Ks component
    const double eks = log((Ko + Nao * pNaK) / (K_i + pNaK * Na_i)) / FoRT;
    const double xsss = 1.0 / (1. + 0.76228973079 * exp(-V_m / 14.));
    const double tauxs = 990. / (1. + 0.842460441617 * exp(-V_m / 14.));
    const double I_Ks = g_Ks * (x_Ks * x_Ks) * (-eks + V_m);
    const double dx_Ks_dt = (-x_Ks + xsss) / tauxs;
    const double dx_Ks_dt_linearized = -1. / tauxs;
    states[STATE_x_Ks] =
            (fabs(dx_Ks_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dx_Ks_dt_linearized)) * dx_Ks_dt / dx_Ks_dt_linearized
                     : dt * dx_Ks_dt)
            + x_Ks;

    // Expressions for the i_to component
    const double q_inf = 1.0 / (1.0 + 58.9637634804 * exp(0.0769230769231 * V_m));
    const double tau_q =
            6. + 39. / (0.0168716780457 * exp(-0.08 * V_m) + 6.46648051673 * exp(0.1 * V_m));
    const double r_inf = 1.0 / (1.0 + 3.28489055021 * exp(-0.0533333333333 * V_m));
    const double tau_r =
            2.75 + 14.4 / (0.0207698622486 * exp(-0.12 * V_m) + 15.7194688773 * exp(0.09 * V_m));
    const double I_to = g_to * (-ek + V_m) * q * r;
    const double dq_dt = (-q + q_inf) / tau_q;
    const double dq_dt_linearized = -1. / tau_q;
    states[STATE_q] = (fabs(dq_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dq_dt_linearized)) * dq_dt / dq_dt_linearized
                               : dt * dq_dt)
                      + q;
    const double dr_dt = (-r + r_inf) / tau_r;
    const double dr_dt_linearized = -1. / tau_r;
    states[STATE_r] = (fabs(dr_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dr_dt_linearized)) * dr_dt / dr_dt_linearized
                               : dt * dr_dt)
                      + r;

    // Expressions for the I_K1 component
    const double aK1 = 1.0 / (1. + 7.50455791508e-6 * exp(0.2 * V_m - 0.2 * ek));
    const double bK1 = (0.745912348821 * exp(0.08 * V_m - 0.08 * ek)
                        + 3.32464030033e-16 * exp(0.06 * V_m - 0.06 * ek))
                       / (1. + 0.0820849986239 * exp(0.5 * ek - 0.5 * V_m));
    const double K1ss = aK1 / (aK1 + bK1);
    const double I_K1 = 0.430331482912 * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m) * K1ss;

    // Expressions for the I_bCl component
    const double I_bCl = g_bCl * (-ecl + V_m);

    // Expressions for the I_Ca component
    const double fss = 1.0 / (1. + 48.8565712749872 * exp(V_m / 9.))
                       + 0.6 / (1. + 12.1824939607035 * exp(-V_m / 20.));
    const double dss = 1.0 / (1. + exp(-5. / 6. - V_m / 6.));
    const double taud = (fabs(5. + V_m) < 0.02 ? 2.38095238095238
                                               : (1. - 0.434598208507078 * exp(-V_m / 6.)) * dss
                                                         / (0.175 + 0.035 * V_m));
    const double tauf = 1.0 / (0.02 + 0.02 * exp(-((0.493 + 0.034 * V_m) * (0.493 + 0.034 * V_m))));
    const double ibarca_j = 4. * Frdy * g_CaL * (-0.34 * ce + 0.34 * cd * exp(2. * FoRT * V_m))
                            * FoRT * V_m / (-1. + exp(2. * FoRT * V_m));
    const double I_CaL = lambda_c_d * pow(Q10CaL, Qpow) * (1. - f_Ca_B) * d * f * ibarca_j;
    const double dd_dt = (-d + dss) / taud;
    const double dd_dt_linearized = -1. / taud;
    states[STATE_d] = (fabs(dd_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized
                               : dt * dd_dt)
                      + d;
    const double df_dt = (-f + fss) / tauf;
    const double df_dt_linearized = -1. / tauf;
    states[STATE_f] = (fabs(df_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * df_dt_linearized)) * df_dt / df_dt_linearized
                               : dt * df_dt)
                      + f;
    const double df_Ca_B_dt = -0.012 * f_Ca_B + (1.7 - 1.7 * f_Ca_B) * cd;
    const double df_Ca_B_dt_linearized = -0.012 - 1.7 * cd;
    states[STATE_f_Ca_B] =
            (fabs(df_Ca_B_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * df_Ca_B_dt_linearized)) * df_Ca_B_dt / df_Ca_B_dt_linearized
                     : dt * df_Ca_B_dt)
            + f_Ca_B;

    // Expressions for the I_NCX component
    const double Ka_sl = 1.0 / (1. + (Kdact * Kdact) / (csl * csl));
    const double s1_sl = ce * (Na_i * Na_i * Na_i) * exp(nu * FoRT * V_m);
    const double s2_sl = (Nao * Nao * Nao) * csl * exp((-1. + nu) * FoRT * V_m);
    const double s3_sl =
            KmCao * (Na_i * Na_i * Na_i) + ce * (Na_i * Na_i * Na_i) + (Nao * Nao * Nao) * csl
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_i * Na_i * Na_i) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + csl / KmCai) * csl;
    const double I_NaCa = g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl
                          / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);

    // Expressions for the I_PCa component
    const double I_pCa = g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl)
                         / ((KmPCa * KmPCa) + (csl * csl));

    // Expressions for the I_CaBK component
    const double I_bCa = g_bCa * lambda_c_e * (-eca_sl + V_m);

    // Expressions for the I_f component
    const double xf_inf = 1.0 / (1. + 5956538.01318461 * exp(V_m / 5.));
    const double tau_xf = 1900. / (1. + 4.48168907033806 * exp(V_m / 10.));
    const double I_f = g_f * (-E_f + V_m) * xf;
    const double dxf_dt = (-xf + xf_inf) / tau_xf;
    const double dxf_dt_linearized = -1. / tau_xf;
    states[STATE_xf] = (fabs(dxf_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dxf_dt_linearized)) * dxf_dt / dxf_dt_linearized
                                : dt * dxf_dt)
                       + xf;

    // Expressions for the I_KATP component
    const double I_KATP =
            0.60295079490657 * g_KATP * pow(Ko, 0.3) * (-ek + V_m) / (40. + 0.0875 * V_m);

    // Expressions for the Ca Fluxes component
    const double J_CaL = -Cm * chi * I_CaL / (2. * Frdy);
    const double J_pCa = -Cm * chi * I_pCa / (2. * Frdy);
    const double J_bCa = -Cm * chi * I_bCa / (2. * Frdy);
    const double J_NaCa = Cm * chi * I_NaCa / Frdy;
    const double J_e_sl = J_NaCa + J_bCa + J_pCa;
    const double Q10SERCA = 2.60000000000000;
    const double J_SERCA = J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                           * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n))
                           / (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n));
    const double J_n_s = alpha_n_s * lambda_c_i * lambda_diff * (-cs + cn);
    const double J_sl_c = alpha_sl_c * lambda_c_i * lambda_diff * (-cc + csl);
    const double J_d_c = alpha_d_c * lambda_c_d * lambda_diff * (-cc + cd);

    // Expressions for the RyRs component
    const double p = 1.0 / (1. + (K_RyR * K_RyR * K_RyR) / (cd * cd * cd));
    const double J_RyR_active = alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p * r_RyR;
    const double J_leak = alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i * (-csl + cs);
    const double J_RyR = J_RyR_active + J_leak;
    const double dr_RyR_dt =
            eta_RyR * (1. - r_RyR) / p - J_RyR_active / (beta_RyR * lambda_RyR * lambda_c_i);
    const double dJ_RyR_active_dr_RyR = alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p;
    const double dr_RyR_dt_linearized =
            -eta_RyR / p - dJ_RyR_active_dr_RyR / (beta_RyR * lambda_RyR * lambda_c_i);
    states[STATE_r_RyR] =
            (fabs(dr_RyR_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dr_RyR_dt_linearized)) * dr_RyR_dt / dr_RyR_dt_linearized
                     : dt * dr_RyR_dt)
            + r_RyR;

    // Expressions for the Ca Buffers component
    const double J_c_b =
            Vc * (-k_off_c * bc + k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c) * cc);
    const double J_d_b =
            Vd * (-k_off_d * bd + k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c) * cd);
    const double J_s_b = Vs * (-k_off_s * bs + k_on_s * (-bs + B_tot_s * lambda_B) * cs);
    const double J_sl_b =
            Vsl * (-k_off_sl * bsl + k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c) * csl);

    // Expressions for the Ca Concentrations component
    const double dcn_dt = (1.0 * J_SERCA - 1.0 * J_n_s) / Vn;
    const double dJ_SERCA_dcn =
            -2. * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow) * cn
                    / ((K_n * K_n) * (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n)))
            - 2. * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                      * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n)) * cn
                      / ((K_n * K_n)
                         * ((1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n))
                            * (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n))));
    const double dJ_n_s_dcn = alpha_n_s * lambda_c_i * lambda_diff;
    const double dcn_dt_linearized = (1.0 * dJ_SERCA_dcn - 1.0 * dJ_n_s_dcn) / Vn;
    states[STATE_cn] = (fabs(dcn_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dcn_dt_linearized)) * dcn_dt / dcn_dt_linearized
                                : dt * dcn_dt)
                       + cn;
    const double dcc_dt = (1.0 * J_d_c + 1.0 * J_sl_c - 1.0 * J_SERCA - 1.0 * J_c_b) / Vc;
    const double dJ_SERCA_dcc =
            2. * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow) * cc
                    / ((K_c * K_c) * (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n)))
            - 2. * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                      * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n)) * cc
                      / ((K_c * K_c)
                         * ((1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n))
                            * (1. + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n))));
    const double dJ_c_b_dcc = Vc * k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c);
    const double dJ_d_c_dcc = -alpha_d_c * lambda_c_d * lambda_diff;
    const double dJ_sl_c_dcc = -alpha_sl_c * lambda_c_i * lambda_diff;
    const double dcc_dt_linearized =
            (1.0 * dJ_d_c_dcc + 1.0 * dJ_sl_c_dcc - 1.0 * dJ_SERCA_dcc - 1.0 * dJ_c_b_dcc) / Vc;
    states[STATE_cc] = (fabs(dcc_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dcc_dt_linearized)) * dcc_dt / dcc_dt_linearized
                                : dt * dcc_dt)
                       + cc;
    const double dcd_dt = (1.0 * J_CaL - 1.0 * J_d_b - 1.0 * J_d_c) / Vd;
    const double dI_CaL_dibarca_j = lambda_c_d * pow(Q10CaL, Qpow) * (1. - f_Ca_B) * d * f;
    const double dJ_CaL_dI_CaL = -Cm * chi / (2. * Frdy);
    const double dJ_d_b_dcd = Vd * k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c);
    const double dJ_d_c_dcd = alpha_d_c * lambda_c_d * lambda_diff;
    const double dibarca_j_dcd =
            1.36 * Frdy * g_CaL * FoRT * V_m * exp(2. * FoRT * V_m) / (-1. + exp(2. * FoRT * V_m));
    const double dcd_dt_linearized = (-1.0 * dJ_d_b_dcd - 1.0 * dJ_d_c_dcd
                                      + 1.0 * dI_CaL_dibarca_j * dJ_CaL_dI_CaL * dibarca_j_dcd)
                                     / Vd;
    states[STATE_cd] = (fabs(dcd_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dcd_dt_linearized)) * dcd_dt / dcd_dt_linearized
                                : dt * dcd_dt)
                       + cd;
    const double dcsl_dt = (1.0 * J_RyR + 1.0 * J_e_sl - 1.0 * J_sl_b - 1.0 * J_sl_c) / Vsl;
    const double dI_NaCa_dKa_sl = g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl)
                                  / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_NaCa_ds2_sl = -g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * Ka_sl
                                  / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_NaCa_ds3_sl = -g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl)
                                  * Ka_sl
                                  / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * (s3_sl * s3_sl));
    const double dI_bCa_deca_sl = -g_bCa * lambda_c_e;
    const double dI_pCa_dcsl =
            -2. * g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl * csl)
                    / (((KmPCa * KmPCa) + (csl * csl)) * ((KmPCa * KmPCa) + (csl * csl)))
            + 2. * g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * csl / ((KmPCa * KmPCa) + (csl * csl));
    const double dJ_NaCa_dI_NaCa = Cm * chi / Frdy;
    const double dJ_RyR_active_dcsl = -alpha_RyR * lambda_RyR * lambda_c_i * p * r_RyR;
    const double dJ_bCa_dI_bCa = -Cm * chi / (2. * Frdy);
    const double dJ_leak_dcsl = -alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i;
    const double dJ_pCa_dI_pCa = -Cm * chi / (2. * Frdy);
    const double dJ_sl_b_dcsl = Vsl * k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c);
    const double dJ_sl_c_dcsl = alpha_sl_c * lambda_c_i * lambda_diff;
    const double dKa_sl_dcsl =
            2.0 * (Kdact * Kdact)
            / (((1. + (Kdact * Kdact) / (csl * csl)) * (1. + (Kdact * Kdact) / (csl * csl)))
               * (csl * csl * csl));
    const double deca_sl_dcsl = -1. / (2. * FoRT * csl);
    const double ds2_sl_dcsl = (Nao * Nao * Nao) * exp((-1. + nu) * FoRT * V_m);
    const double ds3_sl_dcsl = (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (1. + csl / KmCai)
                               + (KmNao * KmNao * KmNao) * csl / KmCai;
    const double dcsl_dt_linearized =
            (1.0 * dJ_RyR_active_dcsl + 1.0 * dJ_leak_dcsl - 1.0 * dJ_sl_b_dcsl - 1.0 * dJ_sl_c_dcsl
             + 1.0
                       * (dI_NaCa_dKa_sl * dKa_sl_dcsl + dI_NaCa_ds2_sl * ds2_sl_dcsl
                          + dI_NaCa_ds3_sl * ds3_sl_dcsl)
                       * dJ_NaCa_dI_NaCa
             + 1.0 * dI_pCa_dcsl * dJ_pCa_dI_pCa
             + 1.0 * dI_bCa_deca_sl * dJ_bCa_dI_bCa * deca_sl_dcsl)
            / Vsl;
    states[STATE_csl] = (fabs(dcsl_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dcsl_dt_linearized))
                                                                     * dcsl_dt / dcsl_dt_linearized
                                                           : dt * dcsl_dt)
                        + csl;
    const double dcs_dt = (1.0 * J_n_s - 1.0 * J_RyR - 1.0 * J_s_b) / Vs;
    const double dJ_RyR_active_dcs = alpha_RyR * lambda_RyR * lambda_c_i * p * r_RyR;
    const double dJ_leak_dcs = alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i;
    const double dJ_n_s_dcs = -alpha_n_s * lambda_c_i * lambda_diff;
    const double dJ_s_b_dcs = Vs * k_on_s * (-bs + B_tot_s * lambda_B);
    const double dcs_dt_linearized =
            (1.0 * dJ_n_s_dcs - 1.0 * dJ_RyR_active_dcs - 1.0 * dJ_leak_dcs - 1.0 * dJ_s_b_dcs)
            / Vs;
    states[STATE_cs] = (fabs(dcs_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dcs_dt_linearized)) * dcs_dt / dcs_dt_linearized
                                : dt * dcs_dt)
                       + cs;

    // Expressions for the Ca Buffer Concentrations component
    const double dbc_dt = 1.0 * J_c_b / Vc;
    const double dJ_c_b_dbc = Vc * (-k_off_c - k_on_c * cc);
    const double dbc_dt_linearized = 1.0 * dJ_c_b_dbc / Vc;
    states[STATE_bc] = (fabs(dbc_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dbc_dt_linearized)) * dbc_dt / dbc_dt_linearized
                                : dt * dbc_dt)
                       + bc;
    const double dbd_dt = 1.0 * J_d_b / Vd;
    const double dJ_d_b_dbd = Vd * (-k_off_d - k_on_d * cd);
    const double dbd_dt_linearized = 1.0 * dJ_d_b_dbd / Vd;
    states[STATE_bd] = (fabs(dbd_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dbd_dt_linearized)) * dbd_dt / dbd_dt_linearized
                                : dt * dbd_dt)
                       + bd;
    const double dbs_dt = 1.0 * J_s_b / Vs;
    const double dJ_s_b_dbs = Vs * (-k_off_s - k_on_s * cs);
    const double dbs_dt_linearized = 1.0 * dJ_s_b_dbs / Vs;
    states[STATE_bs] = (fabs(dbs_dt_linearized) > 1.0e-8
                                ? (-1.0 + exp(dt * dbs_dt_linearized)) * dbs_dt / dbs_dt_linearized
                                : dt * dbs_dt)
                       + bs;
    const double dbsl_dt = 1.0 * J_sl_b / Vsl;
    const double dJ_sl_b_dbsl = Vsl * (-k_off_sl - k_on_sl * csl);
    const double dbsl_dt_linearized = 1.0 * dJ_sl_b_dbsl / Vsl;
    states[STATE_bsl] = (fabs(dbsl_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dbsl_dt_linearized))
                                                                     * dbsl_dt / dbsl_dt_linearized
                                                           : dt * dbsl_dt)
                        + bsl;

    // Expressions for the Membrane potential component
    const double i_Stim =
            (V_m < -40. ? 1. : 0.)
            * (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                               && t - stim_period * floor(t / stim_period) >= stim_start
                       ? -stim_amplitude
                       : 0.);
    const double I_tot = I_CaL + I_K1 + I_KATP + I_Kr + I_Ks + I_Na + I_NaCa + I_NaK + I_NaL + I_bCa
                         + I_bCl + I_f + I_pCa + I_to;
    const double dV_m_dt = -I_tot - i_Stim;
    const double dI_K1_dK1ss = 0.430331482912 * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m);
    const double dK1ss_daK1 = 1.0 / (aK1 + bK1) - aK1 / ((aK1 + bK1) * (aK1 + bK1));
    const double dK1ss_dbK1 = -aK1 / ((aK1 + bK1) * (aK1 + bK1));
    const double daK1_dV_m = -1.500911583016e-6 * exp(0.2 * V_m - 0.2 * ek)
                             / ((1. + 7.50455791508e-6 * exp(0.2 * V_m - 0.2 * ek))
                                * (1. + 7.50455791508e-6 * exp(0.2 * V_m - 0.2 * ek)));
    const double dbK1_dV_m = (1.994784180198e-17 * exp(0.06 * V_m - 0.06 * ek)
                              + 0.05967298790568 * exp(0.08 * V_m - 0.08 * ek))
                                     / (1. + 0.0820849986239 * exp(0.5 * ek - 0.5 * V_m))
                             + 0.04104249931195
                                       * (0.745912348821 * exp(0.08 * V_m - 0.08 * ek)
                                          + 3.32464030033e-16 * exp(0.06 * V_m - 0.06 * ek))
                                       * exp(0.5 * ek - 0.5 * V_m)
                                       / ((1. + 0.0820849986239 * exp(0.5 * ek - 0.5 * V_m))
                                          * (1. + 0.0820849986239 * exp(0.5 * ek - 0.5 * V_m)));
    const double dI_K1_dV_m = 0.430331482912 * g_K1 * lambda_K * sqrt(Ko) * K1ss
                              + 0.430331482912 * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m)
                                        * (dK1ss_daK1 * daK1_dV_m + dK1ss_dbK1 * dbK1_dV_m);
    const double dI_KATP_dV_m = 0.60295079490657 * g_KATP * pow(Ko, 0.3) / (40. + 0.0875 * V_m)
                                - 0.0527581945543249 * g_KATP * pow(Ko, 0.3) * (-ek + V_m)
                                          / ((40. + 0.0875 * V_m) * (40. + 0.0875 * V_m));
    const double dI_Kr_dV_m = 0.430331482912 * g_Kr * sqrt(Ko) * Xr1 * Xr2;
    const double dI_Ks_dV_m = g_Ks * (x_Ks * x_Ks);
    const double dI_Na_dV_m = g_Na * lambda_Na * (m * m * m) * j;
    const double ds1_sl_dV_m = ce * nu * (Na_i * Na_i * Na_i) * FoRT * exp(nu * FoRT * V_m);
    const double ds2_sl_dV_m =
            (Nao * Nao * Nao) * (-1. + nu) * FoRT * csl * exp((-1. + nu) * FoRT * V_m);
    const double dI_NaCa_dV_m =
            g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-ds2_sl_dV_m + ds1_sl_dV_m) * Ka_sl
                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl)
            - g_NaCa * ksat * lambda_c_e * pow(Q10NCX, Qpow) * (-1. + nu) * (-s2_sl + s1_sl) * FoRT
                      * Ka_sl * exp((-1. + nu) * FoRT * V_m)
                      / (((1. + ksat * exp((-1. + nu) * FoRT * V_m))
                          * (1. + ksat * exp((-1. + nu) * FoRT * V_m)))
                         * s3_sl);
    const double dI_NaCa_ds1_sl = g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * Ka_sl
                                  / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_NaK_dfNaK = Ko * g_NaK / ((1. + pow(KmNaip, 4.) / pow(Na_i, 4.)) * (KmKo + Ko));
    const double dI_NaL_dV_m = lambda_Na * GNaL * hL * mL;
    const double dI_to_dV_m = g_to * q * r;
    const double dfNaK_dV_m =
            1.0 * (0.012 * FoRT * exp(-0.1 * FoRT * V_m) + 0.037 * FoRT * exp(-FoRT * V_m) * sigma)
            / ((1. + 0.12 * exp(-0.1 * FoRT * V_m) + 0.037 * exp(-FoRT * V_m) * sigma)
               * (1. + 0.12 * exp(-0.1 * FoRT * V_m) + 0.037 * exp(-FoRT * V_m) * sigma));
    const double dibarca_j_dV_m =
            4. * Frdy * g_CaL * (-0.34 * ce + 0.34 * cd * exp(2. * FoRT * V_m)) * FoRT
                    / (-1. + exp(2. * FoRT * V_m))
            - 8. * Frdy * g_CaL * (FoRT * FoRT) * (-0.34 * ce + 0.34 * cd * exp(2. * FoRT * V_m))
                      * V_m * exp(2. * FoRT * V_m)
                      / ((-1. + exp(2. * FoRT * V_m)) * (-1. + exp(2. * FoRT * V_m)))
            + 2.72 * Frdy * g_CaL * (FoRT * FoRT) * V_m * cd * exp(2. * FoRT * V_m)
                      / (-1. + exp(2. * FoRT * V_m));
    const double dV_m_dt_linearized =
            -g_bCl - dI_K1_dV_m - dI_KATP_dV_m - dI_Kr_dV_m - dI_Ks_dV_m - dI_NaCa_dV_m
            - dI_NaL_dV_m - dI_Na_dV_m - dI_to_dV_m - g_bCa * lambda_c_e - g_f * xf
            - (dK1ss_daK1 * daK1_dV_m + dK1ss_dbK1 * dbK1_dV_m) * dI_K1_dK1ss
            - dI_CaL_dibarca_j * dibarca_j_dV_m - dI_NaCa_ds1_sl * ds1_sl_dV_m
            - dI_NaCa_ds2_sl * ds2_sl_dV_m - dI_NaK_dfNaK * dfNaK_dV_m;
    states[STATE_V_m] = (fabs(dV_m_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dV_m_dt_linearized))
                                                                     * dV_m_dt / dV_m_dt_linearized
                                                           : dt * dV_m_dt)
                        + V_m;

    // Expressions for the Sodium concentration component
    const double I_Na_tot = 3. * I_NaCa + 3. * I_NaK + 0.3293 * I_f + I_Na + I_NaL;
    const double J_Na = -Cm * chi * I_Na_tot / Frdy;
    const double dNa_i_dt = J_Na;
    const double dI_Na_dena = -g_Na * lambda_Na * (m * m * m) * j;
    const double dI_NaK_dNa_i =
            4. * Ko * g_NaK * pow(KmNaip, 4.) * fNaK
            / (((1. + pow(KmNaip, 4.) / pow(Na_i, 4.)) * (1. + pow(KmNaip, 4.) / pow(Na_i, 4.)))
               * (KmKo + Ko) * pow(Na_i, 5.));
    const double dI_NaL_dena = -lambda_Na * GNaL * hL * mL;
    const double dJ_Na_dI_Na_tot = -Cm * chi / Frdy;
    const double dena_dNa_i = -1. / (FoRT * Na_i);
    const double ds1_sl_dNa_i = 3. * ce * (Na_i * Na_i) * exp(nu * FoRT * V_m);
    const double ds3_sl_dNa_i =
            3. * KmCao * (Na_i * Na_i) + 3. * ce * (Na_i * Na_i)
            + 3. * KmCai * (Nao * Nao * Nao) * (Na_i * Na_i) / (KmNai * KmNai * KmNai);
    const double dNa_i_dt_linearized =
            (3. * dI_NaK_dNa_i + dI_NaL_dena * dena_dNa_i + dI_Na_dena * dena_dNa_i
             + 3. * dI_NaCa_ds1_sl * ds1_sl_dNa_i + 3. * dI_NaCa_ds3_sl * ds3_sl_dNa_i)
            * dJ_Na_dI_Na_tot;
    states[STATE_Na_i] =
            Na_i
            + (fabs(dNa_i_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dNa_i_dt_linearized)) * dNa_i_dt / dNa_i_dt_linearized
                       : dt * dNa_i_dt);
}
