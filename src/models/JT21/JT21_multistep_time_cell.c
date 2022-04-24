#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "JT21" model

#include "JT21_simd.h"
#include "cellmodel.h"

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

// Init state values
static void init_state_values(cellmodel_float_t *__restrict states, const long num_cells,
                              const long padded_num_cells)
{
#pragma omp parallel for
    for (long i = 0; i < num_cells; i++) {
        states[STATE_m * padded_num_cells + i] = FP_LITERAL(0.0047593841564);
        states[STATE_j * padded_num_cells + i] = FP_LITERAL(0.59186803670583);
        states[STATE_mL * padded_num_cells + i] = FP_LITERAL(0.00056156997976);
        states[STATE_hL * padded_num_cells + i] = FP_LITERAL(0.19825212801003);
        states[STATE_Xr1 * padded_num_cells + i] = FP_LITERAL(0.02598687936732);
        states[STATE_Xr2 * padded_num_cells + i] = FP_LITERAL(0.46218104728376);
        states[STATE_x_Ks * padded_num_cells + i] = FP_LITERAL(0.00418046147173);
        states[STATE_q * padded_num_cells + i] = FP_LITERAL(0.89186079496414);
        states[STATE_r * padded_num_cells + i] = FP_LITERAL(0.00415791788969);
        states[STATE_d * padded_num_cells + i] = FP_LITERAL(3.47393157e-06);
        states[STATE_f * padded_num_cells + i] = FP_LITERAL(0.99328007454758);
        states[STATE_f_Ca_B * padded_num_cells + i] = FP_LITERAL(0.05601726193501);
        states[STATE_xf * padded_num_cells + i] = FP_LITERAL(0.08014973245989);
        states[STATE_r_RyR * padded_num_cells + i] = FP_LITERAL(0.99999973974165);
        states[STATE_cn * padded_num_cells + i] = FP_LITERAL(0.69858720814252);
        states[STATE_cc * padded_num_cells + i] = FP_LITERAL(0.00011382799663);
        states[STATE_cd * padded_num_cells + i] = FP_LITERAL(0.00020977171788);
        states[STATE_csl * padded_num_cells + i] = FP_LITERAL(0.00011798425266);
        states[STATE_cs * padded_num_cells + i] = FP_LITERAL(0.69792301896013);
        states[STATE_bc * padded_num_cells + i] = FP_LITERAL(0.00839395214247);
        states[STATE_bd * padded_num_cells + i] = FP_LITERAL(0.05570692056887);
        states[STATE_bs * padded_num_cells + i] = FP_LITERAL(31.06660220440676);
        states[STATE_bsl * padded_num_cells + i] = FP_LITERAL(0.10584192890313);
        states[STATE_V_m * padded_num_cells + i] = FP_LITERAL(-80.42101165870085);
        states[STATE_Na_i * padded_num_cells + i] = FP_LITERAL(8.07942377455463);
    }
}

// Default parameter values
static void init_parameters_values(cellmodel_float_t *__restrict parameters)
{
    parameters[PARAM_g_Na] = FP_LITERAL(0.36);
    parameters[PARAM_lambda_Na] = FP_LITERAL(1.0);
    parameters[PARAM_g_NaL] = FP_LITERAL(0.03);
    parameters[PARAM_thL] = FP_LITERAL(200.0);
    parameters[PARAM_KmKo] = FP_LITERAL(1.5);
    parameters[PARAM_KmNaip] = FP_LITERAL(11.0);
    parameters[PARAM_Q10KmNai] = FP_LITERAL(1.4);
    parameters[PARAM_Q10NaK] = FP_LITERAL(1.6);
    parameters[PARAM_g_NaK] = FP_LITERAL(2.3939999999999997);
    parameters[PARAM_epi] = FP_LITERAL(0.0);
    parameters[PARAM_g_Ks] = FP_LITERAL(0.0127);
    parameters[PARAM_pNaK] = FP_LITERAL(0.018);
    parameters[PARAM_g_KATP] = FP_LITERAL(0.01);
    parameters[PARAM_g_Kr] = FP_LITERAL(0.029399999999999996);
    parameters[PARAM_lambda_K] = FP_LITERAL(1.0);
    parameters[PARAM_g_to] = FP_LITERAL(0.105);
    parameters[PARAM_g_K1] = FP_LITERAL(0.27);
    parameters[PARAM_g_bCl] = FP_LITERAL(0.003);
    parameters[PARAM_Q10CaL] = FP_LITERAL(1.8);
    parameters[PARAM_g_CaL] = FP_LITERAL(0.8316000000000001);
    parameters[PARAM_Kdact] = FP_LITERAL(0.00015);
    parameters[PARAM_KmCai] = FP_LITERAL(0.0036);
    parameters[PARAM_KmCao] = FP_LITERAL(1.3);
    parameters[PARAM_KmNai] = FP_LITERAL(12.3);
    parameters[PARAM_KmNao] = FP_LITERAL(87.5);
    parameters[PARAM_Q10NCX] = FP_LITERAL(1.6);
    parameters[PARAM_g_NaCa] = FP_LITERAL(14.1);
    parameters[PARAM_ksat] = FP_LITERAL(0.3);
    parameters[PARAM_nu] = FP_LITERAL(0.3);
    parameters[PARAM_KmPCa] = FP_LITERAL(0.0005);
    parameters[PARAM_Q10SLCaP] = FP_LITERAL(2.35);
    parameters[PARAM_g_pCa] = FP_LITERAL(0.12);
    parameters[PARAM_g_bCa] = FP_LITERAL(0.0021);
    parameters[PARAM_E_f] = FP_LITERAL(-17.0);
    parameters[PARAM_g_f] = FP_LITERAL(0.01);
    parameters[PARAM_Na_sl] = FP_LITERAL(8.0);
    parameters[PARAM_Nao] = FP_LITERAL(140.0);
    parameters[PARAM_K_i] = FP_LITERAL(120.0);
    parameters[PARAM_Ko] = FP_LITERAL(5.0);
    parameters[PARAM_ce] = FP_LITERAL(0.42);
    parameters[PARAM_K_RyR] = FP_LITERAL(0.015);
    parameters[PARAM_alpha_RyR] = FP_LITERAL(0.02075);
    parameters[PARAM_beta_RyR] = FP_LITERAL(0.042);
    parameters[PARAM_eta_RyR] = 1e-05;
    parameters[PARAM_gamma_RyR] = FP_LITERAL(0.001);
    parameters[PARAM_lambda_RyR] = FP_LITERAL(0.63);
    parameters[PARAM_Vc] = FP_LITERAL(0.917);
    parameters[PARAM_Vd] = FP_LITERAL(0.001);
    parameters[PARAM_Vn] = FP_LITERAL(0.05);
    parameters[PARAM_Vs] = FP_LITERAL(0.004);
    parameters[PARAM_Vsl] = FP_LITERAL(0.028);
    parameters[PARAM_J_SERCA_bar] = FP_LITERAL(0.00016);
    parameters[PARAM_K_c] = FP_LITERAL(0.00025);
    parameters[PARAM_K_n] = FP_LITERAL(1.7);
    parameters[PARAM_B_tot_c] = FP_LITERAL(0.063);
    parameters[PARAM_B_tot_d] = FP_LITERAL(2.7);
    parameters[PARAM_B_tot_s] = FP_LITERAL(60.0);
    parameters[PARAM_B_tot_sl] = FP_LITERAL(1.45);
    parameters[PARAM_k_off_c] = FP_LITERAL(0.03);
    parameters[PARAM_k_off_d] = FP_LITERAL(1.0);
    parameters[PARAM_k_off_s] = FP_LITERAL(65.0);
    parameters[PARAM_k_off_sl] = FP_LITERAL(0.15);
    parameters[PARAM_k_on_c] = FP_LITERAL(40.0);
    parameters[PARAM_k_on_d] = FP_LITERAL(100.0);
    parameters[PARAM_k_on_s] = FP_LITERAL(100.0);
    parameters[PARAM_k_on_sl] = FP_LITERAL(100.0);
    parameters[PARAM_lambda_B] = FP_LITERAL(1.0);
    parameters[PARAM_lambda_B_c] = FP_LITERAL(1.0);
    parameters[PARAM_alpha_d_c] = FP_LITERAL(0.0027);
    parameters[PARAM_alpha_n_s] = FP_LITERAL(0.0093);
    parameters[PARAM_alpha_sl_c] = FP_LITERAL(0.3);
    parameters[PARAM_lambda_c_d] = FP_LITERAL(1.0);
    parameters[PARAM_lambda_c_i] = FP_LITERAL(1.0);
    parameters[PARAM_lambda_diff] = FP_LITERAL(1.0);
    parameters[PARAM_Cli] = FP_LITERAL(15.0);
    parameters[PARAM_Clo] = FP_LITERAL(150.0);
    parameters[PARAM_Cm] = FP_LITERAL(0.01);
    parameters[PARAM_Frdy] = FP_LITERAL(96.485);
    parameters[PARAM_R] = FP_LITERAL(8.314);
    parameters[PARAM_Temp] = FP_LITERAL(310.0);
    parameters[PARAM_chi] = FP_LITERAL(0.9);
    parameters[PARAM_lambda_c_e] = FP_LITERAL(1.0);
    parameters[PARAM_stim_amplitude] = FP_LITERAL(5.0);
    parameters[PARAM_stim_duration] = FP_LITERAL(20.0);
    parameters[PARAM_stim_period] = FP_LITERAL(10000.0);
    parameters[PARAM_stim_start] = FP_LITERAL(0.0);
}

// State index
static int state_index(const char name[])
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
static int parameter_index(const char name[])
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


// Compute multiple forward steps using the explicit Euler algorithm to the JT21 ODE
static void multistep_FE(double *__restrict states, const double t_start, const double dt,
                         const double *__restrict parameters, const long num_cells,
                         const long padded_num_cells, const int num_timesteps)
{
    double t = t_start;
    for (int ti = 0; ti < num_timesteps; ti++) {
#if defined(HINT_CLANG_SIMD)
#pragma omp parallel for
#pragma clang loop vectorize(assume_safety)
#elif defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp parallel for simd aligned(states                                                       \
                                      : CELLMODEL_STATES_ALIGNMENT_BYTES) simdlen(VECTOR_LENGTH)
#else
#pragma omp parallel for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES)
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp parallel for
#endif
        for (long i = 0; i < num_cells; i++) {
            // Assign parameters
            const cellmodel_float_t g_Na = parameters[PARAM_g_Na * padded_num_cells + i];
            const cellmodel_float_t lambda_Na = parameters[PARAM_lambda_Na * padded_num_cells + i];
            const cellmodel_float_t g_NaL = parameters[PARAM_g_NaL * padded_num_cells + i];
            const cellmodel_float_t thL = parameters[PARAM_thL * padded_num_cells + i];
            const cellmodel_float_t KmKo = parameters[PARAM_KmKo * padded_num_cells + i];
            const cellmodel_float_t KmNaip = parameters[PARAM_KmNaip * padded_num_cells + i];
            const cellmodel_float_t g_NaK = parameters[PARAM_g_NaK * padded_num_cells + i];
            const cellmodel_float_t epi = parameters[PARAM_epi * padded_num_cells + i];
            const cellmodel_float_t g_Ks = parameters[PARAM_g_Ks * padded_num_cells + i];
            const cellmodel_float_t pNaK = parameters[PARAM_pNaK * padded_num_cells + i];
            const cellmodel_float_t g_KATP = parameters[PARAM_g_KATP * padded_num_cells + i];
            const cellmodel_float_t g_Kr = parameters[PARAM_g_Kr * padded_num_cells + i];
            const cellmodel_float_t lambda_K = parameters[PARAM_lambda_K * padded_num_cells + i];
            const cellmodel_float_t g_to = parameters[PARAM_g_to * padded_num_cells + i];
            const cellmodel_float_t g_K1 = parameters[PARAM_g_K1 * padded_num_cells + i];
            const cellmodel_float_t g_bCl = parameters[PARAM_g_bCl * padded_num_cells + i];
            const cellmodel_float_t Q10CaL = parameters[PARAM_Q10CaL * padded_num_cells + i];
            const cellmodel_float_t g_CaL = parameters[PARAM_g_CaL * padded_num_cells + i];
            const cellmodel_float_t Kdact = parameters[PARAM_Kdact * padded_num_cells + i];
            const cellmodel_float_t KmCai = parameters[PARAM_KmCai * padded_num_cells + i];
            const cellmodel_float_t KmCao = parameters[PARAM_KmCao * padded_num_cells + i];
            const cellmodel_float_t KmNai = parameters[PARAM_KmNai * padded_num_cells + i];
            const cellmodel_float_t KmNao = parameters[PARAM_KmNao * padded_num_cells + i];
            const cellmodel_float_t Q10NCX = parameters[PARAM_Q10NCX * padded_num_cells + i];
            const cellmodel_float_t g_NaCa = parameters[PARAM_g_NaCa * padded_num_cells + i];
            const cellmodel_float_t ksat = parameters[PARAM_ksat * padded_num_cells + i];
            const cellmodel_float_t nu = parameters[PARAM_nu * padded_num_cells + i];
            const cellmodel_float_t KmPCa = parameters[PARAM_KmPCa * padded_num_cells + i];
            const cellmodel_float_t Q10SLCaP = parameters[PARAM_Q10SLCaP * padded_num_cells + i];
            const cellmodel_float_t g_pCa = parameters[PARAM_g_pCa * padded_num_cells + i];
            const cellmodel_float_t g_bCa = parameters[PARAM_g_bCa * padded_num_cells + i];
            const cellmodel_float_t E_f = parameters[PARAM_E_f * padded_num_cells + i];
            const cellmodel_float_t g_f = parameters[PARAM_g_f * padded_num_cells + i];
            const cellmodel_float_t Nao = parameters[PARAM_Nao * padded_num_cells + i];
            const cellmodel_float_t K_i = parameters[PARAM_K_i * padded_num_cells + i];
            const cellmodel_float_t Ko = parameters[PARAM_Ko * padded_num_cells + i];
            const cellmodel_float_t ce = parameters[PARAM_ce * padded_num_cells + i];
            const cellmodel_float_t K_RyR = parameters[PARAM_K_RyR * padded_num_cells + i];
            const cellmodel_float_t alpha_RyR = parameters[PARAM_alpha_RyR * padded_num_cells + i];
            const cellmodel_float_t beta_RyR = parameters[PARAM_beta_RyR * padded_num_cells + i];
            const cellmodel_float_t eta_RyR = parameters[PARAM_eta_RyR * padded_num_cells + i];
            const cellmodel_float_t gamma_RyR = parameters[PARAM_gamma_RyR * padded_num_cells + i];
            const cellmodel_float_t lambda_RyR =
                    parameters[PARAM_lambda_RyR * padded_num_cells + i];
            const cellmodel_float_t Vc = parameters[PARAM_Vc * padded_num_cells + i];
            const cellmodel_float_t Vd = parameters[PARAM_Vd * padded_num_cells + i];
            const cellmodel_float_t Vn = parameters[PARAM_Vn * padded_num_cells + i];
            const cellmodel_float_t Vs = parameters[PARAM_Vs * padded_num_cells + i];
            const cellmodel_float_t Vsl = parameters[PARAM_Vsl * padded_num_cells + i];
            const cellmodel_float_t J_SERCA_bar =
                    parameters[PARAM_J_SERCA_bar * padded_num_cells + i];
            const cellmodel_float_t K_c = parameters[PARAM_K_c * padded_num_cells + i];
            const cellmodel_float_t K_n = parameters[PARAM_K_n * padded_num_cells + i];
            const cellmodel_float_t B_tot_c = parameters[PARAM_B_tot_c * padded_num_cells + i];
            const cellmodel_float_t B_tot_d = parameters[PARAM_B_tot_d * padded_num_cells + i];
            const cellmodel_float_t B_tot_s = parameters[PARAM_B_tot_s * padded_num_cells + i];
            const cellmodel_float_t B_tot_sl = parameters[PARAM_B_tot_sl * padded_num_cells + i];
            const cellmodel_float_t k_off_c = parameters[PARAM_k_off_c * padded_num_cells + i];
            const cellmodel_float_t k_off_d = parameters[PARAM_k_off_d * padded_num_cells + i];
            const cellmodel_float_t k_off_s = parameters[PARAM_k_off_s * padded_num_cells + i];
            const cellmodel_float_t k_off_sl = parameters[PARAM_k_off_sl * padded_num_cells + i];
            const cellmodel_float_t k_on_c = parameters[PARAM_k_on_c * padded_num_cells + i];
            const cellmodel_float_t k_on_d = parameters[PARAM_k_on_d * padded_num_cells + i];
            const cellmodel_float_t k_on_s = parameters[PARAM_k_on_s * padded_num_cells + i];
            const cellmodel_float_t k_on_sl = parameters[PARAM_k_on_sl * padded_num_cells + i];
            const cellmodel_float_t lambda_B = parameters[PARAM_lambda_B * padded_num_cells + i];
            const cellmodel_float_t lambda_B_c =
                    parameters[PARAM_lambda_B_c * padded_num_cells + i];
            const cellmodel_float_t alpha_d_c = parameters[PARAM_alpha_d_c * padded_num_cells + i];
            const cellmodel_float_t alpha_n_s = parameters[PARAM_alpha_n_s * padded_num_cells + i];
            const cellmodel_float_t alpha_sl_c =
                    parameters[PARAM_alpha_sl_c * padded_num_cells + i];
            const cellmodel_float_t lambda_c_d =
                    parameters[PARAM_lambda_c_d * padded_num_cells + i];
            const cellmodel_float_t lambda_c_i =
                    parameters[PARAM_lambda_c_i * padded_num_cells + i];
            const cellmodel_float_t lambda_diff =
                    parameters[PARAM_lambda_diff * padded_num_cells + i];
            const cellmodel_float_t Cli = parameters[PARAM_Cli * padded_num_cells + i];
            const cellmodel_float_t Clo = parameters[PARAM_Clo * padded_num_cells + i];
            const cellmodel_float_t Cm = parameters[PARAM_Cm * padded_num_cells + i];
            const cellmodel_float_t Frdy = parameters[PARAM_Frdy * padded_num_cells + i];
            const cellmodel_float_t R = parameters[PARAM_R * padded_num_cells + i];
            const cellmodel_float_t Temp = parameters[PARAM_Temp * padded_num_cells + i];
            const cellmodel_float_t chi = parameters[PARAM_chi * padded_num_cells + i];
            const cellmodel_float_t lambda_c_e =
                    parameters[PARAM_lambda_c_e * padded_num_cells + i];
            const cellmodel_float_t stim_amplitude =
                    parameters[PARAM_stim_amplitude * padded_num_cells + i];
            const cellmodel_float_t stim_duration =
                    parameters[PARAM_stim_duration * padded_num_cells + i];
            const cellmodel_float_t stim_period =
                    parameters[PARAM_stim_period * padded_num_cells + i];
            const cellmodel_float_t stim_start =
                    parameters[PARAM_stim_start * padded_num_cells + i];

            // Assign states
            const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
            const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
            const cellmodel_float_t mL = states[STATE_mL * padded_num_cells + i];
            const cellmodel_float_t hL = states[STATE_hL * padded_num_cells + i];
            const cellmodel_float_t Xr1 = states[STATE_Xr1 * padded_num_cells + i];
            const cellmodel_float_t Xr2 = states[STATE_Xr2 * padded_num_cells + i];
            const cellmodel_float_t x_Ks = states[STATE_x_Ks * padded_num_cells + i];
            const cellmodel_float_t q = states[STATE_q * padded_num_cells + i];
            const cellmodel_float_t r = states[STATE_r * padded_num_cells + i];
            const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
            const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
            const cellmodel_float_t f_Ca_B = states[STATE_f_Ca_B * padded_num_cells + i];
            const cellmodel_float_t xf = states[STATE_xf * padded_num_cells + i];
            const cellmodel_float_t r_RyR = states[STATE_r_RyR * padded_num_cells + i];
            const cellmodel_float_t cn = states[STATE_cn * padded_num_cells + i];
            const cellmodel_float_t cc = states[STATE_cc * padded_num_cells + i];
            const cellmodel_float_t cd = states[STATE_cd * padded_num_cells + i];
            const cellmodel_float_t csl = states[STATE_csl * padded_num_cells + i];
            const cellmodel_float_t cs = states[STATE_cs * padded_num_cells + i];
            const cellmodel_float_t bc = states[STATE_bc * padded_num_cells + i];
            const cellmodel_float_t bd = states[STATE_bd * padded_num_cells + i];
            const cellmodel_float_t bs = states[STATE_bs * padded_num_cells + i];
            const cellmodel_float_t bsl = states[STATE_bsl * padded_num_cells + i];
            const cellmodel_float_t V_m = states[STATE_V_m * padded_num_cells + i];
            const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];

            // Expressions for the Reversal potentials component
            const cellmodel_float_t FoRT = Frdy / (R * Temp);
            const cellmodel_float_t ena = Log(Nao / Na_i) / FoRT;
            const cellmodel_float_t ek = Log(Ko / K_i) / FoRT;
            const cellmodel_float_t eca_sl = Log(ce / csl) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t ecl = Log(Cli / Clo) / FoRT;
            const cellmodel_float_t Qpow = FP_LITERAL(-31.) + Temp / FP_LITERAL(10.);

            // Expressions for the I_Na component
            const cellmodel_float_t mss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(0.00177610354573438) * Exp(-V_m / FP_LITERAL(9.)))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(0.00177610354573438) * Exp(-V_m / FP_LITERAL(9.))));
            const cellmodel_float_t taum =
                    FP_LITERAL(0.06)
                            * Exp(-((FP_LITERAL(-0.0980392156862745) + V_m / FP_LITERAL(51.))
                                    * (FP_LITERAL(-0.0980392156862745) + V_m / FP_LITERAL(51.))))
                    + FP_LITERAL(0.13)
                              * Exp(-((FP_LITERAL(2.875) + V_m / FP_LITERAL(16.))
                                      * (FP_LITERAL(2.875) + V_m / FP_LITERAL(16.))));
            const cellmodel_float_t aj =
                    (V_m >= FP_LITERAL(-38.)
                             ? FP_LITERAL(0.)
                             : (FP_LITERAL(38.) + V_m)
                                       * (FP_LITERAL(-25000.0) * Exp(FP_LITERAL(0.2) * V_m)
                                          - FP_LITERAL(7.0e-6) * Exp(FP_LITERAL(-0.04) * V_m))
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(19623624323.6513)
                                                    * Exp(FP_LITERAL(0.3) * V_m)));
            const cellmodel_float_t bj =
                    (V_m >= FP_LITERAL(-38.218)
                             ? FP_LITERAL(0.6) * Exp(FP_LITERAL(0.09) * V_m)
                                       / (FP_LITERAL(1.) + Exp(FP_LITERAL(-40.) - V_m))
                             : FP_LITERAL(0.02) * Exp(FP_LITERAL(-0.01) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.00369786371648293)
                                                    * Exp(FP_LITERAL(-0.14) * V_m)));
            const cellmodel_float_t tauj = FP_LITERAL(1.0) / (aj + bj);
            const cellmodel_float_t jss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.) + FP_LITERAL(29310.8866998062) * Exp(V_m / FP_LITERAL(7.)))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(29310.8866998062) * Exp(V_m / FP_LITERAL(7.))));
            const cellmodel_float_t I_Na = g_Na * lambda_Na * (m * m * m) * (-ena + V_m) * j;
            const cellmodel_float_t dm_dt = (-m + mss) / taum;
            states[STATE_m * padded_num_cells + i] = dt * dm_dt + m;
            const cellmodel_float_t dj_dt = (-j + jss) / tauj;
            states[STATE_j * padded_num_cells + i] = dt * dj_dt + j;

            // Expressions for the I_NaL component
            const cellmodel_float_t mLss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.000184105793667579) * Exp(-V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tm =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(9.58097877207697) * Exp(V_m / FP_LITERAL(35.))
                       + FP_LITERAL(2.29642676604141e-5) * Exp(-V_m / FP_LITERAL(6.)));
            const cellmodel_float_t tmL = tm;
            const cellmodel_float_t hLss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(124658.506952) * Exp(FP_LITERAL(2.) / FP_LITERAL(15.) * V_m));
            const cellmodel_float_t GNaL =
                    (epi == FP_LITERAL(1.) ? g_NaL : FP_LITERAL(0.6) * g_NaL);
            const cellmodel_float_t I_NaL = lambda_Na * (-ena + V_m) * GNaL * hL * mL;
            const cellmodel_float_t dmL_dt = (-mL + mLss) / tmL;
            states[STATE_mL * padded_num_cells + i] = dt * dmL_dt + mL;
            const cellmodel_float_t dhL_dt = (-hL + hLss) / thL;
            states[STATE_hL * padded_num_cells + i] = dt * dhL_dt + hL;

            // Expressions for the I_NaK component
            const cellmodel_float_t sigma =
                    FP_LITERAL(-0.142857142857143) + Exp(Nao / FP_LITERAL(67.)) / FP_LITERAL(7.);
            const cellmodel_float_t fNaK =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.12) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                       + FP_LITERAL(0.037) * Exp(-FoRT * V_m) * sigma);
            const cellmodel_float_t I_NaK = Ko * g_NaK * fNaK
                                            / ((1.
                                                + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                                          / (((Na_i) * (Na_i)) * ((Na_i) * (Na_i))))
                                               * (KmKo + Ko));

            // Expressions for the I_Kr component
            const cellmodel_float_t Xr1_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(0.0146327985189294)
                                 * Exp(FP_LITERAL(-0.204081632653061) * V_m));
            const cellmodel_float_t alpha_Xr1 =
                    FP_LITERAL(450.0)
                    / (FP_LITERAL(1.0) + FP_LITERAL(0.0111089965382) * Exp(FP_LITERAL(-0.1) * V_m));
            const cellmodel_float_t beta_Xr1 =
                    FP_LITERAL(6.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(13.5813245226) * Exp(FP_LITERAL(0.0869565217391) * V_m));
            const cellmodel_float_t tau_Xr1 = FP_LITERAL(1.0) * alpha_Xr1 * beta_Xr1;
            const cellmodel_float_t Xr2_infinity =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + Exp(FP_LITERAL(44.) / FP_LITERAL(25.) + V_m / FP_LITERAL(50.)));
            const cellmodel_float_t alpha_Xr2 =
                    FP_LITERAL(3.0)
                    / (FP_LITERAL(1.0) + Exp(FP_LITERAL(-3.) - V_m / FP_LITERAL(20.)));
            const cellmodel_float_t beta_Xr2 =
                    FP_LITERAL(1.12)
                    / (FP_LITERAL(1.0) + Exp(FP_LITERAL(-3.) + V_m / FP_LITERAL(20.)));
            const cellmodel_float_t tau_Xr2 = FP_LITERAL(1.0) * alpha_Xr2 * beta_Xr2;
            const cellmodel_float_t I_Kr =
                    FP_LITERAL(0.430331482912) * g_Kr * sqrt(Ko) * (-ek + V_m) * Xr1 * Xr2;
            const cellmodel_float_t dXr1_dt = (-Xr1 + Xr1_inf) / tau_Xr1;
            states[STATE_Xr1 * padded_num_cells + i] = dt * dXr1_dt + Xr1;
            const cellmodel_float_t dXr2_dt = (-Xr2 + Xr2_infinity) / tau_Xr2;
            states[STATE_Xr2 * padded_num_cells + i] = dt * dXr2_dt + Xr2;

            // Expressions for the I_Ks component
            const cellmodel_float_t eks = Log((Ko + Nao * pNaK) / (K_i + pNaK * Na_i)) / FoRT;
            const cellmodel_float_t xsss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.76228973079) * Exp(-V_m / FP_LITERAL(14.)));
            const cellmodel_float_t tauxs =
                    FP_LITERAL(990.)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.842460441617) * Exp(-V_m / FP_LITERAL(14.)));
            const cellmodel_float_t I_Ks = g_Ks * (x_Ks * x_Ks) * (-eks + V_m);
            const cellmodel_float_t dx_Ks_dt = (-x_Ks + xsss) / tauxs;
            states[STATE_x_Ks * padded_num_cells + i] = dt * dx_Ks_dt + x_Ks;

            // Expressions for the i_to component
            const cellmodel_float_t q_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(58.9637634804) * Exp(FP_LITERAL(0.0769230769231) * V_m));
            const cellmodel_float_t tau_q =
                    FP_LITERAL(6.)
                    + FP_LITERAL(39.)
                              / (FP_LITERAL(0.0168716780457) * Exp(FP_LITERAL(-0.08) * V_m)
                                 + FP_LITERAL(6.46648051673) * Exp(FP_LITERAL(0.1) * V_m));
            const cellmodel_float_t r_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(3.28489055021) * Exp(FP_LITERAL(-0.0533333333333) * V_m));
            const cellmodel_float_t tau_r =
                    FP_LITERAL(2.75)
                    + FP_LITERAL(14.4)
                              / (FP_LITERAL(0.0207698622486) * Exp(FP_LITERAL(-0.12) * V_m)
                                 + FP_LITERAL(15.7194688773) * Exp(FP_LITERAL(0.09) * V_m));
            const cellmodel_float_t I_to = g_to * (-ek + V_m) * q * r;
            const cellmodel_float_t dq_dt = (-q + q_inf) / tau_q;
            states[STATE_q * padded_num_cells + i] = dt * dq_dt + q;
            const cellmodel_float_t dr_dt = (-r + r_inf) / tau_r;
            states[STATE_r * padded_num_cells + i] = dt * dr_dt + r;

            // Expressions for the I_K1 component
            const cellmodel_float_t aK1 =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(7.50455791508e-6)
                                 * Exp(FP_LITERAL(0.2) * V_m - FP_LITERAL(0.2) * ek));
            const cellmodel_float_t bK1 =
                    (FP_LITERAL(0.745912348821)
                             * Exp(FP_LITERAL(0.08) * V_m - FP_LITERAL(0.08) * ek)
                     + FP_LITERAL(3.32464030033e-16)
                               * Exp(FP_LITERAL(0.06) * V_m - FP_LITERAL(0.06) * ek))
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.0820849986239)
                                 * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m));
            const cellmodel_float_t K1ss = aK1 / (aK1 + bK1);
            const cellmodel_float_t I_K1 =
                    FP_LITERAL(0.430331482912) * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m) * K1ss;

            // Expressions for the I_bCl component
            const cellmodel_float_t I_bCl = g_bCl * (-ecl + V_m);

            // Expressions for the I_Ca component
            const cellmodel_float_t fss =
                    FP_LITERAL(1.0)
                            / (FP_LITERAL(1.)
                               + FP_LITERAL(48.8565712749872) * Exp(V_m / FP_LITERAL(9.)))
                    + FP_LITERAL(0.6)
                              / (FP_LITERAL(1.)
                                 + FP_LITERAL(12.1824939607035) * Exp(-V_m / FP_LITERAL(20.)));
            const cellmodel_float_t dss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)));
            const cellmodel_float_t taud =
                    (fabs(FP_LITERAL(5.) + V_m) < FP_LITERAL(0.02)
                             ? FP_LITERAL(2.38095238095238)
                             : (FP_LITERAL(1.)
                                - FP_LITERAL(0.434598208507078) * Exp(-V_m / FP_LITERAL(6.)))
                                       * dss / (FP_LITERAL(0.175) + FP_LITERAL(0.035) * V_m));
            const cellmodel_float_t tauf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(0.02)
                       + FP_LITERAL(0.02)
                                 * Exp(-((FP_LITERAL(0.493) + FP_LITERAL(0.034) * V_m)
                                         * (FP_LITERAL(0.493) + FP_LITERAL(0.034) * V_m))));
            const cellmodel_float_t ibarca_j =
                    FP_LITERAL(4.) * Frdy * g_CaL
                    * (FP_LITERAL(-0.34) * ce
                       + FP_LITERAL(0.34) * cd * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t I_CaL =
                    lambda_c_d * pow(Q10CaL, Qpow) * (FP_LITERAL(1.) - f_Ca_B) * d * f * ibarca_j;
            const cellmodel_float_t dd_dt = (-d + dss) / taud;
            states[STATE_d * padded_num_cells + i] = dt * dd_dt + d;
            const cellmodel_float_t df_dt = (-f + fss) / tauf;
            states[STATE_f * padded_num_cells + i] = dt * df_dt + f;
            const cellmodel_float_t df_Ca_B_dt =
                    FP_LITERAL(-0.012) * f_Ca_B + (FP_LITERAL(1.7) - FP_LITERAL(1.7) * f_Ca_B) * cd;
            states[STATE_f_Ca_B * padded_num_cells + i] = dt * df_Ca_B_dt + f_Ca_B;

            // Expressions for the I_NCX component
            const cellmodel_float_t Ka_sl =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (csl * csl));
            const cellmodel_float_t s1_sl = ce * (Na_i * Na_i * Na_i) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s2_sl =
                    (Nao * Nao * Nao) * csl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_sl =
                    KmCao * (Na_i * Na_i * Na_i) + ce * (Na_i * Na_i * Na_i)
                    + (Nao * Nao * Nao) * csl
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_i * Na_i * Na_i) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + csl / KmCai) * csl;
            const cellmodel_float_t I_NaCa =
                    g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);

            // Expressions for the I_PCa component
            const cellmodel_float_t I_pCa = g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl)
                                            / ((KmPCa * KmPCa) + (csl * csl));

            // Expressions for the I_CaBK component
            const cellmodel_float_t I_bCa = g_bCa * lambda_c_e * (-eca_sl + V_m);

            // Expressions for the I_f component
            const cellmodel_float_t xf_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(5956538.01318461) * Exp(V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tau_xf =
                    FP_LITERAL(1900.)
                    / (FP_LITERAL(1.) + FP_LITERAL(4.48168907033806) * Exp(V_m / FP_LITERAL(10.)));
            const cellmodel_float_t I_f = g_f * (-E_f + V_m) * xf;
            const cellmodel_float_t dxf_dt = (-xf + xf_inf) / tau_xf;
            states[STATE_xf * padded_num_cells + i] = dt * dxf_dt + xf;

            // Expressions for the I_KATP component
            const cellmodel_float_t I_KATP = FP_LITERAL(0.60295079490657) * g_KATP
                                             * pow(Ko, FP_LITERAL(0.3)) * (-ek + V_m)
                                             / (FP_LITERAL(40.) + FP_LITERAL(0.0875) * V_m);

            // Expressions for the Ca Fluxes component
            const cellmodel_float_t J_CaL = -Cm * chi * I_CaL / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_pCa = -Cm * chi * I_pCa / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_bCa = -Cm * chi * I_bCa / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_NaCa = Cm * chi * I_NaCa / Frdy;
            const cellmodel_float_t J_e_sl = J_NaCa + J_bCa + J_pCa;
            const cellmodel_float_t Q10SERCA = FP_LITERAL(2.60000000000000);
            const cellmodel_float_t J_SERCA =
                    J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                    * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n))
                    / (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n));
            const cellmodel_float_t J_n_s = alpha_n_s * lambda_c_i * lambda_diff * (-cs + cn);
            const cellmodel_float_t J_sl_c = alpha_sl_c * lambda_c_i * lambda_diff * (-cc + csl);
            const cellmodel_float_t J_d_c = alpha_d_c * lambda_c_d * lambda_diff * (-cc + cd);

            // Expressions for the RyRs component
            const cellmodel_float_t p =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (K_RyR * K_RyR * K_RyR) / (cd * cd * cd));
            const cellmodel_float_t J_RyR_active =
                    alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p * r_RyR;
            const cellmodel_float_t J_leak =
                    alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i * (-csl + cs);
            const cellmodel_float_t J_RyR = J_RyR_active + J_leak;
            const cellmodel_float_t dr_RyR_dt =
                    eta_RyR * (FP_LITERAL(1.) - r_RyR) / p
                    - J_RyR_active / (beta_RyR * lambda_RyR * lambda_c_i);
            states[STATE_r_RyR * padded_num_cells + i] = dt * dr_RyR_dt + r_RyR;

            // Expressions for the Ca Buffers component
            const cellmodel_float_t J_c_b =
                    Vc * (-k_off_c * bc + k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c) * cc);
            const cellmodel_float_t J_d_b =
                    Vd * (-k_off_d * bd + k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c) * cd);
            const cellmodel_float_t J_s_b =
                    Vs * (-k_off_s * bs + k_on_s * (-bs + B_tot_s * lambda_B) * cs);
            const cellmodel_float_t J_sl_b =
                    Vsl
                    * (-k_off_sl * bsl + k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c) * csl);

            // Expressions for the Ca Concentrations component
            const cellmodel_float_t dcn_dt =
                    (FP_LITERAL(1.0) * J_SERCA - FP_LITERAL(1.0) * J_n_s) / Vn;
            states[STATE_cn * padded_num_cells + i] = dt * dcn_dt + cn;
            const cellmodel_float_t dcc_dt = (FP_LITERAL(1.0) * J_d_c + FP_LITERAL(1.0) * J_sl_c
                                              - FP_LITERAL(1.0) * J_SERCA - FP_LITERAL(1.0) * J_c_b)
                                             / Vc;
            states[STATE_cc * padded_num_cells + i] = dt * dcc_dt + cc;
            const cellmodel_float_t dcd_dt =
                    (FP_LITERAL(1.0) * J_CaL - FP_LITERAL(1.0) * J_d_b - FP_LITERAL(1.0) * J_d_c)
                    / Vd;
            states[STATE_cd * padded_num_cells + i] = dt * dcd_dt + cd;
            const cellmodel_float_t dcsl_dt =
                    (FP_LITERAL(1.0) * J_RyR + FP_LITERAL(1.0) * J_e_sl - FP_LITERAL(1.0) * J_sl_b
                     - FP_LITERAL(1.0) * J_sl_c)
                    / Vsl;
            states[STATE_csl * padded_num_cells + i] = dt * dcsl_dt + csl;
            const cellmodel_float_t dcs_dt =
                    (FP_LITERAL(1.0) * J_n_s - FP_LITERAL(1.0) * J_RyR - FP_LITERAL(1.0) * J_s_b)
                    / Vs;
            states[STATE_cs * padded_num_cells + i] = dt * dcs_dt + cs;

            // Expressions for the Ca Buffer Concentrations component
            const cellmodel_float_t dbc_dt = FP_LITERAL(1.0) * J_c_b / Vc;
            states[STATE_bc * padded_num_cells + i] = dt * dbc_dt + bc;
            const cellmodel_float_t dbd_dt = FP_LITERAL(1.0) * J_d_b / Vd;
            states[STATE_bd * padded_num_cells + i] = dt * dbd_dt + bd;
            const cellmodel_float_t dbs_dt = FP_LITERAL(1.0) * J_s_b / Vs;
            states[STATE_bs * padded_num_cells + i] = dt * dbs_dt + bs;
            const cellmodel_float_t dbsl_dt = FP_LITERAL(1.0) * J_sl_b / Vsl;
            states[STATE_bsl * padded_num_cells + i] = dt * dbsl_dt + bsl;

            // Expressions for the Membrane potential component
            const cellmodel_float_t i_Stim =
                    (V_m < FP_LITERAL(-40.) ? FP_LITERAL(1.) : FP_LITERAL(0.))
                    * (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                       && t - stim_period * floor(t / stim_period) >= stim_start
                               ? -stim_amplitude
                               : FP_LITERAL(0.));
            const cellmodel_float_t I_tot = I_CaL + I_K1 + I_KATP + I_Kr + I_Ks + I_Na + I_NaCa
                                            + I_NaK + I_NaL + I_bCa + I_bCl + I_f + I_pCa + I_to;
            const cellmodel_float_t dV_m_dt = -I_tot - i_Stim;
            states[STATE_V_m * padded_num_cells + i] = dt * dV_m_dt + V_m;

            // Expressions for the Sodium concentration component
            const cellmodel_float_t I_Na_tot = FP_LITERAL(3.) * I_NaCa + FP_LITERAL(3.) * I_NaK
                                               + FP_LITERAL(0.3293) * I_f + I_Na + I_NaL;
            const cellmodel_float_t J_Na = -Cm * chi * I_Na_tot / Frdy;
            const cellmodel_float_t dNa_i_dt = J_Na;
            states[STATE_Na_i * padded_num_cells + i] = dt * dNa_i_dt + Na_i;
        }
        t += dt;
    }
}

// Compute multiple forward steps using the GRL1 scheme to the JT21 ODE
static void multistep_GRL1(double *__restrict states, const double t_start, const double dt,
                           const double *__restrict parameters, const long num_cells,
                           const long padded_num_cells, const int num_timesteps)
{
    double t = t_start;
    for (int ti = 0; ti < num_timesteps; ti++) {
#if defined(HINT_CLANG_SIMD)
#pragma omp parallel for
#pragma clang loop vectorize(assume_safety)
#elif defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp parallel for simd aligned(states                                                       \
                                      : CELLMODEL_STATES_ALIGNMENT_BYTES) simdlen(VECTOR_LENGTH)
#else
#pragma omp parallel for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES)
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp parallel for
#endif
        for (long i = 0; i < num_cells; i++) {
            // Assign parameters
            const cellmodel_float_t g_Na = parameters[PARAM_g_Na * padded_num_cells + i];
            const cellmodel_float_t lambda_Na = parameters[PARAM_lambda_Na * padded_num_cells + i];
            const cellmodel_float_t g_NaL = parameters[PARAM_g_NaL * padded_num_cells + i];
            const cellmodel_float_t thL = parameters[PARAM_thL * padded_num_cells + i];
            const cellmodel_float_t KmKo = parameters[PARAM_KmKo * padded_num_cells + i];
            const cellmodel_float_t KmNaip = parameters[PARAM_KmNaip * padded_num_cells + i];
            const cellmodel_float_t g_NaK = parameters[PARAM_g_NaK * padded_num_cells + i];
            const cellmodel_float_t epi = parameters[PARAM_epi * padded_num_cells + i];
            const cellmodel_float_t g_Ks = parameters[PARAM_g_Ks * padded_num_cells + i];
            const cellmodel_float_t pNaK = parameters[PARAM_pNaK * padded_num_cells + i];
            const cellmodel_float_t g_KATP = parameters[PARAM_g_KATP * padded_num_cells + i];
            const cellmodel_float_t g_Kr = parameters[PARAM_g_Kr * padded_num_cells + i];
            const cellmodel_float_t lambda_K = parameters[PARAM_lambda_K * padded_num_cells + i];
            const cellmodel_float_t g_to = parameters[PARAM_g_to * padded_num_cells + i];
            const cellmodel_float_t g_K1 = parameters[PARAM_g_K1 * padded_num_cells + i];
            const cellmodel_float_t g_bCl = parameters[PARAM_g_bCl * padded_num_cells + i];
            const cellmodel_float_t Q10CaL = parameters[PARAM_Q10CaL * padded_num_cells + i];
            const cellmodel_float_t g_CaL = parameters[PARAM_g_CaL * padded_num_cells + i];
            const cellmodel_float_t Kdact = parameters[PARAM_Kdact * padded_num_cells + i];
            const cellmodel_float_t KmCai = parameters[PARAM_KmCai * padded_num_cells + i];
            const cellmodel_float_t KmCao = parameters[PARAM_KmCao * padded_num_cells + i];
            const cellmodel_float_t KmNai = parameters[PARAM_KmNai * padded_num_cells + i];
            const cellmodel_float_t KmNao = parameters[PARAM_KmNao * padded_num_cells + i];
            const cellmodel_float_t Q10NCX = parameters[PARAM_Q10NCX * padded_num_cells + i];
            const cellmodel_float_t g_NaCa = parameters[PARAM_g_NaCa * padded_num_cells + i];
            const cellmodel_float_t ksat = parameters[PARAM_ksat * padded_num_cells + i];
            const cellmodel_float_t nu = parameters[PARAM_nu * padded_num_cells + i];
            const cellmodel_float_t KmPCa = parameters[PARAM_KmPCa * padded_num_cells + i];
            const cellmodel_float_t Q10SLCaP = parameters[PARAM_Q10SLCaP * padded_num_cells + i];
            const cellmodel_float_t g_pCa = parameters[PARAM_g_pCa * padded_num_cells + i];
            const cellmodel_float_t g_bCa = parameters[PARAM_g_bCa * padded_num_cells + i];
            const cellmodel_float_t E_f = parameters[PARAM_E_f * padded_num_cells + i];
            const cellmodel_float_t g_f = parameters[PARAM_g_f * padded_num_cells + i];
            const cellmodel_float_t Nao = parameters[PARAM_Nao * padded_num_cells + i];
            const cellmodel_float_t K_i = parameters[PARAM_K_i * padded_num_cells + i];
            const cellmodel_float_t Ko = parameters[PARAM_Ko * padded_num_cells + i];
            const cellmodel_float_t ce = parameters[PARAM_ce * padded_num_cells + i];
            const cellmodel_float_t K_RyR = parameters[PARAM_K_RyR * padded_num_cells + i];
            const cellmodel_float_t alpha_RyR = parameters[PARAM_alpha_RyR * padded_num_cells + i];
            const cellmodel_float_t beta_RyR = parameters[PARAM_beta_RyR * padded_num_cells + i];
            const cellmodel_float_t eta_RyR = parameters[PARAM_eta_RyR * padded_num_cells + i];
            const cellmodel_float_t gamma_RyR = parameters[PARAM_gamma_RyR * padded_num_cells + i];
            const cellmodel_float_t lambda_RyR =
                    parameters[PARAM_lambda_RyR * padded_num_cells + i];
            const cellmodel_float_t Vc = parameters[PARAM_Vc * padded_num_cells + i];
            const cellmodel_float_t Vd = parameters[PARAM_Vd * padded_num_cells + i];
            const cellmodel_float_t Vn = parameters[PARAM_Vn * padded_num_cells + i];
            const cellmodel_float_t Vs = parameters[PARAM_Vs * padded_num_cells + i];
            const cellmodel_float_t Vsl = parameters[PARAM_Vsl * padded_num_cells + i];
            const cellmodel_float_t J_SERCA_bar =
                    parameters[PARAM_J_SERCA_bar * padded_num_cells + i];
            const cellmodel_float_t K_c = parameters[PARAM_K_c * padded_num_cells + i];
            const cellmodel_float_t K_n = parameters[PARAM_K_n * padded_num_cells + i];
            const cellmodel_float_t B_tot_c = parameters[PARAM_B_tot_c * padded_num_cells + i];
            const cellmodel_float_t B_tot_d = parameters[PARAM_B_tot_d * padded_num_cells + i];
            const cellmodel_float_t B_tot_s = parameters[PARAM_B_tot_s * padded_num_cells + i];
            const cellmodel_float_t B_tot_sl = parameters[PARAM_B_tot_sl * padded_num_cells + i];
            const cellmodel_float_t k_off_c = parameters[PARAM_k_off_c * padded_num_cells + i];
            const cellmodel_float_t k_off_d = parameters[PARAM_k_off_d * padded_num_cells + i];
            const cellmodel_float_t k_off_s = parameters[PARAM_k_off_s * padded_num_cells + i];
            const cellmodel_float_t k_off_sl = parameters[PARAM_k_off_sl * padded_num_cells + i];
            const cellmodel_float_t k_on_c = parameters[PARAM_k_on_c * padded_num_cells + i];
            const cellmodel_float_t k_on_d = parameters[PARAM_k_on_d * padded_num_cells + i];
            const cellmodel_float_t k_on_s = parameters[PARAM_k_on_s * padded_num_cells + i];
            const cellmodel_float_t k_on_sl = parameters[PARAM_k_on_sl * padded_num_cells + i];
            const cellmodel_float_t lambda_B = parameters[PARAM_lambda_B * padded_num_cells + i];
            const cellmodel_float_t lambda_B_c =
                    parameters[PARAM_lambda_B_c * padded_num_cells + i];
            const cellmodel_float_t alpha_d_c = parameters[PARAM_alpha_d_c * padded_num_cells + i];
            const cellmodel_float_t alpha_n_s = parameters[PARAM_alpha_n_s * padded_num_cells + i];
            const cellmodel_float_t alpha_sl_c =
                    parameters[PARAM_alpha_sl_c * padded_num_cells + i];
            const cellmodel_float_t lambda_c_d =
                    parameters[PARAM_lambda_c_d * padded_num_cells + i];
            const cellmodel_float_t lambda_c_i =
                    parameters[PARAM_lambda_c_i * padded_num_cells + i];
            const cellmodel_float_t lambda_diff =
                    parameters[PARAM_lambda_diff * padded_num_cells + i];
            const cellmodel_float_t Cli = parameters[PARAM_Cli * padded_num_cells + i];
            const cellmodel_float_t Clo = parameters[PARAM_Clo * padded_num_cells + i];
            const cellmodel_float_t Cm = parameters[PARAM_Cm * padded_num_cells + i];
            const cellmodel_float_t Frdy = parameters[PARAM_Frdy * padded_num_cells + i];
            const cellmodel_float_t R = parameters[PARAM_R * padded_num_cells + i];
            const cellmodel_float_t Temp = parameters[PARAM_Temp * padded_num_cells + i];
            const cellmodel_float_t chi = parameters[PARAM_chi * padded_num_cells + i];
            const cellmodel_float_t lambda_c_e =
                    parameters[PARAM_lambda_c_e * padded_num_cells + i];
            const cellmodel_float_t stim_amplitude =
                    parameters[PARAM_stim_amplitude * padded_num_cells + i];
            const cellmodel_float_t stim_duration =
                    parameters[PARAM_stim_duration * padded_num_cells + i];
            const cellmodel_float_t stim_period =
                    parameters[PARAM_stim_period * padded_num_cells + i];
            const cellmodel_float_t stim_start =
                    parameters[PARAM_stim_start * padded_num_cells + i];

            // Assign states
            const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
            const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
            const cellmodel_float_t mL = states[STATE_mL * padded_num_cells + i];
            const cellmodel_float_t hL = states[STATE_hL * padded_num_cells + i];
            const cellmodel_float_t Xr1 = states[STATE_Xr1 * padded_num_cells + i];
            const cellmodel_float_t Xr2 = states[STATE_Xr2 * padded_num_cells + i];
            const cellmodel_float_t x_Ks = states[STATE_x_Ks * padded_num_cells + i];
            const cellmodel_float_t q = states[STATE_q * padded_num_cells + i];
            const cellmodel_float_t r = states[STATE_r * padded_num_cells + i];
            const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
            const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
            const cellmodel_float_t f_Ca_B = states[STATE_f_Ca_B * padded_num_cells + i];
            const cellmodel_float_t xf = states[STATE_xf * padded_num_cells + i];
            const cellmodel_float_t r_RyR = states[STATE_r_RyR * padded_num_cells + i];
            const cellmodel_float_t cn = states[STATE_cn * padded_num_cells + i];
            const cellmodel_float_t cc = states[STATE_cc * padded_num_cells + i];
            const cellmodel_float_t cd = states[STATE_cd * padded_num_cells + i];
            const cellmodel_float_t csl = states[STATE_csl * padded_num_cells + i];
            const cellmodel_float_t cs = states[STATE_cs * padded_num_cells + i];
            const cellmodel_float_t bc = states[STATE_bc * padded_num_cells + i];
            const cellmodel_float_t bd = states[STATE_bd * padded_num_cells + i];
            const cellmodel_float_t bs = states[STATE_bs * padded_num_cells + i];
            const cellmodel_float_t bsl = states[STATE_bsl * padded_num_cells + i];
            const cellmodel_float_t V_m = states[STATE_V_m * padded_num_cells + i];
            const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];

            // Expressions for the Reversal potentials component
            const cellmodel_float_t FoRT = Frdy / (R * Temp);
            const cellmodel_float_t ena = Log(Nao / Na_i) / FoRT;
            const cellmodel_float_t ek = Log(Ko / K_i) / FoRT;
            const cellmodel_float_t eca_sl = Log(ce / csl) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t ecl = Log(Cli / Clo) / FoRT;
            const cellmodel_float_t Qpow = FP_LITERAL(-31.) + Temp / FP_LITERAL(10.);

            // Expressions for the I_Na component
            const cellmodel_float_t mss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(0.00177610354573438) * Exp(-V_m / FP_LITERAL(9.)))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(0.00177610354573438) * Exp(-V_m / FP_LITERAL(9.))));
            const cellmodel_float_t taum =
                    FP_LITERAL(0.06)
                            * Exp(-((FP_LITERAL(-0.0980392156862745) + V_m / FP_LITERAL(51.))
                                    * (FP_LITERAL(-0.0980392156862745) + V_m / FP_LITERAL(51.))))
                    + FP_LITERAL(0.13)
                              * Exp(-((FP_LITERAL(2.875) + V_m / FP_LITERAL(16.))
                                      * (FP_LITERAL(2.875) + V_m / FP_LITERAL(16.))));
            const cellmodel_float_t aj =
                    (V_m >= FP_LITERAL(-38.)
                             ? FP_LITERAL(0.)
                             : (FP_LITERAL(38.) + V_m)
                                       * (FP_LITERAL(-25000.0) * Exp(FP_LITERAL(0.2) * V_m)
                                          - FP_LITERAL(7.0e-6) * Exp(FP_LITERAL(-0.04) * V_m))
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(19623624323.6513)
                                                    * Exp(FP_LITERAL(0.3) * V_m)));
            const cellmodel_float_t bj =
                    (V_m >= FP_LITERAL(-38.218)
                             ? FP_LITERAL(0.6) * Exp(FP_LITERAL(0.09) * V_m)
                                       / (FP_LITERAL(1.) + Exp(FP_LITERAL(-40.) - V_m))
                             : FP_LITERAL(0.02) * Exp(FP_LITERAL(-0.01) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.00369786371648293)
                                                    * Exp(FP_LITERAL(-0.14) * V_m)));
            const cellmodel_float_t tauj = FP_LITERAL(1.0) / (aj + bj);
            const cellmodel_float_t jss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.) + FP_LITERAL(29310.8866998062) * Exp(V_m / FP_LITERAL(7.)))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(29310.8866998062) * Exp(V_m / FP_LITERAL(7.))));
            const cellmodel_float_t I_Na = g_Na * lambda_Na * (m * m * m) * (-ena + V_m) * j;
            const cellmodel_float_t dm_dt = (-m + mss) / taum;
            const cellmodel_float_t dm_dt_linearized = FP_LITERAL(-1.) / taum;
            states[STATE_m * padded_num_cells + i] =
                    (fabs(dm_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dm_dt_linearized)) * dm_dt
                                       / dm_dt_linearized
                             : dt * dm_dt)
                    + m;
            const cellmodel_float_t dj_dt = (-j + jss) / tauj;
            const cellmodel_float_t dj_dt_linearized = FP_LITERAL(-1.) / tauj;
            states[STATE_j * padded_num_cells + i] =
                    (fabs(dj_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dj_dt_linearized)) * dj_dt
                                       / dj_dt_linearized
                             : dt * dj_dt)
                    + j;

            // Expressions for the I_NaL component
            const cellmodel_float_t mLss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.000184105793667579) * Exp(-V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tm =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(9.58097877207697) * Exp(V_m / FP_LITERAL(35.))
                       + FP_LITERAL(2.29642676604141e-5) * Exp(-V_m / FP_LITERAL(6.)));
            const cellmodel_float_t tmL = tm;
            const cellmodel_float_t hLss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(124658.506952) * Exp(FP_LITERAL(2.) / FP_LITERAL(15.) * V_m));
            const cellmodel_float_t GNaL =
                    (epi == FP_LITERAL(1.) ? g_NaL : FP_LITERAL(0.6) * g_NaL);
            const cellmodel_float_t I_NaL = lambda_Na * (-ena + V_m) * GNaL * hL * mL;
            const cellmodel_float_t dmL_dt = (-mL + mLss) / tmL;
            const cellmodel_float_t dmL_dt_linearized = FP_LITERAL(-1.) / tmL;
            states[STATE_mL * padded_num_cells + i] =
                    (fabs(dmL_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dmL_dt_linearized)) * dmL_dt
                                       / dmL_dt_linearized
                             : dt * dmL_dt)
                    + mL;
            const cellmodel_float_t dhL_dt = (-hL + hLss) / thL;
            const cellmodel_float_t dhL_dt_linearized = FP_LITERAL(-1.) / thL;
            states[STATE_hL * padded_num_cells + i] =
                    (fabs(dhL_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dhL_dt_linearized)) * dhL_dt
                                       / dhL_dt_linearized
                             : dt * dhL_dt)
                    + hL;

            // Expressions for the I_NaK component
            const cellmodel_float_t sigma =
                    FP_LITERAL(-0.142857142857143) + Exp(Nao / FP_LITERAL(67.)) / FP_LITERAL(7.);
            const cellmodel_float_t fNaK =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.12) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                       + FP_LITERAL(0.037) * Exp(-FoRT * V_m) * sigma);
            const cellmodel_float_t I_NaK = Ko * g_NaK * fNaK
                                            / ((1.
                                                + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                                          / (((Na_i) * (Na_i)) * ((Na_i) * (Na_i))))
                                               * (KmKo + Ko));

            // Expressions for the I_Kr component
            const cellmodel_float_t Xr1_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(0.0146327985189294)
                                 * Exp(FP_LITERAL(-0.204081632653061) * V_m));
            const cellmodel_float_t alpha_Xr1 =
                    FP_LITERAL(450.0)
                    / (FP_LITERAL(1.0) + FP_LITERAL(0.0111089965382) * Exp(FP_LITERAL(-0.1) * V_m));
            const cellmodel_float_t beta_Xr1 =
                    FP_LITERAL(6.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(13.5813245226) * Exp(FP_LITERAL(0.0869565217391) * V_m));
            const cellmodel_float_t tau_Xr1 = FP_LITERAL(1.0) * alpha_Xr1 * beta_Xr1;
            const cellmodel_float_t Xr2_infinity =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + Exp(FP_LITERAL(44.) / FP_LITERAL(25.) + V_m / FP_LITERAL(50.)));
            const cellmodel_float_t alpha_Xr2 =
                    FP_LITERAL(3.0)
                    / (FP_LITERAL(1.0) + Exp(FP_LITERAL(-3.) - V_m / FP_LITERAL(20.)));
            const cellmodel_float_t beta_Xr2 =
                    FP_LITERAL(1.12)
                    / (FP_LITERAL(1.0) + Exp(FP_LITERAL(-3.) + V_m / FP_LITERAL(20.)));
            const cellmodel_float_t tau_Xr2 = FP_LITERAL(1.0) * alpha_Xr2 * beta_Xr2;
            const cellmodel_float_t I_Kr =
                    FP_LITERAL(0.430331482912) * g_Kr * sqrt(Ko) * (-ek + V_m) * Xr1 * Xr2;
            const cellmodel_float_t dXr1_dt = (-Xr1 + Xr1_inf) / tau_Xr1;
            const cellmodel_float_t dXr1_dt_linearized = FP_LITERAL(-1.) / tau_Xr1;
            states[STATE_Xr1 * padded_num_cells + i] =
                    (fabs(dXr1_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dXr1_dt_linearized)) * dXr1_dt
                                       / dXr1_dt_linearized
                             : dt * dXr1_dt)
                    + Xr1;
            const cellmodel_float_t dXr2_dt = (-Xr2 + Xr2_infinity) / tau_Xr2;
            const cellmodel_float_t dXr2_dt_linearized = FP_LITERAL(-1.) / tau_Xr2;
            states[STATE_Xr2 * padded_num_cells + i] =
                    (fabs(dXr2_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dXr2_dt_linearized)) * dXr2_dt
                                       / dXr2_dt_linearized
                             : dt * dXr2_dt)
                    + Xr2;

            // Expressions for the I_Ks component
            const cellmodel_float_t eks = Log((Ko + Nao * pNaK) / (K_i + pNaK * Na_i)) / FoRT;
            const cellmodel_float_t xsss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.76228973079) * Exp(-V_m / FP_LITERAL(14.)));
            const cellmodel_float_t tauxs =
                    FP_LITERAL(990.)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.842460441617) * Exp(-V_m / FP_LITERAL(14.)));
            const cellmodel_float_t I_Ks = g_Ks * (x_Ks * x_Ks) * (-eks + V_m);
            const cellmodel_float_t dx_Ks_dt = (-x_Ks + xsss) / tauxs;
            const cellmodel_float_t dx_Ks_dt_linearized = FP_LITERAL(-1.) / tauxs;
            states[STATE_x_Ks * padded_num_cells + i] =
                    (fabs(dx_Ks_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dx_Ks_dt_linearized)) * dx_Ks_dt
                                       / dx_Ks_dt_linearized
                             : dt * dx_Ks_dt)
                    + x_Ks;

            // Expressions for the i_to component
            const cellmodel_float_t q_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(58.9637634804) * Exp(FP_LITERAL(0.0769230769231) * V_m));
            const cellmodel_float_t tau_q =
                    FP_LITERAL(6.)
                    + FP_LITERAL(39.)
                              / (FP_LITERAL(0.0168716780457) * Exp(FP_LITERAL(-0.08) * V_m)
                                 + FP_LITERAL(6.46648051673) * Exp(FP_LITERAL(0.1) * V_m));
            const cellmodel_float_t r_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.0)
                       + FP_LITERAL(3.28489055021) * Exp(FP_LITERAL(-0.0533333333333) * V_m));
            const cellmodel_float_t tau_r =
                    FP_LITERAL(2.75)
                    + FP_LITERAL(14.4)
                              / (FP_LITERAL(0.0207698622486) * Exp(FP_LITERAL(-0.12) * V_m)
                                 + FP_LITERAL(15.7194688773) * Exp(FP_LITERAL(0.09) * V_m));
            const cellmodel_float_t I_to = g_to * (-ek + V_m) * q * r;
            const cellmodel_float_t dq_dt = (-q + q_inf) / tau_q;
            const cellmodel_float_t dq_dt_linearized = FP_LITERAL(-1.) / tau_q;
            states[STATE_q * padded_num_cells + i] =
                    (fabs(dq_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dq_dt_linearized)) * dq_dt
                                       / dq_dt_linearized
                             : dt * dq_dt)
                    + q;
            const cellmodel_float_t dr_dt = (-r + r_inf) / tau_r;
            const cellmodel_float_t dr_dt_linearized = FP_LITERAL(-1.) / tau_r;
            states[STATE_r * padded_num_cells + i] =
                    (fabs(dr_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dr_dt_linearized)) * dr_dt
                                       / dr_dt_linearized
                             : dt * dr_dt)
                    + r;

            // Expressions for the I_K1 component
            const cellmodel_float_t aK1 =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(7.50455791508e-6)
                                 * Exp(FP_LITERAL(0.2) * V_m - FP_LITERAL(0.2) * ek));
            const cellmodel_float_t bK1 =
                    (FP_LITERAL(0.745912348821)
                             * Exp(FP_LITERAL(0.08) * V_m - FP_LITERAL(0.08) * ek)
                     + FP_LITERAL(3.32464030033e-16)
                               * Exp(FP_LITERAL(0.06) * V_m - FP_LITERAL(0.06) * ek))
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.0820849986239)
                                 * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m));
            const cellmodel_float_t K1ss = aK1 / (aK1 + bK1);
            const cellmodel_float_t I_K1 =
                    FP_LITERAL(0.430331482912) * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m) * K1ss;

            // Expressions for the I_bCl component
            const cellmodel_float_t I_bCl = g_bCl * (-ecl + V_m);

            // Expressions for the I_Ca component
            const cellmodel_float_t fss =
                    FP_LITERAL(1.0)
                            / (FP_LITERAL(1.)
                               + FP_LITERAL(48.8565712749872) * Exp(V_m / FP_LITERAL(9.)))
                    + FP_LITERAL(0.6)
                              / (FP_LITERAL(1.)
                                 + FP_LITERAL(12.1824939607035) * Exp(-V_m / FP_LITERAL(20.)));
            const cellmodel_float_t dss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)));
            const cellmodel_float_t taud =
                    (fabs(FP_LITERAL(5.) + V_m) < FP_LITERAL(0.02)
                             ? FP_LITERAL(2.38095238095238)
                             : (FP_LITERAL(1.)
                                - FP_LITERAL(0.434598208507078) * Exp(-V_m / FP_LITERAL(6.)))
                                       * dss / (FP_LITERAL(0.175) + FP_LITERAL(0.035) * V_m));
            const cellmodel_float_t tauf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(0.02)
                       + FP_LITERAL(0.02)
                                 * Exp(-((FP_LITERAL(0.493) + FP_LITERAL(0.034) * V_m)
                                         * (FP_LITERAL(0.493) + FP_LITERAL(0.034) * V_m))));
            const cellmodel_float_t ibarca_j =
                    FP_LITERAL(4.) * Frdy * g_CaL
                    * (FP_LITERAL(-0.34) * ce
                       + FP_LITERAL(0.34) * cd * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t I_CaL =
                    lambda_c_d * pow(Q10CaL, Qpow) * (FP_LITERAL(1.) - f_Ca_B) * d * f * ibarca_j;
            const cellmodel_float_t dd_dt = (-d + dss) / taud;
            const cellmodel_float_t dd_dt_linearized = FP_LITERAL(-1.) / taud;
            states[STATE_d * padded_num_cells + i] =
                    (fabs(dd_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dd_dt_linearized)) * dd_dt
                                       / dd_dt_linearized
                             : dt * dd_dt)
                    + d;
            const cellmodel_float_t df_dt = (-f + fss) / tauf;
            const cellmodel_float_t df_dt_linearized = FP_LITERAL(-1.) / tauf;
            states[STATE_f * padded_num_cells + i] =
                    (fabs(df_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * df_dt_linearized)) * df_dt
                                       / df_dt_linearized
                             : dt * df_dt)
                    + f;
            const cellmodel_float_t df_Ca_B_dt =
                    FP_LITERAL(-0.012) * f_Ca_B + (FP_LITERAL(1.7) - FP_LITERAL(1.7) * f_Ca_B) * cd;
            const cellmodel_float_t df_Ca_B_dt_linearized =
                    FP_LITERAL(-0.012) - FP_LITERAL(1.7) * cd;
            states[STATE_f_Ca_B * padded_num_cells + i] =
                    (fabs(df_Ca_B_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * df_Ca_B_dt_linearized)) * df_Ca_B_dt
                                       / df_Ca_B_dt_linearized
                             : dt * df_Ca_B_dt)
                    + f_Ca_B;

            // Expressions for the I_NCX component
            const cellmodel_float_t Ka_sl =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (csl * csl));
            const cellmodel_float_t s1_sl = ce * (Na_i * Na_i * Na_i) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s2_sl =
                    (Nao * Nao * Nao) * csl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_sl =
                    KmCao * (Na_i * Na_i * Na_i) + ce * (Na_i * Na_i * Na_i)
                    + (Nao * Nao * Nao) * csl
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_i * Na_i * Na_i) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + csl / KmCai) * csl;
            const cellmodel_float_t I_NaCa =
                    g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);

            // Expressions for the I_PCa component
            const cellmodel_float_t I_pCa = g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl)
                                            / ((KmPCa * KmPCa) + (csl * csl));

            // Expressions for the I_CaBK component
            const cellmodel_float_t I_bCa = g_bCa * lambda_c_e * (-eca_sl + V_m);

            // Expressions for the I_f component
            const cellmodel_float_t xf_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(5956538.01318461) * Exp(V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tau_xf =
                    FP_LITERAL(1900.)
                    / (FP_LITERAL(1.) + FP_LITERAL(4.48168907033806) * Exp(V_m / FP_LITERAL(10.)));
            const cellmodel_float_t I_f = g_f * (-E_f + V_m) * xf;
            const cellmodel_float_t dxf_dt = (-xf + xf_inf) / tau_xf;
            const cellmodel_float_t dxf_dt_linearized = FP_LITERAL(-1.) / tau_xf;
            states[STATE_xf * padded_num_cells + i] =
                    (fabs(dxf_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dxf_dt_linearized)) * dxf_dt
                                       / dxf_dt_linearized
                             : dt * dxf_dt)
                    + xf;

            // Expressions for the I_KATP component
            const cellmodel_float_t I_KATP = FP_LITERAL(0.60295079490657) * g_KATP
                                             * pow(Ko, FP_LITERAL(0.3)) * (-ek + V_m)
                                             / (FP_LITERAL(40.) + FP_LITERAL(0.0875) * V_m);

            // Expressions for the Ca Fluxes component
            const cellmodel_float_t J_CaL = -Cm * chi * I_CaL / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_pCa = -Cm * chi * I_pCa / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_bCa = -Cm * chi * I_bCa / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t J_NaCa = Cm * chi * I_NaCa / Frdy;
            const cellmodel_float_t J_e_sl = J_NaCa + J_bCa + J_pCa;
            const cellmodel_float_t Q10SERCA = FP_LITERAL(2.60000000000000);
            const cellmodel_float_t J_SERCA =
                    J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                    * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n))
                    / (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c) + (cn * cn) / (K_n * K_n));
            const cellmodel_float_t J_n_s = alpha_n_s * lambda_c_i * lambda_diff * (-cs + cn);
            const cellmodel_float_t J_sl_c = alpha_sl_c * lambda_c_i * lambda_diff * (-cc + csl);
            const cellmodel_float_t J_d_c = alpha_d_c * lambda_c_d * lambda_diff * (-cc + cd);

            // Expressions for the RyRs component
            const cellmodel_float_t p =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (K_RyR * K_RyR * K_RyR) / (cd * cd * cd));
            const cellmodel_float_t J_RyR_active =
                    alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p * r_RyR;
            const cellmodel_float_t J_leak =
                    alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i * (-csl + cs);
            const cellmodel_float_t J_RyR = J_RyR_active + J_leak;
            const cellmodel_float_t dr_RyR_dt =
                    eta_RyR * (FP_LITERAL(1.) - r_RyR) / p
                    - J_RyR_active / (beta_RyR * lambda_RyR * lambda_c_i);
            const cellmodel_float_t dJ_RyR_active_dr_RyR =
                    alpha_RyR * lambda_RyR * lambda_c_i * (-csl + cs) * p;
            const cellmodel_float_t dr_RyR_dt_linearized =
                    -eta_RyR / p - dJ_RyR_active_dr_RyR / (beta_RyR * lambda_RyR * lambda_c_i);
            states[STATE_r_RyR * padded_num_cells + i] =
                    (fabs(dr_RyR_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dr_RyR_dt_linearized)) * dr_RyR_dt
                                       / dr_RyR_dt_linearized
                             : dt * dr_RyR_dt)
                    + r_RyR;

            // Expressions for the Ca Buffers component
            const cellmodel_float_t J_c_b =
                    Vc * (-k_off_c * bc + k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c) * cc);
            const cellmodel_float_t J_d_b =
                    Vd * (-k_off_d * bd + k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c) * cd);
            const cellmodel_float_t J_s_b =
                    Vs * (-k_off_s * bs + k_on_s * (-bs + B_tot_s * lambda_B) * cs);
            const cellmodel_float_t J_sl_b =
                    Vsl
                    * (-k_off_sl * bsl + k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c) * csl);

            // Expressions for the Ca Concentrations component
            const cellmodel_float_t dcn_dt =
                    (FP_LITERAL(1.0) * J_SERCA - FP_LITERAL(1.0) * J_n_s) / Vn;
            const cellmodel_float_t dJ_SERCA_dcn =
                    FP_LITERAL(-2.) * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow) * cn
                            / ((K_n * K_n)
                               * (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                  + (cn * cn) / (K_n * K_n)))
                    - FP_LITERAL(2.) * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                              * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n)) * cn
                              / ((K_n * K_n)
                                 * ((FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                     + (cn * cn) / (K_n * K_n))
                                    * (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                       + (cn * cn) / (K_n * K_n))));
            const cellmodel_float_t dJ_n_s_dcn = alpha_n_s * lambda_c_i * lambda_diff;
            const cellmodel_float_t dcn_dt_linearized =
                    (FP_LITERAL(1.0) * dJ_SERCA_dcn - FP_LITERAL(1.0) * dJ_n_s_dcn) / Vn;
            states[STATE_cn * padded_num_cells + i] =
                    (fabs(dcn_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dcn_dt_linearized)) * dcn_dt
                                       / dcn_dt_linearized
                             : dt * dcn_dt)
                    + cn;
            const cellmodel_float_t dcc_dt = (FP_LITERAL(1.0) * J_d_c + FP_LITERAL(1.0) * J_sl_c
                                              - FP_LITERAL(1.0) * J_SERCA - FP_LITERAL(1.0) * J_c_b)
                                             / Vc;
            const cellmodel_float_t dJ_SERCA_dcc =
                    FP_LITERAL(2.) * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow) * cc
                            / ((K_c * K_c)
                               * (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                  + (cn * cn) / (K_n * K_n)))
                    - FP_LITERAL(2.) * J_SERCA_bar * lambda_c_i * pow(Q10SERCA, Qpow)
                              * ((cc * cc) / (K_c * K_c) - (cn * cn) / (K_n * K_n)) * cc
                              / ((K_c * K_c)
                                 * ((FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                     + (cn * cn) / (K_n * K_n))
                                    * (FP_LITERAL(1.) + (cc * cc) / (K_c * K_c)
                                       + (cn * cn) / (K_n * K_n))));
            const cellmodel_float_t dJ_c_b_dcc =
                    Vc * k_on_c * (-bc + B_tot_c * lambda_B * lambda_B_c);
            const cellmodel_float_t dJ_d_c_dcc = -alpha_d_c * lambda_c_d * lambda_diff;
            const cellmodel_float_t dJ_sl_c_dcc = -alpha_sl_c * lambda_c_i * lambda_diff;
            const cellmodel_float_t dcc_dt_linearized =
                    (FP_LITERAL(1.0) * dJ_d_c_dcc + FP_LITERAL(1.0) * dJ_sl_c_dcc
                     - FP_LITERAL(1.0) * dJ_SERCA_dcc - FP_LITERAL(1.0) * dJ_c_b_dcc)
                    / Vc;
            states[STATE_cc * padded_num_cells + i] =
                    (fabs(dcc_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dcc_dt_linearized)) * dcc_dt
                                       / dcc_dt_linearized
                             : dt * dcc_dt)
                    + cc;
            const cellmodel_float_t dcd_dt =
                    (FP_LITERAL(1.0) * J_CaL - FP_LITERAL(1.0) * J_d_b - FP_LITERAL(1.0) * J_d_c)
                    / Vd;
            const cellmodel_float_t dI_CaL_dibarca_j =
                    lambda_c_d * pow(Q10CaL, Qpow) * (FP_LITERAL(1.) - f_Ca_B) * d * f;
            const cellmodel_float_t dJ_CaL_dI_CaL = -Cm * chi / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t dJ_d_b_dcd =
                    Vd * k_on_d * (-bd + B_tot_d * lambda_B * lambda_B_c);
            const cellmodel_float_t dJ_d_c_dcd = alpha_d_c * lambda_c_d * lambda_diff;
            const cellmodel_float_t dibarca_j_dcd =
                    FP_LITERAL(1.36) * Frdy * g_CaL * FoRT * V_m * Exp(FP_LITERAL(2.) * FoRT * V_m)
                    / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t dcd_dt_linearized =
                    (FP_LITERAL(-1.0) * dJ_d_b_dcd - FP_LITERAL(1.0) * dJ_d_c_dcd
                     + FP_LITERAL(1.0) * dI_CaL_dibarca_j * dJ_CaL_dI_CaL * dibarca_j_dcd)
                    / Vd;
            states[STATE_cd * padded_num_cells + i] =
                    (fabs(dcd_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dcd_dt_linearized)) * dcd_dt
                                       / dcd_dt_linearized
                             : dt * dcd_dt)
                    + cd;
            const cellmodel_float_t dcsl_dt =
                    (FP_LITERAL(1.0) * J_RyR + FP_LITERAL(1.0) * J_e_sl - FP_LITERAL(1.0) * J_sl_b
                     - FP_LITERAL(1.0) * J_sl_c)
                    / Vsl;
            const cellmodel_float_t dI_NaCa_dKa_sl =
                    g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl)
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t dI_NaCa_ds2_sl =
                    -g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t dI_NaCa_ds3_sl =
                    -g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * (s3_sl * s3_sl));
            const cellmodel_float_t dI_bCa_deca_sl = -g_bCa * lambda_c_e;
            const cellmodel_float_t dI_pCa_dcsl =
                    FP_LITERAL(-2.) * g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * (csl * csl * csl)
                            / (((KmPCa * KmPCa) + (csl * csl)) * ((KmPCa * KmPCa) + (csl * csl)))
                    + FP_LITERAL(2.) * g_pCa * lambda_c_e * pow(Q10SLCaP, Qpow) * csl
                              / ((KmPCa * KmPCa) + (csl * csl));
            const cellmodel_float_t dJ_NaCa_dI_NaCa = Cm * chi / Frdy;
            const cellmodel_float_t dJ_RyR_active_dcsl =
                    -alpha_RyR * lambda_RyR * lambda_c_i * p * r_RyR;
            const cellmodel_float_t dJ_bCa_dI_bCa = -Cm * chi / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t dJ_leak_dcsl = -alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i;
            const cellmodel_float_t dJ_pCa_dI_pCa = -Cm * chi / (FP_LITERAL(2.) * Frdy);
            const cellmodel_float_t dJ_sl_b_dcsl =
                    Vsl * k_on_sl * (-bsl + B_tot_sl * lambda_B * lambda_B_c);
            const cellmodel_float_t dJ_sl_c_dcsl = alpha_sl_c * lambda_c_i * lambda_diff;
            const cellmodel_float_t dKa_sl_dcsl =
                    FP_LITERAL(2.0) * (Kdact * Kdact)
                    / (((FP_LITERAL(1.) + (Kdact * Kdact) / (csl * csl))
                        * (FP_LITERAL(1.) + (Kdact * Kdact) / (csl * csl)))
                       * (csl * csl * csl));
            const cellmodel_float_t deca_sl_dcsl = FP_LITERAL(-1.) / (FP_LITERAL(2.) * FoRT * csl);
            const cellmodel_float_t ds2_sl_dcsl =
                    (Nao * Nao * Nao) * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t ds3_sl_dcsl =
                    (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + csl / KmCai)
                    + (KmNao * KmNao * KmNao) * csl / KmCai;
            const cellmodel_float_t dcsl_dt_linearized =
                    (FP_LITERAL(1.0) * dJ_RyR_active_dcsl + FP_LITERAL(1.0) * dJ_leak_dcsl
                     - FP_LITERAL(1.0) * dJ_sl_b_dcsl - FP_LITERAL(1.0) * dJ_sl_c_dcsl
                     + FP_LITERAL(1.0)
                               * (dI_NaCa_dKa_sl * dKa_sl_dcsl + dI_NaCa_ds2_sl * ds2_sl_dcsl
                                  + dI_NaCa_ds3_sl * ds3_sl_dcsl)
                               * dJ_NaCa_dI_NaCa
                     + FP_LITERAL(1.0) * dI_pCa_dcsl * dJ_pCa_dI_pCa
                     + FP_LITERAL(1.0) * dI_bCa_deca_sl * dJ_bCa_dI_bCa * deca_sl_dcsl)
                    / Vsl;
            states[STATE_csl * padded_num_cells + i] =
                    (fabs(dcsl_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dcsl_dt_linearized)) * dcsl_dt
                                       / dcsl_dt_linearized
                             : dt * dcsl_dt)
                    + csl;
            const cellmodel_float_t dcs_dt =
                    (FP_LITERAL(1.0) * J_n_s - FP_LITERAL(1.0) * J_RyR - FP_LITERAL(1.0) * J_s_b)
                    / Vs;
            const cellmodel_float_t dJ_RyR_active_dcs =
                    alpha_RyR * lambda_RyR * lambda_c_i * p * r_RyR;
            const cellmodel_float_t dJ_leak_dcs = alpha_RyR * gamma_RyR * lambda_RyR * lambda_c_i;
            const cellmodel_float_t dJ_n_s_dcs = -alpha_n_s * lambda_c_i * lambda_diff;
            const cellmodel_float_t dJ_s_b_dcs = Vs * k_on_s * (-bs + B_tot_s * lambda_B);
            const cellmodel_float_t dcs_dt_linearized =
                    (FP_LITERAL(1.0) * dJ_n_s_dcs - FP_LITERAL(1.0) * dJ_RyR_active_dcs
                     - FP_LITERAL(1.0) * dJ_leak_dcs - FP_LITERAL(1.0) * dJ_s_b_dcs)
                    / Vs;
            states[STATE_cs * padded_num_cells + i] =
                    (fabs(dcs_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dcs_dt_linearized)) * dcs_dt
                                       / dcs_dt_linearized
                             : dt * dcs_dt)
                    + cs;

            // Expressions for the Ca Buffer Concentrations component
            const cellmodel_float_t dbc_dt = FP_LITERAL(1.0) * J_c_b / Vc;
            const cellmodel_float_t dJ_c_b_dbc = Vc * (-k_off_c - k_on_c * cc);
            const cellmodel_float_t dbc_dt_linearized = FP_LITERAL(1.0) * dJ_c_b_dbc / Vc;
            states[STATE_bc * padded_num_cells + i] =
                    (fabs(dbc_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dbc_dt_linearized)) * dbc_dt
                                       / dbc_dt_linearized
                             : dt * dbc_dt)
                    + bc;
            const cellmodel_float_t dbd_dt = FP_LITERAL(1.0) * J_d_b / Vd;
            const cellmodel_float_t dJ_d_b_dbd = Vd * (-k_off_d - k_on_d * cd);
            const cellmodel_float_t dbd_dt_linearized = FP_LITERAL(1.0) * dJ_d_b_dbd / Vd;
            states[STATE_bd * padded_num_cells + i] =
                    (fabs(dbd_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dbd_dt_linearized)) * dbd_dt
                                       / dbd_dt_linearized
                             : dt * dbd_dt)
                    + bd;
            const cellmodel_float_t dbs_dt = FP_LITERAL(1.0) * J_s_b / Vs;
            const cellmodel_float_t dJ_s_b_dbs = Vs * (-k_off_s - k_on_s * cs);
            const cellmodel_float_t dbs_dt_linearized = FP_LITERAL(1.0) * dJ_s_b_dbs / Vs;
            states[STATE_bs * padded_num_cells + i] =
                    (fabs(dbs_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dbs_dt_linearized)) * dbs_dt
                                       / dbs_dt_linearized
                             : dt * dbs_dt)
                    + bs;
            const cellmodel_float_t dbsl_dt = FP_LITERAL(1.0) * J_sl_b / Vsl;
            const cellmodel_float_t dJ_sl_b_dbsl = Vsl * (-k_off_sl - k_on_sl * csl);
            const cellmodel_float_t dbsl_dt_linearized = FP_LITERAL(1.0) * dJ_sl_b_dbsl / Vsl;
            states[STATE_bsl * padded_num_cells + i] =
                    (fabs(dbsl_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dbsl_dt_linearized)) * dbsl_dt
                                       / dbsl_dt_linearized
                             : dt * dbsl_dt)
                    + bsl;

            // Expressions for the Membrane potential component
            const cellmodel_float_t i_Stim =
                    (V_m < FP_LITERAL(-40.) ? FP_LITERAL(1.) : FP_LITERAL(0.))
                    * (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                       && t - stim_period * floor(t / stim_period) >= stim_start
                               ? -stim_amplitude
                               : FP_LITERAL(0.));
            const cellmodel_float_t I_tot = I_CaL + I_K1 + I_KATP + I_Kr + I_Ks + I_Na + I_NaCa
                                            + I_NaK + I_NaL + I_bCa + I_bCl + I_f + I_pCa + I_to;
            const cellmodel_float_t dV_m_dt = -I_tot - i_Stim;
            const cellmodel_float_t dI_K1_dK1ss =
                    FP_LITERAL(0.430331482912) * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m);
            const cellmodel_float_t dK1ss_daK1 =
                    FP_LITERAL(1.0) / (aK1 + bK1) - aK1 / ((aK1 + bK1) * (aK1 + bK1));
            const cellmodel_float_t dK1ss_dbK1 = -aK1 / ((aK1 + bK1) * (aK1 + bK1));
            const cellmodel_float_t daK1_dV_m =
                    FP_LITERAL(-1.500911583016e-6)
                    * Exp(FP_LITERAL(0.2) * V_m - FP_LITERAL(0.2) * ek)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(7.50455791508e-6)
                                  * Exp(FP_LITERAL(0.2) * V_m - FP_LITERAL(0.2) * ek))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(7.50455791508e-6)
                                    * Exp(FP_LITERAL(0.2) * V_m - FP_LITERAL(0.2) * ek)));
            const cellmodel_float_t dbK1_dV_m =
                    (FP_LITERAL(1.994784180198e-17)
                             * Exp(FP_LITERAL(0.06) * V_m - FP_LITERAL(0.06) * ek)
                     + FP_LITERAL(0.05967298790568)
                               * Exp(FP_LITERAL(0.08) * V_m - FP_LITERAL(0.08) * ek))
                            / (FP_LITERAL(1.)
                               + FP_LITERAL(0.0820849986239)
                                         * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m))
                    + FP_LITERAL(0.04104249931195)
                              * (FP_LITERAL(0.745912348821)
                                         * Exp(FP_LITERAL(0.08) * V_m - FP_LITERAL(0.08) * ek)
                                 + FP_LITERAL(3.32464030033e-16)
                                           * Exp(FP_LITERAL(0.06) * V_m - FP_LITERAL(0.06) * ek))
                              * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m)
                              / ((FP_LITERAL(1.)
                                  + FP_LITERAL(0.0820849986239)
                                            * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m))
                                 * (FP_LITERAL(1.)
                                    + FP_LITERAL(0.0820849986239)
                                              * Exp(FP_LITERAL(0.5) * ek - FP_LITERAL(0.5) * V_m)));
            const cellmodel_float_t dI_K1_dV_m =
                    FP_LITERAL(0.430331482912) * g_K1 * lambda_K * sqrt(Ko) * K1ss
                    + FP_LITERAL(0.430331482912) * g_K1 * lambda_K * sqrt(Ko) * (-ek + V_m)
                              * (dK1ss_daK1 * daK1_dV_m + dK1ss_dbK1 * dbK1_dV_m);
            const cellmodel_float_t dI_KATP_dV_m =
                    FP_LITERAL(0.60295079490657) * g_KATP * pow(Ko, FP_LITERAL(0.3))
                            / (FP_LITERAL(40.) + FP_LITERAL(0.0875) * V_m)
                    - FP_LITERAL(0.0527581945543249) * g_KATP * pow(Ko, FP_LITERAL(0.3))
                              * (-ek + V_m)
                              / ((FP_LITERAL(40.) + FP_LITERAL(0.0875) * V_m)
                                 * (FP_LITERAL(40.) + FP_LITERAL(0.0875) * V_m));
            const cellmodel_float_t dI_Kr_dV_m =
                    FP_LITERAL(0.430331482912) * g_Kr * sqrt(Ko) * Xr1 * Xr2;
            const cellmodel_float_t dI_Ks_dV_m = g_Ks * (x_Ks * x_Ks);
            const cellmodel_float_t dI_Na_dV_m = g_Na * lambda_Na * (m * m * m) * j;
            const cellmodel_float_t ds1_sl_dV_m =
                    ce * nu * (Na_i * Na_i * Na_i) * FoRT * Exp(nu * FoRT * V_m);
            const cellmodel_float_t ds2_sl_dV_m = (Nao * Nao * Nao) * (FP_LITERAL(-1.) + nu) * FoRT
                                                  * csl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t dI_NaCa_dV_m =
                    g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * (-ds2_sl_dV_m + ds1_sl_dV_m) * Ka_sl
                            / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                               * s3_sl)
                    - g_NaCa * ksat * lambda_c_e * pow(Q10NCX, Qpow) * (FP_LITERAL(-1.) + nu)
                              * (-s2_sl + s1_sl) * FoRT * Ka_sl
                              * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)
                              / (((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                                  * (FP_LITERAL(1.)
                                     + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)))
                                 * s3_sl);
            const cellmodel_float_t dI_NaCa_ds1_sl =
                    g_NaCa * lambda_c_e * pow(Q10NCX, Qpow) * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t dI_NaK_dfNaK =
                    Ko * g_NaK
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_i) * (Na_i)) * ((Na_i) * (Na_i))))
                       * (KmKo + Ko));
            const cellmodel_float_t dI_NaL_dV_m = lambda_Na * GNaL * hL * mL;
            const cellmodel_float_t dI_to_dV_m = g_to * q * r;
            const cellmodel_float_t dfNaK_dV_m =
                    FP_LITERAL(1.0)
                    * (FP_LITERAL(0.012) * FoRT * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                       + FP_LITERAL(0.037) * FoRT * Exp(-FoRT * V_m) * sigma)
                    / ((FP_LITERAL(1.) + FP_LITERAL(0.12) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                        + FP_LITERAL(0.037) * Exp(-FoRT * V_m) * sigma)
                       * (FP_LITERAL(1.) + FP_LITERAL(0.12) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                          + FP_LITERAL(0.037) * Exp(-FoRT * V_m) * sigma));
            const cellmodel_float_t dibarca_j_dV_m =
                    FP_LITERAL(4.) * Frdy * g_CaL
                            * (FP_LITERAL(-0.34) * ce
                               + FP_LITERAL(0.34) * cd * Exp(FP_LITERAL(2.) * FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                    - FP_LITERAL(8.) * Frdy * g_CaL * (FoRT * FoRT)
                              * (FP_LITERAL(-0.34) * ce
                                 + FP_LITERAL(0.34) * cd * Exp(FP_LITERAL(2.) * FoRT * V_m))
                              * V_m * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m)))
                    + FP_LITERAL(2.72) * Frdy * g_CaL * (FoRT * FoRT) * V_m * cd
                              * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t dV_m_dt_linearized =
                    -g_bCl - dI_K1_dV_m - dI_KATP_dV_m - dI_Kr_dV_m - dI_Ks_dV_m - dI_NaCa_dV_m
                    - dI_NaL_dV_m - dI_Na_dV_m - dI_to_dV_m - g_bCa * lambda_c_e - g_f * xf
                    - (dK1ss_daK1 * daK1_dV_m + dK1ss_dbK1 * dbK1_dV_m) * dI_K1_dK1ss
                    - dI_CaL_dibarca_j * dibarca_j_dV_m - dI_NaCa_ds1_sl * ds1_sl_dV_m
                    - dI_NaCa_ds2_sl * ds2_sl_dV_m - dI_NaK_dfNaK * dfNaK_dV_m;
            states[STATE_V_m * padded_num_cells + i] =
                    (fabs(dV_m_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dV_m_dt_linearized)) * dV_m_dt
                                       / dV_m_dt_linearized
                             : dt * dV_m_dt)
                    + V_m;

            // Expressions for the Sodium concentration component
            const cellmodel_float_t I_Na_tot = FP_LITERAL(3.) * I_NaCa + FP_LITERAL(3.) * I_NaK
                                               + FP_LITERAL(0.3293) * I_f + I_Na + I_NaL;
            const cellmodel_float_t J_Na = -Cm * chi * I_Na_tot / Frdy;
            const cellmodel_float_t dNa_i_dt = J_Na;
            const cellmodel_float_t dI_Na_dena = -g_Na * lambda_Na * (m * m * m) * j;
            const cellmodel_float_t dI_NaK_dNa_i =
                    FP_LITERAL(4.) * Ko * g_NaK * (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                    * fNaK
                    / (((1.
                         + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                   / (((Na_i) * (Na_i)) * ((Na_i) * (Na_i))))
                        * (1.
                           + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                     / (((Na_i) * (Na_i)) * ((Na_i) * (Na_i)))))
                       * (KmKo + Ko) * pow(Na_i, FP_LITERAL(5.)));
            const cellmodel_float_t dI_NaL_dena = -lambda_Na * GNaL * hL * mL;
            const cellmodel_float_t dJ_Na_dI_Na_tot = -Cm * chi / Frdy;
            const cellmodel_float_t dena_dNa_i = FP_LITERAL(-1.) / (FoRT * Na_i);
            const cellmodel_float_t ds1_sl_dNa_i =
                    FP_LITERAL(3.) * ce * (Na_i * Na_i) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t ds3_sl_dNa_i =
                    FP_LITERAL(3.) * KmCao * (Na_i * Na_i) + FP_LITERAL(3.) * ce * (Na_i * Na_i)
                    + FP_LITERAL(3.) * KmCai * (Nao * Nao * Nao) * (Na_i * Na_i)
                              / (KmNai * KmNai * KmNai);
            const cellmodel_float_t dNa_i_dt_linearized =
                    (FP_LITERAL(3.) * dI_NaK_dNa_i + dI_NaL_dena * dena_dNa_i
                     + dI_Na_dena * dena_dNa_i + FP_LITERAL(3.) * dI_NaCa_ds1_sl * ds1_sl_dNa_i
                     + FP_LITERAL(3.) * dI_NaCa_ds3_sl * ds3_sl_dNa_i)
                    * dJ_Na_dI_Na_tot;
            states[STATE_Na_i * padded_num_cells + i] =
                    Na_i
                    + (fabs(dNa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_i_dt_linearized)) * dNa_i_dt
                                         / dNa_i_dt_linearized
                               : dt * dNa_i_dt);
        }
        t += dt;
    }
}

//static const int num_colour_sets = 4;
#define NUM_COLOUR_SETS 4
static const uint8_t colour_set_map[] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                         0, 1, 0, 1, 2, 3, 2, 2, 2, 0, 0, 1};

const struct state_colour_sets colour_sets = {NUM_COLOUR_SETS, colour_set_map};

// clang-format off
const struct cellmodel model_JT21_multistep_time_cell = {
        .init_states = &init_state_values,
        .init_parameters = &init_parameters_values,
        .state_index = &state_index,
        .parameter_index = &parameter_index,
        .step_FE = NULL,
        .step_RL = NULL,
        .step_GRL1 = NULL,
        .multistep_FE = &multistep_FE,
        .multistep_GRL1 = &multistep_GRL1,
        .num_states = NUM_STATES,
        .num_parameters = NUM_PARAMS,
        .layout = LAYOUT_STRUCT_OF_ARRAYS,
        .colour_sets = &colour_sets,
};
// clang-format on
