#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "base_model" model

#include <array>

//#include "cellmodel_lut.hpp"
#include "TP06_lut.hpp"


enum state {
    STATE_Xr1,
    STATE_Xr2,
    STATE_Xs,
    STATE_m,
    STATE_h,
    STATE_j,
    STATE_d,
    STATE_f,
    STATE_f2,
    STATE_fCass,
    STATE_s,
    STATE_r,
    STATE_Ca_i,
    STATE_R_prime,
    STATE_Ca_SR,
    STATE_Ca_ss,
    STATE_Na_i,
    STATE_V,
    STATE_K_i,
    NUM_STATES,
};

enum parameter {
    PARAM_celltype,
    PARAM_P_kna,
    PARAM_g_K1,
    PARAM_g_Kr,
    PARAM_g_Ks,
    PARAM_g_Na,
    PARAM_g_bna,
    PARAM_g_CaL,
    PARAM_i_CaL_lim_delta,
    PARAM_g_bca,
    PARAM_g_to,
    PARAM_K_mNa,
    PARAM_K_mk,
    PARAM_P_NaK,
    PARAM_K_NaCa,
    PARAM_K_sat,
    PARAM_Km_Ca,
    PARAM_Km_Nai,
    PARAM_alpha,
    PARAM_gamma,
    PARAM_K_pCa,
    PARAM_g_pCa,
    PARAM_g_pK,
    PARAM_Buf_c,
    PARAM_Buf_sr,
    PARAM_Buf_ss,
    PARAM_Ca_o,
    PARAM_EC,
    PARAM_K_buf_c,
    PARAM_K_buf_sr,
    PARAM_K_buf_ss,
    PARAM_K_up,
    PARAM_V_leak,
    PARAM_V_rel,
    PARAM_V_sr,
    PARAM_V_ss,
    PARAM_V_xfer,
    PARAM_Vmax_up,
    PARAM_k1_prime,
    PARAM_k2_prime,
    PARAM_k3,
    PARAM_k4,
    PARAM_max_sr,
    PARAM_min_sr,
    PARAM_Na_o,
    PARAM_Cm,
    PARAM_F,
    PARAM_R,
    PARAM_T,
    PARAM_V_c,
    PARAM_is_stimulated,
    PARAM_stim_amplitude,
    PARAM_K_o,
    NUM_PARAMS,
};

// Init state values
static void init_state_values(cellmodel_float_t *__restrict states, const long num_cells,
                              const long padded_num_cells)
{
    #pragma omp parallel for
    for (long i = 0; i < num_cells; i++) {
        states[STATE_Xr1 * padded_num_cells + i] = FP_LITERAL(0.0165);
        states[STATE_Xr2 * padded_num_cells + i] = FP_LITERAL(0.473);
        states[STATE_Xs * padded_num_cells + i] = FP_LITERAL(0.0174);
        states[STATE_m * padded_num_cells + i] = FP_LITERAL(0.00165);
        states[STATE_h * padded_num_cells + i] = FP_LITERAL(0.749);
        states[STATE_j * padded_num_cells + i] = FP_LITERAL(0.6788);
        states[STATE_d * padded_num_cells + i] = FP_LITERAL(3.288e-05);
        states[STATE_f * padded_num_cells + i] = FP_LITERAL(0.7026);
        states[STATE_f2 * padded_num_cells + i] = FP_LITERAL(0.9526);
        states[STATE_fCass * padded_num_cells + i] = FP_LITERAL(0.9942);
        states[STATE_s * padded_num_cells + i] = FP_LITERAL(0.999998);
        states[STATE_r * padded_num_cells + i] = FP_LITERAL(2.347e-08);
        states[STATE_Ca_i * padded_num_cells + i] = FP_LITERAL(0.000153);
        states[STATE_R_prime * padded_num_cells + i] = FP_LITERAL(0.8978);
        states[STATE_Ca_SR * padded_num_cells + i] = FP_LITERAL(4.272);
        states[STATE_Ca_ss * padded_num_cells + i] = FP_LITERAL(0.00042);
        states[STATE_Na_i * padded_num_cells + i] = FP_LITERAL(10.132);
        states[STATE_V * padded_num_cells + i] = -FP_LITERAL(85.423);
        states[STATE_K_i * padded_num_cells + i] = FP_LITERAL(138.52);
    }
}

// Default parameter values
static void init_parameters_values(cellmodel_float_t *__restrict parameters)
{
    parameters[PARAM_celltype] = 0;
    parameters[PARAM_P_kna] = FP_LITERAL(0.03);
    parameters[PARAM_g_K1] = FP_LITERAL(5.405);
    parameters[PARAM_g_Kr] = FP_LITERAL(0.153);
    parameters[PARAM_g_Ks] = FP_LITERAL(0.098);
    parameters[PARAM_g_Na] = FP_LITERAL(14.838);
    parameters[PARAM_g_bna] = FP_LITERAL(0.00029);
    parameters[PARAM_g_CaL] = FP_LITERAL(3.98e-05);
    parameters[PARAM_i_CaL_lim_delta] = 1e-07;
    parameters[PARAM_g_bca] = FP_LITERAL(0.000592);
    parameters[PARAM_g_to] = FP_LITERAL(0.294);
    parameters[PARAM_K_mNa] = 40;
    parameters[PARAM_K_mk] = 1;
    parameters[PARAM_P_NaK] = FP_LITERAL(2.724);
    parameters[PARAM_K_NaCa] = 1000;
    parameters[PARAM_K_sat] = FP_LITERAL(0.1);
    parameters[PARAM_Km_Ca] = FP_LITERAL(1.38);
    parameters[PARAM_Km_Nai] = FP_LITERAL(87.5);
    parameters[PARAM_alpha] = FP_LITERAL(2.5);
    parameters[PARAM_gamma] = FP_LITERAL(0.35);
    parameters[PARAM_K_pCa] = FP_LITERAL(0.0005);
    parameters[PARAM_g_pCa] = FP_LITERAL(0.1238);
    parameters[PARAM_g_pK] = FP_LITERAL(0.0146);
    parameters[PARAM_Buf_c] = FP_LITERAL(0.2);
    parameters[PARAM_Buf_sr] = 10;
    parameters[PARAM_Buf_ss] = FP_LITERAL(0.4);
    parameters[PARAM_Ca_o] = 2;
    parameters[PARAM_EC] = FP_LITERAL(1.5);
    parameters[PARAM_K_buf_c] = FP_LITERAL(0.001);
    parameters[PARAM_K_buf_sr] = FP_LITERAL(0.3);
    parameters[PARAM_K_buf_ss] = FP_LITERAL(0.00025);
    parameters[PARAM_K_up] = FP_LITERAL(0.00025);
    parameters[PARAM_V_leak] = FP_LITERAL(0.00036);
    parameters[PARAM_V_rel] = FP_LITERAL(0.102);
    parameters[PARAM_V_sr] = FP_LITERAL(0.001094);
    parameters[PARAM_V_ss] = FP_LITERAL(5.468e-05);
    parameters[PARAM_V_xfer] = FP_LITERAL(0.0038);
    parameters[PARAM_Vmax_up] = FP_LITERAL(0.006375);
    parameters[PARAM_k1_prime] = FP_LITERAL(0.15);
    parameters[PARAM_k2_prime] = FP_LITERAL(0.045);
    parameters[PARAM_k3] = FP_LITERAL(0.06);
    parameters[PARAM_k4] = FP_LITERAL(0.005);
    parameters[PARAM_max_sr] = FP_LITERAL(2.5);
    parameters[PARAM_min_sr] = FP_LITERAL(1.0);
    parameters[PARAM_Na_o] = 140;
    parameters[PARAM_Cm] = FP_LITERAL(0.185);
    parameters[PARAM_F] = FP_LITERAL(96485.3415);
    parameters[PARAM_R] = FP_LITERAL(8314.472);
    parameters[PARAM_T] = 310;
    parameters[PARAM_V_c] = FP_LITERAL(0.016404);
    parameters[PARAM_is_stimulated] = 0;
    parameters[PARAM_stim_amplitude] = 52;
    parameters[PARAM_K_o] = FP_LITERAL(5.4);
}

// State index
static int state_index(const char name[])
{
    if (strcmp(name, "Xr1") == 0) {
        return STATE_Xr1;
    } else if (strcmp(name, "Xr2") == 0) {
        return STATE_Xr2;
    } else if (strcmp(name, "Xs") == 0) {
        return STATE_Xs;
    } else if (strcmp(name, "m") == 0) {
        return STATE_m;
    } else if (strcmp(name, "h") == 0) {
        return STATE_h;
    } else if (strcmp(name, "j") == 0) {
        return STATE_j;
    } else if (strcmp(name, "d") == 0) {
        return STATE_d;
    } else if (strcmp(name, "f") == 0) {
        return STATE_f;
    } else if (strcmp(name, "f2") == 0) {
        return STATE_f2;
    } else if (strcmp(name, "fCass") == 0) {
        return STATE_fCass;
    } else if (strcmp(name, "s") == 0) {
        return STATE_s;
    } else if (strcmp(name, "r") == 0) {
        return STATE_r;
    } else if (strcmp(name, "Ca_i") == 0) {
        return STATE_Ca_i;
    } else if (strcmp(name, "R_prime") == 0) {
        return STATE_R_prime;
    } else if (strcmp(name, "Ca_SR") == 0) {
        return STATE_Ca_SR;
    } else if (strcmp(name, "Ca_ss") == 0) {
        return STATE_Ca_ss;
    } else if (strcmp(name, "Na_i") == 0) {
        return STATE_Na_i;
    } else if (strcmp(name, "V") == 0) {
        return STATE_V;
    } else if (strcmp(name, "K_i") == 0) {
        return STATE_K_i;
    }
    return -1;
}

// Parameter index
static int parameter_index(const char name[])
{
    if (strcmp(name, "celltype") == 0) {
        return PARAM_celltype;
    } else if (strcmp(name, "P_kna") == 0) {
        return PARAM_P_kna;
    } else if (strcmp(name, "g_K1") == 0) {
        return PARAM_g_K1;
    } else if (strcmp(name, "g_Kr") == 0) {
        return PARAM_g_Kr;
    } else if (strcmp(name, "g_Ks") == 0) {
        return PARAM_g_Ks;
    } else if (strcmp(name, "g_Na") == 0) {
        return PARAM_g_Na;
    } else if (strcmp(name, "g_bna") == 0) {
        return PARAM_g_bna;
    } else if (strcmp(name, "g_CaL") == 0) {
        return PARAM_g_CaL;
    } else if (strcmp(name, "i_CaL_lim_delta") == 0) {
        return PARAM_i_CaL_lim_delta;
    } else if (strcmp(name, "g_bca") == 0) {
        return PARAM_g_bca;
    } else if (strcmp(name, "g_to") == 0) {
        return PARAM_g_to;
    } else if (strcmp(name, "K_mNa") == 0) {
        return PARAM_K_mNa;
    } else if (strcmp(name, "K_mk") == 0) {
        return PARAM_K_mk;
    } else if (strcmp(name, "P_NaK") == 0) {
        return PARAM_P_NaK;
    } else if (strcmp(name, "K_NaCa") == 0) {
        return PARAM_K_NaCa;
    } else if (strcmp(name, "K_sat") == 0) {
        return PARAM_K_sat;
    } else if (strcmp(name, "Km_Ca") == 0) {
        return PARAM_Km_Ca;
    } else if (strcmp(name, "Km_Nai") == 0) {
        return PARAM_Km_Nai;
    } else if (strcmp(name, "alpha") == 0) {
        return PARAM_alpha;
    } else if (strcmp(name, "gamma") == 0) {
        return PARAM_gamma;
    } else if (strcmp(name, "K_pCa") == 0) {
        return PARAM_K_pCa;
    } else if (strcmp(name, "g_pCa") == 0) {
        return PARAM_g_pCa;
    } else if (strcmp(name, "g_pK") == 0) {
        return PARAM_g_pK;
    } else if (strcmp(name, "Buf_c") == 0) {
        return PARAM_Buf_c;
    } else if (strcmp(name, "Buf_sr") == 0) {
        return PARAM_Buf_sr;
    } else if (strcmp(name, "Buf_ss") == 0) {
        return PARAM_Buf_ss;
    } else if (strcmp(name, "Ca_o") == 0) {
        return PARAM_Ca_o;
    } else if (strcmp(name, "EC") == 0) {
        return PARAM_EC;
    } else if (strcmp(name, "K_buf_c") == 0) {
        return PARAM_K_buf_c;
    } else if (strcmp(name, "K_buf_sr") == 0) {
        return PARAM_K_buf_sr;
    } else if (strcmp(name, "K_buf_ss") == 0) {
        return PARAM_K_buf_ss;
    } else if (strcmp(name, "K_up") == 0) {
        return PARAM_K_up;
    } else if (strcmp(name, "V_leak") == 0) {
        return PARAM_V_leak;
    } else if (strcmp(name, "V_rel") == 0) {
        return PARAM_V_rel;
    } else if (strcmp(name, "V_sr") == 0) {
        return PARAM_V_sr;
    } else if (strcmp(name, "V_ss") == 0) {
        return PARAM_V_ss;
    } else if (strcmp(name, "V_xfer") == 0) {
        return PARAM_V_xfer;
    } else if (strcmp(name, "Vmax_up") == 0) {
        return PARAM_Vmax_up;
    } else if (strcmp(name, "k1_prime") == 0) {
        return PARAM_k1_prime;
    } else if (strcmp(name, "k2_prime") == 0) {
        return PARAM_k2_prime;
    } else if (strcmp(name, "k3") == 0) {
        return PARAM_k3;
    } else if (strcmp(name, "k4") == 0) {
        return PARAM_k4;
    } else if (strcmp(name, "max_sr") == 0) {
        return PARAM_max_sr;
    } else if (strcmp(name, "min_sr") == 0) {
        return PARAM_min_sr;
    } else if (strcmp(name, "Na_o") == 0) {
        return PARAM_Na_o;
    } else if (strcmp(name, "Cm") == 0) {
        return PARAM_Cm;
    } else if (strcmp(name, "F") == 0) {
        return PARAM_F;
    } else if (strcmp(name, "R") == 0) {
        return PARAM_R;
    } else if (strcmp(name, "T") == 0) {
        return PARAM_T;
    } else if (strcmp(name, "V_c") == 0) {
        return PARAM_V_c;
    } else if (strcmp(name, "is_stimulated") == 0) {
        return PARAM_is_stimulated;
    } else if (strcmp(name, "stim_amplitude") == 0) {
        return PARAM_stim_amplitude;
    } else if (strcmp(name, "K_o") == 0) {
        return PARAM_K_o;
    }
    return -1;
}

constexpr std::array<const univariate_func_tuple, 28> expressions_V = {
        univariate_func_tuple{"Xr1_RL_A",
                              [](double V, double dt, double *param) {
                                  double Xr1_inf = (1. / (1. + (exp(((-26. - (V)) / 7.)))));
                                  double aa_Xr1 = (450. / (1. + (exp(((-45. - (V)) / 10.)))));
                                  double bb_Xr1 = (6. / (1. + (exp(((V - (-30.)) / 11.5)))));
                                  double tau_Xr1 = (aa_Xr1 * bb_Xr1);
                                  double Xr1_rush_larsen_C = expm1(((-dt) / tau_Xr1));
                                  return -Xr1_inf * Xr1_rush_larsen_C;
                              }},
        univariate_func_tuple{"Xr1_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_Xr1 = (450. / (1. + (exp(((-45. - (V)) / 10.)))));
                                  double bb_Xr1 = (6. / (1. + (exp(((V - (-30.)) / 11.5)))));
                                  double tau_Xr1 = (aa_Xr1 * bb_Xr1);
                                  return exp((-dt) / tau_Xr1);
                              }},

        univariate_func_tuple{"Xr2_RL_A",
                              [](double V, double dt, double *param) {
                                  double Xr2_inf = (1. / (1. + (exp(((V - ((-88.))) / 24.)))));
                                  double aa_Xr2 = (3. / (1. + (exp(((-60. - (V)) / 20.)))));
                                  double bb_Xr2 = (1.12 / (1. + (exp(((V - (60.)) / 20.)))));
                                  double tau_Xr2 = (aa_Xr2 * bb_Xr2);
                                  double Xr2_rush_larsen_C = expm1(((-dt) / tau_Xr2));
                                  return -Xr2_inf * Xr2_rush_larsen_C;
                              }},
        univariate_func_tuple{"Xr2_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_Xr2 = (3. / (1. + (exp(((-60. - (V)) / 20.)))));
                                  double bb_Xr2 = (1.12 / (1. + (exp(((V - (60.)) / 20.)))));
                                  double tau_Xr2 = (aa_Xr2 * bb_Xr2);
                                  return exp((-dt) / tau_Xr2);
                              }},
        univariate_func_tuple{"Xs_RL_A",
                              [](double V, double dt, double *param) {
                                  double Xs_inf = (1. / (1. + (exp(((-5. - (V)) / 14.)))));
                                  double aa_Xs = (1400. / (sqrt((1. + (exp(((5. - (V)) / 6.)))))));
                                  double bb_Xs = (1. / (1. + (exp(((V - (35.)) / 15.)))));
                                  double tau_Xs = ((aa_Xs * bb_Xs) + 80.);
                                  double Xs_rush_larsen_C = (expm1(((-dt) / tau_Xs)));
                                  return (-Xs_inf) * Xs_rush_larsen_C;
                              }},
        univariate_func_tuple{"Xs_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_Xs = (1400. / (sqrt((1. + (exp(((5. - (V)) / 6.)))))));
                                  double bb_Xs = (1. / (1. + (exp(((V - (35.)) / 15.)))));
                                  double tau_Xs = ((aa_Xs * bb_Xs) + 80.);
                                  return exp((-dt) / tau_Xs);
                              }},
        univariate_func_tuple{"m_RL_A",
                              [](double V, double dt, double *param) {
                                  double M_inf = (1.
                                                  / ((1. + (exp(((-56.86 - (V)) / 9.03))))
                                                     * (1. + (exp(((-56.86 - (V)) / 9.03))))));
                                  double aa_M = (1. / (1. + (exp(((-60. - (V)) / 5.)))));
                                  double bb_M = ((0.1 / (1. + (exp(((V + 35.) / 5.)))))
                                                 + (0.10 / (1. + (exp(((V - (50.)) / 200.))))));
                                  double tau_M = (aa_M * bb_M);
                                  double M_rush_larsen_C = (expm1(((-dt) / tau_M)));
                                  return (-M_inf) * M_rush_larsen_C;
                              }},
        univariate_func_tuple{"m_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_M = (1. / (1. + (exp(((-60. - (V)) / 5.)))));
                                  double bb_M = ((0.1 / (1. + (exp(((V + 35.) / 5.)))))
                                                 + (0.10 / (1. + (exp(((V - (50.)) / 200.))))));
                                  double tau_M = (aa_M * bb_M);
                                  return exp((-dt) / tau_M);
                              }},

        univariate_func_tuple{
                "h_RL_A",
                [](double V, double dt, double *param) {
                    double H_inf = (1.
                                    / ((1. + (exp(((V + 71.55) / 7.43))))
                                       * (1. + (exp(((V + 71.55) / 7.43))))));
                    double aa_H = ((V >= -40.) ? 0. : (0.057 * (exp(((-(V + 80.)) / 6.8)))));
                    double bb_H =
                            ((V >= -40.) ? (0.77 / (0.13 * (1. + (exp(((-(V + 10.66)) / 11.1))))))
                                         : ((2.7 * (exp((0.079 * V))))
                                            + (3.1e5 * (exp((0.3485 * V))))));
                    double tau_H = (1.0 / (aa_H + bb_H));
                    double H_rush_larsen_C = (expm1(((-dt) / tau_H)));
                    return (-H_inf) * H_rush_larsen_C;
                }},
        univariate_func_tuple{
                "h_RL_B",
                [](double V, double dt, double *param) {
                    double aa_H = ((V >= -40.) ? 0. : (0.057 * (exp(((-(V + 80.)) / 6.8)))));
                    double bb_H =
                            ((V >= -40.) ? (0.77 / (0.13 * (1. + (exp(((-(V + 10.66)) / 11.1))))))
                                         : ((2.7 * (exp((0.079 * V))))
                                            + (3.1e5 * (exp((0.3485 * V))))));
                    double tau_H = (1.0 / (aa_H + bb_H));
                    return exp((-dt) / tau_H);
                }},
        univariate_func_tuple{"j_RL_A",
                              [](double V, double dt, double *param) {
                                  double J_inf = (1.
                                                  / ((1. + (exp(((V + 71.55) / 7.43))))
                                                     * (1. + (exp(((V + 71.55) / 7.43))))));
                                  double aa_J =
                                          ((V >= -40.) ? 0.
                                                       : ((((-2.5428e4 * (exp((0.2444 * V))))
                                                            - ((6.948e-6 * (exp((-0.04391 * V))))))
                                                           * (V + 37.78))
                                                          / (1. + (exp((0.311 * (V + 79.23)))))));
                                  double bb_J =
                                          ((V >= -40.) ? ((0.6 * (exp((0.057 * V))))
                                                          / (1. + (exp((-0.1 * (V + 32.))))))
                                                       : ((0.02424 * (exp((-0.01052 * V))))
                                                          / (1. + (exp((-0.1378 * (V + 40.14)))))));
                                  double tau_J = (1.0 / (aa_J + bb_J));
                                  double J_rush_larsen_C = (expm1(((-dt) / tau_J)));
                                  return (-J_inf) * J_rush_larsen_C;
                              }},
        univariate_func_tuple{"j_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_J =
                                          ((V >= -40.) ? 0.
                                                       : ((((-2.5428e4 * (exp((0.2444 * V))))
                                                            - ((6.948e-6 * (exp((-0.04391 * V))))))
                                                           * (V + 37.78))
                                                          / (1. + (exp((0.311 * (V + 79.23)))))));
                                  double bb_J =
                                          ((V >= -40.) ? ((0.6 * (exp((0.057 * V))))
                                                          / (1. + (exp((-0.1 * (V + 32.))))))
                                                       : ((0.02424 * (exp((-0.01052 * V))))
                                                          / (1. + (exp((-0.1378 * (V + 40.14)))))));
                                  double tau_J = (1.0 / (aa_J + bb_J));
                                  return exp((-dt) / tau_J);
                              }},
        univariate_func_tuple{"d_RL_A",
                              [](double V, double dt, double *param) {
                                  double D_inf = (1. / (1. + (exp((((-8.) - (V)) / 7.5)))));
                                  double aa_D = ((1.4 / (1. + (exp(((-35. - (V)) / 13.))))) + 0.25);
                                  double bb_D = (1.4 / (1. + (exp(((V + 5.) / 5.)))));
                                  double cc_D = (1. / (1. + (exp(((50. - (V)) / 20.)))));
                                  double tau_D = ((aa_D * bb_D) + cc_D);
                                  double D_rush_larsen_C = (expm1(((-dt) / tau_D)));
                                  return (-D_inf) * D_rush_larsen_C;
                              }},
        univariate_func_tuple{"d_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_D = ((1.4 / (1. + (exp(((-35. - (V)) / 13.))))) + 0.25);
                                  double bb_D = (1.4 / (1. + (exp(((V + 5.) / 5.)))));
                                  double cc_D = (1. / (1. + (exp(((50. - (V)) / 20.)))));
                                  double tau_D = ((aa_D * bb_D) + cc_D);
                                  return exp((-dt) / tau_D);
                              }},
        univariate_func_tuple{"f_RL_A",
                              [](double V, double dt, double *param) {
                                  double F_inf = (1. / (1. + (exp(((V + 20.) / 7.)))));
                                  double aa_F =
                                          (1102.5 * (exp((((-(V + 27.)) * (V + 27.)) / 225.))));
                                  double bb_F = (200. / (1. + (exp(((13. - (V)) / 10.)))));

                                  double cc_F = ((180. / (1. + (exp(((V + 30.) / 10.))))) + 20.);
                                  double tau_F = ((aa_F + bb_F) + cc_F);
                                  double F_rush_larsen_C = (expm1(((-dt) / tau_F)));
                                  return (-F_inf) * F_rush_larsen_C;
                              }},
        univariate_func_tuple{"f_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_F =
                                          (1102.5 * (exp((((-(V + 27.)) * (V + 27.)) / 225.))));
                                  double bb_F = (200. / (1. + (exp(((13. - (V)) / 10.)))));
                                  double cc_F = ((180. / (1. + (exp(((V + 30.) / 10.))))) + 20.);
                                  double tau_F = ((aa_F + bb_F) + cc_F);
                                  return exp((-dt) / tau_F);
                              }},
        univariate_func_tuple{"f2_RL_A",
                              [](double V, double dt, double *param) {
                                  double F2_inf = ((0.67 / (1. + (exp(((V + 35.) / 7.))))) + 0.33);
                                  double aa_F2 =
                                          (562. * (exp((((-(V + 27.)) * (V + 27.)) / 240.))));
                                  double bb_F2 = (31. / (1. + (exp(((25. - (V)) / 10.)))));
                                  double cc_F2 = (80. / (1. + (exp(((V + 30.) / 10.)))));
                                  double tau_F2 = ((aa_F2 + bb_F2) + cc_F2);
                                  double F2_rush_larsen_C = (expm1(((-dt) / tau_F2)));
                                  return (-F2_inf) * F2_rush_larsen_C;
                              }},
        univariate_func_tuple{"f2_RL_B",
                              [](double V, double dt, double *param) {
                                  double aa_F2 =
                                          (562. * (exp((((-(V + 27.)) * (V + 27.)) / 240.))));
                                  double bb_F2 = (31. / (1. + (exp(((25. - (V)) / 10.)))));
                                  double cc_F2 = (80. / (1. + (exp(((V + 30.) / 10.)))));
                                  double tau_F2 = ((aa_F2 + bb_F2) + cc_F2);
                                  return exp((-dt) / tau_F2);
                              }},
        univariate_func_tuple{
                "r_RL_A",
                [](double V, double dt, double *param) {
                    double R_inf = (1. / (1. + (exp(((20. - (V)) / 6.)))));
                    double tau_R = ((9.5 * (exp((((-(V + 40.)) * (V + 40.)) / 1800.)))) + 0.8);
                    double R_rush_larsen_C = (expm1(((-dt) / tau_R)));
                    return (-R_inf) * R_rush_larsen_C;
                }},
        univariate_func_tuple{
                "r_RL_B",
                [](double V, double dt, double *param) {
                    double tau_R = ((9.5 * (exp((((-(V + 40.)) * (V + 40.)) / 1800.)))) + 0.8);
                    return exp((-dt) / tau_R);
                }},
        univariate_func_tuple{
                "s_RL_A",
                [](double V, double dt, double *param) {
                    bool is_endo = param[PARAM_celltype] == 2;
                    double S_inf = (1. / (1. + (exp(((V + 20.) / 5.)))));
                    double tau_S =
                            ((is_endo)
                                     ? ((1000. * (exp((((-(V + 67.)) * (V + 67.)) / 1000.)))) + 8.)
                                     : (((85. * (exp((((-(V + 45.)) * (V + 45.)) / 320.))))
                                         + (5. / (1. + (exp(((V - (20.)) / 5.))))))
                                        + 3.));
                    double S_rush_larsen_C = (expm1(((-dt) / tau_S)));
                    return (-S_inf) * S_rush_larsen_C;
                }},
        univariate_func_tuple{
                "s_RL_B",
                [](double V, double dt, double *param) {
                    bool is_endo = param[PARAM_celltype] == 2;
                    double tau_S =
                            ((is_endo)
                                     ? ((1000. * (exp((((-(V + 67.)) * (V + 67.)) / 1000.)))) + 8.)
                                     : (((85. * (exp((((-(V + 45.)) * (V + 45.)) / 320.))))
                                         + (5. / (1. + (exp(((V - (20.)) / 5.))))))
                                        + 3.));
                    return exp((-dt) / tau_S);
                }},

        univariate_func_tuple{
                "I_NaCa_A",
                [](double V, double dt, double *param) {
                        double F_RT = param[PARAM_F]
                                / (param[PARAM_R]
                                        * param[PARAM_T]);
                        double gamma = param[PARAM_gamma];
                        double Km_Nai = param[PARAM_Km_Nai];
                        double Na_o = param[PARAM_Na_o];
                        double invKmNai3_Nao3 = 1. / ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o));
                        double invKmCa_Cao = 1./(param[PARAM_Km_Ca]+param[PARAM_Ca_o]);
                        double pmf_INaCa = param[PARAM_K_NaCa] * invKmNai3_Nao3 * invKmCa_Cao;
                        double den =
                                (pmf_INaCa
                                / (1.
                                + (param[PARAM_K_sat]
                                        * (exp((((gamma - 1.) * V)
                                                * F_RT))))));
                        return (den * param[PARAM_Ca_o])
                                * (exp(((gamma * V)
                                        * F_RT)));
                }},
        univariate_func_tuple{
                "I_NaCa_B",
                [](double V, double dt, double *param) {
                        double F_RT = param[PARAM_F]
                                / (param[PARAM_R]
                                        * param[PARAM_T]);
                        double gamma = param[PARAM_gamma];
                        double Km_Nai = param[PARAM_Km_Nai];
                        double Na_o = param[PARAM_Na_o];
                        double Nao3 = Na_o * Na_o * Na_o;
                        double invKmNai3_Nao3 = 1. / ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o));
                        double invKmCa_Cao = 1./(param[PARAM_Km_Ca]+param[PARAM_Ca_o]);
                        double pmf_INaCa = param[PARAM_K_NaCa] * invKmNai3_Nao3 * invKmCa_Cao;
                        double den =
                                (pmf_INaCa
                                / (1.
                                + (param[PARAM_K_sat]
                                        * (exp((((gamma - 1.) * V)
                                                * F_RT))))));
                        return ((den*(exp((((gamma-(1.))*V)*F_RT))))*Nao3)*2.5;
                }},
        univariate_func_tuple{
                "i_CaL_factors_e1",
                [](double V, double dt, double *param) {
                        double F_RT = param[PARAM_F]
                                / (param[PARAM_R] * param[PARAM_T]);
                        double V_eff = V - 15;
                        return 0.25*exp(2.*V_eff*F_RT);
                }},
        univariate_func_tuple{
                "i_CaL_fraction",
                [](double V, double dt, double *param) {
                        double F = param[PARAM_F];
                        double R = param[PARAM_R];
                        double T = param[PARAM_T];
                        double i_CaL_lim_delta = param[PARAM_i_CaL_lim_delta];
                        double V_eff = V - 15;
                        return (fabs(V_eff) < i_CaL_lim_delta
                                ? 0.5
                                : F * V_eff / (R * T * (expm1(2. * F * V_eff / (R * T)))));
                }},
        univariate_func_tuple{
                "rec_iNaK",
                [](double V, double dt, double *param) {
                        double F_RT = param[PARAM_F]
                                / (param[PARAM_R]
                                        * param[PARAM_T]);
                        return 1./((1.+(0.1245*(exp((-0.1*V)*F_RT))))+(0.0353*(exp(((-V)*F_RT)))));
                }},
        univariate_func_tuple{
                "rec_ipK",
                [](double V, double dt, double *param) {
                        return 1./(1.+exp((25.-V)/5.98));
                }},
};

enum {
    LUT_INDEX_Xr1_RL_A,
    LUT_INDEX_Xr1_RL_B,
    LUT_INDEX_Xr2_RL_A,
    LUT_INDEX_Xr2_RL_B,
    LUT_INDEX_Xs_RL_A,
    LUT_INDEX_Xs_RL_B,
    LUT_INDEX_m_RL_A,
    LUT_INDEX_m_RL_B,
    LUT_INDEX_h_RL_A,
    LUT_INDEX_h_RL_B,
    LUT_INDEX_j_RL_A,
    LUT_INDEX_j_RL_B,
    LUT_INDEX_d_RL_A,
    LUT_INDEX_d_RL_B,
    LUT_INDEX_f_RL_A,
    LUT_INDEX_f_RL_B,
    LUT_INDEX_f2_RL_A,
    LUT_INDEX_f2_RL_B,
    LUT_INDEX_r_RL_A,
    LUT_INDEX_r_RL_B,
    LUT_INDEX_s_RL_A,
    LUT_INDEX_s_RL_B,
    LUT_INDEX_I_NaCa_A,
    LUT_INDEX_I_NaCa_B,
    LUT_INDEX_i_CaL_factors_e1,
    LUT_INDEX_i_CaL_fraction,
    LUT_INDEX_rec_iNaK,
    LUT_INDEX_rec_ipK,
};

constexpr std::array<const univariate_func_tuple, 2> expressions_Ca = {
        univariate_func_tuple{
                "fCass_RL_A",
                [](double CaSS, double dt, double *param) {
                        double FCaSS_inf = ((0.6/(1.+((CaSS/0.05)*(CaSS/0.05))))+0.4);
                        double tau_FCaSS = ((80./(1.+((CaSS/0.05)*(CaSS/0.05))))+2.);
                        double FCaSS_rush_larsen_C = (expm1(((-dt)/tau_FCaSS)));
                        return -FCaSS_inf*FCaSS_rush_larsen_C;
                }},
        univariate_func_tuple{
                "fCass_RL_B",
                [](double CaSS, double dt, double *param) {
                        double tau_FCaSS = ((80./(1.+((CaSS/0.05)*(CaSS/0.05))))+2.);
                        return exp(-dt/tau_FCaSS);
                }},
};

enum {
    LUT_Ca_INDEX_fCass_RL_A,
    LUT_Ca_INDEX_fCass_RL_B,
};

#define USE_LUT

// Compute a forward step using the explicit Euler algorithm to the TP06 ODE
template <class LUT_type>
void step_FE(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
             const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
             const long num_cells, const long padded_num_cells, LUT_type &lut_V, LUT_type &lut_Ca)
{
    // Assign parameters
    const cellmodel_float_t celltype = parameters[PARAM_celltype];
    const cellmodel_float_t P_kna = parameters[PARAM_P_kna];
    const cellmodel_float_t g_K1 = parameters[PARAM_g_K1];
    const cellmodel_float_t g_Kr = parameters[PARAM_g_Kr];
    const cellmodel_float_t g_Ks = parameters[PARAM_g_Ks];
    const cellmodel_float_t g_Na = parameters[PARAM_g_Na];
    const cellmodel_float_t g_bna = parameters[PARAM_g_bna];
    const cellmodel_float_t g_CaL = parameters[PARAM_g_CaL];
    const cellmodel_float_t i_CaL_lim_delta = parameters[PARAM_i_CaL_lim_delta];
    const cellmodel_float_t g_bca = parameters[PARAM_g_bca];
    const cellmodel_float_t g_to = parameters[PARAM_g_to];
    const cellmodel_float_t K_mNa = parameters[PARAM_K_mNa];
    const cellmodel_float_t K_mk = parameters[PARAM_K_mk];
    const cellmodel_float_t P_NaK = parameters[PARAM_P_NaK];
    const cellmodel_float_t K_NaCa = parameters[PARAM_K_NaCa];
    const cellmodel_float_t K_sat = parameters[PARAM_K_sat];
    const cellmodel_float_t Km_Ca = parameters[PARAM_Km_Ca];
    const cellmodel_float_t Km_Nai = parameters[PARAM_Km_Nai];
    const cellmodel_float_t alpha = parameters[PARAM_alpha];
    const cellmodel_float_t gamma = parameters[PARAM_gamma];
    const cellmodel_float_t K_pCa = parameters[PARAM_K_pCa];
    const cellmodel_float_t g_pCa = parameters[PARAM_g_pCa];
    const cellmodel_float_t g_pK = parameters[PARAM_g_pK];
    const cellmodel_float_t Buf_c = parameters[PARAM_Buf_c];
    const cellmodel_float_t Buf_sr = parameters[PARAM_Buf_sr];
    const cellmodel_float_t Buf_ss = parameters[PARAM_Buf_ss];
    const cellmodel_float_t Ca_o = parameters[PARAM_Ca_o];
    const cellmodel_float_t EC = parameters[PARAM_EC];
    const cellmodel_float_t K_buf_c = parameters[PARAM_K_buf_c];
    const cellmodel_float_t K_buf_sr = parameters[PARAM_K_buf_sr];
    const cellmodel_float_t K_buf_ss = parameters[PARAM_K_buf_ss];
    const cellmodel_float_t K_up = parameters[PARAM_K_up];
    const cellmodel_float_t V_leak = parameters[PARAM_V_leak];
    const cellmodel_float_t V_rel = parameters[PARAM_V_rel];
    const cellmodel_float_t V_sr = parameters[PARAM_V_sr];
    const cellmodel_float_t V_ss = parameters[PARAM_V_ss];
    const cellmodel_float_t V_xfer = parameters[PARAM_V_xfer];
    const cellmodel_float_t Vmax_up = parameters[PARAM_Vmax_up];
    const cellmodel_float_t k1_prime = parameters[PARAM_k1_prime];
    const cellmodel_float_t k2_prime = parameters[PARAM_k2_prime];
    const cellmodel_float_t k3 = parameters[PARAM_k3];
    const cellmodel_float_t k4 = parameters[PARAM_k4];
    const cellmodel_float_t max_sr = parameters[PARAM_max_sr];
    const cellmodel_float_t min_sr = parameters[PARAM_min_sr];
    const cellmodel_float_t Na_o = parameters[PARAM_Na_o];
    const cellmodel_float_t Cm = parameters[PARAM_Cm];
    const cellmodel_float_t F = parameters[PARAM_F];
    const cellmodel_float_t R = parameters[PARAM_R];
    const cellmodel_float_t T = parameters[PARAM_T];
    const cellmodel_float_t V_c = parameters[PARAM_V_c];
    const cellmodel_float_t is_stimulated = parameters[PARAM_is_stimulated];
    const cellmodel_float_t stim_amplitude = parameters[PARAM_stim_amplitude];
    const cellmodel_float_t K_o = parameters[PARAM_K_o];

    #pragma omp parallel for
    for (long i = 0; i < num_cells; i++) {
        // Assign states
        const cellmodel_float_t Xr1 = states[STATE_Xr1 * padded_num_cells + i];
        const cellmodel_float_t Xr2 = states[STATE_Xr2 * padded_num_cells + i];
        const cellmodel_float_t Xs = states[STATE_Xs * padded_num_cells + i];
        const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
        const cellmodel_float_t h = states[STATE_h * padded_num_cells + i];
        const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
        const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
        const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
        const cellmodel_float_t f2 = states[STATE_f2 * padded_num_cells + i];
        const cellmodel_float_t fCass = states[STATE_fCass * padded_num_cells + i];
        const cellmodel_float_t s = states[STATE_s * padded_num_cells + i];
        const cellmodel_float_t r = states[STATE_r * padded_num_cells + i];
        const cellmodel_float_t Ca_i = states[STATE_Ca_i * padded_num_cells + i];
        const cellmodel_float_t R_prime = states[STATE_R_prime * padded_num_cells + i];
        const cellmodel_float_t Ca_SR = states[STATE_Ca_SR * padded_num_cells + i];
        const cellmodel_float_t Ca_ss = states[STATE_Ca_ss * padded_num_cells + i];
        const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];
        const cellmodel_float_t V = states[STATE_V * padded_num_cells + i];
        const cellmodel_float_t K_i = states[STATE_K_i * padded_num_cells + i];

        const auto lut_V_state = lut_V.compute_input_state(V);

        // Expressions for the Reversal potentials component
        const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
        const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
        const cellmodel_float_t E_Ks = R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
        const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

        // Expressions for the Inward rectifier potassium current component
        const cellmodel_float_t alpha_K1 =
                FP_LITERAL(0.1)
                / (FP_LITERAL(1.)
                   + FP_LITERAL(6.14421235332821e-6)
                             * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
        const cellmodel_float_t beta_K1 =
                (FP_LITERAL(0.367879441171442) * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
                 + FP_LITERAL(3.06060402008027)
                           * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K))
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V));
        const cellmodel_float_t xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
        const cellmodel_float_t i_K1 =
                FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V) * xK1_inf;

        // Expressions for the Rapid time dependent potassium current component
        const cellmodel_float_t i_Kr =
                FP_LITERAL(0.430331482911935) * g_Kr * sqrt(K_o) * (-E_K + V) * Xr1 * Xr2;

        // Expressions for the Xr1 gate component
        const cellmodel_float_t xr1_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
        const cellmodel_float_t alpha_xr1 =
                FP_LITERAL(450.)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
        const cellmodel_float_t beta_xr1 =
                FP_LITERAL(6.)
                / (FP_LITERAL(1.)
                   + Exp(FP_LITERAL(60.) / FP_LITERAL(23.) + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
        const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
        const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
        states[STATE_Xr1 * padded_num_cells + i] = dt * dXr1_dt + Xr1;

        // Expressions for the Xr2 gate component
        const cellmodel_float_t xr2_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
        const cellmodel_float_t alpha_xr2 =
                FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t beta_xr2 =
                FP_LITERAL(1.12) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
        const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
        states[STATE_Xr2 * padded_num_cells + i] = dt * dXr2_dt + Xr2;

        // Expressions for the Slow time dependent potassium current component
        const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

        // Expressions for the Xs gate component
        const cellmodel_float_t xs_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
        const cellmodel_float_t alpha_xs =
                FP_LITERAL(1400.)
                / sqrt(FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t beta_xs =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
        const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
        const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
        states[STATE_Xs * padded_num_cells + i] = dt * dXs_dt + Xs;

        // Expressions for the Fast sodium current component
        const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

        // Expressions for the m gate component
        const cellmodel_float_t m_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                  - FP_LITERAL(100.) * V / FP_LITERAL(903.)))
                                           * (1.
                                              + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                    - FP_LITERAL(100.) * V / FP_LITERAL(903.))));
        const cellmodel_float_t alpha_m =
                FP_LITERAL(1.0) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-12.) - V / FP_LITERAL(5.)));
        const cellmodel_float_t beta_m =
                FP_LITERAL(0.1) / (FP_LITERAL(1.) + Exp(FP_LITERAL(7.) + V / FP_LITERAL(5.)))
                + FP_LITERAL(0.1)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-1.) / FP_LITERAL(4.) + V / FP_LITERAL(200.)));
        const cellmodel_float_t tau_m = alpha_m * beta_m;
        const cellmodel_float_t dm_dt = (-m + m_inf) / tau_m;
        states[STATE_m * padded_num_cells + i] = dt * dm_dt + m;

        // Expressions for the h gate component
        const cellmodel_float_t h_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                  + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                                           * (1.
                                              + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                    + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
        const cellmodel_float_t alpha_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(4.43126792958051e-7) * Exp(FP_LITERAL(-0.147058823529412) * V)
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(310000.) * Exp(FP_LITERAL(0.3485) * V)
                                   + FP_LITERAL(2.7) * Exp(FP_LITERAL(0.079) * V)
                         : FP_LITERAL(0.77)
                                   / (0.13
                                      + FP_LITERAL(0.0497581410839387)
                                                * Exp(FP_LITERAL(-0.0900900900900901) * V)));
        const cellmodel_float_t tau_h = FP_LITERAL(1.0) / (alpha_h + beta_h);
        const cellmodel_float_t dh_dt = (-h + h_inf) / tau_h;
        states[STATE_h * padded_num_cells + i] = dt * dh_dt + h;

        // Expressions for the j gate component
        const cellmodel_float_t j_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                  + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                                           * (1.
                                              + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                    + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
        const cellmodel_float_t alpha_j =
                (V < FP_LITERAL(-40.)
                         ? (FP_LITERAL(37.78) + V)
                                   * (FP_LITERAL(-25428.) * Exp(FP_LITERAL(0.2444) * V)
                                      - FP_LITERAL(6.948e-6) * Exp(FP_LITERAL(-0.04391) * V))
                                   / (FP_LITERAL(1.)
                                      + FP_LITERAL(50262745825.954) * Exp(FP_LITERAL(0.311) * V))
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_j =
                (V < FP_LITERAL(-40.) ? FP_LITERAL(0.02424) * Exp(FP_LITERAL(-0.01052) * V)
                                                / (FP_LITERAL(1.)
                                                   + FP_LITERAL(0.00396086833990426)
                                                             * Exp(FP_LITERAL(-0.1378) * V))
                                      : FP_LITERAL(0.6) * Exp(FP_LITERAL(0.057) * V)
                                                / (FP_LITERAL(1.)
                                                   + Exp(FP_LITERAL(-16.) / FP_LITERAL(5.)
                                                         - V / FP_LITERAL(10.))));
        const cellmodel_float_t tau_j = FP_LITERAL(1.0) / (alpha_j + beta_j);
        const cellmodel_float_t dj_dt = (-j + j_inf) / tau_j;
        states[STATE_j * padded_num_cells + i] = dt * dj_dt + j;

        // Expressions for the Sodium background current component
        const cellmodel_float_t i_b_Na = g_bna * (-E_Na + V);

        // Expressions for the L_type Ca current component
#if 1 && defined(USE_LUT)
        const cellmodel_float_t i_CaL_fraction = lut_V.lookup(LUT_INDEX_i_CaL_fraction, lut_V_state);
        const cellmodel_float_t i_CaL_factors_e1 = lut_V.lookup(LUT_INDEX_i_CaL_factors_e1, lut_V_state);
        const cellmodel_float_t i_CaL_factors = 4. * F * g_CaL
                                 * (-Ca_o + Ca_ss * i_CaL_factors_e1) * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#else
        const cellmodel_float_t V_eff = FP_LITERAL(-15.) + V;
        const cellmodel_float_t i_CaL_factors =
                FP_LITERAL(4.) * F * g_CaL
                * (-Ca_o + Ca_ss * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) / FP_LITERAL(4.))
                * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL_fraction =
                    (fabs(V_eff) < i_CaL_lim_delta
                             ? FP_LITERAL(0.5)
                             : F * V_eff / (R * T * (Expm1(FP_LITERAL(2.) * F * V_eff / (R * T)))));
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#endif

        // Expressions for the d gate component
        const cellmodel_float_t d_inf = FP_LITERAL(1.0)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(-16.) / FP_LITERAL(15.)
                                                 - FP_LITERAL(2.) * V / FP_LITERAL(15.)));
        const cellmodel_float_t alpha_d =
                FP_LITERAL(0.25)
                + FP_LITERAL(1.4)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-35.) / FP_LITERAL(13.) - V / FP_LITERAL(13.)));
        const cellmodel_float_t beta_d =
                FP_LITERAL(1.4) / (FP_LITERAL(1.) + Exp(FP_LITERAL(1.) + V / FP_LITERAL(5.)));
        const cellmodel_float_t gamma_d =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_d = alpha_d * beta_d + gamma_d;
        const cellmodel_float_t dd_dt = (-d + d_inf) / tau_d;
        states[STATE_d * padded_num_cells + i] = dt * dd_dt + d;

        // Expressions for the f gate component
        const cellmodel_float_t f_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(20.) / FP_LITERAL(7.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f =
                FP_LITERAL(20.)
                + FP_LITERAL(180.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(200.)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(13.) / FP_LITERAL(10.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(1102.5)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(225.));
        const cellmodel_float_t df_dt = (-f + f_inf) / tau_f;
        states[STATE_f * padded_num_cells + i] = dt * df_dt + f;

        // Expressions for the F2 gate component
        const cellmodel_float_t f2_inf =
                FP_LITERAL(0.33)
                + FP_LITERAL(0.67) / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f2 =
                FP_LITERAL(31.)
                        / (FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(562.)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(240.));
        const cellmodel_float_t df2_dt = (-f2 + f2_inf) / tau_f2;
        states[STATE_f2 * padded_num_cells + i] = dt * df2_dt + f2;

        // Expressions for the FCass gate component
        const cellmodel_float_t fCass_inf =
                FP_LITERAL(0.4)
                + FP_LITERAL(0.6) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t tau_fCass =
                FP_LITERAL(2.)
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t dfCass_dt = (-fCass + fCass_inf) / tau_fCass;
        states[STATE_fCass * padded_num_cells + i] = dt * dfCass_dt + fCass;

        // Expressions for the Calcium background current component
        const cellmodel_float_t i_b_Ca = g_bca * (-E_Ca + V);

        // Expressions for the Transient outward current component
        const cellmodel_float_t i_to = g_to * (-E_K + V) * r * s;

        // Expressions for the s gate component
        const cellmodel_float_t s_inf =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.)
                                      + Exp(FP_LITERAL(28.) / FP_LITERAL(5.) + V / FP_LITERAL(5.)))
                         : FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.) + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
        const cellmodel_float_t tau_s =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(8.)
                                   + FP_LITERAL(1000.)
                                             * Exp(-((FP_LITERAL(67.) + V) * (FP_LITERAL(67.) + V))
                                                   / FP_LITERAL(1000.))
                         : FP_LITERAL(3.)
                                   + FP_LITERAL(5.)
                                             / (FP_LITERAL(1.)
                                                + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                   + FP_LITERAL(85.)
                                             * Exp(-((FP_LITERAL(45.) + V) * (FP_LITERAL(45.) + V))
                                                   / FP_LITERAL(320.)));
        const cellmodel_float_t ds_dt = (-s + s_inf) / tau_s;
        states[STATE_s * padded_num_cells + i] = dt * ds_dt + s;

        // Expressions for the r gate component
        const cellmodel_float_t r_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(10.) / FP_LITERAL(3.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t tau_r =
                FP_LITERAL(0.8)
                + FP_LITERAL(9.5)
                          * Exp(-((FP_LITERAL(40.) + V) * (FP_LITERAL(40.) + V))
                                / FP_LITERAL(1800.));
        const cellmodel_float_t dr_dt = (-r + r_inf) / tau_r;
        states[STATE_r * padded_num_cells + i] = dt * dr_dt + r;

        // Expressions for the Sodium potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_iNaK = lut_V.lookup(LUT_INDEX_rec_iNaK, lut_V_state);
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i / ((K_mNa + Na_i) * (K_mk + K_o)) * rec_iNaK;
#else
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i
                / ((K_mNa + Na_i) * (K_mk + K_o)
                   * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                      + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));
#endif

        // Expressions for the Sodium calcium exchanger current component
#ifdef USE_LUT
        const cellmodel_float_t i_NaCa_A = lut_V.lookup(LUT_INDEX_I_NaCa_A, lut_V_state);
        const cellmodel_float_t i_NaCa_B = lut_V.lookup(LUT_INDEX_I_NaCa_B, lut_V_state);
        const cellmodel_float_t i_NaCa = (i_NaCa_A * (Na_i * Na_i * Na_i)) - (i_NaCa_B * Ca_i);
#else
        const cellmodel_float_t i_NaCa =
                K_NaCa
                * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                   - alpha * (Na_o * Na_o * Na_o) * Ca_i
                             * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                   * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
#endif

        // Expressions for the Calcium pump current component
        const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

        // Expressions for the Potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_ipK = lut_V.lookup(LUT_INDEX_rec_ipK, lut_V_state);
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V) * rec_ipK;
#else
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                 - FP_LITERAL(50.) * V / FP_LITERAL(299.)));
#endif

        // Expressions for the Calcium dynamics component
        const cellmodel_float_t i_up = Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
        const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
        const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
        const cellmodel_float_t kcasr =
                max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
        const cellmodel_float_t Ca_i_bufc =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
        const cellmodel_float_t Ca_sr_bufsr =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
        const cellmodel_float_t Ca_ss_bufss =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
        const cellmodel_float_t dCa_i_dt =
                (V_sr * (-i_up + i_leak) / V_c
                 - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca) / (FP_LITERAL(2.) * F * V_c)
                 + i_xfer)
                * Ca_i_bufc;
        states[STATE_Ca_i * padded_num_cells + i] = dt * dCa_i_dt + Ca_i;
        const cellmodel_float_t k1 = k1_prime / kcasr;
        const cellmodel_float_t k2 = k2_prime * kcasr;
        const cellmodel_float_t O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
        const cellmodel_float_t dR_prime_dt =
                k4 * (FP_LITERAL(1.) - R_prime) - Ca_ss * R_prime * k2;
        states[STATE_R_prime * padded_num_cells + i] = dt * dR_prime_dt + R_prime;
        const cellmodel_float_t i_rel = V_rel * (-Ca_ss + Ca_SR) * O;
        const cellmodel_float_t dCa_SR_dt = (-i_leak - i_rel + i_up) * Ca_sr_bufsr;
        states[STATE_Ca_SR * padded_num_cells + i] = dt * dCa_SR_dt + Ca_SR;
        const cellmodel_float_t dCa_ss_dt = (V_sr * i_rel / V_ss - V_c * i_xfer / V_ss
                                             - Cm * i_CaL / (FP_LITERAL(2.) * F * V_ss))
                                            * Ca_ss_bufss;
        states[STATE_Ca_ss * padded_num_cells + i] = dt * dCa_ss_dt + Ca_ss;

        // Expressions for the Sodium dynamics component
        const cellmodel_float_t dNa_i_dt =
                Cm * (-i_Na - i_b_Na - FP_LITERAL(3.) * i_NaCa - FP_LITERAL(3.) * i_NaK)
                / (F * V_c);
        states[STATE_Na_i * padded_num_cells + i] = dt * dNa_i_dt + Na_i;

        // Expressions for the Membrane component
        const cellmodel_float_t i_Stim = (is_stimulated ? -stim_amplitude : FP_LITERAL(0.));
        const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK - i_Stim
                                        - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
        states[STATE_V * padded_num_cells + i] = dt * dV_dt + V;

        // Expressions for the Potassium dynamics component
        const cellmodel_float_t dK_i_dt =
                Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + FP_LITERAL(2.) * i_NaK)
                / (F * V_c);
        states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;
    }
}

// Compute a forward step using the GRL1 scheme to the TP06 ODE
template <class LUT_type>
void step_GRL1(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
               const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
               const long num_cells, const long padded_num_cells, LUT_type &lut_V, LUT_type &lut_Ca)
{
    // Assign parameters
    const cellmodel_float_t celltype = parameters[PARAM_celltype];
    const cellmodel_float_t P_kna = parameters[PARAM_P_kna];
    const cellmodel_float_t g_K1 = parameters[PARAM_g_K1];
    const cellmodel_float_t g_Kr = parameters[PARAM_g_Kr];
    const cellmodel_float_t g_Ks = parameters[PARAM_g_Ks];
    const cellmodel_float_t g_Na = parameters[PARAM_g_Na];
    const cellmodel_float_t g_bna = parameters[PARAM_g_bna];
    const cellmodel_float_t g_CaL = parameters[PARAM_g_CaL];
    const cellmodel_float_t i_CaL_lim_delta = parameters[PARAM_i_CaL_lim_delta];
    const cellmodel_float_t g_bca = parameters[PARAM_g_bca];
    const cellmodel_float_t g_to = parameters[PARAM_g_to];
    const cellmodel_float_t K_mNa = parameters[PARAM_K_mNa];
    const cellmodel_float_t K_mk = parameters[PARAM_K_mk];
    const cellmodel_float_t P_NaK = parameters[PARAM_P_NaK];
    const cellmodel_float_t K_NaCa = parameters[PARAM_K_NaCa];
    const cellmodel_float_t K_sat = parameters[PARAM_K_sat];
    const cellmodel_float_t Km_Ca = parameters[PARAM_Km_Ca];
    const cellmodel_float_t Km_Nai = parameters[PARAM_Km_Nai];
    const cellmodel_float_t alpha = parameters[PARAM_alpha];
    const cellmodel_float_t gamma = parameters[PARAM_gamma];
    const cellmodel_float_t K_pCa = parameters[PARAM_K_pCa];
    const cellmodel_float_t g_pCa = parameters[PARAM_g_pCa];
    const cellmodel_float_t g_pK = parameters[PARAM_g_pK];
    const cellmodel_float_t Buf_c = parameters[PARAM_Buf_c];
    const cellmodel_float_t Buf_sr = parameters[PARAM_Buf_sr];
    const cellmodel_float_t Buf_ss = parameters[PARAM_Buf_ss];
    const cellmodel_float_t Ca_o = parameters[PARAM_Ca_o];
    const cellmodel_float_t EC = parameters[PARAM_EC];
    const cellmodel_float_t K_buf_c = parameters[PARAM_K_buf_c];
    const cellmodel_float_t K_buf_sr = parameters[PARAM_K_buf_sr];
    const cellmodel_float_t K_buf_ss = parameters[PARAM_K_buf_ss];
    const cellmodel_float_t K_up = parameters[PARAM_K_up];
    const cellmodel_float_t V_leak = parameters[PARAM_V_leak];
    const cellmodel_float_t V_rel = parameters[PARAM_V_rel];
    const cellmodel_float_t V_sr = parameters[PARAM_V_sr];
    const cellmodel_float_t V_ss = parameters[PARAM_V_ss];
    const cellmodel_float_t V_xfer = parameters[PARAM_V_xfer];
    const cellmodel_float_t Vmax_up = parameters[PARAM_Vmax_up];
    const cellmodel_float_t k1_prime = parameters[PARAM_k1_prime];
    const cellmodel_float_t k2_prime = parameters[PARAM_k2_prime];
    const cellmodel_float_t k3 = parameters[PARAM_k3];
    const cellmodel_float_t k4 = parameters[PARAM_k4];
    const cellmodel_float_t max_sr = parameters[PARAM_max_sr];
    const cellmodel_float_t min_sr = parameters[PARAM_min_sr];
    const cellmodel_float_t Na_o = parameters[PARAM_Na_o];
    const cellmodel_float_t Cm = parameters[PARAM_Cm];
    const cellmodel_float_t F = parameters[PARAM_F];
    const cellmodel_float_t R = parameters[PARAM_R];
    const cellmodel_float_t T = parameters[PARAM_T];
    const cellmodel_float_t V_c = parameters[PARAM_V_c];
    const cellmodel_float_t is_stimulated = parameters[PARAM_is_stimulated];
    const cellmodel_float_t stim_amplitude = parameters[PARAM_stim_amplitude];
    const cellmodel_float_t K_o = parameters[PARAM_K_o];

#ifdef HINT_OMP_SIMD
#ifdef VECTOR_LENGTH
#pragma omp parallel for simd simdlen(VECTOR_LENGTH)
#else
#pragma omp parallel for simd
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp parallel for
#endif
    for (long i = 0; i < num_cells; i++) {
        // Assign states
        const cellmodel_float_t Xr1 = states[STATE_Xr1 * padded_num_cells + i];
        const cellmodel_float_t Xr2 = states[STATE_Xr2 * padded_num_cells + i];
        const cellmodel_float_t Xs = states[STATE_Xs * padded_num_cells + i];
        const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
        const cellmodel_float_t h = states[STATE_h * padded_num_cells + i];
        const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
        const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
        const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
        const cellmodel_float_t f2 = states[STATE_f2 * padded_num_cells + i];
        const cellmodel_float_t fCass = states[STATE_fCass * padded_num_cells + i];
        const cellmodel_float_t s = states[STATE_s * padded_num_cells + i];
        const cellmodel_float_t r = states[STATE_r * padded_num_cells + i];
        const cellmodel_float_t Ca_i = states[STATE_Ca_i * padded_num_cells + i];
        const cellmodel_float_t R_prime = states[STATE_R_prime * padded_num_cells + i];
        const cellmodel_float_t Ca_SR = states[STATE_Ca_SR * padded_num_cells + i];
        const cellmodel_float_t Ca_ss = states[STATE_Ca_ss * padded_num_cells + i];
        const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];
        const cellmodel_float_t V = states[STATE_V * padded_num_cells + i];
        const cellmodel_float_t K_i = states[STATE_K_i * padded_num_cells + i];

        const auto lut_V_state = lut_V.compute_input_state(V);
        const auto lut_Ca_state = lut_Ca.compute_input_state(Ca_ss);

        // Expressions for the Reversal potentials component
        const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
        const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
        const cellmodel_float_t E_Ks = R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
        const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

        // Expressions for the Inward rectifier potassium current component
        const cellmodel_float_t alpha_K1 =
                FP_LITERAL(0.1)
                / (FP_LITERAL(1.)
                   + FP_LITERAL(6.14421235332821e-6)
                             * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
        const cellmodel_float_t beta_K1 =
                (FP_LITERAL(0.367879441171442) * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
                 + FP_LITERAL(3.06060402008027)
                           * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K))
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V));
        const cellmodel_float_t xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
        const cellmodel_float_t i_K1 =
                FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V) * xK1_inf;

        // Expressions for the Rapid time dependent potassium current component
        const cellmodel_float_t i_Kr =
                FP_LITERAL(0.430331482911935) * g_Kr * sqrt(K_o) * (-E_K + V) * Xr1 * Xr2;

        // Expressions for the Xr1 gate component
#ifdef USE_LUT
        const cellmodel_float_t Xr1_RL_A = lut_V.lookup(LUT_INDEX_Xr1_RL_A, lut_V_state);
        const cellmodel_float_t Xr1_RL_B = lut_V.lookup(LUT_INDEX_Xr1_RL_B, lut_V_state);
        states[STATE_Xr1 * padded_num_cells + i] = Xr1_RL_A + Xr1_RL_B * Xr1;
#else
        const cellmodel_float_t xr1_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
        const cellmodel_float_t alpha_xr1 =
                FP_LITERAL(450.)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
        const cellmodel_float_t beta_xr1 =
                FP_LITERAL(6.)
                / (FP_LITERAL(1.)
                   + Exp(FP_LITERAL(60.) / FP_LITERAL(23.) + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
        const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
        const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
        const cellmodel_float_t dXr1_dt_linearized = FP_LITERAL(-1.) / tau_xr1;
        states[STATE_Xr1 * padded_num_cells + i] =
                (fabs(dXr1_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXr1_dt_linearized)) * dXr1_dt
                                   / dXr1_dt_linearized
                         : dt * dXr1_dt)
                + Xr1;
#endif

        // Expressions for the Xr2 gate component
#ifdef USE_LUT
        const cellmodel_float_t Xr2_RL_A = lut_V.lookup(LUT_INDEX_Xr2_RL_A, lut_V_state);
        const cellmodel_float_t Xr2_RL_B = lut_V.lookup(LUT_INDEX_Xr2_RL_B, lut_V_state);
        states[STATE_Xr2 * padded_num_cells + i] = Xr2_RL_A + Xr2_RL_B * Xr2;
#else
        const cellmodel_float_t xr2_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
        const cellmodel_float_t alpha_xr2 =
                FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t beta_xr2 =
                FP_LITERAL(1.12) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
        const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
        const cellmodel_float_t dXr2_dt_linearized = FP_LITERAL(-1.) / tau_xr2;
        states[STATE_Xr2 * padded_num_cells + i] =
                (fabs(dXr2_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXr2_dt_linearized)) * dXr2_dt
                                   / dXr2_dt_linearized
                         : dt * dXr2_dt)
                + Xr2;
#endif

        // Expressions for the Slow time dependent potassium current component
        const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

        // Expressions for the Xs gate component
#ifdef USE_LUT
        const cellmodel_float_t Xs_RL_A = lut_V.lookup(LUT_INDEX_Xs_RL_A, lut_V_state);
        const cellmodel_float_t Xs_RL_B = lut_V.lookup(LUT_INDEX_Xs_RL_B, lut_V_state);
        states[STATE_Xs * padded_num_cells + i] = Xs_RL_A + Xs_RL_B * Xs;
#else
        const cellmodel_float_t xs_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
        const cellmodel_float_t alpha_xs =
                FP_LITERAL(1400.)
                / sqrt(FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t beta_xs =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
        const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
        const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
        const cellmodel_float_t dXs_dt_linearized = FP_LITERAL(-1.) / tau_xs;
        states[STATE_Xs * padded_num_cells + i] =
                (fabs(dXs_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXs_dt_linearized)) * dXs_dt
                                   / dXs_dt_linearized
                         : dt * dXs_dt)
                + Xs;
#endif

        // Expressions for the Fast sodium current component
        const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

        // Expressions for the m gate component
#ifdef USE_LUT
        const cellmodel_float_t m_RL_A = lut_V.lookup(LUT_INDEX_m_RL_A, lut_V_state);
        const cellmodel_float_t m_RL_B = lut_V.lookup(LUT_INDEX_m_RL_B, lut_V_state);
        states[STATE_m * padded_num_cells + i] = m_RL_A + m_RL_B * m;
#else
        const cellmodel_float_t m_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                  - FP_LITERAL(100.) * V / FP_LITERAL(903.)))
                                           * (1.
                                              + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                    - FP_LITERAL(100.) * V / FP_LITERAL(903.))));
        const cellmodel_float_t alpha_m =
                FP_LITERAL(1.0) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-12.) - V / FP_LITERAL(5.)));
        const cellmodel_float_t beta_m =
                FP_LITERAL(0.1) / (FP_LITERAL(1.) + Exp(FP_LITERAL(7.) + V / FP_LITERAL(5.)))
                + FP_LITERAL(0.1)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-1.) / FP_LITERAL(4.) + V / FP_LITERAL(200.)));
        const cellmodel_float_t tau_m = alpha_m * beta_m;
        const cellmodel_float_t dm_dt = (-m + m_inf) / tau_m;
        const cellmodel_float_t dm_dt_linearized = FP_LITERAL(-1.) / tau_m;
        states[STATE_m * padded_num_cells + i] =
                (fabs(dm_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dm_dt_linearized)) * dm_dt
                                   / dm_dt_linearized
                         : dt * dm_dt)
                + m;
#endif

        // Expressions for the h gate component
#ifdef USE_LUT
        const cellmodel_float_t h_RL_A = lut_V.lookup(LUT_INDEX_h_RL_A, lut_V_state);
        const cellmodel_float_t h_RL_B = lut_V.lookup(LUT_INDEX_h_RL_B, lut_V_state);
        states[STATE_h * padded_num_cells + i] = h_RL_A + h_RL_B * h;
#else
        const cellmodel_float_t h_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                  + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                                           * (1.
                                              + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                    + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
        const cellmodel_float_t alpha_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(4.43126792958051e-7) * Exp(FP_LITERAL(-0.147058823529412) * V)
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(310000.) * Exp(FP_LITERAL(0.3485) * V)
                                   + FP_LITERAL(2.7) * Exp(FP_LITERAL(0.079) * V)
                         : FP_LITERAL(0.77)
                                   / (0.13
                                      + FP_LITERAL(0.0497581410839387)
                                                * Exp(FP_LITERAL(-0.0900900900900901) * V)));
        const cellmodel_float_t tau_h = FP_LITERAL(1.0) / (alpha_h + beta_h);
        const cellmodel_float_t dh_dt = (-h + h_inf) / tau_h;
        const cellmodel_float_t dh_dt_linearized = FP_LITERAL(-1.) / tau_h;
        states[STATE_h * padded_num_cells + i] =
                (fabs(dh_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dh_dt_linearized)) * dh_dt
                                   / dh_dt_linearized
                         : dt * dh_dt)
                + h;
#endif

        // Expressions for the j gate component
#ifdef USE_LUT
        const cellmodel_float_t j_RL_A = lut_V.lookup(LUT_INDEX_j_RL_A, lut_V_state);
        const cellmodel_float_t j_RL_B = lut_V.lookup(LUT_INDEX_j_RL_B, lut_V_state);
        states[STATE_j * padded_num_cells + i] = j_RL_A + j_RL_B * j;
#else
        const cellmodel_float_t j_inf =
                FP_LITERAL(1.0)
                / ((FP_LITERAL(1.)
                    + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V))
                   * (FP_LITERAL(1.)
                      + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V)));
        const cellmodel_float_t alpha_j =
                (V < FP_LITERAL(-40.)
                         ? (FP_LITERAL(37.78) + V)
                                   * (FP_LITERAL(-25428.) * Exp(FP_LITERAL(0.2444) * V)
                                      - FP_LITERAL(6.948e-6) * Exp(FP_LITERAL(-0.04391) * V))
                                   / (FP_LITERAL(1.)
                                      + FP_LITERAL(50262745825.954) * Exp(FP_LITERAL(0.311) * V))
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_j =
                (V < FP_LITERAL(-40.) ? FP_LITERAL(0.02424) * Exp(FP_LITERAL(-0.01052) * V)
                                                / (FP_LITERAL(1.)
                                                   + FP_LITERAL(0.00396086833990426)
                                                             * Exp(FP_LITERAL(-0.1378) * V))
                                      : FP_LITERAL(0.6) * Exp(FP_LITERAL(0.057) * V)
                                                / (FP_LITERAL(1.)
                                                   + Exp(FP_LITERAL(-16.) / FP_LITERAL(5.)
                                                         - V / FP_LITERAL(10.))));
        const cellmodel_float_t tau_j = FP_LITERAL(1.0) / (alpha_j + beta_j);
        const cellmodel_float_t dj_dt = (-j + j_inf) / tau_j;
        const cellmodel_float_t dj_dt_linearized = FP_LITERAL(-1.) / tau_j;
        states[STATE_j * padded_num_cells + i] =
                (fabs(dj_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dj_dt_linearized)) * dj_dt
                                   / dj_dt_linearized
                         : dt * dj_dt)
                + j;
#endif

        // Expressions for the Sodium background current component
        const cellmodel_float_t i_b_Na = g_bna * (-E_Na + V);

        // Expressions for the L_type Ca current component
        const cellmodel_float_t V_eff = FP_LITERAL(-15.) + V;
#if 1 && defined(USE_LUT)
        const cellmodel_float_t i_CaL_fraction = lut_V.lookup(LUT_INDEX_i_CaL_fraction, lut_V_state);
        const cellmodel_float_t i_CaL_factors_e1 = lut_V.lookup(LUT_INDEX_i_CaL_factors_e1, lut_V_state);
        const cellmodel_float_t i_CaL_factors = 4. * F * g_CaL
                                 * (-Ca_o + Ca_ss * i_CaL_factors_e1) * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#else
        const cellmodel_float_t i_CaL_factors =
                FP_LITERAL(4.) * F * g_CaL
                * (-Ca_o + Ca_ss * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) / FP_LITERAL(4.))
                * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL_fraction =
                    (fabs(V_eff) < i_CaL_lim_delta
                             ? FP_LITERAL(0.5)
                             : F * V_eff / (R * T * (Expm1(FP_LITERAL(2.) * F * V_eff / (R * T)))));
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#endif

        // Expressions for the d gate component
#ifdef USE_LUT
        const cellmodel_float_t d_RL_A = lut_V.lookup(LUT_INDEX_d_RL_A, lut_V_state);
        const cellmodel_float_t d_RL_B = lut_V.lookup(LUT_INDEX_d_RL_B, lut_V_state);
        states[STATE_d * padded_num_cells + i] = d_RL_A + d_RL_B * d;
#else
        const cellmodel_float_t d_inf = FP_LITERAL(1.0)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(-16.) / FP_LITERAL(15.)
                                                 - FP_LITERAL(2.) * V / FP_LITERAL(15.)));
        const cellmodel_float_t alpha_d =
                FP_LITERAL(0.25)
                + FP_LITERAL(1.4)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-35.) / FP_LITERAL(13.) - V / FP_LITERAL(13.)));
        const cellmodel_float_t beta_d =
                FP_LITERAL(1.4) / (FP_LITERAL(1.) + Exp(FP_LITERAL(1.) + V / FP_LITERAL(5.)));
        const cellmodel_float_t gamma_d =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_d = alpha_d * beta_d + gamma_d;
        const cellmodel_float_t dd_dt = (-d + d_inf) / tau_d;
        const cellmodel_float_t dd_dt_linearized = FP_LITERAL(-1.) / tau_d;
        states[STATE_d * padded_num_cells + i] =
                (fabs(dd_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dd_dt_linearized)) * dd_dt
                                   / dd_dt_linearized
                         : dt * dd_dt)
                + d;
#endif

        // Expressions for the f gate component
#ifdef USE_LUT
        const cellmodel_float_t f_RL_A = lut_V.lookup(LUT_INDEX_f_RL_A, lut_V_state);
        const cellmodel_float_t f_RL_B = lut_V.lookup(LUT_INDEX_f_RL_B, lut_V_state);
        states[STATE_f * padded_num_cells + i] = f_RL_A + f_RL_B * f;
#else
        const cellmodel_float_t f_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(20.) / FP_LITERAL(7.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f =
                FP_LITERAL(20.)
                + FP_LITERAL(180.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(200.)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(13.) / FP_LITERAL(10.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(1102.5)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(225.));
        const cellmodel_float_t df_dt = (-f + f_inf) / tau_f;
        const cellmodel_float_t df_dt_linearized = FP_LITERAL(-1.) / tau_f;
        states[STATE_f * padded_num_cells + i] =
                (fabs(df_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * df_dt_linearized)) * df_dt
                                   / df_dt_linearized
                         : dt * df_dt)
                + f;
#endif

        // Expressions for the F2 gate component
#ifdef USE_LUT
        const cellmodel_float_t f2_RL_A = lut_V.lookup(LUT_INDEX_f2_RL_A, lut_V_state);
        const cellmodel_float_t f2_RL_B = lut_V.lookup(LUT_INDEX_f2_RL_B, lut_V_state);
        states[STATE_f2 * padded_num_cells + i] = f2_RL_A + f2_RL_B * f2;
#else
        const cellmodel_float_t f2_inf =
                FP_LITERAL(0.33)
                + FP_LITERAL(0.67) / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f2 =
                FP_LITERAL(31.)
                        / (FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(562.)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(240.));
        const cellmodel_float_t df2_dt = (-f2 + f2_inf) / tau_f2;
        const cellmodel_float_t df2_dt_linearized = FP_LITERAL(-1.) / tau_f2;
        states[STATE_f2 * padded_num_cells + i] =
                (fabs(df2_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * df2_dt_linearized)) * df2_dt
                                   / df2_dt_linearized
                         : dt * df2_dt)
                + f2;
#endif

        // Expressions for the FCass gate component
#if 1 && defined(USE_LUT)
        const cellmodel_float_t fCass_RL_A = lut_Ca.lookup(LUT_Ca_INDEX_fCass_RL_A, lut_Ca_state);
        const cellmodel_float_t fCass_RL_B = lut_Ca.lookup(LUT_Ca_INDEX_fCass_RL_B, lut_Ca_state);
        states[STATE_fCass * padded_num_cells + i] = fCass_RL_A + fCass_RL_B * fCass;
#else
        const cellmodel_float_t fCass_inf =
                FP_LITERAL(0.4)
                + FP_LITERAL(0.6) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t tau_fCass =
                FP_LITERAL(2.)
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t dfCass_dt = (-fCass + fCass_inf) / tau_fCass;
        const cellmodel_float_t dfCass_dt_linearized = FP_LITERAL(-1.) / tau_fCass;
        states[STATE_fCass * padded_num_cells + i] =
                (fabs(dfCass_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dfCass_dt_linearized)) * dfCass_dt
                                   / dfCass_dt_linearized
                         : dt * dfCass_dt)
                + fCass;
#endif

        // Expressions for the Calcium background current component
        const cellmodel_float_t i_b_Ca = g_bca * (-E_Ca + V);

        // Expressions for the Transient outward current component
        const cellmodel_float_t i_to = g_to * (-E_K + V) * r * s;

        // Expressions for the s gate component
#ifdef USE_LUT
        const cellmodel_float_t s_RL_A = lut_V.lookup(LUT_INDEX_s_RL_A, lut_V_state);
        const cellmodel_float_t s_RL_B = lut_V.lookup(LUT_INDEX_s_RL_B, lut_V_state);
        states[STATE_s * padded_num_cells + i] = s_RL_A + s_RL_B * s;
#else
        const cellmodel_float_t s_inf =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.)
                                      + Exp(FP_LITERAL(28.) / FP_LITERAL(5.) + V / FP_LITERAL(5.)))
                         : FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.) + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
        const cellmodel_float_t tau_s =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(8.)
                                   + FP_LITERAL(1000.)
                                             * Exp(-((FP_LITERAL(67.) + V) * (FP_LITERAL(67.) + V))
                                                   / FP_LITERAL(1000.))
                         : FP_LITERAL(3.)
                                   + FP_LITERAL(5.)
                                             / (FP_LITERAL(1.)
                                                + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                   + FP_LITERAL(85.)
                                             * Exp(-((FP_LITERAL(45.) + V) * (FP_LITERAL(45.) + V))
                                                   / FP_LITERAL(320.)));
        const cellmodel_float_t ds_dt = (-s + s_inf) / tau_s;
        const cellmodel_float_t ds_dt_linearized = FP_LITERAL(-1.) / tau_s;
        states[STATE_s * padded_num_cells + i] =
                (fabs(ds_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * ds_dt_linearized)) * ds_dt
                                   / ds_dt_linearized
                         : dt * ds_dt)
                + s;
#endif

        // Expressions for the r gate component
#ifdef USE_LUT
        const cellmodel_float_t r_RL_A = lut_V.lookup(LUT_INDEX_r_RL_A, lut_V_state);
        const cellmodel_float_t r_RL_B = lut_V.lookup(LUT_INDEX_r_RL_B, lut_V_state);
        states[STATE_r * padded_num_cells + i] = r_RL_A + r_RL_B * r;
#else
        const cellmodel_float_t r_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(10.) / FP_LITERAL(3.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t tau_r =
                FP_LITERAL(0.8)
                + FP_LITERAL(9.5)
                          * Exp(-((FP_LITERAL(40.) + V) * (FP_LITERAL(40.) + V))
                                / FP_LITERAL(1800.));
        const cellmodel_float_t dr_dt = (-r + r_inf) / tau_r;
        const cellmodel_float_t dr_dt_linearized = FP_LITERAL(-1.) / tau_r;
        states[STATE_r * padded_num_cells + i] =
                (fabs(dr_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dr_dt_linearized)) * dr_dt
                                   / dr_dt_linearized
                         : dt * dr_dt)
                + r;
#endif

        // Expressions for the Sodium potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_iNaK = lut_V.lookup(LUT_INDEX_rec_iNaK, lut_V_state);
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i / ((K_mNa + Na_i) * (K_mk + K_o)) * rec_iNaK;
#else
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i
                / ((K_mNa + Na_i) * (K_mk + K_o)
                   * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                      + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));
#endif

        // Expressions for the Sodium calcium exchanger current component
#ifdef USE_LUT
        const cellmodel_float_t i_NaCa_A = lut_V.lookup(LUT_INDEX_I_NaCa_A, lut_V_state);
        const cellmodel_float_t i_NaCa_B = lut_V.lookup(LUT_INDEX_I_NaCa_B, lut_V_state);
        const cellmodel_float_t i_NaCa = (i_NaCa_A * (Na_i * Na_i * Na_i)) - (i_NaCa_B * Ca_i);
#else
        const cellmodel_float_t i_NaCa =
                K_NaCa
                * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                   - alpha * (Na_o * Na_o * Na_o) * Ca_i
                             * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                   * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
#endif

        // Expressions for the Calcium pump current component
        const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

        // Expressions for the Potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_ipK = lut_V.lookup(LUT_INDEX_rec_ipK, lut_V_state);
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V) * rec_ipK;
#else
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                 - FP_LITERAL(50.) * V / FP_LITERAL(299.)));
#endif

        // Expressions for the Calcium dynamics component
        const cellmodel_float_t i_up = Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
        const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
        const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
        const cellmodel_float_t kcasr =
                max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
        const cellmodel_float_t Ca_i_bufc =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
        const cellmodel_float_t Ca_sr_bufsr =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
        const cellmodel_float_t Ca_ss_bufss =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
        const cellmodel_float_t dCa_i_dt =
                (V_sr * (-i_up + i_leak) / V_c
                 - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca) / (FP_LITERAL(2.) * F * V_c)
                 + i_xfer)
                * Ca_i_bufc;
        const cellmodel_float_t di_p_Ca_dCa_i =
                g_pCa / (K_pCa + Ca_i) - g_pCa * Ca_i / ((K_pCa + Ca_i) * (K_pCa + Ca_i));
        const cellmodel_float_t dE_Ca_dCa_i = FP_LITERAL(-0.5) * R * T / (F * Ca_i);
        const cellmodel_float_t dCa_i_bufc_dCa_i =
                FP_LITERAL(2.) * Buf_c * K_buf_c
                / (((FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
                    * (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i))))
                   * ((K_buf_c + Ca_i) * (K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
        const cellmodel_float_t di_NaCa_dCa_i =
                -K_NaCa * alpha * (Na_o * Na_o * Na_o)
                * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T))
                / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                   * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
        const cellmodel_float_t di_up_dCa_i =
                FP_LITERAL(2.) * Vmax_up * (K_up * K_up)
                / (((FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i))
                    * (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i)))
                   * (Ca_i * Ca_i * Ca_i));
        const cellmodel_float_t dCa_i_dt_linearized =
                (-V_xfer + V_sr * (-V_leak - di_up_dCa_i) / V_c
                 - Cm * (FP_LITERAL(-2.) * di_NaCa_dCa_i - g_bca * dE_Ca_dCa_i + di_p_Ca_dCa_i)
                           / (FP_LITERAL(2.) * F * V_c))
                        * Ca_i_bufc
                + (V_sr * (-i_up + i_leak) / V_c
                   - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca) / (FP_LITERAL(2.) * F * V_c)
                   + i_xfer)
                          * dCa_i_bufc_dCa_i;
        states[STATE_Ca_i * padded_num_cells + i] =
                Ca_i
                + (fabs(dCa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                           ? (FP_LITERAL(-1.0) + Exp(dt * dCa_i_dt_linearized)) * dCa_i_dt
                                     / dCa_i_dt_linearized
                           : dt * dCa_i_dt);
        const cellmodel_float_t k1 = k1_prime / kcasr;
        const cellmodel_float_t k2 = k2_prime * kcasr;
        const cellmodel_float_t O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
        const cellmodel_float_t dR_prime_dt =
                k4 * (FP_LITERAL(1.) - R_prime) - Ca_ss * R_prime * k2;
        const cellmodel_float_t dR_prime_dt_linearized = -k4 - Ca_ss * k2;
        states[STATE_R_prime * padded_num_cells + i] =
                (fabs(dR_prime_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dR_prime_dt_linearized)) * dR_prime_dt
                                   / dR_prime_dt_linearized
                         : dt * dR_prime_dt)
                + R_prime;
        const cellmodel_float_t i_rel = V_rel * (-Ca_ss + Ca_SR) * O;
        const cellmodel_float_t dCa_SR_dt = (-i_leak - i_rel + i_up) * Ca_sr_bufsr;
        const cellmodel_float_t dO_dk1 =
                (Ca_ss * Ca_ss) * R_prime / (k3 + (Ca_ss * Ca_ss) * k1)
                - (((Ca_ss) * (Ca_ss)) * ((Ca_ss) * (Ca_ss))) * R_prime * k1
                          / ((k3 + (Ca_ss * Ca_ss) * k1) * (k3 + (Ca_ss * Ca_ss) * k1));
        const cellmodel_float_t dCa_sr_bufsr_dCa_SR =
                FP_LITERAL(2.) * Buf_sr * K_buf_sr
                / (((FP_LITERAL(1.) + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
                    * (FP_LITERAL(1.)
                       + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR))))
                   * ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
        const cellmodel_float_t dkcasr_dCa_SR =
                FP_LITERAL(-2.) * (EC * EC) * (max_sr - min_sr)
                / (((FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR))
                    * (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR)))
                   * (Ca_SR * Ca_SR * Ca_SR));
        const cellmodel_float_t dk1_dkcasr = -k1_prime / (kcasr * kcasr);
        const cellmodel_float_t di_rel_dCa_SR =
                V_rel * O + V_rel * (-Ca_ss + Ca_SR) * dO_dk1 * dk1_dkcasr * dkcasr_dCa_SR;
        const cellmodel_float_t di_rel_dO = V_rel * (-Ca_ss + Ca_SR);
        const cellmodel_float_t dCa_SR_dt_linearized =
                (-V_leak - di_rel_dCa_SR - dO_dk1 * di_rel_dO * dk1_dkcasr * dkcasr_dCa_SR)
                        * Ca_sr_bufsr
                + (-i_leak - i_rel + i_up) * dCa_sr_bufsr_dCa_SR;
        states[STATE_Ca_SR * padded_num_cells + i] =
                Ca_SR
                + (fabs(dCa_SR_dt_linearized) > FP_LITERAL(1.0e-8)
                           ? (FP_LITERAL(-1.0) + Exp(dt * dCa_SR_dt_linearized)) * dCa_SR_dt
                                     / dCa_SR_dt_linearized
                           : dt * dCa_SR_dt);
        const cellmodel_float_t dCa_ss_dt = (V_sr * i_rel / V_ss - V_c * i_xfer / V_ss
                                             - Cm * i_CaL / (FP_LITERAL(2.) * F * V_ss))
                                            * Ca_ss_bufss;

        const cellmodel_float_t dCa_ss_bufss_dCa_ss =
                FP_LITERAL(2.) * Buf_ss * K_buf_ss
                / (((FP_LITERAL(1.) + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
                    * (FP_LITERAL(1.)
                       + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss))))
                   * ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
        const cellmodel_float_t dO_dCa_ss =
                FP_LITERAL(-2.) * (Ca_ss * Ca_ss * Ca_ss) * (k1 * k1) * R_prime
                        / ((k3 + (Ca_ss * Ca_ss) * k1) * (k3 + (Ca_ss * Ca_ss) * k1))
                + FP_LITERAL(2.) * Ca_ss * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
        const cellmodel_float_t di_CaL_factors_dCa_ss =
                    F * g_CaL * d * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) * f * f2 * fCass;
        const cellmodel_float_t di_rel_dCa_ss = -V_rel * O + V_rel * (-Ca_ss + Ca_SR) * dO_dCa_ss;
        const cellmodel_float_t dCa_ss_dt_linearized =
                (V_sr * (dO_dCa_ss * di_rel_dO + di_rel_dCa_ss) / V_ss - V_c * V_xfer / V_ss
                 - Cm * di_CaL_factors_dCa_ss * i_CaL_fraction / (FP_LITERAL(2.) * F * V_ss))
                        * Ca_ss_bufss
                + (V_sr * i_rel / V_ss - V_c * i_xfer / V_ss
                   - Cm * i_CaL / (FP_LITERAL(2.) * F * V_ss))
                          * dCa_ss_bufss_dCa_ss;
        states[STATE_Ca_ss * padded_num_cells + i] =
                Ca_ss
                + (fabs(dCa_ss_dt_linearized) > FP_LITERAL(1.0e-8)
                           ? (FP_LITERAL(-1.0) + Exp(dt * dCa_ss_dt_linearized)) * dCa_ss_dt
                                     / dCa_ss_dt_linearized
                           : dt * dCa_ss_dt);

        // Expressions for the Sodium dynamics component
        const cellmodel_float_t dNa_i_dt =
                Cm * (-i_Na - i_b_Na - FP_LITERAL(3.) * i_NaCa - FP_LITERAL(3.) * i_NaK)
                / (F * V_c);
        const cellmodel_float_t di_NaCa_dNa_i =
                FP_LITERAL(3.) * Ca_o * K_NaCa * (Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                   * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
        const cellmodel_float_t di_Na_dE_Na = -g_Na * (m * m * m) * h * j;
        const cellmodel_float_t dE_Na_dNa_i = -R * T / (F * Na_i);
        const cellmodel_float_t di_NaK_dNa_i =
                K_o * P_NaK
                        / ((K_mNa + Na_i) * (K_mk + K_o)
                           * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                              + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))))
                - K_o * P_NaK * Na_i
                          / (((K_mNa + Na_i) * (K_mNa + Na_i)) * (K_mk + K_o)
                             * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                                + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));
        const cellmodel_float_t dNa_i_dt_linearized =
                Cm
                * (FP_LITERAL(-3.) * di_NaCa_dNa_i - FP_LITERAL(3.) * di_NaK_dNa_i
                   + g_bna * dE_Na_dNa_i - dE_Na_dNa_i * di_Na_dE_Na)
                / (F * V_c);
        states[STATE_Na_i * padded_num_cells + i] =
                Na_i
                + (fabs(dNa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                           ? (FP_LITERAL(-1.0) + Exp(dt * dNa_i_dt_linearized)) * dNa_i_dt
                                     / dNa_i_dt_linearized
                           : dt * dNa_i_dt);

        // Expressions for the Membrane component
        const cellmodel_float_t i_Stim = (is_stimulated ? -stim_amplitude : FP_LITERAL(0.));
        const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK - i_Stim
                                        - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
        const cellmodel_float_t dbeta_K1_dV =
                (FP_LITERAL(0.000612120804016053)
                         * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K)
                 + FP_LITERAL(0.0367879441171442)
                           * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K))
                        / (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V))
                + FP_LITERAL(0.5)
                          * (FP_LITERAL(0.367879441171442)
                                     * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
                             + FP_LITERAL(3.06060402008027)
                                       * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K))
                          * Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)
                          / ((FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V))
                             * (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)));
        const cellmodel_float_t di_CaL_factors_dV_eff =
                    FP_LITERAL(2.) * g_CaL * (F * F) * Ca_ss * d
                    * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) * f * f2 * fCass / (R * T);
        const cellmodel_float_t di_CaL_fraction_dV_eff =
                (fabs(V_eff) < i_CaL_lim_delta
                        ? 0.
                        : F / (R * T * (Expm1(FP_LITERAL(2.) * F * V_eff / (R * T))))
                                - FP_LITERAL(2.) * (F * F) * V_eff
                                                * Exp(FP_LITERAL(2.) * F * V_eff / (R * T))
                                                / ((R * R) * (T * T)
                                                * ((Expm1(FP_LITERAL(2.) * F * V_eff / (R * T)))
                                                * (Expm1(FP_LITERAL(2.) * F * V_eff
                                                        / (R * T))))));
        const cellmodel_float_t di_to_dV = g_to * r * s;
        const cellmodel_float_t di_Kr_dV =
                FP_LITERAL(0.430331482911935) * g_Kr * sqrt(K_o) * Xr1 * Xr2;
        const cellmodel_float_t di_NaK_dV =
                K_o * P_NaK
                * (FP_LITERAL(0.0353) * F * Exp(-F * V / (R * T)) / (R * T)
                   + FP_LITERAL(0.01245) * F * Exp(FP_LITERAL(-0.1) * F * V / (R * T)) / (R * T))
                * Na_i
                / ((K_mNa + Na_i) * (K_mk + K_o)
                   * ((FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                       + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T)))
                      * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                         + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T)))));
        const cellmodel_float_t di_Ks_dV = g_Ks * (Xs * Xs);
        const cellmodel_float_t di_K1_dxK1_inf =
                FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V);
        const cellmodel_float_t di_NaCa_dV =
                K_NaCa
                        * (Ca_o * F * gamma * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                                   / (R * T)
                           - F * alpha * (Na_o * Na_o * Na_o) * (FP_LITERAL(-1.) + gamma) * Ca_i
                                     * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)) / (R * T))
                        / ((FP_LITERAL(1.)
                            + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                           * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)))
                - F * K_NaCa * K_sat * (FP_LITERAL(-1.) + gamma)
                          * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                             - alpha * (Na_o * Na_o * Na_o) * Ca_i
                                       * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                          * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T))
                          / (R * T
                             * ((FP_LITERAL(1.)
                                 + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                                * (FP_LITERAL(1.)
                                   + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T))))
                             * (Ca_o + Km_Ca)
                             * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
        const cellmodel_float_t dalpha_K1_dV =
                FP_LITERAL(-3.68652741199693e-8)
                * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K)
                / ((FP_LITERAL(1.)
                    + FP_LITERAL(6.14421235332821e-6)
                              * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K))
                   * (1.
                      + FP_LITERAL(6.14421235332821e-6)
                                * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K)));
        const cellmodel_float_t dxK1_inf_dbeta_K1 =
                -alpha_K1 / ((alpha_K1 + beta_K1) * (alpha_K1 + beta_K1));
        const cellmodel_float_t dxK1_inf_dalpha_K1 =
                FP_LITERAL(1.0) / (alpha_K1 + beta_K1)
                - alpha_K1 / ((alpha_K1 + beta_K1) * (alpha_K1 + beta_K1));
        const cellmodel_float_t di_K1_dV =
                FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * xK1_inf
                + FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V)
                          * (dalpha_K1_dV * dxK1_inf_dalpha_K1 + dbeta_K1_dV * dxK1_inf_dbeta_K1);
        const cellmodel_float_t di_Na_dV = g_Na * (m * m * m) * h * j;
        const cellmodel_float_t di_p_K_dV =
                g_pK
                        / (FP_LITERAL(1.)
                           + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                 - FP_LITERAL(50.) * V / FP_LITERAL(299.)))
                + FP_LITERAL(50.) * g_pK * (-E_K + V)
                          * Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                - FP_LITERAL(50.) * V / FP_LITERAL(299.))
                          / (299.
                             * ((1.
                                 + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                       - FP_LITERAL(50.) * V / FP_LITERAL(299.)))
                                * (1.
                                   + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                         - FP_LITERAL(50.) * V / FP_LITERAL(299.)))));
        const cellmodel_float_t dV_dt_linearized =
                -g_bca - g_bna - di_K1_dV - di_Kr_dV - di_Ks_dV - di_NaCa_dV - di_NaK_dV - di_Na_dV
                - di_p_K_dV - di_to_dV
                - (dalpha_K1_dV * dxK1_inf_dalpha_K1 + dbeta_K1_dV * dxK1_inf_dbeta_K1)
                          * di_K1_dxK1_inf
                - di_CaL_factors_dV_eff * i_CaL_fraction
                    - di_CaL_fraction_dV_eff * i_CaL_factors;
        states[STATE_V * padded_num_cells + i] =
                (fabs(dV_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dV_dt_linearized)) * dV_dt
                                   / dV_dt_linearized
                         : dt * dV_dt)
                + V;

        // Expressions for the Potassium dynamics component
        const cellmodel_float_t dK_i_dt =
                Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + FP_LITERAL(2.) * i_NaK)
                / (F * V_c);
        const cellmodel_float_t di_Kr_dE_K =
                FP_LITERAL(-0.430331482911935) * g_Kr * sqrt(K_o) * Xr1 * Xr2;
        const cellmodel_float_t dE_K_dK_i = -R * T / (F * K_i);
        const cellmodel_float_t dbeta_K1_dE_K =
                (FP_LITERAL(-0.000612120804016053)
                         * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K)
                 - FP_LITERAL(0.0367879441171442)
                           * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K))
                        / (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V))
                - FP_LITERAL(0.5)
                          * (FP_LITERAL(0.367879441171442)
                                     * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
                             + FP_LITERAL(3.06060402008027)
                                       * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K))
                          * Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)
                          / ((FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V))
                             * (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)));
        const cellmodel_float_t dalpha_K1_dE_K =
                FP_LITERAL(3.68652741199693e-8) * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K)
                / ((FP_LITERAL(1.)
                    + FP_LITERAL(6.14421235332821e-6)
                              * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K))
                   * (1.
                      + FP_LITERAL(6.14421235332821e-6)
                                * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K)));
        const cellmodel_float_t di_K1_dE_K =
                FP_LITERAL(-0.430331482911935) * g_K1 * sqrt(K_o) * xK1_inf
                + FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V)
                          * (dalpha_K1_dE_K * dxK1_inf_dalpha_K1
                             + dbeta_K1_dE_K * dxK1_inf_dbeta_K1);
        const cellmodel_float_t di_p_K_dE_K = -g_pK
                                              / (FP_LITERAL(1.)
                                                 + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                       - FP_LITERAL(50.) * V / FP_LITERAL(299.)));
        const cellmodel_float_t dE_Ks_dK_i = -R * T / (F * (P_kna * Na_i + K_i));
        const cellmodel_float_t di_to_dE_K = -g_to * r * s;
        const cellmodel_float_t di_Ks_dE_Ks = -g_Ks * (Xs * Xs);
        const cellmodel_float_t dK_i_dt_linearized =
                Cm
                * (-(dE_K_dK_i * dalpha_K1_dE_K * dxK1_inf_dalpha_K1
                     + dE_K_dK_i * dbeta_K1_dE_K * dxK1_inf_dbeta_K1)
                           * di_K1_dxK1_inf
                   - dE_K_dK_i * di_K1_dE_K - dE_K_dK_i * di_Kr_dE_K - dE_K_dK_i * di_p_K_dE_K
                   - dE_K_dK_i * di_to_dE_K - dE_Ks_dK_i * di_Ks_dE_Ks)
                / (F * V_c);
        states[STATE_K_i * padded_num_cells + i] =
                K_i
                + (fabs(dK_i_dt_linearized) > FP_LITERAL(1.0e-8)
                           ? (FP_LITERAL(-1.0) + Exp(dt * dK_i_dt_linearized)) * dK_i_dt
                                     / dK_i_dt_linearized
                           : dt * dK_i_dt);
    }
}

// Compute a forward step using the Rush-Larsen scheme to the TP06 ODE
template <class LUT_type>
void step_RL(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
             const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
             const long num_cells, const long padded_num_cells, LUT_type &lut_V, LUT_type &lut_Ca)
{
    // Assign parameters
    const cellmodel_float_t celltype = parameters[PARAM_celltype];
    const cellmodel_float_t P_kna = parameters[PARAM_P_kna];
    const cellmodel_float_t g_K1 = parameters[PARAM_g_K1];
    const cellmodel_float_t g_Kr = parameters[PARAM_g_Kr];
    const cellmodel_float_t g_Ks = parameters[PARAM_g_Ks];
    const cellmodel_float_t g_Na = parameters[PARAM_g_Na];
    const cellmodel_float_t g_bna = parameters[PARAM_g_bna];
    const cellmodel_float_t g_CaL = parameters[PARAM_g_CaL];
    const cellmodel_float_t i_CaL_lim_delta = parameters[PARAM_i_CaL_lim_delta];
    const cellmodel_float_t g_bca = parameters[PARAM_g_bca];
    const cellmodel_float_t g_to = parameters[PARAM_g_to];
    const cellmodel_float_t K_mNa = parameters[PARAM_K_mNa];
    const cellmodel_float_t K_mk = parameters[PARAM_K_mk];
    const cellmodel_float_t P_NaK = parameters[PARAM_P_NaK];
    const cellmodel_float_t K_NaCa = parameters[PARAM_K_NaCa];
    const cellmodel_float_t K_sat = parameters[PARAM_K_sat];
    const cellmodel_float_t Km_Ca = parameters[PARAM_Km_Ca];
    const cellmodel_float_t Km_Nai = parameters[PARAM_Km_Nai];
    const cellmodel_float_t alpha = parameters[PARAM_alpha];
    const cellmodel_float_t gamma = parameters[PARAM_gamma];
    const cellmodel_float_t K_pCa = parameters[PARAM_K_pCa];
    const cellmodel_float_t g_pCa = parameters[PARAM_g_pCa];
    const cellmodel_float_t g_pK = parameters[PARAM_g_pK];
    const cellmodel_float_t Buf_c = parameters[PARAM_Buf_c];
    const cellmodel_float_t Buf_sr = parameters[PARAM_Buf_sr];
    const cellmodel_float_t Buf_ss = parameters[PARAM_Buf_ss];
    const cellmodel_float_t Ca_o = parameters[PARAM_Ca_o];
    const cellmodel_float_t EC = parameters[PARAM_EC];
    const cellmodel_float_t K_buf_c = parameters[PARAM_K_buf_c];
    const cellmodel_float_t K_buf_sr = parameters[PARAM_K_buf_sr];
    const cellmodel_float_t K_buf_ss = parameters[PARAM_K_buf_ss];
    const cellmodel_float_t K_up = parameters[PARAM_K_up];
    const cellmodel_float_t V_leak = parameters[PARAM_V_leak];
    const cellmodel_float_t V_rel = parameters[PARAM_V_rel];
    const cellmodel_float_t V_sr = parameters[PARAM_V_sr];
    const cellmodel_float_t V_ss = parameters[PARAM_V_ss];
    const cellmodel_float_t V_xfer = parameters[PARAM_V_xfer];
    const cellmodel_float_t Vmax_up = parameters[PARAM_Vmax_up];
    const cellmodel_float_t k1_prime = parameters[PARAM_k1_prime];
    const cellmodel_float_t k2_prime = parameters[PARAM_k2_prime];
    const cellmodel_float_t k3 = parameters[PARAM_k3];
    const cellmodel_float_t k4 = parameters[PARAM_k4];
    const cellmodel_float_t max_sr = parameters[PARAM_max_sr];
    const cellmodel_float_t min_sr = parameters[PARAM_min_sr];
    const cellmodel_float_t Na_o = parameters[PARAM_Na_o];
    const cellmodel_float_t Cm = parameters[PARAM_Cm];
    const cellmodel_float_t F = parameters[PARAM_F];
    const cellmodel_float_t R = parameters[PARAM_R];
    const cellmodel_float_t T = parameters[PARAM_T];
    const cellmodel_float_t V_c = parameters[PARAM_V_c];
    const cellmodel_float_t is_stimulated = parameters[PARAM_is_stimulated];
    const cellmodel_float_t stim_amplitude = parameters[PARAM_stim_amplitude];
    const cellmodel_float_t K_o = parameters[PARAM_K_o];

#if defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp parallel for simd simdlen(VECTOR_LENGTH)
#else
#pragma omp parallel for simd
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp parallel for
#endif
    for (long i = 0; i < num_cells; i++) {
        // Assign states
        const cellmodel_float_t Xr1 = states[STATE_Xr1 * padded_num_cells + i];
        const cellmodel_float_t Xr2 = states[STATE_Xr2 * padded_num_cells + i];
        const cellmodel_float_t Xs = states[STATE_Xs * padded_num_cells + i];
        const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
        const cellmodel_float_t h = states[STATE_h * padded_num_cells + i];
        const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
        const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
        const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
        const cellmodel_float_t f2 = states[STATE_f2 * padded_num_cells + i];
        const cellmodel_float_t fCass = states[STATE_fCass * padded_num_cells + i];
        const cellmodel_float_t s = states[STATE_s * padded_num_cells + i];
        const cellmodel_float_t r = states[STATE_r * padded_num_cells + i];
        const cellmodel_float_t Ca_i = states[STATE_Ca_i * padded_num_cells + i];
        const cellmodel_float_t R_prime = states[STATE_R_prime * padded_num_cells + i];
        const cellmodel_float_t Ca_SR = states[STATE_Ca_SR * padded_num_cells + i];
        const cellmodel_float_t Ca_ss = states[STATE_Ca_ss * padded_num_cells + i];
        const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];
        const cellmodel_float_t V = states[STATE_V * padded_num_cells + i];
        const cellmodel_float_t K_i = states[STATE_K_i * padded_num_cells + i];

        const auto lut_V_state = lut_V.compute_input_state(V);
        const auto lut_Ca_state = lut_Ca.compute_input_state(Ca_ss);

        // Expressions for the Reversal potentials component
        const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
        const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
        const cellmodel_float_t E_Ks = R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
        const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

        // Expressions for the Inward rectifier potassium current component
        const cellmodel_float_t alpha_K1 =
                FP_LITERAL(0.1)
                / (FP_LITERAL(1.)
                   + FP_LITERAL(6.14421235332821e-6)
                             * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
        const cellmodel_float_t beta_K1 =
                (FP_LITERAL(0.367879441171442) * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
                 + FP_LITERAL(3.06060402008027)
                           * Exp(FP_LITERAL(0.0002) * V - FP_LITERAL(0.0002) * E_K))
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V));
        const cellmodel_float_t xK1_inf = alpha_K1 / (alpha_K1 + beta_K1);
        const cellmodel_float_t i_K1 =
                FP_LITERAL(0.430331482911935) * g_K1 * sqrt(K_o) * (-E_K + V) * xK1_inf;

        // Expressions for the Rapid time dependent potassium current component
        const cellmodel_float_t i_Kr =
                FP_LITERAL(0.430331482911935) * g_Kr * sqrt(K_o) * (-E_K + V) * Xr1 * Xr2;

        // Expressions for the Xr1 gate component
#ifdef USE_LUT
        const cellmodel_float_t Xr1_RL_A = lut_V.lookup(LUT_INDEX_Xr1_RL_A, lut_V_state);
        const cellmodel_float_t Xr1_RL_B = lut_V.lookup(LUT_INDEX_Xr1_RL_B, lut_V_state);
        states[STATE_Xr1 * padded_num_cells + i] = Xr1_RL_A + Xr1_RL_B * Xr1;
#else
        const cellmodel_float_t xr1_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
        const cellmodel_float_t alpha_xr1 =
                FP_LITERAL(450.)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
        const cellmodel_float_t beta_xr1 =
                FP_LITERAL(6.)
                / (FP_LITERAL(1.)
                   + Exp(FP_LITERAL(60.) / FP_LITERAL(23.) + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
        const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
        const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
        const cellmodel_float_t dXr1_dt_linearized = FP_LITERAL(-1.) / tau_xr1;
        states[STATE_Xr1 * padded_num_cells + i] =
                (fabs(dXr1_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXr1_dt_linearized)) * dXr1_dt
                                   / dXr1_dt_linearized
                         : dt * dXr1_dt)
                + Xr1;
#endif

        // Expressions for the Xr2 gate component
#ifdef USE_LUT
        const cellmodel_float_t Xr2_RL_A = lut_V.lookup(LUT_INDEX_Xr2_RL_A, lut_V_state);
        const cellmodel_float_t Xr2_RL_B = lut_V.lookup(LUT_INDEX_Xr2_RL_B, lut_V_state);
        states[STATE_Xr2 * padded_num_cells + i] = Xr2_RL_A + Xr2_RL_B * Xr2;
#else
        const cellmodel_float_t xr2_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
        const cellmodel_float_t alpha_xr2 =
                FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t beta_xr2 =
                FP_LITERAL(1.12) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
        const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
        const cellmodel_float_t dXr2_dt_linearized = FP_LITERAL(-1.) / tau_xr2;
        states[STATE_Xr2 * padded_num_cells + i] =
                (fabs(dXr2_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXr2_dt_linearized)) * dXr2_dt
                                   / dXr2_dt_linearized
                         : dt * dXr2_dt)
                + Xr2;
#endif

        // Expressions for the Slow time dependent potassium current component
        const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

        // Expressions for the Xs gate component
#ifdef USE_LUT
        const cellmodel_float_t Xs_RL_A = lut_V.lookup(LUT_INDEX_Xs_RL_A, lut_V_state);
        const cellmodel_float_t Xs_RL_B = lut_V.lookup(LUT_INDEX_Xs_RL_B, lut_V_state);
        states[STATE_Xs * padded_num_cells + i] = Xs_RL_A + Xs_RL_B * Xs;
#else
        const cellmodel_float_t xs_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
        const cellmodel_float_t alpha_xs =
                FP_LITERAL(1400.)
                / sqrt(FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t beta_xs =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
        const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
        const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
        const cellmodel_float_t dXs_dt_linearized = FP_LITERAL(-1.) / tau_xs;
        states[STATE_Xs * padded_num_cells + i] =
                (fabs(dXs_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dXs_dt_linearized)) * dXs_dt
                                   / dXs_dt_linearized
                         : dt * dXs_dt)
                + Xs;
#endif

        // Expressions for the Fast sodium current component
        const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

        // Expressions for the m gate component
#ifdef USE_LUT
        const cellmodel_float_t m_RL_A = lut_V.lookup(LUT_INDEX_m_RL_A, lut_V_state);
        const cellmodel_float_t m_RL_B = lut_V.lookup(LUT_INDEX_m_RL_B, lut_V_state);
        states[STATE_m * padded_num_cells + i] = m_RL_A + m_RL_B * m;
#else
        const cellmodel_float_t m_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                  - FP_LITERAL(100.) * V / FP_LITERAL(903.)))
                                           * (1.
                                              + Exp(FP_LITERAL(-5686.) / FP_LITERAL(903.)
                                                    - FP_LITERAL(100.) * V / FP_LITERAL(903.))));
        const cellmodel_float_t alpha_m =
                FP_LITERAL(1.0) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-12.) - V / FP_LITERAL(5.)));
        const cellmodel_float_t beta_m =
                FP_LITERAL(0.1) / (FP_LITERAL(1.) + Exp(FP_LITERAL(7.) + V / FP_LITERAL(5.)))
                + FP_LITERAL(0.1)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-1.) / FP_LITERAL(4.) + V / FP_LITERAL(200.)));
        const cellmodel_float_t tau_m = alpha_m * beta_m;
        const cellmodel_float_t dm_dt = (-m + m_inf) / tau_m;
        const cellmodel_float_t dm_dt_linearized = FP_LITERAL(-1.) / tau_m;
        states[STATE_m * padded_num_cells + i] =
                (fabs(dm_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dm_dt_linearized)) * dm_dt
                                   / dm_dt_linearized
                         : dt * dm_dt)
                + m;
#endif

        // Expressions for the h gate component
#ifdef USE_LUT
        const cellmodel_float_t h_RL_A = lut_V.lookup(LUT_INDEX_h_RL_A, lut_V_state);
        const cellmodel_float_t h_RL_B = lut_V.lookup(LUT_INDEX_h_RL_B, lut_V_state);
        states[STATE_h * padded_num_cells + i] = h_RL_A + h_RL_B * h;
#else
        const cellmodel_float_t h_inf = FP_LITERAL(1.0)
                                        / ((FP_LITERAL(1.)
                                            + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                  + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                                           * (1.
                                              + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                                    + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
        const cellmodel_float_t alpha_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(4.43126792958051e-7) * Exp(FP_LITERAL(-0.147058823529412) * V)
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_h =
                (V < FP_LITERAL(-40.)
                         ? FP_LITERAL(310000.) * Exp(FP_LITERAL(0.3485) * V)
                                   + FP_LITERAL(2.7) * Exp(FP_LITERAL(0.079) * V)
                         : FP_LITERAL(0.77)
                                   / (0.13
                                      + FP_LITERAL(0.0497581410839387)
                                                * Exp(FP_LITERAL(-0.0900900900900901) * V)));
        const cellmodel_float_t tau_h = FP_LITERAL(1.0) / (alpha_h + beta_h);
        const cellmodel_float_t dh_dt = (-h + h_inf) / tau_h;
        const cellmodel_float_t dh_dt_linearized = FP_LITERAL(-1.) / tau_h;
        states[STATE_h * padded_num_cells + i] =
                (fabs(dh_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dh_dt_linearized)) * dh_dt
                                   / dh_dt_linearized
                         : dt * dh_dt)
                + h;
#endif

        // Expressions for the j gate component
#ifdef USE_LUT
        const cellmodel_float_t j_RL_A = lut_V.lookup(LUT_INDEX_j_RL_A, lut_V_state);
        const cellmodel_float_t j_RL_B = lut_V.lookup(LUT_INDEX_j_RL_B, lut_V_state);
        states[STATE_j * padded_num_cells + i] = j_RL_A + j_RL_B * j;
#else
        const cellmodel_float_t j_inf =
                FP_LITERAL(1.0)
                / ((FP_LITERAL(1.)
                    + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V))
                   * (FP_LITERAL(1.)
                      + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V)));
        const cellmodel_float_t alpha_j =
                (V < FP_LITERAL(-40.)
                         ? (FP_LITERAL(37.78) + V)
                                   * (FP_LITERAL(-25428.) * Exp(FP_LITERAL(0.2444) * V)
                                      - FP_LITERAL(6.948e-6) * Exp(FP_LITERAL(-0.04391) * V))
                                   / (FP_LITERAL(1.)
                                      + FP_LITERAL(50262745825.954) * Exp(FP_LITERAL(0.311) * V))
                         : FP_LITERAL(0.));
        const cellmodel_float_t beta_j =
                (V < FP_LITERAL(-40.) ? FP_LITERAL(0.02424) * Exp(FP_LITERAL(-0.01052) * V)
                                                / (FP_LITERAL(1.)
                                                   + FP_LITERAL(0.00396086833990426)
                                                             * Exp(FP_LITERAL(-0.1378) * V))
                                      : FP_LITERAL(0.6) * Exp(FP_LITERAL(0.057) * V)
                                                / (FP_LITERAL(1.)
                                                   + Exp(FP_LITERAL(-16.) / FP_LITERAL(5.)
                                                         - V / FP_LITERAL(10.))));
        const cellmodel_float_t tau_j = FP_LITERAL(1.0) / (alpha_j + beta_j);
        const cellmodel_float_t dj_dt = (-j + j_inf) / tau_j;
        const cellmodel_float_t dj_dt_linearized = FP_LITERAL(-1.) / tau_j;
        states[STATE_j * padded_num_cells + i] =
                (fabs(dj_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dj_dt_linearized)) * dj_dt
                                   / dj_dt_linearized
                         : dt * dj_dt)
                + j;
#endif

        // Expressions for the Sodium background current component
        const cellmodel_float_t i_b_Na = g_bna * (-E_Na + V);

        // Expressions for the L_type Ca current component
#if 1 && defined(USE_LUT)
        const cellmodel_float_t i_CaL_fraction = lut_V.lookup(LUT_INDEX_i_CaL_fraction, lut_V_state);
        const cellmodel_float_t i_CaL_factors_e1 = lut_V.lookup(LUT_INDEX_i_CaL_factors_e1, lut_V_state);
        const cellmodel_float_t i_CaL_factors = 4. * F * g_CaL
                                 * (-Ca_o + Ca_ss * i_CaL_factors_e1) * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#else
        const cellmodel_float_t V_eff = FP_LITERAL(-15.) + V;
        const cellmodel_float_t i_CaL_factors =
                FP_LITERAL(4.) * F * g_CaL
                * (-Ca_o + Ca_ss * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) / FP_LITERAL(4.))
                * d * f * f2 * fCass;
        const cellmodel_float_t i_CaL_fraction =
                    (fabs(V_eff) < i_CaL_lim_delta
                             ? FP_LITERAL(0.5)
                             : F * V_eff / (R * T * (Expm1(FP_LITERAL(2.) * F * V_eff / (R * T)))));
        const cellmodel_float_t i_CaL = i_CaL_factors * i_CaL_fraction;
#endif

        // Expressions for the d gate component
#ifdef USE_LUT
        const cellmodel_float_t d_RL_A = lut_V.lookup(LUT_INDEX_d_RL_A, lut_V_state);
        const cellmodel_float_t d_RL_B = lut_V.lookup(LUT_INDEX_d_RL_B, lut_V_state);
        states[STATE_d * padded_num_cells + i] = d_RL_A + d_RL_B * d;
#else
        const cellmodel_float_t d_inf = FP_LITERAL(1.0)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(-16.) / FP_LITERAL(15.)
                                                 - FP_LITERAL(2.) * V / FP_LITERAL(15.)));
        const cellmodel_float_t alpha_d =
                FP_LITERAL(0.25)
                + FP_LITERAL(1.4)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(-35.) / FP_LITERAL(13.) - V / FP_LITERAL(13.)));
        const cellmodel_float_t beta_d =
                FP_LITERAL(1.4) / (FP_LITERAL(1.) + Exp(FP_LITERAL(1.) + V / FP_LITERAL(5.)));
        const cellmodel_float_t gamma_d =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(20.)));
        const cellmodel_float_t tau_d = alpha_d * beta_d + gamma_d;
        const cellmodel_float_t dd_dt = (-d + d_inf) / tau_d;
        const cellmodel_float_t dd_dt_linearized = FP_LITERAL(-1.) / tau_d;
        states[STATE_d * padded_num_cells + i] =
                (fabs(dd_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dd_dt_linearized)) * dd_dt
                                   / dd_dt_linearized
                         : dt * dd_dt)
                + d;
#endif

        // Expressions for the f gate component
#ifdef USE_LUT
        const cellmodel_float_t f_RL_A = lut_V.lookup(LUT_INDEX_f_RL_A, lut_V_state);
        const cellmodel_float_t f_RL_B = lut_V.lookup(LUT_INDEX_f_RL_B, lut_V_state);
        states[STATE_f * padded_num_cells + i] = f_RL_A + f_RL_B * f;
#else
        const cellmodel_float_t f_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(20.) / FP_LITERAL(7.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f =
                FP_LITERAL(20.)
                + FP_LITERAL(180.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(200.)
                          / (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(13.) / FP_LITERAL(10.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(1102.5)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(225.));
        const cellmodel_float_t df_dt = (-f + f_inf) / tau_f;
        const cellmodel_float_t df_dt_linearized = FP_LITERAL(-1.) / tau_f;
        states[STATE_f * padded_num_cells + i] =
                (fabs(df_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * df_dt_linearized)) * df_dt
                                   / df_dt_linearized
                         : dt * df_dt)
                + f;
#endif

        // Expressions for the F2 gate component
#ifdef USE_LUT
        const cellmodel_float_t f2_RL_A = lut_V.lookup(LUT_INDEX_f2_RL_A, lut_V_state);
        const cellmodel_float_t f2_RL_B = lut_V.lookup(LUT_INDEX_f2_RL_B, lut_V_state);
        states[STATE_f2 * padded_num_cells + i] = f2_RL_A + f2_RL_B * f2;
#else
        const cellmodel_float_t f2_inf =
                FP_LITERAL(0.33)
                + FP_LITERAL(0.67) / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
        const cellmodel_float_t tau_f2 =
                FP_LITERAL(31.)
                        / (FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)))
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                + FP_LITERAL(562.)
                          * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                / FP_LITERAL(240.));
        const cellmodel_float_t df2_dt = (-f2 + f2_inf) / tau_f2;
        const cellmodel_float_t df2_dt_linearized = FP_LITERAL(-1.) / tau_f2;
        states[STATE_f2 * padded_num_cells + i] =
                (fabs(df2_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * df2_dt_linearized)) * df2_dt
                                   / df2_dt_linearized
                         : dt * df2_dt)
                + f2;
#endif

        // Expressions for the FCass gate component
#if 1 && defined(USE_LUT)
        const cellmodel_float_t fCass_RL_A = lut_Ca.lookup(LUT_Ca_INDEX_fCass_RL_A, lut_Ca_state);
        const cellmodel_float_t fCass_RL_B = lut_Ca.lookup(LUT_Ca_INDEX_fCass_RL_B, lut_Ca_state);
        states[STATE_fCass * padded_num_cells + i] = fCass_RL_A + fCass_RL_B * fCass;
#else
        const cellmodel_float_t fCass_inf =
                FP_LITERAL(0.4)
                + FP_LITERAL(0.6) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t tau_fCass =
                FP_LITERAL(2.)
                + FP_LITERAL(80.) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
        const cellmodel_float_t dfCass_dt = (-fCass + fCass_inf) / tau_fCass;
        const cellmodel_float_t dfCass_dt_linearized = FP_LITERAL(-1.) / tau_fCass;
        states[STATE_fCass * padded_num_cells + i] =
                (fabs(dfCass_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dfCass_dt_linearized)) * dfCass_dt
                                   / dfCass_dt_linearized
                         : dt * dfCass_dt)
                + fCass;
#endif

        // Expressions for the Calcium background current component
        const cellmodel_float_t i_b_Ca = g_bca * (-E_Ca + V);

        // Expressions for the Transient outward current component
        const cellmodel_float_t i_to = g_to * (-E_K + V) * r * s;

        // Expressions for the s gate component
#ifdef USE_LUT
        const cellmodel_float_t s_RL_A = lut_V.lookup(LUT_INDEX_s_RL_A, lut_V_state);
        const cellmodel_float_t s_RL_B = lut_V.lookup(LUT_INDEX_s_RL_B, lut_V_state);
        states[STATE_s * padded_num_cells + i] = s_RL_A + s_RL_B * s;
#else
        const cellmodel_float_t s_inf =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.)
                                      + Exp(FP_LITERAL(28.) / FP_LITERAL(5.) + V / FP_LITERAL(5.)))
                         : FP_LITERAL(1.0)
                                   / (FP_LITERAL(1.) + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
        const cellmodel_float_t tau_s =
                (celltype == FP_LITERAL(2.)
                         ? FP_LITERAL(8.)
                                   + FP_LITERAL(1000.)
                                             * Exp(-((FP_LITERAL(67.) + V) * (FP_LITERAL(67.) + V))
                                                   / FP_LITERAL(1000.))
                         : FP_LITERAL(3.)
                                   + FP_LITERAL(5.)
                                             / (FP_LITERAL(1.)
                                                + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                   + FP_LITERAL(85.)
                                             * Exp(-((FP_LITERAL(45.) + V) * (FP_LITERAL(45.) + V))
                                                   / FP_LITERAL(320.)));
        const cellmodel_float_t ds_dt = (-s + s_inf) / tau_s;
        const cellmodel_float_t ds_dt_linearized = FP_LITERAL(-1.) / tau_s;
        states[STATE_s * padded_num_cells + i] =
                (fabs(ds_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * ds_dt_linearized)) * ds_dt
                                   / ds_dt_linearized
                         : dt * ds_dt)
                + s;
#endif

        // Expressions for the r gate component
#ifdef USE_LUT
        const cellmodel_float_t r_RL_A = lut_V.lookup(LUT_INDEX_r_RL_A, lut_V_state);
        const cellmodel_float_t r_RL_B = lut_V.lookup(LUT_INDEX_r_RL_B, lut_V_state);
        states[STATE_r * padded_num_cells + i] = r_RL_A + r_RL_B * r;
#else
        const cellmodel_float_t r_inf =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Exp(FP_LITERAL(10.) / FP_LITERAL(3.) - V / FP_LITERAL(6.)));
        const cellmodel_float_t tau_r =
                FP_LITERAL(0.8)
                + FP_LITERAL(9.5)
                          * Exp(-((FP_LITERAL(40.) + V) * (FP_LITERAL(40.) + V))
                                / FP_LITERAL(1800.));
        const cellmodel_float_t dr_dt = (-r + r_inf) / tau_r;
        const cellmodel_float_t dr_dt_linearized = FP_LITERAL(-1.) / tau_r;
        states[STATE_r * padded_num_cells + i] =
                (fabs(dr_dt_linearized) > FP_LITERAL(1.0e-8)
                         ? (FP_LITERAL(-1.0) + Exp(dt * dr_dt_linearized)) * dr_dt
                                   / dr_dt_linearized
                         : dt * dr_dt)
                + r;
#endif

        // Expressions for the Sodium potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_iNaK = lut_V.lookup(LUT_INDEX_rec_iNaK, lut_V_state);
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i / ((K_mNa + Na_i) * (K_mk + K_o)) * rec_iNaK;
#else
        const cellmodel_float_t i_NaK =
                K_o * P_NaK * Na_i
                / ((K_mNa + Na_i) * (K_mk + K_o)
                   * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                      + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));
#endif

        // Expressions for the Sodium calcium exchanger current component
#ifdef USE_LUT
        const cellmodel_float_t i_NaCa_A = lut_V.lookup(LUT_INDEX_I_NaCa_A, lut_V_state);
        const cellmodel_float_t i_NaCa_B = lut_V.lookup(LUT_INDEX_I_NaCa_B, lut_V_state);
        const cellmodel_float_t i_NaCa = (i_NaCa_A * (Na_i * Na_i * Na_i)) - (i_NaCa_B * Ca_i);
#else
        const cellmodel_float_t i_NaCa =
                K_NaCa
                * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                   - alpha * (Na_o * Na_o * Na_o) * Ca_i
                             * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                   * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));
#endif

        // Expressions for the Calcium pump current component
        const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

        // Expressions for the Potassium pump current component
#ifdef USE_LUT
        const cellmodel_float_t rec_ipK = lut_V.lookup(LUT_INDEX_rec_ipK, lut_V_state);
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V) * rec_ipK;
#else
        const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                        / (FP_LITERAL(1.)
                                           + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                 - FP_LITERAL(50.) * V / FP_LITERAL(299.)));
#endif

        // Expressions for the Calcium dynamics component
        const cellmodel_float_t i_up = Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
        const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
        const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
        const cellmodel_float_t kcasr =
                max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
        const cellmodel_float_t Ca_i_bufc =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
        const cellmodel_float_t Ca_sr_bufsr =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
        const cellmodel_float_t Ca_ss_bufss =
                FP_LITERAL(1.0)
                / (FP_LITERAL(1.) + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
        const cellmodel_float_t dCa_i_dt =
                (V_sr * (-i_up + i_leak) / V_c
                 - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca) / (FP_LITERAL(2.) * F * V_c)
                 + i_xfer)
                * Ca_i_bufc;
        states[STATE_Ca_i * padded_num_cells + i] = dt * dCa_i_dt + Ca_i;
        const cellmodel_float_t k1 = k1_prime / kcasr;
        const cellmodel_float_t k2 = k2_prime * kcasr;
        const cellmodel_float_t O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
        const cellmodel_float_t dR_prime_dt =
                k4 * (FP_LITERAL(1.) - R_prime) - Ca_ss * R_prime * k2;
        states[STATE_R_prime * padded_num_cells + i] = dt * dR_prime_dt + R_prime;
        const cellmodel_float_t i_rel = V_rel * (-Ca_ss + Ca_SR) * O;
        const cellmodel_float_t dCa_SR_dt = (-i_leak - i_rel + i_up) * Ca_sr_bufsr;
        states[STATE_Ca_SR * padded_num_cells + i] = dt * dCa_SR_dt + Ca_SR;
        const cellmodel_float_t dCa_ss_dt = (V_sr * i_rel / V_ss - V_c * i_xfer / V_ss
                                             - Cm * i_CaL / (FP_LITERAL(2.) * F * V_ss))
                                            * Ca_ss_bufss;
        states[STATE_Ca_ss * padded_num_cells + i] = dt * dCa_ss_dt + Ca_ss;

        // Expressions for the Sodium dynamics component
        const cellmodel_float_t dNa_i_dt =
                Cm * (-i_Na - i_b_Na - FP_LITERAL(3.) * i_NaCa - FP_LITERAL(3.) * i_NaK)
                / (F * V_c);
        states[STATE_Na_i * padded_num_cells + i] = dt * dNa_i_dt + Na_i;

        // Expressions for the Membrane component
        const cellmodel_float_t i_Stim = (is_stimulated ? -stim_amplitude : FP_LITERAL(0.));
        const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK - i_Stim
                                        - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
        states[STATE_V * padded_num_cells + i] = dt * dV_dt + V;

        // Expressions for the Potassium dynamics component
        const cellmodel_float_t dK_i_dt =
                Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + FP_LITERAL(2.) * i_NaK)
                / (F * V_c);
        states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;
    }
}

std::vector<univariate_func>
expressions_tuple_to_func_vector(std::vector<univariate_func_tuple> e_tuple_vec)
{
    std::vector<univariate_func> e_func_vec(e_tuple_vec.size());
    for (size_t i = 0; i < e_tuple_vec.size(); i++) {
        e_func_vec[i] = e_tuple_vec[i].f;
    }
    return e_func_vec;
}

const std::vector<univariate_func_tuple> expressions_V_tuple_vec(expressions_V.begin(),
                                                                 expressions_V.end());
const std::vector<univariate_func_tuple> expressions_Ca_ss_tuple_vec(expressions_Ca.begin(),
                                                                     expressions_Ca.end());

const std::vector<univariate_func> expressions_V_vec =
        expressions_tuple_to_func_vector(expressions_V_tuple_vec);
const std::vector<univariate_func> expressions_Ca_ss_vec =
        expressions_tuple_to_func_vector(expressions_Ca_ss_tuple_vec);

const struct cellmodel_lut model_TP06_lut = {
        &init_state_values,
        &init_parameters_values,
        &state_index,
        &parameter_index,
        &step_FE<default_LUT_type>,
        &step_RL<default_LUT_type>,
        &step_GRL1<default_LUT_type>,
        NUM_STATES,
        NUM_PARAMS,
        &expressions_V_vec,
        &expressions_Ca_ss_vec,
};
