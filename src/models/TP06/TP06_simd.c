#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "TP06" model

#include "TP06_simd.h"
#include "cellmodel.h"

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

// Compute a forward step using the explicit Euler algorithm to the TP06 ODE
static void step_FE(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
                    const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
                    const long num_cells, long padded_num_cells)
{
    #pragma omp parallel
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

#if defined(HINT_CLANG_SIMD)
#pragma omp for
#pragma clang loop vectorize(assume_safety)
#elif defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES) simdlen(VECTOR_LENGTH)
#else
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES)
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp for
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

            // Expressions for the Reversal potentials component
            const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
            const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
            const cellmodel_float_t E_Ks =
                    R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
            const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

            // Expressions for the Inward rectifier potassium current component
            const cellmodel_float_t alpha_K1 =
                    FP_LITERAL(0.1)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(6.14421235332821e-6)
                                 * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
            const cellmodel_float_t beta_K1 =
                    (FP_LITERAL(0.367879441171442)
                             * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
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
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
            const cellmodel_float_t alpha_xr1 =
                    FP_LITERAL(450.)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
            const cellmodel_float_t beta_xr1 = FP_LITERAL(6.)
                                               / (FP_LITERAL(1.)
                                                  + Exp(FP_LITERAL(60.) / FP_LITERAL(23.)
                                                        + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
            const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
            const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
            states[STATE_Xr1 * padded_num_cells + i] = dt * dXr1_dt + Xr1;

            // Expressions for the Xr2 gate component
            const cellmodel_float_t xr2_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
            const cellmodel_float_t alpha_xr2 =
                    FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
            const cellmodel_float_t beta_xr2 =
                    FP_LITERAL(1.12)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
            const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
            const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
            states[STATE_Xr2 * padded_num_cells + i] = dt * dXr2_dt + Xr2;

            // Expressions for the Slow time dependent potassium current component
            const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

            // Expressions for the Xs gate component
            const cellmodel_float_t xs_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
            const cellmodel_float_t alpha_xs =
                    FP_LITERAL(1400.)
                    / sqrt(FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
            const cellmodel_float_t beta_xs =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
            const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
            const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
            states[STATE_Xs * padded_num_cells + i] = dt * dXs_dt + Xs;

            // Expressions for the Fast sodium current component
            const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

            // Expressions for the m gate component
            const cellmodel_float_t m_inf =
                    FP_LITERAL(1.0)
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
            const cellmodel_float_t h_inf =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                              + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                       * (1.
                          + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
            const cellmodel_float_t alpha_h =
                    (V < FP_LITERAL(-40.) ? FP_LITERAL(4.43126792958051e-7)
                                                    * Exp(FP_LITERAL(-0.147058823529412) * V)
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
            const cellmodel_float_t j_inf =
                    FP_LITERAL(1.0)
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
                                          + FP_LITERAL(50262745825.954)
                                                    * Exp(FP_LITERAL(0.311) * V))
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
                    + FP_LITERAL(180.)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
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
                    + FP_LITERAL(0.67)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
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
                                          + Exp(FP_LITERAL(28.) / FP_LITERAL(5.)
                                                + V / FP_LITERAL(5.)))
                             : FP_LITERAL(1.0)
                                       / (FP_LITERAL(1.)
                                          + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
            const cellmodel_float_t tau_s =
                    (celltype == FP_LITERAL(2.)
                             ? FP_LITERAL(8.)
                                       + FP_LITERAL(1000.)
                                                 * Exp(-((FP_LITERAL(67.) + V)
                                                         * (FP_LITERAL(67.) + V))
                                                       / FP_LITERAL(1000.))
                             : FP_LITERAL(3.)
                                       + FP_LITERAL(5.)
                                                 / (FP_LITERAL(1.)
                                                    + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                       + FP_LITERAL(85.)
                                                 * Exp(-((FP_LITERAL(45.) + V)
                                                         * (FP_LITERAL(45.) + V))
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
            const cellmodel_float_t i_NaK =
                    K_o * P_NaK * Na_i
                    / ((K_mNa + Na_i) * (K_mk + K_o)
                       * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                          + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));

            // Expressions for the Sodium calcium exchanger current component
            const cellmodel_float_t i_NaCa =
                    K_NaCa
                    * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                       - alpha * (Na_o * Na_o * Na_o) * Ca_i
                                 * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                    / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                       * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));

            // Expressions for the Calcium pump current component
            const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

            // Expressions for the Potassium pump current component
            const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                            / (FP_LITERAL(1.)
                                               + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                     - FP_LITERAL(50.) * V / FP_LITERAL(299.)));

            // Expressions for the Calcium dynamics component
            const cellmodel_float_t i_up =
                    Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
            const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
            const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
            const cellmodel_float_t kcasr =
                    max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
            const cellmodel_float_t Ca_i_bufc =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
            const cellmodel_float_t Ca_sr_bufsr =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
            const cellmodel_float_t Ca_ss_bufss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
            const cellmodel_float_t dCa_i_dt = (V_sr * (-i_up + i_leak) / V_c
                                                - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca)
                                                          / (FP_LITERAL(2.) * F * V_c)
                                                + i_xfer)
                                               * Ca_i_bufc;
            states[STATE_Ca_i * padded_num_cells + i] = dt * dCa_i_dt + Ca_i;
            const cellmodel_float_t k1 = k1_prime / kcasr;
            const cellmodel_float_t k2 = k2_prime * kcasr;
            const cellmodel_float_t O =
                    (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
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
            const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK
                                            - i_Stim - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
            states[STATE_V * padded_num_cells + i] = dt * dV_dt + V;

            // Expressions for the Potassium dynamics component
            const cellmodel_float_t dK_i_dt =
                    Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + FP_LITERAL(2.) * i_NaK)
                    / (F * V_c);
            states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;
        }
    }
}

// Compute a forward step using the GRL1 scheme to the TP06 ODE
static void step_GRL1(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
                      const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
                      const long num_cells, long padded_num_cells)
{
    #pragma omp parallel
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


#if defined(HINT_CLANG_SIMD)
#pragma omp for
#pragma clang loop vectorize(assume_safety)
#elif defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES) simdlen(VECTOR_LENGTH)
#else
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES)
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp for
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

            // Expressions for the Reversal potentials component
            const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
            const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
            const cellmodel_float_t E_Ks =
                    R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
            const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

            // Expressions for the Inward rectifier potassium current component
            const cellmodel_float_t alpha_K1 =
                    FP_LITERAL(0.1)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(6.14421235332821e-6)
                                 * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
            const cellmodel_float_t beta_K1 =
                    (FP_LITERAL(0.367879441171442)
                             * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
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
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
            const cellmodel_float_t alpha_xr1 =
                    FP_LITERAL(450.)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
            const cellmodel_float_t beta_xr1 = FP_LITERAL(6.)
                                               / (FP_LITERAL(1.)
                                                  + Exp(FP_LITERAL(60.) / FP_LITERAL(23.)
                                                        + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
            const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
            const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
            const cellmodel_float_t dXr1_dt_linearized = FP_LITERAL(-1.) / tau_xr1;
            states[STATE_Xr1 * padded_num_cells + i] =
                    (Expm1(dt * dXr1_dt_linearized)) * dXr1_dt / dXr1_dt_linearized + Xr1;

            // Expressions for the Xr2 gate component
            const cellmodel_float_t xr2_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
            const cellmodel_float_t alpha_xr2 =
                    FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
            const cellmodel_float_t beta_xr2 =
                    FP_LITERAL(1.12)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
            const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
            const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
            const cellmodel_float_t dXr2_dt_linearized = FP_LITERAL(-1.) / tau_xr2;
            states[STATE_Xr2 * padded_num_cells + i] =
                    (Expm1(dt * dXr2_dt_linearized)) * dXr2_dt / dXr2_dt_linearized + Xr2;

            // Expressions for the Slow time dependent potassium current component
            const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

            // Expressions for the Xs gate component
            const cellmodel_float_t xs_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
            const cellmodel_float_t alpha_xs =
                    FP_LITERAL(1400.)
                    / sqrt(FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
            const cellmodel_float_t beta_xs =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
            const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
            const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
            const cellmodel_float_t dXs_dt_linearized = FP_LITERAL(-1.) / tau_xs;
            states[STATE_Xs * padded_num_cells + i] =
                    (Expm1(dt * dXs_dt_linearized)) * dXs_dt / dXs_dt_linearized + Xs;

            // Expressions for the Fast sodium current component
            const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

            // Expressions for the m gate component
            const cellmodel_float_t m_inf =
                    FP_LITERAL(1.0)
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
                    (Expm1(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized + m;

            // Expressions for the h gate component
            const cellmodel_float_t h_inf =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                              + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                       * (1.
                          + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
            const cellmodel_float_t alpha_h =
                    (V < FP_LITERAL(-40.) ? FP_LITERAL(4.43126792958051e-7)
                                                    * Exp(FP_LITERAL(-0.147058823529412) * V)
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
                    (Expm1(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized + h;

            // Expressions for the j gate component
            const cellmodel_float_t j_inf =
                    FP_LITERAL(1.0)
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
                                          + FP_LITERAL(50262745825.954)
                                                    * Exp(FP_LITERAL(0.311) * V))
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
                    (Expm1(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized + j;

            // Expressions for the Sodium background current component
            const cellmodel_float_t i_b_Na = g_bna * (-E_Na + V);

            // Expressions for the L_type Ca current component
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
            const cellmodel_float_t dd_dt_linearized = FP_LITERAL(-1.) / tau_d;
            states[STATE_d * padded_num_cells + i] =
                    (Expm1(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized + d;

            // Expressions for the f gate component
            const cellmodel_float_t f_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(20.) / FP_LITERAL(7.) + V / FP_LITERAL(7.)));
            const cellmodel_float_t tau_f =
                    FP_LITERAL(20.)
                    + FP_LITERAL(180.)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                    + FP_LITERAL(200.)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(13.) / FP_LITERAL(10.) - V / FP_LITERAL(10.)))
                    + FP_LITERAL(1102.5)
                              * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                    / FP_LITERAL(225.));
            const cellmodel_float_t df_dt = (-f + f_inf) / tau_f;
            const cellmodel_float_t df_dt_linearized = FP_LITERAL(-1.) / tau_f;
            states[STATE_f * padded_num_cells + i] =
                    (Expm1(dt * df_dt_linearized)) * df_dt / df_dt_linearized + f;

            // Expressions for the F2 gate component
            const cellmodel_float_t f2_inf =
                    FP_LITERAL(0.33)
                    + FP_LITERAL(0.67)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
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
                    (Expm1(dt * df2_dt_linearized)) * df2_dt / df2_dt_linearized + f2;

            // Expressions for the FCass gate component
            const cellmodel_float_t fCass_inf =
                    FP_LITERAL(0.4)
                    + FP_LITERAL(0.6) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
            const cellmodel_float_t tau_fCass =
                    FP_LITERAL(2.)
                    + FP_LITERAL(80.) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
            const cellmodel_float_t dfCass_dt = (-fCass + fCass_inf) / tau_fCass;
            const cellmodel_float_t dfCass_dt_linearized = FP_LITERAL(-1.) / tau_fCass;
            states[STATE_fCass * padded_num_cells + i] =
                    (Expm1(dt * dfCass_dt_linearized)) * dfCass_dt / dfCass_dt_linearized + fCass;

            // Expressions for the Calcium background current component
            const cellmodel_float_t i_b_Ca = g_bca * (-E_Ca + V);

            // Expressions for the Transient outward current component
            const cellmodel_float_t i_to = g_to * (-E_K + V) * r * s;

            // Expressions for the s gate component
            const cellmodel_float_t s_inf =
                    (celltype == FP_LITERAL(2.)
                             ? FP_LITERAL(1.0)
                                       / (FP_LITERAL(1.)
                                          + Exp(FP_LITERAL(28.) / FP_LITERAL(5.)
                                                + V / FP_LITERAL(5.)))
                             : FP_LITERAL(1.0)
                                       / (FP_LITERAL(1.)
                                          + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
            const cellmodel_float_t tau_s =
                    (celltype == FP_LITERAL(2.)
                             ? FP_LITERAL(8.)
                                       + FP_LITERAL(1000.)
                                                 * Exp(-((FP_LITERAL(67.) + V)
                                                         * (FP_LITERAL(67.) + V))
                                                       / FP_LITERAL(1000.))
                             : FP_LITERAL(3.)
                                       + FP_LITERAL(5.)
                                                 / (FP_LITERAL(1.)
                                                    + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                       + FP_LITERAL(85.)
                                                 * Exp(-((FP_LITERAL(45.) + V)
                                                         * (FP_LITERAL(45.) + V))
                                                       / FP_LITERAL(320.)));
            const cellmodel_float_t ds_dt = (-s + s_inf) / tau_s;
            const cellmodel_float_t ds_dt_linearized = FP_LITERAL(-1.) / tau_s;
            states[STATE_s * padded_num_cells + i] =
                    (Expm1(dt * ds_dt_linearized)) * ds_dt / ds_dt_linearized + s;

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
            const cellmodel_float_t dr_dt_linearized = FP_LITERAL(-1.) / tau_r;
            states[STATE_r * padded_num_cells + i] =
                    (Expm1(dt * dr_dt_linearized)) * dr_dt / dr_dt_linearized + r;

            // Expressions for the Sodium potassium pump current component
            const cellmodel_float_t i_NaK =
                    K_o * P_NaK * Na_i
                    / ((K_mNa + Na_i) * (K_mk + K_o)
                       * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                          + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));

            // Expressions for the Sodium calcium exchanger current component
            const cellmodel_float_t i_NaCa =
                    K_NaCa
                    * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                       - alpha * (Na_o * Na_o * Na_o) * Ca_i
                                 * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                    / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                       * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));

            // Expressions for the Calcium pump current component
            const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

            // Expressions for the Potassium pump current component
            const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                            / (FP_LITERAL(1.)
                                               + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                     - FP_LITERAL(50.) * V / FP_LITERAL(299.)));

            // Expressions for the Calcium dynamics component
            const cellmodel_float_t i_up =
                    Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
            const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
            const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
            const cellmodel_float_t kcasr =
                    max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
            const cellmodel_float_t Ca_i_bufc =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
            const cellmodel_float_t Ca_sr_bufsr =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
            const cellmodel_float_t Ca_ss_bufss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
            const cellmodel_float_t dCa_i_dt = (V_sr * (-i_up + i_leak) / V_c
                                                - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca)
                                                          / (FP_LITERAL(2.) * F * V_c)
                                                + i_xfer)
                                               * Ca_i_bufc;
            const cellmodel_float_t di_p_Ca_dCa_i =
                    g_pCa / (K_pCa + Ca_i) - g_pCa * Ca_i / ((K_pCa + Ca_i) * (K_pCa + Ca_i));
            const cellmodel_float_t dE_Ca_dCa_i = FP_LITERAL(-0.5) * R * T / (F * Ca_i);
            const cellmodel_float_t dCa_i_bufc_dCa_i =
                    FP_LITERAL(2.) * Buf_c * K_buf_c
                    / (((FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
                        * (FP_LITERAL(1.)
                           + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i))))
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
                       - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca)
                                 / (FP_LITERAL(2.) * F * V_c)
                       + i_xfer)
                              * dCa_i_bufc_dCa_i;
            states[STATE_Ca_i * padded_num_cells + i] =
                    Ca_i
                    + (fabs(dCa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (Expm1(dt * dCa_i_dt_linearized)) * dCa_i_dt / dCa_i_dt_linearized
                               : dt * dCa_i_dt);
            const cellmodel_float_t k1 = k1_prime / kcasr;
            const cellmodel_float_t k2 = k2_prime * kcasr;
            const cellmodel_float_t O =
                    (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
            const cellmodel_float_t dR_prime_dt =
                    k4 * (FP_LITERAL(1.) - R_prime) - Ca_ss * R_prime * k2;
            const cellmodel_float_t dR_prime_dt_linearized = -k4 - Ca_ss * k2;
            states[STATE_R_prime * padded_num_cells + i] =
                    (fabs(dR_prime_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (Expm1(dt * dR_prime_dt_linearized)) * dR_prime_dt
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
                    / (((FP_LITERAL(1.)
                         + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
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
                               ? (Expm1(dt * dCa_SR_dt_linearized)) * dCa_SR_dt
                                         / dCa_SR_dt_linearized
                               : dt * dCa_SR_dt);
            const cellmodel_float_t dCa_ss_dt = (V_sr * i_rel / V_ss - V_c * i_xfer / V_ss
                                                 - Cm * i_CaL / (FP_LITERAL(2.) * F * V_ss))
                                                * Ca_ss_bufss;
            const cellmodel_float_t dCa_ss_bufss_dCa_ss =
                    FP_LITERAL(2.) * Buf_ss * K_buf_ss
                    / (((FP_LITERAL(1.)
                         + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
                        * (FP_LITERAL(1.)
                           + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss))))
                       * ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
            const cellmodel_float_t dO_dCa_ss =
                    FP_LITERAL(-2.) * (Ca_ss * Ca_ss * Ca_ss) * (k1 * k1) * R_prime
                            / ((k3 + (Ca_ss * Ca_ss) * k1) * (k3 + (Ca_ss * Ca_ss) * k1))
                    + FP_LITERAL(2.) * Ca_ss * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);

            const cellmodel_float_t di_CaL_factors_dCa_ss =
                    F * g_CaL * d * Exp(FP_LITERAL(2.) * F * V_eff / (R * T)) * f * f2 * fCass;
            const cellmodel_float_t di_rel_dCa_ss =
                    -V_rel * O + V_rel * (-Ca_ss + Ca_SR) * dO_dCa_ss;
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
                               ? (Expm1(dt * dCa_ss_dt_linearized)) * dCa_ss_dt
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
                                    + FP_LITERAL(0.1245)
                                              * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));
            const cellmodel_float_t dNa_i_dt_linearized =
                    Cm
                    * (FP_LITERAL(-3.) * di_NaCa_dNa_i - FP_LITERAL(3.) * di_NaK_dNa_i
                       + g_bna * dE_Na_dNa_i - dE_Na_dNa_i * di_Na_dE_Na)
                    / (F * V_c);
            states[STATE_Na_i * padded_num_cells + i] =
                    Na_i
                    + (fabs(dNa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (Expm1(dt * dNa_i_dt_linearized)) * dNa_i_dt / dNa_i_dt_linearized
                               : dt * dNa_i_dt);

            // Expressions for the Membrane component
            const cellmodel_float_t i_Stim = (is_stimulated ? -stim_amplitude : FP_LITERAL(0.));
            const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK
                                            - i_Stim - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
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
                                 * (FP_LITERAL(1.)
                                    + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)));
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
                       + FP_LITERAL(0.01245) * F * Exp(FP_LITERAL(-0.1) * F * V / (R * T))
                                 / (R * T))
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
                            * (Ca_o * F * gamma * (Na_i * Na_i * Na_i)
                                       * Exp(F * gamma * V / (R * T)) / (R * T)
                               - F * alpha * (Na_o * Na_o * Na_o) * (FP_LITERAL(-1.) + gamma) * Ca_i
                                         * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T))
                                         / (R * T))
                            / ((FP_LITERAL(1.)
                                + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                               * (Ca_o + Km_Ca)
                               * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)))
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
                              * (dalpha_K1_dV * dxK1_inf_dalpha_K1
                                 + dbeta_K1_dV * dxK1_inf_dbeta_K1);
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
                    -g_bca - g_bna - di_K1_dV - di_Kr_dV - di_Ks_dV - di_NaCa_dV - di_NaK_dV
                    - di_Na_dV - di_p_K_dV - di_to_dV
                    - (dalpha_K1_dV * dxK1_inf_dalpha_K1 + dbeta_K1_dV * dxK1_inf_dbeta_K1)
                              * di_K1_dxK1_inf
                    - di_CaL_factors_dV_eff * i_CaL_fraction
                    - di_CaL_fraction_dV_eff * i_CaL_factors;
            states[STATE_V * padded_num_cells + i] =
                    (fabs(dV_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (Expm1(dt * dV_dt_linearized)) * dV_dt / dV_dt_linearized
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
                                 * (FP_LITERAL(1.)
                                    + Exp(FP_LITERAL(0.5) * E_K - FP_LITERAL(0.5) * V)));
            const cellmodel_float_t dalpha_K1_dE_K =
                    FP_LITERAL(3.68652741199693e-8)
                    * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K)
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
            const cellmodel_float_t di_p_K_dE_K =
                    -g_pK
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
                               ? (Expm1(dt * dK_i_dt_linearized)) * dK_i_dt / dK_i_dt_linearized
                               : dt * dK_i_dt);
        }
    }
}

// Compute a forward step using the Rush-Larsen scheme to the TP06 ODE
static void step_RL(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
                    const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
                    const long num_cells, long padded_num_cells)
{
    #pragma omp parallel
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

#if defined(HINT_CLANG_SIMD)
#pragma omp for
#pragma clang loop vectorize(assume_safety)
#elif defined(HINT_OMP_SIMD)
#ifdef VECTOR_LENGTH
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES) simdlen(VECTOR_LENGTH)
#else
#pragma omp for simd aligned(states : CELLMODEL_STATES_ALIGNMENT_BYTES)
#endif // defined(VECTOR_LENGTH)
#else
#pragma omp for
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

            // Expressions for the Reversal potentials component
            const cellmodel_float_t E_Na = R * T * Log(Na_o / Na_i) / F;
            const cellmodel_float_t E_K = R * T * Log(K_o / K_i) / F;
            const cellmodel_float_t E_Ks =
                    R * T * Log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F;
            const cellmodel_float_t E_Ca = FP_LITERAL(0.5) * R * T * Log(Ca_o / Ca_i) / F;

            // Expressions for the Inward rectifier potassium current component
            const cellmodel_float_t alpha_K1 =
                    FP_LITERAL(0.1)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(6.14421235332821e-6)
                                 * Exp(FP_LITERAL(0.06) * V - FP_LITERAL(0.06) * E_K));
            const cellmodel_float_t beta_K1 =
                    (FP_LITERAL(0.367879441171442)
                             * Exp(FP_LITERAL(0.1) * V - FP_LITERAL(0.1) * E_K)
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
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-26.) / FP_LITERAL(7.) - V / FP_LITERAL(7.)));
            const cellmodel_float_t alpha_xr1 =
                    FP_LITERAL(450.)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-9.) / FP_LITERAL(2.) - V / FP_LITERAL(10.)));
            const cellmodel_float_t beta_xr1 = FP_LITERAL(6.)
                                               / (FP_LITERAL(1.)
                                                  + Exp(FP_LITERAL(60.) / FP_LITERAL(23.)
                                                        + FP_LITERAL(2.) * V / FP_LITERAL(23.)));
            const cellmodel_float_t tau_xr1 = alpha_xr1 * beta_xr1;
            const cellmodel_float_t dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1;
            const cellmodel_float_t dXr1_dt_linearized = FP_LITERAL(-1.) / tau_xr1;
            states[STATE_Xr1 * padded_num_cells + i] =
                    (Expm1(dt * dXr1_dt_linearized)) * dXr1_dt / dXr1_dt_linearized + Xr1;

            // Expressions for the Xr2 gate component
            const cellmodel_float_t xr2_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(11.) / FP_LITERAL(3.) + V / FP_LITERAL(24.)));
            const cellmodel_float_t alpha_xr2 =
                    FP_LITERAL(3.) / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) - V / FP_LITERAL(20.)));
            const cellmodel_float_t beta_xr2 =
                    FP_LITERAL(1.12)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(-3.) + V / FP_LITERAL(20.)));
            const cellmodel_float_t tau_xr2 = alpha_xr2 * beta_xr2;
            const cellmodel_float_t dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2;
            const cellmodel_float_t dXr2_dt_linearized = FP_LITERAL(-1.) / tau_xr2;
            states[STATE_Xr2 * padded_num_cells + i] =
                    (Expm1(dt * dXr2_dt_linearized)) * dXr2_dt / dXr2_dt_linearized + Xr2;

            // Expressions for the Slow time dependent potassium current component
            const cellmodel_float_t i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V);

            // Expressions for the Xs gate component
            const cellmodel_float_t xs_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(14.) - V / FP_LITERAL(14.)));
            const cellmodel_float_t alpha_xs =
                    FP_LITERAL(1400.)
                    / sqrt(FP_LITERAL(1.)
                           + Exp(FP_LITERAL(5.) / FP_LITERAL(6.) - V / FP_LITERAL(6.)));
            const cellmodel_float_t beta_xs =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-7.) / FP_LITERAL(3.) + V / FP_LITERAL(15.)));
            const cellmodel_float_t tau_xs = FP_LITERAL(80.) + alpha_xs * beta_xs;
            const cellmodel_float_t dXs_dt = (-Xs + xs_inf) / tau_xs;
            const cellmodel_float_t dXs_dt_linearized = FP_LITERAL(-1.) / tau_xs;
            states[STATE_Xs * padded_num_cells + i] =
                    (Expm1(dt * dXs_dt_linearized)) * dXs_dt / dXs_dt_linearized + Xs;

            // Expressions for the Fast sodium current component
            const cellmodel_float_t i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j;

            // Expressions for the m gate component
            const cellmodel_float_t m_inf =
                    FP_LITERAL(1.0)
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
                    (Expm1(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized + m;

            // Expressions for the h gate component
            const cellmodel_float_t h_inf =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                              + FP_LITERAL(100.) * V / FP_LITERAL(743.)))
                       * (1.
                          + Exp(FP_LITERAL(7155.) / FP_LITERAL(743.)
                                + FP_LITERAL(100.) * V / FP_LITERAL(743.))));
            const cellmodel_float_t alpha_h =
                    (V < FP_LITERAL(-40.) ? FP_LITERAL(4.43126792958051e-7)
                                                    * Exp(FP_LITERAL(-0.147058823529412) * V)
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
                    (Expm1(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized + h;

            // Expressions for the j gate component
            const cellmodel_float_t j_inf =
                    FP_LITERAL(1.0)
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
                                          + FP_LITERAL(50262745825.954)
                                                    * Exp(FP_LITERAL(0.311) * V))
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
                    (Expm1(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized + j;

            // Expressions for the Sodium background current component
            const cellmodel_float_t i_b_Na = g_bna * (-E_Na + V);

            // Expressions for the L_type Ca current component
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
            const cellmodel_float_t dd_dt_linearized = FP_LITERAL(-1.) / tau_d;
            states[STATE_d * padded_num_cells + i] =
                    (Expm1(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized + d;

            // Expressions for the f gate component
            const cellmodel_float_t f_inf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(20.) / FP_LITERAL(7.) + V / FP_LITERAL(7.)));
            const cellmodel_float_t tau_f =
                    FP_LITERAL(20.)
                    + FP_LITERAL(180.)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(3.) + V / FP_LITERAL(10.)))
                    + FP_LITERAL(200.)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(13.) / FP_LITERAL(10.) - V / FP_LITERAL(10.)))
                    + FP_LITERAL(1102.5)
                              * Exp(-((FP_LITERAL(27.) + V) * (FP_LITERAL(27.) + V))
                                    / FP_LITERAL(225.));
            const cellmodel_float_t df_dt = (-f + f_inf) / tau_f;
            const cellmodel_float_t df_dt_linearized = FP_LITERAL(-1.) / tau_f;
            states[STATE_f * padded_num_cells + i] =
                    (Expm1(dt * df_dt_linearized)) * df_dt / df_dt_linearized + f;

            // Expressions for the F2 gate component
            const cellmodel_float_t f2_inf =
                    FP_LITERAL(0.33)
                    + FP_LITERAL(0.67)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(5.) + V / FP_LITERAL(7.)));
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
                    (Expm1(dt * df2_dt_linearized)) * df2_dt / df2_dt_linearized + f2;

            // Expressions for the FCass gate component
            const cellmodel_float_t fCass_inf =
                    FP_LITERAL(0.4)
                    + FP_LITERAL(0.6) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
            const cellmodel_float_t tau_fCass =
                    FP_LITERAL(2.)
                    + FP_LITERAL(80.) / (FP_LITERAL(1.) + FP_LITERAL(400.) * (Ca_ss * Ca_ss));
            const cellmodel_float_t dfCass_dt = (-fCass + fCass_inf) / tau_fCass;
            const cellmodel_float_t dfCass_dt_linearized = FP_LITERAL(-1.) / tau_fCass;
            states[STATE_fCass * padded_num_cells + i] =
                    (Expm1(dt * dfCass_dt_linearized)) * dfCass_dt / dfCass_dt_linearized + fCass;

            // Expressions for the Calcium background current component
            const cellmodel_float_t i_b_Ca = g_bca * (-E_Ca + V);

            // Expressions for the Transient outward current component
            const cellmodel_float_t i_to = g_to * (-E_K + V) * r * s;

            // Expressions for the s gate component
            const cellmodel_float_t s_inf =
                    (celltype == FP_LITERAL(2.)
                             ? FP_LITERAL(1.0)
                                       / (FP_LITERAL(1.)
                                          + Exp(FP_LITERAL(28.) / FP_LITERAL(5.)
                                                + V / FP_LITERAL(5.)))
                             : FP_LITERAL(1.0)
                                       / (FP_LITERAL(1.)
                                          + Exp(FP_LITERAL(4.) + V / FP_LITERAL(5.))));
            const cellmodel_float_t tau_s =
                    (celltype == FP_LITERAL(2.)
                             ? FP_LITERAL(8.)
                                       + FP_LITERAL(1000.)
                                                 * Exp(-((FP_LITERAL(67.) + V)
                                                         * (FP_LITERAL(67.) + V))
                                                       / FP_LITERAL(1000.))
                             : FP_LITERAL(3.)
                                       + FP_LITERAL(5.)
                                                 / (FP_LITERAL(1.)
                                                    + Exp(FP_LITERAL(-4.) + V / FP_LITERAL(5.)))
                                       + FP_LITERAL(85.)
                                                 * Exp(-((FP_LITERAL(45.) + V)
                                                         * (FP_LITERAL(45.) + V))
                                                       / FP_LITERAL(320.)));
            const cellmodel_float_t ds_dt = (-s + s_inf) / tau_s;
            const cellmodel_float_t ds_dt_linearized = FP_LITERAL(-1.) / tau_s;
            states[STATE_s * padded_num_cells + i] =
                    (Expm1(dt * ds_dt_linearized)) * ds_dt / ds_dt_linearized + s;

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
            const cellmodel_float_t dr_dt_linearized = FP_LITERAL(-1.) / tau_r;
            states[STATE_r * padded_num_cells + i] =
                    (Expm1(dt * dr_dt_linearized)) * dr_dt / dr_dt_linearized + r;

            // Expressions for the Sodium potassium pump current component
            const cellmodel_float_t i_NaK =
                    K_o * P_NaK * Na_i
                    / ((K_mNa + Na_i) * (K_mk + K_o)
                       * (FP_LITERAL(1.) + FP_LITERAL(0.0353) * Exp(-F * V / (R * T))
                          + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * F * V / (R * T))));

            // Expressions for the Sodium calcium exchanger current component
            const cellmodel_float_t i_NaCa =
                    K_NaCa
                    * (Ca_o * (Na_i * Na_i * Na_i) * Exp(F * gamma * V / (R * T))
                       - alpha * (Na_o * Na_o * Na_o) * Ca_i
                                 * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                    / ((FP_LITERAL(1.) + K_sat * Exp(F * (FP_LITERAL(-1.) + gamma) * V / (R * T)))
                       * (Ca_o + Km_Ca) * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o)));

            // Expressions for the Calcium pump current component
            const cellmodel_float_t i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i);

            // Expressions for the Potassium pump current component
            const cellmodel_float_t i_p_K = g_pK * (-E_K + V)
                                            / (FP_LITERAL(1.)
                                               + Exp(FP_LITERAL(1250.) / FP_LITERAL(299.)
                                                     - FP_LITERAL(50.) * V / FP_LITERAL(299.)));

            // Expressions for the Calcium dynamics component
            const cellmodel_float_t i_up =
                    Vmax_up / (FP_LITERAL(1.) + (K_up * K_up) / (Ca_i * Ca_i));
            const cellmodel_float_t i_leak = V_leak * (-Ca_i + Ca_SR);
            const cellmodel_float_t i_xfer = V_xfer * (-Ca_i + Ca_ss);
            const cellmodel_float_t kcasr =
                    max_sr - (max_sr - min_sr) / (FP_LITERAL(1.) + (EC * EC) / (Ca_SR * Ca_SR));
            const cellmodel_float_t Ca_i_bufc =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)));
            const cellmodel_float_t Ca_sr_bufsr =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)));
            const cellmodel_float_t Ca_ss_bufss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)));
            const cellmodel_float_t dCa_i_dt = (V_sr * (-i_up + i_leak) / V_c
                                                - Cm * (FP_LITERAL(-2.) * i_NaCa + i_b_Ca + i_p_Ca)
                                                          / (FP_LITERAL(2.) * F * V_c)
                                                + i_xfer)
                                               * Ca_i_bufc;
            states[STATE_Ca_i * padded_num_cells + i] = dt * dCa_i_dt + Ca_i;
            const cellmodel_float_t k1 = k1_prime / kcasr;
            const cellmodel_float_t k2 = k2_prime * kcasr;
            const cellmodel_float_t O =
                    (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1);
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
            const cellmodel_float_t dV_dt = -i_CaL - i_K1 - i_Kr - i_Ks - i_Na - i_NaCa - i_NaK
                                            - i_Stim - i_b_Ca - i_b_Na - i_p_Ca - i_p_K - i_to;
            states[STATE_V * padded_num_cells + i] = dt * dV_dt + V;

            // Expressions for the Potassium dynamics component
            const cellmodel_float_t dK_i_dt =
                    Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + FP_LITERAL(2.) * i_NaK)
                    / (F * V_c);
            states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;
        }
    }
}

const char *state_names[] = {
    "Xr1",
    "Xr2",
    "Xs",
    "m",
    "h",
    "j",
    "d",
    "f",
    "f2",
    "fCass",
    "s",
    "r",
    "Ca_i",
    "R_prime",
    "Ca_SR",
    "Ca_ss",
    "Na_i",
    "V",
    "K_i",
    NULL
};

// clang-format off
const struct cellmodel model_TP06_simd = {
        .init_states = &init_state_values,
        .init_parameters = &init_parameters_values,
        .state_index = &state_index,
        .parameter_index = &parameter_index,
        .step_FE = &step_FE,
        .step_RL = &step_RL,
        .step_GRL1 = &step_GRL1,
        .num_states = NUM_STATES,
        .num_parameters = NUM_PARAMS,
        .layout = LAYOUT_STRUCT_OF_ARRAYS,
        .colour_sets = NULL,
        .state_names = state_names,
};
// clang-format on
