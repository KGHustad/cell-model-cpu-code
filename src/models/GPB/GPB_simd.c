#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "GPB" model

#include "GPB_simd.h"

enum state {
    STATE_m,
    STATE_h,
    STATE_j,
    STATE_x_kr,
    STATE_x_ks,
    STATE_x_to_s,
    STATE_y_to_s,
    STATE_x_to_f,
    STATE_y_to_f,
    STATE_d,
    STATE_f,
    STATE_f_Ca_Bj,
    STATE_f_Ca_Bsl,
    STATE_Ry_Rr,
    STATE_Ry_Ro,
    STATE_Ry_Ri,
    STATE_Na_Bj,
    STATE_Na_Bsl,
    STATE_Tn_CL,
    STATE_Tn_CHc,
    STATE_Tn_CHm,
    STATE_CaM,
    STATE_Myo_c,
    STATE_Myo_m,
    STATE_SRB,
    STATE_SLL_j,
    STATE_SLL_sl,
    STATE_SLH_j,
    STATE_SLH_sl,
    STATE_Csqn_b,
    STATE_Ca_sr,
    STATE_Na_j,
    STATE_Na_sl,
    STATE_Na_i,
    STATE_K_i,
    STATE_Ca_j,
    STATE_Ca_sl,
    STATE_Ca_i,
    STATE_V_m,
    NUM_STATES,
};

enum parameter {
    PARAM_Fjunc,
    PARAM_Fjunc_CaL,
    PARAM_cellLength,
    PARAM_cellRadius,
    PARAM_distJuncSL,
    PARAM_distSLcyto,
    PARAM_junctionLength,
    PARAM_junctionRadius,
    PARAM_GNa,
    PARAM_GNaB,
    PARAM_IbarNaK,
    PARAM_KmKo,
    PARAM_KmNaip,
    PARAM_Q10KmNai,
    PARAM_Q10NaK,
    PARAM_GKr,
    PARAM_GKp,
    PARAM_GKs,
    PARAM_pNaK,
    PARAM_GK1,
    PARAM_Gto,
    PARAM_epi,
    PARAM_GClB,
    PARAM_GClCa,
    PARAM_KdClCa,
    PARAM_GCaL,
    PARAM_Q10CaL,
    PARAM_pCa,
    PARAM_pK,
    PARAM_pNa,
    PARAM_IbarNCX,
    PARAM_Kdact,
    PARAM_KmCai,
    PARAM_KmCao,
    PARAM_KmNai,
    PARAM_KmNao,
    PARAM_Q10NCX,
    PARAM_ksat,
    PARAM_nu,
    PARAM_IbarSLCaP,
    PARAM_KmPCa,
    PARAM_Q10SLCaP,
    PARAM_GCaB,
    PARAM_Kmf,
    PARAM_Kmr,
    PARAM_MaxSR,
    PARAM_MinSR,
    PARAM_Q10SRCaP,
    PARAM_Vmax_SRCaP,
    PARAM_ec50SR,
    PARAM_hillSRCaP,
    PARAM_kiCa,
    PARAM_kim,
    PARAM_koCa,
    PARAM_kom,
    PARAM_ks,
    PARAM_Bmax_Naj,
    PARAM_Bmax_Nasl,
    PARAM_koff_na,
    PARAM_kon_na,
    PARAM_Bmax_CaM,
    PARAM_Bmax_SR,
    PARAM_Bmax_TnChigh,
    PARAM_Bmax_TnClow,
    PARAM_Bmax_myosin,
    PARAM_koff_cam,
    PARAM_koff_myoca,
    PARAM_koff_myomg,
    PARAM_koff_sr,
    PARAM_koff_tnchca,
    PARAM_koff_tnchmg,
    PARAM_koff_tncl,
    PARAM_kon_cam,
    PARAM_kon_myoca,
    PARAM_kon_myomg,
    PARAM_kon_sr,
    PARAM_kon_tnchca,
    PARAM_kon_tnchmg,
    PARAM_kon_tncl,
    PARAM_Bmax_SLhighj0,
    PARAM_Bmax_SLhighsl0,
    PARAM_Bmax_SLlowj0,
    PARAM_Bmax_SLlowsl0,
    PARAM_koff_slh,
    PARAM_koff_sll,
    PARAM_kon_slh,
    PARAM_kon_sll,
    PARAM_Bmax_Csqn0,
    PARAM_DcaJuncSL,
    PARAM_DcaSLcyto,
    PARAM_J_ca_juncsl,
    PARAM_J_ca_slmyo,
    PARAM_koff_csqn,
    PARAM_kon_csqn,
    PARAM_DnaJuncSL,
    PARAM_DnaSLcyto,
    PARAM_J_na_juncsl,
    PARAM_J_na_slmyo,
    PARAM_Nao,
    PARAM_Ko,
    PARAM_Cao,
    PARAM_Cli,
    PARAM_Clo,
    PARAM_Mgi,
    PARAM_Cmem,
    PARAM_Frdy,
    PARAM_R,
    PARAM_Temp,
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
        states[STATE_m * padded_num_cells + i] = FP_LITERAL(0.003793087414436);
        states[STATE_h * padded_num_cells + i] = FP_LITERAL(0.626221949492493);
        states[STATE_j * padded_num_cells + i] = FP_LITERAL(0.624553572490432);
        states[STATE_x_kr * padded_num_cells + i] = FP_LITERAL(0.0210022533039071);
        states[STATE_x_ks * padded_num_cells + i] = FP_LITERAL(0.00428016666258923);
        states[STATE_x_to_s * padded_num_cells + i] = FP_LITERAL(0.000440445885642567);
        states[STATE_y_to_s * padded_num_cells + i] = FP_LITERAL(0.785115828275182);
        states[STATE_x_to_f * padded_num_cells + i] = FP_LITERAL(0.000440438103758954);
        states[STATE_y_to_f * padded_num_cells + i] = FP_LITERAL(0.999995844038706);
        states[STATE_d * padded_num_cells + i] = FP_LITERAL(2.92407183949469e-06);
        states[STATE_f * padded_num_cells + i] = FP_LITERAL(0.995135796703515);
        states[STATE_f_Ca_Bj * padded_num_cells + i] = FP_LITERAL(0.0246760872105795);
        states[STATE_f_Ca_Bsl * padded_num_cells + i] = FP_LITERAL(0.0152723084239416);
        states[STATE_Ry_Rr * padded_num_cells + i] = FP_LITERAL(0.890806040818203);
        states[STATE_Ry_Ro * padded_num_cells + i] = FP_LITERAL(7.40481128853622e-07);
        states[STATE_Ry_Ri * padded_num_cells + i] = FP_LITERAL(9.07666168960848e-08);
        states[STATE_Na_Bj * padded_num_cells + i] = FP_LITERAL(3.4543773303328);
        states[STATE_Na_Bsl * padded_num_cells + i] = FP_LITERAL(0.753740951477775);
        states[STATE_Tn_CL * padded_num_cells + i] = FP_LITERAL(0.00893455096919132);
        states[STATE_Tn_CHc * padded_num_cells + i] = FP_LITERAL(0.117412025936615);
        states[STATE_Tn_CHm * padded_num_cells + i] = FP_LITERAL(0.0106160166692932);
        states[STATE_CaM * padded_num_cells + i] = FP_LITERAL(0.000295573424135051);
        states[STATE_Myo_c * padded_num_cells + i] = FP_LITERAL(0.00192322252438022);
        states[STATE_Myo_m * padded_num_cells + i] = FP_LITERAL(0.137560495022823);
        states[STATE_SRB * padded_num_cells + i] = FP_LITERAL(0.00217360235649355);
        states[STATE_SLL_j * padded_num_cells + i] = FP_LITERAL(0.00740524521680039);
        states[STATE_SLL_sl * padded_num_cells + i] = FP_LITERAL(0.00990339304377132);
        states[STATE_SLH_j * padded_num_cells + i] = FP_LITERAL(0.0735890020284214);
        states[STATE_SLH_sl * padded_num_cells + i] = FP_LITERAL(0.114583623436917);
        states[STATE_Csqn_b * padded_num_cells + i] = FP_LITERAL(1.19723145924432);
        states[STATE_Ca_sr * padded_num_cells + i] = FP_LITERAL(0.554760499828172);
        states[STATE_Na_j * padded_num_cells + i] = FP_LITERAL(8.40537012592918);
        states[STATE_Na_sl * padded_num_cells + i] = FP_LITERAL(8.40491910001025);
        states[STATE_Na_i * padded_num_cells + i] = FP_LITERAL(8.40513364344858);
        states[STATE_K_i * padded_num_cells + i] = FP_LITERAL(120.0);
        states[STATE_Ca_j * padded_num_cells + i] = FP_LITERAL(0.000175882395147342);
        states[STATE_Ca_sl * padded_num_cells + i] = FP_LITERAL(0.000106779509977354);
        states[STATE_Ca_i * padded_num_cells + i] = FP_LITERAL(8.72509677797499e-05);
        states[STATE_V_m * padded_num_cells + i] = FP_LITERAL(-81.4552030512661);
    }
}

// Default parameter values
static void init_parameters_values(cellmodel_float_t *__restrict parameters)
{
    parameters[PARAM_Fjunc] = FP_LITERAL(0.11);
    parameters[PARAM_Fjunc_CaL] = FP_LITERAL(0.9);
    parameters[PARAM_cellLength] = FP_LITERAL(100.0);
    parameters[PARAM_cellRadius] = FP_LITERAL(10.25);
    parameters[PARAM_distJuncSL] = FP_LITERAL(0.5);
    parameters[PARAM_distSLcyto] = FP_LITERAL(0.45);
    parameters[PARAM_junctionLength] = FP_LITERAL(0.16);
    parameters[PARAM_junctionRadius] = FP_LITERAL(0.015);
    parameters[PARAM_GNa] = FP_LITERAL(23.0);
    parameters[PARAM_GNaB] = FP_LITERAL(0.000597);
    parameters[PARAM_IbarNaK] = FP_LITERAL(1.8);
    parameters[PARAM_KmKo] = FP_LITERAL(1.5);
    parameters[PARAM_KmNaip] = FP_LITERAL(11.0);
    parameters[PARAM_Q10KmNai] = FP_LITERAL(1.39);
    parameters[PARAM_Q10NaK] = FP_LITERAL(1.63);
    parameters[PARAM_GKr] = FP_LITERAL(0.035);
    parameters[PARAM_GKp] = FP_LITERAL(0.002);
    parameters[PARAM_GKs] = FP_LITERAL(0.0035);
    parameters[PARAM_pNaK] = FP_LITERAL(0.01833);
    parameters[PARAM_GK1] = FP_LITERAL(0.35);
    parameters[PARAM_Gto] = FP_LITERAL(0.13);
    parameters[PARAM_epi] = FP_LITERAL(1.0);
    parameters[PARAM_GClB] = FP_LITERAL(0.009);
    parameters[PARAM_GClCa] = FP_LITERAL(0.0548125);
    parameters[PARAM_KdClCa] = FP_LITERAL(0.1);
    parameters[PARAM_GCaL] = FP_LITERAL(0.5);
    parameters[PARAM_Q10CaL] = FP_LITERAL(1.8);
    parameters[PARAM_pCa] = FP_LITERAL(0.00054);
    parameters[PARAM_pK] = FP_LITERAL(2.7e-07);
    parameters[PARAM_pNa] = FP_LITERAL(1.5e-08);
    parameters[PARAM_IbarNCX] = FP_LITERAL(4.5);
    parameters[PARAM_Kdact] = FP_LITERAL(0.00015);
    parameters[PARAM_KmCai] = FP_LITERAL(0.00359);
    parameters[PARAM_KmCao] = FP_LITERAL(1.3);
    parameters[PARAM_KmNai] = FP_LITERAL(12.29);
    parameters[PARAM_KmNao] = FP_LITERAL(87.5);
    parameters[PARAM_Q10NCX] = FP_LITERAL(1.57);
    parameters[PARAM_ksat] = FP_LITERAL(0.32);
    parameters[PARAM_nu] = FP_LITERAL(0.27);
    parameters[PARAM_IbarSLCaP] = FP_LITERAL(0.0673);
    parameters[PARAM_KmPCa] = FP_LITERAL(0.0005);
    parameters[PARAM_Q10SLCaP] = FP_LITERAL(2.35);
    parameters[PARAM_GCaB] = FP_LITERAL(0.0005513);
    parameters[PARAM_Kmf] = FP_LITERAL(0.000246);
    parameters[PARAM_Kmr] = FP_LITERAL(1.7);
    parameters[PARAM_MaxSR] = FP_LITERAL(15.0);
    parameters[PARAM_MinSR] = FP_LITERAL(1.0);
    parameters[PARAM_Q10SRCaP] = FP_LITERAL(2.6);
    parameters[PARAM_Vmax_SRCaP] = FP_LITERAL(0.0053114);
    parameters[PARAM_ec50SR] = FP_LITERAL(0.45);
    parameters[PARAM_hillSRCaP] = FP_LITERAL(1.787);
    parameters[PARAM_kiCa] = FP_LITERAL(0.5);
    parameters[PARAM_kim] = FP_LITERAL(0.005);
    parameters[PARAM_koCa] = FP_LITERAL(10.0);
    parameters[PARAM_kom] = FP_LITERAL(0.06);
    parameters[PARAM_ks] = FP_LITERAL(25.0);
    parameters[PARAM_Bmax_Naj] = FP_LITERAL(7.561);
    parameters[PARAM_Bmax_Nasl] = FP_LITERAL(1.65);
    parameters[PARAM_koff_na] = FP_LITERAL(0.001);
    parameters[PARAM_kon_na] = FP_LITERAL(0.0001);
    parameters[PARAM_Bmax_CaM] = FP_LITERAL(0.024);
    parameters[PARAM_Bmax_SR] = FP_LITERAL(0.0171);
    parameters[PARAM_Bmax_TnChigh] = FP_LITERAL(0.14);
    parameters[PARAM_Bmax_TnClow] = FP_LITERAL(0.07);
    parameters[PARAM_Bmax_myosin] = FP_LITERAL(0.14);
    parameters[PARAM_koff_cam] = FP_LITERAL(0.238);
    parameters[PARAM_koff_myoca] = FP_LITERAL(0.00046);
    parameters[PARAM_koff_myomg] = FP_LITERAL(5.7e-05);
    parameters[PARAM_koff_sr] = FP_LITERAL(0.06);
    parameters[PARAM_koff_tnchca] = FP_LITERAL(3.2e-05);
    parameters[PARAM_koff_tnchmg] = FP_LITERAL(0.00333);
    parameters[PARAM_koff_tncl] = FP_LITERAL(0.0196);
    parameters[PARAM_kon_cam] = FP_LITERAL(34.0);
    parameters[PARAM_kon_myoca] = FP_LITERAL(13.8);
    parameters[PARAM_kon_myomg] = FP_LITERAL(0.0157);
    parameters[PARAM_kon_sr] = FP_LITERAL(100.0);
    parameters[PARAM_kon_tnchca] = FP_LITERAL(2.37);
    parameters[PARAM_kon_tnchmg] = FP_LITERAL(0.003);
    parameters[PARAM_kon_tncl] = FP_LITERAL(32.7);
    parameters[PARAM_Bmax_SLhighj0] = FP_LITERAL(0.000165);
    parameters[PARAM_Bmax_SLhighsl0] = FP_LITERAL(0.0134);
    parameters[PARAM_Bmax_SLlowj0] = FP_LITERAL(0.00046);
    parameters[PARAM_Bmax_SLlowsl0] = FP_LITERAL(0.0374);
    parameters[PARAM_koff_slh] = FP_LITERAL(0.03);
    parameters[PARAM_koff_sll] = FP_LITERAL(1.3);
    parameters[PARAM_kon_slh] = FP_LITERAL(100.0);
    parameters[PARAM_kon_sll] = FP_LITERAL(100.0);
    parameters[PARAM_Bmax_Csqn0] = FP_LITERAL(0.14);
    parameters[PARAM_DcaJuncSL] = FP_LITERAL(1.64e-06);
    parameters[PARAM_DcaSLcyto] = FP_LITERAL(1.22e-06);
    parameters[PARAM_J_ca_juncsl] = FP_LITERAL(8.2413e-13);
    parameters[PARAM_J_ca_slmyo] = FP_LITERAL(3.7243e-12);
    parameters[PARAM_koff_csqn] = FP_LITERAL(65.0);
    parameters[PARAM_kon_csqn] = FP_LITERAL(100.0);
    parameters[PARAM_DnaJuncSL] = FP_LITERAL(1.09e-05);
    parameters[PARAM_DnaSLcyto] = FP_LITERAL(1.79e-05);
    parameters[PARAM_J_na_juncsl] = FP_LITERAL(1.8313e-14);
    parameters[PARAM_J_na_slmyo] = FP_LITERAL(1.6386e-12);
    parameters[PARAM_Nao] = FP_LITERAL(140.0);
    parameters[PARAM_Ko] = FP_LITERAL(5.4);
    parameters[PARAM_Cao] = FP_LITERAL(1.8);
    parameters[PARAM_Cli] = FP_LITERAL(15.0);
    parameters[PARAM_Clo] = FP_LITERAL(150.0);
    parameters[PARAM_Mgi] = FP_LITERAL(1.0);
    parameters[PARAM_Cmem] = FP_LITERAL(1.381e-10);
    parameters[PARAM_Frdy] = FP_LITERAL(96485.0);
    parameters[PARAM_R] = FP_LITERAL(8314.0);
    parameters[PARAM_Temp] = FP_LITERAL(310.0);
    parameters[PARAM_stim_amplitude] = FP_LITERAL(40.0);
    parameters[PARAM_stim_duration] = FP_LITERAL(1.0);
    parameters[PARAM_stim_period] = FP_LITERAL(1000.0);
    parameters[PARAM_stim_start] = FP_LITERAL(0.0);
}

// State index
static int state_index(const char name[])
{
    if (strcmp(name, "m") == 0) {
        return STATE_m;
    } else if (strcmp(name, "h") == 0) {
        return STATE_h;
    } else if (strcmp(name, "j") == 0) {
        return STATE_j;
    } else if (strcmp(name, "x_kr") == 0) {
        return STATE_x_kr;
    } else if (strcmp(name, "x_ks") == 0) {
        return STATE_x_ks;
    } else if (strcmp(name, "x_to_s") == 0) {
        return STATE_x_to_s;
    } else if (strcmp(name, "y_to_s") == 0) {
        return STATE_y_to_s;
    } else if (strcmp(name, "x_to_f") == 0) {
        return STATE_x_to_f;
    } else if (strcmp(name, "y_to_f") == 0) {
        return STATE_y_to_f;
    } else if (strcmp(name, "d") == 0) {
        return STATE_d;
    } else if (strcmp(name, "f") == 0) {
        return STATE_f;
    } else if (strcmp(name, "f_Ca_Bj") == 0) {
        return STATE_f_Ca_Bj;
    } else if (strcmp(name, "f_Ca_Bsl") == 0) {
        return STATE_f_Ca_Bsl;
    } else if (strcmp(name, "Ry_Rr") == 0) {
        return STATE_Ry_Rr;
    } else if (strcmp(name, "Ry_Ro") == 0) {
        return STATE_Ry_Ro;
    } else if (strcmp(name, "Ry_Ri") == 0) {
        return STATE_Ry_Ri;
    } else if (strcmp(name, "Na_Bj") == 0) {
        return STATE_Na_Bj;
    } else if (strcmp(name, "Na_Bsl") == 0) {
        return STATE_Na_Bsl;
    } else if (strcmp(name, "Tn_CL") == 0) {
        return STATE_Tn_CL;
    } else if (strcmp(name, "Tn_CHc") == 0) {
        return STATE_Tn_CHc;
    } else if (strcmp(name, "Tn_CHm") == 0) {
        return STATE_Tn_CHm;
    } else if (strcmp(name, "CaM") == 0) {
        return STATE_CaM;
    } else if (strcmp(name, "Myo_c") == 0) {
        return STATE_Myo_c;
    } else if (strcmp(name, "Myo_m") == 0) {
        return STATE_Myo_m;
    } else if (strcmp(name, "SRB") == 0) {
        return STATE_SRB;
    } else if (strcmp(name, "SLL_j") == 0) {
        return STATE_SLL_j;
    } else if (strcmp(name, "SLL_sl") == 0) {
        return STATE_SLL_sl;
    } else if (strcmp(name, "SLH_j") == 0) {
        return STATE_SLH_j;
    } else if (strcmp(name, "SLH_sl") == 0) {
        return STATE_SLH_sl;
    } else if (strcmp(name, "Csqn_b") == 0) {
        return STATE_Csqn_b;
    } else if (strcmp(name, "Ca_sr") == 0) {
        return STATE_Ca_sr;
    } else if (strcmp(name, "Na_j") == 0) {
        return STATE_Na_j;
    } else if (strcmp(name, "Na_sl") == 0) {
        return STATE_Na_sl;
    } else if (strcmp(name, "Na_i") == 0) {
        return STATE_Na_i;
    } else if (strcmp(name, "K_i") == 0) {
        return STATE_K_i;
    } else if (strcmp(name, "Ca_j") == 0) {
        return STATE_Ca_j;
    } else if (strcmp(name, "Ca_sl") == 0) {
        return STATE_Ca_sl;
    } else if (strcmp(name, "Ca_i") == 0) {
        return STATE_Ca_i;
    } else if (strcmp(name, "V_m") == 0) {
        return STATE_V_m;
    }
    return -1;
}

// Parameter index
static int parameter_index(const char name[])
{
    if (strcmp(name, "Fjunc") == 0) {
        return PARAM_Fjunc;
    } else if (strcmp(name, "Fjunc_CaL") == 0) {
        return PARAM_Fjunc_CaL;
    } else if (strcmp(name, "cellLength") == 0) {
        return PARAM_cellLength;
    } else if (strcmp(name, "cellRadius") == 0) {
        return PARAM_cellRadius;
    } else if (strcmp(name, "distJuncSL") == 0) {
        return PARAM_distJuncSL;
    } else if (strcmp(name, "distSLcyto") == 0) {
        return PARAM_distSLcyto;
    } else if (strcmp(name, "junctionLength") == 0) {
        return PARAM_junctionLength;
    } else if (strcmp(name, "junctionRadius") == 0) {
        return PARAM_junctionRadius;
    } else if (strcmp(name, "GNa") == 0) {
        return PARAM_GNa;
    } else if (strcmp(name, "GNaB") == 0) {
        return PARAM_GNaB;
    } else if (strcmp(name, "IbarNaK") == 0) {
        return PARAM_IbarNaK;
    } else if (strcmp(name, "KmKo") == 0) {
        return PARAM_KmKo;
    } else if (strcmp(name, "KmNaip") == 0) {
        return PARAM_KmNaip;
    } else if (strcmp(name, "Q10KmNai") == 0) {
        return PARAM_Q10KmNai;
    } else if (strcmp(name, "Q10NaK") == 0) {
        return PARAM_Q10NaK;
    } else if (strcmp(name, "GKr") == 0) {
        return PARAM_GKr;
    } else if (strcmp(name, "GKp") == 0) {
        return PARAM_GKp;
    } else if (strcmp(name, "GKs") == 0) {
        return PARAM_GKs;
    } else if (strcmp(name, "pNaK") == 0) {
        return PARAM_pNaK;
    } else if (strcmp(name, "GK1") == 0) {
        return PARAM_GK1;
    } else if (strcmp(name, "Gto") == 0) {
        return PARAM_Gto;
    } else if (strcmp(name, "epi") == 0) {
        return PARAM_epi;
    } else if (strcmp(name, "GClB") == 0) {
        return PARAM_GClB;
    } else if (strcmp(name, "GClCa") == 0) {
        return PARAM_GClCa;
    } else if (strcmp(name, "KdClCa") == 0) {
        return PARAM_KdClCa;
    } else if (strcmp(name, "GCaL") == 0) {
        return PARAM_GCaL;
    } else if (strcmp(name, "Q10CaL") == 0) {
        return PARAM_Q10CaL;
    } else if (strcmp(name, "pCa") == 0) {
        return PARAM_pCa;
    } else if (strcmp(name, "pK") == 0) {
        return PARAM_pK;
    } else if (strcmp(name, "pNa") == 0) {
        return PARAM_pNa;
    } else if (strcmp(name, "IbarNCX") == 0) {
        return PARAM_IbarNCX;
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
    } else if (strcmp(name, "ksat") == 0) {
        return PARAM_ksat;
    } else if (strcmp(name, "nu") == 0) {
        return PARAM_nu;
    } else if (strcmp(name, "IbarSLCaP") == 0) {
        return PARAM_IbarSLCaP;
    } else if (strcmp(name, "KmPCa") == 0) {
        return PARAM_KmPCa;
    } else if (strcmp(name, "Q10SLCaP") == 0) {
        return PARAM_Q10SLCaP;
    } else if (strcmp(name, "GCaB") == 0) {
        return PARAM_GCaB;
    } else if (strcmp(name, "Kmf") == 0) {
        return PARAM_Kmf;
    } else if (strcmp(name, "Kmr") == 0) {
        return PARAM_Kmr;
    } else if (strcmp(name, "MaxSR") == 0) {
        return PARAM_MaxSR;
    } else if (strcmp(name, "MinSR") == 0) {
        return PARAM_MinSR;
    } else if (strcmp(name, "Q10SRCaP") == 0) {
        return PARAM_Q10SRCaP;
    } else if (strcmp(name, "Vmax_SRCaP") == 0) {
        return PARAM_Vmax_SRCaP;
    } else if (strcmp(name, "ec50SR") == 0) {
        return PARAM_ec50SR;
    } else if (strcmp(name, "hillSRCaP") == 0) {
        return PARAM_hillSRCaP;
    } else if (strcmp(name, "kiCa") == 0) {
        return PARAM_kiCa;
    } else if (strcmp(name, "kim") == 0) {
        return PARAM_kim;
    } else if (strcmp(name, "koCa") == 0) {
        return PARAM_koCa;
    } else if (strcmp(name, "kom") == 0) {
        return PARAM_kom;
    } else if (strcmp(name, "ks") == 0) {
        return PARAM_ks;
    } else if (strcmp(name, "Bmax_Naj") == 0) {
        return PARAM_Bmax_Naj;
    } else if (strcmp(name, "Bmax_Nasl") == 0) {
        return PARAM_Bmax_Nasl;
    } else if (strcmp(name, "koff_na") == 0) {
        return PARAM_koff_na;
    } else if (strcmp(name, "kon_na") == 0) {
        return PARAM_kon_na;
    } else if (strcmp(name, "Bmax_CaM") == 0) {
        return PARAM_Bmax_CaM;
    } else if (strcmp(name, "Bmax_SR") == 0) {
        return PARAM_Bmax_SR;
    } else if (strcmp(name, "Bmax_TnChigh") == 0) {
        return PARAM_Bmax_TnChigh;
    } else if (strcmp(name, "Bmax_TnClow") == 0) {
        return PARAM_Bmax_TnClow;
    } else if (strcmp(name, "Bmax_myosin") == 0) {
        return PARAM_Bmax_myosin;
    } else if (strcmp(name, "koff_cam") == 0) {
        return PARAM_koff_cam;
    } else if (strcmp(name, "koff_myoca") == 0) {
        return PARAM_koff_myoca;
    } else if (strcmp(name, "koff_myomg") == 0) {
        return PARAM_koff_myomg;
    } else if (strcmp(name, "koff_sr") == 0) {
        return PARAM_koff_sr;
    } else if (strcmp(name, "koff_tnchca") == 0) {
        return PARAM_koff_tnchca;
    } else if (strcmp(name, "koff_tnchmg") == 0) {
        return PARAM_koff_tnchmg;
    } else if (strcmp(name, "koff_tncl") == 0) {
        return PARAM_koff_tncl;
    } else if (strcmp(name, "kon_cam") == 0) {
        return PARAM_kon_cam;
    } else if (strcmp(name, "kon_myoca") == 0) {
        return PARAM_kon_myoca;
    } else if (strcmp(name, "kon_myomg") == 0) {
        return PARAM_kon_myomg;
    } else if (strcmp(name, "kon_sr") == 0) {
        return PARAM_kon_sr;
    } else if (strcmp(name, "kon_tnchca") == 0) {
        return PARAM_kon_tnchca;
    } else if (strcmp(name, "kon_tnchmg") == 0) {
        return PARAM_kon_tnchmg;
    } else if (strcmp(name, "kon_tncl") == 0) {
        return PARAM_kon_tncl;
    } else if (strcmp(name, "Bmax_SLhighj0") == 0) {
        return PARAM_Bmax_SLhighj0;
    } else if (strcmp(name, "Bmax_SLhighsl0") == 0) {
        return PARAM_Bmax_SLhighsl0;
    } else if (strcmp(name, "Bmax_SLlowj0") == 0) {
        return PARAM_Bmax_SLlowj0;
    } else if (strcmp(name, "Bmax_SLlowsl0") == 0) {
        return PARAM_Bmax_SLlowsl0;
    } else if (strcmp(name, "koff_slh") == 0) {
        return PARAM_koff_slh;
    } else if (strcmp(name, "koff_sll") == 0) {
        return PARAM_koff_sll;
    } else if (strcmp(name, "kon_slh") == 0) {
        return PARAM_kon_slh;
    } else if (strcmp(name, "kon_sll") == 0) {
        return PARAM_kon_sll;
    } else if (strcmp(name, "Bmax_Csqn0") == 0) {
        return PARAM_Bmax_Csqn0;
    } else if (strcmp(name, "DcaJuncSL") == 0) {
        return PARAM_DcaJuncSL;
    } else if (strcmp(name, "DcaSLcyto") == 0) {
        return PARAM_DcaSLcyto;
    } else if (strcmp(name, "J_ca_juncsl") == 0) {
        return PARAM_J_ca_juncsl;
    } else if (strcmp(name, "J_ca_slmyo") == 0) {
        return PARAM_J_ca_slmyo;
    } else if (strcmp(name, "koff_csqn") == 0) {
        return PARAM_koff_csqn;
    } else if (strcmp(name, "kon_csqn") == 0) {
        return PARAM_kon_csqn;
    } else if (strcmp(name, "DnaJuncSL") == 0) {
        return PARAM_DnaJuncSL;
    } else if (strcmp(name, "DnaSLcyto") == 0) {
        return PARAM_DnaSLcyto;
    } else if (strcmp(name, "J_na_juncsl") == 0) {
        return PARAM_J_na_juncsl;
    } else if (strcmp(name, "J_na_slmyo") == 0) {
        return PARAM_J_na_slmyo;
    } else if (strcmp(name, "Nao") == 0) {
        return PARAM_Nao;
    } else if (strcmp(name, "Ko") == 0) {
        return PARAM_Ko;
    } else if (strcmp(name, "Cao") == 0) {
        return PARAM_Cao;
    } else if (strcmp(name, "Cli") == 0) {
        return PARAM_Cli;
    } else if (strcmp(name, "Clo") == 0) {
        return PARAM_Clo;
    } else if (strcmp(name, "Mgi") == 0) {
        return PARAM_Mgi;
    } else if (strcmp(name, "Cmem") == 0) {
        return PARAM_Cmem;
    } else if (strcmp(name, "Frdy") == 0) {
        return PARAM_Frdy;
    } else if (strcmp(name, "R") == 0) {
        return PARAM_R;
    } else if (strcmp(name, "Temp") == 0) {
        return PARAM_Temp;
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


// Compute a forward step using the explicit Euler algorithm to the GPB ODE
static void step_FE(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
                    const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
                    const long num_cells, long padded_num_cells)
{
    #pragma omp parallel
    {
        // Assign parameters
        const cellmodel_float_t Fjunc = parameters[PARAM_Fjunc];
        const cellmodel_float_t Fjunc_CaL = parameters[PARAM_Fjunc_CaL];
        const cellmodel_float_t cellLength = parameters[PARAM_cellLength];
        const cellmodel_float_t cellRadius = parameters[PARAM_cellRadius];
        const cellmodel_float_t GNa = parameters[PARAM_GNa];
        const cellmodel_float_t GNaB = parameters[PARAM_GNaB];
        const cellmodel_float_t IbarNaK = parameters[PARAM_IbarNaK];
        const cellmodel_float_t KmKo = parameters[PARAM_KmKo];
        const cellmodel_float_t KmNaip = parameters[PARAM_KmNaip];
        const cellmodel_float_t GKr = parameters[PARAM_GKr];
        const cellmodel_float_t GKp = parameters[PARAM_GKp];
        const cellmodel_float_t GKs = parameters[PARAM_GKs];
        const cellmodel_float_t pNaK = parameters[PARAM_pNaK];
        const cellmodel_float_t GK1 = parameters[PARAM_GK1];
        const cellmodel_float_t Gto = parameters[PARAM_Gto];
        const cellmodel_float_t epi = parameters[PARAM_epi];
        const cellmodel_float_t GClB = parameters[PARAM_GClB];
        const cellmodel_float_t GClCa = parameters[PARAM_GClCa];
        const cellmodel_float_t KdClCa = parameters[PARAM_KdClCa];
        const cellmodel_float_t GCaL = parameters[PARAM_GCaL];
        const cellmodel_float_t Q10CaL = parameters[PARAM_Q10CaL];
        const cellmodel_float_t pCa = parameters[PARAM_pCa];
        const cellmodel_float_t pK = parameters[PARAM_pK];
        const cellmodel_float_t pNa = parameters[PARAM_pNa];
        const cellmodel_float_t IbarNCX = parameters[PARAM_IbarNCX];
        const cellmodel_float_t Kdact = parameters[PARAM_Kdact];
        const cellmodel_float_t KmCai = parameters[PARAM_KmCai];
        const cellmodel_float_t KmCao = parameters[PARAM_KmCao];
        const cellmodel_float_t KmNai = parameters[PARAM_KmNai];
        const cellmodel_float_t KmNao = parameters[PARAM_KmNao];
        const cellmodel_float_t Q10NCX = parameters[PARAM_Q10NCX];
        const cellmodel_float_t ksat = parameters[PARAM_ksat];
        const cellmodel_float_t nu = parameters[PARAM_nu];
        const cellmodel_float_t IbarSLCaP = parameters[PARAM_IbarSLCaP];
        const cellmodel_float_t KmPCa = parameters[PARAM_KmPCa];
        const cellmodel_float_t Q10SLCaP = parameters[PARAM_Q10SLCaP];
        const cellmodel_float_t GCaB = parameters[PARAM_GCaB];
        const cellmodel_float_t Kmf = parameters[PARAM_Kmf];
        const cellmodel_float_t Kmr = parameters[PARAM_Kmr];
        const cellmodel_float_t MaxSR = parameters[PARAM_MaxSR];
        const cellmodel_float_t MinSR = parameters[PARAM_MinSR];
        const cellmodel_float_t Q10SRCaP = parameters[PARAM_Q10SRCaP];
        const cellmodel_float_t Vmax_SRCaP = parameters[PARAM_Vmax_SRCaP];
        const cellmodel_float_t ec50SR = parameters[PARAM_ec50SR];
        const cellmodel_float_t hillSRCaP = parameters[PARAM_hillSRCaP];
        const cellmodel_float_t kiCa = parameters[PARAM_kiCa];
        const cellmodel_float_t kim = parameters[PARAM_kim];
        const cellmodel_float_t koCa = parameters[PARAM_koCa];
        const cellmodel_float_t kom = parameters[PARAM_kom];
        const cellmodel_float_t ks = parameters[PARAM_ks];
        const cellmodel_float_t Bmax_Naj = parameters[PARAM_Bmax_Naj];
        const cellmodel_float_t Bmax_Nasl = parameters[PARAM_Bmax_Nasl];
        const cellmodel_float_t koff_na = parameters[PARAM_koff_na];
        const cellmodel_float_t kon_na = parameters[PARAM_kon_na];
        const cellmodel_float_t Bmax_CaM = parameters[PARAM_Bmax_CaM];
        const cellmodel_float_t Bmax_SR = parameters[PARAM_Bmax_SR];
        const cellmodel_float_t Bmax_TnChigh = parameters[PARAM_Bmax_TnChigh];
        const cellmodel_float_t Bmax_TnClow = parameters[PARAM_Bmax_TnClow];
        const cellmodel_float_t Bmax_myosin = parameters[PARAM_Bmax_myosin];
        const cellmodel_float_t koff_cam = parameters[PARAM_koff_cam];
        const cellmodel_float_t koff_myoca = parameters[PARAM_koff_myoca];
        const cellmodel_float_t koff_myomg = parameters[PARAM_koff_myomg];
        const cellmodel_float_t koff_sr = parameters[PARAM_koff_sr];
        const cellmodel_float_t koff_tnchca = parameters[PARAM_koff_tnchca];
        const cellmodel_float_t koff_tnchmg = parameters[PARAM_koff_tnchmg];
        const cellmodel_float_t koff_tncl = parameters[PARAM_koff_tncl];
        const cellmodel_float_t kon_cam = parameters[PARAM_kon_cam];
        const cellmodel_float_t kon_myoca = parameters[PARAM_kon_myoca];
        const cellmodel_float_t kon_myomg = parameters[PARAM_kon_myomg];
        const cellmodel_float_t kon_sr = parameters[PARAM_kon_sr];
        const cellmodel_float_t kon_tnchca = parameters[PARAM_kon_tnchca];
        const cellmodel_float_t kon_tnchmg = parameters[PARAM_kon_tnchmg];
        const cellmodel_float_t kon_tncl = parameters[PARAM_kon_tncl];
        const cellmodel_float_t Bmax_SLhighj0 = parameters[PARAM_Bmax_SLhighj0];
        const cellmodel_float_t Bmax_SLhighsl0 = parameters[PARAM_Bmax_SLhighsl0];
        const cellmodel_float_t Bmax_SLlowj0 = parameters[PARAM_Bmax_SLlowj0];
        const cellmodel_float_t Bmax_SLlowsl0 = parameters[PARAM_Bmax_SLlowsl0];
        const cellmodel_float_t koff_slh = parameters[PARAM_koff_slh];
        const cellmodel_float_t koff_sll = parameters[PARAM_koff_sll];
        const cellmodel_float_t kon_slh = parameters[PARAM_kon_slh];
        const cellmodel_float_t kon_sll = parameters[PARAM_kon_sll];
        const cellmodel_float_t Bmax_Csqn0 = parameters[PARAM_Bmax_Csqn0];
        const cellmodel_float_t J_ca_juncsl = parameters[PARAM_J_ca_juncsl];
        const cellmodel_float_t J_ca_slmyo = parameters[PARAM_J_ca_slmyo];
        const cellmodel_float_t koff_csqn = parameters[PARAM_koff_csqn];
        const cellmodel_float_t kon_csqn = parameters[PARAM_kon_csqn];
        const cellmodel_float_t J_na_juncsl = parameters[PARAM_J_na_juncsl];
        const cellmodel_float_t J_na_slmyo = parameters[PARAM_J_na_slmyo];
        const cellmodel_float_t Nao = parameters[PARAM_Nao];
        const cellmodel_float_t Ko = parameters[PARAM_Ko];
        const cellmodel_float_t Cao = parameters[PARAM_Cao];
        const cellmodel_float_t Cli = parameters[PARAM_Cli];
        const cellmodel_float_t Clo = parameters[PARAM_Clo];
        const cellmodel_float_t Mgi = parameters[PARAM_Mgi];
        const cellmodel_float_t Cmem = parameters[PARAM_Cmem];
        const cellmodel_float_t Frdy = parameters[PARAM_Frdy];
        const cellmodel_float_t R = parameters[PARAM_R];
        const cellmodel_float_t Temp = parameters[PARAM_Temp];
        const cellmodel_float_t stim_amplitude = parameters[PARAM_stim_amplitude];
        const cellmodel_float_t stim_duration = parameters[PARAM_stim_duration];
        const cellmodel_float_t stim_period = parameters[PARAM_stim_period];
        const cellmodel_float_t stim_start = parameters[PARAM_stim_start];

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
            const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
            const cellmodel_float_t h = states[STATE_h * padded_num_cells + i];
            const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
            const cellmodel_float_t x_kr = states[STATE_x_kr * padded_num_cells + i];
            const cellmodel_float_t x_ks = states[STATE_x_ks * padded_num_cells + i];
            const cellmodel_float_t x_to_s = states[STATE_x_to_s * padded_num_cells + i];
            const cellmodel_float_t y_to_s = states[STATE_y_to_s * padded_num_cells + i];
            const cellmodel_float_t x_to_f = states[STATE_x_to_f * padded_num_cells + i];
            const cellmodel_float_t y_to_f = states[STATE_y_to_f * padded_num_cells + i];
            const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
            const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
            const cellmodel_float_t f_Ca_Bj = states[STATE_f_Ca_Bj * padded_num_cells + i];
            const cellmodel_float_t f_Ca_Bsl = states[STATE_f_Ca_Bsl * padded_num_cells + i];
            const cellmodel_float_t Ry_Rr = states[STATE_Ry_Rr * padded_num_cells + i];
            const cellmodel_float_t Ry_Ro = states[STATE_Ry_Ro * padded_num_cells + i];
            const cellmodel_float_t Ry_Ri = states[STATE_Ry_Ri * padded_num_cells + i];
            const cellmodel_float_t Na_Bj = states[STATE_Na_Bj * padded_num_cells + i];
            const cellmodel_float_t Na_Bsl = states[STATE_Na_Bsl * padded_num_cells + i];
            const cellmodel_float_t Tn_CL = states[STATE_Tn_CL * padded_num_cells + i];
            const cellmodel_float_t Tn_CHc = states[STATE_Tn_CHc * padded_num_cells + i];
            const cellmodel_float_t Tn_CHm = states[STATE_Tn_CHm * padded_num_cells + i];
            const cellmodel_float_t CaM = states[STATE_CaM * padded_num_cells + i];
            const cellmodel_float_t Myo_c = states[STATE_Myo_c * padded_num_cells + i];
            const cellmodel_float_t Myo_m = states[STATE_Myo_m * padded_num_cells + i];
            const cellmodel_float_t SRB = states[STATE_SRB * padded_num_cells + i];
            const cellmodel_float_t SLL_j = states[STATE_SLL_j * padded_num_cells + i];
            const cellmodel_float_t SLL_sl = states[STATE_SLL_sl * padded_num_cells + i];
            const cellmodel_float_t SLH_j = states[STATE_SLH_j * padded_num_cells + i];
            const cellmodel_float_t SLH_sl = states[STATE_SLH_sl * padded_num_cells + i];
            const cellmodel_float_t Csqn_b = states[STATE_Csqn_b * padded_num_cells + i];
            const cellmodel_float_t Ca_sr = states[STATE_Ca_sr * padded_num_cells + i];
            const cellmodel_float_t Na_j = states[STATE_Na_j * padded_num_cells + i];
            const cellmodel_float_t Na_sl = states[STATE_Na_sl * padded_num_cells + i];
            const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];
            const cellmodel_float_t K_i = states[STATE_K_i * padded_num_cells + i];
            const cellmodel_float_t Ca_j = states[STATE_Ca_j * padded_num_cells + i];
            const cellmodel_float_t Ca_sl = states[STATE_Ca_sl * padded_num_cells + i];
            const cellmodel_float_t Ca_i = states[STATE_Ca_i * padded_num_cells + i];
            const cellmodel_float_t V_m = states[STATE_V_m * padded_num_cells + i];

            // Expressions for the Geometry component
            const cellmodel_float_t Vcell =
                    FP_LITERAL(1.0e-15) * M_PI * cellLength * (cellRadius * cellRadius);
            const cellmodel_float_t Vmyo = FP_LITERAL(0.65) * Vcell;
            const cellmodel_float_t Vsr = FP_LITERAL(0.035) * Vcell;
            const cellmodel_float_t Vsl = FP_LITERAL(0.02) * Vcell;
            const cellmodel_float_t Vjunc = FP_LITERAL(0.000539) * Vcell;
            const cellmodel_float_t Fsl = FP_LITERAL(1.) - Fjunc;
            const cellmodel_float_t Fsl_CaL = FP_LITERAL(1.) - Fjunc_CaL;

            // Expressions for the Reversal potentials component
            const cellmodel_float_t FoRT = Frdy / (R * Temp);
            const cellmodel_float_t ena_junc = Log(Nao / Na_j) / FoRT;
            const cellmodel_float_t ena_sl = Log(Nao / Na_sl) / FoRT;
            const cellmodel_float_t ek = Log(Ko / K_i) / FoRT;
            const cellmodel_float_t eca_junc = Log(Cao / Ca_j) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t eca_sl = Log(Cao / Ca_sl) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t ecl = Log(Cli / Clo) / FoRT;
            const cellmodel_float_t Qpow = FP_LITERAL(-31.) + Temp / FP_LITERAL(10.);

            // Expressions for the I_Na component
            const cellmodel_float_t mss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(0.00184221158116513)
                                  * Exp(FP_LITERAL(-0.110741971207087) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(0.00184221158116513)
                                    * Exp(FP_LITERAL(-0.110741971207087) * V_m)));
            const cellmodel_float_t taum =
                    FP_LITERAL(0.1292)
                            * Exp(-((FP_LITERAL(2.94658944658945)
                                     + FP_LITERAL(0.0643500643500644) * V_m)
                                    * (FP_LITERAL(2.94658944658945)
                                       + FP_LITERAL(0.0643500643500644) * V_m)))
                    + FP_LITERAL(0.06487)
                              * Exp(-((FP_LITERAL(-0.0943466353677621)
                                       + FP_LITERAL(0.0195618153364632) * V_m)
                                      * (FP_LITERAL(-0.0943466353677621)
                                         + FP_LITERAL(0.0195618153364632) * V_m)));
            const cellmodel_float_t ah =
                    (V_m >= FP_LITERAL(-40.) ? FP_LITERAL(0.)
                                             : FP_LITERAL(4.43126792958051e-7)
                                                       * Exp(FP_LITERAL(-0.147058823529412) * V_m));
            const cellmodel_float_t bh =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.77)
                                       / (FP_LITERAL(0.13)
                                          + FP_LITERAL(0.0497581410839387)
                                                    * Exp(FP_LITERAL(-0.0900900900900901) * V_m))
                             : FP_LITERAL(310000.0) * Exp(FP_LITERAL(0.3485) * V_m)
                                       + FP_LITERAL(2.7) * Exp(FP_LITERAL(0.079) * V_m));
            const cellmodel_float_t tauh = FP_LITERAL(1.0) / (ah + bh);
            const cellmodel_float_t hss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(15212.5932856544)
                                    * Exp(FP_LITERAL(0.134589502018843) * V_m)));
            const cellmodel_float_t aj =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.)
                             : (FP_LITERAL(37.78) + V_m)
                                       * (FP_LITERAL(-25428.0) * Exp(FP_LITERAL(0.2444) * V_m)
                                          - FP_LITERAL(6.948e-6) * Exp(FP_LITERAL(-0.04391) * V_m))
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(50262745825.954)
                                                    * Exp(FP_LITERAL(0.311) * V_m)));
            const cellmodel_float_t bj =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.6) * Exp(FP_LITERAL(0.057) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.0407622039783662)
                                                    * Exp(FP_LITERAL(-0.1) * V_m))
                             : FP_LITERAL(0.02424) * Exp(FP_LITERAL(-0.01052) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.00396086833990426)
                                                    * Exp(FP_LITERAL(-0.1378) * V_m)));
            const cellmodel_float_t tauj = FP_LITERAL(1.0) / (aj + bj);
            const cellmodel_float_t jss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(15212.5932856544)
                                    * Exp(FP_LITERAL(0.134589502018843) * V_m)));
            const cellmodel_float_t dm_dt = (-m + mss) / taum;
            states[STATE_m * padded_num_cells + i] = dt * dm_dt + m;
            const cellmodel_float_t dh_dt = (-h + hss) / tauh;
            states[STATE_h * padded_num_cells + i] = dt * dh_dt + h;
            const cellmodel_float_t dj_dt = (-j + jss) / tauj;
            states[STATE_j * padded_num_cells + i] = dt * dj_dt + j;
            const cellmodel_float_t I_Na_junc =
                    Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j;
            const cellmodel_float_t I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j;

            // Expressions for the I_NaBK component
            const cellmodel_float_t I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m);
            const cellmodel_float_t I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl;

            // Expressions for the I_NaK component
            const cellmodel_float_t sigma =
                    FP_LITERAL(-1.) / FP_LITERAL(7.)
                    + Exp(FP_LITERAL(0.0148588410104012) * Nao) / FP_LITERAL(7.);
            const cellmodel_float_t fnak =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                       + FP_LITERAL(0.0365) * Exp(-FoRT * V_m) * sigma);
            const cellmodel_float_t I_nak_junc =
                    Fjunc * IbarNaK * Ko * fnak
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_j) * (Na_j)) * ((Na_j) * (Na_j))))
                       * (KmKo + Ko));
            const cellmodel_float_t I_nak_sl =
                    IbarNaK * Ko * Fsl * fnak
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_sl) * (Na_sl)) * ((Na_sl) * (Na_sl))))
                       * (KmKo + Ko));
            const cellmodel_float_t I_nak = I_nak_junc + I_nak_sl;

            // Expressions for the I_Kr component
            const cellmodel_float_t gkr = FP_LITERAL(0.430331482911935) * GKr * sqrt(Ko);
            const cellmodel_float_t xrss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(-2.) - V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tauxr =
                    FP_LITERAL(230.)
                            / (FP_LITERAL(1.) + Exp(FP_LITERAL(2.) + V_m / FP_LITERAL(20.)))
                    + FP_LITERAL(3300.)
                              / ((FP_LITERAL(1.)
                                  + Exp(FP_LITERAL(-22.) / FP_LITERAL(9.) - V_m / FP_LITERAL(9.)))
                                 * (FP_LITERAL(1.)
                                    + Exp(FP_LITERAL(11.) / FP_LITERAL(9.)
                                          + V_m / FP_LITERAL(9.))));
            const cellmodel_float_t dx_kr_dt = (-x_kr + xrss) / tauxr;
            states[STATE_x_kr * padded_num_cells + i] = dt * dx_kr_dt + x_kr;
            const cellmodel_float_t rkr =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(37.) / FP_LITERAL(12.) + V_m / FP_LITERAL(24.)));
            const cellmodel_float_t I_kr = (-ek + V_m) * gkr * rkr * x_kr;

            // Expressions for the I_Kp component
            const cellmodel_float_t kp_kp =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(1786.47556537862) * Exp(FP_LITERAL(-0.167224080267559) * V_m));
            const cellmodel_float_t I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp;
            const cellmodel_float_t I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp;
            const cellmodel_float_t I_kp = I_kp_junc + I_kp_sl;

            // Expressions for the I_Ks component
            const cellmodel_float_t eks = Log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT;
            const cellmodel_float_t gks_junc = GKs;
            const cellmodel_float_t gks_sl = GKs;
            const cellmodel_float_t xsss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.765928338364649)
                                 * Exp(FP_LITERAL(-0.0701754385964912) * V_m));
            const cellmodel_float_t tauxs =
                    FP_LITERAL(990.1)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.841540408868102)
                                 * Exp(FP_LITERAL(-0.0708215297450425) * V_m));
            const cellmodel_float_t dx_ks_dt = (-x_ks + xsss) / tauxs;
            states[STATE_x_ks * padded_num_cells + i] = dt * dx_ks_dt + x_ks;
            const cellmodel_float_t I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc;
            const cellmodel_float_t I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl;
            const cellmodel_float_t I_ks = I_ks_junc + I_ks_sl;

            // Expressions for the I_to component
            const cellmodel_float_t GtoSlow =
                    (epi == FP_LITERAL(1.) ? FP_LITERAL(0.12) * Gto : FP_LITERAL(0.2892) * Gto);
            const cellmodel_float_t GtoFast =
                    (epi == FP_LITERAL(1.) ? FP_LITERAL(0.88) * Gto : FP_LITERAL(0.0108) * Gto);
            const cellmodel_float_t xtoss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(19.) / FP_LITERAL(13.) - V_m / FP_LITERAL(13.)));
            const cellmodel_float_t ytoss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(49.4024491055302) * Exp(V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tauxtos =
                    FP_LITERAL(0.5)
                    + FP_LITERAL(9.)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(1.) / FP_LITERAL(5.) + V_m / FP_LITERAL(15.)));
            const cellmodel_float_t tauytos =
                    FP_LITERAL(30.)
                    + FP_LITERAL(800.)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(6.) + V_m / FP_LITERAL(10.)));
            const cellmodel_float_t dx_to_s_dt = (-x_to_s + xtoss) / tauxtos;
            states[STATE_x_to_s * padded_num_cells + i] = dt * dx_to_s_dt + x_to_s;
            const cellmodel_float_t dy_to_s_dt = (-y_to_s + ytoss) / tauytos;
            states[STATE_y_to_s * padded_num_cells + i] = dt * dy_to_s_dt + y_to_s;
            const cellmodel_float_t I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s;
            const cellmodel_float_t tauxtof =
                    FP_LITERAL(0.5)
                    + FP_LITERAL(8.5)
                              * Exp(-((FP_LITERAL(9.) / FP_LITERAL(10.) + V_m / FP_LITERAL(50.))
                                      * (FP_LITERAL(9.) / FP_LITERAL(10.)
                                         + V_m / FP_LITERAL(50.))));
            const cellmodel_float_t tauytof =
                    FP_LITERAL(7.)
                    + FP_LITERAL(85.)
                              * Exp(-((FP_LITERAL(40.) + V_m) * (FP_LITERAL(40.) + V_m))
                                    / FP_LITERAL(220.));
            const cellmodel_float_t dx_to_f_dt = (-x_to_f + xtoss) / tauxtof;
            states[STATE_x_to_f * padded_num_cells + i] = dt * dx_to_f_dt + x_to_f;
            const cellmodel_float_t dy_to_f_dt = (-y_to_f + ytoss) / tauytof;
            states[STATE_y_to_f * padded_num_cells + i] = dt * dy_to_f_dt + y_to_f;
            const cellmodel_float_t I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f;
            const cellmodel_float_t I_to = I_tof + I_tos;

            // Expressions for the I_K1 component
            const cellmodel_float_t aki =
                    FP_LITERAL(1.02)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(7.35454251046446e-7)
                                 * Exp(FP_LITERAL(0.2385) * V_m - FP_LITERAL(0.2385) * ek));
            const cellmodel_float_t bki =
                    (FP_LITERAL(0.762624006506308)
                             * Exp(FP_LITERAL(0.08032) * V_m - FP_LITERAL(0.08032) * ek)
                     + FP_LITERAL(1.15340563518656e-16)
                               * Exp(FP_LITERAL(0.06175) * V_m - FP_LITERAL(0.06175) * ek))
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.0867722941576933)
                                 * Exp(FP_LITERAL(0.5143) * ek - FP_LITERAL(0.5143) * V_m));
            const cellmodel_float_t kiss = aki / (aki + bki);
            const cellmodel_float_t I_K1 =
                    FP_LITERAL(0.430331482911935) * GK1 * sqrt(Ko) * (-ek + V_m) * kiss;

            // Expressions for the I_ClCa component
            const cellmodel_float_t I_ClCa_junc =
                    Fjunc * GClCa * (-ecl + V_m) / (FP_LITERAL(1.) + KdClCa / Ca_j);
            const cellmodel_float_t I_ClCa_sl =
                    GClCa * (-ecl + V_m) * Fsl / (FP_LITERAL(1.) + KdClCa / Ca_sl);
            const cellmodel_float_t I_ClCa = I_ClCa_junc + I_ClCa_sl;
            const cellmodel_float_t I_Clbk = GClB * (-ecl + V_m);

            // Expressions for the I_Ca component
            const cellmodel_float_t fss =
                    FP_LITERAL(1.0)
                            / (FP_LITERAL(1.)
                               + Exp(FP_LITERAL(35.) / FP_LITERAL(9.) + V_m / FP_LITERAL(9.)))
                    + FP_LITERAL(0.6)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V_m / FP_LITERAL(20.)));
            const cellmodel_float_t dss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)));
            const cellmodel_float_t taud =
                    (FP_LITERAL(1.) - Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)))
                    * dss / (FP_LITERAL(0.175) + FP_LITERAL(0.035) * V_m);
            const cellmodel_float_t tauf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(0.02)
                       + FP_LITERAL(0.0197)
                                 * Exp(-((FP_LITERAL(0.48865) + FP_LITERAL(0.0337) * V_m)
                                         * (FP_LITERAL(0.48865) + FP_LITERAL(0.0337) * V_m))));
            const cellmodel_float_t dd_dt = (-d + dss) / taud;
            states[STATE_d * padded_num_cells + i] = dt * dd_dt + d;
            const cellmodel_float_t df_dt = (-f + fss) / tauf;
            states[STATE_f * padded_num_cells + i] = dt * df_dt + f;
            const cellmodel_float_t df_Ca_Bj_dt =
                    FP_LITERAL(-0.0119) * f_Ca_Bj
                    + FP_LITERAL(1.7) * (FP_LITERAL(1.) - f_Ca_Bj) * Ca_j;
            states[STATE_f_Ca_Bj * padded_num_cells + i] = dt * df_Ca_Bj_dt + f_Ca_Bj;
            const cellmodel_float_t df_Ca_Bsl_dt =
                    FP_LITERAL(-0.0119) * f_Ca_Bsl
                    + FP_LITERAL(1.7) * (FP_LITERAL(1.) - f_Ca_Bsl) * Ca_sl;
            states[STATE_f_Ca_Bsl * padded_num_cells + i] = dt * df_Ca_Bsl_dt + f_Ca_Bsl;
            const cellmodel_float_t fcaCaMSL = FP_LITERAL(0.);
            const cellmodel_float_t fcaCaj = FP_LITERAL(0.);
            const cellmodel_float_t ibarca_j =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                    * (FP_LITERAL(-0.341) * Cao
                       + FP_LITERAL(0.341) * Ca_j * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t ibarca_sl =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                    * (FP_LITERAL(-0.341) * Cao
                       + FP_LITERAL(0.341) * Ca_sl * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t ibark =
                    Frdy * GCaL * pK
                    * (FP_LITERAL(-0.75) * Ko + FP_LITERAL(0.75) * K_i * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t ibarna_j =
                    Frdy * GCaL * pNa
                    * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_j * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t ibarna_sl =
                    Frdy * GCaL * pNa
                    * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_sl * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t I_Ca_junc = FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                                                * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f
                                                * ibarca_j;
            const cellmodel_float_t I_Ca_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                              * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d
                                              * f * ibarca_sl;
            const cellmodel_float_t I_CaK = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                            * (Fjunc_CaL * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj)
                                               + (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
                                            * d * f * ibark;
            const cellmodel_float_t I_CaNa_junc = FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                                                  * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f
                                                  * ibarna_j;
            const cellmodel_float_t I_CaNa_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                                * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL
                                                * d * f * ibarna_sl;

            // Expressions for the I_NCX component
            const cellmodel_float_t Ka_junc =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_j * Ca_j));
            const cellmodel_float_t Ka_sl =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_sl * Ca_sl));
            const cellmodel_float_t s1_junc = Cao * (Na_j * Na_j * Na_j) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s2_junc =
                    (Nao * Nao * Nao) * Ca_j * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_junc =
                    Cao * (Na_j * Na_j * Na_j) + KmCao * (Na_j * Na_j * Na_j)
                    + (Nao * Nao * Nao) * Ca_j
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_j / KmCai) * Ca_j;
            const cellmodel_float_t s2_sl =
                    (Nao * Nao * Nao) * Ca_sl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_sl =
                    Cao * (Na_sl * Na_sl * Na_sl) + KmCao * (Na_sl * Na_sl * Na_sl)
                    + (Nao * Nao * Nao) * Ca_sl
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_sl / KmCai) * Ca_sl;
            const cellmodel_float_t I_ncx_junc =
                    Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * s3_junc);
            const cellmodel_float_t I_ncx_sl =
                    IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);

            // Expressions for the I_PCa component
            const cellmodel_float_t I_pca_junc =
                    Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, FP_LITERAL(1.6))
                    / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_j, FP_LITERAL(1.6)));
            const cellmodel_float_t I_pca_sl =
                    IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, FP_LITERAL(1.6)) * Fsl
                    / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_sl, FP_LITERAL(1.6)));

            // Expressions for the I_CaBK component
            const cellmodel_float_t I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m);
            const cellmodel_float_t I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl;

            // Expressions for the SR Fluxes component
            const cellmodel_float_t kCaSR =
                    MaxSR
                    - (MaxSR - MinSR) / (FP_LITERAL(1.) + pow(ec50SR / Ca_sr, FP_LITERAL(2.5)));
            const cellmodel_float_t koSRCa = koCa / kCaSR;
            const cellmodel_float_t kiSRCa = kiCa * kCaSR;
            const cellmodel_float_t RI = FP_LITERAL(1.) - Ry_Ri - Ry_Ro - Ry_Rr;
            const cellmodel_float_t dRy_Rr_dt =
                    kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa;
            states[STATE_Ry_Rr * padded_num_cells + i] = dt * dRy_Rr_dt + Ry_Rr;
            const cellmodel_float_t dRy_Ro_dt = kim * Ry_Ri - kom * Ry_Ro
                                                + (Ca_j * Ca_j) * Ry_Rr * koSRCa
                                                - Ca_j * Ry_Ro * kiSRCa;
            states[STATE_Ry_Ro * padded_num_cells + i] = dt * dRy_Ro_dt + Ry_Ro;
            const cellmodel_float_t dRy_Ri_dt = -kim * Ry_Ri - kom * Ry_Ri
                                                + (Ca_j * Ca_j) * RI * koSRCa
                                                + Ca_j * Ry_Ro * kiSRCa;
            states[STATE_Ry_Ri * padded_num_cells + i] = dt * dRy_Ri_dt + Ry_Ri;
            const cellmodel_float_t J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro;
            const cellmodel_float_t J_serca =
                    Vmax_SRCaP * pow(Q10SRCaP, Qpow)
                    * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                    / (FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP));
            const cellmodel_float_t J_SRleak =
                    FP_LITERAL(5.348e-6) * Ca_sr - FP_LITERAL(5.348e-6) * Ca_j;

            // Expressions for the Na Buffers component
            const cellmodel_float_t dNa_Bj_dt =
                    -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j;
            states[STATE_Na_Bj * padded_num_cells + i] = dt * dNa_Bj_dt + Na_Bj;
            const cellmodel_float_t dNa_Bsl_dt =
                    -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl;
            states[STATE_Na_Bsl * padded_num_cells + i] = dt * dNa_Bsl_dt + Na_Bsl;

            // Expressions for the Cytosolic Ca Buffers component
            const cellmodel_float_t dTn_CL_dt =
                    -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;
            states[STATE_Tn_CL * padded_num_cells + i] = dt * dTn_CL_dt + Tn_CL;
            const cellmodel_float_t dTn_CHc_dt =
                    -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i;
            states[STATE_Tn_CHc * padded_num_cells + i] = dt * dTn_CHc_dt + Tn_CHc;
            const cellmodel_float_t dTn_CHm_dt =
                    -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm);
            states[STATE_Tn_CHm * padded_num_cells + i] = dt * dTn_CHm_dt + Tn_CHm;
            const cellmodel_float_t dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i;
            states[STATE_CaM * padded_num_cells + i] = dt * dCaM_dt + CaM;
            const cellmodel_float_t dMyo_c_dt =
                    -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i;
            states[STATE_Myo_c * padded_num_cells + i] = dt * dMyo_c_dt + Myo_c;
            const cellmodel_float_t dMyo_m_dt =
                    -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m);
            states[STATE_Myo_m * padded_num_cells + i] = dt * dMyo_m_dt + Myo_m;
            const cellmodel_float_t dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i;
            states[STATE_SRB * padded_num_cells + i] = dt * dSRB_dt + SRB;
            const cellmodel_float_t J_CaB_cytosol =
                    -koff_cam * CaM - koff_myoca * Myo_c - koff_myomg * Myo_m - koff_sr * SRB
                    - koff_tnchca * Tn_CHc - koff_tnchmg * Tn_CHm - koff_tncl * Tn_CL
                    + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
                    + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
                    + kon_cam * (Bmax_CaM - CaM) * Ca_i
                    + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
                    + kon_sr * (Bmax_SR - SRB) * Ca_i
                    + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
                    + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;

            // Expressions for the Junctional and SL Ca Buffers component
            const cellmodel_float_t Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl;
            const cellmodel_float_t Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc;
            const cellmodel_float_t Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl;
            const cellmodel_float_t Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc;
            const cellmodel_float_t dSLL_j_dt =
                    -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
            states[STATE_SLL_j * padded_num_cells + i] = dt * dSLL_j_dt + SLL_j;
            const cellmodel_float_t dSLL_sl_dt =
                    -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;
            states[STATE_SLL_sl * padded_num_cells + i] = dt * dSLL_sl_dt + SLL_sl;
            const cellmodel_float_t dSLH_j_dt =
                    -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j;
            states[STATE_SLH_j * padded_num_cells + i] = dt * dSLH_j_dt + SLH_j;
            const cellmodel_float_t dSLH_sl_dt =
                    -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl;
            states[STATE_SLH_sl * padded_num_cells + i] = dt * dSLH_sl_dt + SLH_sl;
            const cellmodel_float_t J_CaB_junction = -koff_slh * SLH_j - koff_sll * SLL_j
                                                     + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
                                                     + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
            const cellmodel_float_t J_CaB_sl = -koff_slh * SLH_sl - koff_sll * SLL_sl
                                               + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
                                               + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;

            // Expressions for the SR Ca Concentrations component
            const cellmodel_float_t Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr;
            const cellmodel_float_t dCsqn_b_dt =
                    -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr;
            states[STATE_Csqn_b * padded_num_cells + i] = dt * dCsqn_b_dt + Csqn_b;
            const cellmodel_float_t dCa_sr_dt = -J_SRCarel + koff_csqn * Csqn_b
                                                - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
                                                - J_SRleak * Vmyo / Vsr + J_serca;
            states[STATE_Ca_sr * padded_num_cells + i] = dt * dCa_sr_dt + Ca_sr;

            // Expressions for the Na Concentrations component
            const cellmodel_float_t I_Na_tot_junc = FP_LITERAL(3.) * I_nak_junc
                                                    + FP_LITERAL(3.) * I_ncx_junc + I_CaNa_junc
                                                    + I_Na_junc + I_nabk_junc;
            const cellmodel_float_t I_Na_tot_sl = FP_LITERAL(3.) * I_nak_sl
                                                  + FP_LITERAL(3.) * I_ncx_sl + I_CaNa_sl + I_Na_sl
                                                  + I_nabk_sl;
            const cellmodel_float_t dNa_j_dt = -dNa_Bj_dt + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
                                               - Cmem * I_Na_tot_junc / (Frdy * Vjunc);
            states[STATE_Na_j * padded_num_cells + i] = dt * dNa_j_dt + Na_j;
            const cellmodel_float_t dNa_sl_dt = -dNa_Bsl_dt + J_na_juncsl * (-Na_sl + Na_j) / Vsl
                                                + J_na_slmyo * (-Na_sl + Na_i) / Vsl
                                                - Cmem * I_Na_tot_sl / (Frdy * Vsl);
            states[STATE_Na_sl * padded_num_cells + i] = dt * dNa_sl_dt + Na_sl;
            const cellmodel_float_t dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo;
            states[STATE_Na_i * padded_num_cells + i] = dt * dNa_i_dt + Na_i;

            // Expressions for the K Concentration component
            const cellmodel_float_t I_K_tot =
                    FP_LITERAL(-2.) * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to;
            const cellmodel_float_t dK_i_dt = FP_LITERAL(0.);
            states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;

            // Expressions for the Ca Concentrations component
            const cellmodel_float_t I_Ca_tot_junc =
                    FP_LITERAL(-2.) * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
            const cellmodel_float_t I_Ca_tot_sl =
                    FP_LITERAL(-2.) * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
            const cellmodel_float_t dCa_j_dt =
                    -J_CaB_junction + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
                    + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc
                    - Cmem * I_Ca_tot_junc / (FP_LITERAL(2.) * Frdy * Vjunc);
            states[STATE_Ca_j * padded_num_cells + i] = dt * dCa_j_dt + Ca_j;
            const cellmodel_float_t dCa_sl_dt =
                    -J_CaB_sl + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
                    + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
                    - Cmem * I_Ca_tot_sl / (FP_LITERAL(2.) * Frdy * Vsl);
            states[STATE_Ca_sl * padded_num_cells + i] = dt * dCa_sl_dt + Ca_sl;
            const cellmodel_float_t dCa_i_dt =
                    -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo;
            states[STATE_Ca_i * padded_num_cells + i] = dt * dCa_i_dt + Ca_i;

            // Expressions for the Membrane potential component
            const cellmodel_float_t i_Stim =
                    (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                     && t - stim_period * floor(t / stim_period) >= stim_start
                             ? -stim_amplitude
                             : FP_LITERAL(0.));
            const cellmodel_float_t I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
            const cellmodel_float_t I_Cl_tot = I_ClCa + I_Clbk;
            const cellmodel_float_t I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
            const cellmodel_float_t I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot;
            const cellmodel_float_t dV_m_dt = -I_tot - i_Stim;
            states[STATE_V_m * padded_num_cells + i] = dt * dV_m_dt + V_m;
        }
    }
}

// Compute a forward step using the rush larsen algorithm to the GPB ODE
static void step_GRL1(cellmodel_float_t *__restrict states, const cellmodel_float_t t,
                      const cellmodel_float_t dt, const cellmodel_float_t *__restrict parameters,
                      const long num_cells, long padded_num_cells)
{
    #pragma omp parallel
    {
        // Assign parameters
        const cellmodel_float_t Fjunc = parameters[PARAM_Fjunc];
        const cellmodel_float_t Fjunc_CaL = parameters[PARAM_Fjunc_CaL];
        const cellmodel_float_t cellLength = parameters[PARAM_cellLength];
        const cellmodel_float_t cellRadius = parameters[PARAM_cellRadius];
        const cellmodel_float_t GNa = parameters[PARAM_GNa];
        const cellmodel_float_t GNaB = parameters[PARAM_GNaB];
        const cellmodel_float_t IbarNaK = parameters[PARAM_IbarNaK];
        const cellmodel_float_t KmKo = parameters[PARAM_KmKo];
        const cellmodel_float_t KmNaip = parameters[PARAM_KmNaip];
        const cellmodel_float_t GKr = parameters[PARAM_GKr];
        const cellmodel_float_t GKp = parameters[PARAM_GKp];
        const cellmodel_float_t GKs = parameters[PARAM_GKs];
        const cellmodel_float_t pNaK = parameters[PARAM_pNaK];
        const cellmodel_float_t GK1 = parameters[PARAM_GK1];
        const cellmodel_float_t Gto = parameters[PARAM_Gto];
        const cellmodel_float_t epi = parameters[PARAM_epi];
        const cellmodel_float_t GClB = parameters[PARAM_GClB];
        const cellmodel_float_t GClCa = parameters[PARAM_GClCa];
        const cellmodel_float_t KdClCa = parameters[PARAM_KdClCa];
        const cellmodel_float_t GCaL = parameters[PARAM_GCaL];
        const cellmodel_float_t Q10CaL = parameters[PARAM_Q10CaL];
        const cellmodel_float_t pCa = parameters[PARAM_pCa];
        const cellmodel_float_t pK = parameters[PARAM_pK];
        const cellmodel_float_t pNa = parameters[PARAM_pNa];
        const cellmodel_float_t IbarNCX = parameters[PARAM_IbarNCX];
        const cellmodel_float_t Kdact = parameters[PARAM_Kdact];
        const cellmodel_float_t KmCai = parameters[PARAM_KmCai];
        const cellmodel_float_t KmCao = parameters[PARAM_KmCao];
        const cellmodel_float_t KmNai = parameters[PARAM_KmNai];
        const cellmodel_float_t KmNao = parameters[PARAM_KmNao];
        const cellmodel_float_t Q10NCX = parameters[PARAM_Q10NCX];
        const cellmodel_float_t ksat = parameters[PARAM_ksat];
        const cellmodel_float_t nu = parameters[PARAM_nu];
        const cellmodel_float_t IbarSLCaP = parameters[PARAM_IbarSLCaP];
        const cellmodel_float_t KmPCa = parameters[PARAM_KmPCa];
        const cellmodel_float_t Q10SLCaP = parameters[PARAM_Q10SLCaP];
        const cellmodel_float_t GCaB = parameters[PARAM_GCaB];
        const cellmodel_float_t Kmf = parameters[PARAM_Kmf];
        const cellmodel_float_t Kmr = parameters[PARAM_Kmr];
        const cellmodel_float_t MaxSR = parameters[PARAM_MaxSR];
        const cellmodel_float_t MinSR = parameters[PARAM_MinSR];
        const cellmodel_float_t Q10SRCaP = parameters[PARAM_Q10SRCaP];
        const cellmodel_float_t Vmax_SRCaP = parameters[PARAM_Vmax_SRCaP];
        const cellmodel_float_t ec50SR = parameters[PARAM_ec50SR];
        const cellmodel_float_t hillSRCaP = parameters[PARAM_hillSRCaP];
        const cellmodel_float_t kiCa = parameters[PARAM_kiCa];
        const cellmodel_float_t kim = parameters[PARAM_kim];
        const cellmodel_float_t koCa = parameters[PARAM_koCa];
        const cellmodel_float_t kom = parameters[PARAM_kom];
        const cellmodel_float_t ks = parameters[PARAM_ks];
        const cellmodel_float_t Bmax_Naj = parameters[PARAM_Bmax_Naj];
        const cellmodel_float_t Bmax_Nasl = parameters[PARAM_Bmax_Nasl];
        const cellmodel_float_t koff_na = parameters[PARAM_koff_na];
        const cellmodel_float_t kon_na = parameters[PARAM_kon_na];
        const cellmodel_float_t Bmax_CaM = parameters[PARAM_Bmax_CaM];
        const cellmodel_float_t Bmax_SR = parameters[PARAM_Bmax_SR];
        const cellmodel_float_t Bmax_TnChigh = parameters[PARAM_Bmax_TnChigh];
        const cellmodel_float_t Bmax_TnClow = parameters[PARAM_Bmax_TnClow];
        const cellmodel_float_t Bmax_myosin = parameters[PARAM_Bmax_myosin];
        const cellmodel_float_t koff_cam = parameters[PARAM_koff_cam];
        const cellmodel_float_t koff_myoca = parameters[PARAM_koff_myoca];
        const cellmodel_float_t koff_myomg = parameters[PARAM_koff_myomg];
        const cellmodel_float_t koff_sr = parameters[PARAM_koff_sr];
        const cellmodel_float_t koff_tnchca = parameters[PARAM_koff_tnchca];
        const cellmodel_float_t koff_tnchmg = parameters[PARAM_koff_tnchmg];
        const cellmodel_float_t koff_tncl = parameters[PARAM_koff_tncl];
        const cellmodel_float_t kon_cam = parameters[PARAM_kon_cam];
        const cellmodel_float_t kon_myoca = parameters[PARAM_kon_myoca];
        const cellmodel_float_t kon_myomg = parameters[PARAM_kon_myomg];
        const cellmodel_float_t kon_sr = parameters[PARAM_kon_sr];
        const cellmodel_float_t kon_tnchca = parameters[PARAM_kon_tnchca];
        const cellmodel_float_t kon_tnchmg = parameters[PARAM_kon_tnchmg];
        const cellmodel_float_t kon_tncl = parameters[PARAM_kon_tncl];
        const cellmodel_float_t Bmax_SLhighj0 = parameters[PARAM_Bmax_SLhighj0];
        const cellmodel_float_t Bmax_SLhighsl0 = parameters[PARAM_Bmax_SLhighsl0];
        const cellmodel_float_t Bmax_SLlowj0 = parameters[PARAM_Bmax_SLlowj0];
        const cellmodel_float_t Bmax_SLlowsl0 = parameters[PARAM_Bmax_SLlowsl0];
        const cellmodel_float_t koff_slh = parameters[PARAM_koff_slh];
        const cellmodel_float_t koff_sll = parameters[PARAM_koff_sll];
        const cellmodel_float_t kon_slh = parameters[PARAM_kon_slh];
        const cellmodel_float_t kon_sll = parameters[PARAM_kon_sll];
        const cellmodel_float_t Bmax_Csqn0 = parameters[PARAM_Bmax_Csqn0];
        const cellmodel_float_t J_ca_juncsl = parameters[PARAM_J_ca_juncsl];
        const cellmodel_float_t J_ca_slmyo = parameters[PARAM_J_ca_slmyo];
        const cellmodel_float_t koff_csqn = parameters[PARAM_koff_csqn];
        const cellmodel_float_t kon_csqn = parameters[PARAM_kon_csqn];
        const cellmodel_float_t J_na_juncsl = parameters[PARAM_J_na_juncsl];
        const cellmodel_float_t J_na_slmyo = parameters[PARAM_J_na_slmyo];
        const cellmodel_float_t Nao = parameters[PARAM_Nao];
        const cellmodel_float_t Ko = parameters[PARAM_Ko];
        const cellmodel_float_t Cao = parameters[PARAM_Cao];
        const cellmodel_float_t Cli = parameters[PARAM_Cli];
        const cellmodel_float_t Clo = parameters[PARAM_Clo];
        const cellmodel_float_t Mgi = parameters[PARAM_Mgi];
        const cellmodel_float_t Cmem = parameters[PARAM_Cmem];
        const cellmodel_float_t Frdy = parameters[PARAM_Frdy];
        const cellmodel_float_t R = parameters[PARAM_R];
        const cellmodel_float_t Temp = parameters[PARAM_Temp];
        const cellmodel_float_t stim_amplitude = parameters[PARAM_stim_amplitude];
        const cellmodel_float_t stim_duration = parameters[PARAM_stim_duration];
        const cellmodel_float_t stim_period = parameters[PARAM_stim_period];
        const cellmodel_float_t stim_start = parameters[PARAM_stim_start];

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
            const cellmodel_float_t m = states[STATE_m * padded_num_cells + i];
            const cellmodel_float_t h = states[STATE_h * padded_num_cells + i];
            const cellmodel_float_t j = states[STATE_j * padded_num_cells + i];
            const cellmodel_float_t x_kr = states[STATE_x_kr * padded_num_cells + i];
            const cellmodel_float_t x_ks = states[STATE_x_ks * padded_num_cells + i];
            const cellmodel_float_t x_to_s = states[STATE_x_to_s * padded_num_cells + i];
            const cellmodel_float_t y_to_s = states[STATE_y_to_s * padded_num_cells + i];
            const cellmodel_float_t x_to_f = states[STATE_x_to_f * padded_num_cells + i];
            const cellmodel_float_t y_to_f = states[STATE_y_to_f * padded_num_cells + i];
            const cellmodel_float_t d = states[STATE_d * padded_num_cells + i];
            const cellmodel_float_t f = states[STATE_f * padded_num_cells + i];
            const cellmodel_float_t f_Ca_Bj = states[STATE_f_Ca_Bj * padded_num_cells + i];
            const cellmodel_float_t f_Ca_Bsl = states[STATE_f_Ca_Bsl * padded_num_cells + i];
            const cellmodel_float_t Ry_Rr = states[STATE_Ry_Rr * padded_num_cells + i];
            const cellmodel_float_t Ry_Ro = states[STATE_Ry_Ro * padded_num_cells + i];
            const cellmodel_float_t Ry_Ri = states[STATE_Ry_Ri * padded_num_cells + i];
            const cellmodel_float_t Na_Bj = states[STATE_Na_Bj * padded_num_cells + i];
            const cellmodel_float_t Na_Bsl = states[STATE_Na_Bsl * padded_num_cells + i];
            const cellmodel_float_t Tn_CL = states[STATE_Tn_CL * padded_num_cells + i];
            const cellmodel_float_t Tn_CHc = states[STATE_Tn_CHc * padded_num_cells + i];
            const cellmodel_float_t Tn_CHm = states[STATE_Tn_CHm * padded_num_cells + i];
            const cellmodel_float_t CaM = states[STATE_CaM * padded_num_cells + i];
            const cellmodel_float_t Myo_c = states[STATE_Myo_c * padded_num_cells + i];
            const cellmodel_float_t Myo_m = states[STATE_Myo_m * padded_num_cells + i];
            const cellmodel_float_t SRB = states[STATE_SRB * padded_num_cells + i];
            const cellmodel_float_t SLL_j = states[STATE_SLL_j * padded_num_cells + i];
            const cellmodel_float_t SLL_sl = states[STATE_SLL_sl * padded_num_cells + i];
            const cellmodel_float_t SLH_j = states[STATE_SLH_j * padded_num_cells + i];
            const cellmodel_float_t SLH_sl = states[STATE_SLH_sl * padded_num_cells + i];
            const cellmodel_float_t Csqn_b = states[STATE_Csqn_b * padded_num_cells + i];
            const cellmodel_float_t Ca_sr = states[STATE_Ca_sr * padded_num_cells + i];
            const cellmodel_float_t Na_j = states[STATE_Na_j * padded_num_cells + i];
            const cellmodel_float_t Na_sl = states[STATE_Na_sl * padded_num_cells + i];
            const cellmodel_float_t Na_i = states[STATE_Na_i * padded_num_cells + i];
            const cellmodel_float_t K_i = states[STATE_K_i * padded_num_cells + i];
            const cellmodel_float_t Ca_j = states[STATE_Ca_j * padded_num_cells + i];
            const cellmodel_float_t Ca_sl = states[STATE_Ca_sl * padded_num_cells + i];
            const cellmodel_float_t Ca_i = states[STATE_Ca_i * padded_num_cells + i];
            const cellmodel_float_t V_m = states[STATE_V_m * padded_num_cells + i];

            // Expressions for the Geometry component
            const cellmodel_float_t Vcell =
                    FP_LITERAL(1.0e-15) * M_PI * cellLength * (cellRadius * cellRadius);
            const cellmodel_float_t Vmyo = FP_LITERAL(0.65) * Vcell;
            const cellmodel_float_t Vsr = FP_LITERAL(0.035) * Vcell;
            const cellmodel_float_t Vsl = FP_LITERAL(0.02) * Vcell;
            const cellmodel_float_t Vjunc = FP_LITERAL(0.000539) * Vcell;
            const cellmodel_float_t Fsl = FP_LITERAL(1.) - Fjunc;
            const cellmodel_float_t Fsl_CaL = FP_LITERAL(1.) - Fjunc_CaL;

            // Expressions for the Reversal potentials component
            const cellmodel_float_t FoRT = Frdy / (R * Temp);
            const cellmodel_float_t ena_junc = Log(Nao / Na_j) / FoRT;
            const cellmodel_float_t ena_sl = Log(Nao / Na_sl) / FoRT;
            const cellmodel_float_t ek = Log(Ko / K_i) / FoRT;
            const cellmodel_float_t eca_junc = Log(Cao / Ca_j) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t eca_sl = Log(Cao / Ca_sl) / (FP_LITERAL(2.) * FoRT);
            const cellmodel_float_t ecl = Log(Cli / Clo) / FoRT;
            const cellmodel_float_t Qpow = FP_LITERAL(-31.) + Temp / FP_LITERAL(10.);

            // Expressions for the I_Na component
            const cellmodel_float_t mss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(0.00184221158116513)
                                  * Exp(FP_LITERAL(-0.110741971207087) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(0.00184221158116513)
                                    * Exp(FP_LITERAL(-0.110741971207087) * V_m)));
            const cellmodel_float_t taum =
                    FP_LITERAL(0.1292)
                            * Exp(-((FP_LITERAL(2.94658944658945)
                                     + FP_LITERAL(0.0643500643500644) * V_m)
                                    * (FP_LITERAL(2.94658944658945)
                                       + FP_LITERAL(0.0643500643500644) * V_m)))
                    + FP_LITERAL(0.06487)
                              * Exp(-((FP_LITERAL(-0.0943466353677621)
                                       + FP_LITERAL(0.0195618153364632) * V_m)
                                      * (FP_LITERAL(-0.0943466353677621)
                                         + FP_LITERAL(0.0195618153364632) * V_m)));
            const cellmodel_float_t ah =
                    (V_m >= FP_LITERAL(-40.) ? FP_LITERAL(0.)
                                             : FP_LITERAL(4.43126792958051e-7)
                                                       * Exp(FP_LITERAL(-0.147058823529412) * V_m));
            const cellmodel_float_t bh =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.77)
                                       / (FP_LITERAL(0.13)
                                          + FP_LITERAL(0.0497581410839387)
                                                    * Exp(FP_LITERAL(-0.0900900900900901) * V_m))
                             : FP_LITERAL(310000.0) * Exp(FP_LITERAL(0.3485) * V_m)
                                       + FP_LITERAL(2.7) * Exp(FP_LITERAL(0.079) * V_m));
            const cellmodel_float_t tauh = FP_LITERAL(1.0) / (ah + bh);
            const cellmodel_float_t hss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(15212.5932856544)
                                    * Exp(FP_LITERAL(0.134589502018843) * V_m)));
            const cellmodel_float_t aj =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.)
                             : (FP_LITERAL(37.78) + V_m)
                                       * (FP_LITERAL(-25428.0) * Exp(FP_LITERAL(0.2444) * V_m)
                                          - FP_LITERAL(6.948e-6) * Exp(FP_LITERAL(-0.04391) * V_m))
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(50262745825.954)
                                                    * Exp(FP_LITERAL(0.311) * V_m)));
            const cellmodel_float_t bj =
                    (V_m >= FP_LITERAL(-40.)
                             ? FP_LITERAL(0.6) * Exp(FP_LITERAL(0.057) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.0407622039783662)
                                                    * Exp(FP_LITERAL(-0.1) * V_m))
                             : FP_LITERAL(0.02424) * Exp(FP_LITERAL(-0.01052) * V_m)
                                       / (FP_LITERAL(1.)
                                          + FP_LITERAL(0.00396086833990426)
                                                    * Exp(FP_LITERAL(-0.1378) * V_m)));
            const cellmodel_float_t tauj = FP_LITERAL(1.0) / (aj + bj);
            const cellmodel_float_t jss =
                    FP_LITERAL(1.0)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(15212.5932856544) * Exp(FP_LITERAL(0.134589502018843) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(15212.5932856544)
                                    * Exp(FP_LITERAL(0.134589502018843) * V_m)));
            const cellmodel_float_t dm_dt = (-m + mss) / taum;
            const cellmodel_float_t dm_dt_linearized = FP_LITERAL(-1.) / taum;
            states[STATE_m * padded_num_cells + i] =
                    (fabs(dm_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dm_dt_linearized)) * dm_dt
                                       / dm_dt_linearized
                             : dt * dm_dt)
                    + m;
            const cellmodel_float_t dh_dt = (-h + hss) / tauh;
            const cellmodel_float_t dh_dt_linearized = FP_LITERAL(-1.) / tauh;
            states[STATE_h * padded_num_cells + i] =
                    (fabs(dh_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dh_dt_linearized)) * dh_dt
                                       / dh_dt_linearized
                             : dt * dh_dt)
                    + h;
            const cellmodel_float_t dj_dt = (-j + jss) / tauj;
            const cellmodel_float_t dj_dt_linearized = FP_LITERAL(-1.) / tauj;
            states[STATE_j * padded_num_cells + i] =
                    (fabs(dj_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dj_dt_linearized)) * dj_dt
                                       / dj_dt_linearized
                             : dt * dj_dt)
                    + j;
            const cellmodel_float_t I_Na_junc =
                    Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j;
            const cellmodel_float_t I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j;

            // Expressions for the I_NaBK component
            const cellmodel_float_t I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m);
            const cellmodel_float_t I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl;

            // Expressions for the I_NaK component
            const cellmodel_float_t sigma =
                    FP_LITERAL(-1.) / FP_LITERAL(7.)
                    + Exp(FP_LITERAL(0.0148588410104012) * Nao) / FP_LITERAL(7.);
            const cellmodel_float_t fnak =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                       + FP_LITERAL(0.0365) * Exp(-FoRT * V_m) * sigma);
            const cellmodel_float_t I_nak_junc =
                    Fjunc * IbarNaK * Ko * fnak
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_j) * (Na_j)) * ((Na_j) * (Na_j))))
                       * (KmKo + Ko));
            const cellmodel_float_t I_nak_sl =
                    IbarNaK * Ko * Fsl * fnak
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_sl) * (Na_sl)) * ((Na_sl) * (Na_sl))))
                       * (KmKo + Ko));
            const cellmodel_float_t I_nak = I_nak_junc + I_nak_sl;

            // Expressions for the I_Kr component
            const cellmodel_float_t gkr = FP_LITERAL(0.430331482911935) * GKr * sqrt(Ko);
            const cellmodel_float_t xrss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + Exp(FP_LITERAL(-2.) - V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tauxr =
                    FP_LITERAL(230.)
                            / (FP_LITERAL(1.) + Exp(FP_LITERAL(2.) + V_m / FP_LITERAL(20.)))
                    + FP_LITERAL(3300.)
                              / ((FP_LITERAL(1.)
                                  + Exp(FP_LITERAL(-22.) / FP_LITERAL(9.) - V_m / FP_LITERAL(9.)))
                                 * (FP_LITERAL(1.)
                                    + Exp(FP_LITERAL(11.) / FP_LITERAL(9.)
                                          + V_m / FP_LITERAL(9.))));
            const cellmodel_float_t dx_kr_dt = (-x_kr + xrss) / tauxr;
            const cellmodel_float_t dx_kr_dt_linearized = FP_LITERAL(-1.) / tauxr;
            states[STATE_x_kr * padded_num_cells + i] =
                    (fabs(dx_kr_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dx_kr_dt_linearized)) * dx_kr_dt
                                       / dx_kr_dt_linearized
                             : dt * dx_kr_dt)
                    + x_kr;
            const cellmodel_float_t rkr =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(37.) / FP_LITERAL(12.) + V_m / FP_LITERAL(24.)));
            const cellmodel_float_t I_kr = (-ek + V_m) * gkr * rkr * x_kr;

            // Expressions for the I_Kp component
            const cellmodel_float_t kp_kp =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(1786.47556537862) * Exp(FP_LITERAL(-0.167224080267559) * V_m));
            const cellmodel_float_t I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp;
            const cellmodel_float_t I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp;
            const cellmodel_float_t I_kp = I_kp_junc + I_kp_sl;

            // Expressions for the I_Ks component
            const cellmodel_float_t eks = Log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT;
            const cellmodel_float_t gks_junc = GKs;
            const cellmodel_float_t gks_sl = GKs;
            const cellmodel_float_t xsss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.765928338364649)
                                 * Exp(FP_LITERAL(-0.0701754385964912) * V_m));
            const cellmodel_float_t tauxs =
                    FP_LITERAL(990.1)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.841540408868102)
                                 * Exp(FP_LITERAL(-0.0708215297450425) * V_m));
            const cellmodel_float_t dx_ks_dt = (-x_ks + xsss) / tauxs;
            const cellmodel_float_t dx_ks_dt_linearized = FP_LITERAL(-1.) / tauxs;
            states[STATE_x_ks * padded_num_cells + i] =
                    (fabs(dx_ks_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dx_ks_dt_linearized)) * dx_ks_dt
                                       / dx_ks_dt_linearized
                             : dt * dx_ks_dt)
                    + x_ks;
            const cellmodel_float_t I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc;
            const cellmodel_float_t I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl;
            const cellmodel_float_t I_ks = I_ks_junc + I_ks_sl;

            // Expressions for the I_to component
            const cellmodel_float_t GtoSlow =
                    (epi == FP_LITERAL(1.) ? FP_LITERAL(0.12) * Gto : FP_LITERAL(0.2892) * Gto);
            const cellmodel_float_t GtoFast =
                    (epi == FP_LITERAL(1.) ? FP_LITERAL(0.88) * Gto : FP_LITERAL(0.0108) * Gto);
            const cellmodel_float_t xtoss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(19.) / FP_LITERAL(13.) - V_m / FP_LITERAL(13.)));
            const cellmodel_float_t ytoss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.) + FP_LITERAL(49.4024491055302) * Exp(V_m / FP_LITERAL(5.)));
            const cellmodel_float_t tauxtos =
                    FP_LITERAL(0.5)
                    + FP_LITERAL(9.)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(1.) / FP_LITERAL(5.) + V_m / FP_LITERAL(15.)));
            const cellmodel_float_t tauytos =
                    FP_LITERAL(30.)
                    + FP_LITERAL(800.)
                              / (FP_LITERAL(1.) + Exp(FP_LITERAL(6.) + V_m / FP_LITERAL(10.)));
            const cellmodel_float_t dx_to_s_dt = (-x_to_s + xtoss) / tauxtos;
            const cellmodel_float_t dx_to_s_dt_linearized = FP_LITERAL(-1.) / tauxtos;
            states[STATE_x_to_s * padded_num_cells + i] =
                    (fabs(dx_to_s_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dx_to_s_dt_linearized)) * dx_to_s_dt
                                       / dx_to_s_dt_linearized
                             : dt * dx_to_s_dt)
                    + x_to_s;
            const cellmodel_float_t dy_to_s_dt = (-y_to_s + ytoss) / tauytos;
            const cellmodel_float_t dy_to_s_dt_linearized = FP_LITERAL(-1.) / tauytos;
            states[STATE_y_to_s * padded_num_cells + i] =
                    (fabs(dy_to_s_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dy_to_s_dt_linearized)) * dy_to_s_dt
                                       / dy_to_s_dt_linearized
                             : dt * dy_to_s_dt)
                    + y_to_s;
            const cellmodel_float_t I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s;
            const cellmodel_float_t tauxtof =
                    FP_LITERAL(0.5)
                    + FP_LITERAL(8.5)
                              * Exp(-((FP_LITERAL(9.) / FP_LITERAL(10.) + V_m / FP_LITERAL(50.))
                                      * (FP_LITERAL(9.) / FP_LITERAL(10.)
                                         + V_m / FP_LITERAL(50.))));
            const cellmodel_float_t tauytof =
                    FP_LITERAL(7.)
                    + FP_LITERAL(85.)
                              * Exp(-((FP_LITERAL(40.) + V_m) * (FP_LITERAL(40.) + V_m))
                                    / FP_LITERAL(220.));
            const cellmodel_float_t dx_to_f_dt = (-x_to_f + xtoss) / tauxtof;
            const cellmodel_float_t dx_to_f_dt_linearized = FP_LITERAL(-1.) / tauxtof;
            states[STATE_x_to_f * padded_num_cells + i] =
                    (fabs(dx_to_f_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dx_to_f_dt_linearized)) * dx_to_f_dt
                                       / dx_to_f_dt_linearized
                             : dt * dx_to_f_dt)
                    + x_to_f;
            const cellmodel_float_t dy_to_f_dt = (-y_to_f + ytoss) / tauytof;
            const cellmodel_float_t dy_to_f_dt_linearized = FP_LITERAL(-1.) / tauytof;
            states[STATE_y_to_f * padded_num_cells + i] =
                    (fabs(dy_to_f_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dy_to_f_dt_linearized)) * dy_to_f_dt
                                       / dy_to_f_dt_linearized
                             : dt * dy_to_f_dt)
                    + y_to_f;
            const cellmodel_float_t I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f;
            const cellmodel_float_t I_to = I_tof + I_tos;

            // Expressions for the I_K1 component
            const cellmodel_float_t aki =
                    FP_LITERAL(1.02)
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(7.35454251046446e-7)
                                 * Exp(FP_LITERAL(0.2385) * V_m - FP_LITERAL(0.2385) * ek));
            const cellmodel_float_t bki =
                    (FP_LITERAL(0.762624006506308)
                             * Exp(FP_LITERAL(0.08032) * V_m - FP_LITERAL(0.08032) * ek)
                     + FP_LITERAL(1.15340563518656e-16)
                               * Exp(FP_LITERAL(0.06175) * V_m - FP_LITERAL(0.06175) * ek))
                    / (FP_LITERAL(1.)
                       + FP_LITERAL(0.0867722941576933)
                                 * Exp(FP_LITERAL(0.5143) * ek - FP_LITERAL(0.5143) * V_m));
            const cellmodel_float_t kiss = aki / (aki + bki);
            const cellmodel_float_t I_K1 =
                    FP_LITERAL(0.430331482911935) * GK1 * sqrt(Ko) * (-ek + V_m) * kiss;

            // Expressions for the I_ClCa component
            const cellmodel_float_t I_ClCa_junc =
                    Fjunc * GClCa * (-ecl + V_m) / (FP_LITERAL(1.) + KdClCa / Ca_j);
            const cellmodel_float_t I_ClCa_sl =
                    GClCa * (-ecl + V_m) * Fsl / (FP_LITERAL(1.) + KdClCa / Ca_sl);
            const cellmodel_float_t I_ClCa = I_ClCa_junc + I_ClCa_sl;
            const cellmodel_float_t I_Clbk = GClB * (-ecl + V_m);

            // Expressions for the I_Ca component
            const cellmodel_float_t fss =
                    FP_LITERAL(1.0)
                            / (FP_LITERAL(1.)
                               + Exp(FP_LITERAL(35.) / FP_LITERAL(9.) + V_m / FP_LITERAL(9.)))
                    + FP_LITERAL(0.6)
                              / (FP_LITERAL(1.)
                                 + Exp(FP_LITERAL(5.) / FP_LITERAL(2.) - V_m / FP_LITERAL(20.)));
            const cellmodel_float_t dss =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(1.)
                       + Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)));
            const cellmodel_float_t taud =
                    (FP_LITERAL(1.) - Exp(FP_LITERAL(-5.) / FP_LITERAL(6.) - V_m / FP_LITERAL(6.)))
                    * dss / (FP_LITERAL(0.175) + FP_LITERAL(0.035) * V_m);
            const cellmodel_float_t tauf =
                    FP_LITERAL(1.0)
                    / (FP_LITERAL(0.02)
                       + FP_LITERAL(0.0197)
                                 * Exp(-((FP_LITERAL(0.48865) + FP_LITERAL(0.0337) * V_m)
                                         * (FP_LITERAL(0.48865) + FP_LITERAL(0.0337) * V_m))));
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
            const cellmodel_float_t df_Ca_Bj_dt =
                    FP_LITERAL(-0.0119) * f_Ca_Bj
                    + FP_LITERAL(1.7) * (FP_LITERAL(1.) - f_Ca_Bj) * Ca_j;
            const cellmodel_float_t df_Ca_Bj_dt_linearized =
                    FP_LITERAL(-0.0119) - FP_LITERAL(1.7) * Ca_j;
            states[STATE_f_Ca_Bj * padded_num_cells + i] =
                    (fabs(df_Ca_Bj_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * df_Ca_Bj_dt_linearized)) * df_Ca_Bj_dt
                                       / df_Ca_Bj_dt_linearized
                             : dt * df_Ca_Bj_dt)
                    + f_Ca_Bj;
            const cellmodel_float_t df_Ca_Bsl_dt =
                    FP_LITERAL(-0.0119) * f_Ca_Bsl
                    + FP_LITERAL(1.7) * (FP_LITERAL(1.) - f_Ca_Bsl) * Ca_sl;
            const cellmodel_float_t df_Ca_Bsl_dt_linearized =
                    FP_LITERAL(-0.0119) - FP_LITERAL(1.7) * Ca_sl;
            states[STATE_f_Ca_Bsl * padded_num_cells + i] =
                    (fabs(df_Ca_Bsl_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * df_Ca_Bsl_dt_linearized)) * df_Ca_Bsl_dt
                                       / df_Ca_Bsl_dt_linearized
                             : dt * df_Ca_Bsl_dt)
                    + f_Ca_Bsl;
            const cellmodel_float_t fcaCaMSL = FP_LITERAL(0.);
            const cellmodel_float_t fcaCaj = FP_LITERAL(0.);
            const cellmodel_float_t ibarca_j =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                    * (FP_LITERAL(-0.341) * Cao
                       + FP_LITERAL(0.341) * Ca_j * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t ibarca_sl =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                    * (FP_LITERAL(-0.341) * Cao
                       + FP_LITERAL(0.341) * Ca_sl * Exp(FP_LITERAL(2.) * FoRT * V_m))
                    * FoRT * V_m / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t ibark =
                    Frdy * GCaL * pK
                    * (FP_LITERAL(-0.75) * Ko + FP_LITERAL(0.75) * K_i * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t ibarna_j =
                    Frdy * GCaL * pNa
                    * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_j * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t ibarna_sl =
                    Frdy * GCaL * pNa
                    * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_sl * Exp(FoRT * V_m)) * FoRT
                    * V_m / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t I_Ca_junc = FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                                                * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f
                                                * ibarca_j;
            const cellmodel_float_t I_Ca_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                              * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d
                                              * f * ibarca_sl;
            const cellmodel_float_t I_CaK = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                            * (Fjunc_CaL * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj)
                                               + (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
                                            * d * f * ibark;
            const cellmodel_float_t I_CaNa_junc = FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                                                  * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f
                                                  * ibarna_j;
            const cellmodel_float_t I_CaNa_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                                * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL
                                                * d * f * ibarna_sl;

            // Expressions for the I_NCX component
            const cellmodel_float_t Ka_junc =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_j * Ca_j));
            const cellmodel_float_t Ka_sl =
                    FP_LITERAL(1.0) / (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_sl * Ca_sl));
            const cellmodel_float_t s1_junc = Cao * (Na_j * Na_j * Na_j) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t s2_junc =
                    (Nao * Nao * Nao) * Ca_j * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_junc =
                    Cao * (Na_j * Na_j * Na_j) + KmCao * (Na_j * Na_j * Na_j)
                    + (Nao * Nao * Nao) * Ca_j
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_j / KmCai) * Ca_j;
            const cellmodel_float_t s2_sl =
                    (Nao * Nao * Nao) * Ca_sl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t s3_sl =
                    Cao * (Na_sl * Na_sl * Na_sl) + KmCao * (Na_sl * Na_sl * Na_sl)
                    + (Nao * Nao * Nao) * Ca_sl
                    + KmCai * (Nao * Nao * Nao)
                              * (FP_LITERAL(1.) + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
                    + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_sl / KmCai) * Ca_sl;
            const cellmodel_float_t I_ncx_junc =
                    Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * s3_junc);
            const cellmodel_float_t I_ncx_sl =
                    IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);

            // Expressions for the I_PCa component
            const cellmodel_float_t I_pca_junc =
                    Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, FP_LITERAL(1.6))
                    / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_j, FP_LITERAL(1.6)));
            const cellmodel_float_t I_pca_sl =
                    IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, FP_LITERAL(1.6)) * Fsl
                    / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_sl, FP_LITERAL(1.6)));

            // Expressions for the I_CaBK component
            const cellmodel_float_t I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m);
            const cellmodel_float_t I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl;

            // Expressions for the SR Fluxes component
            const cellmodel_float_t kCaSR =
                    MaxSR
                    - (MaxSR - MinSR) / (FP_LITERAL(1.) + pow(ec50SR / Ca_sr, FP_LITERAL(2.5)));
            const cellmodel_float_t koSRCa = koCa / kCaSR;
            const cellmodel_float_t kiSRCa = kiCa * kCaSR;
            const cellmodel_float_t RI = FP_LITERAL(1.) - Ry_Ri - Ry_Ro - Ry_Rr;
            const cellmodel_float_t dRy_Rr_dt =
                    kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa;
            const cellmodel_float_t dRy_Rr_dt_linearized =
                    -kim - (Ca_j * Ca_j) * koSRCa - Ca_j * kiSRCa;
            states[STATE_Ry_Rr * padded_num_cells + i] =
                    (fabs(dRy_Rr_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dRy_Rr_dt_linearized)) * dRy_Rr_dt
                                       / dRy_Rr_dt_linearized
                             : dt * dRy_Rr_dt)
                    + Ry_Rr;
            const cellmodel_float_t dRy_Ro_dt = kim * Ry_Ri - kom * Ry_Ro
                                                + (Ca_j * Ca_j) * Ry_Rr * koSRCa
                                                - Ca_j * Ry_Ro * kiSRCa;
            const cellmodel_float_t dRy_Ro_dt_linearized = -kom - Ca_j * kiSRCa;
            states[STATE_Ry_Ro * padded_num_cells + i] =
                    (fabs(dRy_Ro_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dRy_Ro_dt_linearized)) * dRy_Ro_dt
                                       / dRy_Ro_dt_linearized
                             : dt * dRy_Ro_dt)
                    + Ry_Ro;
            const cellmodel_float_t dRy_Ri_dt = -kim * Ry_Ri - kom * Ry_Ri
                                                + (Ca_j * Ca_j) * RI * koSRCa
                                                + Ca_j * Ry_Ro * kiSRCa;
            const cellmodel_float_t dRy_Ri_dt_linearized = -kim - kom - (Ca_j * Ca_j) * koSRCa;
            states[STATE_Ry_Ri * padded_num_cells + i] =
                    (fabs(dRy_Ri_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dRy_Ri_dt_linearized)) * dRy_Ri_dt
                                       / dRy_Ri_dt_linearized
                             : dt * dRy_Ri_dt)
                    + Ry_Ri;
            const cellmodel_float_t J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro;
            const cellmodel_float_t J_serca =
                    Vmax_SRCaP * pow(Q10SRCaP, Qpow)
                    * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                    / (FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP));
            const cellmodel_float_t J_SRleak =
                    FP_LITERAL(5.348e-6) * Ca_sr - FP_LITERAL(5.348e-6) * Ca_j;

            // Expressions for the Na Buffers component
            const cellmodel_float_t dNa_Bj_dt =
                    -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j;
            const cellmodel_float_t dNa_Bj_dt_linearized = -koff_na - kon_na * Na_j;
            states[STATE_Na_Bj * padded_num_cells + i] =
                    Na_Bj
                    + (fabs(dNa_Bj_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_Bj_dt_linearized)) * dNa_Bj_dt
                                         / dNa_Bj_dt_linearized
                               : dt * dNa_Bj_dt);
            const cellmodel_float_t dNa_Bsl_dt =
                    -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl;
            const cellmodel_float_t dNa_Bsl_dt_linearized = -koff_na - kon_na * Na_sl;
            states[STATE_Na_Bsl * padded_num_cells + i] =
                    Na_Bsl
                    + (fabs(dNa_Bsl_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_Bsl_dt_linearized)) * dNa_Bsl_dt
                                         / dNa_Bsl_dt_linearized
                               : dt * dNa_Bsl_dt);

            // Expressions for the Cytosolic Ca Buffers component
            const cellmodel_float_t dTn_CL_dt =
                    -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;
            const cellmodel_float_t dTn_CL_dt_linearized = -koff_tncl - kon_tncl * Ca_i;
            states[STATE_Tn_CL * padded_num_cells + i] =
                    (fabs(dTn_CL_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dTn_CL_dt_linearized)) * dTn_CL_dt
                                       / dTn_CL_dt_linearized
                             : dt * dTn_CL_dt)
                    + Tn_CL;
            const cellmodel_float_t dTn_CHc_dt =
                    -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i;
            const cellmodel_float_t dTn_CHc_dt_linearized = -koff_tnchca - kon_tnchca * Ca_i;
            states[STATE_Tn_CHc * padded_num_cells + i] =
                    (fabs(dTn_CHc_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dTn_CHc_dt_linearized)) * dTn_CHc_dt
                                       / dTn_CHc_dt_linearized
                             : dt * dTn_CHc_dt)
                    + Tn_CHc;
            const cellmodel_float_t dTn_CHm_dt =
                    -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm);
            const cellmodel_float_t dTn_CHm_dt_linearized = -koff_tnchmg - Mgi * kon_tnchmg;
            states[STATE_Tn_CHm * padded_num_cells + i] =
                    (fabs(dTn_CHm_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dTn_CHm_dt_linearized)) * dTn_CHm_dt
                                       / dTn_CHm_dt_linearized
                             : dt * dTn_CHm_dt)
                    + Tn_CHm;
            const cellmodel_float_t dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i;
            const cellmodel_float_t dCaM_dt_linearized = -koff_cam - kon_cam * Ca_i;
            states[STATE_CaM * padded_num_cells + i] =
                    CaM
                    + (fabs(dCaM_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCaM_dt_linearized)) * dCaM_dt
                                         / dCaM_dt_linearized
                               : dt * dCaM_dt);
            const cellmodel_float_t dMyo_c_dt =
                    -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i;
            const cellmodel_float_t dMyo_c_dt_linearized = -koff_myoca - kon_myoca * Ca_i;
            states[STATE_Myo_c * padded_num_cells + i] =
                    Myo_c
                    + (fabs(dMyo_c_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dMyo_c_dt_linearized)) * dMyo_c_dt
                                         / dMyo_c_dt_linearized
                               : dt * dMyo_c_dt);
            const cellmodel_float_t dMyo_m_dt =
                    -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m);
            const cellmodel_float_t dMyo_m_dt_linearized = -koff_myomg - Mgi * kon_myomg;
            states[STATE_Myo_m * padded_num_cells + i] =
                    Myo_m
                    + (fabs(dMyo_m_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dMyo_m_dt_linearized)) * dMyo_m_dt
                                         / dMyo_m_dt_linearized
                               : dt * dMyo_m_dt);
            const cellmodel_float_t dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i;
            const cellmodel_float_t dSRB_dt_linearized = -koff_sr - kon_sr * Ca_i;
            states[STATE_SRB * padded_num_cells + i] =
                    (fabs(dSRB_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dSRB_dt_linearized)) * dSRB_dt
                                       / dSRB_dt_linearized
                             : dt * dSRB_dt)
                    + SRB;
            const cellmodel_float_t J_CaB_cytosol =
                    -koff_cam * CaM - koff_myoca * Myo_c - koff_myomg * Myo_m - koff_sr * SRB
                    - koff_tnchca * Tn_CHc - koff_tnchmg * Tn_CHm - koff_tncl * Tn_CL
                    + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
                    + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
                    + kon_cam * (Bmax_CaM - CaM) * Ca_i
                    + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
                    + kon_sr * (Bmax_SR - SRB) * Ca_i
                    + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
                    + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;

            // Expressions for the Junctional and SL Ca Buffers component
            const cellmodel_float_t Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl;
            const cellmodel_float_t Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc;
            const cellmodel_float_t Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl;
            const cellmodel_float_t Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc;
            const cellmodel_float_t dSLL_j_dt =
                    -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
            const cellmodel_float_t dSLL_j_dt_linearized = -koff_sll - kon_sll * Ca_j;
            states[STATE_SLL_j * padded_num_cells + i] =
                    (fabs(dSLL_j_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dSLL_j_dt_linearized)) * dSLL_j_dt
                                       / dSLL_j_dt_linearized
                             : dt * dSLL_j_dt)
                    + SLL_j;
            const cellmodel_float_t dSLL_sl_dt =
                    -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;
            const cellmodel_float_t dSLL_sl_dt_linearized = -koff_sll - kon_sll * Ca_sl;
            states[STATE_SLL_sl * padded_num_cells + i] =
                    (fabs(dSLL_sl_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dSLL_sl_dt_linearized)) * dSLL_sl_dt
                                       / dSLL_sl_dt_linearized
                             : dt * dSLL_sl_dt)
                    + SLL_sl;
            const cellmodel_float_t dSLH_j_dt =
                    -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j;
            const cellmodel_float_t dSLH_j_dt_linearized = -koff_slh - kon_slh * Ca_j;
            states[STATE_SLH_j * padded_num_cells + i] =
                    (fabs(dSLH_j_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dSLH_j_dt_linearized)) * dSLH_j_dt
                                       / dSLH_j_dt_linearized
                             : dt * dSLH_j_dt)
                    + SLH_j;
            const cellmodel_float_t dSLH_sl_dt =
                    -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl;
            const cellmodel_float_t dSLH_sl_dt_linearized = -koff_slh - kon_slh * Ca_sl;
            states[STATE_SLH_sl * padded_num_cells + i] =
                    (fabs(dSLH_sl_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dSLH_sl_dt_linearized)) * dSLH_sl_dt
                                       / dSLH_sl_dt_linearized
                             : dt * dSLH_sl_dt)
                    + SLH_sl;
            const cellmodel_float_t J_CaB_junction = -koff_slh * SLH_j - koff_sll * SLL_j
                                                     + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
                                                     + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
            const cellmodel_float_t J_CaB_sl = -koff_slh * SLH_sl - koff_sll * SLL_sl
                                               + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
                                               + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;

            // Expressions for the SR Ca Concentrations component
            const cellmodel_float_t Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr;
            const cellmodel_float_t dCsqn_b_dt =
                    -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr;
            const cellmodel_float_t dCsqn_b_dt_linearized = -koff_csqn - kon_csqn * Ca_sr;
            states[STATE_Csqn_b * padded_num_cells + i] =
                    Csqn_b
                    + (fabs(dCsqn_b_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCsqn_b_dt_linearized)) * dCsqn_b_dt
                                         / dCsqn_b_dt_linearized
                               : dt * dCsqn_b_dt);
            const cellmodel_float_t dCa_sr_dt = -J_SRCarel + koff_csqn * Csqn_b
                                                - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr
                                                - J_SRleak * Vmyo / Vsr + J_serca;
            const cellmodel_float_t dJ_serca_dCa_sr =
                    -Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_sr / Kmr, hillSRCaP)
                            / ((FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                + pow(Ca_sr / Kmr, hillSRCaP))
                               * Ca_sr)
                    - Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_sr / Kmr, hillSRCaP)
                              * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                              / (((FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                   + pow(Ca_sr / Kmr, hillSRCaP))
                                  * (FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                     + pow(Ca_sr / Kmr, hillSRCaP)))
                                 * Ca_sr);
            const cellmodel_float_t dCa_sr_dt_linearized =
                    -kon_csqn * (-Csqn_b + Bmax_Csqn) - ks * Ry_Ro
                    - FP_LITERAL(5.348e-6) * Vmyo / Vsr + dJ_serca_dCa_sr;
            states[STATE_Ca_sr * padded_num_cells + i] =
                    Ca_sr
                    + (fabs(dCa_sr_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCa_sr_dt_linearized)) * dCa_sr_dt
                                         / dCa_sr_dt_linearized
                               : dt * dCa_sr_dt);

            // Expressions for the Na Concentrations component
            const cellmodel_float_t I_Na_tot_junc = FP_LITERAL(3.) * I_nak_junc
                                                    + FP_LITERAL(3.) * I_ncx_junc + I_CaNa_junc
                                                    + I_Na_junc + I_nabk_junc;
            const cellmodel_float_t I_Na_tot_sl = FP_LITERAL(3.) * I_nak_sl
                                                  + FP_LITERAL(3.) * I_ncx_sl + I_CaNa_sl + I_Na_sl
                                                  + I_nabk_sl;
            const cellmodel_float_t dNa_j_dt = -dNa_Bj_dt + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
                                               - Cmem * I_Na_tot_junc / (Frdy * Vjunc);
            const cellmodel_float_t dI_ncx_junc_ds3_junc =
                    -Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * (s3_junc * s3_junc));
            const cellmodel_float_t dibarna_j_dNa_j = FP_LITERAL(0.75) * Frdy * GCaL * pNa * FoRT
                                                      * V_m * Exp(FoRT * V_m)
                                                      / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t dI_Na_junc_dena_junc = -Fjunc * GNa * (m * m * m) * h * j;
            const cellmodel_float_t dI_ncx_junc_ds1_junc =
                    Fjunc * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * s3_junc);
            const cellmodel_float_t ds1_junc_dNa_j =
                    FP_LITERAL(3.) * Cao * (Na_j * Na_j) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t dI_nak_junc_dNa_j =
                    FP_LITERAL(4.) * Fjunc * IbarNaK * Ko
                    * (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip))) * fnak
                    / (((1.
                         + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                   / (((Na_j) * (Na_j)) * ((Na_j) * (Na_j))))
                        * (1.
                           + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                     / (((Na_j) * (Na_j)) * ((Na_j) * (Na_j)))))
                       * (KmKo + Ko) * pow(Na_j, FP_LITERAL(5.)));
            const cellmodel_float_t ds3_junc_dNa_j =
                    FP_LITERAL(3.) * Cao * (Na_j * Na_j) + FP_LITERAL(3.) * KmCao * (Na_j * Na_j)
                    + FP_LITERAL(3.) * KmCai * (Nao * Nao * Nao) * (Na_j * Na_j)
                              / (KmNai * KmNai * KmNai);
            const cellmodel_float_t dena_junc_dNa_j = FP_LITERAL(-1.) / (FoRT * Na_j);
            const cellmodel_float_t dI_nabk_junc_dena_junc = -Fjunc * GNaB;
            const cellmodel_float_t dI_CaNa_junc_dibarna_j =
                    FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                    * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f;
            const cellmodel_float_t dNa_j_dt_linearized =
                    -J_na_juncsl / Vjunc
                    - Cmem
                              * (FP_LITERAL(3.) * dI_nak_junc_dNa_j
                                 + dI_CaNa_junc_dibarna_j * dibarna_j_dNa_j
                                 + dI_Na_junc_dena_junc * dena_junc_dNa_j
                                 + dI_nabk_junc_dena_junc * dena_junc_dNa_j
                                 + FP_LITERAL(3.) * dI_ncx_junc_ds1_junc * ds1_junc_dNa_j
                                 + FP_LITERAL(3.) * dI_ncx_junc_ds3_junc * ds3_junc_dNa_j)
                              / (Frdy * Vjunc);
            states[STATE_Na_j * padded_num_cells + i] =
                    Na_j
                    + (fabs(dNa_j_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_j_dt_linearized)) * dNa_j_dt
                                         / dNa_j_dt_linearized
                               : dt * dNa_j_dt);
            const cellmodel_float_t dNa_sl_dt = -dNa_Bsl_dt + J_na_juncsl * (-Na_sl + Na_j) / Vsl
                                                + J_na_slmyo * (-Na_sl + Na_i) / Vsl
                                                - Cmem * I_Na_tot_sl / (Frdy * Vsl);
            const cellmodel_float_t dI_CaNa_sl_dibarna_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                                            * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl)
                                                            * Fsl_CaL * d * f;
            const cellmodel_float_t dI_nabk_sl_dena_sl = -GNaB * Fsl;
            const cellmodel_float_t ds3_sl_dNa_sl = FP_LITERAL(3.) * Cao * (Na_sl * Na_sl)
                                                    + FP_LITERAL(3.) * KmCao * (Na_sl * Na_sl)
                                                    + FP_LITERAL(3.) * KmCai * (Nao * Nao * Nao)
                                                              * (Na_sl * Na_sl)
                                                              / (KmNai * KmNai * KmNai);
            const cellmodel_float_t ds1_sl_dNa_sl =
                    FP_LITERAL(3.) * Cao * (Na_sl * Na_sl) * Exp(nu * FoRT * V_m);
            const cellmodel_float_t dI_ncx_sl_ds3_sl =
                    -IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * (s3_sl * s3_sl));
            const cellmodel_float_t dena_sl_dNa_sl = FP_LITERAL(-1.) / (FoRT * Na_sl);
            const cellmodel_float_t dI_ncx_sl_ds1_sl =
                    IbarNCX * pow(Q10NCX, Qpow) * Fsl * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t dI_Na_sl_dena_sl = -GNa * (m * m * m) * Fsl * h * j;
            const cellmodel_float_t dibarna_sl_dNa_sl = FP_LITERAL(0.75) * Frdy * GCaL * pNa * FoRT
                                                        * V_m * Exp(FoRT * V_m)
                                                        / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t dI_nak_sl_dNa_sl =
                    FP_LITERAL(4.) * IbarNaK * Ko * (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                    * Fsl * fnak
                    / (((1.
                         + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                   / (((Na_sl) * (Na_sl)) * ((Na_sl) * (Na_sl))))
                        * (1.
                           + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                     / (((Na_sl) * (Na_sl)) * ((Na_sl) * (Na_sl)))))
                       * (KmKo + Ko) * pow(Na_sl, FP_LITERAL(5.)));
            const cellmodel_float_t dNa_sl_dt_linearized =
                    -J_na_juncsl / Vsl - J_na_slmyo / Vsl
                    - Cmem
                              * (FP_LITERAL(3.) * dI_nak_sl_dNa_sl
                                 + dI_CaNa_sl_dibarna_sl * dibarna_sl_dNa_sl
                                 + dI_Na_sl_dena_sl * dena_sl_dNa_sl
                                 + dI_nabk_sl_dena_sl * dena_sl_dNa_sl
                                 + FP_LITERAL(3.) * dI_ncx_sl_ds1_sl * ds1_sl_dNa_sl
                                 + FP_LITERAL(3.) * dI_ncx_sl_ds3_sl * ds3_sl_dNa_sl)
                              / (Frdy * Vsl);
            states[STATE_Na_sl * padded_num_cells + i] =
                    Na_sl
                    + (fabs(dNa_sl_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_sl_dt_linearized)) * dNa_sl_dt
                                         / dNa_sl_dt_linearized
                               : dt * dNa_sl_dt);
            const cellmodel_float_t dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo;
            const cellmodel_float_t dNa_i_dt_linearized = -J_na_slmyo / Vmyo;
            states[STATE_Na_i * padded_num_cells + i] =
                    Na_i
                    + (fabs(dNa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dNa_i_dt_linearized)) * dNa_i_dt
                                         / dNa_i_dt_linearized
                               : dt * dNa_i_dt);

            // Expressions for the K Concentration component
            const cellmodel_float_t I_K_tot =
                    FP_LITERAL(-2.) * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to;
            const cellmodel_float_t dK_i_dt = FP_LITERAL(0.);
            states[STATE_K_i * padded_num_cells + i] = dt * dK_i_dt + K_i;

            // Expressions for the Ca Concentrations component
            const cellmodel_float_t I_Ca_tot_junc =
                    FP_LITERAL(-2.) * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
            const cellmodel_float_t I_Ca_tot_sl =
                    FP_LITERAL(-2.) * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
            const cellmodel_float_t dCa_j_dt =
                    -J_CaB_junction + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
                    + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc
                    - Cmem * I_Ca_tot_junc / (FP_LITERAL(2.) * Frdy * Vjunc);
            const cellmodel_float_t ds3_junc_dCa_j =
                    (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_j / KmCai)
                    + (KmNao * KmNao * KmNao) * Ca_j / KmCai;
            const cellmodel_float_t dJ_SRCarel_dCa_j = -ks * Ry_Ro;
            const cellmodel_float_t dI_cabk_junc_deca_junc = -Fjunc * GCaB;
            const cellmodel_float_t dI_pca_junc_dCa_j =
                    FP_LITERAL(1.6) * Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow)
                            * pow(Ca_j, FP_LITERAL(0.6))
                            / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_j, FP_LITERAL(1.6)))
                    - FP_LITERAL(1.6) * Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow)
                              * pow(Ca_j, FP_LITERAL(2.2))
                              / ((pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_j, FP_LITERAL(1.6)))
                                 * (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_j, FP_LITERAL(1.6))));
            const cellmodel_float_t dI_Ca_junc_dibarca_j =
                    FP_LITERAL(0.45) * Fjunc_CaL * pow(Q10CaL, Qpow)
                    * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj) * d * f;
            const cellmodel_float_t dI_ncx_junc_ds2_junc =
                    -Fjunc * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * s3_junc);
            const cellmodel_float_t dJ_CaB_junction_dCa_j =
                    kon_slh * (-SLH_j + Bmax_SLhighj) + kon_sll * (-SLL_j + Bmax_SLlowj);
            const cellmodel_float_t dI_ncx_junc_dKa_junc =
                    Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc)
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                       * s3_junc);
            const cellmodel_float_t dKa_junc_dCa_j =
                    FP_LITERAL(2.) * (Kdact * Kdact)
                    / (((FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_j * Ca_j))
                        * (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_j * Ca_j)))
                       * (Ca_j * Ca_j * Ca_j));
            const cellmodel_float_t dibarca_j_dCa_j =
                    FP_LITERAL(1.364) * Frdy * GCaL * pCa * FoRT * V_m
                    * Exp(FP_LITERAL(2.) * FoRT * V_m)
                    / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t deca_junc_dCa_j =
                    FP_LITERAL(-1.) / (FP_LITERAL(2.) * Ca_j * FoRT);
            const cellmodel_float_t ds2_junc_dCa_j =
                    (Nao * Nao * Nao) * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t dCa_j_dt_linearized =
                    -dJ_CaB_junction_dCa_j - J_ca_juncsl / Vjunc
                    - FP_LITERAL(5.348e-6) * Vmyo / Vjunc + Vsr * dJ_SRCarel_dCa_j / Vjunc
                    - Cmem
                              * (dI_Ca_junc_dibarca_j * dibarca_j_dCa_j
                                 + dI_cabk_junc_deca_junc * deca_junc_dCa_j
                                 - FP_LITERAL(2.) * dI_ncx_junc_dKa_junc * dKa_junc_dCa_j
                                 - FP_LITERAL(2.) * dI_ncx_junc_ds2_junc * ds2_junc_dCa_j
                                 - FP_LITERAL(2.) * dI_ncx_junc_ds3_junc * ds3_junc_dCa_j
                                 + dI_pca_junc_dCa_j)
                              / (FP_LITERAL(2.) * Frdy * Vjunc);
            states[STATE_Ca_j * padded_num_cells + i] =
                    Ca_j
                    + (fabs(dCa_j_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCa_j_dt_linearized)) * dCa_j_dt
                                         / dCa_j_dt_linearized
                               : dt * dCa_j_dt);
            const cellmodel_float_t dCa_sl_dt =
                    -J_CaB_sl + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
                    + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
                    - Cmem * I_Ca_tot_sl / (FP_LITERAL(2.) * Frdy * Vsl);
            const cellmodel_float_t dKa_sl_dCa_sl =
                    FP_LITERAL(2.) * (Kdact * Kdact)
                    / (((FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_sl * Ca_sl))
                        * (FP_LITERAL(1.) + (Kdact * Kdact) / (Ca_sl * Ca_sl)))
                       * (Ca_sl * Ca_sl * Ca_sl));
            const cellmodel_float_t dI_ncx_sl_dKa_sl =
                    IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t dibarca_sl_dCa_sl =
                    FP_LITERAL(1.364) * Frdy * GCaL * pCa * FoRT * V_m
                    * Exp(FP_LITERAL(2.) * FoRT * V_m)
                    / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t dI_ncx_sl_ds2_sl =
                    -IbarNCX * pow(Q10NCX, Qpow) * Fsl * Ka_sl
                    / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)) * s3_sl);
            const cellmodel_float_t ds2_sl_dCa_sl =
                    (Nao * Nao * Nao) * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t dI_cabk_sl_deca_sl = -GCaB * Fsl;
            const cellmodel_float_t deca_sl_dCa_sl =
                    FP_LITERAL(-1.) / (FP_LITERAL(2.) * Ca_sl * FoRT);
            const cellmodel_float_t dJ_CaB_sl_dCa_sl =
                    kon_slh * (-SLH_sl + Bmax_SLhighsl) + kon_sll * (-SLL_sl + Bmax_SLlowsl);
            const cellmodel_float_t dI_Ca_sl_dibarca_sl = FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                                                          * (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl)
                                                          * Fsl_CaL * d * f;
            const cellmodel_float_t dI_pca_sl_dCa_sl =
                    FP_LITERAL(1.6) * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, FP_LITERAL(0.6))
                            * Fsl / (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_sl, FP_LITERAL(1.6)))
                    - FP_LITERAL(1.6) * IbarSLCaP * pow(Q10SLCaP, Qpow)
                              * pow(Ca_sl, FP_LITERAL(2.2)) * Fsl
                              / ((pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_sl, FP_LITERAL(1.6)))
                                 * (pow(KmPCa, FP_LITERAL(1.6)) + pow(Ca_sl, FP_LITERAL(1.6))));
            const cellmodel_float_t ds3_sl_dCa_sl =
                    (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (FP_LITERAL(1.) + Ca_sl / KmCai)
                    + (KmNao * KmNao * KmNao) * Ca_sl / KmCai;
            const cellmodel_float_t dCa_sl_dt_linearized =
                    -dJ_CaB_sl_dCa_sl - J_ca_juncsl / Vsl - J_ca_slmyo / Vsl
                    - Cmem
                              * (dI_Ca_sl_dibarca_sl * dibarca_sl_dCa_sl
                                 + dI_cabk_sl_deca_sl * deca_sl_dCa_sl
                                 - FP_LITERAL(2.) * dI_ncx_sl_dKa_sl * dKa_sl_dCa_sl
                                 - FP_LITERAL(2.) * dI_ncx_sl_ds2_sl * ds2_sl_dCa_sl
                                 - FP_LITERAL(2.) * dI_ncx_sl_ds3_sl * ds3_sl_dCa_sl
                                 + dI_pca_sl_dCa_sl)
                              / (FP_LITERAL(2.) * Frdy * Vsl);
            states[STATE_Ca_sl * padded_num_cells + i] =
                    Ca_sl
                    + (fabs(dCa_sl_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCa_sl_dt_linearized)) * dCa_sl_dt
                                         / dCa_sl_dt_linearized
                               : dt * dCa_sl_dt);
            const cellmodel_float_t dCa_i_dt =
                    -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo;
            const cellmodel_float_t dJ_CaB_cytosol_dCa_i =
                    kon_cam * (Bmax_CaM - CaM) + kon_myoca * (Bmax_myosin - Myo_c - Myo_m)
                    + kon_sr * (Bmax_SR - SRB) + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
                    + kon_tncl * (Bmax_TnClow - Tn_CL);
            const cellmodel_float_t dJ_serca_dCa_i =
                    Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_i / Kmf, hillSRCaP)
                            / ((FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                + pow(Ca_sr / Kmr, hillSRCaP))
                               * Ca_i)
                    - Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_i / Kmf, hillSRCaP)
                              * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                              / (((FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                   + pow(Ca_sr / Kmr, hillSRCaP))
                                  * (FP_LITERAL(1.) + pow(Ca_i / Kmf, hillSRCaP)
                                     + pow(Ca_sr / Kmr, hillSRCaP)))
                                 * Ca_i);
            const cellmodel_float_t dCa_i_dt_linearized =
                    -dJ_CaB_cytosol_dCa_i - J_ca_slmyo / Vmyo - Vsr * dJ_serca_dCa_i / Vmyo;
            states[STATE_Ca_i * padded_num_cells + i] =
                    Ca_i
                    + (fabs(dCa_i_dt_linearized) > FP_LITERAL(1.0e-8)
                               ? (FP_LITERAL(-1.0) + Exp(dt * dCa_i_dt_linearized)) * dCa_i_dt
                                         / dCa_i_dt_linearized
                               : dt * dCa_i_dt);

            // Expressions for the Membrane potential component
            const cellmodel_float_t i_Stim =
                    (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                     && t - stim_period * floor(t / stim_period) >= stim_start
                             ? -stim_amplitude
                             : FP_LITERAL(0.));
            const cellmodel_float_t I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
            const cellmodel_float_t I_Cl_tot = I_ClCa + I_Clbk;
            const cellmodel_float_t I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
            const cellmodel_float_t I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot;
            const cellmodel_float_t dV_m_dt = -I_tot - i_Stim;
            const cellmodel_float_t daki_dV_m =
                    FP_LITERAL(-1.78913955652069e-7)
                    * Exp(FP_LITERAL(0.2385) * V_m - FP_LITERAL(0.2385) * ek)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(7.35454251046446e-7)
                                  * Exp(FP_LITERAL(0.2385) * V_m - FP_LITERAL(0.2385) * ek))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(7.35454251046446e-7)
                                    * Exp(FP_LITERAL(0.2385) * V_m - FP_LITERAL(0.2385) * ek)));
            const cellmodel_float_t dI_nak_junc_dfnak =
                    Fjunc * IbarNaK * Ko
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_j) * (Na_j)) * ((Na_j) * (Na_j))))
                       * (KmKo + Ko));
            const cellmodel_float_t dibark_dV_m =
                    Frdy * GCaL * pK
                            * (FP_LITERAL(-0.75) * Ko + FP_LITERAL(0.75) * K_i * Exp(FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FoRT * V_m))
                    - Frdy * GCaL * pK * (FoRT * FoRT)
                              * (FP_LITERAL(-0.75) * Ko + FP_LITERAL(0.75) * K_i * Exp(FoRT * V_m))
                              * V_m * Exp(FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FoRT * V_m)))
                    + FP_LITERAL(0.75) * Frdy * GCaL * pK * (FoRT * FoRT) * K_i * V_m
                              * Exp(FoRT * V_m) / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t dI_CaK_dibark =
                    FP_LITERAL(0.45) * pow(Q10CaL, Qpow)
                    * (Fjunc_CaL * (FP_LITERAL(1.) + fcaCaj - f_Ca_Bj)
                       + (FP_LITERAL(1.) + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL)
                    * d * f;
            const cellmodel_float_t ds2_sl_dV_m = (Nao * Nao * Nao) * (FP_LITERAL(-1.) + nu) * Ca_sl
                                                  * FoRT * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t dI_kr_drkr = (-ek + V_m) * gkr * x_kr;
            const cellmodel_float_t dbki_dV_m =
                    (FP_LITERAL(7.12227979727698e-18)
                             * Exp(FP_LITERAL(0.06175) * V_m - FP_LITERAL(0.06175) * ek)
                     + FP_LITERAL(0.0612539602025867)
                               * Exp(FP_LITERAL(0.08032) * V_m - FP_LITERAL(0.08032) * ek))
                            / (FP_LITERAL(1.)
                               + FP_LITERAL(0.0867722941576933)
                                         * Exp(FP_LITERAL(0.5143) * ek - FP_LITERAL(0.5143) * V_m))
                    + FP_LITERAL(0.0446269908853017)
                              * (FP_LITERAL(0.762624006506308)
                                         * Exp(FP_LITERAL(0.08032) * V_m - FP_LITERAL(0.08032) * ek)
                                 + FP_LITERAL(1.15340563518656e-16)
                                           * Exp(FP_LITERAL(0.06175) * V_m
                                                 - FP_LITERAL(0.06175) * ek))
                              * Exp(FP_LITERAL(0.5143) * ek - FP_LITERAL(0.5143) * V_m)
                              / ((FP_LITERAL(1.)
                                  + FP_LITERAL(0.0867722941576933)
                                            * Exp(FP_LITERAL(0.5143) * ek
                                                  - FP_LITERAL(0.5143) * V_m))
                                 * (FP_LITERAL(1.)
                                    + FP_LITERAL(0.0867722941576933)
                                              * Exp(FP_LITERAL(0.5143) * ek
                                                    - FP_LITERAL(0.5143) * V_m)));
            const cellmodel_float_t dkiss_daki =
                    FP_LITERAL(1.0) / (aki + bki) - aki / ((aki + bki) * (aki + bki));
            const cellmodel_float_t dkiss_dbki = -aki / ((aki + bki) * (aki + bki));
            const cellmodel_float_t dI_K1_dV_m =
                    FP_LITERAL(0.430331482911935) * GK1 * sqrt(Ko) * kiss
                    + FP_LITERAL(0.430331482911935) * GK1 * sqrt(Ko) * (-ek + V_m)
                              * (daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki);
            const cellmodel_float_t dI_ks_sl_dV_m = (x_ks * x_ks) * Fsl * gks_sl;
            const cellmodel_float_t dibarca_sl_dV_m =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                            * (FP_LITERAL(-0.341) * Cao
                               + FP_LITERAL(0.341) * Ca_sl * Exp(FP_LITERAL(2.) * FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                    - FP_LITERAL(8.) * Frdy * GCaL * pCa * (FoRT * FoRT)
                              * (FP_LITERAL(-0.341) * Cao
                                 + FP_LITERAL(0.341) * Ca_sl * Exp(FP_LITERAL(2.) * FoRT * V_m))
                              * V_m * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m)))
                    + FP_LITERAL(2.728) * Frdy * GCaL * pCa * (FoRT * FoRT) * Ca_sl * V_m
                              * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t dibarca_j_dV_m =
                    FP_LITERAL(4.) * Frdy * GCaL * pCa
                            * (FP_LITERAL(-0.341) * Cao
                               + FP_LITERAL(0.341) * Ca_j * Exp(FP_LITERAL(2.) * FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                    - FP_LITERAL(8.) * Frdy * GCaL * pCa * (FoRT * FoRT)
                              * (FP_LITERAL(-0.341) * Cao
                                 + FP_LITERAL(0.341) * Ca_j * Exp(FP_LITERAL(2.) * FoRT * V_m))
                              * V_m * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m)))
                    + FP_LITERAL(2.728) * Frdy * GCaL * pCa * (FoRT * FoRT) * Ca_j * V_m
                              * Exp(FP_LITERAL(2.) * FoRT * V_m)
                              / (FP_LITERAL(-1.) + Exp(FP_LITERAL(2.) * FoRT * V_m));
            const cellmodel_float_t dI_ks_junc_dV_m = Fjunc * (x_ks * x_ks) * gks_junc;
            const cellmodel_float_t ds2_junc_dV_m = (Nao * Nao * Nao) * (FP_LITERAL(-1.) + nu)
                                                    * Ca_j * FoRT
                                                    * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m);
            const cellmodel_float_t ds1_sl_dV_m =
                    Cao * nu * (Na_sl * Na_sl * Na_sl) * FoRT * Exp(nu * FoRT * V_m);
            const cellmodel_float_t dI_kp_junc_dkp_kp = Fjunc * GKp * (-ek + V_m);
            const cellmodel_float_t dI_ClCa_junc_dV_m =
                    Fjunc * GClCa / (FP_LITERAL(1.) + KdClCa / Ca_j);
            const cellmodel_float_t dI_kp_sl_dkp_kp = GKp * (-ek + V_m) * Fsl;
            const cellmodel_float_t dI_nak_sl_dfnak =
                    IbarNaK * Ko * Fsl
                    / ((1.
                        + (((KmNaip) * (KmNaip)) * ((KmNaip) * (KmNaip)))
                                  / (((Na_sl) * (Na_sl)) * ((Na_sl) * (Na_sl))))
                       * (KmKo + Ko));
            const cellmodel_float_t dfnak_dV_m =
                    (FP_LITERAL(0.01245) * FoRT * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                     + FP_LITERAL(0.0365) * FoRT * Exp(-FoRT * V_m) * sigma)
                    / ((FP_LITERAL(1.) + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                        + FP_LITERAL(0.0365) * Exp(-FoRT * V_m) * sigma)
                       * (FP_LITERAL(1.) + FP_LITERAL(0.1245) * Exp(FP_LITERAL(-0.1) * FoRT * V_m)
                          + FP_LITERAL(0.0365) * Exp(-FoRT * V_m) * sigma));
            const cellmodel_float_t dI_K1_dkiss =
                    FP_LITERAL(0.430331482911935) * GK1 * sqrt(Ko) * (-ek + V_m);
            const cellmodel_float_t dibarna_j_dV_m =
                    Frdy * GCaL * pNa
                            * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_j * Exp(FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FoRT * V_m))
                    - Frdy * GCaL * pNa * (FoRT * FoRT)
                              * (FP_LITERAL(-0.75) * Nao
                                 + FP_LITERAL(0.75) * Na_j * Exp(FoRT * V_m))
                              * V_m * Exp(FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FoRT * V_m)))
                    + FP_LITERAL(0.75) * Frdy * GCaL * pNa * (FoRT * FoRT) * Na_j * V_m
                              * Exp(FoRT * V_m) / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t dkp_kp_dV_m =
                    FP_LITERAL(298.741733340907) * Exp(FP_LITERAL(-0.167224080267559) * V_m)
                    / ((FP_LITERAL(1.)
                        + FP_LITERAL(1786.47556537862) * Exp(FP_LITERAL(-0.167224080267559) * V_m))
                       * (FP_LITERAL(1.)
                          + FP_LITERAL(1786.47556537862)
                                    * Exp(FP_LITERAL(-0.167224080267559) * V_m)));
            const cellmodel_float_t dI_kp_sl_dV_m =
                    GKp * Fsl * kp_kp + GKp * (-ek + V_m) * Fsl * dkp_kp_dV_m;
            const cellmodel_float_t dI_kp_junc_dV_m =
                    Fjunc * GKp * kp_kp + Fjunc * GKp * (-ek + V_m) * dkp_kp_dV_m;
            const cellmodel_float_t dI_Na_junc_dV_m = Fjunc * GNa * (m * m * m) * h * j;
            const cellmodel_float_t drkr_dV_m =
                    -Exp(FP_LITERAL(37.) / FP_LITERAL(12.) + V_m / FP_LITERAL(24.))
                    / (FP_LITERAL(24.)
                       * ((FP_LITERAL(1.)
                           + Exp(FP_LITERAL(37.) / FP_LITERAL(12.) + V_m / FP_LITERAL(24.)))
                          * (FP_LITERAL(1.)
                             + Exp(FP_LITERAL(37.) / FP_LITERAL(12.) + V_m / FP_LITERAL(24.)))));
            const cellmodel_float_t ds1_junc_dV_m =
                    Cao * nu * (Na_j * Na_j * Na_j) * FoRT * Exp(nu * FoRT * V_m);
            const cellmodel_float_t dI_Na_sl_dV_m = GNa * (m * m * m) * Fsl * h * j;
            const cellmodel_float_t dI_ncx_junc_dV_m =
                    Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-ds2_junc_dV_m + ds1_junc_dV_m) * Ka_junc
                            / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                               * s3_junc)
                    - Fjunc * IbarNCX * ksat * pow(Q10NCX, Qpow) * (FP_LITERAL(-1.) + nu)
                              * (-s2_junc + s1_junc) * FoRT * Ka_junc
                              * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)
                              / (((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                                  * (FP_LITERAL(1.)
                                     + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)))
                                 * s3_junc);
            const cellmodel_float_t dI_kr_dV_m =
                    gkr * rkr * x_kr + (-ek + V_m) * drkr_dV_m * gkr * x_kr;
            const cellmodel_float_t dI_tof_dV_m = GtoFast * x_to_f * y_to_f;
            const cellmodel_float_t dI_ClCa_sl_dV_m =
                    GClCa * Fsl / (FP_LITERAL(1.) + KdClCa / Ca_sl);
            const cellmodel_float_t dI_ncx_sl_dV_m =
                    IbarNCX * pow(Q10NCX, Qpow) * (-ds2_sl_dV_m + ds1_sl_dV_m) * Fsl * Ka_sl
                            / ((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                               * s3_sl)
                    - IbarNCX * ksat * pow(Q10NCX, Qpow) * (FP_LITERAL(-1.) + nu) * (-s2_sl + s1_sl)
                              * FoRT * Fsl * Ka_sl * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)
                              / (((FP_LITERAL(1.) + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m))
                                  * (FP_LITERAL(1.)
                                     + ksat * Exp((FP_LITERAL(-1.) + nu) * FoRT * V_m)))
                                 * s3_sl);
            const cellmodel_float_t dI_tos_dV_m = GtoSlow * x_to_s * y_to_s;
            const cellmodel_float_t dibarna_sl_dV_m =
                    Frdy * GCaL * pNa
                            * (FP_LITERAL(-0.75) * Nao + FP_LITERAL(0.75) * Na_sl * Exp(FoRT * V_m))
                            * FoRT / (FP_LITERAL(-1.) + Exp(FoRT * V_m))
                    - Frdy * GCaL * pNa * (FoRT * FoRT)
                              * (FP_LITERAL(-0.75) * Nao
                                 + FP_LITERAL(0.75) * Na_sl * Exp(FoRT * V_m))
                              * V_m * Exp(FoRT * V_m)
                              / ((FP_LITERAL(-1.) + Exp(FoRT * V_m))
                                 * (FP_LITERAL(-1.) + Exp(FoRT * V_m)))
                    + FP_LITERAL(0.75) * Frdy * GCaL * pNa * (FoRT * FoRT) * Na_sl * V_m
                              * Exp(FoRT * V_m) / (FP_LITERAL(-1.) + Exp(FoRT * V_m));
            const cellmodel_float_t dV_m_dt_linearized =
                    -GClB - dI_ClCa_junc_dV_m - dI_ClCa_sl_dV_m - dI_K1_dV_m - dI_Na_junc_dV_m
                    - dI_Na_sl_dV_m - dI_kp_junc_dV_m - dI_kp_sl_dV_m - dI_kr_dV_m - dI_ks_junc_dV_m
                    - dI_ks_sl_dV_m - dI_ncx_junc_dV_m - dI_ncx_sl_dV_m - dI_tof_dV_m - dI_tos_dV_m
                    - Fjunc * GCaB - Fjunc * GNaB - GCaB * Fsl - GNaB * Fsl
                    - (daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki) * dI_K1_dkiss
                    - dI_CaK_dibark * dibark_dV_m - dI_CaNa_junc_dibarna_j * dibarna_j_dV_m
                    - dI_CaNa_sl_dibarna_sl * dibarna_sl_dV_m
                    - dI_Ca_junc_dibarca_j * dibarca_j_dV_m - dI_Ca_sl_dibarca_sl * dibarca_sl_dV_m
                    - dI_kp_junc_dkp_kp * dkp_kp_dV_m - dI_kp_sl_dkp_kp * dkp_kp_dV_m
                    - dI_kr_drkr * drkr_dV_m - dI_nak_junc_dfnak * dfnak_dV_m
                    - dI_nak_sl_dfnak * dfnak_dV_m - dI_ncx_junc_ds1_junc * ds1_junc_dV_m
                    - dI_ncx_junc_ds2_junc * ds2_junc_dV_m - dI_ncx_sl_ds1_sl * ds1_sl_dV_m
                    - dI_ncx_sl_ds2_sl * ds2_sl_dV_m;
            states[STATE_V_m * padded_num_cells + i] =
                    (fabs(dV_m_dt_linearized) > FP_LITERAL(1.0e-8)
                             ? (FP_LITERAL(-1.0) + Exp(dt * dV_m_dt_linearized)) * dV_m_dt
                                       / dV_m_dt_linearized
                             : dt * dV_m_dt)
                    + V_m;
        }
    }
}

const char *state_names[] = {
    "m",
    "h",
    "j",
    "x_kr",
    "x_ks",
    "x_to_s",
    "y_to_s",
    "x_to_f",
    "y_to_f",
    "d",
    "f",
    "f_Ca_Bj",
    "f_Ca_Bsl",
    "Ry_Rr",
    "Ry_Ro",
    "Ry_Ri",
    "Na_Bj",
    "Na_Bsl",
    "Tn_CL",
    "Tn_CHc",
    "Tn_CHm",
    "CaM",
    "Myo_c",
    "Myo_m",
    "SRB",
    "SLL_j",
    "SLL_sl",
    "SLH_j",
    "SLH_sl",
    "Csqn_b",
    "Ca_sr",
    "Na_j",
    "Na_sl",
    "Na_i",
    "K_i",
    "Ca_j",
    "Ca_sl",
    "Ca_i",
    "V_m",
    NULL
};

// clang-format off
const struct cellmodel model_GPB_simd = {
        .init_states = &init_state_values,
        .init_parameters = &init_parameters_values,
        .state_index = &state_index,
        .parameter_index = &parameter_index,
        .step_FE = &step_FE,
        .step_RL = NULL,
        .step_GRL1 = &step_GRL1,
        .num_states = NUM_STATES,
        .num_parameters = NUM_PARAMS,
        .layout = LAYOUT_STRUCT_OF_ARRAYS,
        .colour_sets = NULL,
        .state_names = state_names,
};
// clang-format on
