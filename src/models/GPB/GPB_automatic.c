#include <math.h>
#include <string.h>
// Gotran generated C/C++ code for the "GPB" model

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

enum monitored {
    NUM_MONITORED,
};

// Init state values
void init_state_values(double *states)
{
    states[STATE_m] = 0.003793087414436;
    states[STATE_h] = 0.626221949492493;
    states[STATE_j] = 0.624553572490432;
    states[STATE_x_kr] = 0.0210022533039071;
    states[STATE_x_ks] = 0.00428016666258923;
    states[STATE_x_to_s] = 0.000440445885642567;
    states[STATE_y_to_s] = 0.785115828275182;
    states[STATE_x_to_f] = 0.000440438103758954;
    states[STATE_y_to_f] = 0.999995844038706;
    states[STATE_d] = 2.92407183949469e-06;
    states[STATE_f] = 0.995135796703515;
    states[STATE_f_Ca_Bj] = 0.0246760872105795;
    states[STATE_f_Ca_Bsl] = 0.0152723084239416;
    states[STATE_Ry_Rr] = 0.890806040818203;
    states[STATE_Ry_Ro] = 7.40481128853622e-07;
    states[STATE_Ry_Ri] = 9.07666168960848e-08;
    states[STATE_Na_Bj] = 3.4543773303328;
    states[STATE_Na_Bsl] = 0.753740951477775;
    states[STATE_Tn_CL] = 0.00893455096919132;
    states[STATE_Tn_CHc] = 0.117412025936615;
    states[STATE_Tn_CHm] = 0.0106160166692932;
    states[STATE_CaM] = 0.000295573424135051;
    states[STATE_Myo_c] = 0.00192322252438022;
    states[STATE_Myo_m] = 0.137560495022823;
    states[STATE_SRB] = 0.00217360235649355;
    states[STATE_SLL_j] = 0.00740524521680039;
    states[STATE_SLL_sl] = 0.00990339304377132;
    states[STATE_SLH_j] = 0.0735890020284214;
    states[STATE_SLH_sl] = 0.114583623436917;
    states[STATE_Csqn_b] = 1.19723145924432;
    states[STATE_Ca_sr] = 0.554760499828172;
    states[STATE_Na_j] = 8.40537012592918;
    states[STATE_Na_sl] = 8.40491910001025;
    states[STATE_Na_i] = 8.40513364344858;
    states[STATE_K_i] = 120.0;
    states[STATE_Ca_j] = 0.000175882395147342;
    states[STATE_Ca_sl] = 0.000106779509977354;
    states[STATE_Ca_i] = 8.72509677797499e-05;
    states[STATE_V_m] = -81.4552030512661;
}

// Default parameter values
void init_parameters_values(double *parameters)
{
    parameters[PARAM_Fjunc] = 0.11;
    parameters[PARAM_Fjunc_CaL] = 0.9;
    parameters[PARAM_cellLength] = 100.0;
    parameters[PARAM_cellRadius] = 10.25;
    parameters[PARAM_distJuncSL] = 0.5;
    parameters[PARAM_distSLcyto] = 0.45;
    parameters[PARAM_junctionLength] = 0.16;
    parameters[PARAM_junctionRadius] = 0.015;
    parameters[PARAM_GNa] = 23.0;
    parameters[PARAM_GNaB] = 0.000597;
    parameters[PARAM_IbarNaK] = 1.8;
    parameters[PARAM_KmKo] = 1.5;
    parameters[PARAM_KmNaip] = 11.0;
    parameters[PARAM_Q10KmNai] = 1.39;
    parameters[PARAM_Q10NaK] = 1.63;
    parameters[PARAM_GKr] = 0.035;
    parameters[PARAM_GKp] = 0.002;
    parameters[PARAM_GKs] = 0.0035;
    parameters[PARAM_pNaK] = 0.01833;
    parameters[PARAM_GK1] = 0.35;
    parameters[PARAM_Gto] = 0.13;
    parameters[PARAM_epi] = 1.0;
    parameters[PARAM_GClB] = 0.009;
    parameters[PARAM_GClCa] = 0.0548125;
    parameters[PARAM_KdClCa] = 0.1;
    parameters[PARAM_GCaL] = 0.5;
    parameters[PARAM_Q10CaL] = 1.8;
    parameters[PARAM_pCa] = 0.00054;
    parameters[PARAM_pK] = 2.7e-07;
    parameters[PARAM_pNa] = 1.5e-08;
    parameters[PARAM_IbarNCX] = 4.5;
    parameters[PARAM_Kdact] = 0.00015;
    parameters[PARAM_KmCai] = 0.00359;
    parameters[PARAM_KmCao] = 1.3;
    parameters[PARAM_KmNai] = 12.29;
    parameters[PARAM_KmNao] = 87.5;
    parameters[PARAM_Q10NCX] = 1.57;
    parameters[PARAM_ksat] = 0.32;
    parameters[PARAM_nu] = 0.27;
    parameters[PARAM_IbarSLCaP] = 0.0673;
    parameters[PARAM_KmPCa] = 0.0005;
    parameters[PARAM_Q10SLCaP] = 2.35;
    parameters[PARAM_GCaB] = 0.0005513;
    parameters[PARAM_Kmf] = 0.000246;
    parameters[PARAM_Kmr] = 1.7;
    parameters[PARAM_MaxSR] = 15.0;
    parameters[PARAM_MinSR] = 1.0;
    parameters[PARAM_Q10SRCaP] = 2.6;
    parameters[PARAM_Vmax_SRCaP] = 0.0053114;
    parameters[PARAM_ec50SR] = 0.45;
    parameters[PARAM_hillSRCaP] = 1.787;
    parameters[PARAM_kiCa] = 0.5;
    parameters[PARAM_kim] = 0.005;
    parameters[PARAM_koCa] = 10.0;
    parameters[PARAM_kom] = 0.06;
    parameters[PARAM_ks] = 25.0;
    parameters[PARAM_Bmax_Naj] = 7.561;
    parameters[PARAM_Bmax_Nasl] = 1.65;
    parameters[PARAM_koff_na] = 0.001;
    parameters[PARAM_kon_na] = 0.0001;
    parameters[PARAM_Bmax_CaM] = 0.024;
    parameters[PARAM_Bmax_SR] = 0.0171;
    parameters[PARAM_Bmax_TnChigh] = 0.14;
    parameters[PARAM_Bmax_TnClow] = 0.07;
    parameters[PARAM_Bmax_myosin] = 0.14;
    parameters[PARAM_koff_cam] = 0.238;
    parameters[PARAM_koff_myoca] = 0.00046;
    parameters[PARAM_koff_myomg] = 5.7e-05;
    parameters[PARAM_koff_sr] = 0.06;
    parameters[PARAM_koff_tnchca] = 3.2e-05;
    parameters[PARAM_koff_tnchmg] = 0.00333;
    parameters[PARAM_koff_tncl] = 0.0196;
    parameters[PARAM_kon_cam] = 34.0;
    parameters[PARAM_kon_myoca] = 13.8;
    parameters[PARAM_kon_myomg] = 0.0157;
    parameters[PARAM_kon_sr] = 100.0;
    parameters[PARAM_kon_tnchca] = 2.37;
    parameters[PARAM_kon_tnchmg] = 0.003;
    parameters[PARAM_kon_tncl] = 32.7;
    parameters[PARAM_Bmax_SLhighj0] = 0.000165;
    parameters[PARAM_Bmax_SLhighsl0] = 0.0134;
    parameters[PARAM_Bmax_SLlowj0] = 0.00046;
    parameters[PARAM_Bmax_SLlowsl0] = 0.0374;
    parameters[PARAM_koff_slh] = 0.03;
    parameters[PARAM_koff_sll] = 1.3;
    parameters[PARAM_kon_slh] = 100.0;
    parameters[PARAM_kon_sll] = 100.0;
    parameters[PARAM_Bmax_Csqn0] = 0.14;
    parameters[PARAM_DcaJuncSL] = 1.64e-06;
    parameters[PARAM_DcaSLcyto] = 1.22e-06;
    parameters[PARAM_J_ca_juncsl] = 8.2413e-13;
    parameters[PARAM_J_ca_slmyo] = 3.7243e-12;
    parameters[PARAM_koff_csqn] = 65.0;
    parameters[PARAM_kon_csqn] = 100.0;
    parameters[PARAM_DnaJuncSL] = 1.09e-05;
    parameters[PARAM_DnaSLcyto] = 1.79e-05;
    parameters[PARAM_J_na_juncsl] = 1.8313e-14;
    parameters[PARAM_J_na_slmyo] = 1.6386e-12;
    parameters[PARAM_Nao] = 140.0;
    parameters[PARAM_Ko] = 5.4;
    parameters[PARAM_Cao] = 1.8;
    parameters[PARAM_Cli] = 15.0;
    parameters[PARAM_Clo] = 150.0;
    parameters[PARAM_Mgi] = 1.0;
    parameters[PARAM_Cmem] = 1.381e-10;
    parameters[PARAM_Frdy] = 96485.0;
    parameters[PARAM_R] = 8314.0;
    parameters[PARAM_Temp] = 310.0;
    parameters[PARAM_stim_amplitude] = 40.0;
    parameters[PARAM_stim_duration] = 1.0;
    parameters[PARAM_stim_period] = 1000.0;
    parameters[PARAM_stim_start] = 0.0;
}

// State index
int state_index(const char name[])
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
int parameter_index(const char name[])
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
void step_FE_singlecell(double *__restrict states, const double t, const double dt,
                        const double *__restrict parameters)
{

    // Assign states
    const double m = states[STATE_m];
    const double h = states[STATE_h];
    const double j = states[STATE_j];
    const double x_kr = states[STATE_x_kr];
    const double x_ks = states[STATE_x_ks];
    const double x_to_s = states[STATE_x_to_s];
    const double y_to_s = states[STATE_y_to_s];
    const double x_to_f = states[STATE_x_to_f];
    const double y_to_f = states[STATE_y_to_f];
    const double d = states[STATE_d];
    const double f = states[STATE_f];
    const double f_Ca_Bj = states[STATE_f_Ca_Bj];
    const double f_Ca_Bsl = states[STATE_f_Ca_Bsl];
    const double Ry_Rr = states[STATE_Ry_Rr];
    const double Ry_Ro = states[STATE_Ry_Ro];
    const double Ry_Ri = states[STATE_Ry_Ri];
    const double Na_Bj = states[STATE_Na_Bj];
    const double Na_Bsl = states[STATE_Na_Bsl];
    const double Tn_CL = states[STATE_Tn_CL];
    const double Tn_CHc = states[STATE_Tn_CHc];
    const double Tn_CHm = states[STATE_Tn_CHm];
    const double CaM = states[STATE_CaM];
    const double Myo_c = states[STATE_Myo_c];
    const double Myo_m = states[STATE_Myo_m];
    const double SRB = states[STATE_SRB];
    const double SLL_j = states[STATE_SLL_j];
    const double SLL_sl = states[STATE_SLL_sl];
    const double SLH_j = states[STATE_SLH_j];
    const double SLH_sl = states[STATE_SLH_sl];
    const double Csqn_b = states[STATE_Csqn_b];
    const double Ca_sr = states[STATE_Ca_sr];
    const double Na_j = states[STATE_Na_j];
    const double Na_sl = states[STATE_Na_sl];
    const double Na_i = states[STATE_Na_i];
    const double K_i = states[STATE_K_i];
    const double Ca_j = states[STATE_Ca_j];
    const double Ca_sl = states[STATE_Ca_sl];
    const double Ca_i = states[STATE_Ca_i];
    const double V_m = states[STATE_V_m];

    // Assign parameters
    const double Fjunc = parameters[PARAM_Fjunc];
    const double Fjunc_CaL = parameters[PARAM_Fjunc_CaL];
    const double cellLength = parameters[PARAM_cellLength];
    const double cellRadius = parameters[PARAM_cellRadius];
    const double GNa = parameters[PARAM_GNa];
    const double GNaB = parameters[PARAM_GNaB];
    const double IbarNaK = parameters[PARAM_IbarNaK];
    const double KmKo = parameters[PARAM_KmKo];
    const double KmNaip = parameters[PARAM_KmNaip];
    const double GKr = parameters[PARAM_GKr];
    const double GKp = parameters[PARAM_GKp];
    const double GKs = parameters[PARAM_GKs];
    const double pNaK = parameters[PARAM_pNaK];
    const double GK1 = parameters[PARAM_GK1];
    const double Gto = parameters[PARAM_Gto];
    const double epi = parameters[PARAM_epi];
    const double GClB = parameters[PARAM_GClB];
    const double GClCa = parameters[PARAM_GClCa];
    const double KdClCa = parameters[PARAM_KdClCa];
    const double GCaL = parameters[PARAM_GCaL];
    const double Q10CaL = parameters[PARAM_Q10CaL];
    const double pCa = parameters[PARAM_pCa];
    const double pK = parameters[PARAM_pK];
    const double pNa = parameters[PARAM_pNa];
    const double IbarNCX = parameters[PARAM_IbarNCX];
    const double Kdact = parameters[PARAM_Kdact];
    const double KmCai = parameters[PARAM_KmCai];
    const double KmCao = parameters[PARAM_KmCao];
    const double KmNai = parameters[PARAM_KmNai];
    const double KmNao = parameters[PARAM_KmNao];
    const double Q10NCX = parameters[PARAM_Q10NCX];
    const double ksat = parameters[PARAM_ksat];
    const double nu = parameters[PARAM_nu];
    const double IbarSLCaP = parameters[PARAM_IbarSLCaP];
    const double KmPCa = parameters[PARAM_KmPCa];
    const double Q10SLCaP = parameters[PARAM_Q10SLCaP];
    const double GCaB = parameters[PARAM_GCaB];
    const double Kmf = parameters[PARAM_Kmf];
    const double Kmr = parameters[PARAM_Kmr];
    const double MaxSR = parameters[PARAM_MaxSR];
    const double MinSR = parameters[PARAM_MinSR];
    const double Q10SRCaP = parameters[PARAM_Q10SRCaP];
    const double Vmax_SRCaP = parameters[PARAM_Vmax_SRCaP];
    const double ec50SR = parameters[PARAM_ec50SR];
    const double hillSRCaP = parameters[PARAM_hillSRCaP];
    const double kiCa = parameters[PARAM_kiCa];
    const double kim = parameters[PARAM_kim];
    const double koCa = parameters[PARAM_koCa];
    const double kom = parameters[PARAM_kom];
    const double ks = parameters[PARAM_ks];
    const double Bmax_Naj = parameters[PARAM_Bmax_Naj];
    const double Bmax_Nasl = parameters[PARAM_Bmax_Nasl];
    const double koff_na = parameters[PARAM_koff_na];
    const double kon_na = parameters[PARAM_kon_na];
    const double Bmax_CaM = parameters[PARAM_Bmax_CaM];
    const double Bmax_SR = parameters[PARAM_Bmax_SR];
    const double Bmax_TnChigh = parameters[PARAM_Bmax_TnChigh];
    const double Bmax_TnClow = parameters[PARAM_Bmax_TnClow];
    const double Bmax_myosin = parameters[PARAM_Bmax_myosin];
    const double koff_cam = parameters[PARAM_koff_cam];
    const double koff_myoca = parameters[PARAM_koff_myoca];
    const double koff_myomg = parameters[PARAM_koff_myomg];
    const double koff_sr = parameters[PARAM_koff_sr];
    const double koff_tnchca = parameters[PARAM_koff_tnchca];
    const double koff_tnchmg = parameters[PARAM_koff_tnchmg];
    const double koff_tncl = parameters[PARAM_koff_tncl];
    const double kon_cam = parameters[PARAM_kon_cam];
    const double kon_myoca = parameters[PARAM_kon_myoca];
    const double kon_myomg = parameters[PARAM_kon_myomg];
    const double kon_sr = parameters[PARAM_kon_sr];
    const double kon_tnchca = parameters[PARAM_kon_tnchca];
    const double kon_tnchmg = parameters[PARAM_kon_tnchmg];
    const double kon_tncl = parameters[PARAM_kon_tncl];
    const double Bmax_SLhighj0 = parameters[PARAM_Bmax_SLhighj0];
    const double Bmax_SLhighsl0 = parameters[PARAM_Bmax_SLhighsl0];
    const double Bmax_SLlowj0 = parameters[PARAM_Bmax_SLlowj0];
    const double Bmax_SLlowsl0 = parameters[PARAM_Bmax_SLlowsl0];
    const double koff_slh = parameters[PARAM_koff_slh];
    const double koff_sll = parameters[PARAM_koff_sll];
    const double kon_slh = parameters[PARAM_kon_slh];
    const double kon_sll = parameters[PARAM_kon_sll];
    const double Bmax_Csqn0 = parameters[PARAM_Bmax_Csqn0];
    const double J_ca_juncsl = parameters[PARAM_J_ca_juncsl];
    const double J_ca_slmyo = parameters[PARAM_J_ca_slmyo];
    const double koff_csqn = parameters[PARAM_koff_csqn];
    const double kon_csqn = parameters[PARAM_kon_csqn];
    const double J_na_juncsl = parameters[PARAM_J_na_juncsl];
    const double J_na_slmyo = parameters[PARAM_J_na_slmyo];
    const double Nao = parameters[PARAM_Nao];
    const double Ko = parameters[PARAM_Ko];
    const double Cao = parameters[PARAM_Cao];
    const double Cli = parameters[PARAM_Cli];
    const double Clo = parameters[PARAM_Clo];
    const double Mgi = parameters[PARAM_Mgi];
    const double Cmem = parameters[PARAM_Cmem];
    const double Frdy = parameters[PARAM_Frdy];
    const double R = parameters[PARAM_R];
    const double Temp = parameters[PARAM_Temp];
    const double stim_amplitude = parameters[PARAM_stim_amplitude];
    const double stim_duration = parameters[PARAM_stim_duration];
    const double stim_period = parameters[PARAM_stim_period];
    const double stim_start = parameters[PARAM_stim_start];

    // Expressions for the Geometry component
    const double Vcell = 1.0e-15 * M_PI * cellLength * (cellRadius * cellRadius);
    const double Vmyo = 0.65 * Vcell;
    const double Vsr = 0.035 * Vcell;
    const double Vsl = 0.02 * Vcell;
    const double Vjunc = 0.000539 * Vcell;
    const double Fsl = 1. - Fjunc;
    const double Fsl_CaL = 1. - Fjunc_CaL;

    // Expressions for the Reversal potentials component
    const double FoRT = Frdy / (R * Temp);
    const double ena_junc = log(Nao / Na_j) / FoRT;
    const double ena_sl = log(Nao / Na_sl) / FoRT;
    const double ek = log(Ko / K_i) / FoRT;
    const double eca_junc = log(Cao / Ca_j) / (2. * FoRT);
    const double eca_sl = log(Cao / Ca_sl) / (2. * FoRT);
    const double ecl = log(Cli / Clo) / FoRT;
    const double Qpow = -31. + Temp / 10.;

    // Expressions for the I_Na component
    const double mss = 1.0
                       / ((1. + 0.00184221158116513 * exp(-0.110741971207087 * V_m))
                          * (1. + 0.00184221158116513 * exp(-0.110741971207087 * V_m)));
    const double taum = 0.1292
                                * exp(-((2.94658944658945 + 0.0643500643500644 * V_m)
                                        * (2.94658944658945 + 0.0643500643500644 * V_m)))
                        + 0.06487
                                  * exp(-((-0.0943466353677621 + 0.0195618153364632 * V_m)
                                          * (-0.0943466353677621 + 0.0195618153364632 * V_m)));
    const double ah = (V_m >= -40. ? 0. : 4.43126792958051e-7 * exp(-0.147058823529412 * V_m));
    const double bh =
            (V_m >= -40. ? 0.77 / (0.13 + 0.0497581410839387 * exp(-0.0900900900900901 * V_m))
                         : 310000.0 * exp(0.3485 * V_m) + 2.7 * exp(0.079 * V_m));
    const double tauh = 1.0 / (ah + bh);
    const double hss = 1.0
                       / ((1. + 15212.5932856544 * exp(0.134589502018843 * V_m))
                          * (1. + 15212.5932856544 * exp(0.134589502018843 * V_m)));
    const double aj =
            (V_m >= -40. ? 0.
                         : (37.78 + V_m)
                                   * (-25428.0 * exp(0.2444 * V_m) - 6.948e-6 * exp(-0.04391 * V_m))
                                   / (1. + 50262745825.954 * exp(0.311 * V_m)));
    const double bj =
            (V_m >= -40. ? 0.6 * exp(0.057 * V_m) / (1. + 0.0407622039783662 * exp(-0.1 * V_m))
                         : 0.02424 * exp(-0.01052 * V_m)
                                   / (1. + 0.00396086833990426 * exp(-0.1378 * V_m)));
    const double tauj = 1.0 / (aj + bj);
    const double jss = 1.0
                       / ((1. + 15212.5932856544 * exp(0.134589502018843 * V_m))
                          * (1. + 15212.5932856544 * exp(0.134589502018843 * V_m)));
    const double dm_dt = (-m + mss) / taum;
    states[STATE_m] = dt * dm_dt + m;
    const double dh_dt = (-h + hss) / tauh;
    states[STATE_h] = dt * dh_dt + h;
    const double dj_dt = (-j + jss) / tauj;
    states[STATE_j] = dt * dj_dt + j;
    const double I_Na_junc = Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j;
    const double I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j;

    // Expressions for the I_NaBK component
    const double I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m);
    const double I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl;

    // Expressions for the I_NaK component
    const double sigma = -1. / 7. + exp(0.0148588410104012 * Nao) / 7.;
    const double fnak =
            1.0 / (1. + 0.1245 * exp(-0.1 * FoRT * V_m) + 0.0365 * exp(-FoRT * V_m) * sigma);
    const double I_nak_junc =
            Fjunc * IbarNaK * Ko * fnak / ((1. + pow(KmNaip, 4.) / pow(Na_j, 4.)) * (KmKo + Ko));
    const double I_nak_sl =
            IbarNaK * Ko * Fsl * fnak / ((1. + pow(KmNaip, 4.) / pow(Na_sl, 4.)) * (KmKo + Ko));
    const double I_nak = I_nak_junc + I_nak_sl;

    // Expressions for the I_Kr component
    const double gkr = 0.430331482911935 * GKr * sqrt(Ko);
    const double xrss = 1.0 / (1. + exp(-2. - V_m / 5.));
    const double tauxr =
            230. / (1. + exp(2. + V_m / 20.))
            + 3300. / ((1. + exp(-22. / 9. - V_m / 9.)) * (1. + exp(11. / 9. + V_m / 9.)));
    const double dx_kr_dt = (-x_kr + xrss) / tauxr;
    states[STATE_x_kr] = dt * dx_kr_dt + x_kr;
    const double rkr = 1.0 / (1. + exp(37. / 12. + V_m / 24.));
    const double I_kr = (-ek + V_m) * gkr * rkr * x_kr;

    // Expressions for the I_Kp component
    const double kp_kp = 1.0 / (1. + 1786.47556537862 * exp(-0.167224080267559 * V_m));
    const double I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp;
    const double I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp;
    const double I_kp = I_kp_junc + I_kp_sl;

    // Expressions for the I_Ks component
    const double eks = log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT;
    const double gks_junc = GKs;
    const double gks_sl = GKs;
    const double xsss = 1.0 / (1. + 0.765928338364649 * exp(-0.0701754385964912 * V_m));
    const double tauxs = 990.1 / (1. + 0.841540408868102 * exp(-0.0708215297450425 * V_m));
    const double dx_ks_dt = (-x_ks + xsss) / tauxs;
    states[STATE_x_ks] = dt * dx_ks_dt + x_ks;
    const double I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc;
    const double I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl;
    const double I_ks = I_ks_junc + I_ks_sl;

    // Expressions for the I_to component
    const double GtoSlow = (epi == 1. ? 0.12 * Gto : 0.2892 * Gto);
    const double GtoFast = (epi == 1. ? 0.88 * Gto : 0.0108 * Gto);
    const double xtoss = 1.0 / (1. + exp(19. / 13. - V_m / 13.));
    const double ytoss = 1.0 / (1. + 49.4024491055302 * exp(V_m / 5.));
    const double tauxtos = 0.5 + 9. / (1. + exp(1. / 5. + V_m / 15.));
    const double tauytos = 30. + 800. / (1. + exp(6. + V_m / 10.));
    const double dx_to_s_dt = (-x_to_s + xtoss) / tauxtos;
    states[STATE_x_to_s] = dt * dx_to_s_dt + x_to_s;
    const double dy_to_s_dt = (-y_to_s + ytoss) / tauytos;
    states[STATE_y_to_s] = dt * dy_to_s_dt + y_to_s;
    const double I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s;
    const double tauxtof = 0.5 + 8.5 * exp(-((9. / 10. + V_m / 50.) * (9. / 10. + V_m / 50.)));
    const double tauytof = 7. + 85. * exp(-((40. + V_m) * (40. + V_m)) / 220.);
    const double dx_to_f_dt = (-x_to_f + xtoss) / tauxtof;
    states[STATE_x_to_f] = dt * dx_to_f_dt + x_to_f;
    const double dy_to_f_dt = (-y_to_f + ytoss) / tauytof;
    states[STATE_y_to_f] = dt * dy_to_f_dt + y_to_f;
    const double I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f;
    const double I_to = I_tof + I_tos;

    // Expressions for the I_K1 component
    const double aki = 1.02 / (1. + 7.35454251046446e-7 * exp(0.2385 * V_m - 0.2385 * ek));
    const double bki = (0.762624006506308 * exp(0.08032 * V_m - 0.08032 * ek)
                        + 1.15340563518656e-16 * exp(0.06175 * V_m - 0.06175 * ek))
                       / (1. + 0.0867722941576933 * exp(0.5143 * ek - 0.5143 * V_m));
    const double kiss = aki / (aki + bki);
    const double I_K1 = 0.430331482911935 * GK1 * sqrt(Ko) * (-ek + V_m) * kiss;

    // Expressions for the I_ClCa component
    const double I_ClCa_junc = Fjunc * GClCa * (-ecl + V_m) / (1. + KdClCa / Ca_j);
    const double I_ClCa_sl = GClCa * (-ecl + V_m) * Fsl / (1. + KdClCa / Ca_sl);
    const double I_ClCa = I_ClCa_junc + I_ClCa_sl;
    const double I_Clbk = GClB * (-ecl + V_m);

    // Expressions for the I_Ca component
    const double fss =
            1.0 / (1. + exp(35. / 9. + V_m / 9.)) + 0.6 / (1. + exp(5. / 2. - V_m / 20.));
    const double dss = 1.0 / (1. + exp(-5. / 6. - V_m / 6.));
    const double taud = (1. - exp(-5. / 6. - V_m / 6.)) * dss / (0.175 + 0.035 * V_m);
    const double tauf =
            1.0 / (0.02 + 0.0197 * exp(-((0.48865 + 0.0337 * V_m) * (0.48865 + 0.0337 * V_m))));
    const double dd_dt = (-d + dss) / taud;
    states[STATE_d] = dt * dd_dt + d;
    const double df_dt = (-f + fss) / tauf;
    states[STATE_f] = dt * df_dt + f;
    const double df_Ca_Bj_dt = -0.0119 * f_Ca_Bj + 1.7 * (1. - f_Ca_Bj) * Ca_j;
    states[STATE_f_Ca_Bj] = dt * df_Ca_Bj_dt + f_Ca_Bj;
    const double df_Ca_Bsl_dt = -0.0119 * f_Ca_Bsl + 1.7 * (1. - f_Ca_Bsl) * Ca_sl;
    states[STATE_f_Ca_Bsl] = dt * df_Ca_Bsl_dt + f_Ca_Bsl;
    const double fcaCaMSL = 0.;
    const double fcaCaj = 0.;
    const double ibarca_j = 4. * Frdy * GCaL * pCa
                            * (-0.341 * Cao + 0.341 * Ca_j * exp(2. * FoRT * V_m)) * FoRT * V_m
                            / (-1. + exp(2. * FoRT * V_m));
    const double ibarca_sl = 4. * Frdy * GCaL * pCa
                             * (-0.341 * Cao + 0.341 * Ca_sl * exp(2. * FoRT * V_m)) * FoRT * V_m
                             / (-1. + exp(2. * FoRT * V_m));
    const double ibark = Frdy * GCaL * pK * (-0.75 * Ko + 0.75 * K_i * exp(FoRT * V_m)) * FoRT * V_m
                         / (-1. + exp(FoRT * V_m));
    const double ibarna_j = Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_j * exp(FoRT * V_m)) * FoRT
                            * V_m / (-1. + exp(FoRT * V_m));
    const double ibarna_sl = Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_sl * exp(FoRT * V_m))
                             * FoRT * V_m / (-1. + exp(FoRT * V_m));
    const double I_Ca_junc =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f * ibarca_j;
    const double I_Ca_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f * ibarca_sl;
    const double I_CaK =
            0.45 * pow(Q10CaL, Qpow)
            * (Fjunc_CaL * (1. + fcaCaj - f_Ca_Bj) + (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL) * d * f
            * ibark;
    const double I_CaNa_junc =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f * ibarna_j;
    const double I_CaNa_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f * ibarna_sl;

    // Expressions for the I_NCX component
    const double Ka_junc = 1.0 / (1. + (Kdact * Kdact) / (Ca_j * Ca_j));
    const double Ka_sl = 1.0 / (1. + (Kdact * Kdact) / (Ca_sl * Ca_sl));
    const double s1_junc = Cao * (Na_j * Na_j * Na_j) * exp(nu * FoRT * V_m);
    const double s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * exp(nu * FoRT * V_m);
    const double s2_junc = (Nao * Nao * Nao) * Ca_j * exp((-1. + nu) * FoRT * V_m);
    const double s3_junc =
            Cao * (Na_j * Na_j * Na_j) + KmCao * (Na_j * Na_j * Na_j) + (Nao * Nao * Nao) * Ca_j
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + Ca_j / KmCai) * Ca_j;
    const double s2_sl = (Nao * Nao * Nao) * Ca_sl * exp((-1. + nu) * FoRT * V_m);
    const double s3_sl =
            Cao * (Na_sl * Na_sl * Na_sl) + KmCao * (Na_sl * Na_sl * Na_sl)
            + (Nao * Nao * Nao) * Ca_sl
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + Ca_sl / KmCai) * Ca_sl;
    const double I_ncx_junc = Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
                              / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc);
    const double I_ncx_sl = IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
                            / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);

    // Expressions for the I_PCa component
    const double I_pca_junc = Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, 1.6)
                              / (pow(KmPCa, 1.6) + pow(Ca_j, 1.6));
    const double I_pca_sl = IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, 1.6) * Fsl
                            / (pow(KmPCa, 1.6) + pow(Ca_sl, 1.6));

    // Expressions for the I_CaBK component
    const double I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m);
    const double I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl;

    // Expressions for the SR Fluxes component
    const double kCaSR = MaxSR - (MaxSR - MinSR) / (1. + pow(ec50SR / Ca_sr, 2.5));
    const double koSRCa = koCa / kCaSR;
    const double kiSRCa = kiCa * kCaSR;
    const double RI = 1. - Ry_Ri - Ry_Ro - Ry_Rr;
    const double dRy_Rr_dt =
            kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa;
    states[STATE_Ry_Rr] = dt * dRy_Rr_dt + Ry_Rr;
    const double dRy_Ro_dt =
            kim * Ry_Ri - kom * Ry_Ro + (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Ro * kiSRCa;
    states[STATE_Ry_Ro] = dt * dRy_Ro_dt + Ry_Ro;
    const double dRy_Ri_dt =
            -kim * Ry_Ri - kom * Ry_Ri + (Ca_j * Ca_j) * RI * koSRCa + Ca_j * Ry_Ro * kiSRCa;
    states[STATE_Ry_Ri] = dt * dRy_Ri_dt + Ry_Ri;
    const double J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro;
    const double J_serca = Vmax_SRCaP * pow(Q10SRCaP, Qpow)
                           * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                           / (1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP));
    const double J_SRleak = 5.348e-6 * Ca_sr - 5.348e-6 * Ca_j;

    // Expressions for the Na Buffers component
    const double dNa_Bj_dt = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j;
    states[STATE_Na_Bj] = dt * dNa_Bj_dt + Na_Bj;
    const double dNa_Bsl_dt = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl;
    states[STATE_Na_Bsl] = dt * dNa_Bsl_dt + Na_Bsl;

    // Expressions for the Cytosolic Ca Buffers component
    const double dTn_CL_dt = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;
    states[STATE_Tn_CL] = dt * dTn_CL_dt + Tn_CL;
    const double dTn_CHc_dt =
            -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i;
    states[STATE_Tn_CHc] = dt * dTn_CHc_dt + Tn_CHc;
    const double dTn_CHm_dt =
            -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm);
    states[STATE_Tn_CHm] = dt * dTn_CHm_dt + Tn_CHm;
    const double dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i;
    states[STATE_CaM] = dt * dCaM_dt + CaM;
    const double dMyo_c_dt = -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i;
    states[STATE_Myo_c] = dt * dMyo_c_dt + Myo_c;
    const double dMyo_m_dt = -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m);
    states[STATE_Myo_m] = dt * dMyo_m_dt + Myo_m;
    const double dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i;
    states[STATE_SRB] = dt * dSRB_dt + SRB;
    const double J_CaB_cytosol =
            -koff_cam * CaM - koff_myoca * Myo_c - koff_myomg * Myo_m - koff_sr * SRB
            - koff_tnchca * Tn_CHc - koff_tnchmg * Tn_CHm - koff_tncl * Tn_CL
            + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
            + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
            + kon_cam * (Bmax_CaM - CaM) * Ca_i + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
            + kon_sr * (Bmax_SR - SRB) * Ca_i + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
            + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;

    // Expressions for the Junctional and SL Ca Buffers component
    const double Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl;
    const double Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc;
    const double Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl;
    const double Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc;
    const double dSLL_j_dt = -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
    states[STATE_SLL_j] = dt * dSLL_j_dt + SLL_j;
    const double dSLL_sl_dt = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;
    states[STATE_SLL_sl] = dt * dSLL_sl_dt + SLL_sl;
    const double dSLH_j_dt = -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j;
    states[STATE_SLH_j] = dt * dSLH_j_dt + SLH_j;
    const double dSLH_sl_dt = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl;
    states[STATE_SLH_sl] = dt * dSLH_sl_dt + SLH_sl;
    const double J_CaB_junction = -koff_slh * SLH_j - koff_sll * SLL_j
                                  + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
                                  + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
    const double J_CaB_sl = -koff_slh * SLH_sl - koff_sll * SLL_sl
                            + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
                            + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;

    // Expressions for the SR Ca Concentrations component
    const double Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr;
    const double dCsqn_b_dt = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr;
    states[STATE_Csqn_b] = dt * dCsqn_b_dt + Csqn_b;
    const double dCa_sr_dt = -J_SRCarel + koff_csqn * Csqn_b
                             - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr - J_SRleak * Vmyo / Vsr
                             + J_serca;
    states[STATE_Ca_sr] = dt * dCa_sr_dt + Ca_sr;

    // Expressions for the Na Concentrations component
    const double I_Na_tot_junc =
            3. * I_nak_junc + 3. * I_ncx_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc;
    const double I_Na_tot_sl = 3. * I_nak_sl + 3. * I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl;
    const double dNa_j_dt = -dNa_Bj_dt + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
                            - Cmem * I_Na_tot_junc / (Frdy * Vjunc);
    states[STATE_Na_j] = dt * dNa_j_dt + Na_j;
    const double dNa_sl_dt = -dNa_Bsl_dt + J_na_juncsl * (-Na_sl + Na_j) / Vsl
                             + J_na_slmyo * (-Na_sl + Na_i) / Vsl
                             - Cmem * I_Na_tot_sl / (Frdy * Vsl);
    states[STATE_Na_sl] = dt * dNa_sl_dt + Na_sl;
    const double dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo;
    states[STATE_Na_i] = dt * dNa_i_dt + Na_i;

    // Expressions for the K Concentration component
    const double I_K_tot = -2. * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to;
    const double dK_i_dt = 0.;
    states[STATE_K_i] = dt * dK_i_dt + K_i;

    // Expressions for the Ca Concentrations component
    const double I_Ca_tot_junc = -2. * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
    const double I_Ca_tot_sl = -2. * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
    const double dCa_j_dt = -J_CaB_junction + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
                            + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc
                            - Cmem * I_Ca_tot_junc / (2. * Frdy * Vjunc);
    states[STATE_Ca_j] = dt * dCa_j_dt + Ca_j;
    const double dCa_sl_dt = -J_CaB_sl + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
                             + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
                             - Cmem * I_Ca_tot_sl / (2. * Frdy * Vsl);
    states[STATE_Ca_sl] = dt * dCa_sl_dt + Ca_sl;
    const double dCa_i_dt =
            -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo;
    states[STATE_Ca_i] = dt * dCa_i_dt + Ca_i;

    // Expressions for the Membrane potential component
    const double i_Stim = (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                           && t - stim_period * floor(t / stim_period) >= stim_start
                                   ? -stim_amplitude
                                   : 0.);
    const double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
    const double I_Cl_tot = I_ClCa + I_Clbk;
    const double I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
    const double I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot;
    const double dV_m_dt = -I_tot - i_Stim;
    states[STATE_V_m] = dt * dV_m_dt + V_m;
}

// Compute a forward step using the rush larsen algorithm to the GPB ODE
void step_GRL1_singlecell(double *__restrict states, const double t, const double dt,
                          const double *__restrict parameters)
{

    // Assign states
    const double m = states[STATE_m];
    const double h = states[STATE_h];
    const double j = states[STATE_j];
    const double x_kr = states[STATE_x_kr];
    const double x_ks = states[STATE_x_ks];
    const double x_to_s = states[STATE_x_to_s];
    const double y_to_s = states[STATE_y_to_s];
    const double x_to_f = states[STATE_x_to_f];
    const double y_to_f = states[STATE_y_to_f];
    const double d = states[STATE_d];
    const double f = states[STATE_f];
    const double f_Ca_Bj = states[STATE_f_Ca_Bj];
    const double f_Ca_Bsl = states[STATE_f_Ca_Bsl];
    const double Ry_Rr = states[STATE_Ry_Rr];
    const double Ry_Ro = states[STATE_Ry_Ro];
    const double Ry_Ri = states[STATE_Ry_Ri];
    const double Na_Bj = states[STATE_Na_Bj];
    const double Na_Bsl = states[STATE_Na_Bsl];
    const double Tn_CL = states[STATE_Tn_CL];
    const double Tn_CHc = states[STATE_Tn_CHc];
    const double Tn_CHm = states[STATE_Tn_CHm];
    const double CaM = states[STATE_CaM];
    const double Myo_c = states[STATE_Myo_c];
    const double Myo_m = states[STATE_Myo_m];
    const double SRB = states[STATE_SRB];
    const double SLL_j = states[STATE_SLL_j];
    const double SLL_sl = states[STATE_SLL_sl];
    const double SLH_j = states[STATE_SLH_j];
    const double SLH_sl = states[STATE_SLH_sl];
    const double Csqn_b = states[STATE_Csqn_b];
    const double Ca_sr = states[STATE_Ca_sr];
    const double Na_j = states[STATE_Na_j];
    const double Na_sl = states[STATE_Na_sl];
    const double Na_i = states[STATE_Na_i];
    const double K_i = states[STATE_K_i];
    const double Ca_j = states[STATE_Ca_j];
    const double Ca_sl = states[STATE_Ca_sl];
    const double Ca_i = states[STATE_Ca_i];
    const double V_m = states[STATE_V_m];

    // Assign parameters
    const double Fjunc = parameters[PARAM_Fjunc];
    const double Fjunc_CaL = parameters[PARAM_Fjunc_CaL];
    const double cellLength = parameters[PARAM_cellLength];
    const double cellRadius = parameters[PARAM_cellRadius];
    const double GNa = parameters[PARAM_GNa];
    const double GNaB = parameters[PARAM_GNaB];
    const double IbarNaK = parameters[PARAM_IbarNaK];
    const double KmKo = parameters[PARAM_KmKo];
    const double KmNaip = parameters[PARAM_KmNaip];
    const double GKr = parameters[PARAM_GKr];
    const double GKp = parameters[PARAM_GKp];
    const double GKs = parameters[PARAM_GKs];
    const double pNaK = parameters[PARAM_pNaK];
    const double GK1 = parameters[PARAM_GK1];
    const double Gto = parameters[PARAM_Gto];
    const double epi = parameters[PARAM_epi];
    const double GClB = parameters[PARAM_GClB];
    const double GClCa = parameters[PARAM_GClCa];
    const double KdClCa = parameters[PARAM_KdClCa];
    const double GCaL = parameters[PARAM_GCaL];
    const double Q10CaL = parameters[PARAM_Q10CaL];
    const double pCa = parameters[PARAM_pCa];
    const double pK = parameters[PARAM_pK];
    const double pNa = parameters[PARAM_pNa];
    const double IbarNCX = parameters[PARAM_IbarNCX];
    const double Kdact = parameters[PARAM_Kdact];
    const double KmCai = parameters[PARAM_KmCai];
    const double KmCao = parameters[PARAM_KmCao];
    const double KmNai = parameters[PARAM_KmNai];
    const double KmNao = parameters[PARAM_KmNao];
    const double Q10NCX = parameters[PARAM_Q10NCX];
    const double ksat = parameters[PARAM_ksat];
    const double nu = parameters[PARAM_nu];
    const double IbarSLCaP = parameters[PARAM_IbarSLCaP];
    const double KmPCa = parameters[PARAM_KmPCa];
    const double Q10SLCaP = parameters[PARAM_Q10SLCaP];
    const double GCaB = parameters[PARAM_GCaB];
    const double Kmf = parameters[PARAM_Kmf];
    const double Kmr = parameters[PARAM_Kmr];
    const double MaxSR = parameters[PARAM_MaxSR];
    const double MinSR = parameters[PARAM_MinSR];
    const double Q10SRCaP = parameters[PARAM_Q10SRCaP];
    const double Vmax_SRCaP = parameters[PARAM_Vmax_SRCaP];
    const double ec50SR = parameters[PARAM_ec50SR];
    const double hillSRCaP = parameters[PARAM_hillSRCaP];
    const double kiCa = parameters[PARAM_kiCa];
    const double kim = parameters[PARAM_kim];
    const double koCa = parameters[PARAM_koCa];
    const double kom = parameters[PARAM_kom];
    const double ks = parameters[PARAM_ks];
    const double Bmax_Naj = parameters[PARAM_Bmax_Naj];
    const double Bmax_Nasl = parameters[PARAM_Bmax_Nasl];
    const double koff_na = parameters[PARAM_koff_na];
    const double kon_na = parameters[PARAM_kon_na];
    const double Bmax_CaM = parameters[PARAM_Bmax_CaM];
    const double Bmax_SR = parameters[PARAM_Bmax_SR];
    const double Bmax_TnChigh = parameters[PARAM_Bmax_TnChigh];
    const double Bmax_TnClow = parameters[PARAM_Bmax_TnClow];
    const double Bmax_myosin = parameters[PARAM_Bmax_myosin];
    const double koff_cam = parameters[PARAM_koff_cam];
    const double koff_myoca = parameters[PARAM_koff_myoca];
    const double koff_myomg = parameters[PARAM_koff_myomg];
    const double koff_sr = parameters[PARAM_koff_sr];
    const double koff_tnchca = parameters[PARAM_koff_tnchca];
    const double koff_tnchmg = parameters[PARAM_koff_tnchmg];
    const double koff_tncl = parameters[PARAM_koff_tncl];
    const double kon_cam = parameters[PARAM_kon_cam];
    const double kon_myoca = parameters[PARAM_kon_myoca];
    const double kon_myomg = parameters[PARAM_kon_myomg];
    const double kon_sr = parameters[PARAM_kon_sr];
    const double kon_tnchca = parameters[PARAM_kon_tnchca];
    const double kon_tnchmg = parameters[PARAM_kon_tnchmg];
    const double kon_tncl = parameters[PARAM_kon_tncl];
    const double Bmax_SLhighj0 = parameters[PARAM_Bmax_SLhighj0];
    const double Bmax_SLhighsl0 = parameters[PARAM_Bmax_SLhighsl0];
    const double Bmax_SLlowj0 = parameters[PARAM_Bmax_SLlowj0];
    const double Bmax_SLlowsl0 = parameters[PARAM_Bmax_SLlowsl0];
    const double koff_slh = parameters[PARAM_koff_slh];
    const double koff_sll = parameters[PARAM_koff_sll];
    const double kon_slh = parameters[PARAM_kon_slh];
    const double kon_sll = parameters[PARAM_kon_sll];
    const double Bmax_Csqn0 = parameters[PARAM_Bmax_Csqn0];
    const double J_ca_juncsl = parameters[PARAM_J_ca_juncsl];
    const double J_ca_slmyo = parameters[PARAM_J_ca_slmyo];
    const double koff_csqn = parameters[PARAM_koff_csqn];
    const double kon_csqn = parameters[PARAM_kon_csqn];
    const double J_na_juncsl = parameters[PARAM_J_na_juncsl];
    const double J_na_slmyo = parameters[PARAM_J_na_slmyo];
    const double Nao = parameters[PARAM_Nao];
    const double Ko = parameters[PARAM_Ko];
    const double Cao = parameters[PARAM_Cao];
    const double Cli = parameters[PARAM_Cli];
    const double Clo = parameters[PARAM_Clo];
    const double Mgi = parameters[PARAM_Mgi];
    const double Cmem = parameters[PARAM_Cmem];
    const double Frdy = parameters[PARAM_Frdy];
    const double R = parameters[PARAM_R];
    const double Temp = parameters[PARAM_Temp];
    const double stim_amplitude = parameters[PARAM_stim_amplitude];
    const double stim_duration = parameters[PARAM_stim_duration];
    const double stim_period = parameters[PARAM_stim_period];
    const double stim_start = parameters[PARAM_stim_start];

    // Expressions for the Geometry component
    const double Vcell = 1.0e-15 * M_PI * cellLength * (cellRadius * cellRadius);
    const double Vmyo = 0.65 * Vcell;
    const double Vsr = 0.035 * Vcell;
    const double Vsl = 0.02 * Vcell;
    const double Vjunc = 0.000539 * Vcell;
    const double Fsl = 1. - Fjunc;
    const double Fsl_CaL = 1. - Fjunc_CaL;

    // Expressions for the Reversal potentials component
    const double FoRT = Frdy / (R * Temp);
    const double ena_junc = log(Nao / Na_j) / FoRT;
    const double ena_sl = log(Nao / Na_sl) / FoRT;
    const double ek = log(Ko / K_i) / FoRT;
    const double eca_junc = log(Cao / Ca_j) / (2. * FoRT);
    const double eca_sl = log(Cao / Ca_sl) / (2. * FoRT);
    const double ecl = log(Cli / Clo) / FoRT;
    const double Qpow = -31. + Temp / 10.;

    // Expressions for the I_Na component
    const double mss = 1.0
                       / ((1. + 0.00184221158116513 * exp(-0.110741971207087 * V_m))
                          * (1. + 0.00184221158116513 * exp(-0.110741971207087 * V_m)));
    const double taum = 0.1292
                                * exp(-((2.94658944658945 + 0.0643500643500644 * V_m)
                                        * (2.94658944658945 + 0.0643500643500644 * V_m)))
                        + 0.06487
                                  * exp(-((-0.0943466353677621 + 0.0195618153364632 * V_m)
                                          * (-0.0943466353677621 + 0.0195618153364632 * V_m)));
    const double ah = (V_m >= -40. ? 0. : 4.43126792958051e-7 * exp(-0.147058823529412 * V_m));
    const double bh =
            (V_m >= -40. ? 0.77 / (0.13 + 0.0497581410839387 * exp(-0.0900900900900901 * V_m))
                         : 310000.0 * exp(0.3485 * V_m) + 2.7 * exp(0.079 * V_m));
    const double tauh = 1.0 / (ah + bh);
    const double hss = 1.0
                       / ((1. + 15212.5932856544 * exp(0.134589502018843 * V_m))
                          * (1. + 15212.5932856544 * exp(0.134589502018843 * V_m)));
    const double aj =
            (V_m >= -40. ? 0.
                         : (37.78 + V_m)
                                   * (-25428.0 * exp(0.2444 * V_m) - 6.948e-6 * exp(-0.04391 * V_m))
                                   / (1. + 50262745825.954 * exp(0.311 * V_m)));
    const double bj =
            (V_m >= -40. ? 0.6 * exp(0.057 * V_m) / (1. + 0.0407622039783662 * exp(-0.1 * V_m))
                         : 0.02424 * exp(-0.01052 * V_m)
                                   / (1. + 0.00396086833990426 * exp(-0.1378 * V_m)));
    const double tauj = 1.0 / (aj + bj);
    const double jss = 1.0
                       / ((1. + 15212.5932856544 * exp(0.134589502018843 * V_m))
                          * (1. + 15212.5932856544 * exp(0.134589502018843 * V_m)));
    const double dm_dt = (-m + mss) / taum;
    const double dm_dt_linearized = -1. / taum;
    states[STATE_m] = (fabs(dm_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized
                               : dt * dm_dt)
                      + m;
    const double dh_dt = (-h + hss) / tauh;
    const double dh_dt_linearized = -1. / tauh;
    states[STATE_h] = (fabs(dh_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized
                               : dt * dh_dt)
                      + h;
    const double dj_dt = (-j + jss) / tauj;
    const double dj_dt_linearized = -1. / tauj;
    states[STATE_j] = (fabs(dj_dt_linearized) > 1.0e-8
                               ? (-1.0 + exp(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized
                               : dt * dj_dt)
                      + j;
    const double I_Na_junc = Fjunc * GNa * (m * m * m) * (-ena_junc + V_m) * h * j;
    const double I_Na_sl = GNa * (m * m * m) * (-ena_sl + V_m) * Fsl * h * j;

    // Expressions for the I_NaBK component
    const double I_nabk_junc = Fjunc * GNaB * (-ena_junc + V_m);
    const double I_nabk_sl = GNaB * (-ena_sl + V_m) * Fsl;

    // Expressions for the I_NaK component
    const double sigma = -1. / 7. + exp(0.0148588410104012 * Nao) / 7.;
    const double fnak =
            1.0 / (1. + 0.1245 * exp(-0.1 * FoRT * V_m) + 0.0365 * exp(-FoRT * V_m) * sigma);
    const double I_nak_junc =
            Fjunc * IbarNaK * Ko * fnak / ((1. + pow(KmNaip, 4.) / pow(Na_j, 4.)) * (KmKo + Ko));
    const double I_nak_sl =
            IbarNaK * Ko * Fsl * fnak / ((1. + pow(KmNaip, 4.) / pow(Na_sl, 4.)) * (KmKo + Ko));
    const double I_nak = I_nak_junc + I_nak_sl;

    // Expressions for the I_Kr component
    const double gkr = 0.430331482911935 * GKr * sqrt(Ko);
    const double xrss = 1.0 / (1. + exp(-2. - V_m / 5.));
    const double tauxr =
            230. / (1. + exp(2. + V_m / 20.))
            + 3300. / ((1. + exp(-22. / 9. - V_m / 9.)) * (1. + exp(11. / 9. + V_m / 9.)));
    const double dx_kr_dt = (-x_kr + xrss) / tauxr;
    const double dx_kr_dt_linearized = -1. / tauxr;
    states[STATE_x_kr] =
            (fabs(dx_kr_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dx_kr_dt_linearized)) * dx_kr_dt / dx_kr_dt_linearized
                     : dt * dx_kr_dt)
            + x_kr;
    const double rkr = 1.0 / (1. + exp(37. / 12. + V_m / 24.));
    const double I_kr = (-ek + V_m) * gkr * rkr * x_kr;

    // Expressions for the I_Kp component
    const double kp_kp = 1.0 / (1. + 1786.47556537862 * exp(-0.167224080267559 * V_m));
    const double I_kp_junc = Fjunc * GKp * (-ek + V_m) * kp_kp;
    const double I_kp_sl = GKp * (-ek + V_m) * Fsl * kp_kp;
    const double I_kp = I_kp_junc + I_kp_sl;

    // Expressions for the I_Ks component
    const double eks = log((Ko + Nao * pNaK) / (pNaK * Na_i + K_i)) / FoRT;
    const double gks_junc = GKs;
    const double gks_sl = GKs;
    const double xsss = 1.0 / (1. + 0.765928338364649 * exp(-0.0701754385964912 * V_m));
    const double tauxs = 990.1 / (1. + 0.841540408868102 * exp(-0.0708215297450425 * V_m));
    const double dx_ks_dt = (-x_ks + xsss) / tauxs;
    const double dx_ks_dt_linearized = -1. / tauxs;
    states[STATE_x_ks] =
            (fabs(dx_ks_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dx_ks_dt_linearized)) * dx_ks_dt / dx_ks_dt_linearized
                     : dt * dx_ks_dt)
            + x_ks;
    const double I_ks_junc = Fjunc * (x_ks * x_ks) * (-eks + V_m) * gks_junc;
    const double I_ks_sl = (x_ks * x_ks) * (-eks + V_m) * Fsl * gks_sl;
    const double I_ks = I_ks_junc + I_ks_sl;

    // Expressions for the I_to component
    const double GtoSlow = (epi == 1. ? 0.12 * Gto : 0.2892 * Gto);
    const double GtoFast = (epi == 1. ? 0.88 * Gto : 0.0108 * Gto);
    const double xtoss = 1.0 / (1. + exp(19. / 13. - V_m / 13.));
    const double ytoss = 1.0 / (1. + 49.4024491055302 * exp(V_m / 5.));
    const double tauxtos = 0.5 + 9. / (1. + exp(1. / 5. + V_m / 15.));
    const double tauytos = 30. + 800. / (1. + exp(6. + V_m / 10.));
    const double dx_to_s_dt = (-x_to_s + xtoss) / tauxtos;
    const double dx_to_s_dt_linearized = -1. / tauxtos;
    states[STATE_x_to_s] =
            (fabs(dx_to_s_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dx_to_s_dt_linearized)) * dx_to_s_dt / dx_to_s_dt_linearized
                     : dt * dx_to_s_dt)
            + x_to_s;
    const double dy_to_s_dt = (-y_to_s + ytoss) / tauytos;
    const double dy_to_s_dt_linearized = -1. / tauytos;
    states[STATE_y_to_s] =
            (fabs(dy_to_s_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dy_to_s_dt_linearized)) * dy_to_s_dt / dy_to_s_dt_linearized
                     : dt * dy_to_s_dt)
            + y_to_s;
    const double I_tos = (-ek + V_m) * GtoSlow * x_to_s * y_to_s;
    const double tauxtof = 0.5 + 8.5 * exp(-((9. / 10. + V_m / 50.) * (9. / 10. + V_m / 50.)));
    const double tauytof = 7. + 85. * exp(-((40. + V_m) * (40. + V_m)) / 220.);
    const double dx_to_f_dt = (-x_to_f + xtoss) / tauxtof;
    const double dx_to_f_dt_linearized = -1. / tauxtof;
    states[STATE_x_to_f] =
            (fabs(dx_to_f_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dx_to_f_dt_linearized)) * dx_to_f_dt / dx_to_f_dt_linearized
                     : dt * dx_to_f_dt)
            + x_to_f;
    const double dy_to_f_dt = (-y_to_f + ytoss) / tauytof;
    const double dy_to_f_dt_linearized = -1. / tauytof;
    states[STATE_y_to_f] =
            (fabs(dy_to_f_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dy_to_f_dt_linearized)) * dy_to_f_dt / dy_to_f_dt_linearized
                     : dt * dy_to_f_dt)
            + y_to_f;
    const double I_tof = (-ek + V_m) * GtoFast * x_to_f * y_to_f;
    const double I_to = I_tof + I_tos;

    // Expressions for the I_K1 component
    const double aki = 1.02 / (1. + 7.35454251046446e-7 * exp(0.2385 * V_m - 0.2385 * ek));
    const double bki = (0.762624006506308 * exp(0.08032 * V_m - 0.08032 * ek)
                        + 1.15340563518656e-16 * exp(0.06175 * V_m - 0.06175 * ek))
                       / (1. + 0.0867722941576933 * exp(0.5143 * ek - 0.5143 * V_m));
    const double kiss = aki / (aki + bki);
    const double I_K1 = 0.430331482911935 * GK1 * sqrt(Ko) * (-ek + V_m) * kiss;

    // Expressions for the I_ClCa component
    const double I_ClCa_junc = Fjunc * GClCa * (-ecl + V_m) / (1. + KdClCa / Ca_j);
    const double I_ClCa_sl = GClCa * (-ecl + V_m) * Fsl / (1. + KdClCa / Ca_sl);
    const double I_ClCa = I_ClCa_junc + I_ClCa_sl;
    const double I_Clbk = GClB * (-ecl + V_m);

    // Expressions for the I_Ca component
    const double fss =
            1.0 / (1. + exp(35. / 9. + V_m / 9.)) + 0.6 / (1. + exp(5. / 2. - V_m / 20.));
    const double dss = 1.0 / (1. + exp(-5. / 6. - V_m / 6.));
    const double taud = (1. - exp(-5. / 6. - V_m / 6.)) * dss / (0.175 + 0.035 * V_m);
    const double tauf =
            1.0 / (0.02 + 0.0197 * exp(-((0.48865 + 0.0337 * V_m) * (0.48865 + 0.0337 * V_m))));
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
    const double df_Ca_Bj_dt = -0.0119 * f_Ca_Bj + 1.7 * (1. - f_Ca_Bj) * Ca_j;
    const double df_Ca_Bj_dt_linearized = -0.0119 - 1.7 * Ca_j;
    states[STATE_f_Ca_Bj] =
            (fabs(df_Ca_Bj_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * df_Ca_Bj_dt_linearized))
                                                             * df_Ca_Bj_dt / df_Ca_Bj_dt_linearized
                                                   : dt * df_Ca_Bj_dt)
            + f_Ca_Bj;
    const double df_Ca_Bsl_dt = -0.0119 * f_Ca_Bsl + 1.7 * (1. - f_Ca_Bsl) * Ca_sl;
    const double df_Ca_Bsl_dt_linearized = -0.0119 - 1.7 * Ca_sl;
    states[STATE_f_Ca_Bsl] = (fabs(df_Ca_Bsl_dt_linearized) > 1.0e-8
                                      ? (-1.0 + exp(dt * df_Ca_Bsl_dt_linearized)) * df_Ca_Bsl_dt
                                                / df_Ca_Bsl_dt_linearized
                                      : dt * df_Ca_Bsl_dt)
                             + f_Ca_Bsl;
    const double fcaCaMSL = 0.;
    const double fcaCaj = 0.;
    const double ibarca_j = 4. * Frdy * GCaL * pCa
                            * (-0.341 * Cao + 0.341 * Ca_j * exp(2. * FoRT * V_m)) * FoRT * V_m
                            / (-1. + exp(2. * FoRT * V_m));
    const double ibarca_sl = 4. * Frdy * GCaL * pCa
                             * (-0.341 * Cao + 0.341 * Ca_sl * exp(2. * FoRT * V_m)) * FoRT * V_m
                             / (-1. + exp(2. * FoRT * V_m));
    const double ibark = Frdy * GCaL * pK * (-0.75 * Ko + 0.75 * K_i * exp(FoRT * V_m)) * FoRT * V_m
                         / (-1. + exp(FoRT * V_m));
    const double ibarna_j = Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_j * exp(FoRT * V_m)) * FoRT
                            * V_m / (-1. + exp(FoRT * V_m));
    const double ibarna_sl = Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_sl * exp(FoRT * V_m))
                             * FoRT * V_m / (-1. + exp(FoRT * V_m));
    const double I_Ca_junc =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f * ibarca_j;
    const double I_Ca_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f * ibarca_sl;
    const double I_CaK =
            0.45 * pow(Q10CaL, Qpow)
            * (Fjunc_CaL * (1. + fcaCaj - f_Ca_Bj) + (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL) * d * f
            * ibark;
    const double I_CaNa_junc =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f * ibarna_j;
    const double I_CaNa_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f * ibarna_sl;

    // Expressions for the I_NCX component
    const double Ka_junc = 1.0 / (1. + (Kdact * Kdact) / (Ca_j * Ca_j));
    const double Ka_sl = 1.0 / (1. + (Kdact * Kdact) / (Ca_sl * Ca_sl));
    const double s1_junc = Cao * (Na_j * Na_j * Na_j) * exp(nu * FoRT * V_m);
    const double s1_sl = Cao * (Na_sl * Na_sl * Na_sl) * exp(nu * FoRT * V_m);
    const double s2_junc = (Nao * Nao * Nao) * Ca_j * exp((-1. + nu) * FoRT * V_m);
    const double s3_junc =
            Cao * (Na_j * Na_j * Na_j) + KmCao * (Na_j * Na_j * Na_j) + (Nao * Nao * Nao) * Ca_j
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_j * Na_j * Na_j) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + Ca_j / KmCai) * Ca_j;
    const double s2_sl = (Nao * Nao * Nao) * Ca_sl * exp((-1. + nu) * FoRT * V_m);
    const double s3_sl =
            Cao * (Na_sl * Na_sl * Na_sl) + KmCao * (Na_sl * Na_sl * Na_sl)
            + (Nao * Nao * Nao) * Ca_sl
            + KmCai * (Nao * Nao * Nao) * (1. + (Na_sl * Na_sl * Na_sl) / (KmNai * KmNai * KmNai))
            + (KmNao * KmNao * KmNao) * (1. + Ca_sl / KmCai) * Ca_sl;
    const double I_ncx_junc = Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
                              / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc);
    const double I_ncx_sl = IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
                            / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);

    // Expressions for the I_PCa component
    const double I_pca_junc = Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, 1.6)
                              / (pow(KmPCa, 1.6) + pow(Ca_j, 1.6));
    const double I_pca_sl = IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, 1.6) * Fsl
                            / (pow(KmPCa, 1.6) + pow(Ca_sl, 1.6));

    // Expressions for the I_CaBK component
    const double I_cabk_junc = Fjunc * GCaB * (-eca_junc + V_m);
    const double I_cabk_sl = GCaB * (-eca_sl + V_m) * Fsl;

    // Expressions for the SR Fluxes component
    const double kCaSR = MaxSR - (MaxSR - MinSR) / (1. + pow(ec50SR / Ca_sr, 2.5));
    const double koSRCa = koCa / kCaSR;
    const double kiSRCa = kiCa * kCaSR;
    const double RI = 1. - Ry_Ri - Ry_Ro - Ry_Rr;
    const double dRy_Rr_dt =
            kim * RI + kom * Ry_Ro - (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Rr * kiSRCa;
    const double dRy_Rr_dt_linearized = -kim - (Ca_j * Ca_j) * koSRCa - Ca_j * kiSRCa;
    states[STATE_Ry_Rr] =
            (fabs(dRy_Rr_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dRy_Rr_dt_linearized)) * dRy_Rr_dt / dRy_Rr_dt_linearized
                     : dt * dRy_Rr_dt)
            + Ry_Rr;
    const double dRy_Ro_dt =
            kim * Ry_Ri - kom * Ry_Ro + (Ca_j * Ca_j) * Ry_Rr * koSRCa - Ca_j * Ry_Ro * kiSRCa;
    const double dRy_Ro_dt_linearized = -kom - Ca_j * kiSRCa;
    states[STATE_Ry_Ro] =
            (fabs(dRy_Ro_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dRy_Ro_dt_linearized)) * dRy_Ro_dt / dRy_Ro_dt_linearized
                     : dt * dRy_Ro_dt)
            + Ry_Ro;
    const double dRy_Ri_dt =
            -kim * Ry_Ri - kom * Ry_Ri + (Ca_j * Ca_j) * RI * koSRCa + Ca_j * Ry_Ro * kiSRCa;
    const double dRy_Ri_dt_linearized = -kim - kom - (Ca_j * Ca_j) * koSRCa;
    states[STATE_Ry_Ri] =
            (fabs(dRy_Ri_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dRy_Ri_dt_linearized)) * dRy_Ri_dt / dRy_Ri_dt_linearized
                     : dt * dRy_Ri_dt)
            + Ry_Ri;
    const double J_SRCarel = ks * (-Ca_j + Ca_sr) * Ry_Ro;
    const double J_serca = Vmax_SRCaP * pow(Q10SRCaP, Qpow)
                           * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                           / (1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP));
    const double J_SRleak = 5.348e-6 * Ca_sr - 5.348e-6 * Ca_j;

    // Expressions for the Na Buffers component
    const double dNa_Bj_dt = -koff_na * Na_Bj + kon_na * (Bmax_Naj - Na_Bj) * Na_j;
    const double dNa_Bj_dt_linearized = -koff_na - kon_na * Na_j;
    states[STATE_Na_Bj] =
            Na_Bj
            + (fabs(dNa_Bj_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dNa_Bj_dt_linearized)) * dNa_Bj_dt / dNa_Bj_dt_linearized
                       : dt * dNa_Bj_dt);
    const double dNa_Bsl_dt = -koff_na * Na_Bsl + kon_na * (Bmax_Nasl - Na_Bsl) * Na_sl;
    const double dNa_Bsl_dt_linearized = -koff_na - kon_na * Na_sl;
    states[STATE_Na_Bsl] =
            Na_Bsl
            + (fabs(dNa_Bsl_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dNa_Bsl_dt_linearized))
                                                              * dNa_Bsl_dt / dNa_Bsl_dt_linearized
                                                    : dt * dNa_Bsl_dt);

    // Expressions for the Cytosolic Ca Buffers component
    const double dTn_CL_dt = -koff_tncl * Tn_CL + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;
    const double dTn_CL_dt_linearized = -koff_tncl - kon_tncl * Ca_i;
    states[STATE_Tn_CL] =
            (fabs(dTn_CL_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dTn_CL_dt_linearized)) * dTn_CL_dt / dTn_CL_dt_linearized
                     : dt * dTn_CL_dt)
            + Tn_CL;
    const double dTn_CHc_dt =
            -koff_tnchca * Tn_CHc + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i;
    const double dTn_CHc_dt_linearized = -koff_tnchca - kon_tnchca * Ca_i;
    states[STATE_Tn_CHc] =
            (fabs(dTn_CHc_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dTn_CHc_dt_linearized)) * dTn_CHc_dt / dTn_CHc_dt_linearized
                     : dt * dTn_CHc_dt)
            + Tn_CHc;
    const double dTn_CHm_dt =
            -koff_tnchmg * Tn_CHm + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm);
    const double dTn_CHm_dt_linearized = -koff_tnchmg - Mgi * kon_tnchmg;
    states[STATE_Tn_CHm] =
            (fabs(dTn_CHm_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dTn_CHm_dt_linearized)) * dTn_CHm_dt / dTn_CHm_dt_linearized
                     : dt * dTn_CHm_dt)
            + Tn_CHm;
    const double dCaM_dt = -koff_cam * CaM + kon_cam * (Bmax_CaM - CaM) * Ca_i;
    const double dCaM_dt_linearized = -koff_cam - kon_cam * Ca_i;
    states[STATE_CaM] =
            CaM
            + (fabs(dCaM_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dCaM_dt_linearized)) * dCaM_dt / dCaM_dt_linearized
                       : dt * dCaM_dt);
    const double dMyo_c_dt = -koff_myoca * Myo_c + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i;
    const double dMyo_c_dt_linearized = -koff_myoca - kon_myoca * Ca_i;
    states[STATE_Myo_c] =
            Myo_c
            + (fabs(dMyo_c_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dMyo_c_dt_linearized)) * dMyo_c_dt / dMyo_c_dt_linearized
                       : dt * dMyo_c_dt);
    const double dMyo_m_dt = -koff_myomg * Myo_m + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m);
    const double dMyo_m_dt_linearized = -koff_myomg - Mgi * kon_myomg;
    states[STATE_Myo_m] =
            Myo_m
            + (fabs(dMyo_m_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dMyo_m_dt_linearized)) * dMyo_m_dt / dMyo_m_dt_linearized
                       : dt * dMyo_m_dt);
    const double dSRB_dt = -koff_sr * SRB + kon_sr * (Bmax_SR - SRB) * Ca_i;
    const double dSRB_dt_linearized = -koff_sr - kon_sr * Ca_i;
    states[STATE_SRB] = (fabs(dSRB_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dSRB_dt_linearized))
                                                                     * dSRB_dt / dSRB_dt_linearized
                                                           : dt * dSRB_dt)
                        + SRB;
    const double J_CaB_cytosol =
            -koff_cam * CaM - koff_myoca * Myo_c - koff_myomg * Myo_m - koff_sr * SRB
            - koff_tnchca * Tn_CHc - koff_tnchmg * Tn_CHm - koff_tncl * Tn_CL
            + Mgi * kon_myomg * (Bmax_myosin - Myo_c - Myo_m)
            + Mgi * kon_tnchmg * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
            + kon_cam * (Bmax_CaM - CaM) * Ca_i + kon_myoca * (Bmax_myosin - Myo_c - Myo_m) * Ca_i
            + kon_sr * (Bmax_SR - SRB) * Ca_i + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm) * Ca_i
            + kon_tncl * (Bmax_TnClow - Tn_CL) * Ca_i;

    // Expressions for the Junctional and SL Ca Buffers component
    const double Bmax_SLlowsl = Bmax_SLlowsl0 * Vmyo / Vsl;
    const double Bmax_SLlowj = Bmax_SLlowj0 * Vmyo / Vjunc;
    const double Bmax_SLhighsl = Bmax_SLhighsl0 * Vmyo / Vsl;
    const double Bmax_SLhighj = Bmax_SLhighj0 * Vmyo / Vjunc;
    const double dSLL_j_dt = -koff_sll * SLL_j + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
    const double dSLL_j_dt_linearized = -koff_sll - kon_sll * Ca_j;
    states[STATE_SLL_j] =
            (fabs(dSLL_j_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dSLL_j_dt_linearized)) * dSLL_j_dt / dSLL_j_dt_linearized
                     : dt * dSLL_j_dt)
            + SLL_j;
    const double dSLL_sl_dt = -koff_sll * SLL_sl + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;
    const double dSLL_sl_dt_linearized = -koff_sll - kon_sll * Ca_sl;
    states[STATE_SLL_sl] =
            (fabs(dSLL_sl_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dSLL_sl_dt_linearized)) * dSLL_sl_dt / dSLL_sl_dt_linearized
                     : dt * dSLL_sl_dt)
            + SLL_sl;
    const double dSLH_j_dt = -koff_slh * SLH_j + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j;
    const double dSLH_j_dt_linearized = -koff_slh - kon_slh * Ca_j;
    states[STATE_SLH_j] =
            (fabs(dSLH_j_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dSLH_j_dt_linearized)) * dSLH_j_dt / dSLH_j_dt_linearized
                     : dt * dSLH_j_dt)
            + SLH_j;
    const double dSLH_sl_dt = -koff_slh * SLH_sl + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl;
    const double dSLH_sl_dt_linearized = -koff_slh - kon_slh * Ca_sl;
    states[STATE_SLH_sl] =
            (fabs(dSLH_sl_dt_linearized) > 1.0e-8
                     ? (-1.0 + exp(dt * dSLH_sl_dt_linearized)) * dSLH_sl_dt / dSLH_sl_dt_linearized
                     : dt * dSLH_sl_dt)
            + SLH_sl;
    const double J_CaB_junction = -koff_slh * SLH_j - koff_sll * SLL_j
                                  + kon_slh * (-SLH_j + Bmax_SLhighj) * Ca_j
                                  + kon_sll * (-SLL_j + Bmax_SLlowj) * Ca_j;
    const double J_CaB_sl = -koff_slh * SLH_sl - koff_sll * SLL_sl
                            + kon_slh * (-SLH_sl + Bmax_SLhighsl) * Ca_sl
                            + kon_sll * (-SLL_sl + Bmax_SLlowsl) * Ca_sl;

    // Expressions for the SR Ca Concentrations component
    const double Bmax_Csqn = Bmax_Csqn0 * Vmyo / Vsr;
    const double dCsqn_b_dt = -koff_csqn * Csqn_b + kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr;
    const double dCsqn_b_dt_linearized = -koff_csqn - kon_csqn * Ca_sr;
    states[STATE_Csqn_b] =
            Csqn_b
            + (fabs(dCsqn_b_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dCsqn_b_dt_linearized))
                                                              * dCsqn_b_dt / dCsqn_b_dt_linearized
                                                    : dt * dCsqn_b_dt);
    const double dCa_sr_dt = -J_SRCarel + koff_csqn * Csqn_b
                             - kon_csqn * (-Csqn_b + Bmax_Csqn) * Ca_sr - J_SRleak * Vmyo / Vsr
                             + J_serca;
    const double dJ_serca_dCa_sr =
            -Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_sr / Kmr, hillSRCaP)
                    / ((1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP)) * Ca_sr)
            - Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_sr / Kmr, hillSRCaP)
                      * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                      / (((1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP))
                          * (1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP)))
                         * Ca_sr);
    const double dCa_sr_dt_linearized = -kon_csqn * (-Csqn_b + Bmax_Csqn) - ks * Ry_Ro
                                        - 5.348e-6 * Vmyo / Vsr + dJ_serca_dCa_sr;
    states[STATE_Ca_sr] =
            Ca_sr
            + (fabs(dCa_sr_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dCa_sr_dt_linearized)) * dCa_sr_dt / dCa_sr_dt_linearized
                       : dt * dCa_sr_dt);

    // Expressions for the Na Concentrations component
    const double I_Na_tot_junc =
            3. * I_nak_junc + 3. * I_ncx_junc + I_CaNa_junc + I_Na_junc + I_nabk_junc;
    const double I_Na_tot_sl = 3. * I_nak_sl + 3. * I_ncx_sl + I_CaNa_sl + I_Na_sl + I_nabk_sl;
    const double dNa_j_dt = -dNa_Bj_dt + J_na_juncsl * (-Na_j + Na_sl) / Vjunc
                            - Cmem * I_Na_tot_junc / (Frdy * Vjunc);
    const double dI_CaNa_junc_dibarna_j =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f;
    const double dI_Na_junc_dena_junc = -Fjunc * GNa * (m * m * m) * h * j;
    const double dI_nabk_junc_dena_junc = -Fjunc * GNaB;
    const double dI_nak_junc_dNa_j =
            4. * Fjunc * IbarNaK * Ko * pow(KmNaip, 4.) * fnak
            / (((1. + pow(KmNaip, 4.) / pow(Na_j, 4.)) * (1. + pow(KmNaip, 4.) / pow(Na_j, 4.)))
               * (KmKo + Ko) * pow(Na_j, 5.));
    const double dI_ncx_junc_ds1_junc = Fjunc * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc
                                        / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc);
    const double dI_ncx_junc_ds3_junc =
            -Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc) * Ka_junc
            / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * (s3_junc * s3_junc));
    const double dena_junc_dNa_j = -1. / (FoRT * Na_j);
    const double dibarna_j_dNa_j =
            0.75 * Frdy * GCaL * pNa * FoRT * V_m * exp(FoRT * V_m) / (-1. + exp(FoRT * V_m));
    const double ds1_junc_dNa_j = 3. * Cao * (Na_j * Na_j) * exp(nu * FoRT * V_m);
    const double ds3_junc_dNa_j =
            3. * Cao * (Na_j * Na_j) + 3. * KmCao * (Na_j * Na_j)
            + 3. * KmCai * (Nao * Nao * Nao) * (Na_j * Na_j) / (KmNai * KmNai * KmNai);
    const double dNa_j_dt_linearized =
            -J_na_juncsl / Vjunc
            - Cmem
                      * (3. * dI_nak_junc_dNa_j + dI_CaNa_junc_dibarna_j * dibarna_j_dNa_j
                         + dI_Na_junc_dena_junc * dena_junc_dNa_j
                         + dI_nabk_junc_dena_junc * dena_junc_dNa_j
                         + 3. * dI_ncx_junc_ds1_junc * ds1_junc_dNa_j
                         + 3. * dI_ncx_junc_ds3_junc * ds3_junc_dNa_j)
                      / (Frdy * Vjunc);
    states[STATE_Na_j] =
            Na_j
            + (fabs(dNa_j_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dNa_j_dt_linearized)) * dNa_j_dt / dNa_j_dt_linearized
                       : dt * dNa_j_dt);
    const double dNa_sl_dt = -dNa_Bsl_dt + J_na_juncsl * (-Na_sl + Na_j) / Vsl
                             + J_na_slmyo * (-Na_sl + Na_i) / Vsl
                             - Cmem * I_Na_tot_sl / (Frdy * Vsl);
    const double dI_CaNa_sl_dibarna_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f;
    const double dI_Na_sl_dena_sl = -GNa * (m * m * m) * Fsl * h * j;
    const double dI_nabk_sl_dena_sl = -GNaB * Fsl;
    const double dI_nak_sl_dNa_sl =
            4. * IbarNaK * Ko * pow(KmNaip, 4.) * Fsl * fnak
            / (((1. + pow(KmNaip, 4.) / pow(Na_sl, 4.)) * (1. + pow(KmNaip, 4.) / pow(Na_sl, 4.)))
               * (KmKo + Ko) * pow(Na_sl, 5.));
    const double dI_ncx_sl_ds1_sl = IbarNCX * pow(Q10NCX, Qpow) * Fsl * Ka_sl
                                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_ncx_sl_ds3_sl =
            -IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl * Ka_sl
            / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * (s3_sl * s3_sl));
    const double dena_sl_dNa_sl = -1. / (FoRT * Na_sl);
    const double dibarna_sl_dNa_sl =
            0.75 * Frdy * GCaL * pNa * FoRT * V_m * exp(FoRT * V_m) / (-1. + exp(FoRT * V_m));
    const double ds1_sl_dNa_sl = 3. * Cao * (Na_sl * Na_sl) * exp(nu * FoRT * V_m);
    const double ds3_sl_dNa_sl =
            3. * Cao * (Na_sl * Na_sl) + 3. * KmCao * (Na_sl * Na_sl)
            + 3. * KmCai * (Nao * Nao * Nao) * (Na_sl * Na_sl) / (KmNai * KmNai * KmNai);
    const double dNa_sl_dt_linearized =
            -J_na_juncsl / Vsl - J_na_slmyo / Vsl
            - Cmem
                      * (3. * dI_nak_sl_dNa_sl + dI_CaNa_sl_dibarna_sl * dibarna_sl_dNa_sl
                         + dI_Na_sl_dena_sl * dena_sl_dNa_sl + dI_nabk_sl_dena_sl * dena_sl_dNa_sl
                         + 3. * dI_ncx_sl_ds1_sl * ds1_sl_dNa_sl
                         + 3. * dI_ncx_sl_ds3_sl * ds3_sl_dNa_sl)
                      / (Frdy * Vsl);
    states[STATE_Na_sl] =
            Na_sl
            + (fabs(dNa_sl_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dNa_sl_dt_linearized)) * dNa_sl_dt / dNa_sl_dt_linearized
                       : dt * dNa_sl_dt);
    const double dNa_i_dt = J_na_slmyo * (-Na_i + Na_sl) / Vmyo;
    const double dNa_i_dt_linearized = -J_na_slmyo / Vmyo;
    states[STATE_Na_i] =
            Na_i
            + (fabs(dNa_i_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dNa_i_dt_linearized)) * dNa_i_dt / dNa_i_dt_linearized
                       : dt * dNa_i_dt);

    // Expressions for the K Concentration component
    const double I_K_tot = -2. * I_nak + I_CaK + I_K1 + I_kp + I_kr + I_ks + I_to;
    const double dK_i_dt = 0.;
    states[STATE_K_i] = dt * dK_i_dt + K_i;

    // Expressions for the Ca Concentrations component
    const double I_Ca_tot_junc = -2. * I_ncx_junc + I_Ca_junc + I_cabk_junc + I_pca_junc;
    const double I_Ca_tot_sl = -2. * I_ncx_sl + I_Ca_sl + I_cabk_sl + I_pca_sl;
    const double dCa_j_dt = -J_CaB_junction + J_ca_juncsl * (-Ca_j + Ca_sl) / Vjunc
                            + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc
                            - Cmem * I_Ca_tot_junc / (2. * Frdy * Vjunc);
    const double dI_Ca_junc_dibarca_j =
            0.45 * Fjunc_CaL * pow(Q10CaL, Qpow) * (1. + fcaCaj - f_Ca_Bj) * d * f;
    const double dI_cabk_junc_deca_junc = -Fjunc * GCaB;
    const double dI_ncx_junc_dKa_junc = Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-s2_junc + s1_junc)
                                        / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc);
    const double dI_ncx_junc_ds2_junc = -Fjunc * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc
                                        / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc);
    const double dI_pca_junc_dCa_j =
            1.6 * Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, 0.6)
                    / (pow(KmPCa, 1.6) + pow(Ca_j, 1.6))
            - 1.6 * Fjunc * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_j, 2.2)
                      / ((pow(KmPCa, 1.6) + pow(Ca_j, 1.6)) * (pow(KmPCa, 1.6) + pow(Ca_j, 1.6)));
    const double dJ_CaB_junction_dCa_j =
            kon_slh * (-SLH_j + Bmax_SLhighj) + kon_sll * (-SLL_j + Bmax_SLlowj);
    const double dJ_SRCarel_dCa_j = -ks * Ry_Ro;
    const double dKa_junc_dCa_j =
            2. * (Kdact * Kdact)
            / (((1. + (Kdact * Kdact) / (Ca_j * Ca_j)) * (1. + (Kdact * Kdact) / (Ca_j * Ca_j)))
               * (Ca_j * Ca_j * Ca_j));
    const double deca_junc_dCa_j = -1. / (2. * Ca_j * FoRT);
    const double dibarca_j_dCa_j = 1.364 * Frdy * GCaL * pCa * FoRT * V_m * exp(2. * FoRT * V_m)
                                   / (-1. + exp(2. * FoRT * V_m));
    const double ds2_junc_dCa_j = (Nao * Nao * Nao) * exp((-1. + nu) * FoRT * V_m);
    const double ds3_junc_dCa_j = (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (1. + Ca_j / KmCai)
                                  + (KmNao * KmNao * KmNao) * Ca_j / KmCai;
    const double dCa_j_dt_linearized =
            -dJ_CaB_junction_dCa_j - J_ca_juncsl / Vjunc - 5.348e-6 * Vmyo / Vjunc
            + Vsr * dJ_SRCarel_dCa_j / Vjunc
            - Cmem
                      * (dI_Ca_junc_dibarca_j * dibarca_j_dCa_j
                         + dI_cabk_junc_deca_junc * deca_junc_dCa_j
                         - 2. * dI_ncx_junc_dKa_junc * dKa_junc_dCa_j
                         - 2. * dI_ncx_junc_ds2_junc * ds2_junc_dCa_j
                         - 2. * dI_ncx_junc_ds3_junc * ds3_junc_dCa_j + dI_pca_junc_dCa_j)
                      / (2. * Frdy * Vjunc);
    states[STATE_Ca_j] =
            Ca_j
            + (fabs(dCa_j_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dCa_j_dt_linearized)) * dCa_j_dt / dCa_j_dt_linearized
                       : dt * dCa_j_dt);
    const double dCa_sl_dt = -J_CaB_sl + J_ca_juncsl * (-Ca_sl + Ca_j) / Vsl
                             + J_ca_slmyo * (-Ca_sl + Ca_i) / Vsl
                             - Cmem * I_Ca_tot_sl / (2. * Frdy * Vsl);
    const double dI_Ca_sl_dibarca_sl =
            0.45 * pow(Q10CaL, Qpow) * (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL * d * f;
    const double dI_cabk_sl_deca_sl = -GCaB * Fsl;
    const double dI_ncx_sl_dKa_sl = IbarNCX * pow(Q10NCX, Qpow) * (-s2_sl + s1_sl) * Fsl
                                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_ncx_sl_ds2_sl = -IbarNCX * pow(Q10NCX, Qpow) * Fsl * Ka_sl
                                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl);
    const double dI_pca_sl_dCa_sl =
            1.6 * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, 0.6) * Fsl
                    / (pow(KmPCa, 1.6) + pow(Ca_sl, 1.6))
            - 1.6 * IbarSLCaP * pow(Q10SLCaP, Qpow) * pow(Ca_sl, 2.2) * Fsl
                      / ((pow(KmPCa, 1.6) + pow(Ca_sl, 1.6)) * (pow(KmPCa, 1.6) + pow(Ca_sl, 1.6)));
    const double dJ_CaB_sl_dCa_sl =
            kon_slh * (-SLH_sl + Bmax_SLhighsl) + kon_sll * (-SLL_sl + Bmax_SLlowsl);
    const double dKa_sl_dCa_sl =
            2. * (Kdact * Kdact)
            / (((1. + (Kdact * Kdact) / (Ca_sl * Ca_sl)) * (1. + (Kdact * Kdact) / (Ca_sl * Ca_sl)))
               * (Ca_sl * Ca_sl * Ca_sl));
    const double deca_sl_dCa_sl = -1. / (2. * Ca_sl * FoRT);
    const double dibarca_sl_dCa_sl = 1.364 * Frdy * GCaL * pCa * FoRT * V_m * exp(2. * FoRT * V_m)
                                     / (-1. + exp(2. * FoRT * V_m));
    const double ds2_sl_dCa_sl = (Nao * Nao * Nao) * exp((-1. + nu) * FoRT * V_m);
    const double ds3_sl_dCa_sl = (Nao * Nao * Nao) + (KmNao * KmNao * KmNao) * (1. + Ca_sl / KmCai)
                                 + (KmNao * KmNao * KmNao) * Ca_sl / KmCai;
    const double dCa_sl_dt_linearized =
            -dJ_CaB_sl_dCa_sl - J_ca_juncsl / Vsl - J_ca_slmyo / Vsl
            - Cmem
                      * (dI_Ca_sl_dibarca_sl * dibarca_sl_dCa_sl
                         + dI_cabk_sl_deca_sl * deca_sl_dCa_sl
                         - 2. * dI_ncx_sl_dKa_sl * dKa_sl_dCa_sl
                         - 2. * dI_ncx_sl_ds2_sl * ds2_sl_dCa_sl
                         - 2. * dI_ncx_sl_ds3_sl * ds3_sl_dCa_sl + dI_pca_sl_dCa_sl)
                      / (2. * Frdy * Vsl);
    states[STATE_Ca_sl] =
            Ca_sl
            + (fabs(dCa_sl_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dCa_sl_dt_linearized)) * dCa_sl_dt / dCa_sl_dt_linearized
                       : dt * dCa_sl_dt);
    const double dCa_i_dt =
            -J_CaB_cytosol + J_ca_slmyo * (-Ca_i + Ca_sl) / Vmyo - J_serca * Vsr / Vmyo;
    const double dJ_CaB_cytosol_dCa_i =
            kon_cam * (Bmax_CaM - CaM) + kon_myoca * (Bmax_myosin - Myo_c - Myo_m)
            + kon_sr * (Bmax_SR - SRB) + kon_tnchca * (Bmax_TnChigh - Tn_CHc - Tn_CHm)
            + kon_tncl * (Bmax_TnClow - Tn_CL);
    const double dJ_serca_dCa_i =
            Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_i / Kmf, hillSRCaP)
                    / ((1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP)) * Ca_i)
            - Vmax_SRCaP * hillSRCaP * pow(Q10SRCaP, Qpow) * pow(Ca_i / Kmf, hillSRCaP)
                      * (pow(Ca_i / Kmf, hillSRCaP) - pow(Ca_sr / Kmr, hillSRCaP))
                      / (((1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP))
                          * (1. + pow(Ca_i / Kmf, hillSRCaP) + pow(Ca_sr / Kmr, hillSRCaP)))
                         * Ca_i);
    const double dCa_i_dt_linearized =
            -dJ_CaB_cytosol_dCa_i - J_ca_slmyo / Vmyo - Vsr * dJ_serca_dCa_i / Vmyo;
    states[STATE_Ca_i] =
            Ca_i
            + (fabs(dCa_i_dt_linearized) > 1.0e-8
                       ? (-1.0 + exp(dt * dCa_i_dt_linearized)) * dCa_i_dt / dCa_i_dt_linearized
                       : dt * dCa_i_dt);

    // Expressions for the Membrane potential component
    const double i_Stim = (t - stim_period * floor(t / stim_period) <= stim_duration + stim_start
                                           && t - stim_period * floor(t / stim_period) >= stim_start
                                   ? -stim_amplitude
                                   : 0.);
    const double I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
    const double I_Cl_tot = I_ClCa + I_Clbk;
    const double I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
    const double I_tot = I_Ca_tot + I_Cl_tot + I_K_tot + I_Na_tot;
    const double dV_m_dt = -I_tot - i_Stim;
    const double dI_CaK_dibark =
            0.45 * pow(Q10CaL, Qpow)
            * (Fjunc_CaL * (1. + fcaCaj - f_Ca_Bj) + (1. + fcaCaMSL - f_Ca_Bsl) * Fsl_CaL) * d * f;
    const double dI_ClCa_junc_dV_m = Fjunc * GClCa / (1. + KdClCa / Ca_j);
    const double dI_ClCa_sl_dV_m = GClCa * Fsl / (1. + KdClCa / Ca_sl);
    const double daki_dV_m = -1.78913955652069e-7 * exp(0.2385 * V_m - 0.2385 * ek)
                             / ((1. + 7.35454251046446e-7 * exp(0.2385 * V_m - 0.2385 * ek))
                                * (1. + 7.35454251046446e-7 * exp(0.2385 * V_m - 0.2385 * ek)));
    const double dbki_dV_m =
            (7.12227979727698e-18 * exp(0.06175 * V_m - 0.06175 * ek)
             + 0.0612539602025867 * exp(0.08032 * V_m - 0.08032 * ek))
                    / (1. + 0.0867722941576933 * exp(0.5143 * ek - 0.5143 * V_m))
            + 0.0446269908853017
                      * (0.762624006506308 * exp(0.08032 * V_m - 0.08032 * ek)
                         + 1.15340563518656e-16 * exp(0.06175 * V_m - 0.06175 * ek))
                      * exp(0.5143 * ek - 0.5143 * V_m)
                      / ((1. + 0.0867722941576933 * exp(0.5143 * ek - 0.5143 * V_m))
                         * (1. + 0.0867722941576933 * exp(0.5143 * ek - 0.5143 * V_m)));
    const double dkiss_daki = 1.0 / (aki + bki) - aki / ((aki + bki) * (aki + bki));
    const double dkiss_dbki = -aki / ((aki + bki) * (aki + bki));
    const double dI_K1_dV_m = 0.430331482911935 * GK1 * sqrt(Ko) * kiss
                              + 0.430331482911935 * GK1 * sqrt(Ko) * (-ek + V_m)
                                        * (daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki);
    const double dI_K1_dkiss = 0.430331482911935 * GK1 * sqrt(Ko) * (-ek + V_m);
    const double dI_Na_junc_dV_m = Fjunc * GNa * (m * m * m) * h * j;
    const double dI_Na_sl_dV_m = GNa * (m * m * m) * Fsl * h * j;
    const double dkp_kp_dV_m = 298.741733340907 * exp(-0.167224080267559 * V_m)
                               / ((1. + 1786.47556537862 * exp(-0.167224080267559 * V_m))
                                  * (1. + 1786.47556537862 * exp(-0.167224080267559 * V_m)));
    const double dI_kp_junc_dV_m = Fjunc * GKp * kp_kp + Fjunc * GKp * (-ek + V_m) * dkp_kp_dV_m;
    const double dI_kp_junc_dkp_kp = Fjunc * GKp * (-ek + V_m);
    const double dI_kp_sl_dV_m = GKp * Fsl * kp_kp + GKp * (-ek + V_m) * Fsl * dkp_kp_dV_m;
    const double dI_kp_sl_dkp_kp = GKp * (-ek + V_m) * Fsl;
    const double drkr_dV_m =
            -exp(37. / 12. + V_m / 24.)
            / (24. * ((1. + exp(37. / 12. + V_m / 24.)) * (1. + exp(37. / 12. + V_m / 24.))));
    const double dI_kr_dV_m = gkr * rkr * x_kr + (-ek + V_m) * drkr_dV_m * gkr * x_kr;
    const double dI_kr_drkr = (-ek + V_m) * gkr * x_kr;
    const double dI_ks_junc_dV_m = Fjunc * (x_ks * x_ks) * gks_junc;
    const double dI_ks_sl_dV_m = (x_ks * x_ks) * Fsl * gks_sl;
    const double dI_nak_junc_dfnak =
            Fjunc * IbarNaK * Ko / ((1. + pow(KmNaip, 4.) / pow(Na_j, 4.)) * (KmKo + Ko));
    const double dI_nak_sl_dfnak =
            IbarNaK * Ko * Fsl / ((1. + pow(KmNaip, 4.) / pow(Na_sl, 4.)) * (KmKo + Ko));
    const double ds1_junc_dV_m = Cao * nu * (Na_j * Na_j * Na_j) * FoRT * exp(nu * FoRT * V_m);
    const double ds2_junc_dV_m =
            (Nao * Nao * Nao) * (-1. + nu) * Ca_j * FoRT * exp((-1. + nu) * FoRT * V_m);
    const double dI_ncx_junc_dV_m =
            Fjunc * IbarNCX * pow(Q10NCX, Qpow) * (-ds2_junc_dV_m + ds1_junc_dV_m) * Ka_junc
                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_junc)
            - Fjunc * IbarNCX * ksat * pow(Q10NCX, Qpow) * (-1. + nu) * (-s2_junc + s1_junc) * FoRT
                      * Ka_junc * exp((-1. + nu) * FoRT * V_m)
                      / (((1. + ksat * exp((-1. + nu) * FoRT * V_m))
                          * (1. + ksat * exp((-1. + nu) * FoRT * V_m)))
                         * s3_junc);
    const double ds1_sl_dV_m = Cao * nu * (Na_sl * Na_sl * Na_sl) * FoRT * exp(nu * FoRT * V_m);
    const double ds2_sl_dV_m =
            (Nao * Nao * Nao) * (-1. + nu) * Ca_sl * FoRT * exp((-1. + nu) * FoRT * V_m);
    const double dI_ncx_sl_dV_m =
            IbarNCX * pow(Q10NCX, Qpow) * (-ds2_sl_dV_m + ds1_sl_dV_m) * Fsl * Ka_sl
                    / ((1. + ksat * exp((-1. + nu) * FoRT * V_m)) * s3_sl)
            - IbarNCX * ksat * pow(Q10NCX, Qpow) * (-1. + nu) * (-s2_sl + s1_sl) * FoRT * Fsl
                      * Ka_sl * exp((-1. + nu) * FoRT * V_m)
                      / (((1. + ksat * exp((-1. + nu) * FoRT * V_m))
                          * (1. + ksat * exp((-1. + nu) * FoRT * V_m)))
                         * s3_sl);
    const double dI_tof_dV_m = GtoFast * x_to_f * y_to_f;
    const double dI_tos_dV_m = GtoSlow * x_to_s * y_to_s;
    const double dfnak_dV_m =
            (0.01245 * FoRT * exp(-0.1 * FoRT * V_m) + 0.0365 * FoRT * exp(-FoRT * V_m) * sigma)
            / ((1. + 0.1245 * exp(-0.1 * FoRT * V_m) + 0.0365 * exp(-FoRT * V_m) * sigma)
               * (1. + 0.1245 * exp(-0.1 * FoRT * V_m) + 0.0365 * exp(-FoRT * V_m) * sigma));
    const double dibarca_j_dV_m =
            4. * Frdy * GCaL * pCa * (-0.341 * Cao + 0.341 * Ca_j * exp(2. * FoRT * V_m)) * FoRT
                    / (-1. + exp(2. * FoRT * V_m))
            - 8. * Frdy * GCaL * pCa * (FoRT * FoRT)
                      * (-0.341 * Cao + 0.341 * Ca_j * exp(2. * FoRT * V_m)) * V_m
                      * exp(2. * FoRT * V_m)
                      / ((-1. + exp(2. * FoRT * V_m)) * (-1. + exp(2. * FoRT * V_m)))
            + 2.728 * Frdy * GCaL * pCa * (FoRT * FoRT) * Ca_j * V_m * exp(2. * FoRT * V_m)
                      / (-1. + exp(2. * FoRT * V_m));
    const double dibarca_sl_dV_m =
            4. * Frdy * GCaL * pCa * (-0.341 * Cao + 0.341 * Ca_sl * exp(2. * FoRT * V_m)) * FoRT
                    / (-1. + exp(2. * FoRT * V_m))
            - 8. * Frdy * GCaL * pCa * (FoRT * FoRT)
                      * (-0.341 * Cao + 0.341 * Ca_sl * exp(2. * FoRT * V_m)) * V_m
                      * exp(2. * FoRT * V_m)
                      / ((-1. + exp(2. * FoRT * V_m)) * (-1. + exp(2. * FoRT * V_m)))
            + 2.728 * Frdy * GCaL * pCa * (FoRT * FoRT) * Ca_sl * V_m * exp(2. * FoRT * V_m)
                      / (-1. + exp(2. * FoRT * V_m));
    const double dibark_dV_m =
            Frdy * GCaL * pK * (-0.75 * Ko + 0.75 * K_i * exp(FoRT * V_m)) * FoRT
                    / (-1. + exp(FoRT * V_m))
            - Frdy * GCaL * pK * (FoRT * FoRT) * (-0.75 * Ko + 0.75 * K_i * exp(FoRT * V_m)) * V_m
                      * exp(FoRT * V_m) / ((-1. + exp(FoRT * V_m)) * (-1. + exp(FoRT * V_m)))
            + 0.75 * Frdy * GCaL * pK * (FoRT * FoRT) * K_i * V_m * exp(FoRT * V_m)
                      / (-1. + exp(FoRT * V_m));
    const double dibarna_j_dV_m =
            Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_j * exp(FoRT * V_m)) * FoRT
                    / (-1. + exp(FoRT * V_m))
            - Frdy * GCaL * pNa * (FoRT * FoRT) * (-0.75 * Nao + 0.75 * Na_j * exp(FoRT * V_m))
                      * V_m * exp(FoRT * V_m) / ((-1. + exp(FoRT * V_m)) * (-1. + exp(FoRT * V_m)))
            + 0.75 * Frdy * GCaL * pNa * (FoRT * FoRT) * Na_j * V_m * exp(FoRT * V_m)
                      / (-1. + exp(FoRT * V_m));
    const double dibarna_sl_dV_m =
            Frdy * GCaL * pNa * (-0.75 * Nao + 0.75 * Na_sl * exp(FoRT * V_m)) * FoRT
                    / (-1. + exp(FoRT * V_m))
            - Frdy * GCaL * pNa * (FoRT * FoRT) * (-0.75 * Nao + 0.75 * Na_sl * exp(FoRT * V_m))
                      * V_m * exp(FoRT * V_m) / ((-1. + exp(FoRT * V_m)) * (-1. + exp(FoRT * V_m)))
            + 0.75 * Frdy * GCaL * pNa * (FoRT * FoRT) * Na_sl * V_m * exp(FoRT * V_m)
                      / (-1. + exp(FoRT * V_m));
    const double dV_m_dt_linearized =
            -GClB - dI_ClCa_junc_dV_m - dI_ClCa_sl_dV_m - dI_K1_dV_m - dI_Na_junc_dV_m
            - dI_Na_sl_dV_m - dI_kp_junc_dV_m - dI_kp_sl_dV_m - dI_kr_dV_m - dI_ks_junc_dV_m
            - dI_ks_sl_dV_m - dI_ncx_junc_dV_m - dI_ncx_sl_dV_m - dI_tof_dV_m - dI_tos_dV_m
            - Fjunc * GCaB - Fjunc * GNaB - GCaB * Fsl - GNaB * Fsl
            - (daki_dV_m * dkiss_daki + dbki_dV_m * dkiss_dbki) * dI_K1_dkiss
            - dI_CaK_dibark * dibark_dV_m - dI_CaNa_junc_dibarna_j * dibarna_j_dV_m
            - dI_CaNa_sl_dibarna_sl * dibarna_sl_dV_m - dI_Ca_junc_dibarca_j * dibarca_j_dV_m
            - dI_Ca_sl_dibarca_sl * dibarca_sl_dV_m - dI_kp_junc_dkp_kp * dkp_kp_dV_m
            - dI_kp_sl_dkp_kp * dkp_kp_dV_m - dI_kr_drkr * drkr_dV_m
            - dI_nak_junc_dfnak * dfnak_dV_m - dI_nak_sl_dfnak * dfnak_dV_m
            - dI_ncx_junc_ds1_junc * ds1_junc_dV_m - dI_ncx_junc_ds2_junc * ds2_junc_dV_m
            - dI_ncx_sl_ds1_sl * ds1_sl_dV_m - dI_ncx_sl_ds2_sl * ds2_sl_dV_m;
    states[STATE_V_m] = (fabs(dV_m_dt_linearized) > 1.0e-8 ? (-1.0 + exp(dt * dV_m_dt_linearized))
                                                                     * dV_m_dt / dV_m_dt_linearized
                                                           : dt * dV_m_dt)
                        + V_m;
}
