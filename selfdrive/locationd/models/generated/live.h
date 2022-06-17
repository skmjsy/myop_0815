#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1005564015022223644);
void live_err_fun(double *nom_x, double *delta_x, double *out_7139359416976386283);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3528438406138962813);
void live_H_mod_fun(double *state, double *out_857020396668958774);
void live_f_fun(double *state, double dt, double *out_8594626317175388658);
void live_F_fun(double *state, double dt, double *out_2686999570181155255);
void live_h_4(double *state, double *unused, double *out_1981103364723989574);
void live_H_4(double *state, double *unused, double *out_2811390932214119266);
void live_h_9(double *state, double *unused, double *out_4657013390034409161);
void live_H_9(double *state, double *unused, double *out_8348134206230984880);
void live_h_10(double *state, double *unused, double *out_2164776807616556046);
void live_H_10(double *state, double *unused, double *out_244317375542918480);
void live_h_12(double *state, double *unused, double *out_6579105792362423226);
void live_H_12(double *state, double *unused, double *out_7830847340246081061);
void live_h_31(double *state, double *unused, double *out_3986579323893169117);
void live_H_31(double *state, double *unused, double *out_6178052989586726642);
void live_h_32(double *state, double *unused, double *out_7691855483359838423);
void live_H_32(double *state, double *unused, double *out_4597277548981993471);
void live_h_13(double *state, double *unused, double *out_6574883009380718746);
void live_H_13(double *state, double *unused, double *out_5933706523911350354);
void live_h_14(double *state, double *unused, double *out_4657013390034409161);
void live_H_14(double *state, double *unused, double *out_8348134206230984880);
void live_h_33(double *state, double *unused, double *out_8108416934746195994);
void live_H_33(double *state, double *unused, double *out_9118134079483967370);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}