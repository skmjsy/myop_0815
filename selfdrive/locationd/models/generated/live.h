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
void live_H(double *in_vec, double *out_8023738009988559242);
void live_err_fun(double *nom_x, double *delta_x, double *out_542607325599421320);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2976604130129466173);
void live_H_mod_fun(double *state, double *out_5799357857479944923);
void live_f_fun(double *state, double dt, double *out_4696192010801758683);
void live_F_fun(double *state, double dt, double *out_8028924799786396209);
void live_h_4(double *state, double *unused, double *out_4032574820212496956);
void live_H_4(double *state, double *unused, double *out_1856445282913742426);
void live_h_9(double *state, double *unused, double *out_5969760726097425182);
void live_H_9(double *state, double *unused, double *out_5430773652350705044);
void live_h_10(double *state, double *unused, double *out_452366081758892872);
void live_H_10(double *state, double *unused, double *out_481319251548546161);
void live_h_12(double *state, double *unused, double *out_8247243165387203066);
void live_H_12(double *state, double *unused, double *out_3163011125118219369);
void live_h_31(double *state, double *unused, double *out_785597178875220859);
void live_H_31(double *state, double *unused, double *out_5908574157443233078);
void live_h_32(double *state, double *unused, double *out_3478676145413613168);
void live_H_32(double *state, double *unused, double *out_6736865544058705980);
void live_h_13(double *state, double *unused, double *out_6752979694932810425);
void live_H_13(double *state, double *unused, double *out_4339981498911877711);
void live_h_14(double *state, double *unused, double *out_5969760726097425182);
void live_H_14(double *state, double *unused, double *out_5430773652350705044);
void live_h_33(double *state, double *unused, double *out_1744391063262447564);
void live_H_33(double *state, double *unused, double *out_9059131162082090682);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}