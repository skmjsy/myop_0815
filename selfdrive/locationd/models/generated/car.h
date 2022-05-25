#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5880628066235207371);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1461421936946751991);
void car_H_mod_fun(double *state, double *out_755478003627873076);
void car_f_fun(double *state, double dt, double *out_1213401948881567203);
void car_F_fun(double *state, double dt, double *out_4768455973214330715);
void car_h_25(double *state, double *unused, double *out_9212746903059268727);
void car_H_25(double *state, double *unused, double *out_3862204661873397673);
void car_h_24(double *state, double *unused, double *out_6625563200448412729);
void car_H_24(double *state, double *unused, double *out_5282662048077479489);
void car_h_30(double *state, double *unused, double *out_1319807408174745368);
void car_H_30(double *state, double *unused, double *out_6380537620380646300);
void car_h_26(double *state, double *unused, double *out_2263317261076407803);
void car_H_26(double *state, double *unused, double *out_120701342999341449);
void car_h_27(double *state, double *unused, double *out_5958149357902826300);
void car_H_27(double *state, double *unused, double *out_2840254980054635436);
void car_h_29(double *state, double *unused, double *out_5682955295618320411);
void car_H_29(double *state, double *unused, double *out_155260323939818341);
void car_h_28(double *state, double *unused, double *out_956476090241980655);
void car_H_28(double *state, double *unused, double *out_1808369947625507910);
void car_h_31(double *state, double *unused, double *out_8937552840774762838);
void car_H_31(double *state, double *unused, double *out_505506759234010027);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}