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
void car_err_fun(double *nom_x, double *delta_x, double *out_2128698594186070786);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4181753220305909536);
void car_H_mod_fun(double *state, double *out_6400451532049496356);
void car_f_fun(double *state, double dt, double *out_4626213852693615128);
void car_F_fun(double *state, double dt, double *out_4088273288682318823);
void car_h_25(double *state, double *unused, double *out_1795064971337587010);
void car_H_25(double *state, double *unused, double *out_2881603371350459916);
void car_h_24(double *state, double *unused, double *out_3436917844040241541);
void car_H_24(double *state, double *unused, double *out_6439528970152783052);
void car_h_30(double *state, double *unused, double *out_1842384556916435340);
void car_H_30(double *state, double *unused, double *out_5399936329857708543);
void car_h_26(double *state, double *unused, double *out_5258657846565505356);
void car_H_26(double *state, double *unused, double *out_859899947523596308);
void car_h_27(double *state, double *unused, double *out_2962808772413728747);
void car_H_27(double *state, double *unused, double *out_3820856270577573193);
void car_h_29(double *state, double *unused, double *out_4556656960497413806);
void car_H_29(double *state, double *unused, double *out_1135861614462756098);
void car_h_28(double *state, double *unused, double *out_5775360875683582130);
void car_H_28(double *state, double *unused, double *out_827768657102570153);
void car_h_31(double *state, double *unused, double *out_1193646813724223456);
void car_H_31(double *state, double *unused, double *out_1486108049756947784);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}