#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5880628066235207371) {
   out_5880628066235207371[0] = delta_x[0] + nom_x[0];
   out_5880628066235207371[1] = delta_x[1] + nom_x[1];
   out_5880628066235207371[2] = delta_x[2] + nom_x[2];
   out_5880628066235207371[3] = delta_x[3] + nom_x[3];
   out_5880628066235207371[4] = delta_x[4] + nom_x[4];
   out_5880628066235207371[5] = delta_x[5] + nom_x[5];
   out_5880628066235207371[6] = delta_x[6] + nom_x[6];
   out_5880628066235207371[7] = delta_x[7] + nom_x[7];
   out_5880628066235207371[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1461421936946751991) {
   out_1461421936946751991[0] = -nom_x[0] + true_x[0];
   out_1461421936946751991[1] = -nom_x[1] + true_x[1];
   out_1461421936946751991[2] = -nom_x[2] + true_x[2];
   out_1461421936946751991[3] = -nom_x[3] + true_x[3];
   out_1461421936946751991[4] = -nom_x[4] + true_x[4];
   out_1461421936946751991[5] = -nom_x[5] + true_x[5];
   out_1461421936946751991[6] = -nom_x[6] + true_x[6];
   out_1461421936946751991[7] = -nom_x[7] + true_x[7];
   out_1461421936946751991[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_755478003627873076) {
   out_755478003627873076[0] = 1.0;
   out_755478003627873076[1] = 0;
   out_755478003627873076[2] = 0;
   out_755478003627873076[3] = 0;
   out_755478003627873076[4] = 0;
   out_755478003627873076[5] = 0;
   out_755478003627873076[6] = 0;
   out_755478003627873076[7] = 0;
   out_755478003627873076[8] = 0;
   out_755478003627873076[9] = 0;
   out_755478003627873076[10] = 1.0;
   out_755478003627873076[11] = 0;
   out_755478003627873076[12] = 0;
   out_755478003627873076[13] = 0;
   out_755478003627873076[14] = 0;
   out_755478003627873076[15] = 0;
   out_755478003627873076[16] = 0;
   out_755478003627873076[17] = 0;
   out_755478003627873076[18] = 0;
   out_755478003627873076[19] = 0;
   out_755478003627873076[20] = 1.0;
   out_755478003627873076[21] = 0;
   out_755478003627873076[22] = 0;
   out_755478003627873076[23] = 0;
   out_755478003627873076[24] = 0;
   out_755478003627873076[25] = 0;
   out_755478003627873076[26] = 0;
   out_755478003627873076[27] = 0;
   out_755478003627873076[28] = 0;
   out_755478003627873076[29] = 0;
   out_755478003627873076[30] = 1.0;
   out_755478003627873076[31] = 0;
   out_755478003627873076[32] = 0;
   out_755478003627873076[33] = 0;
   out_755478003627873076[34] = 0;
   out_755478003627873076[35] = 0;
   out_755478003627873076[36] = 0;
   out_755478003627873076[37] = 0;
   out_755478003627873076[38] = 0;
   out_755478003627873076[39] = 0;
   out_755478003627873076[40] = 1.0;
   out_755478003627873076[41] = 0;
   out_755478003627873076[42] = 0;
   out_755478003627873076[43] = 0;
   out_755478003627873076[44] = 0;
   out_755478003627873076[45] = 0;
   out_755478003627873076[46] = 0;
   out_755478003627873076[47] = 0;
   out_755478003627873076[48] = 0;
   out_755478003627873076[49] = 0;
   out_755478003627873076[50] = 1.0;
   out_755478003627873076[51] = 0;
   out_755478003627873076[52] = 0;
   out_755478003627873076[53] = 0;
   out_755478003627873076[54] = 0;
   out_755478003627873076[55] = 0;
   out_755478003627873076[56] = 0;
   out_755478003627873076[57] = 0;
   out_755478003627873076[58] = 0;
   out_755478003627873076[59] = 0;
   out_755478003627873076[60] = 1.0;
   out_755478003627873076[61] = 0;
   out_755478003627873076[62] = 0;
   out_755478003627873076[63] = 0;
   out_755478003627873076[64] = 0;
   out_755478003627873076[65] = 0;
   out_755478003627873076[66] = 0;
   out_755478003627873076[67] = 0;
   out_755478003627873076[68] = 0;
   out_755478003627873076[69] = 0;
   out_755478003627873076[70] = 1.0;
   out_755478003627873076[71] = 0;
   out_755478003627873076[72] = 0;
   out_755478003627873076[73] = 0;
   out_755478003627873076[74] = 0;
   out_755478003627873076[75] = 0;
   out_755478003627873076[76] = 0;
   out_755478003627873076[77] = 0;
   out_755478003627873076[78] = 0;
   out_755478003627873076[79] = 0;
   out_755478003627873076[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_1213401948881567203) {
   out_1213401948881567203[0] = state[0];
   out_1213401948881567203[1] = state[1];
   out_1213401948881567203[2] = state[2];
   out_1213401948881567203[3] = state[3];
   out_1213401948881567203[4] = state[4];
   out_1213401948881567203[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1213401948881567203[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1213401948881567203[7] = state[7];
   out_1213401948881567203[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4768455973214330715) {
   out_4768455973214330715[0] = 1;
   out_4768455973214330715[1] = 0;
   out_4768455973214330715[2] = 0;
   out_4768455973214330715[3] = 0;
   out_4768455973214330715[4] = 0;
   out_4768455973214330715[5] = 0;
   out_4768455973214330715[6] = 0;
   out_4768455973214330715[7] = 0;
   out_4768455973214330715[8] = 0;
   out_4768455973214330715[9] = 0;
   out_4768455973214330715[10] = 1;
   out_4768455973214330715[11] = 0;
   out_4768455973214330715[12] = 0;
   out_4768455973214330715[13] = 0;
   out_4768455973214330715[14] = 0;
   out_4768455973214330715[15] = 0;
   out_4768455973214330715[16] = 0;
   out_4768455973214330715[17] = 0;
   out_4768455973214330715[18] = 0;
   out_4768455973214330715[19] = 0;
   out_4768455973214330715[20] = 1;
   out_4768455973214330715[21] = 0;
   out_4768455973214330715[22] = 0;
   out_4768455973214330715[23] = 0;
   out_4768455973214330715[24] = 0;
   out_4768455973214330715[25] = 0;
   out_4768455973214330715[26] = 0;
   out_4768455973214330715[27] = 0;
   out_4768455973214330715[28] = 0;
   out_4768455973214330715[29] = 0;
   out_4768455973214330715[30] = 1;
   out_4768455973214330715[31] = 0;
   out_4768455973214330715[32] = 0;
   out_4768455973214330715[33] = 0;
   out_4768455973214330715[34] = 0;
   out_4768455973214330715[35] = 0;
   out_4768455973214330715[36] = 0;
   out_4768455973214330715[37] = 0;
   out_4768455973214330715[38] = 0;
   out_4768455973214330715[39] = 0;
   out_4768455973214330715[40] = 1;
   out_4768455973214330715[41] = 0;
   out_4768455973214330715[42] = 0;
   out_4768455973214330715[43] = 0;
   out_4768455973214330715[44] = 0;
   out_4768455973214330715[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4768455973214330715[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4768455973214330715[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4768455973214330715[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4768455973214330715[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4768455973214330715[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4768455973214330715[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4768455973214330715[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4768455973214330715[53] = -9.8000000000000007*dt;
   out_4768455973214330715[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4768455973214330715[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4768455973214330715[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4768455973214330715[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4768455973214330715[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4768455973214330715[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4768455973214330715[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4768455973214330715[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4768455973214330715[62] = 0;
   out_4768455973214330715[63] = 0;
   out_4768455973214330715[64] = 0;
   out_4768455973214330715[65] = 0;
   out_4768455973214330715[66] = 0;
   out_4768455973214330715[67] = 0;
   out_4768455973214330715[68] = 0;
   out_4768455973214330715[69] = 0;
   out_4768455973214330715[70] = 1;
   out_4768455973214330715[71] = 0;
   out_4768455973214330715[72] = 0;
   out_4768455973214330715[73] = 0;
   out_4768455973214330715[74] = 0;
   out_4768455973214330715[75] = 0;
   out_4768455973214330715[76] = 0;
   out_4768455973214330715[77] = 0;
   out_4768455973214330715[78] = 0;
   out_4768455973214330715[79] = 0;
   out_4768455973214330715[80] = 1;
}
void h_25(double *state, double *unused, double *out_9212746903059268727) {
   out_9212746903059268727[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3862204661873397673) {
   out_3862204661873397673[0] = 0;
   out_3862204661873397673[1] = 0;
   out_3862204661873397673[2] = 0;
   out_3862204661873397673[3] = 0;
   out_3862204661873397673[4] = 0;
   out_3862204661873397673[5] = 0;
   out_3862204661873397673[6] = 1;
   out_3862204661873397673[7] = 0;
   out_3862204661873397673[8] = 0;
}
void h_24(double *state, double *unused, double *out_6625563200448412729) {
   out_6625563200448412729[0] = state[4];
   out_6625563200448412729[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5282662048077479489) {
   out_5282662048077479489[0] = 0;
   out_5282662048077479489[1] = 0;
   out_5282662048077479489[2] = 0;
   out_5282662048077479489[3] = 0;
   out_5282662048077479489[4] = 1;
   out_5282662048077479489[5] = 0;
   out_5282662048077479489[6] = 0;
   out_5282662048077479489[7] = 0;
   out_5282662048077479489[8] = 0;
   out_5282662048077479489[9] = 0;
   out_5282662048077479489[10] = 0;
   out_5282662048077479489[11] = 0;
   out_5282662048077479489[12] = 0;
   out_5282662048077479489[13] = 0;
   out_5282662048077479489[14] = 1;
   out_5282662048077479489[15] = 0;
   out_5282662048077479489[16] = 0;
   out_5282662048077479489[17] = 0;
}
void h_30(double *state, double *unused, double *out_1319807408174745368) {
   out_1319807408174745368[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6380537620380646300) {
   out_6380537620380646300[0] = 0;
   out_6380537620380646300[1] = 0;
   out_6380537620380646300[2] = 0;
   out_6380537620380646300[3] = 0;
   out_6380537620380646300[4] = 1;
   out_6380537620380646300[5] = 0;
   out_6380537620380646300[6] = 0;
   out_6380537620380646300[7] = 0;
   out_6380537620380646300[8] = 0;
}
void h_26(double *state, double *unused, double *out_2263317261076407803) {
   out_2263317261076407803[0] = state[7];
}
void H_26(double *state, double *unused, double *out_120701342999341449) {
   out_120701342999341449[0] = 0;
   out_120701342999341449[1] = 0;
   out_120701342999341449[2] = 0;
   out_120701342999341449[3] = 0;
   out_120701342999341449[4] = 0;
   out_120701342999341449[5] = 0;
   out_120701342999341449[6] = 0;
   out_120701342999341449[7] = 1;
   out_120701342999341449[8] = 0;
}
void h_27(double *state, double *unused, double *out_5958149357902826300) {
   out_5958149357902826300[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2840254980054635436) {
   out_2840254980054635436[0] = 0;
   out_2840254980054635436[1] = 0;
   out_2840254980054635436[2] = 0;
   out_2840254980054635436[3] = 1;
   out_2840254980054635436[4] = 0;
   out_2840254980054635436[5] = 0;
   out_2840254980054635436[6] = 0;
   out_2840254980054635436[7] = 0;
   out_2840254980054635436[8] = 0;
}
void h_29(double *state, double *unused, double *out_5682955295618320411) {
   out_5682955295618320411[0] = state[1];
}
void H_29(double *state, double *unused, double *out_155260323939818341) {
   out_155260323939818341[0] = 0;
   out_155260323939818341[1] = 1;
   out_155260323939818341[2] = 0;
   out_155260323939818341[3] = 0;
   out_155260323939818341[4] = 0;
   out_155260323939818341[5] = 0;
   out_155260323939818341[6] = 0;
   out_155260323939818341[7] = 0;
   out_155260323939818341[8] = 0;
}
void h_28(double *state, double *unused, double *out_956476090241980655) {
   out_956476090241980655[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1808369947625507910) {
   out_1808369947625507910[0] = 1;
   out_1808369947625507910[1] = 0;
   out_1808369947625507910[2] = 0;
   out_1808369947625507910[3] = 0;
   out_1808369947625507910[4] = 0;
   out_1808369947625507910[5] = 0;
   out_1808369947625507910[6] = 0;
   out_1808369947625507910[7] = 0;
   out_1808369947625507910[8] = 0;
}
void h_31(double *state, double *unused, double *out_8937552840774762838) {
   out_8937552840774762838[0] = state[8];
}
void H_31(double *state, double *unused, double *out_505506759234010027) {
   out_505506759234010027[0] = 0;
   out_505506759234010027[1] = 0;
   out_505506759234010027[2] = 0;
   out_505506759234010027[3] = 0;
   out_505506759234010027[4] = 0;
   out_505506759234010027[5] = 0;
   out_505506759234010027[6] = 0;
   out_505506759234010027[7] = 0;
   out_505506759234010027[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5880628066235207371) {
  err_fun(nom_x, delta_x, out_5880628066235207371);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1461421936946751991) {
  inv_err_fun(nom_x, true_x, out_1461421936946751991);
}
void car_H_mod_fun(double *state, double *out_755478003627873076) {
  H_mod_fun(state, out_755478003627873076);
}
void car_f_fun(double *state, double dt, double *out_1213401948881567203) {
  f_fun(state,  dt, out_1213401948881567203);
}
void car_F_fun(double *state, double dt, double *out_4768455973214330715) {
  F_fun(state,  dt, out_4768455973214330715);
}
void car_h_25(double *state, double *unused, double *out_9212746903059268727) {
  h_25(state, unused, out_9212746903059268727);
}
void car_H_25(double *state, double *unused, double *out_3862204661873397673) {
  H_25(state, unused, out_3862204661873397673);
}
void car_h_24(double *state, double *unused, double *out_6625563200448412729) {
  h_24(state, unused, out_6625563200448412729);
}
void car_H_24(double *state, double *unused, double *out_5282662048077479489) {
  H_24(state, unused, out_5282662048077479489);
}
void car_h_30(double *state, double *unused, double *out_1319807408174745368) {
  h_30(state, unused, out_1319807408174745368);
}
void car_H_30(double *state, double *unused, double *out_6380537620380646300) {
  H_30(state, unused, out_6380537620380646300);
}
void car_h_26(double *state, double *unused, double *out_2263317261076407803) {
  h_26(state, unused, out_2263317261076407803);
}
void car_H_26(double *state, double *unused, double *out_120701342999341449) {
  H_26(state, unused, out_120701342999341449);
}
void car_h_27(double *state, double *unused, double *out_5958149357902826300) {
  h_27(state, unused, out_5958149357902826300);
}
void car_H_27(double *state, double *unused, double *out_2840254980054635436) {
  H_27(state, unused, out_2840254980054635436);
}
void car_h_29(double *state, double *unused, double *out_5682955295618320411) {
  h_29(state, unused, out_5682955295618320411);
}
void car_H_29(double *state, double *unused, double *out_155260323939818341) {
  H_29(state, unused, out_155260323939818341);
}
void car_h_28(double *state, double *unused, double *out_956476090241980655) {
  h_28(state, unused, out_956476090241980655);
}
void car_H_28(double *state, double *unused, double *out_1808369947625507910) {
  H_28(state, unused, out_1808369947625507910);
}
void car_h_31(double *state, double *unused, double *out_8937552840774762838) {
  h_31(state, unused, out_8937552840774762838);
}
void car_H_31(double *state, double *unused, double *out_505506759234010027) {
  H_31(state, unused, out_505506759234010027);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
