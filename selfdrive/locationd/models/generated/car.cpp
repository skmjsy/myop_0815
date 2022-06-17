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
void err_fun(double *nom_x, double *delta_x, double *out_2128698594186070786) {
   out_2128698594186070786[0] = delta_x[0] + nom_x[0];
   out_2128698594186070786[1] = delta_x[1] + nom_x[1];
   out_2128698594186070786[2] = delta_x[2] + nom_x[2];
   out_2128698594186070786[3] = delta_x[3] + nom_x[3];
   out_2128698594186070786[4] = delta_x[4] + nom_x[4];
   out_2128698594186070786[5] = delta_x[5] + nom_x[5];
   out_2128698594186070786[6] = delta_x[6] + nom_x[6];
   out_2128698594186070786[7] = delta_x[7] + nom_x[7];
   out_2128698594186070786[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4181753220305909536) {
   out_4181753220305909536[0] = -nom_x[0] + true_x[0];
   out_4181753220305909536[1] = -nom_x[1] + true_x[1];
   out_4181753220305909536[2] = -nom_x[2] + true_x[2];
   out_4181753220305909536[3] = -nom_x[3] + true_x[3];
   out_4181753220305909536[4] = -nom_x[4] + true_x[4];
   out_4181753220305909536[5] = -nom_x[5] + true_x[5];
   out_4181753220305909536[6] = -nom_x[6] + true_x[6];
   out_4181753220305909536[7] = -nom_x[7] + true_x[7];
   out_4181753220305909536[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_6400451532049496356) {
   out_6400451532049496356[0] = 1.0;
   out_6400451532049496356[1] = 0;
   out_6400451532049496356[2] = 0;
   out_6400451532049496356[3] = 0;
   out_6400451532049496356[4] = 0;
   out_6400451532049496356[5] = 0;
   out_6400451532049496356[6] = 0;
   out_6400451532049496356[7] = 0;
   out_6400451532049496356[8] = 0;
   out_6400451532049496356[9] = 0;
   out_6400451532049496356[10] = 1.0;
   out_6400451532049496356[11] = 0;
   out_6400451532049496356[12] = 0;
   out_6400451532049496356[13] = 0;
   out_6400451532049496356[14] = 0;
   out_6400451532049496356[15] = 0;
   out_6400451532049496356[16] = 0;
   out_6400451532049496356[17] = 0;
   out_6400451532049496356[18] = 0;
   out_6400451532049496356[19] = 0;
   out_6400451532049496356[20] = 1.0;
   out_6400451532049496356[21] = 0;
   out_6400451532049496356[22] = 0;
   out_6400451532049496356[23] = 0;
   out_6400451532049496356[24] = 0;
   out_6400451532049496356[25] = 0;
   out_6400451532049496356[26] = 0;
   out_6400451532049496356[27] = 0;
   out_6400451532049496356[28] = 0;
   out_6400451532049496356[29] = 0;
   out_6400451532049496356[30] = 1.0;
   out_6400451532049496356[31] = 0;
   out_6400451532049496356[32] = 0;
   out_6400451532049496356[33] = 0;
   out_6400451532049496356[34] = 0;
   out_6400451532049496356[35] = 0;
   out_6400451532049496356[36] = 0;
   out_6400451532049496356[37] = 0;
   out_6400451532049496356[38] = 0;
   out_6400451532049496356[39] = 0;
   out_6400451532049496356[40] = 1.0;
   out_6400451532049496356[41] = 0;
   out_6400451532049496356[42] = 0;
   out_6400451532049496356[43] = 0;
   out_6400451532049496356[44] = 0;
   out_6400451532049496356[45] = 0;
   out_6400451532049496356[46] = 0;
   out_6400451532049496356[47] = 0;
   out_6400451532049496356[48] = 0;
   out_6400451532049496356[49] = 0;
   out_6400451532049496356[50] = 1.0;
   out_6400451532049496356[51] = 0;
   out_6400451532049496356[52] = 0;
   out_6400451532049496356[53] = 0;
   out_6400451532049496356[54] = 0;
   out_6400451532049496356[55] = 0;
   out_6400451532049496356[56] = 0;
   out_6400451532049496356[57] = 0;
   out_6400451532049496356[58] = 0;
   out_6400451532049496356[59] = 0;
   out_6400451532049496356[60] = 1.0;
   out_6400451532049496356[61] = 0;
   out_6400451532049496356[62] = 0;
   out_6400451532049496356[63] = 0;
   out_6400451532049496356[64] = 0;
   out_6400451532049496356[65] = 0;
   out_6400451532049496356[66] = 0;
   out_6400451532049496356[67] = 0;
   out_6400451532049496356[68] = 0;
   out_6400451532049496356[69] = 0;
   out_6400451532049496356[70] = 1.0;
   out_6400451532049496356[71] = 0;
   out_6400451532049496356[72] = 0;
   out_6400451532049496356[73] = 0;
   out_6400451532049496356[74] = 0;
   out_6400451532049496356[75] = 0;
   out_6400451532049496356[76] = 0;
   out_6400451532049496356[77] = 0;
   out_6400451532049496356[78] = 0;
   out_6400451532049496356[79] = 0;
   out_6400451532049496356[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4626213852693615128) {
   out_4626213852693615128[0] = state[0];
   out_4626213852693615128[1] = state[1];
   out_4626213852693615128[2] = state[2];
   out_4626213852693615128[3] = state[3];
   out_4626213852693615128[4] = state[4];
   out_4626213852693615128[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4626213852693615128[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4626213852693615128[7] = state[7];
   out_4626213852693615128[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4088273288682318823) {
   out_4088273288682318823[0] = 1;
   out_4088273288682318823[1] = 0;
   out_4088273288682318823[2] = 0;
   out_4088273288682318823[3] = 0;
   out_4088273288682318823[4] = 0;
   out_4088273288682318823[5] = 0;
   out_4088273288682318823[6] = 0;
   out_4088273288682318823[7] = 0;
   out_4088273288682318823[8] = 0;
   out_4088273288682318823[9] = 0;
   out_4088273288682318823[10] = 1;
   out_4088273288682318823[11] = 0;
   out_4088273288682318823[12] = 0;
   out_4088273288682318823[13] = 0;
   out_4088273288682318823[14] = 0;
   out_4088273288682318823[15] = 0;
   out_4088273288682318823[16] = 0;
   out_4088273288682318823[17] = 0;
   out_4088273288682318823[18] = 0;
   out_4088273288682318823[19] = 0;
   out_4088273288682318823[20] = 1;
   out_4088273288682318823[21] = 0;
   out_4088273288682318823[22] = 0;
   out_4088273288682318823[23] = 0;
   out_4088273288682318823[24] = 0;
   out_4088273288682318823[25] = 0;
   out_4088273288682318823[26] = 0;
   out_4088273288682318823[27] = 0;
   out_4088273288682318823[28] = 0;
   out_4088273288682318823[29] = 0;
   out_4088273288682318823[30] = 1;
   out_4088273288682318823[31] = 0;
   out_4088273288682318823[32] = 0;
   out_4088273288682318823[33] = 0;
   out_4088273288682318823[34] = 0;
   out_4088273288682318823[35] = 0;
   out_4088273288682318823[36] = 0;
   out_4088273288682318823[37] = 0;
   out_4088273288682318823[38] = 0;
   out_4088273288682318823[39] = 0;
   out_4088273288682318823[40] = 1;
   out_4088273288682318823[41] = 0;
   out_4088273288682318823[42] = 0;
   out_4088273288682318823[43] = 0;
   out_4088273288682318823[44] = 0;
   out_4088273288682318823[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4088273288682318823[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4088273288682318823[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4088273288682318823[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4088273288682318823[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4088273288682318823[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4088273288682318823[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4088273288682318823[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4088273288682318823[53] = -9.8000000000000007*dt;
   out_4088273288682318823[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4088273288682318823[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4088273288682318823[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4088273288682318823[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4088273288682318823[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4088273288682318823[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4088273288682318823[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4088273288682318823[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4088273288682318823[62] = 0;
   out_4088273288682318823[63] = 0;
   out_4088273288682318823[64] = 0;
   out_4088273288682318823[65] = 0;
   out_4088273288682318823[66] = 0;
   out_4088273288682318823[67] = 0;
   out_4088273288682318823[68] = 0;
   out_4088273288682318823[69] = 0;
   out_4088273288682318823[70] = 1;
   out_4088273288682318823[71] = 0;
   out_4088273288682318823[72] = 0;
   out_4088273288682318823[73] = 0;
   out_4088273288682318823[74] = 0;
   out_4088273288682318823[75] = 0;
   out_4088273288682318823[76] = 0;
   out_4088273288682318823[77] = 0;
   out_4088273288682318823[78] = 0;
   out_4088273288682318823[79] = 0;
   out_4088273288682318823[80] = 1;
}
void h_25(double *state, double *unused, double *out_1795064971337587010) {
   out_1795064971337587010[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2881603371350459916) {
   out_2881603371350459916[0] = 0;
   out_2881603371350459916[1] = 0;
   out_2881603371350459916[2] = 0;
   out_2881603371350459916[3] = 0;
   out_2881603371350459916[4] = 0;
   out_2881603371350459916[5] = 0;
   out_2881603371350459916[6] = 1;
   out_2881603371350459916[7] = 0;
   out_2881603371350459916[8] = 0;
}
void h_24(double *state, double *unused, double *out_3436917844040241541) {
   out_3436917844040241541[0] = state[4];
   out_3436917844040241541[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6439528970152783052) {
   out_6439528970152783052[0] = 0;
   out_6439528970152783052[1] = 0;
   out_6439528970152783052[2] = 0;
   out_6439528970152783052[3] = 0;
   out_6439528970152783052[4] = 1;
   out_6439528970152783052[5] = 0;
   out_6439528970152783052[6] = 0;
   out_6439528970152783052[7] = 0;
   out_6439528970152783052[8] = 0;
   out_6439528970152783052[9] = 0;
   out_6439528970152783052[10] = 0;
   out_6439528970152783052[11] = 0;
   out_6439528970152783052[12] = 0;
   out_6439528970152783052[13] = 0;
   out_6439528970152783052[14] = 1;
   out_6439528970152783052[15] = 0;
   out_6439528970152783052[16] = 0;
   out_6439528970152783052[17] = 0;
}
void h_30(double *state, double *unused, double *out_1842384556916435340) {
   out_1842384556916435340[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5399936329857708543) {
   out_5399936329857708543[0] = 0;
   out_5399936329857708543[1] = 0;
   out_5399936329857708543[2] = 0;
   out_5399936329857708543[3] = 0;
   out_5399936329857708543[4] = 1;
   out_5399936329857708543[5] = 0;
   out_5399936329857708543[6] = 0;
   out_5399936329857708543[7] = 0;
   out_5399936329857708543[8] = 0;
}
void h_26(double *state, double *unused, double *out_5258657846565505356) {
   out_5258657846565505356[0] = state[7];
}
void H_26(double *state, double *unused, double *out_859899947523596308) {
   out_859899947523596308[0] = 0;
   out_859899947523596308[1] = 0;
   out_859899947523596308[2] = 0;
   out_859899947523596308[3] = 0;
   out_859899947523596308[4] = 0;
   out_859899947523596308[5] = 0;
   out_859899947523596308[6] = 0;
   out_859899947523596308[7] = 1;
   out_859899947523596308[8] = 0;
}
void h_27(double *state, double *unused, double *out_2962808772413728747) {
   out_2962808772413728747[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3820856270577573193) {
   out_3820856270577573193[0] = 0;
   out_3820856270577573193[1] = 0;
   out_3820856270577573193[2] = 0;
   out_3820856270577573193[3] = 1;
   out_3820856270577573193[4] = 0;
   out_3820856270577573193[5] = 0;
   out_3820856270577573193[6] = 0;
   out_3820856270577573193[7] = 0;
   out_3820856270577573193[8] = 0;
}
void h_29(double *state, double *unused, double *out_4556656960497413806) {
   out_4556656960497413806[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1135861614462756098) {
   out_1135861614462756098[0] = 0;
   out_1135861614462756098[1] = 1;
   out_1135861614462756098[2] = 0;
   out_1135861614462756098[3] = 0;
   out_1135861614462756098[4] = 0;
   out_1135861614462756098[5] = 0;
   out_1135861614462756098[6] = 0;
   out_1135861614462756098[7] = 0;
   out_1135861614462756098[8] = 0;
}
void h_28(double *state, double *unused, double *out_5775360875683582130) {
   out_5775360875683582130[0] = state[0];
}
void H_28(double *state, double *unused, double *out_827768657102570153) {
   out_827768657102570153[0] = 1;
   out_827768657102570153[1] = 0;
   out_827768657102570153[2] = 0;
   out_827768657102570153[3] = 0;
   out_827768657102570153[4] = 0;
   out_827768657102570153[5] = 0;
   out_827768657102570153[6] = 0;
   out_827768657102570153[7] = 0;
   out_827768657102570153[8] = 0;
}
void h_31(double *state, double *unused, double *out_1193646813724223456) {
   out_1193646813724223456[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1486108049756947784) {
   out_1486108049756947784[0] = 0;
   out_1486108049756947784[1] = 0;
   out_1486108049756947784[2] = 0;
   out_1486108049756947784[3] = 0;
   out_1486108049756947784[4] = 0;
   out_1486108049756947784[5] = 0;
   out_1486108049756947784[6] = 0;
   out_1486108049756947784[7] = 0;
   out_1486108049756947784[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2128698594186070786) {
  err_fun(nom_x, delta_x, out_2128698594186070786);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4181753220305909536) {
  inv_err_fun(nom_x, true_x, out_4181753220305909536);
}
void car_H_mod_fun(double *state, double *out_6400451532049496356) {
  H_mod_fun(state, out_6400451532049496356);
}
void car_f_fun(double *state, double dt, double *out_4626213852693615128) {
  f_fun(state,  dt, out_4626213852693615128);
}
void car_F_fun(double *state, double dt, double *out_4088273288682318823) {
  F_fun(state,  dt, out_4088273288682318823);
}
void car_h_25(double *state, double *unused, double *out_1795064971337587010) {
  h_25(state, unused, out_1795064971337587010);
}
void car_H_25(double *state, double *unused, double *out_2881603371350459916) {
  H_25(state, unused, out_2881603371350459916);
}
void car_h_24(double *state, double *unused, double *out_3436917844040241541) {
  h_24(state, unused, out_3436917844040241541);
}
void car_H_24(double *state, double *unused, double *out_6439528970152783052) {
  H_24(state, unused, out_6439528970152783052);
}
void car_h_30(double *state, double *unused, double *out_1842384556916435340) {
  h_30(state, unused, out_1842384556916435340);
}
void car_H_30(double *state, double *unused, double *out_5399936329857708543) {
  H_30(state, unused, out_5399936329857708543);
}
void car_h_26(double *state, double *unused, double *out_5258657846565505356) {
  h_26(state, unused, out_5258657846565505356);
}
void car_H_26(double *state, double *unused, double *out_859899947523596308) {
  H_26(state, unused, out_859899947523596308);
}
void car_h_27(double *state, double *unused, double *out_2962808772413728747) {
  h_27(state, unused, out_2962808772413728747);
}
void car_H_27(double *state, double *unused, double *out_3820856270577573193) {
  H_27(state, unused, out_3820856270577573193);
}
void car_h_29(double *state, double *unused, double *out_4556656960497413806) {
  h_29(state, unused, out_4556656960497413806);
}
void car_H_29(double *state, double *unused, double *out_1135861614462756098) {
  H_29(state, unused, out_1135861614462756098);
}
void car_h_28(double *state, double *unused, double *out_5775360875683582130) {
  h_28(state, unused, out_5775360875683582130);
}
void car_H_28(double *state, double *unused, double *out_827768657102570153) {
  H_28(state, unused, out_827768657102570153);
}
void car_h_31(double *state, double *unused, double *out_1193646813724223456) {
  h_31(state, unused, out_1193646813724223456);
}
void car_H_31(double *state, double *unused, double *out_1486108049756947784) {
  H_31(state, unused, out_1486108049756947784);
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
