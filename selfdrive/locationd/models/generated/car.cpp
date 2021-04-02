
extern "C"{

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

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5431386215939578595) {
   out_5431386215939578595[0] = delta_x[0] + nom_x[0];
   out_5431386215939578595[1] = delta_x[1] + nom_x[1];
   out_5431386215939578595[2] = delta_x[2] + nom_x[2];
   out_5431386215939578595[3] = delta_x[3] + nom_x[3];
   out_5431386215939578595[4] = delta_x[4] + nom_x[4];
   out_5431386215939578595[5] = delta_x[5] + nom_x[5];
   out_5431386215939578595[6] = delta_x[6] + nom_x[6];
   out_5431386215939578595[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5340998221229236371) {
   out_5340998221229236371[0] = -nom_x[0] + true_x[0];
   out_5340998221229236371[1] = -nom_x[1] + true_x[1];
   out_5340998221229236371[2] = -nom_x[2] + true_x[2];
   out_5340998221229236371[3] = -nom_x[3] + true_x[3];
   out_5340998221229236371[4] = -nom_x[4] + true_x[4];
   out_5340998221229236371[5] = -nom_x[5] + true_x[5];
   out_5340998221229236371[6] = -nom_x[6] + true_x[6];
   out_5340998221229236371[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_2564987129028209696) {
   out_2564987129028209696[0] = 1.0;
   out_2564987129028209696[1] = 0.0;
   out_2564987129028209696[2] = 0.0;
   out_2564987129028209696[3] = 0.0;
   out_2564987129028209696[4] = 0.0;
   out_2564987129028209696[5] = 0.0;
   out_2564987129028209696[6] = 0.0;
   out_2564987129028209696[7] = 0.0;
   out_2564987129028209696[8] = 0.0;
   out_2564987129028209696[9] = 1.0;
   out_2564987129028209696[10] = 0.0;
   out_2564987129028209696[11] = 0.0;
   out_2564987129028209696[12] = 0.0;
   out_2564987129028209696[13] = 0.0;
   out_2564987129028209696[14] = 0.0;
   out_2564987129028209696[15] = 0.0;
   out_2564987129028209696[16] = 0.0;
   out_2564987129028209696[17] = 0.0;
   out_2564987129028209696[18] = 1.0;
   out_2564987129028209696[19] = 0.0;
   out_2564987129028209696[20] = 0.0;
   out_2564987129028209696[21] = 0.0;
   out_2564987129028209696[22] = 0.0;
   out_2564987129028209696[23] = 0.0;
   out_2564987129028209696[24] = 0.0;
   out_2564987129028209696[25] = 0.0;
   out_2564987129028209696[26] = 0.0;
   out_2564987129028209696[27] = 1.0;
   out_2564987129028209696[28] = 0.0;
   out_2564987129028209696[29] = 0.0;
   out_2564987129028209696[30] = 0.0;
   out_2564987129028209696[31] = 0.0;
   out_2564987129028209696[32] = 0.0;
   out_2564987129028209696[33] = 0.0;
   out_2564987129028209696[34] = 0.0;
   out_2564987129028209696[35] = 0.0;
   out_2564987129028209696[36] = 1.0;
   out_2564987129028209696[37] = 0.0;
   out_2564987129028209696[38] = 0.0;
   out_2564987129028209696[39] = 0.0;
   out_2564987129028209696[40] = 0.0;
   out_2564987129028209696[41] = 0.0;
   out_2564987129028209696[42] = 0.0;
   out_2564987129028209696[43] = 0.0;
   out_2564987129028209696[44] = 0.0;
   out_2564987129028209696[45] = 1.0;
   out_2564987129028209696[46] = 0.0;
   out_2564987129028209696[47] = 0.0;
   out_2564987129028209696[48] = 0.0;
   out_2564987129028209696[49] = 0.0;
   out_2564987129028209696[50] = 0.0;
   out_2564987129028209696[51] = 0.0;
   out_2564987129028209696[52] = 0.0;
   out_2564987129028209696[53] = 0.0;
   out_2564987129028209696[54] = 1.0;
   out_2564987129028209696[55] = 0.0;
   out_2564987129028209696[56] = 0.0;
   out_2564987129028209696[57] = 0.0;
   out_2564987129028209696[58] = 0.0;
   out_2564987129028209696[59] = 0.0;
   out_2564987129028209696[60] = 0.0;
   out_2564987129028209696[61] = 0.0;
   out_2564987129028209696[62] = 0.0;
   out_2564987129028209696[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5775547959059003624) {
   out_5775547959059003624[0] = state[0];
   out_5775547959059003624[1] = state[1];
   out_5775547959059003624[2] = state[2];
   out_5775547959059003624[3] = state[3];
   out_5775547959059003624[4] = state[4];
   out_5775547959059003624[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5775547959059003624[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5775547959059003624[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8799050675779776787) {
   out_8799050675779776787[0] = 1;
   out_8799050675779776787[1] = 0;
   out_8799050675779776787[2] = 0;
   out_8799050675779776787[3] = 0;
   out_8799050675779776787[4] = 0;
   out_8799050675779776787[5] = 0;
   out_8799050675779776787[6] = 0;
   out_8799050675779776787[7] = 0;
   out_8799050675779776787[8] = 0;
   out_8799050675779776787[9] = 1;
   out_8799050675779776787[10] = 0;
   out_8799050675779776787[11] = 0;
   out_8799050675779776787[12] = 0;
   out_8799050675779776787[13] = 0;
   out_8799050675779776787[14] = 0;
   out_8799050675779776787[15] = 0;
   out_8799050675779776787[16] = 0;
   out_8799050675779776787[17] = 0;
   out_8799050675779776787[18] = 1;
   out_8799050675779776787[19] = 0;
   out_8799050675779776787[20] = 0;
   out_8799050675779776787[21] = 0;
   out_8799050675779776787[22] = 0;
   out_8799050675779776787[23] = 0;
   out_8799050675779776787[24] = 0;
   out_8799050675779776787[25] = 0;
   out_8799050675779776787[26] = 0;
   out_8799050675779776787[27] = 1;
   out_8799050675779776787[28] = 0;
   out_8799050675779776787[29] = 0;
   out_8799050675779776787[30] = 0;
   out_8799050675779776787[31] = 0;
   out_8799050675779776787[32] = 0;
   out_8799050675779776787[33] = 0;
   out_8799050675779776787[34] = 0;
   out_8799050675779776787[35] = 0;
   out_8799050675779776787[36] = 1;
   out_8799050675779776787[37] = 0;
   out_8799050675779776787[38] = 0;
   out_8799050675779776787[39] = 0;
   out_8799050675779776787[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8799050675779776787[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8799050675779776787[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8799050675779776787[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8799050675779776787[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8799050675779776787[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8799050675779776787[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8799050675779776787[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8799050675779776787[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8799050675779776787[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8799050675779776787[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8799050675779776787[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8799050675779776787[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8799050675779776787[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8799050675779776787[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8799050675779776787[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8799050675779776787[56] = 0;
   out_8799050675779776787[57] = 0;
   out_8799050675779776787[58] = 0;
   out_8799050675779776787[59] = 0;
   out_8799050675779776787[60] = 0;
   out_8799050675779776787[61] = 0;
   out_8799050675779776787[62] = 0;
   out_8799050675779776787[63] = 1;
}
void h_25(double *state, double *unused, double *out_6799911384264016311) {
   out_6799911384264016311[0] = state[6];
}
void H_25(double *state, double *unused, double *out_649369991836513177) {
   out_649369991836513177[0] = 0;
   out_649369991836513177[1] = 0;
   out_649369991836513177[2] = 0;
   out_649369991836513177[3] = 0;
   out_649369991836513177[4] = 0;
   out_649369991836513177[5] = 0;
   out_649369991836513177[6] = 1;
   out_649369991836513177[7] = 0;
}
void h_24(double *state, double *unused, double *out_7862081742138175572) {
   out_7862081742138175572[0] = state[4];
   out_7862081742138175572[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6927229677840864079) {
   out_6927229677840864079[0] = 0;
   out_6927229677840864079[1] = 0;
   out_6927229677840864079[2] = 0;
   out_6927229677840864079[3] = 0;
   out_6927229677840864079[4] = 1;
   out_6927229677840864079[5] = 0;
   out_6927229677840864079[6] = 0;
   out_6927229677840864079[7] = 0;
   out_6927229677840864079[8] = 0;
   out_6927229677840864079[9] = 0;
   out_6927229677840864079[10] = 0;
   out_6927229677840864079[11] = 0;
   out_6927229677840864079[12] = 0;
   out_6927229677840864079[13] = 1;
   out_6927229677840864079[14] = 0;
   out_6927229677840864079[15] = 0;
}
void h_30(double *state, double *unused, double *out_5349279991635133144) {
   out_5349279991635133144[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8183475170923855159) {
   out_8183475170923855159[0] = 0;
   out_8183475170923855159[1] = 0;
   out_8183475170923855159[2] = 0;
   out_8183475170923855159[3] = 0;
   out_8183475170923855159[4] = 1;
   out_8183475170923855159[5] = 0;
   out_8183475170923855159[6] = 0;
   out_8183475170923855159[7] = 0;
}
void h_26(double *state, double *unused, double *out_6841639360115371218) {
   out_6841639360115371218[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5055223163102696096) {
   out_5055223163102696096[0] = 0;
   out_5055223163102696096[1] = 0;
   out_5055223163102696096[2] = 0;
   out_5055223163102696096[3] = 0;
   out_5055223163102696096[4] = 0;
   out_5055223163102696096[5] = 0;
   out_5055223163102696096[6] = 0;
   out_5055223163102696096[7] = 1;
}
void h_27(double *state, double *unused, double *out_410998101990205200) {
   out_410998101990205200[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6895893183087229847) {
   out_6895893183087229847[0] = 0;
   out_6895893183087229847[1] = 0;
   out_6895893183087229847[2] = 0;
   out_6895893183087229847[3] = 1;
   out_6895893183087229847[4] = 0;
   out_6895893183087229847[5] = 0;
   out_6895893183087229847[6] = 0;
   out_6895893183087229847[7] = 0;
}
void h_29(double *state, double *unused, double *out_8095790814413572827) {
   out_8095790814413572827[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7041038383279671585) {
   out_7041038383279671585[0] = 0;
   out_7041038383279671585[1] = 1;
   out_7041038383279671585[2] = 0;
   out_7041038383279671585[3] = 0;
   out_7041038383279671585[4] = 0;
   out_7041038383279671585[5] = 0;
   out_7041038383279671585[6] = 0;
   out_7041038383279671585[7] = 0;
}
void h_28(double *state, double *unused, double *out_2070050746269667637) {
   out_2070050746269667637[0] = state[5];
   out_2070050746269667637[1] = state[6];
}
void H_28(double *state, double *unused, double *out_90602800395630787) {
   out_90602800395630787[0] = 0;
   out_90602800395630787[1] = 0;
   out_90602800395630787[2] = 0;
   out_90602800395630787[3] = 0;
   out_90602800395630787[4] = 0;
   out_90602800395630787[5] = 1;
   out_90602800395630787[6] = 0;
   out_90602800395630787[7] = 0;
   out_90602800395630787[8] = 0;
   out_90602800395630787[9] = 0;
   out_90602800395630787[10] = 0;
   out_90602800395630787[11] = 0;
   out_90602800395630787[12] = 0;
   out_90602800395630787[13] = 0;
   out_90602800395630787[14] = 1;
   out_90602800395630787[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
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



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
