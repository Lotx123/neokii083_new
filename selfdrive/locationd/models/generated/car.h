/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5431386215939578595);
void inv_err_fun(double *nom_x, double *true_x, double *out_5340998221229236371);
void H_mod_fun(double *state, double *out_2564987129028209696);
void f_fun(double *state, double dt, double *out_5775547959059003624);
void F_fun(double *state, double dt, double *out_8799050675779776787);
void h_25(double *state, double *unused, double *out_6799911384264016311);
void H_25(double *state, double *unused, double *out_649369991836513177);
void h_24(double *state, double *unused, double *out_7862081742138175572);
void H_24(double *state, double *unused, double *out_6927229677840864079);
void h_30(double *state, double *unused, double *out_5349279991635133144);
void H_30(double *state, double *unused, double *out_8183475170923855159);
void h_26(double *state, double *unused, double *out_6841639360115371218);
void H_26(double *state, double *unused, double *out_5055223163102696096);
void h_27(double *state, double *unused, double *out_410998101990205200);
void H_27(double *state, double *unused, double *out_6895893183087229847);
void h_29(double *state, double *unused, double *out_8095790814413572827);
void H_29(double *state, double *unused, double *out_7041038383279671585);
void h_28(double *state, double *unused, double *out_2070050746269667637);
void H_28(double *state, double *unused, double *out_90602800395630787);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
