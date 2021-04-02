/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8431838611904138524);
void inv_err_fun(double *nom_x, double *true_x, double *out_5168483825837128788);
void H_mod_fun(double *state, double *out_6292912079102671410);
void f_fun(double *state, double dt, double *out_2898850656836748995);
void F_fun(double *state, double dt, double *out_2192292858371402905);
void h_3(double *state, double *unused, double *out_1063937345248663084);
void H_3(double *state, double *unused, double *out_8614481489080304765);
void h_4(double *state, double *unused, double *out_8722033022806688192);
void H_4(double *state, double *unused, double *out_2340003312655773576);
void h_9(double *state, double *unused, double *out_5301673063030202187);
void H_9(double *state, double *unused, double *out_2339848545563697117);
void h_10(double *state, double *unused, double *out_7870446520771126277);
void H_10(double *state, double *unused, double *out_4330286694514542411);
void h_12(double *state, double *unused, double *out_84988027843208930);
void H_12(double *state, double *unused, double *out_4446558824394983145);
void h_31(double *state, double *unused, double *out_5748600440082468485);
void H_31(double *state, double *unused, double *out_6169826432235931452);
void h_32(double *state, double *unused, double *out_6706269596097888632);
void H_32(double *state, double *unused, double *out_4219324710588179805);
void h_13(double *state, double *unused, double *out_8521959366806518552);
void H_13(double *state, double *unused, double *out_2714277075622101920);
void h_14(double *state, double *unused, double *out_5301673063030202187);
void H_14(double *state, double *unused, double *out_2339848545563697117);
void h_19(double *state, double *unused, double *out_2873651588677245017);
void H_19(double *state, double *unused, double *out_2092831305538734064);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);