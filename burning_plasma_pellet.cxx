/*******************************************************************************
 *  * Burning_plasma case
 *  *
 *  *
 *  * Solves equations for deuterium and tritium nuclear reaction
 *  * with fluid method
 *  * deuterium ion  density                          Ni
 *  * helium ion density                              Ni_he
 *  * tritium ion density                             Ni_t
 *  * electron density                                Ne=Ni+Ni_t+2*Ni_he
 *  * alpha particle density                          N_alpha
 *  * electron and deuterium ion temperatures         Te, Ti
 *  * tritium ion temperatures                        Ti_t
 *  * helium ion temperatures                         Ti_he
 *  *
 *  * fueling background ion profile
 *  *                                                 SnD      Deuterium atom
 *  *                                                 SnT      Tritium atom
 *  *                                                 S_pellet pellet profile
 *  * RF heating                                      Prfe
 *  * NBI heating                                     Pauxi
 *  *                                                 Pauxe
 *  * Radiation                                       Prad
 *  * The fusion formula is vaild when Ti<25keV
 *  * 1D transport problem for burning loss and fueling
 *  * Author: Lin Weisheng    Email:2018lws@pku.edu.cn
 ********************************************************************************/


#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>

#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <interpolation.hxx>

#include <cmath>
#include <math.h>

// Vectors
Vector2D b0xcv; // Curvature term

//2D Evolving fields from Grid file
Field2D Te_grid, Ti_grid, Ni_grid;
Field2D Te_exp, Ti_exp, Ni_exp;

//2D dx,dy
Field2D dx,dy;

//2D Psi from Grid file
Field2D psixy;

//ShiftAngle
Field2D q_profile;

// 3D evolving fields
Field3D Te, Ne;             // density,temperature,parallel velocity of Deuterium ion
Field3D Ni, Ti;             // density,temperature,parallel velocity of Deuterium ion
Field3D Ni_t, Ti_t;      // density, temperature and Velocity of Tritium ion
Field3D Ni_he, Ti_he;      // density, temperature and Velocity of Helium ion
Field3D temp_Ni, temp_Ti, temp_Te;         // temporary variables for minivalue protection
Field3D temp_Ni_t, temp_Ti_t, temp_Ni_he, temp_Ti_he;         // temporary variables for minivalue protection
Field3D temp_Ne;
Field3D N_alpha;  //density of fast alpha particle


//Test
//Field3D Test, Test1, Test2, Test3, Test4;
// parameters for initial sources of atom D and T
Field3D SnD, SnT;
BoutReal SnD_amp, SnD_x0, SnD_xwidth, SnT_amp, SnT_x0, SnT_xwidth;
bool fueling_atomD, fueling_atomT;
Field3D P_fuel_loss;

//parameters for Pauxe,Pauxi,Pauxhe,Pauxi_t   NBI
Field3D Pauxe,Pauxi,Pauxhe,Pauxi_t;
BoutReal Pauxe_amp, Pauxe_x0, Pauxe_xwidth;
bool auxiliary_heating_e;
BoutReal Pauxi_amp, Pauxi_x0, Pauxi_xwidth;
bool auxiliary_heating_i;
BoutReal Pauxi_t_amp, Pauxi_t_x0, Pauxi_t_xwidth;
bool auxiliary_heating_i_t;
BoutReal Pauxhe_amp, Pauxhe_x0, Pauxhe_xwidth;
bool auxiliary_heating_he;

//parameters for Prfe  ECCD
Field3D Prfe;
bool auxiliary_rfheating_e;
BoutReal Prfe_amp, Prfe_x0, Prfe_xwidth;

//parameters for heating profile grid
Field3D Prfe_grid,Pauxi_grid,Pauxe_grid,Prfi_grid,Prad_grid;

//Diffusion coeffient factor
BoutReal Diff_psi_factor;
BoutReal kappa_psi_factor;//This value is kappa/D*0.1

//pinch factor
BoutReal pinch_factor;

//Scaling thermal conductivites
Field3D kappa_psi;


//Fusion reaction rate between D and T
Field3D sigma_DT_V_th;
Field3D temp_sigma_DT_V_th;

//radiation loss
Field3D P_brem,P_cyclo;
//Field3D P_line_loss;
//Field3D P_rec;
Field3D P_radiation;

//alpha power to D,T,He and electron
Field3D PaD,PaT,Pahe,Pae;

// Thermal Ion and Electron Speed
Field3D V_th_i, V_th_e;

////Collision time between same particles
//Field3D tau_e, tau_i, tau_i_he, tau_i_t;

//Slowing time between alpha and electron
Field3D tau_slowing_se;

//alpha particle slowing down time
Field3D tau_slowing_down;

//alpha particle pressure
Field3D alpha_pressure;

////Collision rate between ions
Field3D col_eD, col_eT, col_ehe;
Field3D col_De, col_DT, col_Dhe;
Field3D col_Te, col_TD, col_The;
Field3D col_hee, col_heD, col_heT, col_hehe;
Field3D lambda_De, lambda_Te, lambda_hee;
Field3D lambda_DT, lambda_The, lambda_Dhe;
Field3D lambda_hehe;

//Classical diffusion coefficient of  plasma
//Field3D Diff_eD, Diff_eT, Diff_ehe;


//Diffusion coefficient depended on psi
Field3D Diff_psi;
Field3D sqrt_psi;
Field3D psi_norm;
Vector3D V_pinch; // pinch velocity
BoutReal Ip_98, B_98, P_98, n_98, M_98, R_98, epsi_98, kappa_98,P_cd;
BoutReal tau_e98;
BoutReal a_radius;
Field3D F_r_pinch;
Field3D F_r_pinch_2;
BoutReal rho_coef,p_coef,A_coef;
BoutReal n_average;
BoutReal S_fu_total;
BoutReal P_fuel_loss_total;
BoutReal P_radiation_total;
BoutReal P_cd_total;
BoutReal P_fusion_total;
BoutReal Prfe_grid_total;
BoutReal Pauxi_grid_total;
BoutReal Pauxe_grid_total;

//update kappa_factor
Field3D Wth;
BoutReal tau_e_real;
BoutReal Wth_total;//total thermal energy

//Alpha particle diffusion coefficent and pinch velocity
Field3D Diff_alpha;
Field3D V_pinch_alpha;
Field3D DDX_Te;

//Slowing down caculate parameter
BoutReal Z_He,Z_D,Z_T,Z_b;
BoutReal A_b,A_D,A_T,A_He;
Field3D Z_average;
Field3D W_c;
Field3D ln_lambda;
BoutReal W_b0;
BoutReal E_alpha;

//fraction to ion and electron
Field3D fi,fe;

//para diffusion coefficent
Field3D Diff_psi_para;

//para thermal conductivites
Field3D kappa_psi_para;

//Diffusion coefficiency of helium
//Field3D Diff_he;


// Fusion source Terms
Field3D S_fu;      //Deuterium and Tritium fusion sources term
bool burning;

//pellet parameter
Field3D pellet_rate;
Field3D pellet_rp_profile;
Field3D pellet_dN_rate;
Field3D pellet_dN_dr;

//ELM section   Note:Need the caculation result from EPED
bool trigger_ELM;
bool ELM_happen;
BoutReal ped_position;
BoutReal ped_width;
BoutReal critical_pressure;
Field3D pressure;
int ped_index; //should find out before running
Field3D Diff_ELM;
BoutReal Diff_ELM_amp;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe, Zxy;
Field2D I;

// Max grid number in X Y Z
int NX, NY, NZ;


// Initial profile parameters
BoutReal Te_core, Te_edge, Ti_core, Ti_edge, Ni_core, Ni_edge;
BoutReal Initfile_x0, Initfile_w_ped;
BoutReal density_unit;                   // Number density [m^-3]

// parameters
BoutReal Te_x, Ti_x, Tm_x, Ni_x, Vi_x, bmag, rho_s, AA, ZZ;
BoutReal minimum_val;
BoutReal Mm, Mi, Mn, Mi_t, Mi_he, Mn_t;
BoutReal W_ionz, W_diss, W_bind, W_rec;
BoutReal Lbar, tbar;

BoutReal wci;

//pellet parameter
Field3D S_pellet;
int mid;
BoutReal vp,Mi_p,np;
BoutReal rp0;
BoutReal S_pellet_total;
int timestep;
BoutReal deposit_depth;
BoutReal pellet_life_time;
BoutReal pellet_frequency;
int pellet_timestep;
Field3D T_pellet;
Field3D r_cloud;
Field3D L_c;
Field3D pellet_displacement;
Field3D add_pellet_displacement;
int t_output;
Field3D pellet_drift_psi;
Field3D R_new;

//radial global x
BoutReal x_rela;

//constant gradient at xin with auto unit
BoutReal dNidx_xin_au, dTidx_xin_au, dTedx_xin_au;
BoutReal psi_xout_y0, psi_axis, psi_bndry;


//logical paramter
bool load_grid_profiles, initial_profile_exp, initial_profile_linear, initial_profile_square;
bool load_experiment_profiles;
bool initial_profile_Hmode;
bool noshear;

const Field3D field_larger(const Field3D &f, const BoutReal limit);
const BoutReal average(const Field2D &f,const Field3D &k);
const BoutReal total(const Field2D &f,const Field3D &k);
const Field3D integral_fe(const BoutReal down, const BoutReal up, int n);
const Field3D integral_fi(const BoutReal down, const BoutReal up, int n);
//const Field3D Pellet_rp_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp,const int mid,\
//                             const Field3D &Ne,const Field3D &Te, const Field3D &R);
const Field3D Pellet_rp_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp,const int mid);
const Field3D Pellet_rp_profile(const Field3D &rate,const BoutReal rp0);
const Field3D Pellet_dN_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp);
const Field3D Pellet_dN_dr(const BoutReal vp);
const Field3D Pellet_dN_profile();
const BoutReal find_rp_deposit_depth();
const Field3D Pellet_displacement();
const Field3D Add_pellet_displacement();
const Field3D fit_dN_profile();
const Field3D Interp_R_new();
//const void test();
const bool test_ELM_happen();
//const void revise_ELM_coef();

const BoutReal PI = 3.14159265;
//const BoutReal MU0 = 4.0e-7 * PI;
//const BoutReal Mi = 2.0*1.6726e-27; // Ion mass
//const BoutReal Me = 0.5 / 1836.2;        // in unit of Mass_deuterium
//const BoutReal KB = 1.38065e-23;     // Boltamann constant
//const BoutReal ee = 1.602e-19;       // ln(Lambda)
//const BoutReal eV_K = 11605.0;         // 1eV = 11605K

//const BoutReal C_fe_sheat = 7.;     // coefficient of parallel heat transmission in sheat BC
//const BoutReal C_fi_sheat = 2.5;
const BoutReal mass_e = 9.1e-28;    //in csg unit
const BoutReal mass_p = 1.672e-24;

int physics_init(bool restarting) {

    t_output = 1;


    /////////////// LOAD DATA FROM GRID FILE //////////////

// Load 2D profiles
    mesh->get(Te_grid, "Te0");    // eV
    mesh->get(Ti_grid, "Ti0");    // eV
    mesh->get(Ni_grid, "Ni0");    // hl2a 10^21/m^3

    mesh->get(Te_exp, "Te_exp");    // eV
    mesh->get(Ti_exp, "Ti_exp");    // eV
    mesh->get(Ni_exp, "Ni_exp");    // hl2a 10^19/m^3

    // Load Psi
    mesh->get(psixy, "psixy");         // unit m^2 T
    mesh->get(psi_axis, "psi_axis");
    mesh->get(psi_bndry, "psi_bndry");
    // Load curvature term
    b0xcv.covariant = false; // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

    // Load metrics
    GRID_LOAD(Rxy);         // Major radius [m]
    GRID_LOAD(Zxy);         // Height [m]
    GRID_LOAD2(Bpxy, Btxy); // Poloidal, Toroidal B field [T]
    GRID_LOAD(hthe);        // Poloidal arc length [m / radian]
    //mesh->get(mesh->dx,   "dpsi");
    //mesh->get(mesh->dx, "dx");  //1D test only
    mesh->get(I, "sinty");// m^-2 T^-1
    mesh->get(NX, "nx");
    mesh->get(NY, "ny");
    mesh->get(dx,"dx");
    mesh->get(dy,"dy");

    if (mesh->get(bmag, "bmag"))
        bmag = 1.0;
    if (mesh->get(Lbar, "rmag"))
        Lbar = 1.0;

    // Load normalisation values
    //GRID_LOAD(Te_x);
    //GRID_LOAD(Ti_x);
    //GRID_LOAD(Ni_x);
    //GRID_LOAD(bmag);

    //load heating profile
    mesh->get(Prfe_grid,"Prfe");
    mesh->get(Pauxi_grid,"Pnbii");
    mesh->get(Pauxe_grid,"Pnbie");
    mesh->get(Prad_grid,"Prad");

    //q_profile
    mesh->get(q_profile,"q");
    /////////////// READ OPTIONS //////////////////////////

    // Read some parameters
    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("nucl-fusion");


    OPTION(options, minimum_val, 1.e-10);    // minimum value limit for densities Ns

    OPTION(options, NZ, 1);        // maximum grid number in Z

    OPTION(options, Te_x, 10.0);    // Read in eV
    OPTION(options, Ti_x, 10.0);    // eV
    OPTION(options, Tm_x, 0.0258);    // eV
    OPTION(options, Ni_x, 1.0);     // in 1.e^19 m^-3
    OPTION(options, density_unit, 1.0e19); // Number density [m^-3]


    OPTION(options, noshear, true);

    OPTION(options, load_grid_profiles, false);     // priority II
    OPTION(options, load_experiment_profiles, false);     // priority III
    OPTION(options, initial_profile_exp, true);     // priority III
    OPTION(options, initial_profile_linear, false); // priority III
    OPTION(options, initial_profile_Hmode, false); // priority III
    OPTION(options, initial_profile_square, false); // priority III

    OPTION(options, AA, 2.0);
    OPTION(options, ZZ, 1.0);

    OPTION(options, Mi, 1.0);        // Read in Mi
    OPTION(options, Mi_t, 1.5);        //Read in Mi
    OPTION(options, Mi_he, 2);         //Read in Mi
    OPTION(options, Mn, 1.0);        // Read in Mi
    OPTION(options, Mn_t, 1.5);        // Read in Mi
    OPTION(options, Mm, 2.0);        // Read in Mi
    OPTION(options, W_ionz, 20.0);    // Read in eV
    OPTION(options, W_diss, 4.5);     // in eV
    OPTION(options, W_rec, 4.5);      // in eV


    OPTION(options, W_bind, 0.5);     // in eV

    OPTION(options, dNidx_xin_au, -65.);    // a.u.
    OPTION(options, dTidx_xin_au, -3500.);  // a.u.
    OPTION(options, dTedx_xin_au, -3500.);  // a.u.
    OPTION(options, psi_xout_y0, 0.229543);  //m^2 T
    OPTION(options, Te_core, 1000.);  //eV
    OPTION(options, Te_edge, 10.);    //eV
    OPTION(options, Ti_core, 1000.);  //eV
    OPTION(options, Ti_edge, 10.);    //eV
    OPTION(options, Ni_core, 2.);     //in 1.e^19 m^-3
    OPTION(options, Ni_edge, 0.1);    //in 1.e^19 m^-3
    OPTION(options, Initfile_x0, 0.01);    //in a.u.
    OPTION(options, Initfile_w_ped, 0.2);    //in a.u.
    // Initial Sources of Deuterium, Tritium

    OPTION(options, fueling_atomD, false);    //in a.u., turn on atom Deuterium fueling
    OPTION(options, SnD_amp, 0.1);    //in a.u., amplitude
    OPTION(options, SnD_x0, 0.4);    //in a.u., central position
    OPTION(options, SnD_xwidth, 0.1);    //in a.u., width in x direction

    OPTION(options, fueling_atomT, false);    //in a.u., turn on atom Tritium fueling
    OPTION(options, SnT_amp, 0.1);    //in a.u., amplitude
    OPTION(options, SnT_x0, 0.4);    //in a.u., central position
    OPTION(options, SnT_xwidth, 0.2);    //in a.u., width in x direction


    //parameter for ECCD
    OPTION(options, auxiliary_rfheating_e, false);
    OPTION(options, Prfe_amp, 0.1);    //in a.u., amplitude
    OPTION(options, Prfe_x0, 0.01);    //in a.u., central position
    OPTION(options, Prfe_xwidth, 0.5);    //in a.u., width in x direction

    //parameter for auxiliary heating electron
    OPTION(options, auxiliary_heating_e, false);
    OPTION(options, Pauxe_amp, 0.1);    //in a.u., amplitude
    OPTION(options, Pauxe_x0, 0.01);    //in a.u., central position
    OPTION(options, Pauxe_xwidth, 0.5);    //in a.u., width in x direction

    //parameter for auxiliary heating i
    OPTION(options, auxiliary_heating_i, false);
    OPTION(options, Pauxi_amp, 0.1);    //in a.u., amplitude
    OPTION(options, Pauxi_x0, 0.01);    //in a.u., central position
    OPTION(options, Pauxi_xwidth, 0.5);    //in a.u., width in x direction

    //parameter for auxiliary heating i_t
    OPTION(options, auxiliary_heating_i_t, false);
    OPTION(options, Pauxi_t_amp, 0.1);    //in a.u., amplitude
    OPTION(options, Pauxi_t_x0, 0.01);    //in a.u., central position
    OPTION(options, Pauxi_t_xwidth, 0.5);    //in a.u., width in x direction

    //parameter for auxiliary heating he
    OPTION(options, auxiliary_heating_he, false);    //in a.u., turn on atom Tritium fueling
    OPTION(options, Pauxhe_amp, 0.1);    //in a.u., amplitude
    OPTION(options, Pauxhe_x0, 0.01);    //in a.u., central position
    OPTION(options, Pauxhe_xwidth, 0.5);    //in a.u., width in x direction

    //parameter for heating profile from grid
//    OPTION(options, Pauxi_grid, 1);
//    OPTION(options, Pauxe_grid, 1);
//    OPTION(options, Prfe,i);

//Diffusion coefficent
    OPTION(options, Diff_psi_factor, 0.1);
    OPTION(options, kappa_psi_factor, 0.8);
    OPTION(options, pinch_factor, 1.);

    //tau_e98 caculate parameter
    OPTION(options, Ip_98, 13.0);//MA
    OPTION(options, B_98, 6.53);//T
    OPTION(options, P_98, 296);//MW
    OPTION(options, P_cd,80);//MW
    OPTION(options, n_98, 11.6);//10^19m^-3
    OPTION(options, M_98, 2.45);
    OPTION(options, R_98, 7.2);
    OPTION(options, epsi_98, 0.31);//inverse of aspect ratio
    OPTION(options, kappa_98, 2.03);
    OPTION(options, tau_e98, 1.);
    OPTION(options, a_radius, 2.22);//m
    OPTION(options, A_coef, 0.15);
    OPTION(options, p_coef, 40);
    OPTION(options, rho_coef, 0.95);
    OPTION(options, Wth_total, 10000);
    OPTION(options, P_radiation_total, 0.1);
    OPTION(options, P_fuel_loss_total, 1);
    OPTION(options, P_cd_total, 1);
    OPTION(options, P_fusion_total, 1);

    //tau_e_real
    OPTION(options, tau_e_real,1.);

    //slowing down caculate parameter
    OPTION(options, Z_He, 2);
    OPTION(options, Z_D, 1);
    OPTION(options, Z_T, 1);
    OPTION(options, Z_b, 2);
    OPTION(options, A_b, 4);
    OPTION(options, A_D, 2);
    OPTION(options, A_T, 2);
    OPTION(options, A_He, 4);
    OPTION(options, W_b0, 3500);//keV
    OPTION(options, E_alpha, 3.5e5);//10eV

    //Pellet parameter
    OPTION(options, vp, 500);//cm/s
    OPTION(options, Mi_p, 2); //D pellet
    OPTION(options, np, 5.96e22);// cm^-3
    OPTION(options, mid, 46);
    OPTION(options, rp0, 0.25);
    OPTION(options, timestep, 1);
    OPTION(options, deposit_depth, 1.0);
    OPTION(options, pellet_life_time, 1.0);
    OPTION(options, pellet_frequency, 16);
    OPTION(options, pellet_timestep, 290);

    //ELM
    OPTION(options, trigger_ELM, true);
    OPTION(options, ELM_happen, false);
    OPTION(options, ped_position, 0.95);
    OPTION(options, ped_width, 0.0456/2);
    OPTION(options, critical_pressure, 80.99);//kPa
    OPTION(options, ped_index, 63);
    OPTION(options, Diff_ELM_amp, 10);

    //burning
    OPTION(options, burning, false);    //in a.u., turn on burning effect

    ////////////// CALCULATE PARAMETERS ///////////////////

    //in order to using formulas of wci and rhos in Gaussian cgs units,
    // bmag in SI unit Tesla multiply by 1.e4 changes to Gaussian unit gauss
    // Te_x in unit eV, no need to transfer

    // !!!Be very careful when using quantities calculated below, units transfer needed

    Ni_x *= 1.0e13;    // in unit cm^-3 now
    bmag *= 1.0e4;     // in unit gauss now

    output.write("Calculating parameters in Gaussian units and then transfer to SI units \n");
    output.write("\tIn Gaussian cgs units:  Ni_x %e cm^-3, Te_x %e eV, bmag %e gauss  \n", Ni_x, Te_x, bmag);
    output.write("\tparameters:  AA=mi/mp %e,  ZZ %e \n", AA, ZZ);

    rho_s = 1.02e2 * sqrt(AA * Te_x) / ZZ / bmag;   // unit cm
    wci = 9.58e3 * ZZ * bmag / AA;
    Vi_x = wci * rho_s;                     // in unit cm/s
    Ni_x /= 1.0e13;                         // back to unit 1.e^19 m^-3 mow
    bmag /= 1.0e4;                           // back to unit tesla now
    Vi_x /= 1.e2;                           // back to unit m/s

    tbar = Lbar / Vi_x;                       // in unit s

    output.write("\tIn SI units:  Ni_x %e e^19 m^-3,  Te_x %e eV, bmag %e tesla \n",
                 Ni_x, Te_x, bmag);

    ///////////// PRINT Z INFORMATION /////////////////////

    BoutReal hthe0;
    if (GRID_LOAD(hthe0) == 0) {
        output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0 / Lbar);
    }

    b0xcv = 0.0;

    if (noshear) {
        mesh->ShiftXderivs = false;
        I = 0.0;
    }

    //////////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    if (mesh->ShiftXderivs) {
        if (mesh->IncIntShear) {
            // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
            mesh->IntShiftTorsion = I;

        } else {
            // Dimits style, using local coordinate system
            I = 0.0;  // I disappears from metric
        }
    }

    ///////////// NORMALISE QUANTITIES ////////////////////

    output.write("\tNormalising to Lbar = %e m, tbar %e s, V_Ti %e m/s \n", Lbar, tbar, Vi_x);

    // Normalise geometry
    Rxy /= Lbar;
    Zxy /= Lbar;
    hthe /= Lbar;
    mesh->dx /= Lbar * Lbar * bmag;
    //mesh->dx /= Lbar;
    I *= Lbar * Lbar * bmag;
    psixy /= Lbar * Lbar * bmag;
    psi_axis /= Lbar * Lbar * bmag;
    psi_bndry /= Lbar * Lbar * bmag;
    psi_xout_y0 /= Lbar * Lbar * bmag;
    // Normalise magnetic field
    b0xcv.x /= bmag;
    b0xcv.y *= Lbar * Lbar;
    b0xcv.z *= Lbar * Lbar;

    Bpxy /= bmag;
    Btxy /= bmag;
    mesh->Bxy /= bmag;

    // Normalise coefficients
    W_ionz /= Te_x;
    W_diss /= Te_x;
    W_bind /= Te_x;
    Tm_x /= Te_x;
    Te_core /= Te_x;
    Te_edge /= Te_x;
    Ti_core /= Te_x;
    Ti_edge /= Te_x;
    Ni_core /= Ni_x;
    Ni_edge /= Ni_x;

//    Te_grid /= Te_x;
//    Ti_grid /= Te_x;
//    Ni_grid = Ni_grid * 100. / Ni_x;  //hl2a grid in unit 10^19/m^3
    Te_grid = Te_grid * 1000/Te_x;               //CFETR grid in unit keV
    Ti_grid = Ti_grid * 1000/Te_x;                //CFETR grid in unit keV
    Ni_grid = Ni_grid / Ni_x;  //CFETR Hybrid3530 grid in unit 10^19/m^3

    Pauxi_grid = Pauxi_grid*tbar/Ti_x/1.602/Ni_x;//CFETR
    Pauxe_grid = Pauxe_grid*tbar/Ti_x/1.602/Ni_x;//CFETR
    Prfe_grid = Prfe_grid*tbar/Ti_x/1.602/Ni_x;//CFETR
    Prad_grid = Prad_grid*tbar/Ti_x/1.602/Ni_x;//CFETR

    Te_exp /= Te_x;
    Ti_exp /= Te_x;
    Ni_exp /= Ni_x;

    /////////////// CALCULATE METRICS /////////////////

    mesh->g11 = (Rxy * Bpxy) ^ 2;
    mesh->g22 = 1.0 / (hthe ^ 2);
    mesh->g33 = (I ^ 2) * mesh->g11 + (mesh->Bxy ^ 2) / mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -I * mesh->g11;
    mesh->g23 = -Btxy / (hthe * Bpxy * Rxy);

    mesh->J = hthe / Bpxy;

    mesh->g_11 = 1.0 / mesh->g11 + ((I * Rxy) ^ 2);
    mesh->g_22 = (mesh->Bxy * hthe / Bpxy) ^ 2;
    mesh->g_33 = Rxy * Rxy;
    mesh->g_12 = Btxy * hthe * I * Rxy / Bpxy;
    mesh->g_13 = I * Rxy * Rxy;
    mesh->g_23 = Btxy * hthe * Rxy / Bpxy;

    mesh->geometry(); // Calculate other metrics

    // SET VARIABLE LOCATIONS *******************//

    Ni.setLocation(CELL_CENTRE);
    Ti.setLocation(CELL_CENTRE);
    Te.setLocation(CELL_CENTRE);

    // Tritium
    Ni_t.setLocation(CELL_CENTRE);
    Ti_t.setLocation(CELL_CENTRE);
    // Helium
    Ni_he.setLocation(CELL_CENTRE);
    Ti_he.setLocation(CELL_CENTRE);
    //Alpha particle
    N_alpha.setLocation(CELL_CENTER);
    //////////////// BOUNDARIES ///////////////////////
    // Set BOUNDARiES first here, and then apply them every time in physics run/////
    //

    Ni.setBoundary("Ni");
    Ti.setBoundary("Ti");
    Te.setBoundary("Te");

    // Tritium
    Ni_t.setBoundary("Ni_t");
    Ti_t.setBoundary("Ti_t");

    //Alpha
    N_alpha.setBoundary("N_alpha");

    // Helium
    Ni_he.setBoundary("Ni_he");
    Ti_he.setBoundary("Ti_he");

    //Set Boundary for other output variables
//    tau_i.setBoundary("Tn");
//    tau_i_t.setBoundary("Tn");
//    tau_i_he.setBoundary("Tn");
//    tau_e.setBoundary("Tn");
    tau_slowing_se.setBoundary("Tn");
    tau_slowing_down.setBoundary("Tn");
    W_c.setBoundary("Tn");
//    S_fu.setBoundary("Tn");
//    col_eD.setBoundary("Tn");
//    col_eT.setBoundary("Tn");
//    col_ehe.setBoundary("Tn");
//    col_De.setBoundary("Tn");
//    col_DT.setBoundary("Tn");
//    col_Dhe.setBoundary("Tn");
//    col_Te.setBoundary("Tn");
//    col_TD.setBoundary("Tn");
//    col_The.setBoundary("Tn");
//    col_hee.setBoundary("Tn");
//    col_heD.setBoundary("Tn");
//    col_heT.setBoundary("Tn");
    sigma_DT_V_th.setBoundary("Tn");
    temp_sigma_DT_V_th.setBoundary("Tn");
//    P_brem.setBoundary("Tn");
//    P_line_loss.setBoundary("Tn");
//    P_rec.setBoundary("Tn");
//    P_radiation.setBoundary("Tn");
    DDX_Te.setBoundary("DDX_Te");

    V_pinch_alpha.setBoundary("Tn");
//    Test.setBoundary("Tn");





    ///////////// SET EVOLVING VARIABLES //////////////
    //
    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object


    Ne = Ni = Te = Ti = 0.0;
    //NB: it is **NECESSARY** to set values for evolving quantities if not read from grid
    Ni_t = Ti_t = 0.0;
    Ni_he = Ti_he = 0.0;
    N_alpha = 0.0;
//    Test = 0.0;
    S_fu = 0.0;
    pressure =0.0;
//    Prfe = minimum_val;
//    Pauxe = minimum_val;
//    Pauxi = minimum_val;
//    Pauxhe = minimum_val;
    P_fuel_loss = minimum_val;
    P_brem = minimum_val;
    P_cyclo = minimum_val;
    P_radiation = minimum_val;
//    Prad_grid = minimum_val;
    PaD = Pae = Pahe = PaT = 0.0;
    V_pinch.x = 0.0;
    V_pinch_alpha=0.0;


    SOLVE_FOR3(Ni, Te, Ti);
    SOLVE_FOR2(Ni_t, Ti_t);
    SOLVE_FOR2(Ni_he, Ti_he);
    SOLVE_FOR(N_alpha);


    if (load_grid_profiles) {
        output.write("\tInitial profiles of Ti Te Ni are loaded from grid file\n");

        Ni = Ni_grid;
        Ni_t = Ni_grid;
        Ni_he = minimum_val;
        Te = Te_grid;
        Ti = Ti_grid;
        Ti_t = Ti_grid;
        Ti_he = Ti_grid;


    }

    if (load_experiment_profiles) {
        output.write("\tInitial experiment profiles of Ti Te Ni are loaded from grid file\n");
        Te = Te_exp;
        Ti = Ti_exp;
        Ni = Ni_exp;

    }

    // Initialization of Profiles
    // ****** NB: profiles should be initialized after "SOLVE_FOR()" ****//

    // Sources of atom Deuterium and Tritium
    SnD = 0.0;
    SnT = 0.0;
    Pauxe = 0.0;
    Pauxi = 0.0;
    Pauxi_t = 0.0;
    Pauxhe = 0.0;
    Prfe = 0.0;
    Diff_ELM= 0.0;

    for (int jx = 0; jx < mesh->ngx; jx++) {
        x_rela = mesh->GlobalX(jx);

        //BoutReal x0=0.3,w_ped=0.1;
        BoutReal temp = exp(2. * (x_rela - Initfile_x0) / Initfile_w_ped);
        BoutReal x0_nn = 1.02, w_nn = 0.05;
        BoutReal temp2 = exp(-(x_rela - x0_nn) * (x_rela - x0_nn) / w_nn / w_nn);
        BoutReal tempSnD = exp(-(x_rela - SnD_x0) * (x_rela - SnD_x0) / SnD_xwidth / SnD_xwidth);
        BoutReal tempSnT = exp(-(x_rela - SnT_x0) * (x_rela - SnT_x0) / SnT_xwidth / SnT_xwidth);
        BoutReal tempPauxe = exp(-(x_rela - Pauxe_x0) * (x_rela - Pauxe_x0) / Pauxe_xwidth / Pauxe_xwidth);
        BoutReal tempPauxi = exp(-(x_rela - Pauxi_x0) * (x_rela - Pauxi_x0) / Pauxi_xwidth / Pauxi_xwidth);
        BoutReal tempPauxi_t = exp(-(x_rela - Pauxi_t_x0) * (x_rela - Pauxi_t_x0) / Pauxi_t_xwidth / Pauxi_t_xwidth);
        BoutReal tempPauxhe = exp(-(x_rela - Pauxhe_x0) * (x_rela - Pauxhe_x0) / Pauxhe_xwidth / Pauxhe_xwidth);
        BoutReal tempPrfe = exp(-(x_rela - Prfe_x0) * (x_rela - Prfe_x0) / Prfe_xwidth / Prfe_xwidth);
        BoutReal temDiff_ELM = exp(-(x_rela - ped_position) * (x_rela - ped_position) / ped_width / ped_width);

        for (int jy = 0; jy < mesh->ngy; jy++) {
            BoutReal x_psi_l = psixy[jx][jy] - psi_xout_y0;
            BoutReal psi_normal = (psixy[jx][jy] - psi_axis) / (psi_bndry - psi_axis);
            //   output.write("\t psi_normal %e \n",psi_normal );
            BoutReal y_rela = mesh->GlobalY(jy);
            int jy_global = mesh->YGLOBAL(jy);
            BoutReal y0_nn = 0.5, wy_nn = 0.05;
            BoutReal temp2_y = exp(-(y_rela - y0_nn) * (y_rela - y0_nn) / wy_nn / wy_nn);
            for (int jz = 0; jz < mesh->ngz; jz++) {
                if (!load_grid_profiles) {

                    if (initial_profile_exp) {
                        Ni[jx][jy][jz] = Ni_edge + 2. * Ni_core / (1. + temp);
//                        Test[jx][jy][jz] = Ni[jx][jy][jz];
                        // Ti[jx][jy][jz]=Te_edge+Te_core/(1.+temp);
                        Ni_t[jx][jy][jz] = Ni_edge + 2. * Ni_core / (1. + temp);
                        Ti[jx][jy][jz] = Ti_edge + 2. * Ti_core / (1. + temp);
                        Ti[jx][jy][jz] = Ti_edge + 2. * Ti_core / (1. + temp);
                        Te[jx][jy][jz] = Ti_edge + 2. * Ti_core / (1. + temp);
                        // Te[jx][jy][jz]=Te_edge+Te_core/(1.+temp);
                    }

                    if (initial_profile_linear) {
                        Ni[jx][jy][jz] = Ni_edge + dNidx_xin_au * x_psi_l;
//                        Test[jx][jy][jz] = Ni[jx][jy][jz];
                        Ni_t[jx][jy][jz] = Ni_edge + dNidx_xin_au * x_psi_l;
                        Ti[jx][jy][jz] = Ti_edge + dTidx_xin_au * x_psi_l;
                        Ti_t[jx][jy][jz] = Ti_edge + dTidx_xin_au * x_psi_l;
                        Te[jx][jy][jz] = Te_edge + dTedx_xin_au * x_psi_l;
                    }


                    if (initial_profile_square) {
                        Ni[jx][jy][jz] = -2.311e-14 * exp(32.34 * psi_normal) + 4.89 * exp(-0.4763 * psi_normal);
                       // Ni[jx][jy][jz] = 5-5*psi_normal;
//                        Test[jx][jy][jz] = Ni[jx][jy][jz];
                        Ni_t[jx][jy][jz] = Ni[jx][jy][jz];

                        Ti[jx][jy][jz] = 100. * (-124.3 * psi_normal * psi_normal * psi_normal * psi_normal * psi_normal
                                                 + 321.1 * psi_normal * psi_normal * psi_normal * psi_normal
                                                 - 344.6 * psi_normal * psi_normal * psi_normal
                                                 + 198.3 * psi_normal * psi_normal
                                                 - 77.87 * psi_normal
                                                 + 28.34);
                        Ti_t[jx][jy][jz] = Ti[jx][jy][jz];
                        Te[jx][jy][jz] = Ti[jx][jy][jz];
//                        Ni_he[jx][jy][jz]=Ti[jx][jy][jz];
                        Ti_he[jx][jy][jz]=Ti[jx][jy][jz];
                    }


                }


                //~~~~~~~~~~~~~~~~~~
                // INITIALIZE
                //__________________
                // Ni_t[jx][jy][jz] = minimum_val;
                // Ti_t[jx][jy][jz] = minimum_val;
                Ni_he[jx][jy][jz] = minimum_val;
                //Ti_he[jx][jy][jz] = minimum_val;


                if (fueling_atomD) SnD[jx][jy][jz] = SnD_amp * tempSnD;
                if (fueling_atomT) SnT[jx][jy][jz] = SnT_amp * tempSnT;
                if (auxiliary_heating_e) Pauxe[jx][jy][jz] = Pauxe_amp * tempPauxe;
                if (auxiliary_heating_i) Pauxi[jx][jy][jz] = Pauxi_amp * tempPauxi;
                if (auxiliary_heating_i_t) Pauxi_t[jx][jy][jz] = Pauxi_t_amp * tempPauxi_t;
                if (auxiliary_heating_he) Pauxhe[jx][jy][jz] = Pauxhe_amp * tempPauxhe;
                if (auxiliary_rfheating_e) Prfe[jx][jy][jz] = Prfe_amp * tempPrfe;
                if (trigger_ELM) Diff_ELM[jx][jy][jz] = Diff_ELM_amp*temDiff_ELM;
            }
        }
    }



    //End of Initialization of Profiles

    // Set step functions of diffusion coefficients
    // //NB: it is **NECESSARY** to set initial Zero values
    F_r_pinch = 0.;
    F_r_pinch_2=0.;

    for (int jx = 0; jx < mesh->ngx; jx++) {
        for (int jy = 0; jy < mesh->ngy; jy++) {
            int jy_global = mesh->YGLOBAL(jy);
            BoutReal psi_normal = (psixy[jx][jy] - psi_axis) / (psi_bndry - psi_axis);
            for (int jz = 0; jz < mesh->ngz; jz++) {
                if (sqrt(psi_normal) < rho_coef-0.01) {
                    // F_r_pinch[jx][jy][jz]=(0.1+(0.9*psi_normal*psi_normal))*((0.9*exp(-psi_normal/0.92/0.92)^100.0)+0.1);
                    F_r_pinch[jx][jy][jz] = (0.25 + (0.75 * psi_normal * psi_normal )) *
                                             ((1.-A_coef) * exp(-pow((psi_normal / rho_coef / rho_coef), p_coef)) + A_coef);
                    F_r_pinch_2[jx][jy][jz]=F_r_pinch[jx][jy][jz];
                } else {
                    F_r_pinch[jx][jy][jz] = (0.25 + (0.75 * rho_coef * rho_coef * rho_coef * rho_coef)) *
                                            (((1.-A_coef) * exp(-pow((rho_coef / rho_coef), (2*p_coef)))) + A_coef) * (1 - sqrt(psi_normal)) *
                                            (1 - sqrt(psi_normal)) / (1-rho_coef) / (1-rho_coef);
                    F_r_pinch_2[jx][jy][jz]=0.0;
                }


            }
        }
    }

    ///////////// ADD OUTPUT VARIABLES ////////////////
    //
    // Add any other variables to be dumped to file

    SAVE_ONCE5(Te_x, Ti_x, Ni_x, Lbar, tbar); // Normalisation factors
    SAVE_ONCE2(bmag, Vi_x);

    // electron density calculation with quasi-neutral condition
    Ne = Ni + Ni_t + 2. * Ni_he + 2. * N_alpha;

    // Set flux limit for kappa
    V_th_e = 4.19e5 * sqrt(Te * Te_x);
    V_th_i = 9.79e3 * sqrt(Ti * Te_x / AA);
    output.write("\tion thermal velocity: %e -> %e [m/s]\n", min(V_th_i), max(V_th_i));
    output.write("\telectron thermal velocity: %e -> %e [m/s]\n", min(V_th_e), max(V_th_e));
    V_th_e /= Lbar / tbar;
    V_th_i /= Lbar / tbar;
    output.write("\tNormalized ion thermal velocity: %e -> %e [Lbar/tbar]\n", min(V_th_i), max(V_th_i));
    output.write("\tNormalized electron thermal velocity: %e -> %e [Lbar/tbar]\n", min(V_th_e), max(V_th_e));



//output.write("\t Test 01_03  By Lin\n");

    dump.add(psixy, "psixy", 0);


    if (fueling_atomD) dump.add(SnD, "SnD", 1);
    if (fueling_atomT) dump.add(SnT, "SnT", 1);
    if (auxiliary_rfheating_e) dump.add(Prfe,"Prfe",1);
    if (auxiliary_heating_e) dump.add(Pauxe,"Pauxe",1);
    if (auxiliary_heating_i) dump.add(Pauxi,"Pauxi",1);
    if (auxiliary_heating_i_t) dump.add(Pauxi_t,"Pauxi_t",1);
    if (auxiliary_heating_he) dump.add(Pauxhe,"Pauxhe",1);
    dump.add(Ne, "Ne", 1);
    dump.add(Diff_psi, "Diff_psi", 1);
    dump.add(V_pinch.x, "V_pinch", 1);
    dump.add(sqrt_psi, "sqrt_psi", 1);
    dump.add(psi_norm,"psi_norm",1);
    dump.add(kappa_psi, "kappa_psi", 1);
//    dump.add(Diff_he, "Diff_he", 1);
    dump.add(Diff_psi_para, "Diff_psi_para", 1);
    dump.add(kappa_psi_para,"kappa_psi_para",1);
    dump.add(S_fu, "S_fu", 1);
    dump.add(PaD,"PaD",1);
    dump.add(Pae,"Pae",1);
    dump.add(Pahe,"Pahe",1);
    dump.add(PaT,"PaT",1);
    dump.add(tau_slowing_se, "tau_slowing_se", 1);
    dump.add(tau_slowing_down, "tau_slowing_down", 1);
    dump.add(Diff_alpha,"Diff_alpha",1);
    dump.add(V_pinch_alpha,"V_pinch_alpha",1);
    dump.add(alpha_pressure,"alpha_pressure",1);

//  dump.add(col_eD,"col_eD",1);
//  dump.add(col_eT,"col_eT",1);
//  dump.add(col_ehe,"col_ehe",1);
    dump.add(col_De,"col_De",1);
    dump.add(col_DT,"col_DT",1);
    dump.add(col_Dhe,"col_Dhe",1);
    dump.add(col_Te,"col_Te",1);
//  dump.add(col_TD,"col_TD",1);
    dump.add(col_The,"col_The",1);
    dump.add(col_hee,"col_hee",1);
    dump.add(col_heT,"col_heT",1);
    dump.add(col_heD,"col_heD",1);
    dump.add(col_hehe,"col_hehe",1);

    dump.add(P_brem, "P_brem", 1);
    dump.add(P_cyclo, "P_cyclo", 1);
//    dump.add(P_line_loss, "P_line_loss", 1);
//    dump.add(P_rec, "P_rec", 1);
    dump.add(P_radiation, "P_radiation", 1);
    dump.add(sigma_DT_V_th, "sigma_DT_V_th", 1);
//    dump.add(tau_i, "tau_i", 1);
//    dump.add(tau_e, "tau_e", 1);
//    dump.add(tau_i_t, "tau_i_t", 1);
//    dump.add(tau_i_he, "tau_i_he", 1);
//dump.add(Test,"Test",1);
//    dump.add(Test1, "Test1", 1);
//    dump.add(Test2, "Test2", 1);
//    dump.add(Test3, "Test3", 1);
//    dump.add(Test4, "Test4", 1);
//    dump.add(lambda_De, "lambda_De", 1);

    dump.add(Z_average,"Z_average",1);
    dump.add(W_c,"W_c",1);
    dump.add(fi,"fi",1);
    dump.add(fe,"fe",1);
    dump.add(Prfe_grid,"Prfe_grid",1);
    dump.add(Prad_grid,"Prad_grid",1);
    dump.add(Pauxi_grid,"Pauxi_grid",1);
    dump.add(Pauxe_grid,"Pauxe_grid",1);
    dump.add(tau_e_real,"tau_e_real",1);
    dump.add(tau_e98,"tau_e98",1);
    dump.add(P_fuel_loss,"P_fuel_loss",1);
    dump.add(P_98,"P_98",1);
    dump.add(P_fuel_loss_total,"P_fuel_loss_total",1);
    dump.add(P_radiation_total,"P_radiation_total",1);
    dump.add(P_cd_total,"P_cd_total",1);
    dump.add(P_fusion_total, "P_fusion_total",1);
    dump.add(Prfe_grid_total,"Prfe_grid_total",1);
    dump.add(Pauxi_grid_total,"Pauxi_grid_total",1);
    dump.add(Pauxe_grid_total,"Pauxe_grid_total",1);

    dump.add(pellet_rate,"pellet_rate",1);
    dump.add(pellet_rp_profile,"pellet_rp_profile",1);
    dump.add(pellet_dN_rate,"pellet_dN_rate",1);
    dump.add(pellet_dN_dr,"pellet_dN_dr",1);
    dump.add(S_pellet,"S_pellet",1);
    dump.add(S_pellet_total,"S_pellet_total",1);
    dump.add(n_average,"n_average",1);
    dump.add(deposit_depth,"deposit_depth",1);
    dump.add(pellet_life_time,"pellet_life_time",1);
    dump.add(pellet_timestep,"pellet_timestep",1);
    dump.add(pressure,"pressure",1);
    dump.add(T_pellet,"T_pellet",1);
    dump.add(r_cloud,"r_cloud",1);
    dump.add(L_c,"L_c",1);
    dump.add(pellet_displacement,"pellet_displacement",1);
    dump.add(add_pellet_displacement,"add_pellet_displacement",1);
    dump.add(pellet_drift_psi,"pellet_drift_psi",1);
    dump.add(R_new,"R_new",1);
}


int physics_run(BoutReal t) {

    // Communicate variables
    mesh->communicate(Ni, Te, Ti);
    mesh->communicate(Ni_t, Ti_t);
    mesh->communicate(Ni_he, Ti_he);
    mesh->communicate(N_alpha);
    // NB: Intermediate variables calculated with Grad operators are all necessary to be communicated
    // after being calculated

    Ni.applyBoundary();
    Ti.applyBoundary();

    Ni_t.applyBoundary();
    Ti_t.applyBoundary();

    Ni_he.applyBoundary();
    Ti_he.applyBoundary();

    Te.applyBoundary();

    N_alpha.applyBoundary();

//    Test.applyBoundary();
    V_pinch_alpha.applyBoundary();

//  output.write("\t Test 11 By Lin\n");
    // electron density calculation with quasi-neutral condition

    Ne = Ni + Ni_t + 2. * Ni_he + 2. * N_alpha;

    //smooth noisies

    //*****@!!@*****
    // NB: Any value re-assignment should be given HERE ONLY!
    //*****@!!@*****

    temp_Ni = field_larger(Ni, minimum_val);
    temp_Ti = field_larger(Ti, minimum_val);

    temp_Ne = field_larger(Ne, minimum_val);
    temp_Te = field_larger(Te, minimum_val);


    temp_Ni_t = field_larger(Ni_t, minimum_val);
    temp_Ti_t = field_larger(Ti_t, minimum_val);

    temp_Ni_he = field_larger(Ni_he, minimum_val);
    temp_Ti_he = field_larger(Ti_he, minimum_val);


    Ni = field_larger(Ni, minimum_val);
    Ti = field_larger(Ti, minimum_val);

    Ne = field_larger(Ne, minimum_val);
    Te = field_larger(Te, minimum_val);

    Ni_t = field_larger(Ni_t, minimum_val);
    Ti_t = field_larger(Ti_t, minimum_val);

    Ni_he = field_larger(Ni_he, minimum_val);
    Ti_he = field_larger(Ti_he, minimum_val);

    N_alpha = field_larger(N_alpha, minimum_val);

    //Critical energy W_c: in order to caculate slowing down time
    Z_average = temp_Ni*Z_D*Z_D*A_b/Ne/A_D+temp_Ni_t*Z_T*Z_T*A_b/Ne/A_T+temp_Ni_he*Z_He*Z_He*A_b/Ne/A_He;
     W_c = 14.8*(Z_average^0.6667)*pow(A_b,0.3333)*(temp_Te*Te_x/1000);//keV
     fi=integral_fi(0,W_b0,1000);
     fe=integral_fe(0,W_b0,1000);
     //W_c= (temp_Te*Te_x/1000);
    //Slowing time between alpha and electron
    ln_lambda = 29.27 + 1.5*log(temp_Te*Ti_x)-0.5*log(temp_Ne*density_unit*Ni_x);
    //tau_slowing_se = 1.17e18 * ((temp_Te * Ti_x / 1000.) ^ 1.5) / (temp_Ne * density_unit) / tbar;
    tau_slowing_se = 0.2 * A_b * ((temp_Te*Te_x/1000)^1.5) / Z_b / Z_b / (temp_Ne * density_unit*Ni_x/1e20)/ln_lambda/tbar;
    //alpha particle slowing down time
   tau_slowing_down = tau_slowing_se/3*log((pow(W_b0,1.5)+W_c^1.5)/((W_c^1.5)+(temp_Ti_he*Te_x/1000.)^1.5));


    // Collision rate calculation between e, D,T and He
    lambda_De = 24. - log((sqrt(temp_Ne * Ni_x * density_unit * 1.e-6)) / (temp_Te * Te_x));
    lambda_Te = lambda_De;
    lambda_hee = lambda_De;
//    // output.write("\t Test 131 By Lin\n");
    lambda_DT = 23. - log(5. / (3.0 * temp_Ti * Te_x + 2.0 * temp_Ti_t * Te_x) *
                          ((temp_Ni * Ni_x * density_unit * 1.e-6 / (temp_Ti * Te_x) +
                            temp_Ni_t * Ni_x * density_unit * 1.e-6 / (temp_Ti_t * Te_x)) ^ 0.5));
    lambda_Dhe = 23. - log(6. / (4.0 * temp_Ti * Te_x + 2.0 * temp_Ti_he * Te_x) *
                           ((temp_Ni * Ni_x * density_unit * 1.e-6 / (temp_Ti * Te_x) +
                             4.0 * temp_Ni_he * Ni_x * density_unit * 1.e-6 / (temp_Ti_he * Te_x)) ^ 0.5));
    lambda_The = 23. - log(7. / (3.0 * temp_Ti_he * Te_x + 4 * temp_Ti_t * Te_x) *
                           ((temp_Ni_t * Ni_x * density_unit * 1.e-6 / (temp_Ti_t * Te_x) +
                             4.0 * temp_Ni_he * Ni_x * density_unit * 1.e-6 / (temp_Ti_he * Te_x)) ^ 0.5));
    lambda_hehe = 23. - log(8. * (sqrt(temp_Ni_he * Ni_x * density_unit * 1.e-6)) / ((temp_Ti_he * Te_x) ^ 1.5));
//note: The formula unit in NRL formula is cgs.
    col_De = tbar * 1.8e-19 * sqrt(mass_e * 2 * mass_p) * (temp_Ne * Ni_x * density_unit * 1.e-6) *
            lambda_De / (2 * mass_p * Te * Te_x + mass_e * Ti * Te_x) /
            sqrt(2 * mass_p * Te * Te_x + mass_e * Ti * Te_x);
    col_DT = tbar * 1.8e-19 * sqrt(3 * mass_p * 2 * mass_p) * (temp_Ni_t * Ni_x * density_unit * 1.e-6) *
            lambda_DT / (2 * mass_p * Ti_t * Te_x + 3 * mass_p * Ti * Te_x) /
            sqrt(2 * mass_p * Ti_t * Te_x + 3 * mass_p * Ti * Te_x);
    col_Dhe = tbar * 1.8e-19 * sqrt(4 * mass_p * 2 * mass_p) * 4 * (temp_Ni_he * Ni_x * density_unit * 1.e-6) *
            lambda_Dhe / (2 * mass_p * Ti_he * Te_x + 4 * mass_p * Ti * Te_x) /
            sqrt(2 * mass_p * Ti_he * Te_x + 4 * mass_p * Ti * Te_x);
    col_Te = tbar * 1.8e-19 * sqrt(mass_e * 3 * mass_p) * (temp_Ne * Ni_x * density_unit * 1.e-6) *
            lambda_Te / (3 * mass_p * Te * Te_x + mass_e * Ti_t * Te_x) /
            sqrt(3 * mass_p * Te * Te_x + mass_e * Ti_t * Te_x);
    col_The = tbar * 1.8e-19 * sqrt(3 * mass_p * 4 * mass_p) * 4 *(temp_Ni_he * Ni_x * density_unit * 1.e-6) *
              lambda_The / (3 * mass_p * Ti_he * Te_x + 4 * mass_p * Ti_t * Te_x) /
              sqrt(3 * mass_p * Ti_t * Te_x + 4 * mass_p * Ti_he * Te_x);
    col_hee =  tbar * 1.8e-19 * sqrt(mass_e * 4 * mass_p) * 4 * (temp_Ne * Ni_x * density_unit * 1.e-6) *
               lambda_hee / (4 * mass_p * Te * Te_x + mass_e * Ti_he * Te_x) /
               sqrt(4 * mass_p * Te * Te_x + mass_e * Ti_he * Te_x);
    col_hehe = tbar * 1.8e-19 * sqrt(4 * mass_p * 4 * mass_p) * 4 * 4 *(temp_Ni_he * Ni_x * density_unit * 1.e-6) *
               lambda_hehe / (4 * mass_p * Ti_he * Te_x + 4 * mass_p * Ti_he * Te_x) /
               sqrt(4 * mass_p * Ti_he * Te_x + 4 * mass_p * Ti_he * Te_x);
//    col_TD = col_DT * Ni * Mi / temp_Ni_t / Mi_t;
    col_heT = col_The * temp_Ni_t / temp_Ni_he ;
    col_heD = col_Dhe * temp_Ni /temp_Ni_he ;
//    col_eD = col_De * Ni * Mi / temp_Ne / Me;
//    col_eT = col_Te * Ni_t * Mi_t / temp_Ne / Me;
//    col_ehe = col_hee * Ni_he * Mi_he / temp_Ne / Me;

//Classical diffusion coefficiency of  plasma
//Diff_eD=Te*Me*col_eD/tbar/tbar/(mesh->Bxy*mesh->Bxy)*3.344e-8/1.602*3.344e-8/1.602;
//Diff_eT=Te*Me*col_eT/tbar/tbar/(mesh->Bxy*mesh->Bxy)*3.344e-8/1.602*3.344e-8/1.602;
//Diff_ehe=Te*Me*col_ehe/tbar/tbar/(mesh->Bxy*mesh->Bxy)*3.344e-8/1.602*3.344e-8/1.602;

//radiation loss   From NRL
//    P_brem = 1.69e-32*(Ne*1.e-6*Ni_x*1e19)*sqrt(temp_Te*Te_x)*(temp_Ni*1.e-6*Ni_x*1e19+temp_Ni_t*1.e-6*Ni_x*1e19+Z_He*Z_He*temp_Ni_he*1e19)
//             *tbar/Ti_x/1.602/Ni_x;  // normalization
//    P_cyclo = 6.21e-28*B_98*B_98*1e8*(Ne*1.e-6*Ni_x*1e19)*(temp_Te*Te_x)
//             *tbar/Ti_x/1.602/Ni_x;  // normalization
//    P_brem = 4.8e1 * (4 * temp_Ni_he + temp_Ni + temp_Ni_t) * temp_Ne * ((temp_Te * Te_x * 1.e-3) ^ 0.5) * tbar / Te_x /
//             1.602;
//    P_line_loss =
//            1.8 * (16 * temp_Ni_he + temp_Ni + temp_Ni_t) * temp_Ne * tbar / ((temp_Te * Te_x * 1.e-3) ^ 0.5) / Te_x /
//            1.602;
//    P_rec = 4.1e-2 * (64 * temp_Ni_he + temp_Ni + temp_Ni_t) * temp_Ne * tbar / ((temp_Te * Te_x * 1.e-3) ^ 1.5) /
//            Te_x / 1.602;
//    P_radiation = Prad_grid;
//    P_radiation = P_radiation*1e6 // back to m^-3
//                  /1e19/Ti_x*tbar/ee;// dimensionless
    //fusion reaction rate
//  sigma_DT_V_th = tbar*(Ni_x*density_unit*1.e-6)*3.68e-12/((temp_Ti*Te_x*0.5*1.e-3+temp_Ti_t*Te_x*0.5*1.e-3)^0.6667)*exp(-19.94/((temp_Ti*Te_x*0.5*1.e-3+temp_Ti_t*Te_x*0.5*1.e-3)^0.3333));
    sigma_DT_V_th = tbar * (Ni_x * density_unit * 1.e-6) *
                    (-2.926e-20 * (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) *
                     (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) *
                     (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) +
                     1.931e-18 * (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) *
                     (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) -
                     6.047e-18 * (temp_Ti * Te_x * 0.5 * 1.e-3 + temp_Ti_t * Te_x * 0.5 * 1.e-3) + 3.48e-18);
    temp_sigma_DT_V_th = field_larger(sigma_DT_V_th, minimum_val);//set the min value of nuclear fusion

    S_fu = temp_sigma_DT_V_th * Ni * Ni_t;// Deuterium and Tritium fusion sources term
// caculate the n_average
    n_average=average(mesh->J,Ne);
    S_fu_total=total(mesh->J,S_fu);
//    P_cd = total(mesh->J,Pauxe)*2*PI*1.602*Ti_x/tbar/1e6 +
//            total(mesh->J,Pauxi)*2*PI*1.602*Ti_x/tbar/1e6 +
//            total(mesh->J,Pauxi_t)*2*PI*1.602*Ti_x/tbar/1e6 +
//            total(mesh->J,Prfe)*2*PI*1.602*Ti_x/tbar/1e6;
    P_cd =    total(mesh->J,Pauxe_grid)*2*PI*1.602*Ti_x/tbar/1e6
            + total(mesh->J,Pauxi_grid)*2*PI*1.602*Ti_x/tbar/1e6
//            total(mesh->J,Pauxi_t)*2*PI*1.602*Ti_x/tbar/1e6 +
            + total(mesh->J,Prfe_grid)*2*PI*1.602*Ti_x/tbar/1e6
            ;
    P_radiation_total = total(mesh->J,Prad_grid)*2*PI*1.602*Ti_x/tbar/1e6;
    P_fusion_total = S_fu_total/tbar*2*PI*3.5*1.602;
//    P_98=P_cd+S_fu_total/tbar*2*PI*3.5*1.602 - total(mesh->J,P_radiation)*2*PI*1.602*Ti_x/tbar/1e6;
    P_98 = P_cd + P_fusion_total
            - P_radiation_total;
    //        - P_fuel_loss_total;
    M_98=(2*total(mesh->J,Ni)+3*total(mesh->J,Ni_t)+4*total(mesh->J,Ni_he)+4*total(mesh->J,N_alpha))/(total(mesh->J,Ni)+total(mesh->J,Ni_t)+total(mesh->J,Ni_he)+total(mesh->J,N_alpha));
    n_98=n_average;
    P_fuel_loss = 1.5*(SnD+SnT)*Te+1.5*SnD*Ti+1.5*SnT*Ti_t;
    P_fuel_loss_total = total(mesh->J,P_fuel_loss)*2*PI*1.602*Ti_x/tbar/1e6;
    P_cd_total = P_cd;

    Prfe_grid_total = total(mesh->J,Prfe_grid)*2*PI*1.602*Ti_x/tbar/1e6;
    Pauxi_grid_total = total(mesh->J,Pauxi_grid)*2*PI*1.602*Ti_x/tbar/1e6;
    Pauxe_grid_total = total(mesh->J,Pauxe_grid)*2*PI*1.602*Ti_x/tbar/1e6;

//  output.write("M_98 is %e\n",M_98);
//  output.write("n_average is %e\n",n_average);
//   output.write("S_fu_average is %e\n",S_fu_total/tbar*2*PI*3.5*1.602);
//   output.write("P_98 is %e\n",P_98);
//   output.write("P_cd is %e\n",P_cd);
//   output.write("P_fusion_total is %e\n",P_fusion_total);
//   output.write("P_radiation_total is %e\n",P_radiation_total);
//caculate tau_e98
    tau_e98 = 0.0562 * pow(Ip_98, 0.93) * pow(B_98, 0.15) * pow(P_98, -0.69) * pow(n_98, 0.41) * pow(M_98, 0.19) *
              pow(R_98, 1.97) * pow(epsi_98, 0.58) * pow(kappa_98, 0.78);
//  tau_e98=0.0562*Ip_98^0.93*B_98^0.15*P_98^-0.69*n_98^0.41*M_98^0.19*R_98^1.97*epsi_98^0.58*kappa_98^0.78;
//output.write("tau_e98 is %e\n",tau_e98);

//Update kappa_factor
    Wth = 1.5*(Ni*Ti+Ni_t*Ti_t+Ni_he*Ti_he+Ne*Te);
    Wth_total = total(mesh->J,Wth)*2*PI*1.602*Ti_x/1e6;
    tau_e_real = Wth_total/P_98;//normalization
//    if(abs(tau_e_real/tau_e98-1)>=0.1){
//        kappa_psi_factor = kappa_psi_factor*sqrt(tau_e_real/tau_e98);
//   }


//diffusion coefficiency depended on psi
    sqrt_psi = sqrt((psixy - psi_axis) / (psi_bndry - psi_axis));
    psi_norm = (psixy - psi_axis) / (psi_bndry - psi_axis);
    //Scaling thermal conductivity
//    kappa_psi = kappa_psi_factor  * 2. * a_radius * a_radius  / tau_e98 *
//                (0.25 + (0.75 *  sqrt_psi * sqrt_psi * sqrt_psi * sqrt_psi)) *
//                ((1. - A_coef) * exp(-((sqrt_psi / rho_coef) ^ (2*p_coef))) + A_coef) * tbar / Lbar / Lbar;
    // Diff_psi=2*a_radius*a_radius/tau_e98*(0.25+0.75*sqrt_psi^4)*((1-0.15)*exp(-(sqrt_psi/0.93)^40)+0.15);
    //Diff_psi=Diff_psi_factor*0.35*2*a_radius*a_radius/tau_e98*(0.25+0.75*sqrt_psi^4)*((1.-0.15)*exp(-((sqrt_psi/0.93)^40.0))+0.15);
//    Diff_psi = Diff_psi_factor  *kappa_psi;
    //Diff_psi= 1.0*tbar/Lbar/Lbar*Diff_psi_factor;
//    Diff_he = Diff_psi;
    Diff_psi_para = 10000 * tbar / Lbar / Lbar ;
    kappa_psi_para = 10000* tbar / Lbar / Lbar;
    //pinch velocity
    // V_pinch.x=-2*0.65*sqrt_psi/a_radius*Diff_psi_factor*0.35*2.*a_radius*a_radius*1.2/tau_e98*tbar/Lbar*F_r_pinch;
    V_pinch.x =
            -2 * 0.65 * sqrt_psi *kappa_psi_factor* Diff_psi_factor  * 2.  / tau_e98 * tbar / Lbar * F_r_pinch * a_radius*pinch_factor;
    V_pinch.y = 0;
    V_pinch.z = 0;

    //alpha pressure
    alpha_pressure = tau_slowing_se * S_fu / 2 * fe * 1.602 ; //kPa

    //trigger ELM
    pressure = (Ni*Ti+Ni_t*Ti_t+Ni_he*Ti_he+Ne*Te)*Te_x*1.602/1000 + alpha_pressure ;
    if(trigger_ELM)
        {
        ELM_happen = test_ELM_happen();
            }
    if(ELM_happen)
        {
	//Scaling thermal conductivity
        kappa_psi = kappa_psi_factor  * 2. * a_radius * a_radius  / tau_e98 *
                (0.25 + (0.75 *  sqrt_psi * sqrt_psi * sqrt_psi * sqrt_psi)) *
                ((1. - A_coef) * exp(-((sqrt_psi / rho_coef) ^ (2*p_coef))) + A_coef) * tbar / Lbar / Lbar + Diff_ELM*tbar/Lbar/Lbar;        
        Diff_psi = Diff_psi_factor  *kappa_psi;
//	Diff_he = Diff_psi + Diff_ELM*tbar/Lbar/Lbar; 
        }
     else
     {
        kappa_psi = kappa_psi_factor  * 2. * a_radius * a_radius  / tau_e98 *
                (0.25 + (0.75 *  sqrt_psi * sqrt_psi * sqrt_psi * sqrt_psi)) *
                ((1. - A_coef) * exp(-((sqrt_psi / rho_coef) ^ (2*p_coef))) + A_coef) * tbar / Lbar / Lbar;
        Diff_psi = Diff_psi_factor  *kappa_psi;
 //       Diff_he = Diff_psi;
     }

//Alpha particle's diffusion coefficient

    //Diff_alpha = Diff_psi*(0.02+4.5*(temp_Te*Te_x/1000./W_b0)+8*(temp_Te*Te_x/1000./W_b0)^2+350*(temp_Te*Te_x/1000./W_b0)^3);
    Diff_alpha = Diff_psi*(0.02+4.5*(temp_Te*Te_x/1000./W_b0)+8.*((temp_Te*Te_x/1000./W_b0)^2.)+350*((temp_Te*Te_x/1000./W_b0)^3.));
    DDX_Te = DDX(temp_Te);
    mesh->communicate(DDX_Te);
    DDX_Te.applyBoundary();

//    V_pinch_alpha = -1.5 * sqrt(mesh->g11) * DDX_Te / temp_Te *(1/((1+((W_c/W_b0)^1.5))*log((((W_b0/W_c)^1.5)+1)))-1)*Diff_alpha*Lbar;

    V_pinch_alpha = -1.5 * sqrt(mesh->g11) * DDX_Te / temp_Te *(1/((1+((W_c/W_b0)^1.5))*log(((W_b0/W_c)^1.5)+1))-1)*Lbar
                    *(0.02+4.5*(temp_Te*Te_x/1000./W_b0)+8.*((temp_Te*Te_x/1000./W_b0)^2.)+350*((temp_Te*Te_x/1000./W_b0)^3.))
                    *Diff_psi_factor  * 2. * a_radius * a_radius / tau_e98*F_r_pinch_2* tbar / Lbar / Lbar;
                    ;



    PaD = S_fu*(E_alpha-Ti_he)*fi/(fi+fe)*col_heD/
            (col_heD+col_heT+col_hehe);
    PaT = S_fu*(E_alpha-Ti_he)*fi/(fi+fe)*col_heT/
            (col_heD+col_heT+col_hehe);
    Pahe = S_fu*(E_alpha-Ti_he)*fi/(fi+fe)*col_hehe/
            (col_heD+col_heT+col_hehe);
    Pae = S_fu*(E_alpha-Ti_he)*fe/(fi+fe);

    //pellet parameter
//    pellet_rate=Pellet_rp_ablate_rate(np,Mi_p,vp,mid,Ne,Te,Rxy);
      pellet_rate=Pellet_rp_ablate_rate(np,Mi_p,vp,mid);
      pellet_rp_profile=Pellet_rp_profile(pellet_rate,rp0);
      T_pellet=1.88e-9*((Mi_p/Te/10.)^0.3333)*((Ne*1e19*1e-6*pellet_rp_profile*log(2*Te*10/7.5))^0.66667);  //eV
      r_cloud=1.54e4*((10*Te/Mi_p)^0.1333333)*(4*1.66667*T_pellet+13.6+2.2)*(pellet_rp_profile^0.66667)
              *((Ne*1e19*1e-6*log(2*Te*10/7.5))^(-0.3333));  //cm
      L_c=(r_cloud*R_98*100)^0.5;  //cm
      pellet_dN_rate=Pellet_dN_ablate_rate(np,Mi_p,vp);
      pellet_dN_dr=Pellet_dN_dr(vp);
      pellet_displacement = Pellet_displacement()/(psi_bndry-psi_axis);
      add_pellet_displacement = Add_pellet_displacement();
      pellet_drift_psi=psi_norm-add_pellet_displacement;
      R_new=Interp_R_new();
      S_pellet=Pellet_dN_profile();  //need normalization
      S_pellet_total=total(mesh->J,S_pellet)*2*PI;
      deposit_depth = find_rp_deposit_depth();
//      output.write("deposit_depth is %f\n",deposit_depth);
      pellet_life_time = deposit_depth/vp*100;
      pellet_timestep = int(1/pellet_frequency/tbar);
//      test();

    //Set Boundary and Apply Boundary for output variables
//    S_fu.applyBoundary();
//    col_eD.applyBoundary();
//    col_eT.applyBoundary();
//    col_ehe.applyBoundary();
//    col_De.applyBoundary();
//    col_DT.applyBoundary();
//    col_Dhe.applyBoundary();
//    col_Te.applyBoundary();
//    col_The.applyBoundary();
//    col_TD.applyBoundary();
//    col_hee.applyBoundary();
//    col_heT.applyBoundary();
//    col_heD.applyBoundary();
//    tau_i.applyBoundary();
//    tau_i_t.applyBoundary();
//    tau_i_he.applyBoundary();
//    tau_e.applyBoundary();
    sigma_DT_V_th.applyBoundary();
    temp_sigma_DT_V_th.applyBoundary();
//    P_brem.applyBoundary();
//    P_line_loss.applyBoundary();
//    P_rec.applyBoundary();
//    P_radiation.applyBoundary();

//    if (t==10.0)
//        {
//        Ni += S_pellet;
//        output.write("Add pellet");
//        }
//    if (t==10.0) Ti -= Ti*S_pellet/temp_Ni;
//    if (t==10.0) Ni_t += S_pellet;
//    if (t==10.0) Ti_t -= S_pellet*Ti_t/temp_Ni_t;
//    if (t==10.0) Te -= S_pellet*Te/temp_Ne;
//    output.write("t is %f\n",t);
//*************************************************
    // DENSITY EQUATION        --Ni--    for deuterium
//*************************************************
    // output.write("Now updata Ni \n");

    //ddt(Ni)= 0.;

    ddt(Ni) =   Diff_psi * Delp2(Ni) + DDX(Ni) * DDX(Diff_psi) * mesh->g11
              - Ni * DDX(V_pinch.x) * sqrt(mesh->g11) - V_pinch.x * DDX(Ni) * sqrt(mesh->g11)
              + Diff_psi_para * Grad2_par2(Ni)
              // - Ni*Grad_par(Vi)
              // - Vpar_Grad_par(Vi,Ni)
              // + Si_p
              + SnD//for no atom case
            ;
    if (burning) ddt(Ni) += -S_fu;
    if (t>=0 && int(t)%pellet_timestep>=0 && int(t)%pellet_timestep<=timestep) ddt(Ni) += S_pellet/2/timestep;
    if (t>=0) ddt(Ni) += -SnD;
//***********************************************
    // ION TEMPERATURE   ---Ti---     for deuterium
//***********************************************
  ddt(Ti) =   0.66667*kappa_psi_para*Grad2_par2(Ti)     //para term
             - Ti*(SnD)/temp_Ni // source term
             - 0.66667*Ti*(-DDX(Ni) * DDX(Diff_psi) * mesh->g11/temp_Ni-Diff_psi * Delp2(Ni)/temp_Ni
             + Diff_psi*DDX(Ni)*DDX(Ni)*mesh->g11/temp_Ni/temp_Ni + DDX(V_pinch.x)*sqrt(mesh->g11))  //0.6667*T*Gard(u)
             + 0.66667*kappa_psi*Delp2(Ti) + 0.66667*DDX(Ti)*DDX(kappa_psi)*mesh->g11
             + 0.66667*kappa_psi*DDX(Ni)*DDX(Ti)*mesh->g11/temp_Ni // conductive term
             + Diff_psi*DDX(Ni)*DDX(Ti)*mesh->g11/temp_Ni - V_pinch.x*DDX(Ti)*sqrt(mesh->g11)
             + col_De*(Te-Ti)
             + col_DT*(Ti_t-Ti)
             + col_Dhe*(Ti_he-Ti)
             + 0.66667*PaD/temp_Ni
             + 0.66667*Pauxi_grid/(Ni+Ni_t+Ni_he)
             ;
    if (auxiliary_heating_i) ddt(Ti) += 0.66667*Pauxi/temp_Ni;
    if (t>=0 && int(t)%pellet_timestep>=0 && int(t)%pellet_timestep<=timestep) ddt(Ti) -= Ti*S_pellet/temp_Ni/2/timestep;
    if (t>=0) ddt(Ti) += Ti*(SnD)/temp_Ni;
//************************************
// Ion Density of T    ---Ni_t---
//************************************
    ddt(Ni_t) =   Diff_psi * Delp2(Ni_t) + DDX(Ni_t) * DDX(Diff_psi) * mesh->g11
                - Ni_t * DDX(V_pinch.x) * sqrt(mesh->g11) - V_pinch.x * DDX(Ni_t) * sqrt(mesh->g11)
                + Diff_psi_para * Grad2_par2(Ni_t)
                + SnT//no atom case
                ;
    if (burning) ddt(Ni_t) += -S_fu;
    if (t>=0 && int(t)%pellet_timestep>=0 && int(t)%pellet_timestep<=timestep) ddt(Ni_t) += S_pellet/2/timestep;
    if (t>=0) ddt(Ni_t) += -SnT;
//    if (fueling_atomT1) ddt(Ni_t) += SnT1;
//******************************************
// Ion Temperature of T         ---Ti_t---
//*******************************************
    ddt(Ti_t) =   0.66667*kappa_psi_para*Grad2_par2(Ti_t)
                - Ti_t*(SnT)/temp_Ni_t // source term
                - 0.66667*Ti_t*(-DDX(Ni_t) * DDX(Diff_psi) * mesh->g11/temp_Ni_t-Diff_psi * Delp2(Ni_t)/temp_Ni_t
                + Diff_psi*DDX(Ni_t)*DDX(Ni_t)*mesh->g11/temp_Ni_t/temp_Ni_t + DDX(V_pinch.x)*sqrt(mesh->g11))  //0.6667*T*Gard(u)
                + 0.66667*kappa_psi*Delp2(Ti_t) + 0.66667*DDX(Ti_t)*DDX(kappa_psi)*mesh->g11
                + 0.66667*kappa_psi*DDX(Ni_t)*DDX(Ti_t)*mesh->g11/temp_Ni_t // conductive term
                + Diff_psi*DDX(Ni_t)*DDX(Ti_t)*mesh->g11/temp_Ni_t - V_pinch.x*DDX(Ti_t)*sqrt(mesh->g11)
                + col_Te*(Te-Ti_t)
                + col_DT*temp_Ni/temp_Ni_t*(Ti-Ti_t)
                + col_The*(Ti_he-Ti_t)
                + 0.66667*PaT/temp_Ni_t
                + 0.66667*Pauxi_grid/(Ni+Ni_t+Ni_he)
                ;
    if (auxiliary_heating_i_t) ddt(Ti_t) += 0.66667*Pauxi_t/temp_Ni_t;
    if (t>=0 && int(t)%pellet_timestep>=0 && int(t)%pellet_timestep<=timestep) ddt(Ti_t) -= Ti_t*S_pellet/temp_Ni_t/2/timestep;
    if (t>=0) ddt(Ti_t) += Ti_t*(SnT)/temp_Ni_t;
//************************************
// Ion Density of He    ---Ni_he---
//************************************
    // output.write("Now updata Ni_he \n");
//output.write("\t Ni_he %e -> %e \n", min(Ni_he), max(Ni_he));

    ddt(Ni_he) = 0.0;

    if (burning) { ddt(Ni_he) += Diff_psi * Delp2(Ni_he) + DDX(Ni_he) * DDX(Diff_psi) * mesh->g11
                                  + N_alpha / tau_slowing_down
                                  - Ni_he * DDX(V_pinch.x) * sqrt(mesh->g11) - V_pinch.x * DDX(Ni_he) * sqrt(mesh->g11)
                                  + Diff_psi_para * Grad2_par2(Ni_he);
                       }
//******************************************
// Ion Temperature of He         ---Ti_he---
//*******************************************
    ddt(Ti_he) =   0.66667*kappa_psi_para*Grad2_par2(Ti_he)
              // - Ti_he*(N_alpha / tau_slowing_down)/temp_Ni_he // source term
                 - 0.66667*Ti_he*(-DDX(Ni_he) * DDX(Diff_psi) * mesh->g11/temp_Ni_he-Diff_psi * Delp2(Ni_he)/temp_Ni_he
                 + Diff_psi*DDX(Ni_he)*DDX(Ni_he)*mesh->g11/temp_Ni_he/temp_Ni_he + DDX(V_pinch.x)*sqrt(mesh->g11))  //0.6667*T*Gard(u)
                 + 0.66667*kappa_psi*Delp2(Ti_he) + 0.66667*DDX(Ti_he)*DDX(kappa_psi)*mesh->g11
                 + 0.66667*kappa_psi*DDX(Ni_he)*DDX(Ti_he)*mesh->g11/temp_Ni_he // conductive term
                 + Diff_psi*DDX(Ni_he)*DDX(Ti_he)*mesh->g11/temp_Ni_he - V_pinch.x*DDX(Ti_he)*sqrt(mesh->g11)
                 + col_Dhe*temp_Ni/temp_Ni_he*(Ti-Ti_he)
                 + col_The*temp_Ni_t/temp_Ni_he*(Ti_t-Ti_he)
                 + col_hee*(Te-Ti_he)
                 + 0.66667*Pahe/temp_Ni_he
                 + 0.66667*Pauxi_grid/(Ni+Ni_t+Ni_he)
                 ;
    if(auxiliary_heating_he) ddt(Ti_he) += 0.66667*Pauxhe/temp_Ni_he;
//************************************
    // ELECTRON TEMPERATURE   ---Te---  temp_Ne=temp_Ni+temp_Ni_t+2*temp_Ni_he+2*N_alpha
//************************************
    ddt(Te) =   0.66667*kappa_psi_para*Grad2_par2(Te)
              - Te*(SnT+SnD)/temp_Ne // source term
              - 0.66667*Te*(-DDX(Ne) * DDX(Diff_psi) * mesh->g11/temp_Ne-Diff_psi * Delp2(Ne)/temp_Ne
              + Diff_psi*DDX(Ne)*DDX(Ne)*mesh->g11/temp_Ne/temp_Ne + DDX(V_pinch.x)*sqrt(mesh->g11))  //0.6667*T*Gard(u)
              + 0.66667*kappa_psi*Delp2(Te) + 0.66667*DDX(Te)*DDX(kappa_psi)*mesh->g11
              + 0.66667*kappa_psi*DDX(Ne)*DDX(Te)*mesh->g11/temp_Ne // conductive term
              + Diff_psi*DDX(Ne)*DDX(Te)*mesh->g11/temp_Ne - V_pinch.x*DDX(Te)*sqrt(mesh->g11)
              + col_De*temp_Ni/temp_Ne*(Ti-Te)
              + col_Te*temp_Ni_t/temp_Ne*(Ti_t-Te)
              + col_hee*temp_Ni_he/temp_Ne*(Ti_he-Te)
              + 0.66667*Pae/temp_Ne
              + 0.66667*Pauxe_grid/temp_Ne
              + 0.66667*Prfe_grid/temp_Ne;
              - 0.66667*Prad_grid/temp_Ne;

    if (auxiliary_heating_e) ddt(Te) += 0.66667*Pauxe/temp_Ne;
    if (auxiliary_rfheating_e) ddt(Te) += 0.66667*Prfe/temp_Ne;
    if (t>=0 && int(t)%pellet_timestep>=0 && int(t)%pellet_timestep<=timestep) ddt(Te) -= S_pellet*Te/temp_Ne/timestep;
    if (t>=0) ddt(Te) += Te*(SnT+SnD)/temp_Ne;
//************************************
    // Alpha density   ---N_alpha---   Its velocity is slowind down distribution now.
//************************************
    ddt(N_alpha) =  0.0;

    if (burning) {
        ddt(N_alpha) += S_fu
                        + Diff_alpha * Delp2(N_alpha) + DDX(N_alpha) * DDX(Diff_alpha) * mesh->g11
                        - N_alpha / tau_slowing_down
                        - N_alpha * DDX(V_pinch_alpha) * sqrt(mesh->g11)
                        -V_pinch_alpha * DDX(N_alpha) * sqrt(mesh->g11)
                        + Diff_psi_para * Grad2_par2(N_alpha);
    }
//********************************************
// For Test Only linws
//********************************************
//    Test1 = DDX_Te / temp_Te;
//
//    Test2 = V_pinch_alpha;
//
//    Test3 = DDX(V_pinch_alpha);
//
//    Test4 = Diff_alpha * Delp2(N_alpha) + DDX(N_alpha) * DDX(Diff_alpha) * mesh->g11;


    return (0);
}


const Field3D field_larger(const Field3D &f, const BoutReal limit) {

    Field3D result;

   result.allocate();

//  #pragma omp parallel for

    for (int jx = 0; jx < mesh->ngx; jx++)

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < mesh->ngz; jz++) {

                if (f[jx][jy][jz] >= limit)

                    result[jx][jy][jz] = f[jx][jy][jz];

                else

                    result[jx][jy][jz] = limit;

            }

    mesh->communicate(result);

    return (result);

}


const BoutReal average(const Field2D &f, const Field3D &k) {

    BoutReal result,S,result2,S2;
    result=0.0;
    S=0;
//    output.write("mesh->ngy is %d\n",mesh->ngy);
    //  result.allocate();
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//  #pragma omp parallel for

    for (int jx = 0; jx <  mesh->ngx; jx++)

        for (int jy = 2; jy <  mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {
                  S=S+dx[jx][jy]*dy[jx][jy]*f[jx][jy];
                  result=result+k[jx][jy][jz]*dx[jx][jy]*dy[jx][jy]*f[jx][jy];
            }

    MPI_Allreduce(&result,&result2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&S,&S2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return (result2/S2);

}

const BoutReal total(const Field2D &f, const Field3D &k) {

    BoutReal result,result2;
    result=0.0;
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // result.allocate();

//  #pragma omp parallel for

    for (int jx = 0; jx < mesh->ngx; jx++)

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

                result=result+k[jx][jy][jz]*dx[jx][jy]*dy[jx][jy]*f[jx][jy];

            }

    MPI_Allreduce(&result,&result2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    return (result2);

}

const Field3D integral_fe(const BoutReal down, const BoutReal up, int n)
{
    BoutReal h;
    Field3D s=0;
    int i;

    h=(up-down)/n;

    s=pow(down,1.5)/(pow(down,1.5)+(W_c^1.5))*h;

    for(i=1;i<n;i++)
    {
        s=s+pow(down+i*h,1.5)/(pow(down+i*h,1.5)+(W_c^1.5))*h;
    }

    return s;
}

const Field3D integral_fi(const BoutReal down, const BoutReal up, int n)
{
    BoutReal h;
    Field3D s=0;
    int i;

    h=(up-down)/n;

    s=(W_c^1.5)/(pow(down,1.5)+(W_c^1.5))*h;

    for(i=1;i<n;i++)
    {
        s=s+(W_c^1.5)/(pow(down+i*h,1.5)+(W_c^1.5))*h;
    }

    return s;
}

//const Field3D Pellet_rp_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp,const int mid,\
//                             const Field3D &Ne,const Field3D &Te, const Field3D &R)
const Field3D Pellet_rp_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp,const int mid)
{
    Field3D p=0;


     for (int jx = 1; jx < mesh->ngx; jx++)   // jx=0 is initial index

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {
//                p[jx][jy][jz] = -1.12e16*(pow(Ne[jx][jy][jz]*1e19*1e-6,0.333))*(pow(Te[jx][jy][jz]*10,1.64))*1.666\
//                *pow(Mi_p,-0.333)*(Rxy[jx][mid]*100*Lbar-Rxy[jx-1][mid]*100*Lbar)/vp/4/PI/np;
                  p[jx][jy][jz] = -1.12e16*(pow(Ne[jx][jy][jz]*1e19*1e-6,0.333))*(pow(Te[jx][jy][jz]*10,1.64))*1.666\
                  *pow(Mi_p,-0.333)*(Rxy[jx][jy]*100*Lbar-Rxy[jx-1][jy]*100*Lbar)/vp/4/PI/np;
            }

    return p;
}

const Field3D Pellet_dN_ablate_rate(const BoutReal np,const BoutReal Mi_p,const BoutReal vp)

{
    Field3D p=0;


     for (int jx = 1; jx < mesh->ngx; jx++)   // jx=0 is initial index

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {
                  p[jx][jy][jz] = 1.12e16*(pow(Ne[jx][jy][jz]*1e19*1e-6,0.333))*(pow(Te[jx][jy][jz]*10,1.64))\
                  *pow(Mi_p,-0.333)*pow(pellet_rp_profile[jx][jy][jz],1.3333)/1e19;
            }

    return p;
}

const Field3D Pellet_rp_profile(const Field3D &rate,const BoutReal rp0)

{
    Field3D rp_profile =0;
    for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

                rp_profile[mesh->ngx-1][jy][jz] = rp0;

            }


    for (int jx = mesh->ngx-2; jx > 0; jx--)   // jx=0 is initial index

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

                if((pow(rp_profile[jx+1][jy][jz],1.66666)+rate[jx][jy][jz])<0)
                    {
                     rp_profile[jx][jy][jz] = 0;
                    }
                else{
                rp_profile[jx][jy][jz] = pow(pow(rp_profile[jx+1][jy][jz],1.66666)+rate[jx][jy][jz],0.6);
                 }
            }

    return rp_profile;
}

const Field3D Pellet_dN_dr(const BoutReal vp)

{
    Field3D p=0;

    for (int jx = 1; jx < mesh->ngx; jx++)   // jx=0 is initial index

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

                  p[jx][jy][jz] = pellet_dN_rate[jx][jy][jz]/vp;
            }

    return p;
}

const Field3D Pellet_displacement()

{   
    int rank;
    Field3D beta;
    beta= pressure*4*PI*1e-7/B_98/B_98*1000;
    Field3D d=0.0;
    BoutReal k;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    BoutReal k2;

   for (int jx = 1; jx < mesh->ngx; jx++)

        { for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {
                       if(rank==(mid/(mesh->ngy-4)))
                      {

                      k=pellet_dN_dr[jx][mid%(mesh->ngy-4)+2][jz]*(Rxy[jx][mid%(mesh->ngy-4)+2]*100*Lbar
                              -Rxy[jx-1][mid%(mesh->ngy-4)+2]*100*Lbar);
		     if(k==0)
                           { k2=0;
                                }
		      else{
/*		     
                      k2=q_profile[jx][mid%(mesh->ngy-4)+2]*q_profile[jx][mid%(mesh->ngy-4)+2]
			*beta[jx][mid%(mesh->ngy-4)+2][jz]/2*R_98/a_radius
                        *(1+q_profile[jx][mid%(mesh->ngy-4)+2]*L_c[jx][mid%(mesh->ngy-4)+2][jz]/a_radius/100)
			*k*100/L_c[jx][mid%(mesh->ngy-4)+2][jz]/a_radius/a_radius
			/(Ne[jx][mid%(mesh->ngy-4)+2][jz]+k*1e6/L_c[jx][mid%(mesh->ngy-4)+2][jz]
							     /r_cloud[jx][mid%(mesh->ngy-4)+2][jz]
							     /r_cloud[jx][mid%(mesh->ngy-4)+2][jz]);
*/
                      k2=q_profile[jx][mid%(mesh->ngy-4)+2]*beta[jx][mid%(mesh->ngy-4)+2][jz]*B_98
			*R_98/a_radius
			*(1+q_profile[jx][mid%(mesh->ngy-4)+2]*L_c[jx][mid%(mesh->ngy-4)+2][jz]/a_radius/100)
			*k*100/L_c[jx][mid%(mesh->ngy-4)+2][jz]
			/(Ne[jx][mid%(mesh->ngy-4)+2][jz]+k*1e6/L_c[jx][mid%(mesh->ngy-4)+2][jz]
                                                             /r_cloud[jx][mid%(mesh->ngy-4)+2][jz]
                                                             /r_cloud[jx][mid%(mesh->ngy-4)+2][jz]);
		      }

		      }
	    }
		       MPI_Bcast(&k2,1,MPI_DOUBLE,mid/(mesh->ngy-4),MPI_COMM_WORLD);
              
          
           for (int jy = 2; jy < mesh->ngy-2; jy++)

              for (int jz = 0; jz < NZ; jz++) {
	    

		      d[jx][jy][jz]=k2;
	    }
     }
   return d;
}

const Field3D Add_pellet_displacement()

{
	Field3D k=0;
	BoutReal f;
        for (int jy = 2; jy < mesh->ngy-2; jy++)

              for (int jz = 0; jz < NZ; jz++) {

                     for (int jx = mesh->ngx-2; jx > 0; jx--)
		     {
			     k[jx][jy][jz]=pellet_displacement[jx][jy][jz]+k[jx+1][jy][jz];
		     }
                     
            }
	return k;

}

const Field3D Pellet_dN_profile()

{
    int rank;
    Field3D p=0;
    BoutReal k;
    BoutReal V=0;
    BoutReal dV=0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int jx = 1; jx < mesh->ngx; jx++)

        { for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

                  if(rank==(mid/(mesh->ngy-4)))
                      {
                      k=pellet_dN_dr[jx][mid%(mesh->ngy-4)+2][jz]*(Rxy[jx][mid%(mesh->ngy-4)+2]*100*Lbar
                              -Rxy[jx-1][mid%(mesh->ngy-4)+2]*100*Lbar);
                      }
              }
            MPI_Bcast(&k,1,MPI_DOUBLE,mid/(mesh->ngy-4),MPI_COMM_WORLD);

            dV=0;
            V=0;

            for (int jy = 2; jy < mesh->ngy-2; jy++)

                  for (int jz = 0; jz < NZ; jz++) {

                  dV=dV+dx[jx][jy]*dy[jx][jy]*mesh->J[jx][jy]*2*PI;

            }
            MPI_Allreduce(&dV,&V,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

            for (int jy = 2; jy < mesh->ngy-2; jy++)

                  for (int jz = 0; jz < NZ; jz++) {

                  p[jx][jy][jz]=k/V;

            }
        }

//
//
//    for (int jx = 1; jx < mesh->ngx-1; jx++)
//        {
//         l=0;
//         dl=0;
//         for (int jy = 2; jy < mesh->ngy-2; jy++)
//
//            for (int jz = 0; jz < NZ; jz++) {
//
//                dl=sqrt(pow(Rxy[jx][jy+1]-Rxy[jx][jy],2)+pow(Zxy[jx][jy+1]-Zxy[jx][jy],2))*Lbar+dl;
//            }
//
//         MPI_Allreduce(&dl,&l,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//
//
//         for (int jy = 2; jy < mesh->ngy-2; jy++)
//
//            for (int jz = 0; jz < NZ; jz++) {
//
//                p[jx][jy][jz] = p[jx][jy][jz]/l/100;
//            }
//
//        }
    return p;
}
//const void test()
//{
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    output.write("rank %d,Rxy[64][32] is %f\n",rank,Rxy[64][32]);
//}
const BoutReal find_rp_deposit_depth()

{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    BoutReal a=0;
    BoutReal b=0;
    BoutReal l=0;
    int jx1=mesh->ngx-1;
    int jx2=0;

    for (int jx = mesh->ngx-1; jx > 0; jx--)   // jx=0 is initial index

        for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

               if(rank==(mid/(mesh->ngy-4)) && pellet_rp_profile[jx][mid%(mesh->ngy-4)+2][jz] == 0)
                      {
                   jx1=jx;
                   if(jx1>jx2) jx2=jx1;
//                      a= Rxy[jx][mid%(mesh->ngy-4)+2];
 //                     if(a>b) b=a;
 //                     if(a>b) jx1=jx;
//                      l= Rxy[mesh->ngx-1][mid%(mesh->ngy-4)+2]-b;
//                      output.write("Rxy[mesh->ngx2] is %f,b is %f,l is %f\n",Lbar*Rxy[mesh->ngx-2][mid%(mesh->ngy-4)+2],Lbar*b,l*Lbar);
                      }
            }
//    output.write("jx1 is %d\n",jx1);
//            output.write("jx2 is %d\n",jx2);
      if(rank==(mid/(mesh->ngy-4))) l=Rxy[mesh->ngx-1][mid%(mesh->ngy-4)+2]*Lbar-Rxy[jx2][mid%(mesh->ngy-4)+2]*Lbar;
//      if(rank==(mid/(mesh->ngy-4)))
//          {
//          output.write("l is %f\n",l);
//          output.write("mid%(mesh->ngy-4) is %d\n",mid%(mesh->ngy-4)+2);
//          output.write("Rxy[mesh->ngx-1][mid%(mesh->ngy-4)+2]*Lbar %f\n",Rxy[mesh->ngx-1][mid%(mesh->ngy-4)+2]*Lbar);
//          output.write("Rxy[jx2][mid%(mesh->ngy-4)+2]*Lbar %f\n",Rxy[jx2][mid%(mesh->ngy-4)+2]*Lbar);
//          output.write("jx2 is %d\n",jx2);
//          output.write("mesh->ngy is %d\n",mesh->ngy);
//          output.write("mid is %d\n",mid);
//          output.write("mesh->ngx-1 is %d\n",mesh->ngx-1);
//          }

    MPI_Bcast(&l,1,MPI_DOUBLE,mid/(mesh->ngy-4),MPI_COMM_WORLD);

    return l;
}

const bool test_ELM_happen()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int ELM = 0;
    bool ELM_happen_mark = false;
    if(rank==(mid/(mesh->ngy-4)) && pressure[ped_index][mid%(mesh->ngy-4)+2][0] > critical_pressure)
        {
        ELM = 1;
        }

    MPI_Bcast(&ELM,1,MPI_INT,mid/(mesh->ngy-4),MPI_COMM_WORLD);

    if(ELM)
        {
        ELM_happen_mark = true;
        }

    return ELM_happen_mark;
}
/*
const void revise_ELM_coef()
{
     for (int jx = 0; jx < mesh->ngx; jx++)

         for (int jy = 2; jy < mesh->ngy-2; jy++)

            for (int jz = 0; jz < NZ; jz++) {

//                  if(psi_norm[jx][jy][jz]> ped_position-ped_width/2)
//                      {
//                      Diff_psi[jx][jy][jz]=50*tbar/Lbar/Lbar;
//                      kappa_psi[jx][jy][jz]=50*tbar/Lbar/Lbar;
//                      }
                    Diff_psi[jx][jy][jz] += Diff_ELM[jx][jy][jz]*tbar/Lbar/Lbar;
                    kappa_psi[jx][jy][jz] += Diff_ELM[jx][jy][jz]*tbar/Lbar/Lbar;
              }

}*/
 
const Field3D Interp_R_new()
{
     int rank;
     Field3D Diff=0;
     Field3D f=0;
     BoutReal newton;
     BoutReal tmp;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==(mid/(mesh->ngy-4)))
        {
        for (int jx = 0; jx < mesh->ngx; jx++)
          {
               Diff[jx][mid%(mesh->ngy-4)+2][0]=Rxy[jx][mid%(mesh->ngy-4)+2];
	  }
        for (int i=0;i<mesh->ngx-1;i++)
	{
		for(int j=mesh->ngx-1;j>i;j--)
		{
		Diff[j][mid%(mesh->ngy-4)+2][0]=(Diff[j][mid%(mesh->ngy-4)+2][0]-Diff[j-1][mid%(mesh->ngy-4)+2][0])
			                       /(psi_norm[j][mid%(mesh->ngy-4)+2][0]-psi_norm[j-1-i][mid%(mesh->ngy-4)+2][0]);
		}
	}
	for (int jx = 0; jx < mesh->ngx; jx++)
         {
		 tmp=1;newton=Diff[jx][mid%(mesh->ngy-4)+2][0];
		 for(int i=0;i<mesh->ngx-1;i++)
		 {
			 tmp=tmp*(pellet_drift_psi[jx][mid%(mesh->ngy-4)+2][0]-psi_norm[i][mid%(mesh->ngy-4)+2][0]);
			 f[jx][mid%(mesh->ngy-4)+2][0]+=tmp*Diff[i+1][mid%(mesh->ngy-4)+2][0];
		 }
	 }
        }

    return f;

}
