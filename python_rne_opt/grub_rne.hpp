/*--------------- MANUAL CHANGES -------------------
	need to update return type in both the function name and function body. This will likely be py::array. Ensure this matches the .c file
	make sure input vars for each function are correct type (e.g. std::vector). Ensure this matches the .c file.
	make sure code doesn't contain invalid c syntax that slipped through
*/

#ifndef GRUB_RNE_HPP
#define GRUB_RNE_HPP

#include <math.h>
#include <vector>
#include <string.h>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <map>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

using namespace Eigen;
using namespace std;


class Grub_rne
{
public:
    // Grub_rne(); //This is a default constructor that doesn't make sense
    // public methods
    Grub_rne(string links, VectorXd discrete_segments, double l);
    void build_robot();
    void update_mat();
    void calc_G();
    void calc_Cqdot();
    void calc_M();
    VectorXd get_xdot(VectorXd uv, VectorXd uv_dot, VectorXd tau);
    void set_uv(VectorXd uv, VectorXd uv_dot);
    VectorXd get_G();
    VectorXd get_Cqdot();
    MatrixXd get_M();

private:

    // All jacobians should be the geometric jacobian found by transforming from spatial jacobian
    // Methods 
    void small_inertia(double h, double mu, double r, double u, double v, double *inertia_h);
    void true_inertia(double h, double mu, double r, double u, double v, double *inertia_h);
    MatrixXd J_dot_udot(double u, double v, double h);
    MatrixXd J_dot_vdot(double u, double v, double h);
    MatrixXd J_end(double u, double v, double h);
    MatrixXd J_dot_udot_com(double u, double v, double h);
    MatrixXd J_dot_vdot_com(double u, double v, double h);
    MatrixXd J_com(double u, double v, double h);
    void get_tau(Vector3d g, VectorXd uvdot, VectorXd uvddot);
    void get_e_rel(int link_num, VectorXd uvdot, VectorXd uvddot);
    void get_c_rel(int link_num, VectorXd uvdot, VectorXd uvddot);
    // MatrixXd get_jacobian(double u, double v, double h);
    void discretize_joint(double h, double mass, double r, int* joint_count, int* link_count, Matrix3d* prev_rot);
    VectorXd get_real_taus();

    // properties
    int _real_joint_num, _total_joint_num, _rigid_num, _link_num;
    int _num_discrete_segments;
    double _l;
    string _links;
    vector<string> _link_type; 
    std::map<std::string, int> _linkmap;
    // types of links are joint4,joint8,jabberwockbase,jabberwockdistal,id
    Vector3d _omega_c_rel, _alpha_c_rel, _omega_e_rel, _alpha_e_rel, _a_c_rel, _a_e_rel, _v_c_rel, _v_e_rel;
    VectorXd _uv, _uv_real, _uvdot_real ,_uvdot, _uvddot, _m_mat, _G, _Cqdot, _qdd, 
            _flexi_link, _h, _report_tau, _discrete_segments;
    MatrixXd _rc_mat, _re_mat, _rot_mat, _rot_ground_mat, _Inertia_mat, 
            _omega_mat, _alpha_mat, _a_c_mat, _a_e_mat, _f_mat,
            _tau_mat, _qdot_mat, _qddot_mat, _M, _K_ps; 
 
};

PYBIND11_MODULE(grub_rne, m) {
    py::class_<Grub_rne>(m, "rne_model")
        .def(py::init<string, VectorXd, double>())
        .def("calc_G", &Grub_rne::calc_G)
        .def("calc_Cqdot", &Grub_rne::calc_Cqdot)
        .def("calc_M", &Grub_rne::calc_M)
        .def("get_xdot",&Grub_rne::get_xdot)
        .def("set_uv",&Grub_rne::set_uv)
        .def("get_G",&Grub_rne::get_G)
        .def("get_Cqdot",&Grub_rne::get_Cqdot)
        .def("get_M",&Grub_rne::get_M);
}

#endif // GRUB_RNE_HPP
