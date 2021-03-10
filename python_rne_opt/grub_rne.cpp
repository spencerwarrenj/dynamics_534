/*------------------------- MANUAL CHANGES ------------------------------
	need to add 'namespace py = pybind11;' after #include statements
	need to update return type in both the function name and function body. This will likely be py::array. Ensure this matches the header file.
	make sure input vars for each function are correct type (e.g. std::vector). These should also match the header file prototypes.
	make sure code doesn't contain invalid c syntax that slipped through
*/

#include <math.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <map>

// #include "bellows_rne.hpp"
#include "grub_rne.hpp"

using namespace Eigen;
// Grub_rne::Grub_rne()
// {
//     std::cout << "[Rectangle] default constructor!!! Don't Do it!" << std::endl;
//     _links = "helloWorld";
// }

// This is the constructor for the grub_rne class
Grub_rne::Grub_rne(string links,VectorXd discrete_segments,double l)
{
    _links = links;
    _discrete_segments = discrete_segments;
    _num_discrete_segments = discrete_segments.size();
    _l = l;
    // std::cout << "l: " << _l << std::endl;
    // std::cout << "number of discrete segments: " << _num_discrete_segments << std::endl;
    // std::cout << "discrete segments: " << _discrete_segments << std::endl;
    // std::cout << "Constructor: link: " << _links << "\n";
    build_robot();
}

void Grub_rne::discretize_joint(double h, double mass, double r, 
                int* joint_count, int* link_count, Matrix3d* prev_rot)
{
    Vector3d rc,re;
    Matrix3d rot, Inertia_tens, Inertia_diag, skewr;
    double *Inertia_h = new double[9];
    double mu, sigma, theta, v_tilde, u_tilde, cphi, sphi, phi, u, v, rho;
    int i;
    int discrete_joint_count = 0;
    VectorXd l_vec, mass_vec;


    l_vec = h*_discrete_segments;
    mass_vec = mass*_discrete_segments;
    for (discrete_joint_count=0;discrete_joint_count<
                    _num_discrete_segments;discrete_joint_count++)
    {
        u = _uv(*joint_count*2);
        v = _uv(*joint_count*2+1);
        phi = sqrt(pow(u,2.)+pow(v,2));
        sphi = sin(phi);
        cphi = cos(phi);
        h = l_vec(discrete_joint_count);
        // std::cout << "link_count: " << *link_count << std::endl;
        
        _h(*link_count) = h;
        // std::cout << "_h: " << _h << std::endl;
        if (phi > pow(10,-6))
        {
            u_tilde = u/phi;
            v_tilde = v/phi;
            sigma = cphi-1;
            theta = atan2(-u,v);
            rho = h/phi;   
            r = rho*(1-cphi);
            rc << (phi-sphi)*v_tilde, (phi-sphi)*-u_tilde, 1-cphi;
            _rc_mat.block<3,1>(0,*link_count) = h/pow(phi,2)*rc;
            re << r*cos(theta), r*sin(theta), rho*sphi;
            rot << sigma*pow(v_tilde,2)+1, -sigma*u_tilde*v_tilde, v_tilde*sphi,
                    -sigma*u_tilde*v_tilde, sigma*pow(u_tilde,2)+1, -u_tilde*sphi,
                    -v_tilde*sphi, u_tilde*sphi, cphi;
        }
        else
        {
            re << 1/2*v*h, -1/2*u*h, h-1/6*h*pow(phi,2);
            rc = re*1/2; // just halved it since I couldn't find anything else should be minimal error
            _rc_mat.block<3,1>(0,*link_count) = rc;
            re << 1/2*v*h, -1/2*u*h, h-1/6*h*pow(phi,2);
            rot << 1-1/2*pow(v,2), 1/2*u*v, v,
                    1/2*u*v, 1-1/2*pow(u,2), -u,
                    -v, u, 1-1/2*pow(phi,2);
        }
        _re_mat.block<3,1>(0,*link_count) = re;
        // std::cout << "rc: " << v << "\n";
        _rot_mat.block<3,3>(0,*link_count*3) = rot;
        
        _m_mat(*link_count) = mass_vec(discrete_joint_count);
        _rot_ground_mat.block<3,3>(0,*link_count*3) = rot**prev_rot;
        *prev_rot = rot**prev_rot;
        mu = mass/h;
        if (abs(phi)<M_PI/24)
        {
            small_inertia(h, mu, r, u, v, Inertia_h);
        }
        else
        {
            true_inertia(h, mu, r, u, v, Inertia_h);
        }
        Inertia_tens << Inertia_h[0],Inertia_h[1],Inertia_h[2],
                        Inertia_h[3],Inertia_h[4],Inertia_h[5],
                        Inertia_h[6],Inertia_h[7],Inertia_h[8];
        _Inertia_mat.block<3,3>(0,*link_count*3) = Inertia_tens;
        _flexi_link(*link_count) = 1; 
        if (discrete_joint_count == 0)
        {
            _report_tau(*link_count) = 1;
        }
        *link_count = *link_count+1; 
        *joint_count = *joint_count+1;
    } 
    // std::cout << "tau_report: " << _report_tau << std::endl;
}

// This function gets all of the important info for the robot
void Grub_rne::build_robot()
{
    Vector3d rc,re;
    Matrix3d rot, Inertia_tens, Inertia_diag, skewr;
    double *Inertia_h = new double[9];
    double r, mass, mu, sigma, theta, v_tilde, u_tilde, cphi, sphi, phi, u, v, h, rho;

    // parse links string
    stringstream s_stream(_links); //create string stream from the string
    while(s_stream.good()) {
      string substr;
      getline(s_stream, substr, ','); //get first string delimited by comma
      _link_type.push_back(substr);
    //   std::cout << "link_type: " << substr << std::endl;
    }

    // std::cout << _link_type[0] << std::endl;

    // create map from joint type to int
    
    _linkmap.insert(std::make_pair("joint8", 1));
    _linkmap.insert(std::make_pair("joint4", 2));
    _linkmap.insert(std::make_pair("jabberwockbase", 3));
    _linkmap.insert(std::make_pair("jabberwockdistal", 4));
    _linkmap.insert(std::make_pair("id", 5));

    _real_joint_num = 0;
    _rigid_num = 0;
    // loop through every joint
    for (int i=0; i<_link_type.size(); i++)
    {
        switch (_linkmap[_link_type[i]])
        {
            case 1 :
                _real_joint_num++;
                break;
            case 2 :
                _real_joint_num++;
                break;
            case 3 :
                _rigid_num++;
                break;
            case 4 :
                _rigid_num++;
                break;
            case 5 :
                _rigid_num++; 
                break;               
            default :
                std::cout << "You need to specify real joint names!!" << std::endl;
                break;
        }
    }
    _total_joint_num = _real_joint_num*_num_discrete_segments;
    _link_num = _total_joint_num+_rigid_num;
    // resize matrices appropriately
    // parameters
    _rc_mat.resize(3,_total_joint_num+_rigid_num);
    _re_mat.resize(3,_total_joint_num+_rigid_num);
    _rot_mat.resize(3,(_total_joint_num+_rigid_num)*3);
    _rot_ground_mat.resize(3,(_total_joint_num+_rigid_num)*3);
    _Inertia_mat.resize(3,(_total_joint_num+_rigid_num)*3);
    _uv.resize(_total_joint_num*2);
    _uvdot.resize(_total_joint_num*2);
    _uvddot.resize(_total_joint_num*2);
    _K_ps.resize(_real_joint_num*2,_real_joint_num*2);
    _m_mat.resize(_rigid_num+_total_joint_num);
    // forward recursion
    _omega_mat.resize(3,_link_num);
    _alpha_mat.resize(3,_link_num);
    _a_c_mat.resize(3,_link_num);
    _a_e_mat.resize(3,_link_num);
    // backward recursion
    _f_mat.resize(3,_link_num+1);
    _tau_mat.resize(3,_link_num+1);
    // _qdot_mat.resize(3,_link_num);
    // _qddot_mat.resize(3,_link_num);
    _flexi_link.resize(_total_joint_num+_rigid_num);
    _report_tau.resize(_total_joint_num+_rigid_num);
    _h.resize(_total_joint_num+_rigid_num);

    _report_tau.setZero();
    _flexi_link.setZero();
    _h.setZero();
    // _uv.setZero();

    VectorXd vec(_real_joint_num*2);
    vec.setOnes();
    vec = vec*.1;
    set_uv(vec,vec);
    // _uv.setOnes();
    // _uv = _uv*.1;
    // _uvdot.setZero();
    _uvddot.setZero();

    int joint_count = 0;
    int link_count = 0;
    Matrix3d prev_rot =  MatrixXd::Identity(3,3);
    // loop through every joint
    for (int i=0; i<_link_type.size(); i++)
    {
        switch (_linkmap[_link_type[i]])
        {
            case 1 : // joint 8
                h = 0.212;
                mass = 2.65273333333333;
                r = 0.25;
                discretize_joint(h,mass,r,&joint_count, &link_count, &prev_rot);
                break;
            case 2 : //joint 4
                h = .202;
                mass = 1.32636666666667;
                r = 0.16;
                discretize_joint(h,mass,r,&joint_count, &link_count, &prev_rot);
                break;
            case 3 : //jabberwockbase
                rc << 0.0, 0.0, 0.115;
                _rc_mat.block<3,1>(0,link_count) = rc;
                re << -.02234224, .02234214, .2819785;
                _re_mat.block<3,1>(0,link_count) = re;
                rot <<  0.8536,  0.1464, -0.5000,
                        0.1464,  0.8536,  0.5000,
                        0.5000, -0.5000, 0.7071;
                _rot_mat.block<3,3>(0,link_count*3) = rot;
                _m_mat(link_count) = 6.534-2.65273333333333;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                Inertia_diag << 0.108, 0.0, 0.0,
                                0.0, 0.108, 0.0,
                                0.0, 0.0, 0.023;
                skewr <<  0.0,-rc(2), rc(1),
                            rc(2), 0.0, -rc(0),
                            -rc(1), rc(0), 0.0;
                Inertia_tens << Inertia_diag-mass*skewr*skewr; // shifting inertia tensor
                _Inertia_mat.block<3,3>(0,link_count*3) = Inertia_tens;
                _flexi_link(link_count) = 0;
                link_count++;
                break;
            case 4 : //jabberwock distal
                rc << -0.07, 0.07, 0.13;
                _rc_mat.block<3,1>(0,link_count) = rc;
                re << -.10307463,.10307473,0.19021958;
                _re_mat.block<3,1>(0,link_count) = re;
                rot <<  0.8536,  0.1464, -0.5000,
                        0.1464,  0.8536,  0.5000,
                        0.5000, -0.5000, 0.7071;
                _rot_mat.block<3,3>(0,link_count*3) = rot;
                _m_mat(link_count) = 4.8-1.32636666666667;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                Inertia_diag << 0.05, 0.0, 0.0,
                                0.0, 0.05, 0.0,
                                0.0, 0.0, 0.017;
                skewr <<  0.0,-rc(2), rc(1),
                            rc(2), 0.0, -rc(0),
                            -rc(1), rc(0), 0.0;
                Inertia_tens << Inertia_diag-mass*skewr*skewr; // shifting inertia tensor
                _Inertia_mat.block<3,3>(0,link_count*3) = Inertia_tens;
                _flexi_link(link_count) = 0;
                link_count++;
                break;
            case 5 : //simple nothings
                rc << 0.0, 0.0, 0;
                _rc_mat.block<3,1>(0,link_count) = rc;
                re << 0.0,0.0,0.0;
                _re_mat.block<3,1>(0,link_count) = re;
                rot <<  1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;
                _rot_mat.block<3,3>(0,link_count*3) = rot;
                _m_mat(link_count) = .1;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                Inertia_diag << 0.001, 0.0, 0.0,
                                0.0, 0.001, 0.0,
                                0.0, 0.0, 0.001;
                skewr <<  0.0,-rc(2), rc(1),
                            rc(2), 0.0, -rc(0),
                            -rc(1), rc(0), 0.0;
                Inertia_tens << Inertia_diag-mass*skewr*skewr; // shifting inertia tensor
                _Inertia_mat.block<3,3>(0,link_count*3) = Inertia_tens;
                _flexi_link(link_count) = 0;
                link_count++;
                break;
            default :
                std::cout << "You need to specify real joint names!!" << std::endl;
                break;
        }

        // std::cout << i << std::endl;
    }

    // std::cout << "Prev_rot: " << prev_rot << std::endl;
    // std::cout << "rc_mat: " << _rc_mat << "\n";
    // std::cout << "re_mat: " << _re_mat << "\n";
    // std::cout << "rot_mat: " << _rot_mat << "\n";
    // std::cout << "rot_ground_mat: " << _rot_ground_mat << "\n";
    // std::cout << "m_mat: " << _m_mat << "\n";
    // std::cout << "Inertia_mat: " << _Inertia_mat << "\n";

}

// This function updates only the changing elements for the robot
void Grub_rne::update_mat()
{
    Vector3d rc,re;
    Matrix3d rot, Inertia_tens, Inertia_diag, skewr;
    double *Inertia_h = new double[9];
    double r, mass, mu, sigma, theta, v_tilde, u_tilde, cphi, sphi, phi, u, v, h, rho;

    int joint_count = 0;
    int link_count = 0;
    Matrix3d prev_rot =  MatrixXd::Identity(3,3);
    // loop through every joint
    for (int i=0; i<_link_type.size(); i++)
    {
        switch (_linkmap[_link_type[i]])
        {
            case 1 : // joint 8
                h = 0.212;
                mass = 2.65273333333333;
                r = 0.25;
                discretize_joint(h,mass,r,&joint_count, &link_count, &prev_rot);
                break;
            case 2 : //joint 4
                h = .202;
                mass = 1.32636666666667;
                r = 0.16;
                discretize_joint(h,mass,r,&joint_count, &link_count, &prev_rot);
                break;
            case 3 : //jabberwockbase
                rot <<  0.8536,  0.1464, -0.5000,
                        0.1464,  0.8536,  0.5000,
                        0.5000, -0.5000, 0.7071;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                link_count++;
                break;
            case 4 : //jabberwock distal
                rot <<  0.8536,  0.1464, -0.5000,
                        0.1464,  0.8536,  0.5000,
                        0.5000, -0.5000, 0.7071;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                link_count++;
                break;
            case 5 : //simple nothings
                rot <<  1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;
                _rot_ground_mat.block<3,3>(0,link_count*3) =  rot*prev_rot;
                prev_rot = rot*prev_rot;
                link_count++;
                break;
            default :
                std::cout << "You need to specify real joint names!!" << std::endl;
                break;
        }

        // std::cout << i << std::endl;
    }

    // std::cout << "rc_mat: " << _rc_mat << "\n";
    // std::cout << "re_mat: " << _re_mat << "\n";
    // std::cout << "rot_mat: " << _rot_mat << "\n";
    // std::cout << "rot_ground_mat: " << _rot_ground_mat << "\n";
    // std::cout << "m_mat: " << _m_mat << "\n";
    // std::cout << "Inertia_mat: " << _Inertia_mat << "\n";

}

VectorXd Grub_rne::get_xdot(VectorXd uv, VectorXd uv_dot,VectorXd tau)
{
    set_uv(uv,uv_dot);
    _K_ps.setIdentity();
    _K_ps = _K_ps*9;
    calc_G();
    calc_Cqdot();
    calc_M();
    MatrixXd minv = _M.inverse();
    // std::cout << "minv: " << minv << std::endl;
    // std::cout << "tau: " << tau << std::endl;
    // std::cout << "K_ps: " << _K_ps << std::endl;
    _qdd = minv*(-_Cqdot-_G+tau+_K_ps*_uv_real);
    

    return _qdd;
}

VectorXd Grub_rne::get_real_taus()
{
    Matrix3d rot;
    VectorXd tau_block(_real_joint_num*2);

    int tau_report_count = 0;

    for (int link_count = 0; link_count < _link_num; link_count++)
    {
    // std::cout << "i: " << i << "\n";
    // std::cout << "set G" << "\n";
        if (_report_tau(link_count)==1)
        {
            rot = _rot_mat.block<3,3>(0,link_count*3);
            
            Vector3d dum = rot*_tau_mat.block<3,1>(0,link_count);
            
            // std::cout << "G: " << _G << "\n";
            tau_block.block<2,1>(link_count*2,0) = dum.block<2,1>(0,0);
            tau_report_count++;
        }
    }
    return tau_block;
}

// This function calculates the G vector from RNE equations
void Grub_rne::calc_G()
{
    // Inputs:
    // rc_mat/re_mat is a 3x6 matrix with rc/re [rc1,rc2,rc3,...]
    // rot_mat is a 3x18 matrix with the rotations from base to tip [R01, R12, R23,...]
    // rot_ground_mat is a 3x18 with the rotations from the ground frame to the link frame [R01, R02, R03,...]
    // m_mat is a 6x1 vector of the masses of each rc_link
    // Inertia_mat is a 3x18 matrix of the inertia tensors about the center of mass
    // of each link [I1, I2, I3, ...]
    // G is a pointer to a 1x6 vector to be filled

    //regular Stuff
    VectorXd uvdot(_total_joint_num*2), uvddot(_total_joint_num*2), tau_block;
    Matrix3d rot;
    // std::cout << "omega_mat: " << omega_mat << "\n";

    //Special Stuff to find G
    uvdot.setZero();
    uvddot.setZero();
    Vector3d g;
    g << 0,0,-9.81;
    _G.resize(2*_real_joint_num,1);

    // std::cout << "uvdot: " << uvdot <<"\n";
    // std::cout << "uvddot: " << uvddot <<"\n";

    get_tau(g,uvdot,uvddot);

    // std::cout << "tau: " << _tau_mat <<"\n";
    // R01,R12,R23,R34,R45,R56
    _G = get_real_taus();
    
    // std::cout << "G: " << _G << "\n";

}

// This function calculates the Cqdot vector from RNE equations
void Grub_rne::calc_Cqdot()
{
    // Inputs:
    // rc_mat/re_mat is a 3x6 matrix with rc/re [rc1,rc2,rc3,...]
    // rot_mat is a 3x18 matrix with the rotations from base to tip [R01, R12, R23,...]
    // rot_ground_mat is a 3x18 with the rotations from the ground frame to the link frame [R01, R02, R03,...]
    // m_mat is a 6x1 vector of the masses of each rc_link
    // Inertia_mat is a 3x18 matrix of the inertia tensors about the center of mass
    // of each link [I1, I2, I3, ...]
    // G is a pointer to a 1x6 vector to be filled

    //regular Stuff
    VectorXd uvdot(_total_joint_num*2), uvddot(_total_joint_num*2);
    Matrix3d rot;
    // std::cout << "omega_mat: " << omega_mat << "\n";

    //Special Stuff to find Cqdot
    uvdot = _uvdot;
    uvddot.setZero();
    Vector3d g;
    g << 0.,0.,0.;
    _Cqdot.resize(2*_real_joint_num,1);

    get_tau(g,uvdot,uvddot);

    // std::cout << "tau: " << _tau_mat <<"\n";
    // R01,R12,R23,R34,R45,R56
    _Cqdot = get_real_taus();
    // std::cout << "Cqdot: " << _Cqdot << "\n";
}

void Grub_rne::calc_M()
{
    // Inputs:
    // rc_mat/re_mat is a 3x6 matrix with rc/re [rc1,rc2,rc3,...]
    // rot_mat is a 3x18 matrix with the rotations from base to tip [R01, R12, R23,...]
    // rot_ground_mat is a 3x18 with the rotations from the ground frame to the link frame [R01, R02, R03,...]
    // m_mat is a 6x1 vector of the masses of each rc_link
    // Inertia_mat is a 3x18 matrix of the inertia tensors about the center of mass
    // of each link [I1, I2, I3, ...]
    // G is a pointer to a 1x6 vector to be filled

    //regular Stuff
    VectorXd uvdot(_total_joint_num*2), uvddot(_total_joint_num*2), num_ones(_num_discrete_segments);
    Matrix3d rot, rott;

    //Special Stuff to find M
    num_ones.setOnes();
    uvdot.setZero();
    uvddot.setZero();
    Vector3d g;
    g << 0.,0.,0.;
    _M.resize(_real_joint_num*2,_real_joint_num*2);

    
    for (int j = 0; j < _real_joint_num*2; j++)
    {
        // set u_ddot to 1
        // std::cout << "_discrete_segments: " << _discrete_segments << "\n";
        uvddot.setZero();
        for (int i = 0; i < _real_joint_num; i++)
        {
            for (int k=(j%2); k < _num_discrete_segments*2; k=k+2)
            {
                uvddot(i*_num_discrete_segments+k) = _discrete_segments(k/2);
            }
        }
        
        // std::cout << "uvddot: " << uvddot << "\n";
        get_tau(g,uvdot,uvddot);
        // std::cout << "tau: " << _tau_mat <<"\n";
        // R01,R12,R23,R34,R45,R56
        _M.col(j) = get_real_taus();
    }
}

// This function runs through the forward and backward recursion to get dynamics
void Grub_rne::get_tau(Vector3d g, VectorXd uvdot, VectorXd uvddot)
{
    // std::cout << "got here " <<"\n";
    Vector3d prev_alpha, prev_omega, prev_a_e, rc, re, omega,
            alpha, a_c, a_e, f, tau, prev_f, prev_tau, qdot, qddot;
    Matrix3d rot, rott, rotf, rotb, rotground, Inertia_tensor, 
            rotgroundt;
    
    double phi, sphi, cphi, u_tilde, v_tilde, sigma, theta, rho, r, m;

    _omega_mat.setZero();
    _alpha_mat.setZero();
    _a_c_mat.setZero();
    _a_e_mat.setZero();
    _f_mat.setZero();
    _tau_mat.setZero();
    // _qdot_mat.setZero();
    // _qddot_mat.setZero();

    // std::cout << "qdot_mat: " << _qdot_mat << "\n";
    // std::cout << "qddot_mat: " << _qddot_mat << "\n";
   //forward recursion
    for (int i = 0; i < _link_num; ++i)
    {
        //get params
        if (i == 0)
        {
            prev_alpha.setZero();
            prev_omega.setZero();
            prev_a_e.setZero();
        }
        else
        {
            // prev_alpha = _alpha_mat.block<3,1>(0,i-1);
            // prev_omega = _omega_mat.block<3,1>(0,i-1);
            prev_a_e = _a_e_mat.block<3,1>(0,i-1);
        }
        // qdot = _qdot_mat.block<3,1>(0,i);
        // qddot = _qddot_mat.block<3,1>(0,i);
        // std::cout << "qdot: " << qdot << "\n";
        // std::cout << "qddot: " << qddot << "\n";
        
        rc = _rc_mat.block<3,1>(0,i);
        re = _re_mat.block<3,1>(0,i);
        rot = _rot_mat.block<3,3>(0,i*3);
        // std::cout << "rot_mat: " << _rot_mat << "\n";
        // std::cout << "Rot" << rot << "\n";
        rott = rot.transpose();
        rot = rott;
        // std::cout << "Rot" << rot << "\n";
        // std::cout << "Start Forward recursion" << "\n";
        
        // forward recursion for kinematics
        // set omega_c_rel, alpha_c_rel, omega_e_rel, alpha_e_rel, 
        //              a_c_rel, v_c_rel, a_e_rel, v_e_rel
        get_e_rel(i,uvdot,uvddot);
        get_c_rel(i,uvdot,uvddot);
        // std::cout << "omega_rel: " << _omega_c_rel << "\n";
        // alpha_rel = get_alpha_rel(i,uvdot,uvddot);
        // std::cout << "alpha_rel: " << _alpha_c_rel << "\n";
        // std::cout << "prev_alpha: " << prev_alpha << "\n";
        // std::cout << "prev_a_e: " << prev_a_e << "\n";
        // std::cout << "a_e_rel: " << _a_e_rel << "\n";
        // std::cout << "re: " << re << "\n";


        a_c = rot*(prev_a_e+_a_c_rel+prev_alpha.cross(rc)+
                prev_omega.cross(prev_omega.cross(rc))+
                2*prev_omega.cross(_v_c_rel));
        a_e = rot*(prev_a_e+_a_e_rel+prev_alpha.cross(re)+
                prev_omega.cross(prev_omega.cross(re))+
                2*prev_omega.cross(_v_e_rel));
        omega = rot*prev_omega+rot*_omega_c_rel;
        alpha = rot*prev_alpha+rot*_alpha_c_rel;

        // std::cout << "End Forward recursion" << "\n";
        // std::cout << "omega: " << omega << "\n";
        // std::cout << "alpha:  " << alpha << "\n";
        // std::cout << "ac:  " << a_c << "\n";
        // std::cout << "ae:  " << a_e << "\n";

        // set kinematics
        _omega_mat.block<3,1>(0,i) = omega;
        _alpha_mat.block<3,1>(0,i) = alpha;
        _a_c_mat.block<3,1>(0,i) = a_c;
        _a_e_mat.block<3,1>(0,i) = a_e;

        prev_omega = rot*prev_omega+rot*_omega_e_rel;
        prev_alpha = rot*prev_alpha+rot*_alpha_e_rel;
        // prev_a_e = a_e;

        // std::cout << "omega_mat: " << _omega_mat << "\n";
        // std::cout << "alpha_mat: " << _alpha_mat << "\n";
        // std::cout << "a_c_mat: " << _a_c_mat << "\n";
        // std::cout << "a_e_mat: " << _a_e_mat << "\n";
    }
    
    // std::cout << "a_e_mat: " << _a_e_mat << "\n";
    // std::cout << "a_c_mat: " << _a_c_mat << "\n";
    
    // backward recursion
    for (int i = _link_num-1; i >= 0; i--)
    {
        // get params
        omega = _omega_mat.block<3,1>(0,i);
        alpha = _alpha_mat.block<3,1>(0,i);
        a_c = _a_c_mat.block<3,1>(0,i);
        // std::cout << "omega: " << omega << "\n";
        // std::cout << "alpha: " << alpha << "\n";
        // std::cout << "a_c: " << a_c << "\n"; 
        // std::cout << "rot_mat: " << _rot_mat << "\n";
        rot = _rot_mat.block<3,3>(0,i*3);
        rott = rot.transpose();
        rotb = rott;
        // std::cout << "rotb: " << rotb << "\n";
        if (i != _link_num-1)
        {
        rotf = _rot_mat.block<3,3>(0,(i+1)*3);
        }
        else
        {
        rotf = MatrixXd::Identity(3,3);
        }
        // std::cout << "rotf: " << rotf << "\n";
        rotground = _rot_ground_mat.block<3,3>(0,i*3);
        rotgroundt = rotground.transpose();
        rotground = rotgroundt;
        m = _m_mat(i);
        // std::cout << "m: " << m << std::endl << "m_mat: " << _m_mat << "\n";
    //  std::cout << "rc_mat: " << rc_mat << "\n";
        rc = _rc_mat.block<3,1>(0,i);
        re = _re_mat.block<3,1>(0,i);
    //  std::cout << "rc" << rc << "\n";
    //  std::cout << "re_mat" << re_mat << "\n";
    //  std::cout << "re: " << re << "\n";
        Inertia_tensor = _Inertia_mat.block<3,3>(0,i*3);
        prev_f = _f_mat.block<3,1>(0,i+1);
        prev_tau = _tau_mat.block<3,1>(0,i+1);

        // std::cout << "Start Backward Recursion" << "\n";
        
        // get forward kinetics
        // std::cout << "prev_f: " << prev_f <<"\n";
        // std::cout << "g: " << g <<"\n";
        // std::cout << "a_c: " << a_c <<"\n";
        f = rotf*prev_f-m*rotground*g+m*a_c;
        // std::cout << "f: " << f << std::endl;
        tau = rotf*prev_tau-f.cross(rotb*rc)+(rotf*prev_f).cross(rotb*(rc-re))+
                        Inertia_tensor*alpha+omega.cross(Inertia_tensor*omega);
        // std::cout << "tau: " << tau << std::endl;
        // set dynamics
        _f_mat.block<3,1>(0,i) = f;
        // std::cout << "f_mat: " << _f_mat << "\n";
        _tau_mat.block<3,1>(0,i) = tau;
    }
        // std::cout << "tau_mat: " << _tau_mat << "\n";
}

// This is the function to the small angle approximation of the inertia to avoid numerical problems
void Grub_rne::small_inertia(double h, double mu, double r, double u, double v, double *inertia_h) {
   const double x0 = pow(h, 3);
   const double x1 = 0.041666666666666671*mu*x0;
   const double x2 = pow(r, 2);
   const double x3 = mu*x2;
   const double x4 = 0.125*x3;
   const double x5 = h*x4;
   const double x6 = pow(v, 3);
   const double x7 = 0.083333333333333329*mu;
   const double x8 = pow(u, 2);
   const double x9 = mu*x8;
   const double x10 = 0.083333333333333329*x9;
   const double x11 = 0.16666666666666666*x0;
   const double x12 = 0.125*x2;
   const double x13 = 0.25*h;
   const double x14 = v*x1 - v*x5 - x11*(v*x10 + x6*x7) + x13*(v*x12*x9 + x4*x6);
   const double x15 = v*x14;
   const double x16 = 0.5*u;
   const double x17 = pow(v, 2);
   const double x18 = 1 - 0.5*x17;
   const double x19 = pow(u, 3);
   const double x20 = x17*x7;
   const double x21 = mu*x17;
   const double x22 = -u*x1 + u*x5 - x11*(-u*x20 - x19*x7) + x13*(-u*x12*x21 - x19*x4);
   const double x23 = 0.25*x9;
   const double x24 = 0.25*x21;
   const double x25 = 0.20000000000000001*x0;
   const double x26 = x2*x23;
   const double x27 = pow(u, 4);
   const double x28 = x27*x3;
   const double x29 = pow(v, 4);
   const double x30 = 0.20000000000000001*h;
   const double x31 = 0.5*h*x3 + 0.33333333333333331*h*(-x2*x24 - x26) - x25*(-x23 - x24) + x30*(x17*x26 + 0.125*x28 + x29*x4);
   const double x32 = -v*x31 + x15*x16 + x18*x22;
   const double x33 = u*v;
   const double x34 = h*x2;
   const double x35 = 0.0625*x3;
   const double x36 = 0.050000000000000003*mu*x0*x33 + x30*(-u*x35*x6 - v*x19*x35) - x33*x34*x7;
   const double x37 = x17*x9;
   const double x38 = 0.0625*x2*x37;
   const double x39 = 0.083333333333333315*mu;
   const double x40 = 0.16666666666666666*mu;
   const double x41 = 0.027777777777777776*mu;
   const double x42 = x0*x39 - 0.25*x0*(-x17*x40 - x40*x8) - 0.14285714285714285*x0*(-x27*x41 - x29*x41 - 0.055555555555555552*x37) + x13*x3;
   const double x43 = x20*x34 - x25*(0.33333333333333331*x21 + x39*x8) + x30*(x29*x35 + x38) + x42;
   const double x44 = u*x43;
   const double x45 = 0.5*v;
   const double x46 = -x15 + x18*x36 + x44*x45;
   const double x47 = u*x46;
   const double x48 = u*x36;
   const double x49 = x45*x48;
   const double x50 = v*x22;
   const double x51 = x10*x34 - x25*(x17*x39 + 0.33333333333333331*x9) + x30*(0.0625*x28 + x38) + x42;
   const double x52 = x18*x51 + x49 - x50;
   const double x53 = v*x52;
   const double x54 = -0.5*x8;
   const double x55 = x54 + 1;
   const double x56 = x18 + x54;
   const double x57 = u*x31 + x14*x55 + x16*x50;
   const double x58 = u*x14;
   const double x59 = x43*x55 + x49 + x58;
   const double x60 = u*x59;
   const double x61 = v*x51;
   const double x62 = u*x22 + x16*x61 + x36*x55;
   const double x63 = v*x62;
   const double x64 = x31*x56 + x50 - x58;
   const double x65 = v*x36 + x14*x56 - x44;
   const double x66 = u*x65;
   const double x67 = x22*x56 - x48 + x61;
   const double x68 = v*x67;

   inertia_h[0] = -v*x32 + x18*x52 + x45*x47;
   inertia_h[1] = u*x32 + x16*x53 + x46*x55;
   inertia_h[2] = x32*x56 - x47 + x53;
   inertia_h[3] = -v*x57 + x18*x62 + x45*x60;
   inertia_h[4] = u*x57 + x16*x63 + x55*x59;
   inertia_h[5] = x56*x57 - x60 + x63;
   inertia_h[6] = -v*x64 + x18*x67 + x45*x66;
   inertia_h[7] = u*x64 + x16*x68 + x55*x65;
   inertia_h[8] = x56*x64 - x66 + x68;

}

// This is the True inertial term. This method runs into problems with small phis close to zero
void Grub_rne::true_inertia(double h, double mu, double r, double u, double v, double *inertia_h) {
   const double x0 = pow(u, 2);
   const double x1 = v*x0;
   const double x2 = pow(v, 2);
   const double x3 = x0 + x2;
   const double x4 = sqrt(x3);
   const double x5 = 2*x4;
   const double x6 = x0*x5;
   const double x7 = x2*x5;
   const double x8 = 1.0/(x6 + x7);
   const double x9 = 1.0/x4;
   const double x10 = (1.0/2.0)*h;
   const double x11 = mu*pow(r, 2);
   const double x12 = x10*x11;
   const double x13 = x12*x9;
   const double x14 = x13*x8;
   const double x15 = pow(h, 3);
   const double x16 = pow(u, 4);
   const double x17 = x16*x5;
   const double x18 = pow(v, 4);
   const double x19 = x18*x5;
   const double x20 = 4*x4;
   const double x21 = x0*x2;
   const double x22 = 1.0/(x17 + x19 + x20*x21);
   const double x23 = x15*x22;
   const double x24 = sin(x4);
   const double x25 = 2*x24;
   const double x26 = mu*x23*x25;
   const double x27 = cos(x4);
   const double x28 = 2*x27;
   const double x29 = h*x8;
   const double x30 = x0*x29;
   const double x31 = pow(x24, 2);
   const double x32 = x2*x31;
   const double x33 = x29*x32;
   const double x34 = (1.0/4.0)*x11;
   const double x35 = x34*x9;
   const double x36 = x31*x4;
   const double x37 = x23*x36;
   const double x38 = x24*x27;
   const double x39 = x23*x38;
   const double x40 = 2*x39;
   const double x41 = 1.0/x3;
   const double x42 = h*x31;
   const double x43 = -v*x34*x41*x42;
   const double x44 = fabs(v);
   const double x45 = u == x44;
   const double x46 = x4/h != 0;
   const double x47 = u == 0;
   const double x48 = v == 0;
   const double x49 = u == -x44;
   const double x50 = x45 && x46 || x46 && x49 || x45 && x46 && x47 || x45 && x46 && x48 || x45 && x46 && x49 || x46 && x47 && x48 || x46 && x47 && x49 || x46 && x48 && x49 || x45 && x46 && x47 && x48 || x45 && x46 && x47 && x49 || x45 && x46 && x48 && x49 || x46 && x47 && x48 && x49 || x45 && x46 && x47 && x48 && x49;
   const double x51 = x45 || x49 || x45 && x47 || x45 && x48 || x45 && x49 || x47 && x48 || x47 && x49 || x48 && x49 || x45 && x47 && x48 || x45 && x47 && x49 || x45 && x48 && x49 || x47 && x48 && x49 || x45 && x47 && x48 && x49;
   const double x52 = x34/pow(x3, 2);
   const double x53 = h*x52;
   const double x54 = pow(x3, -3.0/2.0);
   const double x55 = h*x9;
   const double x56 = x10*x31;
   const double x57 = x34*x54*(-x27*x55 - x56*x9);
   const double x58 = -x1*x53 - x1*x57;
   const double x59 = x46 || x46 && x47 || x46 && x48;
   const double x60 = -mu*(u*x37 + u*x40) + u*x26 + v*x35*(-x28*x30 + x33) + x1*x14 + ((x50) ? (
      x43
   )
   : ((x51) ? (
      0
   )
   : ((x59) ? (
      x43 + x58
   )
   : (
      x58
   ))));
   const double x61 = u*x41;
   const double x62 = x61*(1 - x27);
   const double x63 = v*x62;
   const double x64 = x41*(x27 - 1);
   const double x65 = x2*x64 + 1;
   const double x66 = u*x2;
   const double x67 = x2*x29;
   const double x68 = x0*x31;
   const double x69 = x29*x68;
   const double x70 = x34*x61;
   const double x71 = x42*x70;
   const double x72 = x53*x66 + x57*x66;
   const double x73 = -mu*(-v*x37 - v*x40) - u*x35*(-x28*x67 + x69) - v*x26 - x14*x66 + ((x50) ? (
      x71
   )
   : ((x51) ? (
      0
   )
   : ((x59) ? (
      x71 + x72
   )
   : (
      x72
   ))));
   const double x74 = (1.0/2.0)*x4;
   const double x75 = (1.0/2.0)*x38;
   const double x76 = x54*(x74 - x75);
   const double x77 = x0*x76;
   const double x78 = h*x34;
   const double x79 = x2*x76;
   const double x80 = pow(x27, 2);
   const double x81 = x4*x80;
   const double x82 = 6*x4;
   const double x83 = x15/(pow(u, 6)*x5 + pow(v, 6)*x5 + x0*x18*x82 + x16*x2*x82);
   const double x84 = x16*x83;
   const double x85 = x31*x83;
   const double x86 = x21*x83;
   const double x87 = -x36*x86 - x38*x86 - x81*x86;
   const double x88 = -mu*(-x36*x84 - x38*x84 + x6*x85 - x81*x84 + x87);
   const double x89 = x18*x83;
   const double x90 = -mu*(-x36*x89 - x38*x89 + x7*x85 - x81*x89 + x87);
   const double x91 = x88 + x90 + ((x46) ? (
      x13*(x74 + x75) + x77*x78 + x78*x79
   )
   : (
      x12
   ));
   const double x92 = x24*x9;
   const double x93 = v*x92;
   const double x94 = x60*x63 + x65*x73 + x91*x93;
   const double x95 = u*v;
   const double x96 = x25*x30;
   const double x97 = x25*x67;
   const double x98 = v*x70;
   const double x99 = u*pow(v, 3);
   const double x100 = x83*x99;
   const double x101 = pow(u, 3)*v;
   const double x102 = x101*x83;
   const double x103 = x81*x83;
   const double x104 = -mu*(-x100*x36 - x100*x38 - x101*x103 - x102*x36 - x102*x38 - x103*x99 + x5*x85*x95) + x98*(x29*x6 - x33*x4 - x38*x67 - x67*x81 - x96 + x97) + x98*(x29*x7 - x30*x38 - x30*x81 - x4*x69 + x96 - x97) + ((x46) ? (
      -x12*x76*x95
   )
   : (
      0
   ));
   const double x105 = h*x22;
   const double x106 = 4*x105*x21*x24;
   const double x107 = x105*x18;
   const double x108 = x12*x79;
   const double x109 = x21*x52*(h + x10*x80 - x25*x55 + x55*x75 + x56);
   const double x110 = x20*x23;
   const double x111 = x23*x5;
   const double x112 = x23*x4;
   const double x113 = x23*x81;
   const double x114 = mu*(x110*x27 - x110) - mu*(-x0*x113 + x0*x39 + x111*x80 - x111 - x112*x32 - x112*x68 - x113*x2 + x2*x39);
   const double x115 = x114 + x34*(x105*x17 + x106 + x107*x36 + x107*x38 + x107*x81) + x90 + ((x50) ? (
      x108
   )
   : ((x51) ? (
      0
   )
   : ((x59) ? (
      x108 + x109
   )
   : (
      x109
   ))));
   const double x116 = x104*x65 + x115*x63 + x60*x93;
   const double x117 = v*x73;
   const double x118 = x117*x92;
   const double x119 = x104*x63;
   const double x120 = x105*x16;
   const double x121 = x12*x77;
   const double x122 = x114 + x34*(x105*x19 + x106 + x120*x36 + x120*x38 + x120*x81) + x88 + ((x50) ? (
      x121
   )
   : ((x51) ? (
      0
   )
   : ((x59) ? (
      x109 + x121
   )
   : (
      x109
   ))));
   const double x123 = x118 + x119 + x122*x65;
   const double x124 = u*x92;
   const double x125 = x0*x64 + 1;
   const double x126 = x117*x62 - x124*x91 + x125*x60;
   const double x127 = x124*x60;
   const double x128 = x115*x125 + x119 - x127;
   const double x129 = x104*x125 + x122*x63 - x124*x73;
   const double x130 = -x118 + x127 + x27*x91;
   const double x131 = -x104*x93 + x115*x124 + x27*x60;
   const double x132 = x104*x124 - x122*x93 + x27*x73;

   inertia_h[0] = x116*x63 + x123*x65 + x93*x94;
   inertia_h[1] = x116*x125 + x123*x63 - x124*x94;
   inertia_h[2] = x116*x124 - x123*x93 + x27*x94;
   inertia_h[3] = x126*x93 + x128*x63 + x129*x65;
   inertia_h[4] = -x124*x126 + x125*x128 + x129*x63;
   inertia_h[5] = x124*x128 + x126*x27 - x129*x93;
   inertia_h[6] = x130*x93 + x131*x63 + x132*x65;
   inertia_h[7] = -x124*x130 + x125*x131 + x132*x63;
   inertia_h[8] = x124*x131 + x130*x27 - x132*x93;

}

MatrixXd Grub_rne::J_dot_udot(double u, double v, double h) {
    // updated:02/26/2021
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = cos(x3);
    const double x5 = x4 - 1;
    const double x6 = pow(x5, 2);
    const double x7 = pow(x2, -2);
    const double x8 = h*x7;
    const double x9 = v*x8;
    const double x10 = pow(x2, -3);
    const double x11 = h*x10;
    const double x12 = 4*x11;
    const double x13 = x12*x6;
    const double x14 = sin(x3);
    const double x15 = x0*x14;
    const double x16 = pow(x2, -5.0/2.0);
    const double x17 = h*x5;
    const double x18 = x16*x17;
    const double x19 = v*x18;
    const double x20 = x15*x19;
    const double x21 = -x14 + x3;
    const double x22 = x14*x21;
    const double x23 = x22*x9;
    const double x24 = h*x4;
    const double x25 = x0*x21;
    const double x26 = v*x16;
    const double x27 = x24*x25*x26;
    const double x28 = x15*x21;
    const double x29 = v*x12*x28;
    const double x30 = 1.0/x3;
    const double x31 = u*x30;
    const double x32 = x31*x4;
    const double x33 = x31 - x32;
    const double x34 = u*x14;
    const double x35 = x33*x34*x9;
    const double x36 = pow(x2, -3.0/2.0);
    const double x37 = h*x36;
    const double x38 = x34*x37;
    const double x39 = 2*u;
    const double x40 = x39*x8;
    const double x41 = x40*x5;
    const double x42 = pow(u, 3);
    const double x43 = h*x16;
    const double x44 = 1 - x4;
    const double x45 = x42*x44;
    const double x46 = 4*x10*x17;
    const double x47 = x1*x3 + x15;
    const double x48 = x16*x47;
    const double x49 = u*x24;
    const double x50 = x34*x47;
    const double x51 = x30*x42;
    const double x52 = x1*x31 + x14*x39 + x4*x51;
    const double x53 = x14*x8;
    const double x54 = x1*x14;
    const double x55 = x0*x3 + x54;
    const double x56 = x11*x55;
    const double x57 = x16*x55;
    const double x58 = x1*x32 + x3*x39 + x51;
    const double x59 = v*x43;
    const double x60 = x0*x44;
    const double x61 = x21*x36;
    const double x62 = 3*x25;
    const double x63 = x33*x37;
    const double x64 = x1*x18;
    const double x65 = 5*x17/pow(x2, 7.0/2.0);
    const double x66 = u*x33;
    const double x67 = 3*u;
    const double x68 = v*x11;
    const double x69 = v*x65;
    const double x70 = v*x36;
    const double x71 = v*x61 - x26*x62 + x66*x70;

    double *J_dot_udot = new double[12];

    J_dot_udot[0] = v*x0*x13 + 2*x20 + x23 + x27 - x29 + x35 - x6*x9;
    J_dot_udot[1] = -x12*x50 - x14*x18*x42 + x14*x43*x45 - x38 - x41*x44 - x41 + x45*x46 + x48*x49 + x52*x53;
    J_dot_udot[2] = u*x1*x13 + x18*x39*x54 + 4*x34*x56 + x38 - x40*x44 - x49*x57 - x53*x58;
    J_dot_udot[3] = v*x46*x60 + x15*x44*x59 - x20 - x23 - x27 + x29 - x35 - x44*x5*x9;
    J_dot_udot[4] = h*x61 + u*x18*x58 + u*x63 - x0*x55*x65 - x1*x11*x28 - x1*x25*x65 - x15*x56 + x17*x57 + x21*x64 - x43*x62 + x64*x66;
    J_dot_udot[5] = -u*x47*x69 + v*x63 + x0*x19*x33 + x19*x21*x39 + x19*x52 - x21*x42*x69 - x21*x59*x67 - x22*x42*x68 - x50*x68;
    J_dot_udot[6] = x36*x58 - x57*x67;
    J_dot_udot[7] = x71;
    J_dot_udot[8] = x71;
    J_dot_udot[9] = x36*x52 - x48*x67;
    J_dot_udot[10] = -v*x39*x5*x7 - x34*x70;
    J_dot_udot[11] = x15*x36 - 2*x60*x7 + x44/x2;

    MatrixXd J_dot(6,2);
    J_dot << J_dot_udot[0], J_dot_udot[1], J_dot_udot[2], J_dot_udot[3],
                J_dot_udot[4], J_dot_udot[5], J_dot_udot[6], J_dot_udot[7],
                J_dot_udot[8], J_dot_udot[9], J_dot_udot[10], J_dot_udot[11];

    delete J_dot_udot;

    return J_dot;

}

MatrixXd Grub_rne::J_dot_vdot(double u, double v, double h) {
    // Updated: 02/26/2021
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = cos(x3);
    const double x5 = x4 - 1;
    const double x6 = pow(x5, 2);
    const double x7 = pow(x2, -2);
    const double x8 = h*x7;
    const double x9 = u*x8;
    const double x10 = pow(x2, -3);
    const double x11 = h*x10;
    const double x12 = 4*x11;
    const double x13 = x12*x6;
    const double x14 = sin(x3);
    const double x15 = x1*x14;
    const double x16 = pow(x2, -5.0/2.0);
    const double x17 = h*x5;
    const double x18 = x16*x17;
    const double x19 = u*x18;
    const double x20 = x15*x19;
    const double x21 = -x14 + x3;
    const double x22 = x14*x21;
    const double x23 = x22*x9;
    const double x24 = h*x4;
    const double x25 = x1*x21;
    const double x26 = u*x16;
    const double x27 = x24*x25*x26;
    const double x28 = u*x21;
    const double x29 = x12*x15*x28;
    const double x30 = 1.0/x3;
    const double x31 = v*x30;
    const double x32 = x31*x4;
    const double x33 = x31 - x32;
    const double x34 = v*x14;
    const double x35 = x33*x34*x9;
    const double x36 = pow(x2, -3.0/2.0);
    const double x37 = h*x36;
    const double x38 = x34*x37;
    const double x39 = 2*v;
    const double x40 = x39*x8;
    const double x41 = x0*x14;
    const double x42 = v*x18;
    const double x43 = h*x16;
    const double x44 = 1 - x4;
    const double x45 = v*x44;
    const double x46 = 4*x10;
    const double x47 = x0*x17;
    const double x48 = x1*x3 + x41;
    const double x49 = x16*x48;
    const double x50 = v*x24;
    const double x51 = x12*x34;
    const double x52 = pow(v, 3);
    const double x53 = x30*x52;
    const double x54 = x0*x32 + x3*x39 + x53;
    const double x55 = x14*x8;
    const double x56 = x0*x3 + x15;
    const double x57 = x16*x56;
    const double x58 = x0*x31 + x14*x39 + x4*x53;
    const double x59 = u*x44;
    const double x60 = x1*x17;
    const double x61 = 3*v;
    const double x62 = x33*x37;
    const double x63 = u*x11;
    const double x64 = 5/pow(x2, 7.0/2.0);
    const double x65 = x17*x64;
    const double x66 = x21*x36;
    const double x67 = 3*x25;
    const double x68 = x0*x18;
    const double x69 = v*x33;
    const double x70 = u*x36;
    const double x71 = u*x66 - x26*x67 + x69*x70;

    double *J_dot_vdot = new double[12];

    J_dot_vdot[0] = u*x1*x13 + 2*x20 + x23 + x27 - x29 + x35 - x6*x9;
    J_dot_vdot[1] = -x38 - x40*x5 - x41*x42 + x41*x43*x45 + x45*x46*x47 - x48*x51 + x49*x50 + x54*x55;
    J_dot_vdot[2] = x13*x52 + 2*x14*x18*x52 + x38 - x40*x44 - x40*x6 - x50*x57 + x51*x56 - x55*x58;
    J_dot_vdot[3] = x15*x43*x59 - x20 - x23 - x27 + x29 - x35 - x44*x5*x9 + x46*x59*x60;
    J_dot_vdot[4] = -u*v*x56*x65 + u*x62 + x1*x19*x33 + x19*x21*x39 + x19*x58 - x22*x52*x63 - x28*x43*x61 - x28*x52*x65 - x34*x56*x63;
    J_dot_vdot[5] = h*x66 + v*x62 - x11*x15*x48 - x11*x25*x41 + x17*x49 + x21*x68 - x25*x47*x64 + x42*x54 - x43*x67 - x48*x60*x64 + x68*x69;
    J_dot_vdot[6] = x36*x58 - x57*x61;
    J_dot_vdot[7] = x71;
    J_dot_vdot[8] = x71;
    J_dot_vdot[9] = x36*x54 - x49*x61;
    J_dot_vdot[10] = -2*x1*x5*x7 - x15*x36 + x5/x2;
    J_dot_vdot[11] = x34*x70 - x39*x59*x7;

    MatrixXd J_dot(6,2);
    J_dot << J_dot_vdot[0], J_dot_vdot[1], J_dot_vdot[2], J_dot_vdot[3],
                J_dot_vdot[4], J_dot_vdot[5], J_dot_vdot[6], J_dot_vdot[7],
                J_dot_vdot[8], J_dot_vdot[9], J_dot_vdot[10], J_dot_vdot[11];
    
    delete J_dot_vdot;

    return J_dot;
}

MatrixXd Grub_rne::J_end(double u, double v, double h) {
    // Updated: 02/26/2021
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = cos(x3);
    const double x5 = x4 - 1;
    const double x6 = h/pow(x2, 2);
    const double x7 = pow(x5, 2)*x6;
    const double x8 = u*v;
    const double x9 = sin(x3);
    const double x10 = x3 - x9;
    const double x11 = x6*x9;
    const double x12 = x10*x11*x8;
    const double x13 = 1.0/x2;
    const double x14 = h*x13;
    const double x15 = 1 - x4;
    const double x16 = x0*x9 + x1*x3;
    const double x17 = x0*x3 + x1*x9;
    const double x18 = v*x5;
    const double x19 = u*x15;
    const double x20 = h*u;
    const double x21 = pow(x2, -3.0/2.0);
    const double x22 = x10*x21;
    const double x23 = pow(x2, -5.0/2.0);
    const double x24 = x20*x23*x5;
    const double x25 = v*x22;
    const double x26 = h*x18*x23;
    const double x27 = u*x25;

    double *J = new double[12];

    J[0] = x12 - x7*x8;
    J[1] = -x0*x15*x5*x6 + x11*x16 + x14*x5;
    J[2] = -x1*x7 - x11*x17 + x14*x15;
    J[3] = -x12 - x18*x19*x6;
    J[4] = x1*x10*x24 + x17*x24 + x20*x22;
    J[5] = h*x25 + x0*x10*x26 + x16*x26;
    J[6] = x17*x21;
    J[7] = x27;
    J[8] = x27;
    J[9] = x16*x21;
    J[10] = x13*x18;
    J[11] = x13*x19;

    MatrixXd J_rel(6,2);

    J_rel << J[0], J[1], 
                J[2], J[3],
                J[4], J[5], 
                J[6], J[7],
                J[8], J[9], 
                J[10], J[11];

    delete J;

    return J_rel;
}

MatrixXd Grub_rne::J_dot_udot_com(double u, double v, double h) {
    // updated:03/03/2021
    double l;
    l = _l;
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = l*x3;
    const double x5 = cos(x4);
    const double x6 = x5 - 1;
    const double x7 = sin(x3);
    const double x8 = x3 - x7;
    const double x9 = h*x8;
    const double x10 = pow(x2, -5.0/2.0);
    const double x11 = v*x10;
    const double x12 = x11*x9;
    const double x13 = x0*x7;
    const double x14 = sin(x4);
    const double x15 = -x14 + x4;
    const double x16 = pow(x2, -3);
    const double x17 = x15*x16;
    const double x18 = v*x17;
    const double x19 = h*x13*x18;
    const double x20 = l*x0;
    const double x21 = x16*x9;
    const double x22 = x14*x21;
    const double x23 = v*x20*x22;
    const double x24 = cos(x3);
    const double x25 = h*(1 - x24);
    const double x26 = x11*x15;
    const double x27 = x25*x26;
    const double x28 = v*x6;
    const double x29 = 5/pow(x2, 7.0/2.0);
    const double x30 = x29*x9;
    const double x31 = v*x0;
    const double x32 = x15*x25*x29*x31;
    const double x33 = 1.0/x3;
    const double x34 = u*x33;
    const double x35 = -x24*x34 + x34;
    const double x36 = h*x35;
    const double x37 = u*x36;
    const double x38 = x11*x37;
    const double x39 = l*x5;
    const double x40 = l*x34 - x34*x39;
    const double x41 = u*x40;
    const double x42 = x11*x25*x41;
    const double x43 = l*x1;
    const double x44 = x0*x24;
    const double x45 = x1*x24;
    const double x46 = x20*x24;
    const double x47 = x24*x43;
    const double x48 = x0*x5;
    const double x49 = x1*x5;
    const double x50 = 2*x3;
    const double x51 = x50*x7;
    const double x52 = x14*x3;
    const double x53 = l - 1;
    const double x54 = x3*x53;
    const double x55 = sin(x54);
    const double x56 = x3*x55;
    const double x57 = -x0 - x1 - x20 - x43 - x44 - x45 + x46 + x47 - x48 - x49 + x51 + x52 - x56;
    const double x58 = pow(v, 4);
    const double x59 = 3*x0;
    const double x60 = pow(u, 4);
    const double x61 = 3*x1;
    const double x62 = pow(u, 6) + pow(v, 6) + x58*x59 + x60*x61;
    const double x63 = h/x62;
    const double x64 = v*x63;
    const double x65 = pow(u, 5);
    const double x66 = pow(u, 3);
    const double x67 = h*(-6*u*x58 - 12*x1*x66 - 6*x65)/pow(x62, 2);
    const double x68 = u*v;
    const double x69 = x67*x68;
    const double x70 = 2*u;
    const double x71 = l*x70;
    const double x72 = x24*x71;
    const double x73 = x5*x70;
    const double x74 = u*x39;
    const double x75 = x53*cos(x54);
    const double x76 = u*x75;
    const double x77 = x33*x7;
    const double x78 = x70*x77;
    const double x79 = x14*x34;
    const double x80 = x33*x66;
    const double x81 = x7*x80;
    const double x82 = x1*x7;
    const double x83 = x34*x82;
    const double x84 = l*x81;
    const double x85 = x34*x55;
    const double x86 = x14*x80;
    const double x87 = l*x86;
    const double x88 = x43*x7;
    const double x89 = x34*x88;
    const double x90 = x43*x79;
    const double x91 = u*x64;
    const double x92 = l*x14;
    const double x93 = 1 - x5;
    const double x94 = x9*x93;
    const double x95 = x0*x93;
    const double x96 = x10*x36;
    const double x97 = pow(x2, -3.0/2.0);
    const double x98 = u*x97;
    const double x99 = l*x58;
    const double x100 = x1*x20;
    const double x101 = x0*x52;
    const double x102 = x100 + x101 + x99;
    const double x103 = x0*x1;
    const double x104 = 2*x103 + x58 + x60;
    const double x105 = 1.0/x104;
    const double x106 = h*x105;
    const double x107 = x102*x106;
    const double x108 = pow(x2, -2);
    const double x109 = x108*x70;
    const double x110 = x109*x25;
    const double x111 = 4*x66;
    const double x112 = 4*u;
    const double x113 = (-x1*x112 - x111)/pow(x104, 2);
    const double x114 = x102*x113;
    const double x115 = 1.0/x2;
    const double x116 = x115*x25;
    const double x117 = x43*x70;
    const double x118 = x39*x66;
    const double x119 = x52*x70;
    const double x120 = x105*(x117 + x118 + x119 + x86);
    const double x121 = x1*x44;
    const double x122 = x20*x45;
    const double x123 = x1*x48;
    const double x124 = x1*x70;
    const double x125 = x47*x70;
    const double x126 = x49*x70;
    const double x127 = x34*x7;
    const double x128 = x33*x70;
    const double x129 = x80*x82;
    const double x130 = x43*x81;
    const double x131 = x43*x86;
    const double x132 = u*x43;
    const double x133 = x1*x6;
    const double x134 = x14*x97;
    const double x135 = x1*x134 + x115*x20;
    const double x136 = x135*x97;
    const double x137 = u*x136;
    const double x138 = 2*x66;
    const double x139 = x132*x5;
    const double x140 = -l*x108*x138 - u*x10*x14*x61 + x108*x139 + x115*x71;
    const double x141 = l*x60;
    const double x142 = l*x111;
    const double x143 = x65*x77;
    const double x144 = x2 + x20 + x43 + x44 + x45 - x46 - x47 + x48 + x49 - x51 - x52 + x56;
    const double x145 = x1*x9;
    const double x146 = 6*x15*x9/pow(x2, 4);
    const double x147 = x0*x4 + x1*x4 + x13 - x20*x7 + x24*x50 - x50 + x82 - x88;
    const double x148 = x106*x147;
    const double x149 = h*x113*x147;
    const double x150 = l*x80;
    const double x151 = x106*(x128*x24 - x128 - x150*x24 + x150 + x24*x80 + x34*x43 + x34*x45 - x34*x47 + x4*x70 - x7*x71);
    const double x152 = v*x97;
    const double x153 = x152*x9;
    const double x154 = v*x33;
    const double x155 = x15*x152 + x152*x41 - x26*x59;

    double *J_dot_udot_com = new double[12];

    J_dot_udot_com[0] = -x0*x28*x30 + x12*x6 + x19 - x23 + x27 - x32 + x38*x6 + x42 + x57*x64 + x57*x69 + x91*(-x70 - x71 + x72 - x73 + x74 - x76 + x78 + x79 + x81 + x83 - x84 - x85 + x87 - x89 + x90);
    J_dot_udot_com[1] = x10*x70*x94 - x102*x105*x110 + x107*x7*x98 + x114*x116 + x116*x120 + x21*x66*x92 - x29*x66*x94 + x63*(x111*x5 - x117 - x118 - x119 - x124 + x125 + x126 + x127*x58 - x127*x99 + x128*x82 + x129 - x130 - x131 - x33*x65*x92 + x55*x80 + x56*x70 + x66*x75 - x86) + x67*(x0*x56 - x100 - x101 - x103 - x121 + x122 + x123 - x24*x58 + x24*x99 + x5*x60 + x50*x82 - x58 - x99) + x95*x96;
    J_dot_udot_com[2] = -h*x137*x7 - u*x133*x30 + x110*x135 - x116*x140 - x132*x22 + x133*x96 + x63*(l*x143 - x1*x76 + x1*x79 - x1*x85 + x111 - x112*x3*x7 + x117 + x124 - x125 - x126 - x129 + x130 + x131 + x138*x24 + x139 - x142*x24 + x142 - x143 + x45*x70 + x79*x99 - 2*x81) + x67*(x1*x52 - x1*x56 + x100 + x103 + x121 - x122 - x123 - x13*x50 - x141*x24 + x141 + x24*x60 - x5*x58 + x60);
    J_dot_udot_com[3] = -v*x30*x95 + x12*x93 + x144*x64 + x144*x69 - x19 + x23 - x27 + x32 + x38*x93 - x42 + x91*(x70 + x71 - x72 + x73 - x74 + x76 - x78 - x79 - x81 - x83 + x84 + x85 - x87 + x89 - x90);
    J_dot_udot_com[4] = -x0*x148*x97 - x1*x17*x37 + x10*x135*x59*x9 + x103*x146 - x136*x9 - x137*x36 - x140*x9*x98 - x145*x16*x41 - x145*x17 + x148*x33 + x149*x34 + x151*x34;
    J_dot_udot_com[5] = 3*u*x107*x11*x8 - u*x148*x152 + v*x146*x66 - x0*x18*x36 - x107*x152*x35 - x114*x153 - x120*x153 + x149*x154 + x151*x154 - x18*x70*x9 - x21*x31*x40;
    J_dot_udot_com[6] = x140;
    J_dot_udot_com[7] = x155;
    J_dot_udot_com[8] = x155;
    J_dot_udot_com[9] = x114 + x120;
    J_dot_udot_com[10] = -l*x134*x68 - x109*x28;
    J_dot_udot_com[11] = -2*x108*x95 + x115*x93 + x134*x20;

    MatrixXd J_dot(6,2);
    J_dot << J_dot_udot_com[0], J_dot_udot_com[1], 
                J_dot_udot_com[2], J_dot_udot_com[3],
                J_dot_udot_com[4], J_dot_udot_com[5],
                J_dot_udot_com[6], J_dot_udot_com[7],
                J_dot_udot_com[8], J_dot_udot_com[9],
                J_dot_udot_com[10], J_dot_udot_com[11];

    delete J_dot_udot_com;

    return J_dot;
}

MatrixXd Grub_rne::J_dot_vdot_com(double u, double v, double h) {
    // Updated: 03/03/2021
    double l;
    l = _l;
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = l*x3;
    const double x5 = cos(x4);
    const double x6 = x5 - 1;
    const double x7 = pow(x2, -5.0/2.0);
    const double x8 = u*x7;
    const double x9 = sin(x3);
    const double x10 = x3 - x9;
    const double x11 = h*x10;
    const double x12 = x11*x8;
    const double x13 = x1*x9;
    const double x14 = sin(x4);
    const double x15 = -x14 + x4;
    const double x16 = pow(x2, -3);
    const double x17 = x15*x16;
    const double x18 = h*u;
    const double x19 = x13*x17*x18;
    const double x20 = l*x1;
    const double x21 = x11*x16;
    const double x22 = x14*x21;
    const double x23 = u*x20*x22;
    const double x24 = cos(x3);
    const double x25 = h*(1 - x24);
    const double x26 = x15*x8;
    const double x27 = x25*x26;
    const double x28 = x1*x6;
    const double x29 = 5/pow(x2, 7.0/2.0);
    const double x30 = x11*x29;
    const double x31 = u*x1;
    const double x32 = x15*x25*x29*x31;
    const double x33 = 1.0/x3;
    const double x34 = v*x33;
    const double x35 = -x24*x34 + x34;
    const double x36 = h*x35;
    const double x37 = v*x36;
    const double x38 = x37*x8;
    const double x39 = l*x5;
    const double x40 = l*x34 - x34*x39;
    const double x41 = v*x40;
    const double x42 = x25*x41*x8;
    const double x43 = l*x0;
    const double x44 = x0*x24;
    const double x45 = x1*x24;
    const double x46 = x24*x43;
    const double x47 = x20*x24;
    const double x48 = x0*x5;
    const double x49 = x1*x5;
    const double x50 = 2*x3;
    const double x51 = x50*x9;
    const double x52 = x14*x3;
    const double x53 = l - 1;
    const double x54 = x3*x53;
    const double x55 = sin(x54);
    const double x56 = x3*x55;
    const double x57 = -x0 - x1 - x20 - x43 - x44 - x45 + x46 + x47 - x48 - x49 + x51 + x52 - x56;
    const double x58 = pow(v, 4);
    const double x59 = pow(u, 4);
    const double x60 = 3*x1;
    const double x61 = pow(u, 6) + pow(v, 6) + 3*x0*x58 + x59*x60;
    const double x62 = h/x61;
    const double x63 = u*x62;
    const double x64 = pow(v, 5);
    const double x65 = pow(v, 3);
    const double x66 = h*(-6*v*x59 - 12*x0*x65 - 6*x64)/pow(x61, 2);
    const double x67 = u*v;
    const double x68 = x66*x67;
    const double x69 = 2*v;
    const double x70 = l*x69;
    const double x71 = x24*x70;
    const double x72 = x5*x69;
    const double x73 = v*x39;
    const double x74 = x53*cos(x54);
    const double x75 = v*x74;
    const double x76 = x33*x9;
    const double x77 = x69*x76;
    const double x78 = x14*x34;
    const double x79 = x65*x76;
    const double x80 = x0*x9;
    const double x81 = x34*x80;
    const double x82 = l*x79;
    const double x83 = x34*x55;
    const double x84 = x33*x65;
    const double x85 = x14*x84;
    const double x86 = l*x85;
    const double x87 = x43*x9;
    const double x88 = x34*x87;
    const double x89 = x43*x78;
    const double x90 = v*x63;
    const double x91 = v*x43;
    const double x92 = 1 - x5;
    const double x93 = x0*x92;
    const double x94 = x36*x7;
    const double x95 = pow(x2, -3.0/2.0);
    const double x96 = v*x95;
    const double x97 = x9*x96;
    const double x98 = l*x58;
    const double x99 = x1*x43;
    const double x100 = x0*x52;
    const double x101 = x100 + x98 + x99;
    const double x102 = x0*x1;
    const double x103 = 2*x102 + x58 + x59;
    const double x104 = 1.0/x103;
    const double x105 = h*x104;
    const double x106 = x101*x105;
    const double x107 = pow(x2, -2);
    const double x108 = x107*x69;
    const double x109 = x108*x25;
    const double x110 = 4*x65;
    const double x111 = -x110;
    const double x112 = 4*v;
    const double x113 = (-x0*x112 + x111)/pow(x103, 2);
    const double x114 = x101*x113;
    const double x115 = 1.0/x2;
    const double x116 = x115*x25;
    const double x117 = l*x110;
    const double x118 = x43*x69;
    const double x119 = x5*x91;
    const double x120 = x0*x78;
    const double x121 = x104*(x117 + x118 + x119 + x120);
    const double x122 = x1*x44;
    const double x123 = x43*x45;
    const double x124 = x1*x48;
    const double x125 = x0*x69;
    const double x126 = x24*x65;
    const double x127 = x46*x69;
    const double x128 = x48*x69;
    const double x129 = x64*x76;
    const double x130 = x80*x84;
    const double x131 = x43*x79;
    const double x132 = l*x59;
    const double x133 = x43*x85;
    const double x134 = x11*x6;
    const double x135 = l*x14;
    const double x136 = x14*x95;
    const double x137 = x1*x136 + x115*x43;
    const double x138 = 3*x7;
    const double x139 = x39*x65;
    const double x140 = -x107*x118 + x107*x139 + x136*x69 - x138*x14*x65;
    const double x141 = x33*x69;
    const double x142 = x34*x9;
    const double x143 = u*x92;
    const double x144 = x2 + x20 + x43 + x44 + x45 - x46 - x47 + x48 + x49 - x51 - x52 + x56;
    const double x145 = u*x11;
    const double x146 = 6*x15/pow(x2, 4);
    const double x147 = u*x95;
    const double x148 = x0*x4 + x1*x4 + x13 - x20*x9 + x24*x50 - x50 + x80 - x87;
    const double x149 = x113*x148;
    const double x150 = x126*x33;
    const double x151 = -l*x150 + l*x84 + x141*x24 - x141 + x150 + x34*x43 + x34*x44 - x34*x46 + x4*x69 - x70*x9;
    const double x152 = x105*x33;
    const double x153 = x0*x11;
    const double x154 = x105*x95;
    const double x155 = x10*x101;
    const double x156 = x11*x96;
    const double x157 = x147*x15 + x147*x41 - x26*x60;

    double *J_dot_vdot_com = new double[12];

    J_dot_vdot_com[0] = -u*x28*x30 + x12*x6 + x19 - x23 + x27 - x32 + x38*x6 + x42 + x57*x63 + x57*x68 + x90*(-x69 - x70 + x71 - x72 + x73 - x75 + x77 + x78 + x79 + x81 - x82 - x83 + x86 - x88 + x89);
    J_dot_vdot_com[1] = -v*x30*x93 - x101*x104*x109 + x106*x97 + x114*x116 + x116*x121 + x22*x91 + x62*(-l*x129 + x0*x75 + x0*x83 + x111 + x112*x3*x9 + x117*x24 - x117 - x118 - x119 - x120 - x125 - 2*x126 + x127 + x128 + x129 + x130 - x131 - x132*x78 - x133 - x44*x69 + 2*x79) + x66*(x0*x56 - x100 - x102 - x122 + x123 + x124 + x13*x50 - x24*x58 + x24*x98 + x5*x59 - x58 - x98 - x99) + x93*x94;
    J_dot_vdot_com[2] = -h*x137*x97 + x109*x137 - x116*x140 - x134*x29*x65 + x134*x69*x7 - x135*x21*x65 + x28*x94 + x62*(-x110*x5 + x118 + x125 - x127 - x128 - x130 + x131 + x132*x142 + x133 + x135*x33*x64 + x139 - x141*x80 - x142*x59 + x52*x69 - x55*x84 - x56*x69 - x65*x74 + x85) + x66*(x1*x52 - x1*x56 + x102 + x122 - x123 - x124 - x132*x24 + x132 + x24*x59 - x5*x58 - x50*x80 + x59 + x99);
    J_dot_vdot_com[3] = -x1*x143*x30 + x12*x92 + x144*x63 + x144*x68 - x19 + x23 - x27 + x32 + x38*x92 - x42 + x90*(x69 + x70 - x71 + x72 - x73 + x75 - x77 - x78 - x79 - x81 + x82 + x83 - x86 + x88 - x89);
    J_dot_vdot_com[4] = u*x151*x152 - v*x105*x147*x148 + x11*x137*x138*x67 - x11*x140*x147 - x137*x147*x36 + x145*x146*x65 - x145*x17*x69 + x149*x18*x33 - x17*x31*x36 - x21*x31*x40;
    J_dot_vdot_com[5] = h*x149*x34 - x0*x17*x37 - x1*x148*x154 + x102*x11*x146 + x105*x151*x34 + x105*x155*x60*x7 - x106*x35*x96 - x114*x156 - x121*x156 + x148*x152 - x153*x16*x41 - x153*x17 - x154*x155;
    J_dot_vdot_com[6] = x140;
    J_dot_vdot_com[7] = x157;
    J_dot_vdot_com[8] = x157;
    J_dot_vdot_com[9] = x114 + x121;
    J_dot_vdot_com[10] = -2*x107*x28 + x115*x6 - x136*x20;
    J_dot_vdot_com[11] = l*x136*x67 - x108*x143;

    MatrixXd J_dot(6,2);
    J_dot << J_dot_vdot_com[0], J_dot_vdot_com[1], 
                J_dot_vdot_com[2], J_dot_vdot_com[3],
                J_dot_vdot_com[4], J_dot_vdot_com[5],
                J_dot_vdot_com[6], J_dot_vdot_com[7],
                J_dot_vdot_com[8], J_dot_vdot_com[9],
                J_dot_vdot_com[10], J_dot_vdot_com[11];

    delete J_dot_vdot_com;

    return J_dot;
}

MatrixXd Grub_rne::J_com(double u, double v, double h) {
    // Updated: 02/26/2021
    double l;
    l = _l;
    const double x0 = pow(u, 2);
    const double x1 = pow(v, 2);
    const double x2 = x0 + x1;
    const double x3 = sqrt(x2);
    const double x4 = l*x3;
    const double x5 = cos(x4);
    const double x6 = x5 - 1;
    const double x7 = v*x6;
    const double x8 = sin(x3);
    const double x9 = x3 - x8;
    const double x10 = h/pow(x2, 5.0/2.0);
    const double x11 = x10*x9;
    const double x12 = cos(x3);
    const double x13 = 1 - x12;
    const double x14 = sin(x4);
    const double x15 = -x14 + x4;
    const double x16 = u*v;
    const double x17 = x15*x16;
    const double x18 = x10*x13*x17;
    const double x19 = l*x0;
    const double x20 = l*x1;
    const double x21 = x0*x12;
    const double x22 = x1*x12;
    const double x23 = x12*x19;
    const double x24 = x12*x20;
    const double x25 = x0*x5;
    const double x26 = x1*x5;
    const double x27 = 2*x3;
    const double x28 = x27*x8;
    const double x29 = x14*x3;
    const double x30 = x3*sin(x3*(l - 1));
    const double x31 = pow(v, 4);
    const double x32 = pow(u, 4);
    const double x33 = h/(pow(u, 6) + pow(v, 6) + 3*x0*x31 + 3*x1*x32);
    const double x34 = x16*x33;
    const double x35 = 1 - x5;
    const double x36 = 1.0/x2;
    const double x37 = x13*x36;
    const double x38 = x0*x1;
    const double x39 = 1.0/(x31 + x32 + 2*x38);
    const double x40 = l*x31;
    const double x41 = x1*x19;
    const double x42 = x0*x29;
    const double x43 = x39*(x40 + x41 + x42);
    const double x44 = h*x43;
    const double x45 = x1*x21;
    const double x46 = x19*x22;
    const double x47 = x1*x25;
    const double x48 = x1*x8;
    const double x49 = pow(x2, -3.0/2.0);
    const double x50 = x1*x14*x49 + x19*x36;
    const double x51 = h*x50;
    const double x52 = l*x32;
    const double x53 = x0*x8;
    const double x54 = u*x35;
    const double x55 = h*u;
    const double x56 = x15*x9/pow(x2, 3);
    const double x57 = x49*x9;
    const double x58 = x39*(x0*x4 + x1*x4 + x12*x27 - x19*x8 - x20*x8 - x27 + x48 + x53)/x3;
    const double x59 = h*v;
    const double x60 = x17*x49;

    double *J_com = new double[12];

    J_com[0] = u*x11*x7 + x18 + x34*(-x0 - x1 - x19 - x20 - x21 - x22 + x23 + x24 - x25 - x26 + x28 + x29 - x30);
    J_com[1] = x0*x11*x35 + x33*(x0*x30 - x12*x31 + x12*x40 + x27*x48 - x31 + x32*x5 - x38 - x40 - x41 - x42 - x45 + x46 + x47) + x37*x44;
    J_com[2] = x1*x11*x6 + x33*(x1*x29 - x1*x30 + x12*x32 - x12*x52 - x27*x53 - x31*x5 + x32 + x38 + x41 + x45 - x46 - x47 + x52) - x37*x51;
    J_com[3] = v*x11*x54 - x18 + x34*(x19 + x2 + x20 + x21 + x22 - x23 - x24 + x25 + x26 - x28 - x29 + x30);
    J_com[4] = -u*x51*x57 - x1*x55*x56 + x55*x58;
    J_com[5] = -v*x44*x57 - x0*x56*x59 + x58*x59;
    J_com[6] = x50;
    J_com[7] = x60;
    J_com[8] = x60;
    J_com[9] = x43;
    J_com[10] = x36*x7;
    J_com[11] = x36*x54;

    MatrixXd J_rel(6,2);

    J_rel << J_com[0], J_com[1], 
                J_com[2], J_com[3],
                J_com[4], J_com[5], 
                J_com[6], J_com[7],
                J_com[8], J_com[9], 
                J_com[10], J_com[11];

    delete J_com;

    return J_rel;
}

void Grub_rne::get_e_rel(int link_num, VectorXd uvdot, VectorXd uvddot)
{
    VectorXd qdot(2), xi_dot(6), qddot(2), xi_ddot(6);
    double u, v, h, udot,vdot;
    int i;
    
    // test if it is a flexible joint
    if (_flexi_link(link_num)==1)
    {
        int joint_num=0;
        for (i=0; i<link_num; i++)
        {
            if (_flexi_link(i)==1)
            {
                joint_num++;
            }
        }

        u = _uv(joint_num*2);
        v = _uv(joint_num*2+1);
        h = _h(link_num);
        udot = uvdot(joint_num*2);
        vdot = uvdot(joint_num*2+1);
        qdot << udot,vdot;
        qddot << uvddot(joint_num*2),uvddot(joint_num*2+1);
        // std::cout << "qddot: " << qddot << std::endl; 
        // find omega from jacobian
        MatrixXd J = J_end(u,v,h);
        xi_dot = J*qdot;
        _v_e_rel = xi_dot.block<3,1>(0,0);
        _omega_e_rel = xi_dot.block<3,1>(3,0);
        // std::cout << "v_e_rel: " << _v_e_rel << std::endl; 
        // std::cout << "omega_rel: " << _omega_rel << std::endl; 
        // find omega from jacobian
        MatrixXd Jdot_u = J_dot_udot(u,v,h);
        MatrixXd Jdot_v = J_dot_vdot(u,v,h);
        MatrixXd Jdot = Jdot_u*udot+Jdot_v*vdot;
        // std::cout << "Jdot: " << Jdot << std::endl;  
        xi_ddot = J*qddot+Jdot*qdot;
        _a_e_rel = xi_ddot.block<3,1>(0,0);
        // std::cout << "a_e_rel: " << _a_e_rel << std::endl;  
        _alpha_e_rel = xi_ddot.block<3,1>(3,0);
        // std::cout << "alpha_rel: " << _alpha_rel << std::endl; 
        // omega_rel, alpha_rel, a_c_rel, v_c_rel, a_e_rel, v_e_rel
    }
    else
    {
        _v_e_rel.setZero();
        _omega_e_rel.setZero();
        _a_e_rel.setZero();
        _alpha_e_rel.setZero();
    }
}

void Grub_rne::get_c_rel(int link_num, VectorXd uvdot, VectorXd uvddot)
{
    VectorXd qdot(2), xi_dot(6), qddot(2), xi_ddot(6);
    double u, v, h, udot,vdot;
    int i;
    
    // test if it is a flexible joint
    if (_flexi_link(link_num)==1)
    {
        int joint_num=0;
        for (i=0; i<link_num; i++)
        {
            if (_flexi_link(i)==1)
            {
                joint_num++;
            }
        }

        u = _uv(joint_num*2);
        v = _uv(joint_num*2+1);
        h = _h(link_num);
        udot = uvdot(joint_num*2);
        vdot = uvdot(joint_num*2+1);
        qdot << udot,vdot;
        qddot << uvddot(joint_num*2),uvddot(joint_num*2+1);
        // find omega from jacobian
        MatrixXd J = J_com(u,v,h);
        xi_dot = J*qdot;
        _v_c_rel = xi_dot.block<3,1>(0,0);
        _omega_c_rel = xi_dot.block<3,1>(3,0);
        // std::cout << "v_c_rel: " << _v_c_rel << std::endl;
        // std::cout << "omega_c_rel: " << _omega_c_rel << std::endl;
        // find omega from jacobian
        MatrixXd Jdot_u = J_dot_udot_com(u,v,h);
        MatrixXd Jdot_v = J_dot_vdot_com(u,v,h);
        MatrixXd Jdot = Jdot_u*udot+Jdot_v*vdot;
        // std::cout << "Jdot: " << Jdot << std::endl;  
        xi_ddot = J*qddot+Jdot*qdot;
        _a_c_rel = xi_ddot.block<3,1>(0,0);
        // std::cout << "a_c_rel: " << _a_c_rel << std::endl; 
        _alpha_c_rel = xi_ddot.block<3,1>(3,0);
        // std::cout << "alpha_c_rel: " << _alpha_c_rel << std::endl; 
        // omega_rel, alpha_rel, a_c_rel, v_c_rel, a_e_rel, v_e_rel
    }
    else
    {
        _v_e_rel.setZero();
        _omega_c_rel.setZero();
        _a_e_rel.setZero();
        _alpha_c_rel.setZero();
    }
}

void Grub_rne::set_uv(VectorXd uv, VectorXd uv_dot)
{
    int joint_count, total_joint_count, discrete_joint_count;
    VectorXd phin_vec, l_vec;
    double u,v,phi, h, rho, rhox, rhoy, phin,l,vn,un, u_dot, v_dot;

    // set actual uv for xdot equations
    _uv_real = uv;
    total_joint_count = 0;
    
    // discretize uv for each of the joints
    for (joint_count = 0; joint_count < _real_joint_num; ++joint_count)
    {
        u = uv(joint_count*2);
        v = uv(joint_count*2+1);

        for (discrete_joint_count = 0; discrete_joint_count < 
                            _num_discrete_segments; discrete_joint_count++)
        {
            _uv(total_joint_count*2) = u*_discrete_segments(discrete_joint_count);
            _uv(total_joint_count*2+1) = v*_discrete_segments(discrete_joint_count);
            total_joint_count++;
        }   
    }
    // std::cout << "uv: " << _uv << std::endl;
    // set actual uv_dot
    _uvdot_real = uv_dot;

    total_joint_count = 0;

    // discretize uv_dot for each of the segments
    for (joint_count = 0; joint_count < _real_joint_num; ++joint_count)
    {
        u_dot = uv_dot(joint_count*2);
        v_dot = uv_dot(joint_count*2+1);

        for (discrete_joint_count = 0; discrete_joint_count < 
                            _num_discrete_segments; discrete_joint_count++)
        {
            _uvdot(total_joint_count*2) = u_dot*_discrete_segments(discrete_joint_count);
            _uvdot(total_joint_count*2+1) = v_dot*_discrete_segments(discrete_joint_count);
            total_joint_count++;
        }   
    }
    // std::cout << "uv_dot: " << _uvdot << std::endl;
    update_mat();
}

// getters and setters
MatrixXd Grub_rne::get_M()
{
    return _M;
}

VectorXd Grub_rne::get_Cqdot()
{
    return _Cqdot;
}

VectorXd Grub_rne::get_G()
{
    return _G;
}


// MatrixXd Grub_rne::get_jacobian(double u, double v, double h)
// {
//     MatrixXd Jacobian(6,2);
//     double phi;
//     phi = sqrt(pow(u,2.)+pow(v,2));
//     // test if phi is close to zero
//     if (phi > 10e6)
//     {
//         Jacobian = J(u,v,h);
//     }
//     else
//     {
//         Jacobian = J(u,v,h);
//     }
    
//     return Jacobian;
// }