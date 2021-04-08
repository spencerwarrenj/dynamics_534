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
    // std::cout << "number of discrete segments: " << _num_discrete_segments << std::endl;
    // std::cout << "Constructor: link: " << _links << "\n";
    build_robot();
}

void Grub_rne::discretize_joint(double h, double mass, double r_width, 
                int* joint_count, int* link_count, Matrix3d* prev_rot)
{
    Vector3d rc,re;
    Matrix3d rot, Inertia_tens, Inertia_diag, skewr;
    double *Inertia_h = new double[9];
    double mu, sigma, theta, v_tilde, u_tilde, cphi, sphi, phi, u, v, rho, r;
    int i;
    int discrete_joint_count = 0;
    VectorXd l_vec, mass_vec;

    l_vec = h*_discrete_segments;
    mass_vec = mass*_discrete_segments;
    mu = mass/h;
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
        _mu(*link_count) = mu;
        _r(*link_count) = r_width;
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
        if (abs(phi)<M_PI/24)
        {
            small_inertia(h, mu, r_width, u, v, Inertia_h);
        }
        else
        {
            true_inertia(h, mu, r_width, u, v, Inertia_h);
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
    _K_s.resize(_real_joint_num*2,_real_joint_num*2);
    _K_d.resize(_real_joint_num*2,_real_joint_num*2);
    _m_mat.resize(_rigid_num+_total_joint_num);
    // forward recursion
    _omega_mat.resize(3,_link_num);
    _alpha_mat.resize(3,_link_num);
    _a_c_mat.resize(3,_link_num);
    _a_e_mat.resize(3,_link_num);
    _v_e_mat.resize(3,_link_num);
    // backward recursion
    _f_mat.resize(3,_link_num+1);
    _tau_mat.resize(3,_link_num+1);
    // _qdot_mat.resize(3,_link_num);
    // _qddot_mat.resize(3,_link_num);
    _flexi_link.resize(_total_joint_num+_rigid_num);
    _report_tau.resize(_total_joint_num+_rigid_num);
    _h.resize(_total_joint_num+_rigid_num);
    _mu.resize(_total_joint_num+_rigid_num);
    _r.resize(_total_joint_num+_rigid_num);


    _report_tau.setZero();
    _flexi_link.setZero();
    _h.setZero();
    // _uv.setZero();

    VectorXd vec(_real_joint_num*2);
    vec.setOnes();
    vec = vec;
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
                _mu(link_count) = 0;
                _r(link_count) = 0;
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
                _mu(link_count) = 0;
                _r(link_count) = 0;
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
                _mu(link_count) = 0;
                _r(link_count) = 0;
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
                mass = 2; //1.32636666666667;
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
    _K_s.setIdentity();
    _K_s = _K_s*2;
    _K_d.setIdentity();
    _K_d = _K_d*.08;
    calc_G();
    calc_Cqdot();
    calc_M();
    // std::cout << "G: " << _G << std::endl;
    // std::cout << "Cqdot: " << _Cqdot << std::endl;
    // std::cout << "M: " << _M << std::endl;
    MatrixXd minv = _M.inverse();
    // std::cout << "M_inv: " << minv << std::endl;
    // std::cout << "minv: " << minv << std::endl;
    // std::cout << "tau: " << tau << std::endl;
    // std::cout << "K_s: " << _K_s << std::endl;
    _qdd = minv*(-_Cqdot-_G+tau); //-_K_s*_uv_real-_K_d*_uvdot_real);

    return _qdd;
}

VectorXd Grub_rne::get_real_taus()
{
    Matrix3d rot;
    VectorXd tau_block(_real_joint_num*2);

    int tau_report_count = 0;

    for (int link_count = 0; link_count < _link_num; link_count++)
    {
    // std::cout << "tau_mat: " << _tau_mat << std::endl;
    // std::cout << "i: " << i << "\n";
    // std::cout << "set G" << "\n";
        if (_report_tau(link_count)==1)
        {
            // std::cout << "_rot_mat: " << _rot_mat << "\n";
            rot = _rot_mat.block<3,3>(0,link_count*3);
            // std::cout << "rot: " << rot << "\n";
            Vector3d dum = rot*_tau_mat.block<3,1>(0,link_count);
            
            // std::cout << "link_count: " << link_count << "\n";
            // std::cout << "dum: " << dum << "\n";
            tau_block.block<2,1>(link_count/_num_discrete_segments*2,0) = dum.block<2,1>(0,0);
            // std::cout << "tau_block: " << tau_block << std::endl;
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
    // For EL for single joint
    // _G << -0.19134463,-0.19134463; ///////////////Testing//////////////////
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
    // for EL for single joint
    // _Cqdot << -0.00019864,-0.00019864; ///////////////////Testing

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

    
    for (int j = 0; j < _real_joint_num; j++)
    {
        for (int i = 0; i < 2; i++)
        {
            // set u_ddot to 1
            // std::cout << "_discrete_segments: " << _discrete_segments << "\n";
            uvddot.setZero();
            for (int k=0; k < _num_discrete_segments; k++)
            {
                // std::cout << "uvddot: " << uvddot << "\n";
                uvddot(j*(2*_num_discrete_segments)+i+k*2) = _discrete_segments(k);
            }
            
            // std::cout << "uvddot: " << uvddot << "\n";
            get_tau(g,uvdot,uvddot);
            // std::cout << "tau: " << _tau_mat <<"\n";
            // R01,R12,R23,R34,R45,R56
            _M.col(j*2+i) = get_real_taus();
        }
    }
    // from EL for single joint
    // _M << 5.34455854e-03,-1.55293806e-05,-1.55293806e-05,  5.34455854e-03; %%%%%////
}

// This function runs through the forward and backward recursion to get dynamics
void Grub_rne::get_tau(Vector3d g, VectorXd uvdot, VectorXd uvddot)
{
    // std::cout << "got here " <<"\n";
    Vector3d prev_alpha, prev_omega, prev_a_e, prev_v_e, rc, re, omega,
            alpha, a_c, a_e, v_e, f, tau, prev_f, prev_tau, qdot, qddot, Inertial_term;
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
            prev_v_e.setZero();
        }
        else
        {
            // prev_alpha = _alpha_mat.block<3,1>(0,i-1);
            // prev_omega = _omega_mat.block<3,1>(0,i-1);
            prev_a_e = _a_e_mat.block<3,1>(0,i-1);
            prev_v_e = _v_e_mat.block<3,1>(0,i-1);
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
        // std::cout << "Rot" << rot << "\n";
        // std::cout << "Start Forward recursion" << "\n";
        
        // forward recursion for kinematics
        // set omega_c_rel, alpha_c_rel, omega_e_rel, alpha_e_rel, 
        //              a_c_rel, v_c_rel, a_e_rel, v_e_rel
        get_e_rel(i,uvdot,uvddot);
        get_c_rel(i,uvdot,uvddot);
        // std::cout << "omega_rel: " << omega_rel << "\n";
        // alpha_rel = get_alpha_rel(i,uvdot,uvddot);
        // std::cout << "alpha_rel: " << alpha_rel << "\n";
        // std::cout << "prev_alpha: " << prev_alpha << "\n";
        // std::cout << "prev_a_e: " << prev_a_e << "\n";
        // std::cout << "a_e_rel: " << _a_e_rel << "\n";
        // std::cout << "re: " << re << "\n";
        v_e = rott*(_v_e_rel+prev_v_e+prev_omega.cross(re));
        a_c = rott*(prev_a_e+_a_c_rel+prev_alpha.cross(rc)+
                prev_omega.cross(prev_omega.cross(rc))+
                2*prev_omega.cross(_v_c_rel));
        a_e = rott*(prev_a_e+_a_e_rel+prev_alpha.cross(re)+
                prev_omega.cross(prev_omega.cross(re))+
                2*prev_omega.cross(_v_e_rel));
        omega = rott*prev_omega+rott*_omega_c_rel;
        alpha = rott*prev_alpha+rott*_alpha_c_rel+omega.cross(rott*_omega_c_rel);

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
        _v_e_mat.block<3,1>(0,i) = v_e;

        prev_omega = omega;
        prev_alpha = alpha;
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
        if (i == 0)
        {
            prev_omega.setZero();
            prev_alpha.setZero();
            prev_v_e.setZero();
            prev_a_e.setZero();
        }
        else
        {
            prev_omega = _omega_mat.block<3,1>(0,i-1);
            prev_alpha = _alpha_mat.block<3,1>(0,i-1);
            prev_v_e = _v_e_mat.block<3,1>(0,i-1);
            prev_a_e = _a_e_mat.block<3,1>(0,i-1);
            // std::cout << "prev_omega: " << prev_omega << "\n";
            // std::cout << "prev_alpha: " << prev_alpha << "\n";
            // std::cout << "prev_a_e: " << prev_a_e << "\n"; 
            // std::cout << "prev_v_e: " << prev_v_e << "\n"; 
        }
        omega = _omega_mat.block<3,1>(0,i);
        alpha = _alpha_mat.block<3,1>(0,i);
        a_c = _a_c_mat.block<3,1>(0,i);
        a_e = _a_e_mat.block<3,1>(0,i);
        // std::cout << "omega: " << omega << "\n";
        // std::cout << "alpha: " << alpha << "\n";
        // std::cout << "a_c: " << a_c << "\n"; 
        // std::cout << "a_e: " << a_e << "\n";
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
        // std::cout << "rotf: " << rotf << std::endl;
        // std::cout << "rot: " << rot <<"\n";
        // // std::cout << "prev_f: " << prev_f <<"\n";
        // std::cout << "m: " << m << std::endl;
        // std::cout << "rotground: " << rotground <<"\n";
        // std::cout << "rott: " << rott <<"\n";
        // std::cout << "rotb: " << rott <<"\n";
        // std::cout << "g: " << g <<"\n";
        // std::cout << "g_term: " << m*rotground*g << std::endl;
        // std::cout << "a_c: " << a_c <<"\n";
        f = rotf*prev_f-m*rotground*g+m*a_c;
        // std::cout << "f: " << f << std::endl;
        // std::cout << "rot*f: " << rot*f << std::endl;
        Inertial_term = get_inertial_term(i,uvdot,uvddot,
                                prev_alpha,prev_omega,alpha,omega,
                                rott,Inertia_tensor);
        // std::cout << "prev_tau: " << prev_tau << std::endl;
        // std::cout << "omega: " << omega << std::endl;
        // std::cout << "alpha: " << alpha << std::endl;
        Vector3d rigid_inertia = Inertia_tensor*alpha+omega.cross(Inertia_tensor*omega);
        // std::cout << "Inertial_term: " << Inertial_term << std::endl;
        // std::cout << "rigid_Inertial_term: " << rigid_inertia << std::endl;
        // std::cout << "rc: " << rc << std::endl;
        // std::cout << "rotb: " << rotb << std::endl;
        // std::cout << "rotb*rc: " << rotb*rc << std::endl;
        tau = rotf*prev_tau-f.cross(rotb*rc)+
                    (rotf*prev_f).cross(rotb*(rc-re))+Inertial_term;
        // std::cout << "rot*tau: " << rot*tau << std::endl;
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

Vector3d Grub_rne::get_inertial_term(int link_num, VectorXd uvdot, 
            VectorXd uvddot, Vector3d alpha_prev, Vector3d omega_prev,
            Vector3d alpha, Vector3d omega, Matrix3d rott,
            Matrix3d Inertia_tensor)
{
    Vector3d Inertial;
    double u, v, h, mu, r, udot,vdot, uddot, vddot;
    int i;
    
    // // test if it is a flexible joint
    // if (_flexi_link(link_num)==1)
    // {
    //     int joint_num=0;
    //     for (i=0; i<link_num; i++)
    //     {
    //         if (_flexi_link(i)==1)
    //         {
    //             joint_num++;
    //         }
    //     }

    //     u = _uv(joint_num*2);
    //     v = _uv(joint_num*2+1);
    //     h = _h(link_num);
    //     mu = _mu(link_num);
    //     r = _r(link_num);
    //     udot = uvdot(joint_num*2);
    //     vdot = uvdot(joint_num*2+1);
    //     uddot = uvddot(joint_num*2);
    //     vddot = uvddot(joint_num*2+1);
    //     // std::cout << "_mu: " << _mu << std::endl;
    //     // std::cout << "_r: " << _r << std::endl;
    //     // std::cout << "mu: " << mu << std::endl;
    //     // std::cout << "r: " << r << std::endl;
    //     // std::cout << "u: " << u << std::endl;
    //     // std::cout << "v: " << v << std::endl;
    //     // std::cout << "h: " << h << std::endl;
    //     // std::cout << "uvdot: " << uvdot << std::endl;
    //     // std::cout << "udot: " << udot << std::endl;
    //     // std::cout << "vdot: " << vdot << std::endl;
    //     // std::cout << "uddot: " << uddot << std::endl;
    //     // std::cout << "vddot: " << vddot << std::endl;
    //     // std::cout << "alpha_prev: " << alpha_prev << std::endl;
    //     // std::cout << "omega_prev: " << omega_prev << std::endl;
    //     Inertial = full_deriv_ang_mom(u,v,h,mu,r,udot,vdot,uddot,vddot,
    //                                     alpha_prev,omega_prev);
    //     // std::cout << "Inertial: " << Inertial << std::endl;
    //     Inertial = rott*Inertial;
    // }
    // else
    // {
        Inertial = Inertia_tensor*alpha+omega.cross(Inertia_tensor*omega);
    // }
    return Inertial;
}

Vector3d Grub_rne::full_deriv_ang_mom(double u, double v, double h, double mu, double r,
        double u_dot, double v_dot, double u_ddot, double v_ddot,
        Vector3d alpha_prev, Vector3d omega_prev)
{
    // updated: 04/06/2021
    double alphax, alphay, alphaz, w_x, w_y, w_z;
    alphax = alpha_prev(0);
    alphay = alpha_prev(1);
    alphaz = alpha_prev(2);
    w_x = omega_prev(0);
    w_y = omega_prev(1);
    w_z = omega_prev(2);

    const double x0 = pow(u, 3);
    const double x1 = v_dot*x0;
    const double x2 = pow(h, 2);
    const double x3 = pow(u, 2);
    const double x4 = pow(v, 2);
    const double x5 = x3 + x4;
    const double x6 = sqrt(x5);
    const double x7 = (1.0/2.0)*x6;
    const double x8 = x2*x7;
    const double x9 = x1*x8;
    const double x10 = x3*x6;
    const double x11 = x4*x6;
    const double x12 = x10 + x11;
    const double x13 = 1.0/x12;
    const double x14 = pow(u, 4);
    const double x15 = h*x14;
    const double x16 = pow(v, 4);
    const double x17 = h*x16;
    const double x18 = x3*x4;
    const double x19 = 2*x18;
    const double x20 = h*x19 + x15 + x17;
    const double x21 = pow(r, 2);
    const double x22 = mu*x21;
    const double x23 = x22/x20;
    const double x24 = x13*x23;
    const double x25 = x10*x2;
    const double x26 = u_dot*x25;
    const double x27 = (1.0/2.0)*x23;
    const double x28 = v*x27;
    const double x29 = x26*x28;
    const double x30 = 4*x10;
    const double x31 = 4*x11;
    const double x32 = x30 + x31;
    const double x33 = 1.0/x32;
    const double x34 = x23*x33;
    const double x35 = pow(v, 3);
    const double x36 = u_dot*x35;
    const double x37 = x34*x8;
    const double x38 = x2*x6;
    const double x39 = 9*x11;
    const double x40 = 9*x10 + x39;
    const double x41 = 1.0/x40;
    const double x42 = x23*x41;
    const double x43 = (9.0/2.0)*x42;
    const double x44 = x38*x43;
    const double x45 = x11*x2;
    const double x46 = v_dot*x45;
    const double x47 = u*x27;
    const double x48 = x33*x47;
    const double x49 = u*x43;
    const double x50 = pow(h, 3);
    const double x51 = mu*x50;
    const double x52 = u_dot*x51;
    const double x53 = x52*x6;
    const double x54 = 4*x2;
    const double x55 = x14*x54;
    const double x56 = x16*x54;
    const double x57 = 8*x18*x2;
    const double x58 = x55 + x56 + x57;
    const double x59 = 1.0/x58;
    const double x60 = x21*x33;
    const double x61 = x59*x60;
    const double x62 = x35*x61;
    const double x63 = x51*x6;
    const double x64 = x61*x63;
    const double x65 = x10*x52;
    const double x66 = v*x61;
    const double x67 = v_dot*x51;
    const double x68 = x11*x67;
    const double x69 = u*x61;
    const double x70 = h*x22;
    const double x71 = (1.0/4.0)*x70;
    const double x72 = w_x*x71;
    const double x73 = 1.0/x5;
    const double x74 = u*x73;
    const double x75 = x72*x74;
    const double x76 = w_y*x71;
    const double x77 = x73*x76;
    const double x78 = v*x77;
    const double x79 = (1.0/2.0)*w_z;
    const double x80 = x70*x79;
    const double x81 = sin(x6);
    const double x82 = 1.0/x6;
    const double x83 = x81*x82;
    const double x84 = cos(x6);
    const double x85 = x6/h != 0;
    const double x86 = x2*x3;
    const double x87 = x13*x81;
    const double x88 = x2*x4;
    const double x89 = x13*x84;
    const double x90 = x38*x89 + x86*x87 + x87*x88;
    const double x91 = x27*x90;
    const double x92 = u_dot*x3;
    const double x93 = v*x92;
    const double x94 = pow(x81, 3);
    const double x95 = x41*x94;
    const double x96 = 3*x86;
    const double x97 = x88*x95;
    const double x98 = pow(x84, 3);
    const double x99 = x41*x98;
    const double x100 = 2*x38;
    const double x101 = pow(x81, 2);
    const double x102 = 3*x101;
    const double x103 = x41*x84;
    const double x104 = x103*x38;
    const double x105 = x100*x99 + x102*x104 + x95*x96 + 3*x97;
    const double x106 = x27*x36;
    const double x107 = x4*x47;
    const double x108 = v_dot*x107;
    const double x109 = 6*x86;
    const double x110 = 7*x38;
    const double x111 = 9*x86;
    const double x112 = pow(x84, 2);
    const double x113 = x41*x81;
    const double x114 = x112*x113;
    const double x115 = 9*x88;
    const double x116 = 6*x101;
    const double x117 = x104*x116 + x109*x95 + x110*x99 + x111*x114 + x114*x115 + 6*x97;
    const double x118 = x112*x33;
    const double x119 = x118*x38;
    const double x120 = x81*x84;
    const double x121 = x120*x33;
    const double x122 = 2*x121;
    const double x123 = x122*x86;
    const double x124 = x122*x88;
    const double x125 = x101*x33;
    const double x126 = x118*x25 + x118*x45 + x125*x25 + x125*x45;
    const double x127 = x119 + x123 + x124 + x126;
    const double x128 = x127*x27;
    const double x129 = v_dot*x127;
    const double x130 = -x119 - x123 - x124 + x126;
    const double x131 = x59*x70;
    const double x132 = x130*x131;
    const double x133 = u*x4;
    const double x134 = v_dot*x132;
    const double x135 = -x1*x128 - x1*x132 - x1*x64 + x1*x91 - x105*x106 + x105*x108 - x106*x117 - x107*x129 + x108*x117 + x128*x36 + x128*x93 + x13*x29 + x132*x36 + x132*x93 - x133*x134 - x24*x9 - x29*x33 + x34*x9 - x36*x37 + x36*x44 + x46*x48 - x46*x49 + x53*x62 + x65*x66 - x68*x69 - x91*x93 + ((x85) ? (
        x75*x84 - x75 - x78*x84 + x78 + x80*x83
    )
    : (
        x80
    ));
    const double x136 = x14 + x16 + x19;
    const double x137 = 1.0/x136;
    const double x138 = x137*x51;
    const double x139 = u_dot*x138;
    const double x140 = v*x139;
    const double x141 = v_dot*x138;
    const double x142 = u*x141;
    const double x143 = w_z*x138;
    const double x144 = x143*x3;
    const double x145 = x143*x4;
    const double x146 = pow(u, 6);
    const double x147 = pow(v, 6);
    const double x148 = 3*x16;
    const double x149 = x148*x3;
    const double x150 = 3*x14;
    const double x151 = x150*x4;
    const double x152 = x146 + x147 + x149 + x151;
    const double x153 = x51/x152;
    const double x154 = x1*x153;
    const double x155 = w_z*x153;
    const double x156 = u_dot*x153;
    const double x157 = x156*x35;
    const double x158 = x153*x4;
    const double x159 = v_dot*x158;
    const double x160 = u*x159;
    const double x161 = x156*x3;
    const double x162 = v*x161;
    const double x163 = u*x153;
    const double x164 = x101*x163;
    const double x165 = v_dot*x164;
    const double x166 = v*x101;
    const double x167 = x156*x166;
    const double x168 = x101*x155;
    const double x169 = 2*x6;
    const double x170 = x169*x81;
    const double x171 = v_dot*x163;
    const double x172 = x170*x171;
    const double x173 = v*x170;
    const double x174 = x156*x173;
    const double x175 = x10*x81;
    const double x176 = 2*x155;
    const double x177 = x11*x81;
    const double x178 = x14*x6;
    const double x179 = x16*x6;
    const double x180 = 2*x4;
    const double x181 = x10*x180 + x178 + x179;
    const double x182 = 1.0/x181;
    const double x183 = x182*x81;
    const double x184 = w_x*x51;
    const double x185 = v*x184;
    const double x186 = x183*x185;
    const double x187 = x183*x67;
    const double x188 = 2*u;
    const double x189 = x187*x188;
    const double x190 = x183*x51;
    const double x191 = u*x190;
    const double x192 = w_y*x191;
    const double x193 = u_dot*x190;
    const double x194 = 2*v;
    const double x195 = x193*x194;
    const double x196 = w_z*x190;
    const double x197 = 2*x3;
    const double x198 = x196*x197;
    const double x199 = x180*x196;
    const double x200 = x120*x182;
    const double x201 = x200*x51;
    const double x202 = u*x201;
    const double x203 = x185*x200;
    const double x204 = w_y*x138;
    const double x205 = u*x204;
    const double x206 = 2*x83;
    const double x207 = w_x*x138;
    const double x208 = v*x207;
    const double x209 = x120*x82;
    const double x210 = 2*x209;
    const double x211 = x51*x82;
    const double x212 = w_y*x211;
    const double x213 = (1.0/2.0)*x13;
    const double x214 = x101*x213;
    const double x215 = x212*x214;
    const double x216 = x185*x82;
    const double x217 = (1.0/2.0)*x120;
    const double x218 = x217 + x7;
    const double x219 = x218*x82;
    const double x220 = x101*x67;
    const double x221 = x182*x82;
    const double x222 = x188*x221;
    const double x223 = x220*x222;
    const double x224 = x182*x211;
    const double x225 = x166*x224;
    const double x226 = 2*u_dot;
    const double x227 = x101*x224;
    const double x228 = w_z*x227;
    const double x229 = x208*x84;
    const double x230 = x205*x84;
    const double x231 = w_y*x202 + x14*x155 + x140 - x142 - x144 - x145 + x154 + x155*x16 + x155*x19 - x157 + x160 - x162 + x165 - x167 + x168*x3 + x168*x4 - x172 + x174 - x175*x176 - x176*x177 + x186 + x189 - x192 - x195 + x198 + x199 - x203 + ((x85) ? (
        -u*x215 - x140*x219 + x142*x219 + x144*x219 + x145*x219 - x180*x228 - x197*x228 + x205*x206 - x205*x210 - x206*x208 + x208*x210 + x214*x216 - x223 + x225*x226
    )
    : (
        -x140 + x142 + x144 + x145 - x189 + x195 - x198 - x199 + x205 - x208 + x229 - x230
    ));
    const double x232 = v*x75;
    const double x233 = x4*x77;
    const double x234 = x73*x80;
    const double x235 = v*x234;
    const double x236 = u*v;
    const double x237 = pow(x5, -3.0/2.0);
    const double x238 = x237*x81;
    const double x239 = x236*x238;
    const double x240 = x4*x76;
    const double x241 = v*v_dot;
    const double x242 = u*x241;
    const double x243 = x38*x87 - x86*x89 - x88*x89;
    const double x244 = x25 + x45;
    const double x245 = 1.0/x244;
    const double x246 = x245*x71;
    const double x247 = x243*x246;
    const double x248 = u_dot*x4;
    const double x249 = h*x4;
    const double x250 = x15*x169 + x169*x17 + x249*x30;
    const double x251 = x22/x250;
    const double x252 = x243*x251;
    const double x253 = x0*x241;
    const double x254 = u_dot*x18;
    const double x255 = x178*x54 + x179*x54 + 8*x25*x4;
    const double x256 = 1.0/x255;
    const double x257 = x22*x256;
    const double x258 = x17*x257;
    const double x259 = u_dot*x258;
    const double x260 = v_dot*x35;
    const double x261 = u*x260;
    const double x262 = x256*x70;
    const double x263 = x243*x262;
    const double x264 = x88*x99;
    const double x265 = 3*x112;
    const double x266 = x38*x81;
    const double x267 = x265*x266;
    const double x268 = x100*x95 - 3*x264 + x267*x41 - x96*x99;
    const double x269 = x251*x268;
    const double x270 = u_dot*x16;
    const double x271 = -x118*x86 - x118*x88 + x121*x38 + x125*x86 + x125*x88;
    const double x272 = x251*x270;
    const double x273 = x251*x271;
    const double x274 = u_dot*x271;
    const double x275 = x262*x271;
    const double x276 = x18*x262;
    const double x277 = x101*x103;
    const double x278 = x112*x38;
    const double x279 = 6*x278;
    const double x280 = -x109*x99 + x110*x95 - x111*x277 + x113*x279 - x115*x277 - 6*x264;
    const double x281 = x251*x280;
    const double x282 = 12*x14 + 12*x16 + 24*x18;
    const double x283 = x50/x282;
    const double x284 = x283*x3;
    const double x285 = x265*x284;
    const double x286 = x283*x4;
    const double x287 = x265*x286;
    const double x288 = x102*x284;
    const double x289 = x102*x286;
    const double x290 = x120*x283;
    const double x291 = 3*x6;
    const double x292 = x290*x291;
    const double x293 = 6*x290;
    const double x294 = x10*x293;
    const double x295 = x11*x293;
    const double x296 = 2*x283;
    const double x297 = x14*x296;
    const double x298 = x16*x296;
    const double x299 = 4*x18;
    const double x300 = x112*x283;
    const double x301 = x101*x283;
    const double x302 = x101*x297 + x101*x298 + x112*x297 + x112*x298 + x299*x300 + x299*x301;
    const double x303 = -x285 - x287 + x288 + x289 + x292 - x294 - x295 + x302;
    const double x304 = x22*x245;
    const double x305 = x303*x304;
    const double x306 = (1.0/4.0)*x305;
    const double x307 = x242*x6;
    const double x308 = x285 + x287 - x288 - x289 - x292 + x294 + x295 + x302;
    const double x309 = x304*x308;
    const double x310 = (1.0/4.0)*x309;
    const double x311 = u_dot*x10;
    const double x312 = (1.0/4.0)*x311;
    const double x313 = x232 - x233 - x242*x247 - x243*x259 + x247*x248 - x252*x253 + x252*x254 + x253*x263 + x253*x273 - x253*x275 - x254*x263 - x254*x273 + x258*x274 + x261*x263 - x261*x269 + x261*x273 - x261*x275 - x261*x281 + x269*x270 - x271*x272 + x272*x280 + x274*x276 + x306*x307 + x306*x311 + x307*x310 + x309*x312 + x76 + ((x85) ? (
        x235*x84 - x235 + x238*x240 - x239*x72
    )
    : (
        -x232 + x233
    ));
    const double x314 = x153*x16;
    const double x315 = x204*x4;
    const double x316 = v_dot*x153;
    const double x317 = x180*x316;
    const double x318 = 2*x84;
    const double x319 = x163*x35;
    const double x320 = v*x153;
    const double x321 = x0*x320;
    const double x322 = mu*x137;
    const double x323 = x249*x322;
    const double x324 = x130*x323;
    const double x325 = x187*x4;
    const double x326 = x4*x67;
    const double x327 = x10*x148 + x11*x150 + x146*x6 + x147*x6;
    const double x328 = 1.0/x327;
    const double x329 = x328*x81;
    const double x330 = x200*x67;
    const double x331 = x153*x18;
    const double x332 = x101*x158;
    const double x333 = u*x208;
    const double x334 = v_dot*x243;
    const double x335 = mu*x182;
    const double x336 = x249*x335;
    const double x337 = x334*x336;
    const double x338 = u*u_dot;
    const double x339 = 2*x338;
    const double x340 = x320*x339;
    const double x341 = 2*x328;
    const double x342 = x180*x190;
    const double x343 = w_y*x342;
    const double x344 = h*x322;
    const double x345 = x127*x344;
    const double x346 = v*x338;
    const double x347 = x130*x344;
    const double x348 = v*x190;
    const double x349 = x338*x348;
    const double x350 = x320*x338;
    const double x351 = x329*x51;
    const double x352 = v*x164;
    const double x353 = h*x335;
    const double x354 = x243*x353;
    const double x355 = x346*x354;
    const double x356 = v_dot*x90;
    const double x357 = mu*x183;
    const double x358 = x249*x357;
    const double x359 = x120*x328;
    const double x360 = x318*x328;
    const double x361 = w_y*x153;
    const double x362 = 2*x177;
    const double x363 = x328*x51;
    const double x364 = x339*x6;
    const double x365 = v*x364;
    const double x366 = x363*x365;
    const double x367 = x186*x188;
    const double x368 = x182*x87;
    const double x369 = x357*x90;
    const double x370 = h*x369;
    const double x371 = x6*x81;
    const double x372 = x359*x51;
    const double x373 = v*x163;
    const double x374 = x170*x373;
    const double x375 = x368*x63;
    const double x376 = -x217 + x7;
    const double x377 = x237*x51;
    const double x378 = x376*x377;
    const double x379 = 2*x13;
    const double x380 = x212*x379;
    const double x381 = x221*x326;
    const double x382 = 4*x89;
    const double x383 = x220*x82;
    const double x384 = x213*x383;
    const double x385 = x224*x346;
    const double x386 = 2*x141;
    const double x387 = x182*x212;
    const double x388 = x101*x180;
    const double x389 = u*x143;
    const double x390 = x13*x79;
    const double x391 = x101*x211;
    const double x392 = x390*x391;
    const double x393 = (1.0/2.0)*x101;
    const double x394 = x101*x185;
    const double x395 = (1.0/2.0)*x225;
    const double x396 = x159*x84;
    const double x397 = x350*x84;
    const double x398 = x396 + x397;
    const double x399 = x141*x84;
    const double x400 = v_dot*x324 + w_x*x319 + w_x*x321 + w_x*x352 - w_x*x374 + w_y*x314 + w_y*x331 + w_y*x332 - w_z*x191 + w_z*x202 + x101*x159 + x101*x350 + x112*x159 + x112*x204 + x112*x350 + x129*x323 - x177*x316 - x187 - x204*x318 + x204 - x315 - x317 + x325 + x326*x329 - x326*x359 + x330 - x333 + x337*x84 - x337 - x340 + x341*x68 + x343 + x345*x346 + x346*x347 + x346*x351 - x346*x370 - x346*x372 + x346*x375 + x349 - x350*x371 + x355*x84 - x355 - x356*x358 - x360*x68 - x361*x362 - x366*x84 + x366 + x367 + x368*x68 + x398 + ((x85) ? (
        -u*x392 + w_y*x378 - x112*x380 - x112*x381 - x112*x385 + x141*x206 - x159 + x206*x389 + x209*x317 + x209*x340 - x209*x386 - x210*x389 + x212*x382 + x219*x315 + x219*x333 - x222*x394 - x317*x83 - x338*x395 - x340*x83 - x350 - x380 - x381*x393 + x381 - x384 + x385 - x387*x388 + x398
    )
    : (
        x141 + x315 - x325 + x333 - x343 - x349 - x367 - x389*x84 + x389 - x399
    ));
    const double x401 = alphax*x138;
    const double x402 = alphax*x71;
    const double x403 = x14*x153;
    const double x404 = u_dot*x0;
    const double x405 = 4*x404;
    const double x406 = 4*x260;
    const double x407 = x338*x4;
    const double x408 = 4*x407;
    const double x409 = x241*x3;
    const double x410 = 4*x409;
    const double x411 = (-x405 - x406 - x408 - x410)/pow(x136, 2);
    const double x412 = x184*x411;
    const double x413 = x3*x401;
    const double x414 = u_ddot*x190;
    const double x415 = pow(u_dot, 2);
    const double x416 = x163*x415;
    const double x417 = 4*x416;
    const double x418 = pow(v_dot, 2);
    const double x419 = x163*x418;
    const double x420 = 2*x419;
    const double x421 = u_ddot*x153;
    const double x422 = x3*x421;
    const double x423 = 2*x422;
    const double x424 = alphax*x3;
    const double x425 = x101*x153;
    const double x426 = u*x418;
    const double x427 = x191*x418;
    const double x428 = u_ddot*x3;
    const double x429 = x3*x414;
    const double x430 = 8*x6;
    const double x431 = x241 + x338;
    const double x432 = x431*x82;
    const double x433 = x3*x432;
    const double x434 = x4*x432;
    const double x435 = (-x241*x430 - x338*x430 - 4*x433 - 4*x434)/pow(x32, 2);
    const double x436 = x278*x435;
    const double x437 = x2*x432;
    const double x438 = x118*x437;
    const double x439 = x2*x431;
    const double x440 = x122*x439;
    const double x441 = x2*x433;
    const double x442 = x125*x441;
    const double x443 = x2*x434;
    const double x444 = x125*x443;
    const double x445 = x120*x435;
    const double x446 = 2*x445;
    const double x447 = x446*x86;
    const double x448 = x446*x88;
    const double x449 = x118*x441;
    const double x450 = x118*x443;
    const double x451 = x121*x54;
    const double x452 = x338*x451;
    const double x453 = x241*x451;
    const double x454 = x112*x435;
    const double x455 = x101*x435;
    const double x456 = x339*x38;
    const double x457 = 2*x241;
    const double x458 = x38*x457;
    const double x459 = x119*x339 + x119*x457 + x125*x456 + x125*x458 + x25*x454 + x25*x455 + x45*x454 + x45*x455;
    const double x460 = x436 + x438 - x440 - x442 - x444 + x447 + x448 + 3*x449 + 3*x450 + x452 + x453 + x459;
    const double x461 = x344*x460;
    const double x462 = -x436 - x438 + x440 + 3*x442 + 3*x444 - x447 - x448 - x449 - x450 - x452 - x453 + x459;
    const double x463 = x344*x462;
    const double x464 = x14*x251;
    const double x465 = v_ddot*x243;
    const double x466 = x2*x89;
    const double x467 = x457*x6;
    const double x468 = (-x364 - x433 - x434 - x467)/pow(x12, 2);
    const double x469 = x468*x84;
    const double x470 = x266*x468 - x339*x466 + x437*x87 + x439*x89 + x441*x87 + x443*x87 - x457*x466 - x469*x86 - x469*x88;
    const double x471 = v_dot*x464;
    const double x472 = 6*pow(u, 5)*u_dot;
    const double x473 = 6*pow(v, 5)*v_dot;
    const double x474 = 6*x338;
    const double x475 = 6*x241;
    const double x476 = 12*x260;
    const double x477 = 12*x404;
    const double x478 = (-x14*x475 - x16*x474 - x3*x476 - x4*x477 - x472 - x473)/pow(x152, 2);
    const double x479 = x14*x478;
    const double x480 = x35*x415;
    const double x481 = alphay*x138;
    const double x482 = u*x481;
    const double x483 = v*x482;
    const double x484 = v_dot*x205;
    const double x485 = x243*x426;
    const double x486 = x353*x485;
    const double x487 = x354*x428;
    const double x488 = u_dot*x204;
    const double x489 = v*x488;
    const double x490 = x353*x470;
    const double x491 = x490*x92;
    const double x492 = (-x14*x432 - x16*x432 - x19*x432 - x241*x30 - x31*x338 - x405*x6 - x406*x6)/pow(x181, 2);
    const double x493 = x492*x81;
    const double x494 = x493*x52;
    const double x495 = v_ddot*x271;
    const double x496 = x125*x2;
    const double x497 = x118*x2;
    const double x498 = x118*x439 + x121*x437 - x125*x439 + x339*x496 - x339*x497 + x38*x445 + x433*x451 + x434*x451 - x454*x86 - x454*x88 + x455*x86 + x455*x88 + x457*x496 - x457*x497;
    const double x499 = x3*x412;
    const double x500 = x207*x339;
    const double x501 = v_ddot*x373;
    const double x502 = 2*x501;
    const double x503 = x188*x415;
    const double x504 = x156*x457;
    const double x505 = x478*x51;
    const double x506 = x505*x92;
    const double x507 = 2*x424;
    const double x508 = x190*x507;
    const double x509 = 2*x415;
    const double x510 = x191*x509;
    const double x511 = 2*x416;
    const double x512 = x169*x363;
    const double x513 = x426*x512;
    const double x514 = x10*x341;
    const double x515 = u_ddot*x51;
    const double x516 = u*x415;
    const double x517 = 4*x6;
    const double x518 = x363*x517;
    const double x519 = x516*x518;
    const double x520 = w_x*x153;
    const double x521 = x3*x73;
    const double x522 = x402*x521;
    const double x523 = u*x35;
    const double x524 = u_ddot*x273;
    const double x525 = x251*x498;
    const double x526 = x338*x35;
    const double x527 = v_ddot*x236;
    const double x528 = v_ddot*x190;
    const double x529 = x236*x528;
    const double x530 = v_ddot*x164;
    const double x531 = w_y*x505;
    const double x532 = v*x0;
    const double x533 = u_dot*x241;
    const double x534 = x193*x241;
    const double x535 = x156*x241;
    const double x536 = v*x404;
    const double x537 = h*x92;
    const double x538 = mu*x411;
    const double x539 = x537*x538;
    const double x540 = x51*x92;
    const double x541 = x493*x540;
    const double x542 = (-x10*x476 - x11*x477 - x146*x432 - x147*x432 - x149*x432 - x151*x432 - x178*x475 - x179*x474 - x472*x6 - x473*x6)/pow(x327, 2);
    const double x543 = x542*x81;
    const double x544 = x182*x432;
    const double x545 = x112*x52;
    const double x546 = x120*x492;
    const double x547 = w_z*x51;
    const double x548 = x493*x547;
    const double x549 = x15*x257;
    const double x550 = v_ddot*x18;
    const double x551 = v_dot*x549;
    const double x552 = x2*x99;
    const double x553 = x439*x99;
    const double x554 = x437*x95;
    const double x555 = 9*x433;
    const double x556 = x114*x2;
    const double x557 = 9*x434;
    const double x558 = 18*x6;
    const double x559 = (-x241*x558 - x338*x558 - x555 - x557)/pow(x40, 2);
    const double x560 = x559*x98;
    const double x561 = x560*x88;
    const double x562 = x559*x94;
    const double x563 = x100*x562 + x113*x265*x437 + x267*x559 - x474*x552 - x475*x552 + 3*x553 + 2*x554 + x555*x556 + x556*x557 - x560*x96 - 3*x561;
    const double x564 = x251*x563;
    const double x565 = v_dot*x18;
    const double x566 = 12*x552;
    const double x567 = 18*x2;
    const double x568 = x277*x567;
    const double x569 = 9*x439;
    const double x570 = x2*x95;
    const double x571 = x559*x84;
    const double x572 = x101*x571;
    const double x573 = x559*x81;
    const double x574 = -x109*x560 + x110*x562 - x111*x572 + 6*x114*x437 - x115*x572 - x241*x566 - x241*x568 + x277*x569 + x279*x573 - x338*x566 - x338*x568 + 6*x553 + 7*x554 + x555*x570 + x557*x570 - 6*x561;
    const double x575 = x251*x574;
    const double x576 = x22*x334;
    const double x577 = h*x430;
    const double x578 = 8*h;
    const double x579 = 2*x432;
    const double x580 = h*x432;
    const double x581 = (-x10*x241*x578 - x11*x338*x578 - x15*x579 - x17*x579 - x260*x577 - x299*x580 - x404*x577)/pow(x250, 2);
    const double x582 = x14*x581;
    const double x583 = x184*x478;
    const double x584 = x101*x3;
    const double x585 = x262*x480;
    const double x586 = v*x418;
    const double x587 = x586*x7;
    const double x588 = v*x201;
    const double x589 = u_ddot*x523;
    const double x590 = x236*x465;
    const double x591 = x353*x590;
    const double x592 = x242*x490;
    const double x593 = x411*x51;
    const double x594 = w_y*x593;
    const double x595 = u*x594;
    const double x596 = v*x595;
    const double x597 = u_ddot*x252;
    const double x598 = x354*x533;
    const double x599 = x251*x470;
    const double x600 = x2*x87;
    const double x601 = x468*x81;
    const double x602 = x339*x600 + x38*x469 + x437*x89 - x439*x87 + x441*x89 + x443*x89 + x457*x600 + x601*x86 + x601*x88;
    const double x603 = x357*x602;
    const double x604 = mu*x492;
    const double x605 = x243*x604;
    const double x606 = x537*x605;
    const double x607 = x161*x84;
    const double x608 = x52*x544;
    const double x609 = v_dot*x470;
    const double x610 = x22*x271;
    const double x611 = v_dot*x610;
    const double x612 = 3*x404;
    const double x613 = v_dot*x273;
    const double x614 = v*x415;
    const double x615 = 3*x3*x614;
    const double x616 = 2*x175;
    const double x617 = u*x457;
    const double x618 = x318*x363;
    const double x619 = u_ddot*x10;
    const double x620 = x194*x418;
    const double x621 = x3*x620;
    const double x622 = x3*x520;
    const double x623 = x318*x431;
    const double x624 = alphay*x191;
    const double x625 = x194*x624;
    const double x626 = w_x*x339;
    const double x627 = x512*x527;
    const double x628 = w_y*x189;
    const double x629 = x318*x354;
    const double x630 = x328*x52;
    const double x631 = x467*x630;
    const double x632 = w_y*x193;
    const double x633 = x194*x632;
    const double x634 = 2*x542;
    const double x635 = x433*x52;
    const double x636 = x351*x92;
    const double x637 = 2*x431;
    const double x638 = x184*x197;
    const double x639 = x493*x638;
    const double x640 = x579*x81;
    const double x641 = 3*w_y;
    const double x642 = x334*x612;
    const double x643 = 4*x184;
    const double x644 = x338*x643;
    const double x645 = x183*x644;
    const double x646 = (1.0/2.0)*x70;
    const double x647 = x646*x73;
    const double x648 = w_x*x338;
    const double x649 = x647*x648;
    const double x650 = pow(x5, -2);
    const double x651 = x650*(-x339 - x457);
    const double x652 = x3*x72;
    const double x653 = x651*x652;
    const double x654 = alphay*x71;
    const double x655 = v*x74;
    const double x656 = x654*x655;
    const double x657 = v_dot*x76;
    const double x658 = x657*x74;
    const double x659 = u_dot*x78;
    const double x660 = x6*x614;
    const double x661 = x246*x3;
    const double x662 = v_ddot*x11;
    const double x663 = v_dot*x11;
    const double x664 = x283*x431;
    const double x665 = x265*x664;
    const double x666 = x102*x664;
    const double x667 = x301*x474;
    const double x668 = x301*x475;
    const double x669 = x116*x431;
    const double x670 = x284*x669;
    const double x671 = x286*x669;
    const double x672 = x50*(-48*x260 - 48*x404 - 48*x407 - 48*x409)/pow(x282, 2);
    const double x673 = x102*x672;
    const double x674 = x3*x673;
    const double x675 = x4*x673;
    const double x676 = x265*x672;
    const double x677 = x3*x676;
    const double x678 = x4*x676;
    const double x679 = x300*x474;
    const double x680 = x300*x475;
    const double x681 = x112*x431;
    const double x682 = 6*x681;
    const double x683 = x284*x682;
    const double x684 = x286*x682;
    const double x685 = x120*x672;
    const double x686 = x291*x685;
    const double x687 = 3*x290*x432;
    const double x688 = x293*x433;
    const double x689 = x293*x434;
    const double x690 = 6*x685;
    const double x691 = x10*x690;
    const double x692 = x11*x690;
    const double x693 = 12*x290*x6;
    const double x694 = x338*x693;
    const double x695 = x241*x693;
    const double x696 = 2*x672;
    const double x697 = x14*x696;
    const double x698 = x16*x696;
    const double x699 = 8*x404;
    const double x700 = 8*x260;
    const double x701 = x299*x672;
    const double x702 = 8*x407;
    const double x703 = 8*x409;
    const double x704 = x101*x697 + x101*x698 + x101*x701 + x112*x697 + x112*x698 + x112*x701 + x300*x699 + x300*x700 + x300*x702 + x300*x703 + x301*x699 + x301*x700 + x301*x702 + x301*x703;
    const double x705 = x304*(-x665 + x666 - x667 - x668 - x670 - x671 - x674 - x675 + x677 + x678 + x679 + x680 + x683 + x684 - x686 - x687 - x688 - x689 + x691 + x692 + x694 + x695 + x704);
    const double x706 = (1.0/4.0)*x705;
    const double x707 = x304*(x665 - x666 + x667 + x668 + x670 + x671 + x674 + x675 - x677 - x678 - x679 - x680 - x683 - x684 + x686 + x687 + x688 + x689 - x691 - x692 - x694 - x695 + x704);
    const double x708 = (1.0/4.0)*x663;
    const double x709 = x262*x470;
    const double x710 = x581*x610;
    const double x711 = h*x242;
    const double x712 = x538*x711;
    const double x713 = x242*x51;
    const double x714 = x493*x713;
    const double x715 = x242*x505;
    const double x716 = u_ddot*x532;
    const double x717 = x368*x515;
    const double x718 = x182*x84;
    const double x719 = x635*x718;
    const double x720 = x433*x630;
    const double x721 = x431*x89;
    const double x722 = x182*x721;
    const double x723 = x544*x547;
    const double x724 = v*x723;
    const double x725 = x262*x498;
    const double x726 = 16*x38;
    const double x727 = 16*x338;
    const double x728 = 16*x241;
    const double x729 = (-x25*x728 - x260*x726 - x404*x726 - x432*x55 - x432*x56 - x432*x57 - x45*x727)/pow(x255, 2);
    const double x730 = x15*x729;
    const double x731 = x22*x581;
    const double x732 = x565*x731;
    const double x733 = v_dot*x407;
    const double x734 = x526*x731;
    const double x735 = x605*x711;
    const double x736 = x163*x241;
    const double x737 = x736*x84;
    const double x738 = x243*x731;
    const double x739 = mu*x493;
    const double x740 = x739*x90;
    const double x741 = x478*x52;
    const double x742 = x120*x542;
    const double x743 = x546*x547;
    const double x744 = x517*x81;
    const double x745 = x338*x744;
    const double x746 = x373*x623;
    const double x747 = x318*x542;
    const double x748 = x433*x81;
    const double x749 = x156*x748;
    const double x750 = 2*x520;
    const double x751 = x120*x579;
    const double x752 = u*x467;
    const double x753 = x51*x542;
    const double x754 = x752*x753;
    const double x755 = x363*x432;
    const double x756 = x351*x431;
    const double x757 = w_y*x51;
    const double x758 = x493*x757;
    const double x759 = u*x758;
    const double x760 = x194*x759;
    const double x761 = x169*x51;
    const double x762 = x368*x761;
    const double x763 = x182*x318;
    const double x764 = x433*x763;
    const double x765 = x184*x764;
    const double x766 = u_ddot*x236;
    const double x767 = x246*x470;
    const double x768 = x6*x766;
    const double x769 = x346*x6;
    const double x770 = (1.0/4.0)*x707;
    const double x771 = x334*x338;
    const double x772 = v_dot*x338;
    const double x773 = x6*x772;
    const double x774 = x651*x76;
    const double x775 = x236*x774;
    const double x776 = (-x441 - x443 - x456 - x458)/pow(x244, 2);
    const double x777 = x71*x776;
    const double x778 = x22*x776;
    const double x779 = x708*x778;
    const double x780 = v_dot*x434;
    const double x781 = x243*x729;
    const double x782 = x70*x781;
    const double x783 = x544*x84;
    const double x784 = x713*x783;
    const double x785 = x242*x755;
    const double x786 = x492*x87;
    const double x787 = x70*x729;
    const double x788 = x271*x787;
    const double x789 = u_dot*x433;
    const double x790 = x353*x84;
    const double x791 = x790*x90;
    const double x792 = h*x357;
    const double x793 = x432*x81;
    const double x794 = x163*x457;
    const double x795 = x373*x640;
    const double x796 = x579*x718;
    const double x797 = x236*x757*x796;
    const double x798 = x120*x433;
    const double x799 = x243*x777;
    const double x800 = (1.0/4.0)*x778;
    const double x801 = x769*x800;
    const double x802 = x346*x432;
    const double x803 = x63*x786;
    const double x804 = x191*x241;
    const double x805 = x468*x6;
    const double x806 = x544*x87;
    const double x807 = x242*x432;
    const double x808 = x243*x357*x580;
    const double x809 = x373*x751;
    const double x810 = alphaz*x646;
    const double x811 = x74*x810;
    const double x812 = u_dot*x234;
    const double x813 = x651*x80;
    const double x814 = u*x813;
    const double x815 = x238*x654;
    const double x816 = x238*x646;
    const double x817 = u*x238;
    const double x818 = v*x238;
    const double x819 = x76*x818;
    const double x820 = x238*x402;
    const double x821 = x431*x80;
    const double x822 = x431*x84;
    const double x823 = x650*x822;
    const double x824 = x236*x76;
    const double x825 = (-3*x241 - 3*x338)/pow(x5, 5.0/2.0);
    const double x826 = x81*x825;
    const double x827 = (1.0/2.0)*x432;
    const double x828 = x101*x827;
    const double x829 = x112*x827;
    const double x830 = x377*(x827 + x828 - x829);
    const double x831 = x211*x379;
    const double x832 = alphax*x831;
    const double x833 = x224*x426;
    const double x834 = x224*x428;
    const double x835 = x376*x825;
    const double x836 = -x241 - x338;
    const double x837 = x377*x836;
    const double x838 = x379*x837;
    const double x839 = w_x*x838;
    const double x840 = x184*x82;
    const double x841 = 2*x468;
    const double x842 = x840*x841;
    const double x843 = x224*x503;
    const double x844 = u_ddot*x138;
    const double x845 = x211*x382;
    const double x846 = x211*x214;
    const double x847 = x224*x527;
    const double x848 = x224*x533;
    const double x849 = x182*x837;
    const double x850 = x849*x92;
    const double x851 = x211*x492;
    const double x852 = x851*x92;
    const double x853 = x207*x3;
    const double x854 = x82*(x827 - x828 + x829);
    const double x855 = x166*x211;
    const double x856 = x213*x855;
    const double x857 = x431*x73;
    const double x858 = x857*x87;
    const double x859 = alphaz*x138;
    const double x860 = v*x859;
    const double x861 = x139*x73;
    const double x862 = 2*x681;
    const double x863 = w_z*x141;
    const double x864 = x101*x637;
    const double x865 = x139*x84;
    const double x866 = x637*x73;
    const double x867 = x237*x836;
    const double x868 = x81*x867;
    const double x869 = 2*x139;
    const double x870 = u_dot*x593;
    const double x871 = x382*x837;
    const double x872 = 4*x469;
    const double x873 = x214*x837;
    const double x874 = (1.0/2.0)*x468;
    const double x875 = x242*x849;
    const double x876 = x242*x851;
    const double x877 = v*x205;
    const double x878 = x218*x867;
    const double x879 = x390*x837;
    const double x880 = x468*x79;
    const double x881 = x120*x13;
    const double x882 = x857*x881;
    const double x883 = x101*x221;
    const double x884 = x188*x225;
    const double x885 = x166*x387;
    const double x886 = 2*x161;
    const double x887 = x101*x857;
    const double x888 = x623*x73;
    const double x889 = 2*x92;
    const double x890 = x478*x81;
    const double x891 = x211*x890;
    const double x892 = x120*x867;
    const double x893 = v*x143;
    const double x894 = x73*x893;
    const double x895 = 2*x868;
    const double x896 = w_z*x593;
    const double x897 = v*x896;
    const double x898 = x101*x197;
    const double x899 = w_x*x849;
    const double x900 = x681*x73;
    const double x901 = u*x395;
    const double x902 = x201*x857;
    const double x903 = v*x73;
    const double x904 = x431*x881;
    const double x905 = x547*x904;
    const double x906 = x431*x521;
    const double x907 = x84*x857;
    const double x908 = u*x478;
    const double x909 = x211*x457*x908;
    const double x910 = x166*x188;
    const double x911 = x212*x492;
    const double x912 = x120*x211*x478;
    const double x913 = 2*x892;
    const double x914 = x431*x74;
    const double x915 = 4*w_y;
    const double x916 = x318*x416 + x419*x84 + x422*x84 + x501*x84 + x506*x84 + x535*x84 + x715*x84;
    const double x917 = x655*x76;
    const double x918 = x521*x72;
    const double x919 = x74*x80;
    const double x920 = -x247*x346 - x252*x536 + x263*x526 + x263*x536 - x269*x526 + x269*x565 - x271*x471 + x271*x551 + x273*x526 + x273*x536 - x273*x565 - x275*x526 - x275*x536 + x275*x565 - x276*x334 - x281*x526 + x281*x565 + x306*x663 + x306*x769 + x310*x663 + x310*x769 + x334*x464 - x334*x549 + x334*x661 + x72 + x917 - x918 + ((x85) ? (
        x238*x652 - x239*x76 - x84*x919 + x919
    )
    : (
        -x917 + x918
    ));
    const double x921 = x193*x3;
    const double x922 = u_dot*x201;
    const double x923 = x3*x425;
    const double x924 = x354*x92;
    const double x925 = x183*x638;
    const double x926 = x242*x354;
    const double x927 = x363*x752;
    const double x928 = x192*x194;
    const double x929 = x379*x840;
    const double x930 = x224*x92;
    const double x931 = x224*x242;
    const double x932 = x607 + x737;
    const double x933 = w_x*x331 + w_x*x403 + w_x*x923 + w_y*x319 + w_y*x321 + w_y*x352 - w_y*x374 + w_z*x348 - w_z*x588 + x101*x161 + x112*x161 + x112*x207 + x112*x736 - x156*x175 + x164*x241 - x193 - x207*x318 + x207 + x242*x345 + x242*x347 + x242*x351 - x242*x370 - x242*x372 + x242*x375 + x345*x92 + x347*x92 - x359*x540 - x360*x65 + x368*x65 - x369*x537 - x371*x736 + x514*x52 - x520*x616 + x636 - x794 + x804 + x84*x924 + x84*x926 - x84*x927 - x853 - x877 - x886 + x921 + x922 - x924 + x925 - x926 + x927 + x928 + x932 + ((x85) ? (
        -u_dot*x846 + w_x*x378 - x112*x929 - x112*x930 - x112*x931 + x139*x206 - x139*x210 - x161 - x188*x885 - x206*x893 + x209*x794 + x209*x886 + x210*x893 + x219*x853 + x219*x877 - 1.0/2.0*x227*x242 + x382*x840 + x390*x855 - x393*x930 - x638*x883 - x736 - x794*x83 - x83*x886 - x929 + x930 + x931 + x932
    )
    : (
        x139 - x804 + x84*x893 + x853 - x865 + x877 - x893 - x921 - x925 - x928
    ));
    const double x934 = x4*x481;
    const double x935 = x320*x418;
    const double x936 = 4*x935;
    const double x937 = x320*x415;
    const double x938 = 2*x937;
    const double x939 = v_ddot*x153;
    const double x940 = x180*x939;
    const double x941 = u_ddot*x16;
    const double x942 = x348*x415;
    const double x943 = v_ddot*x127;
    const double x944 = x4*x528;
    const double x945 = v_ddot*x158;
    const double x946 = v_ddot*x4;
    const double x947 = v_dot*x323;
    const double x948 = x16*x505;
    const double x949 = x0*x418;
    const double x950 = v*x401;
    const double x951 = u*x950;
    const double x952 = w_x*x142;
    const double x953 = u_dot*x208;
    const double x954 = x354*x614;
    const double x955 = x336*x465;
    const double x956 = x336*x609;
    const double x957 = x493*x67;
    const double x958 = x4*x594;
    const double x959 = u_ddot*x373;
    const double x960 = 2*x959;
    const double x961 = x316*x339;
    const double x962 = x204*x457;
    const double x963 = x180*x67;
    const double x964 = x478*x963;
    const double x965 = alphay*x342;
    const double x966 = x512*x614;
    const double x967 = 2*x348*x418;
    const double x968 = 2*x935;
    const double x969 = x51*x662;
    const double x970 = x518*x586;
    const double x971 = x654*x73;
    const double x972 = x4*x971;
    const double x973 = x236*x414;
    const double x974 = x187*x338;
    const double x975 = x316*x338;
    const double x976 = x338*x67;
    const double x977 = v_ddot*x273;
    const double x978 = x35*x478;
    const double x979 = u_ddot*x258;
    const double x980 = x270*x581;
    const double x981 = x22*x980;
    const double x982 = x185*x478;
    const double x983 = x249*x538;
    const double x984 = v_dot*x130;
    const double x985 = x326*x493;
    const double x986 = x326*x478;
    const double x987 = x544*x67;
    const double x988 = x101*x4;
    const double x989 = x262*x949;
    const double x990 = x516*x7;
    const double x991 = x354*x766;
    const double x992 = x346*x490;
    const double x993 = x353*x771;
    const double x994 = v*x412;
    const double x995 = u*x994;
    const double x996 = v_ddot*x523;
    const double x997 = x465*x532;
    const double x998 = x249*x334*x604;
    const double x999 = 3*x269;
    const double x1000 = x4*x426;
    const double x1001 = 3*x1000;
    const double x1002 = u_dot*x260;
    const double x1003 = 3*x1002;
    const double x1004 = v*x339;
    const double x1005 = x180*x516;
    const double x1006 = x180*x822;
    const double x1007 = alphax*x191*x194;
    const double x1008 = x512*x766;
    const double x1009 = x328*x67;
    const double x1010 = x1009*x364;
    const double x1011 = w_x*x187;
    const double x1012 = x1011*x188;
    const double x1013 = x186*x226;
    const double x1014 = x434*x67;
    const double x1015 = x180*x758;
    const double x1016 = 3*w_x;
    const double x1017 = w_y*x241;
    const double x1018 = 4*x1017*x190;
    const double x1019 = x1017*x647;
    const double x1020 = x240*x651;
    const double x1021 = x402*x74;
    const double x1022 = v*x1021;
    const double x1023 = v_dot*x75;
    const double x1024 = x426*x6;
    const double x1025 = u_dot*x72;
    const double x1026 = x1025*x903;
    const double x1027 = h*x346;
    const double x1028 = x1027*x538;
    const double x1029 = x346*x51;
    const double x1030 = x1029*x493;
    const double x1031 = x346*x505;
    const double x1032 = x166*x338;
    const double x1033 = u*x101;
    const double x1034 = x262*x523;
    const double x1035 = u*x723;
    const double x1036 = u_dot*x409;
    const double x1037 = x274*x729;
    const double x1038 = x17*x22;
    const double x1039 = x1014*x718;
    const double x1040 = x182*x326;
    const double x1041 = x1027*x605;
    const double x1042 = x261*x731;
    const double x1043 = x1009*x434;
    const double x1044 = x241*x744;
    const double x1045 = x434*x81;
    const double x1046 = x1045*x316;
    const double x1047 = 2*x361;
    const double x1048 = x365*x753;
    const double x1049 = x185*x493;
    const double x1050 = x1049*x188;
    const double x1051 = x434*x763;
    const double x1052 = x1051*x757;
    const double x1053 = 3*x262;
    const double x1054 = x527*x6;
    const double x1055 = x651*x72;
    const double x1056 = x1055*x236;
    const double x1057 = x533*x6;
    const double x1058 = x312*x778;
    const double x1059 = x1029*x783;
    const double x1060 = x346*x755;
    const double x1061 = u*x796;
    const double x1062 = x1061*x185;
    const double x1063 = x120*x434;
    const double x1064 = x307*x800;
    const double x1065 = x810*x903;
    const double x1066 = v_dot*x234;
    const double x1067 = v*x813;
    const double x1068 = x72*x817;
    const double x1069 = x236*x72;
    const double x1070 = alphay*x831;
    const double x1071 = x224*x614;
    const double x1072 = x224*x946;
    const double x1073 = w_y*x838;
    const double x1074 = x212*x841;
    const double x1075 = x224*x620;
    const double x1076 = v_ddot*x138;
    const double x1077 = x224*x766;
    const double x1078 = x221*x976;
    const double x1079 = x1040*x867;
    const double x1080 = x492*x82;
    const double x1081 = x1080*x326;
    const double x1082 = 4*x757;
    const double x1083 = x180*x227;
    const double x1084 = u*x859;
    const double x1085 = w_z*x139;
    const double x1086 = x411*x67;
    const double x1087 = u*x846;
    const double x1088 = x220*x867;
    const double x1089 = x346*x849;
    const double x1090 = x346*x851;
    const double x1091 = 4*x241;
    const double x1092 = x389*x73;
    const double x1093 = x388*x849;
    const double x1094 = u*x896;
    const double x1095 = x4*x902;
    const double x1096 = (1.0/2.0)*x1032;
    const double x1097 = 4*x914;
    const double x1098 = x1031*x84 + x318*x935 + x84*x937 + x84*x945 + x84*x959 + x84*x975 + x84*x986;
    const double x1099 = v*x844;
    const double x1100 = u*x1076;
    const double x1101 = x143*x339;
    const double x1102 = x143*x457;
    const double x1103 = x3*x859;
    const double x1104 = x4*x859;
    const double x1105 = v_dot*x180;
    const double x1106 = alphaz*x153;
    const double x1107 = v*x156;
    const double x1108 = v*x870;
    const double x1109 = u*x1086;
    const double x1110 = x3*x896;
    const double x1111 = x4*x896;
    const double x1112 = x188*x528;
    const double x1113 = x194*x414;
    const double x1114 = 4*x196;
    const double x1115 = x1114*x338;
    const double x1116 = x1114*x241;
    const double x1117 = alphaz*x197;
    const double x1118 = x1117*x190;
    const double x1119 = alphaz*x342;
    const double x1120 = v_ddot*x0;
    const double x1121 = u*x614;
    const double x1122 = x1121*x23;
    const double x1123 = x1122*x38;
    const double x1124 = x13*x28;
    const double x1125 = u_ddot*x25;
    const double x1126 = v_dot*x26;
    const double x1127 = u_ddot*x35;
    const double x1128 = v*x426;
    const double x1129 = x1128*x38;
    const double x1130 = v_ddot*x45;
    const double x1131 = x28*x33;
    const double x1132 = x185*x544;
    const double x1133 = x1061*x67;
    const double x1134 = u*x757;
    const double x1135 = v*x52*x796;
    const double x1136 = x112*x544;
    const double x1137 = u_dot*x441;
    const double x1138 = x547*x764;
    const double x1139 = x1051*x547;
    const double x1140 = x2*x34*x827;
    const double x1141 = w_z*x505;
    const double x1142 = v_dot*x443;
    const double x1143 = x61*x761;
    const double x1144 = x226*x61;
    const double x1145 = x22*(-h*x405 - h*x406 - h*x408 - h*x410)/pow(x20, 2);
    const double x1146 = v*x26;
    const double x1147 = x1145*x33;
    const double x1148 = x36*x8;
    const double x1149 = (9.0/2.0)*x36*x38;
    const double x1150 = x1145*x41;
    const double x1151 = u*x46;
    const double x1152 = (1.0/2.0)*x1145;
    const double x1153 = x1152*x33;
    const double x1154 = (9.0/2.0)*x1151;
    const double x1155 = x35*x53;
    const double x1156 = 16*x2;
    const double x1157 = (-x1156*x260 - x1156*x404 - x727*x88 - x728*x86)/pow(x58, 2);
    const double x1158 = x1157*x60;
    const double x1159 = x1*x63;
    const double x1160 = v*x65;
    const double x1161 = u*x68;
    const double x1162 = x23*x435;
    const double x1163 = x23*x559;
    const double x1164 = x21*x435*x59;
    const double x1165 = x23*x92;
    const double x1166 = v*x428;
    const double x1167 = x188*x957;
    const double x1168 = x194*x494;
    const double x1169 = x197*x548;
    const double x1170 = x180*x548;
    const double x1171 = x1*x1152;
    const double x1172 = x1152*x93;
    const double x1173 = x1127*x27;
    const double x1174 = x105*x23;
    const double x1175 = x47*x946;
    const double x1176 = v_dot*x248;
    const double x1177 = x1152*x36;
    const double x1178 = v_dot*x133;
    const double x1179 = x1152*x1178;
    const double x1180 = x117*x23;
    const double x1181 = x27*x943;
    const double x1182 = x127*x23;
    const double x1183 = x1157*x70;
    const double x1184 = x1183*x130;
    const double x1185 = v*x971;
    const double x1186 = x1025*x73;
    const double x1187 = v_dot*x77;
    const double x1188 = u*x1055;
    const double x1189 = v*x774;
    const double x1190 = x1*x27;
    const double x1191 = x27*x93;
    const double x1192 = x439*x95;
    const double x1193 = x437*x99;
    const double x1194 = x103*x437;
    const double x1195 = x2*x277;
    const double x1196 = x562*x88;
    const double x1197 = x38*x571;
    const double x1198 = x100*x560 + x102*x1194 + x102*x1197 - 3*x1192 + 2*x1193 + x1195*x555 + x1195*x557 + 3*x1196 + x474*x570 + x475*x570 + x562*x96;
    const double x1199 = 12*x570;
    const double x1200 = x114*x567;
    const double x1201 = x112*x573;
    const double x1202 = x109*x562 + x110*x560 + x111*x1201 - x114*x569 + x115*x1201 + x116*x1194 + x116*x1197 - 6*x1192 + 7*x1193 + 6*x1196 + x1199*x241 + x1199*x338 + x1200*x241 + x1200*x338 + x552*x555 + x552*x557;
    const double x1203 = x131*x462;
    const double x1204 = w_x*x141;
    const double x1205 = x205*x73;
    const double x1206 = x208*x73;
    const double x1207 = x166*x226;
    const double x1208 = w_z*x898;
    const double x1209 = 4*w_z;
    const double x1210 = x393*x468;

    double *ang_mom_dot = new double[3];

    ang_mom_dot[0] = -alphax*x153*x616 + alphax*x331 + alphax*x403 + alphay*x319 + alphay*x321 + alphay*x352 - alphay*x374 + alphaz*x348 - alphaz*x588 + u*x166*x531 + u_ddot*x201 + v*x530 + v*x548 - v*x743 + v_dot*x275*x612 + w_y*x135 + w_y*x154 + w_y*x157 + w_y*x165 + w_y*x167 - w_y*x172 - w_y*x174 + w_y*x231 - w_y*x746 - w_y*x795 + w_y*x809 + w_z*x187 - w_z*x313 - w_z*x330 - w_z*x400 + x10*x193*x468 + x10*x717 + x101*x419 + x101*x422 + x101*x506 + x101*x535 - x101*x608 + x101*x715 + x101*x720 + x101*x785 + x112*x401 + x112*x412 + x112*x419 + x112*x422 + x112*x501 + x112*x506 + x112*x511 + x112*x535 + x112*x715 - x112*x724 - x112*x785 + x127*x539 + x127*x712 + x130*x539 + x130*x712 + x158*x626 + x160*x641 + x162*x641 + x164*x509 + x166*x723 - x170*x236*x531 - x170*x416 - x175*x421 - x175*x741 - x18*x334*x787 + x18*x583 + x184*x479 + x207*x640 - x207*x751 - x241*x359*x52 + x241*x368*x53 + x242*x461 + x242*x463 + x242*x803 - x242*x808 + x243*x585 - x243*x789*x792 + x246*x771 - x247*x614 - x247*x766 + x251*x642 - x252*x615 + x262*x334*x407 - x262*x642 + x263*x589 + x263*x615 - x263*x621 + x263*x716 + x268*x732 - x268*x734 - x269*x480 + x269*x550 - x269*x589 + x269*x621 - x269*x733 - x271*x585 + x273*x480 - x273*x550 + x273*x615 - x273*x621 - x275*x589 - x275*x615 + x275*x621 - x275*x716 - x275*x733 - x276*x465 + x276*x495 - x276*x609 + x280*x732 - x280*x734 - x281*x480 + x281*x550 - x281*x589 + x281*x621 - x281*x733 + x3*x334*x777 + x303*x779 + x303*x801 + x305*x587 + x306*x660 + x306*x662 + x306*x768 + x306*x773 + x306*x780 + x306*x802 + x308*x779 + x308*x801 + x309*x587 + x310*x660 + x310*x662 + x310*x768 + x310*x773 + x310*x780 + x310*x802 - x318*x401 - x318*x412 - x328*x433*x545 + x341*x635 + x345*x426 + x345*x428 + x345*x503 + x345*x527 + x345*x533 - x346*x767 - x346*x799 + x347*x426 + x347*x428 + x347*x503 + x347*x527 + x347*x533 + x351*x426 + x351*x428 + x351*x503 + x351*x527 + x351*x533 - x354*x503 + x368*x635 - x370*x426 - x370*x428 - x370*x503 - x370*x527 - x370*x533 - x371*x419 - x371*x501 - x371*x535 - x371*x715 - x372*x426 - x372*x428 - x372*x503 - x372*x527 + x375*x426 + x375*x527 + x401 + x402 + x405*x520 + x407*x613 + x412 - x413 - x414 - x417 - x420 - x423 + x424*x425 + x425*x626 + x427 + x429 - x431*x607 - x431*x737 + x457*x622 + x461*x92 + x463*x92 + x464*x465 - x464*x495 - x465*x549 + x465*x661 + x470*x471 - x471*x498 - x483 - x484 + x486*x84 - x486 + x487*x84 - x487 - x489 + x491*x84 - x491 - x494 + x495*x549 + x498*x551 - x499 - x500 - x502 - x504 - x505*x617 - 2*x506 + x508 + x510 - x513*x84 + x513 + x514*x515 + x516*x629 + x516*x762 - x519*x84 + x519 + x52*x546 - x520*x745 - x522 + x523*x524 + x523*x531 + x524*x532 + x525*x526 + x525*x536 - x525*x565 - x526*x564 - x526*x575 + x526*x709 + x526*x710 - x526*x725 + x526*x782 - x526*x788 + x529 + x531*x532 - x532*x597 + x534 - x536*x599 + x536*x709 + x536*x710 - x536*x725 - x536*x738 + x536*x782 - x536*x788 - x537*x603 - x537*x740 + x540*x543 + x540*x722 - x540*x742 + x541 + x543*x713 + x544*x545 - x549*x609 + x564*x565 + x565*x575 - x565*x710 + x565*x725 + x565*x788 + x576*x582 - x576*x730 - x582*x611 + x583*x584 - x583*x616 + x591*x84 - x591 + x592*x84 - x592 - x596 + x598*x84 - x598 - x603*x711 + x606*x84 - x606 - x608*x84 + x609*x661 + x611*x730 - x612*x613 + x617*x755 + x617*x756 - x618*x619 - x622*x623 + x625 - x627*x84 + x627 + x628 - x631*x84 + x631 + x633 + x634*x65 + x636*x637 + x639 + x645 - x649 - x65*x747 + x65*x786 - x653 + x656 + x658 + x659 + x663*x706 + x706*x769 + x707*x708 - x711*x740 + x713*x722 - x713*x742 + x713*x806 + x714 + x719 - x720*x84 + x724*x84 + x735*x84 - x735 - x748*x750 - 2*x749 + x750*x798 - x754*x84 + x754 + x760 + x765 + x769*x770 + x775 + x784 - x785*x84 - x789*x791 - x791*x807 - x793*x794 + x797 + x804*x805 + x916 + ((x85) ? (
        -u_dot*x819 - x236*x815 + x3*x820 + x648*x816 + x652*x823 + x652*x826 - x657*x817 - x811*x84 + x811 - x812*x84 + x812 - x814*x84 + x814 + x817*x821 - x823*x824 - x824*x826
    )
    : (
        x522 + x649 + x653 - x656 - x658 - x659 - x775
    )) + ((x85) ? (
        alphax*x378 + alphax*x845 - alphay*x884 + alphaz*x856 - u_ddot*x846 - u_dot*x391*x874 - u_dot*x873 - v_ddot*x901 + w_x*x830 + w_x*x871 - w_y*x223 - w_y*x849*x910 - x112*x832 - x112*x833 - x112*x834 - x112*x839 - x112*x842 - x112*x843 - x112*x847 - x112*x848 - x112*x850 - x112*x852 - x112*x875 - x112*x876 + x120*x909 - x161*x888 - x164*x457*x857 + x166*x879 + x184*x835 - x200*x643*x906 + x201*x241*x914 - x206*x501 + x206*x844 - x206*x860 - x206*x863 + x206*x870 - x206*x897 + x209*x417 + x209*x420 + x209*x423 + x209*x502 + x209*x504 - x210*x844 + x210*x860 + x210*x863 - x210*x870 + x210*x897 + x219*x413 + x219*x483 + x219*x484 + x219*x489 + x219*x499 + x219*x500 + x219*x596 - x226*x885 - x227*x507 - x227*x516 + x383*x390 - x393*x833 - x393*x834 - x393*x848 - x393*x850 - x393*x852 - x393*x875 - x393*x876 - x417*x83 - x419 - x420*x83 - x422 - x423*x83 - x492*x840*x898 - x501 - x504*x83 - x506 - x511 - x52*x882 - x535 - x588*x914*x915 - x643*x858 + x643*x882 - x644*x883 - x715 - x736*x793 - x749 - x794*x868 + x794*x892 + x794*x900 - x794*x907 - x81*x909 - x832 + x833 + x834 - x839 + x840*x872 - x842 + x843 + x847 + x848 + x850 + x852 + x853*x854 + x853*x878 + x854*x877 + x855*x880 - x861*x862 + x861*x864 + x862*x894 - x864*x894 + x865*x866 + x868*x869 - x868*x886 - x869*x892 + x875 + x876 + x877*x878 - x886*x887 + x886*x892 + x886*x900 - x888*x893 - x889*x891 + x889*x912 - x893*x895 + x893*x913 - x898*x899 + x902*x92 + x903*x905 - x910*x911 + x916
    )
    : (
        w_z*x399 + x139*x793 + x413 - x427 - x429 + x483 + x484 + x489 + x499 + x500 - x508 - x510 - x529 - x534 - x541 + x596 - x625 - x628 - x633 - x639 - x645 - x714 - x719 - x760 - x765 - x784 - x793*x893 - x797 - x84*x844 + x84*x860 - x84*x870 + x84*x897 + x844 - x860 - x863 + x870 - x897
    ));
    ang_mom_dot[1] = alphax*x319 + alphax*x321 + alphax*x352 - alphax*x374 - alphay*x153*x362 + alphay*x314 + alphay*x331 + alphay*x332 - alphaz*x191 + alphaz*x202 - u*x170*x982 + u*x184*x978 - u*x548 + u*x743 - u_ddot*x18*x263 + u_ddot*x247*x4 + u_ddot*x271*x276 + u_ddot*x352 - u_dot*x1038*x781 + v_ddot*x201 + v_ddot*x324 - v_ddot*x358*x90 - v_dot*x358*x602 - w_x*x135 + w_x*x154 + w_x*x157 + w_x*x165 + w_x*x167 - w_x*x172 - w_x*x174 - w_x*x231 - w_x*x746 - w_x*x795 + w_x*x809 + w_y*x158*x339 + w_y*x425*x457 + w_y*x948 - w_z*x193 + w_z*x920 + w_z*x922 + w_z*x933 + x0*x982 - x1000*x999 + x1001*x273 - x1001*x275 - x1001*x281 + x1002*x999 - x1003*x263 - x1003*x273 + x1003*x281 - x1004*x505 + x1004*x755 + x1004*x756 + x1005*x252 - x1005*x263 - x1005*x273 + x1005*x275 - x1006*x361 + x1007 - x1008*x84 + x1008 + x101*x937 + x101*x945 + x101*x968 + x101*x975 + x101*x986 - x1010*x84 + x1010 + x1012 + x1013 + x1014*x341 + x1014*x368 + x1015 + x1016*x160 + x1016*x162 + x1018 - x1019 - x1020 + x1022 + x1023 + x1024*x306 + x1024*x310 + x1026 - x1027*x603 - x1027*x740 + x1028*x127 + x1028*x130 + x1029*x543 + x1029*x722 - x1029*x742 + x1029*x806 + x1030 + x1031*x112 - x1031*x371 + x1032*x505 + x1032*x755 - x1033*x723 + x1033*x982 + x1034*x465 - x1034*x495 + x1035*x112 - x1035*x84 - x1036*x252 + x1036*x263 + x1036*x273 + x1037*x1038 + x1037*x18*x70 + x1039 + x1040*x721 + x1041*x84 - x1041 - x1042*x268 - x1042*x280 - x1043*x112 - x1043*x84 - x1044*x361 - x1045*x1047 - 2*x1046 + x1047*x1063 - x1048*x84 + x1048 + x1050 + x1052 + x1053*x260*x274 + x1053*x4*x485 + x1054*x306 + x1054*x310 + x1056 + x1057*x306 + x1057*x310 + x1058*x303 + x1058*x308 + x1059 - x1060*x112 - x1060*x84 + x1062 + x1064*x303 + x1064*x308 + x1098 + x11*x187*x468 + x112*x481 + x112*x594 + x112*x937 + x112*x945 + x112*x959 + x112*x968 + x112*x975 + x112*x986 + x112*x987 + x129*x338*x344 + x129*x983 - x170*x935 - x177*x939 - x18*x274*x731 - x18*x524 + x18*x531 + x18*x597 + x204*x640 - x204*x751 + x220*x328*x434 - x220*x544 - x236*x359*x515 + x236*x6*x717 - x242*x767 - x242*x799 - x243*x979 + x243*x989 - x246*x485 - x246*x590 + x247*x533 + x248*x767 + x248*x799 - x249*x356*x739 - x251*x997 - x252*x949 + x253*x525 - x253*x599 + x253*x709 + x253*x710 - x253*x725 - x253*x738 + x253*x782 - x253*x788 - x254*x525 + x254*x599 - x254*x709 + x254*x725 + x254*x738 - x254*x782 - x259*x470 + x259*x498 + x261*x525 - x261*x564 - x261*x575 + x261*x709 + x261*x710 - x261*x725 + x261*x782 - x261*x788 - x262*x274*x409 - x262*x495*x532 + x262*x997 + x268*x981 + x269*x941 - x269*x996 + x271*x979 - x271*x989 - x272*x498 + x272*x563 + x272*x574 - x273*x941 + x273*x949 + x280*x981 + x281*x941 - x281*x996 + x3*x361*x457 + x305*x990 + x306*x619 + x306*x789 + x306*x807 + x307*x706 + x307*x770 + x309*x990 + x310*x619 + x310*x789 + x310*x807 + x312*x705 + x312*x707 - x318*x481 - x318*x594 + x323*x943 + x326*x543 - x326*x742 + x329*x431*x963 + x329*x976 - x334*x434*x792 - x338*x356*x792 - x340*x793 + x341*x969 + x345*x614 + x345*x620 + x345*x766 + x346*x461 + x346*x463 + x346*x803 - x346*x808 + x347*x614 + x347*x620 + x347*x766 + x347*x772 + x349*x805 + x351*x614 + x351*x620 + x351*x766 + x351*x946 - x354*x620 - x356*x434*x790 - x359*x976 + x361*x406 - x362*x531 + x368*x6*x976 + x368*x969 - x370*x614 - x370*x620 - x370*x766 - x371*x937 - x371*x959 - x371*x975 - x372*x614 - x372*x620 - x372*x946 + x375*x614 - x396*x431 - x397*x431 + x460*x947 + x462*x947 + x481 + x523*x977 - x528 + x531*x988 + x532*x977 + x546*x67 + x586*x629 + x586*x762 + x594 - x610*x980 - x618*x662 + x634*x68 + x654 - x68*x747 + x68*x786 - x68*x890 - x791*x802 + x84*x954 + x84*x955 + x84*x956 - x84*x966 - x84*x970 - x84*x987 + x84*x991 + x84*x992 + x84*x993 + x84*x998 - x934 - x936 - x938 - x940 + x942 + x944 - x951 - x952 - x953 - x954 - x955 - x956 - x957 - x958 - x960 - x961 - x962 - x964 + x965 + x966 + x967 + x970 - x972 + x973 + x974 + x983*x984 + x985 - x991 - x992 - x993 - x995 - x998 + ((x85) ? (
        -v_dot*x1068 + x1017*x816 - x1025*x818 + x1065*x84 - x1065 + x1066*x84 - x1066 + x1067*x84 - x1067 - x1069*x823 - x1069*x826 - x236*x820 + x240*x823 + x240*x826 + x4*x815 - x818*x821
    )
    : (
        x1019 + x1020 - x1022 - x1023 - x1026 - x1056 + x972
    )) + ((x85) ? (
        -alphax*x884 - alphay*x1083 + alphay*x378 + alphay*x845 - alphaz*x1087 - u*x391*x880 - u_ddot*x901 - u_dot*x392 - v_ddot*x846 - w_x*x223 - w_y*x1093 + w_y*x830 + w_y*x871 - x1004*x891 + x1004*x912 - x101*x1091*x387 - x1031 - x1033*x879 - x1046 - x1070*x112 - x1070 - x1071*x112 + x1071 - x1072*x112 - x1072*x393 + x1072 - x1073*x112 - x1073 - x1074*x112 - x1074 - x1075*x112 + x1075 + x1076*x206 - x1076*x210 - x1077*x112 + x1077 - x1078*x112 + x1078 - x1079*x112 - x1079*x393 + x1079 - x1080*x188*x394 - x1081*x112 - x1081*x393 + x1081 - x1082*x858 + x1082*x882 + x1084*x206 - x1084*x210 + x1085*x206 - x1085*x210 + x1086*x206 - x1086*x210 - x1088*x213 - x1089*x112 + x1089 - x1090*x112 + x1090 - x1092*x862 + x1092*x864 + x1094*x206 - x1094*x210 - x1095*x915 - x1096*x849 - x1096*x851 - x1097*x203 + x1098 - x185*x226*x883 + x200*x326*x857 - x206*x959 + x209*x936 + x209*x938 + x209*x940 + x209*x960 + x209*x961 + x209*x964 + x212*x872 + x219*x934 + x219*x951 + x219*x952 + x219*x953 + x219*x958 + x219*x962 + x219*x995 - 1.0/2.0*x220*x221*x338 - x225*x418 + x315*x854 + x315*x878 - x317*x868 - x317*x887 + x317*x892 + x317*x900 - x317*x907 + x333*x854 + x333*x878 - x340*x868 - x340*x887 + x340*x892 + x340*x900 - x340*x907 + x346*x902 - x350*x793 - x383*x874 + x386*x868 + x386*x887 - x386*x892 - x386*x900 - x388*x911 + x389*x888 + x389*x895 - x389*x913 - x395*x415 + x399*x866 - x67*x882 - x74*x905 + x757*x835 - x83*x936 - x83*x938 - x83*x940 - x83*x961 - x83*x964 - x899*x910 - x937 - x945 - x959 - x968 - x975 - x986
    )
    : (
        -w_z*x865 - x1007 - x1012 - x1013 - x1015 - x1018 - x1030 - x1039 - x1050 - x1052 - x1059 - x1062 - x1076*x84 + x1076 - x1084*x84 + x1084 + x1085 - x1086*x84 + x1086 - x1094*x84 + x1094 + x141*x793 + x389*x793 + x934 - x942 - x944 + x951 + x952 + x953 + x958 + x962 - x965 - x967 - x973 - x974 - x985 + x995
    ));
    ang_mom_dot[2] = alphax*x348 - alphax*x588 + alphay*x202 + alphaz*x314 + alphaz*x332 + alphaz*x403 + alphaz*x923 - u*x132*x946 + u*x945 + u*x986 + u_dot*v_dot*x2*x39*x42 + u_dot*x1105*x132 - u_dot*x34*x46 + v*x420 - v*x422 - v*x506 - v*x511 - v_ddot*x163*x170 + v_dot*x886 + w_x*x313 - w_x*x330 + w_x*x400 - w_y*x920 + w_y*x922 - w_y*x933 + w_z*x948 - x0*x1181 + x0*x939 + x1*x1140 - x1*x1184 - x1*x1203 - x1*x2*x24*x827 - x1*x432*x51*x61 + x1*x505 - x10*x1144*x67 + x10*x515*x66 - x1006*x155 + x101*x1132 + x1011 - x1014*x69 - x1033*x544*x757 - x1044*x155 - x1045*x176 + x1049 - x105*x1173 + x105*x1175 - x105*x1177 + x105*x1179 - x106*x1198 - x106*x1202 + x106*x460 + x1063*x176 + x108*x1198 + x108*x1202 - x108*x460 + x1099 - x1100 - x1101 - x1102 - x1103 - x1104 - x1105*x156 + x1106*x19 - x1106*x362 - x1106*x616 + x1107*x623 + x1107*x640 - x1107*x751 + x1108 - x1109 - x1110 - x1111 + x1112 - x1113 + x1115 + x1116 + x1118 + x1119 - x1120*x132 - x1120*x24*x8 + x1120*x37 - x1120*x64 + x1120*x91 + x1121*x1143 + x1121*x1182 - x1122*x90 + x1123*x13 - x1123*x33 + x1124*x1125 + x1124*x1137 - x1125*x1131 - x1126*x24 + x1126*x34 + x1127*x128 + x1127*x132 - x1127*x37 + x1127*x44 - x1128*x1143 + x1128*x1174 + x1128*x1180 - x1128*x1182 + x1129*x34 - 9*x1129*x42 + x1130*x48 - x1130*x49 - x1131*x1137 + x1132*x84 + x1133 + x1134*x1136 + x1134*x546 - x1134*x783 - x1135 - x1136*x185 + x1138 + x1139 - x1140*x36 + x1141*x19 - x1141*x362 + x1141*x584 - x1141*x616 + x1141*x988 + x1142*x48 - x1142*x49 + x1144*x68 + x1145*x1146*x213 - x1145*x13*x9 - x1146*x1153 - x1147*x1148 + x1147*x9 - x1148*x1162 + x1149*x1150 + x1149*x1163 - x1150*x1154 + x1151*x1153 - x1152*x129*x133 - x1154*x1163 + x1155*x1158 + x1155*x1164 - x1158*x1159 + x1158*x1160 - x1158*x1161 - x1159*x1164 + x1160*x1164 - x1161*x1164 + x1162*x9 - x1165*x129 + x1165*x356 + x1166*x128 + x1166*x132 - x1166*x91 + x1167 - x1168 + x1169 - x117*x1173 + x117*x1175 - x117*x1177 + x117*x1179 + x1170 - x1171*x127 + x1171*x90 + x1172*x127 - x1172*x90 - x1174*x1176 - x1176*x1180 + x1177*x127 - x1178*x1203 - x1181*x133 - x1183*x133*x984 + x1184*x36 + x1184*x93 - x1190*x460 + x1190*x602 + x1191*x460 - x1191*x602 + x1203*x36 + x1203*x93 + x129*x23*x248 + x132*x188*x614 - x132*x194*x426 - x134*x889 - x155*x3*x623 + x155*x405 + x155*x406 + x155*x408 + x155*x410 - x155*x745 - x166*x421 - x166*x741 + x168*x339 + x168*x457 - x170*x67*x908 - x171*x623 - x171*x640 + x171*x751 + x173*x421 + x173*x741 - x176*x748 + x176*x798 - x185*x546 + x220*x908 - x23*x468*x9 - x29*x435 + x29*x468 - x35*x421 + x36*x43*x437 + x432*x52*x62 + x435*x46*x47 + x479*x547 + x515*x6*x62 - x52*x978 + x530 - x624 - x632 + x635*x66 - x69*x969 - x759 + ((x85) ? (
        x1021*x84 - x1021 - x1068*x431 - x1185*x84 + x1185 + x1186*x84 - x1186 - x1187*x84 + x1187 + x1188*x84 - x1188 - x1189*x84 + x1189 + x234*x822 + x431*x819 + x80*x868 + x810*x83
    )
    : (
        x810
    )) + ((x85) ? (
        alphax*x856 - alphay*x1087 - alphaz*x1083 - u*w_y*x873 - u*x1210*x212 + 2*u_ddot*x225 - u_dot*x215 - v_ddot*x188*x227 + w_x*x166*x213*x837 + w_x*x384 - w_z*x1093 - w_z*x388*x851 - x1088*x182*x188 - x1091*x228 - x1095*x1209 - x1097*x330 - x1099*x219 + x1100*x219 + x1101*x219 + x1102*x219 + x1103*x219 + x1104*x219 - x1108*x219 + x1109*x219 + x1110*x219 + x1111*x219 - x1117*x227 - x1204*x206 + x1204*x210 - x1205*x862 + x1205*x864 + x1206*x862 - x1206*x864 + x1207*x849 + x1207*x851 - x1208*x849 - x1208*x851 - x1209*x201*x906 + x1210*x216 - x140*x854 - x140*x878 + x142*x854 + x142*x878 + x144*x854 + x144*x878 + x145*x854 + x145*x878 + x185*x882 - x188*x383*x492 + x205*x895 - x205*x913 + x206*x482 + x206*x488 + x206*x595 - x206*x950 - x206*x994 - x208*x895 + x208*x913 - x210*x482 - x210*x488 - x210*x595 + x210*x950 + x210*x994 - 4*x228*x338 - x229*x866 + x230*x866 + 4*x431*x903*x922 - x74*x757*x904
    )
    : (
        w_x*x399 - x1099 + x1100 + x1101 + x1102 + x1103 + x1104 - x1108 + x1109 + x1110 + x1111 - x1112 + x1113 - x1115 - x1116 - x1118 - x1119 - x1133 + x1135 - x1138 - x1139 - x1167 + x1168 - x1169 - x1170 - x1204 + x205*x793 - x208*x793 - x482*x84 + x482 - x488*x84 + x488 - x595*x84 + x595 + x84*x950 + x84*x994 - x950 - x994
    ));

    Vector3d H_dot;
    H_dot << ang_mom_dot[0], ang_mom_dot[1], ang_mom_dot[2];

    delete ang_mom_dot;

    return H_dot;
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
        // std::cout << "u: " << u << std::endl;
        // std::cout << "v: " << v << std::endl;
        // std::cout << "h: " << h << std::endl;
        udot = uvdot(joint_num*2);
        vdot = uvdot(joint_num*2+1);
        qdot << udot,vdot;
        qddot << uvddot(joint_num*2),uvddot(joint_num*2+1);
        // std::cout << "qddot: " << qddot << std::endl; 
        // find omega from jacobian
        MatrixXd J = J_end(u,v,h);
        // std::cout << "J: " << J << std::endl; 
        xi_dot = J*qdot;
        _v_e_rel = xi_dot.block<3,1>(0,0);
        _omega_e_rel = xi_dot.block<3,1>(3,0);
        // std::cout << "v_e_rel: " << _v_e_rel << std::endl; 
        // std::cout << "omega_rel: " << _omega_e_rel << std::endl; 
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