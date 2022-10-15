#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <set>
#include <thread>
#include <vector>
#include <Eigen/Dense>

using Vector12d = Eigen::Matrix<double, 12, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Vector3d = Eigen::Matrix<double, 3, 1>;
using Matrix12d = Eigen::Matrix<double, 12, 12>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Matrix3d = Eigen::Matrix<double, 3, 3>;
using Vector4d = Eigen::Matrix<double, 4, 1>;
using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using MatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
 using namespace matplot;
using namespace std;


Vector6d jerk_minimum_path(const Vector3d & start, const Vector3d & end, const double & T);

double duration_of_path(const double & dt,
                                                    const double & arc_length, const double &jerk_limit,
                                                    const double & acc_limit, const double & vel_limit);

Vector4d evaluate_path(const double & t, const Vector6d & coeff);

int main() {

    /**************************************************************************
     * Simulation Variables
     ***************************************************************************/

    std::vector<float> x_vals;
    std::vector<float> y_vals;
    std::vector<float> z_vals;
    std::vector<float> x_dot1_vals;
    std::vector<float> y_dot1_vals;
    std::vector<float> z_dot1_vals;
    std::vector<float> x_dot2_vals;
    std::vector<float> y_dot2_vals;
    std::vector<float> z_dot2_vals;
    std::vector<float> x_dot3_vals;
    std::vector<float> y_dot3_vals;
    std::vector<float> z_dot3_vals;

    std::vector<float> time_vals;
    std::vector<float> w_vals;
    std::vector<float> curvature_vals;

    double Cx = 0.0;
    double Cy = 0.0;
    double Cz = 0.0;

    double R = 5.0;
    double L = 5.0;


    double chi0 =  0.0;
    double chi1 = 6.28318530718;

    double gamma = 0.1;

    double arc_length = R * fabs(chi1 - chi0);

    double w = 0.00;
    double t = 0.0;
    double h = 0.2;
    double kappa = 0.0;



    double max_jerk = 0.00001;
    double max_acc = 0.125;
    double max_vel = 0.5;

    double T_z = duration_of_path(h, L, max_jerk, max_acc, max_vel);
    double T_s= duration_of_path(h, arc_length, max_jerk, max_acc, max_vel);

    Vector3d s_i(0.0, 0.0, 0.0);
    Vector3d s_f(arc_length, 0.0, 0.0);
    Vector6d C_s = jerk_minimum_path(s_i, s_f, T_s);

    Vector3d z_i(0.0, 0.0, 0.0);
    Vector3d z_f(L, 0.0, 0.0);
    Vector6d C_z = jerk_minimum_path(z_i, z_f, T_s);


    Vector4d S;
    Vector4d Z ;

    Vector3d P;
    Vector3d P_dot;
    Vector3d P_dot_dot;

    while (t <=  1.0 * T_s) {

        S = evaluate_path(t, C_s);
        Z = evaluate_path(t, C_z);

        P << Cx + R * cos(S(0)/R) * cos (gamma),
                  Cy + R * sin(S(0)/R) * cos(gamma),
                      Z(0);

        P_dot <<   S(1) * -sin(S(0)/R) * cos (gamma),
                            S(1) * cos(S(0)/R) * cos(gamma),
                            Z(1);


        P_dot_dot << S(2) * -sin(S(0)/R) * cos(gamma) +  (1/R) * S(1) * -cos(S(0)/R) * cos(gamma),
                                   S(2) * cos(S(0)/R) * cos(gamma) +  (1/R)  * S(1) * -sin(S(0)/R) * cos(gamma) ,
                                   Z(2);
        t    += 0.2;
        kappa = (P_dot_dot(1) * P_dot(0) - P_dot_dot(0) * P_dot(1)) /
                         pow((P_dot(0) * P_dot(0) + P_dot(1) * P_dot(1)), (1/3));

        std::cout << " ~ ** line segment position ** ~" << std::endl;
        std::cout << P << std::endl;
        std::cout << " ~ ** line segment velocity ** ~" << std::endl;
        std::cout << P_dot << std::endl;
        std::cout << " ~ ** line segment acceleration  ** ~" << std::endl;
        std::cout << P_dot_dot << std::endl;
        x_vals.push_back(P(0));
        y_vals.push_back(P(1));
        z_vals.push_back(P(2));
        x_dot1_vals.push_back(P_dot(0));
        y_dot1_vals.push_back(P_dot(1));
        z_dot1_vals.push_back(P_dot(2));
        x_dot2_vals.push_back(P_dot_dot(0));
        y_dot2_vals.push_back(P_dot_dot(1));
        z_dot2_vals.push_back(P_dot_dot(2));
        time_vals.push_back(t);
        w_vals.push_back(w);
         curvature_vals.push_back(kappa);
    }
//    plot3(x_vals, y_vals, z_vals, "--+");
//    title("East-North-UP");
//    xlabel("X (m)");
//    ylabel("Y (m)");
//    zlabel("Z (m)");
//    legend("P(t)", "z");


//    plot(time_vals, z_vals, "d");
//    title("Z");
//    ylabel("Z");
//    xlabel("t");

//    plot(time_vals, curvature_vals);
//    title("Curvature");
//    ylabel("κ");
//    xlabel("t");

//    plot(time_vals, w_vals);
//    title("Path Variable");
//    ylabel("ω");
//    xlabel("t");


//    hold(on);
//    plot(time_vals, x_dot1_vals, "o");
//    plot(time_vals, y_dot1_vals, "+");
//    plot(time_vals, z_dot1_vals, "d");
//    hold(off);
//    title("Path Velocity");
//    ylabel("Ẋ Ẏ Ż");
//    xlabel("t");
//    legend({"Ẋ", "Ẏ", "Ż"});
//    grid(on);

        hold(on);
        plot(time_vals, x_dot2_vals, "--o");
        plot(time_vals, y_dot2_vals, "--+");
        plot(time_vals, z_dot2_vals, "--d");
        hold(off);
        title("Path Acceleration");
        ylabel("ẍ Ÿ Z^{\. \. }");
        xlabel("t");
        legend({"ẍ", "Ÿ", "Z^{\. \. }"});

    show();

     return 0;
}
/**************************************************************************************/
Vector6d jerk_minimum_path(const Vector3d & start, const Vector3d & end, const double & T)
{
    /*
    Calculate the Jerk Minimizing path that connects the initial state
    to the final state in time T.
    INPUTS
    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]
    end   - the desired end state for vehicle. Like "start" this is a
        length three array.
    T     - The duration, in seconds, over which this maneuver should occur.
    OUTPUT
    an array of length 6, each value corresponding to a coefficent in the polynomial
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
    */

    MatrixXd A = MatrixXd(3, 3);
        A << T*T*T, T*T*T*T, T*T*T*T*T,
                            3*T*T, 4*T*T*T,5*T*T*T*T,
                            6*T, 12*T*T, 20*T*T*T;

        MatrixXd B = MatrixXd(3,1);
        B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
                            end[1]-(start[1]+start[2]*T),
                            end[2]-start[2];

       MatrixXd Ai = A.inverse();
       Vector3d C0to2(start(0), start(1), 0.5 * start(2));
       Vector3d C3to5 = Ai*B;
       Vector6d  C;
       C << C0to2, C3to5;
       return C;
}
/**************************************************************************************/
double duration_of_path(const double & dt,
                                                    const double & arc_length, const double &jerk_limit,
                                                    const double & acc_limit, const double & vel_limit)
{
    double s = 0.0;
    double t =0.0;
    while (s <= arc_length)
    {
        s += vel_limit * dt +
                (1/2) * acc_limit * dt * dt +
                (1/6) * jerk_limit * dt * dt *dt;
        t += dt;
    }
    return  t;
}
/**************************************************************************************/
Vector4d evaluate_path(const double & t, const Vector6d & coeff)
{
    Vector4d p;
    p <<  coeff(0) + coeff(1) * t + coeff(2) * t*t + coeff(3) * t*t*t + coeff(4) * t*t*t*t + coeff(5) * t*t*t*t*t,
              coeff(1) + 2 * coeff(2) * t + 3 * coeff(3) * t*t + 4 * coeff(4) * t*t*t + 5 * coeff(5) * t*t*t*t,
              2 * coeff(2) + 6 * coeff(3) * t + 12 * coeff(4) * t*t + 20 * coeff(5) * t*t*t,
              6 * coeff(3) + 24 * coeff(4) * t + 60 * coeff(5) * t*t;
    return p;
}
