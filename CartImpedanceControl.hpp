#include "../include/CartesianControl.hpp"
#include "../include/DynamicsControl.hpp"
#include "../include/JointData.hpp"
#include <CosimaUtilities/SkewMatrix.hpp>
#include <CosimaUtilities/TorqueSaturation.hpp>
#include <CosimaUtilities/LowPassFilter.hpp>

class CartImpedanceControl : public RTT::TaskContext {
public:
    CartImpedanceControl(std::string const &name);

    bool configureHook();
    bool startHook();
    void updateHook();
    void stopHook();
    void cleanupHook();
    void compute();

    void setGains(float gainP, float gainD, float orientation_gainP, float orientation_gainD);
    void setGains_anisotropic(const Eigen::VectorXf& Stiffness_Des, const Eigen::VectorXf& Damping_Des);
    void addChain(int dof);
    void preparePorts();
    void displayCurrentState();
    double getSimulationTime();
    
    void quatToEuler(const double qW, const double qX, const double qY, const double qZ, Eigen::VectorXf &EulerAngle);
    void calculate_K_Max_Beta_Pos_Ori(Eigen::VectorXf &F_MAX, Eigen::VectorXf &Max_Disp, Eigen::VectorXf &Tau_MAX, Eigen::VectorXf &Max_Disp_Ori, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Beta ,Eigen::MatrixXf &K_Max_System_Ori, Eigen::VectorXf &Beta_Ori, int i);
    void get_K_Total_Pos_Ori(Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Var_Pos, Eigen::MatrixXf &K_Total_Des_Pos, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Var_Ori, Eigen::MatrixXf &K_Total_Des_Ori, int i);
    void calculate_F_Max (Eigen::VectorXf &Max_Disp, Eigen::VectorXf &F_Max, int i);
    void calculate_Tau_Max(Eigen::VectorXf &Max_Disp_Ori, Eigen::VectorXf &Tau_Max, int i);
    void boxShaping_Online(Eigen::VectorXf &Max_Disp, int i);
    void calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(Eigen::MatrixXf &K1, Eigen::MatrixXf &K2, Eigen::MatrixXf &Kd_Max, Eigen::VectorXf &Pmax, Eigen::VectorXf &Pmid, Eigen::VectorXf &Max_Dist_Vec, int i);
    void calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(Eigen::MatrixXf &K1_Ori, Eigen::MatrixXf &K2_Ori, Eigen::MatrixXf &Kd_Max_Ori, Eigen::VectorXf &Pmax_Ori, Eigen::VectorXf &Pmid_Ori, Eigen::VectorXf &Max_Dist_Vec_Ori, int i);
    void calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Eigen::MatrixXf &Kc, Eigen::MatrixXf &DampConstLinear, Eigen::VectorXf &FMidCtrl, Eigen::VectorXf &FkCtrl, int i);
    void calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Eigen::MatrixXf &Kc_Ori, Eigen::MatrixXf &DampConstAngular, Eigen::VectorXf &FMidCtrl_Ori, Eigen::VectorXf &FkCtrl_Ori, int i);
    void trajectory_execution(Eigen::VectorXf &DesPos, Eigen::VectorXf &EndEff_VelLinear, int i);
    void curPos_on_line_projection (Eigen::VectorXf &P1, Eigen::VectorXf &P2);

    template <typename T> static int sgn(const T &val, const T &eps = T(0)) {
        return (T(eps) < val) - (val < T(-eps));
    }

private:

    // DEFINING THE USE OF CARTESIAN CTRL, DYNAMICS CTRL AND JOINT DATA
    std::unique_ptr<CartesianControl> cc;
    std::unique_ptr<DynamicsControl> dc;
    std::unique_ptr<JointData> jd;

    // DEFINING JACOBIAN PSEUDO INVERSE & NULLSPACE & TORQUE SATURATION
    std::unique_ptr<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> > pinv;
    std::unique_ptr<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> > nullSpace;
    std::unique_ptr<CosimaUtilities::TorqueSaturation<Eigen::VectorXf> > trqSat;

    float thresh = 1e-4;

    // NUMBER OF END-EFFECTOR(S) INITIALIZATION
    int numOfEndEffectors;

    //INITIALIZING AND READING ELEMENTS FROM TXT FILE
    int numTrajAct, numTraj_CircTraj_Act , numBoxShapingAct;
    double numTrajSize, numTrajOmega, numCircTrajSize;
    double numVirtPosBound_x, numVirtPosBound_y, numVirtPosBound_z, numVirtOriBound_x, numVirtOriBound_y, numVirtOriBound_z, numDesPos_x, numDesPos_y, numDesPos_z;

    // SETTING CONSTANT GAINS - ELEMENTS
    double gainP, gainD, orientation_gainP, orientation_gainD;

    // CALCULATION OF DELTA-TIME - ELEMENTS
    double time0, time1, time2, deltaTime;

    // END-EFFECTOR CURRENT AND DESIRED POSE
    Eigen::VectorXf curPos3x1, desPos3x1, curOri3x1, desOri3x1;

    // END-EFFECTOR POSE ERROR
    Eigen::VectorXf position_error, angular_error;
    Eigen::VectorXf posOriError_6x1, posOriVelError_6x1;

    // END-EFFECTOR LINEAR AND ANGULAR VELOCITY ERROR
    Eigen::VectorXf velLinearError, velAngularError;

    // IMPEDANCE CONTROL GAIN: STIFFNESS
    Eigen::MatrixXf kConstPos, kVarPos, kTotalDesPos;
    Eigen::MatrixXf kConstOri, kVarOri, kTotalDesOri;
    Eigen::MatrixXf kTotal6x6;

    // IMPEDANCE CONTROL GAIN: DAMPING
    Eigen::MatrixXf dConstLinear, dVarPos;
    Eigen::MatrixXf dConstAngular, dVarOri;
    Eigen::MatrixXf dTotal6x6, dCritical6x6;

    // DESIRED STIFFNESS AND DAMPING INITIALIZATION
    Eigen::VectorXf Stiffness_Des, Damping_Des;

    // KdMAX - ELEMENTS
    Eigen::MatrixXf KdMax;
    Eigen::MatrixXf KdMax_Ori;

    // NATURAL FREQUENCY OBSERVER
    Eigen::MatrixXf natFreq_observer;

    // NULL-SPACE - ELEMENTS
    Eigen::VectorXf curJntPos7x1, desJntPos7x1;
    Eigen::VectorXf curJntVel7x1, desJntVel7x1;

    // NULL-SPACE - STIFFNESS - ELEMENTS
    Eigen::MatrixXf kConstNull, kVarNull;
    Eigen::MatrixXf Kd_NullSpace;

    // NULL-SPACE - DAMPING - ELEMENTS
    Eigen::MatrixXf dConstNull, dVarNull;

    // CALCULATION OF PREVIOUS VALUES - ELEMENTS
    Eigen::VectorXf preEndEff_PosError3x1, preEndEff_VelError3x1;
    Eigen::VectorXf preEndEff_OriError3x1, preEndEff_VelOriError3x1;

    // NULL-SPACE CTRL PLUS TORQUE SATURATION - ELEMENTS
    Eigen::VectorXf nullSpaceCtrl;
    Eigen::VectorXf torqueSaturation;

    // SYSTEM ENERGY - POTENTIAL AND KINETIC INITIALIZATION
    double potentialEnergy, kineticEnergy;

    // PROPOSED METHOD - POSITION - ELEMENTS
    Eigen::VectorXf maxVel3x1;
    Eigen::VectorXf Pmax3x1;
    Eigen::VectorXf Pmid3x1;
    Eigen::VectorXf maxDistVec_3x1;
    Eigen::MatrixXf K1_3x3;
    Eigen::MatrixXf K2_3x3;
    Eigen::MatrixXf Kc_3x3;
    Eigen::VectorXf F_midPointErrorCtrl;
    Eigen::VectorXf Fk;

    // PROPOSED METHOD - ORIENTATION - ELEMENTS
    Eigen::VectorXf maxVel3x1_Ori;
    Eigen::VectorXf Pmax3x1_Ori;
    Eigen::VectorXf Pmid3x1_Ori;
    Eigen::VectorXf maxDistVec_Ori_3x1;
    Eigen::MatrixXf K1_Ori_3x3;
    Eigen::MatrixXf K2_Ori_3x3;
    Eigen::MatrixXf Kc_Ori_3x3;
    Eigen::VectorXf F_midPointErrorCtrl_Ori;
    Eigen::VectorXf Fk_Ori;

    // DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - REQUIRED ELEMENTS
    Eigen::VectorXf Fd_Pos, Fd_Ori;

    // Fmax & BETA & MAXIMUM DISPLACEMENT (POSITION) & KmaxSystem
    Eigen::VectorXf Fmax;
    Eigen::VectorXf maxDisp;
    Eigen::VectorXf beta;
    Eigen::MatrixXf KmaxSystem;

    // TauMax &  BETA_Ori & MAXIMUM DISPLACEMENT (ORIENTATION) & KmaxSystem_Ori
    Eigen::VectorXf TauMax;
    Eigen::VectorXf maxDisp_Ori;
    Eigen::VectorXf beta_Ori;
    Eigen::MatrixXf KmaxSystem_Ori;

    // Fmax & TauMax SLOPE INITIALIZATION
    double slopeFmax, slopeTauMax_Ori;

    // ONLINE BOX SHAPING SWITCH
    double onlineBoxShapingSwitch;

    // TRAJECTORY EXECUTION - ELEMENTS
    double trajectory_activation;
    double omegaTraj, lineDispTrajCoef, circRadiusTraj;
    double ctrlStartTime, trajStartTime;

	// END-EFFECTOR TASK SPACE INERTIA & ESTIMATED JOINT TORQUES
    Eigen::MatrixXf lambda_c;
    RTT::OutputPort<Eigen::VectorXf> out_estimated_torque_port;
    Eigen::VectorXf estimated_torque_var;

    //JUST TESTING
    double stabilityVerfication;

    //JUST TESTING 2
    Eigen::VectorXf q_projCurPos, p1, p2;

};
