#include "../include/CartImpedanceControl.hpp"
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

#define _USE_MATH_DEFINES

CartImpedanceControl::CartImpedanceControl(std::string const &name) : RTT::TaskContext(name) {

    //OPERATIONS PREPARATION
    addOperation("setGains", &CartImpedanceControl::setGains, this).doc("set Gains");
    addOperation("setGains_anisotropic", &CartImpedanceControl::setGains_anisotropic, this).doc("set Gains anisotropic");
    addOperation("addChain", &CartImpedanceControl::addChain, this).doc("add Chain");
    addOperation("preparePorts", &CartImpedanceControl::preparePorts, this).doc("preparePorts");
    addOperation("displayCurrentState", &CartImpedanceControl::displayCurrentState, this).doc("display current state");
    addProperty("thresh", thresh);

    // DEFINING THE USE OF CARTESIAN CTRL, DYNAMICS CTRL AND JOINT DATA
    cc = std::make_unique<CartesianControl>(this);
    dc = std::make_unique<DynamicsControl>(this);
    jd = std::make_unique<JointData>(this);

    // NUMBER OF END-EFFECTOR(S) INITIALIZATION
    numOfEndEffectors = 0;

    //INITIALIZING AND READING ELEMENTS FROM TXT FILE
    numTrajAct = 0; numTraj_CircTraj_Act = 0; numBoxShapingAct = 0; numTrajSize = 0; numTrajOmega = 0; numCircTrajSize = 0;
    numVirtPosBound_x = 0; numVirtPosBound_y = 0; numVirtPosBound_z = 0;
    numVirtOriBound_x = 0; numVirtOriBound_y = 0; numVirtOriBound_z = 0;
    numDesPos_x = 0; numDesPos_y = 0; numDesPos_z = 0;

    std::fstream inFile;
    inFile.open("/home/kkouhkiloui/keyhan_ws/numbers.txt");
    inFile >> numTrajAct >> numTraj_CircTraj_Act >> numBoxShapingAct >>  numTrajSize >> numTrajOmega >> numCircTrajSize >>
                            numVirtPosBound_x >> numVirtPosBound_y >> numVirtPosBound_z >> numVirtOriBound_x >> numVirtOriBound_y >> numVirtOriBound_z >>
                            numDesPos_x >> numDesPos_y >> numDesPos_z;

    // SETTING CONSTANT GAINS - ELEMENTS
    gainP = 0; gainD = 0;
    orientation_gainP = 0; orientation_gainD = 0;

    // CALCULATION OF DELTA-TIME - ELEMENTS
    time1 = 0; time2= 0;
    deltaTime = 0;

    // END-EFFECTOR CURRENT AND DESIRED POSE
    curPos3x1 = Eigen::VectorXf (3); curPos3x1.setZero();
    desPos3x1 = Eigen::VectorXf (3); desPos3x1.setZero();
    curOri3x1 = Eigen::VectorXf (3); curOri3x1.setZero();
    desOri3x1 = Eigen::VectorXf (3); desOri3x1.setZero();

    // END-EFFECTOR POSE ERROR
    position_error = Eigen::VectorXf (3); position_error.setZero();
    angular_error = Eigen::VectorXf (3); angular_error.setZero();
    posOriError_6x1 = Eigen::VectorXf (6); posOriError_6x1.setZero();
    posOriVelError_6x1 = Eigen::VectorXf (6); posOriVelError_6x1.setZero();

    // END-EFFECTOR LINEAR AND ANGULAR VELOCITY ERROR
    velLinearError = Eigen::VectorXf (3); velLinearError.setZero();
    velAngularError = Eigen::VectorXf (3); velAngularError.setZero();

    // IMPEDANCE CONTROL GAIN: STIFFNESS
    kConstPos = Eigen::MatrixXf(3, 3); kConstPos.setIdentity();
    kVarPos = Eigen::MatrixXf(3, 3); kVarPos.setIdentity();
    kTotalDesPos = Eigen::MatrixXf(3, 3); kTotalDesPos.setIdentity();
    kConstOri = Eigen::MatrixXf(3, 3); kConstOri.setIdentity();
    kVarOri = Eigen::MatrixXf(3, 3); kVarOri.setIdentity();
    kTotalDesOri = Eigen::MatrixXf(3, 3); kTotalDesOri.setIdentity();
    kTotal6x6 = Eigen::MatrixXf (6,6); kTotal6x6.setIdentity();

    // IMPEDANCE CONTROL GAIN: DAMPING
    dConstLinear = Eigen::MatrixXf (3,3); dConstLinear.setIdentity();
    dVarPos = Eigen::MatrixXf(3, 3); dVarPos.setIdentity();
    dConstAngular = Eigen::MatrixXf (3, 3); dConstAngular.setIdentity();
    dVarOri = Eigen::MatrixXf(3, 3); dVarOri.setIdentity();
    dTotal6x6 = Eigen::MatrixXf (6,6); dTotal6x6.setIdentity();
    dCritical6x6 = Eigen::MatrixXf (6,6); dCritical6x6.setIdentity();

    // DESIRED STIFFNESS AND DAMPING INITIALIZATION
    Stiffness_Des = Eigen::VectorXf(6); Stiffness_Des.setZero();
    Damping_Des = Eigen::VectorXf(6); Damping_Des.setZero();

    // KdMAX - ELEMENTS
    KdMax = Eigen::MatrixXf(3, 3); KdMax.setIdentity();
    KdMax_Ori = Eigen::MatrixXf(3, 3); KdMax_Ori.setIdentity();

    // NATURAL FREQUENCY OBSERVER
    natFreq_observer = Eigen::MatrixXf (6,6); natFreq_observer.setIdentity();

    // NULL-SPACE - ELEMENTS
    curJntPos7x1 = Eigen::VectorXf (7); curJntPos7x1.setZero();
    desJntPos7x1 = Eigen::VectorXf (7); desJntPos7x1.setZero();
    curJntVel7x1 = Eigen::VectorXf (7); curJntVel7x1.setZero();
    desJntVel7x1 = Eigen::VectorXf (7); desJntVel7x1.setZero();

    // NULL-SPACE - STIFFNESS - ELEMENTS
    kConstNull = Eigen::MatrixXf (7,7); kConstNull.setIdentity();
    kVarNull = Eigen::MatrixXf (7,7); kVarNull.setIdentity();
    Kd_NullSpace = Eigen::MatrixXf (7,7); Kd_NullSpace.setIdentity();

    // NULL-SPACE - DAMPING - ELEMENTS
    dConstNull = Eigen::MatrixXf (7,7); dConstNull.setIdentity();
    dVarNull = Eigen::MatrixXf (7,7); dVarNull.setIdentity();

    // CALCULATION OF PREVIOUS VALUES - ELEMENTS
    preEndEff_PosError3x1 = Eigen::VectorXf(3); preEndEff_PosError3x1.setZero();
    preEndEff_OriError3x1 = Eigen::VectorXf(3); preEndEff_OriError3x1.setZero();
    preEndEff_VelError3x1 = Eigen::VectorXf(3); preEndEff_VelError3x1.setZero();
    preEndEff_VelOriError3x1 = Eigen::VectorXf(3); preEndEff_VelOriError3x1.setZero();

    // NULL-SPACE CTRL PLUS TORQUE SATURATION - ELEMENTS
    nullSpaceCtrl = Eigen::VectorXf (7); nullSpaceCtrl.setZero();
    torqueSaturation = Eigen::VectorXf (7); torqueSaturation.setZero();

    // SYSTEM ENERGY - POTENTIAL AND KINETIC INITIALIZATION
    potentialEnergy = 0; kineticEnergy = 0;

    // PROPOSED METHOD - POSITION - ELEMENTS
    maxVel3x1 = Eigen::VectorXf(3); maxVel3x1.setZero();
    Pmax3x1 = Eigen::VectorXf(3); Pmax3x1.setZero();
    Pmid3x1 = Eigen::VectorXf(3); Pmid3x1.setZero();
    maxDistVec_3x1 = Eigen::VectorXf(3); maxDistVec_3x1.setZero();
    maxDistVec_3x1.setConstant(0.001);
    K1_3x3 = Eigen::MatrixXf(3,3); K1_3x3.setIdentity();
    K2_3x3 = Eigen::MatrixXf(3,3); K2_3x3.setIdentity();
    Kc_3x3 = Eigen::MatrixXf(3,3); Kc_3x3.setZero();
    F_midPointErrorCtrl = Eigen::VectorXf(3); F_midPointErrorCtrl.setZero();
    Fk = Eigen::VectorXf(3); Fk.setZero();

    // PROPOSED METHOD - ORIENTATION - ELEMENTS
    maxVel3x1_Ori = Eigen::VectorXf(3); maxVel3x1_Ori.setZero();
    Pmax3x1_Ori = Eigen::VectorXf(3);Pmax3x1_Ori.setZero();
    Pmid3x1_Ori = Eigen::VectorXf(3);Pmid3x1_Ori.setZero();
    maxDistVec_Ori_3x1 = Eigen::VectorXf(3); maxDistVec_Ori_3x1.setZero();
    K1_Ori_3x3 = Eigen::MatrixXf(3,3); K1_Ori_3x3.setIdentity();
    K2_Ori_3x3 = Eigen::MatrixXf(3,3); K2_Ori_3x3.setIdentity();
    Kc_Ori_3x3 = Eigen::MatrixXf(3,3); Kc_Ori_3x3.setIdentity();
    F_midPointErrorCtrl_Ori = Eigen::VectorXf(3); F_midPointErrorCtrl_Ori.setZero();
    Fk_Ori = Eigen::VectorXf(3); Fk_Ori.setZero();

    // DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - REQUIRED ELEMENTS
    Fd_Pos = Eigen::VectorXf(3); Fd_Pos.setZero();
    Fd_Ori = Eigen::VectorXf(3); Fd_Ori.setZero();

    // Fmax & BETA & MAXIMUM DISPLACEMENT (POSITION) & KmaxSystem
    Fmax = Eigen::VectorXf(3); Fmax.setZero();
    maxDisp = Eigen::VectorXf(3);
    beta = Eigen::VectorXf(3); beta.setZero();
    KmaxSystem = Eigen::MatrixXf(3,3); KmaxSystem.setIdentity();
    maxDisp[0] = numVirtPosBound_x;
    maxDisp[1] = numVirtPosBound_y;
    maxDisp[2] = numVirtPosBound_z;

    // TauMax &  BETA_Ori & MAXIMUM DISPLACEMENT (ORIENTATION) & KmaxSystem_Ori
    TauMax = Eigen::VectorXf(3); TauMax.setZero();
    maxDisp_Ori = Eigen::VectorXf(3);
    beta_Ori = Eigen::VectorXf(3); beta_Ori.setZero();
    KmaxSystem_Ori = Eigen::MatrixXf(3,3); KmaxSystem_Ori.setIdentity();
    maxDisp_Ori[0] = numVirtOriBound_x;
    maxDisp_Ori[1] = numVirtOriBound_y;
    maxDisp_Ori[2] = numVirtOriBound_z;

    // Fmax & TauMax SLOPE INITIALIZATION
    slopeFmax = 0; slopeTauMax_Ori = 0;

    // ONLINE BOX SHAPING SWITCH
    onlineBoxShapingSwitch = numBoxShapingAct;
    time0 = 0;

    // TRAJECTORY EXECUTION - ELEMENTS
    trajectory_activation = numTrajAct;
    lineDispTrajCoef = numTrajSize;
    circRadiusTraj = numCircTrajSize;
    omegaTraj = numTrajOmega;
    trajStartTime = 10;
    ctrlStartTime = 10;

    //JUST TESTING
    maxVel3x1.setConstant(1);
    stabilityVerfication = 0;

    //JUST TESTING 2
    q_projCurPos = Eigen::VectorXf(3); q_projCurPos.setZero();
    p1 = Eigen::VectorXf(3); p1.setZero();
    p2 = Eigen::VectorXf(3); p2.setZero();

    p1[0] = 0.4; p1[1] = 0; p1[2] = 1;
    p2[0] = 0.4; p2[1] = 0; p2[2] = -1;

}

void CartImpedanceControl::addChain(int dof) {
    numOfEndEffectors++;

    cc->addChain(dof);
    dc->addChain(dof);
    jd->addChain(dof);
    estimated_torque_var.resize(dof + estimated_torque_var.size());
    estimated_torque_var.setZero();
}

void CartImpedanceControl::preparePorts() {
    cc->preparePorts();
    dc->preparePorts();
    jd->preparePorts();

    Eigen::MatrixXf tmp2(6 * numOfEndEffectors, 6 * numOfEndEffectors);
    pinv = std::make_unique<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> >(tmp2);

    Eigen::MatrixXf tmp3(7 * numOfEndEffectors, 7 * numOfEndEffectors);
    nullSpace = std::make_unique<CosimaUtilities::Pseudoinverse<Eigen::MatrixXf> >(tmp3);

    Eigen::VectorXf tmp4(7 * numOfEndEffectors);
    tmp4.setZero();
    trqSat = std::make_unique<CosimaUtilities::TorqueSaturation<Eigen::VectorXf> > (0.9f,tmp4);

    out_estimated_torque_port.setName("out_estimated_torque_port");
    out_estimated_torque_port.doc("Output port for sending estimated force values");
    out_estimated_torque_port.setDataSample(estimated_torque_var);
    ports()->addPort(out_estimated_torque_port);
}

void CartImpedanceControl::displayCurrentState() {
    cc->displayCurrentState();
    dc->displayCurrentState();
    jd->displayCurrentState();
    return;
}

void CartImpedanceControl::setGains_anisotropic(const Eigen::VectorXf& Stiffness_Des, const Eigen::VectorXf& Damping_Des) {
    assert(Stiffness_Des[0] >= 0); assert(Damping_Des[0] >= 0);
    assert(Stiffness_Des[1] >= 0); assert(Damping_Des[1] >= 0);
    assert(Stiffness_Des[2] >= 0); assert(Damping_Des[2] >= 0);
    assert(Stiffness_Des[3] >= 0); assert(Damping_Des[3] >= 0);
    assert(Stiffness_Des[4] >= 0); assert(Damping_Des[4] >= 0);
    assert(Stiffness_Des[5] >= 0); assert(Damping_Des[5] >= 0);

    this->Stiffness_Des[0] = Stiffness_Des[0];
    this->Stiffness_Des[1] = Stiffness_Des[1];
    this->Stiffness_Des[2] = Stiffness_Des[2];
    this->Stiffness_Des[3] = Stiffness_Des[3];
    this->Stiffness_Des[4] = Stiffness_Des[4];
    this->Stiffness_Des[5] = Stiffness_Des[5];

    this->Damping_Des[0] = Damping_Des[0];
    this->Damping_Des[1] = Damping_Des[1];
    this->Damping_Des[2] = Damping_Des[2];
    this->Damping_Des[3] = Damping_Des[3];
    this->Damping_Des[4] = Damping_Des[4];
    this->Damping_Des[5] = Damping_Des[5];
}


void CartImpedanceControl::setGains(float gainP, float gainD, float orientation_gainP, float orientation_gainD) {
    assert(gainP >= 0);
    assert(gainD >= 0);
    assert(orientation_gainP >= 0);
    assert(orientation_gainD >= 0);
    this->gainP = gainP;
    this->gainD = gainD;
    this->orientation_gainP = orientation_gainP;
    this->orientation_gainD = orientation_gainD;
}

double CartImpedanceControl::getSimulationTime() {
    return 1E-9 * RTT::os::TimeService::ticks2nsecs(RTT::os::TimeService::Instance()->getTicks());

}

void CartImpedanceControl::quatToEuler(const double qW, const double qX, const double qY, const double qZ, Eigen::VectorXf &EulerAngle){

    double sinr_cosp, cosr_cosp, sinp, siny_cosp, cosy_cosp;
    double roll, pitch, yaw;

    sinr_cosp = 0; cosr_cosp = 0; sinp = 0; siny_cosp = 0; cosy_cosp = 0;
    roll = 0; pitch = 0; yaw = 0;

    // ROLL
    sinr_cosp = +2.0 * (qW * qX + qY * qZ);
    cosr_cosp = +1.0 - 2.0 * (qX * qX + qY * qY);
    roll = atan2(sinr_cosp, cosr_cosp);

    // PITCH
    sinp = +2.0 * (qW * qY - qZ * qX);
    if (fabs(sinp) >= 1)
        pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
    else
        pitch = asin(sinp);

    // YAW
    siny_cosp = +2.0 * (qW * qZ + qX * qY);
    cosy_cosp = +1.0 - 2.0 * (qY * qY + qZ * qZ);
    yaw = atan2(siny_cosp, cosy_cosp);

    EulerAngle[0] = roll;
    EulerAngle[1] = pitch;
    EulerAngle[2] = yaw;
}

void CartImpedanceControl::calculate_K_Max_Beta_Pos_Ori(Eigen::VectorXf &F_MAX, Eigen::VectorXf &Max_Disp, Eigen::VectorXf &Tau_MAX, Eigen::VectorXf &Max_Disp_Ori, Eigen::MatrixXf &K_Max_System, Eigen::VectorXf &Beta ,Eigen::MatrixXf &K_Max_System_Ori, Eigen::VectorXf &Beta_Ori, int i) {

    K_Max_System(i,i) = F_MAX[i]/Max_Disp[i];
    Beta[i] = sqrt(log(K_Max_System(i,i) - kConstPos(i,i))/pow(Max_Disp[i],2));
    if (std::isnan(Beta[i])){
        Beta[i] = 0;
    }

    K_Max_System_Ori(i,i) = Tau_MAX[i]/Max_Disp_Ori[i];
    Beta_Ori[i] = sqrt(log(K_Max_System_Ori(i,i) - kConstOri(i,i))/pow(Max_Disp_Ori[i],2));
    if (std::isnan(Beta_Ori[i])){
        Beta_Ori[i] = 0;
    }

}

void CartImpedanceControl::get_K_Total_Pos_Ori(Eigen::VectorXf &Beta, Eigen::MatrixXf &K_Var_Pos, Eigen::MatrixXf &K_Total_Des_Pos, Eigen::VectorXf &Beta_Ori, Eigen::MatrixXf &K_Var_Ori, Eigen::MatrixXf &K_Total_Des_Ori, int i) {

    K_Var_Pos(i,i) = exp(pow(Beta[i] * (position_error[i]), 2)) - 1;
    K_Total_Des_Pos(i,i) = kConstPos(i,i) + K_Var_Pos(i,i);

    K_Var_Ori(i, i) = exp(pow(Beta_Ori[i] * (angular_error[i]),2)) - 1;
    K_Total_Des_Ori(i,i) = kConstOri(i,i) + K_Var_Ori(i,i);

}

void CartImpedanceControl::calculate_F_Max(Eigen::VectorXf &Max_Disp, Eigen::VectorXf &F_Max, int i){

    // FRANKA CONTROLLER CALIBRATION
    /*if (i == 0){
        if (Max_Disp[i] >= 0.003){
            F_Max[i] = 30;
        }else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.003) {
            F_Max[i] = 15;
        }else {
            F_Max[i] = 12.5;
        }
    }else if (i == 1){

        if (Max_Disp[i] >= 0.09){
            F_Max[i] = 30;
        }else if(Max_Disp[i] >= 0.05 && Max_Disp[i] < 0.09) {
            F_Max[i] = 25;
        }else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.05) {
            F_Max[i] = 22.5;
        }else if(Max_Disp[i] >= 0.03 && Max_Disp[i] < 0.04) {
            F_Max[i] = 20;
        }else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.03) {
            F_Max[i] = 18;
        }else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.02) {
            F_Max[i] = 12.5;
        }else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007) {
            F_Max[i] = 10;
        }else if(Max_Disp[i] >= 0.004 && Max_Disp[i] < 0.006) {
            F_Max[i] = 6;
        }else if(Max_Disp[i] >= 0.003 && Max_Disp[i] < 0.004) {
            F_Max[i] = 5;
        }else if(Max_Disp[i] >= 0.002 && Max_Disp[i] < 0.003) {
            F_Max[i] = 4;
        }else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.002) {
            F_Max[i] = 3;
        }else {
            F_Max[i] = 2;
        }
    }else{

        if (Max_Disp[i] >= 0.08){
            F_Max[i] = 30;
        }else if(Max_Disp[i] >= 0.05 && Max_Disp[i] < 0.08) {
            F_Max[i] = 28;
        }else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.05) {
            F_Max[i] = 26;
        }else if(Max_Disp[i] >= 0.03 && Max_Disp[i] < 0.04) {
            F_Max[i] = 24;
        }else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.03) {
            F_Max[i] = 20;
        }else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.02) {
            F_Max[i] = 12.5;
        }else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007) {
            F_Max[i] = 10;
        }else if(Max_Disp[i] >= 0.004 && Max_Disp[i] < 0.006) {
            F_Max[i] = 6;
        }else if(Max_Disp[i] >= 0.003 && Max_Disp[i] < 0.004) {
            F_Max[i] = 5;
        }else if(Max_Disp[i] >= 0.002 && Max_Disp[i] < 0.003) {
            F_Max[i] = 4;
        }else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.002) {
            F_Max[i] = 3;
        }else {
            F_Max[i] = 2;
        }
    }*/


    // KUKA CONTROLLER CALIBRATION
    if (maxDisp[i] >= 0.08){
        slopeFmax = 50;
        F_Max[i] = (-slopeFmax * (0.08-Max_Disp[i]) + 11)/1;

    }else if(Max_Disp[i] >= 0.06 && Max_Disp[i] < 0.08){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.06-Max_Disp[i]) + 9)/1;

    }else if(Max_Disp[i] >= 0.04 && Max_Disp[i] < 0.06){
        slopeFmax = 50;
        F_Max[i] = (-slopeFmax * (0.04-Max_Disp[i]) + 8)/1;

    }else if(Max_Disp[i] >= 0.02 && Max_Disp[i] < 0.04){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.02-Max_Disp[i]) + 6)/1;

    }else if(Max_Disp[i] >= 0.01 && Max_Disp[i] < 0.02){
        slopeFmax = 100;
        F_Max[i] = (-slopeFmax * (0.01-Max_Disp[i]) + 5)/1;

    }else if(Max_Disp[i] >= 0.007 && Max_Disp[i] < 0.01){
        slopeFmax = 333;
        F_Max[i] = (-slopeFmax * (0.007-Max_Disp[i]) + 4)/1;

    }else if(Max_Disp[i] >= 0.006 && Max_Disp[i] < 0.007){
        slopeFmax = 1000;
        F_Max[i] = (-slopeFmax * (0.006-Max_Disp[i]) + 3)/1;

    }else if(Max_Disp[i] >= 0.001 && Max_Disp[i] < 0.006){
        slopeFmax = 250;
        F_Max[i] = -slopeFmax * (0.001-Max_Disp[i]) + 1;

    }else
        F_Max[i] = 1;
}

void CartImpedanceControl::calculate_Tau_Max(Eigen::VectorXf &Max_Disp_Ori, Eigen::VectorXf &Tau_Max, int i){

    // FRANKA CONTROLLER CALIBRATION
    /*if (i == 0){
        Tau_Max[i] = 30;
    }
    else{
        if (Max_Disp_Ori[i] >= 0.261799){
            Tau_Max[i] = 20;
        }else if (Max_Disp_Ori[i] >= 0.0872665 && Max_Disp_Ori[i] < 0.261799){
            Tau_Max[i] = 15;
        }else {
            Tau_Max[i] = 5;
        }
    }*/

    // KUKA CONTROLLER CALIBRATION
    if (Max_Disp_Ori[i] >= 0.08){
        slopeTauMax_Ori = 50;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.08-Max_Disp_Ori[i]) + 11)/2;

    }else if(Max_Disp_Ori[i] >= 0.06 && Max_Disp_Ori[i] < 0.08){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.06-Max_Disp_Ori[i]) + 9)/2;

    }else if(Max_Disp_Ori[i] >= 0.04 && Max_Disp_Ori[i] < 0.06){
        slopeTauMax_Ori = 50;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.04-Max_Disp_Ori[i]) + 8)/2;

    }else if(Max_Disp_Ori[i] >= 0.02 && Max_Disp_Ori[i] < 0.04){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.02-Max_Disp_Ori[i]) + 6)/2;

    }else if(Max_Disp_Ori[i] >= 0.01 && Max_Disp_Ori[i] < 0.02){
        slopeTauMax_Ori = 100;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.01-maxDisp_Ori[i]) + 5)/2;

    }else if(maxDisp_Ori[i] >= 0.007 && maxDisp_Ori[i] < 0.01){
        slopeTauMax_Ori = 333;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.007-maxDisp_Ori[i]) + 4)/2;

    }else if(maxDisp_Ori[i] >= 0.006 && maxDisp_Ori[i] < 0.007){
        slopeTauMax_Ori = 1000;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.006-maxDisp_Ori[i]) + 3)/2;

    }else if(maxDisp_Ori[i] >= 0.001 && maxDisp_Ori[i] < 0.006){
        slopeTauMax_Ori = 250;
        Tau_Max[i] = (-slopeTauMax_Ori * (0.001-maxDisp_Ori[i]) + 1)/2;

    }else
        Tau_Max[i] = 1.5;
}

void CartImpedanceControl::boxShaping_Online(Eigen::VectorXf &Max_Disp, int i){

    if (getSimulationTime() >= time0 + 1){
            Max_Disp[i] = Max_Disp[i] - 0.0025;
            if(i == 2){
                time0 = getSimulationTime();
            }
            if(Max_Disp[i] <= 0.01) {
                Max_Disp[i] = 0.01;
            }
    }
}

void CartImpedanceControl::trajectory_execution(Eigen::VectorXf &DesPos, Eigen::VectorXf &EndEff_VelLinear, int i) {

    if (cc->des_poses_var[0] != numDesPos_x || cc->des_poses_var[1] != numDesPos_y ||cc->des_poses_var[2] != numDesPos_z) {
        numDesPos_x = cc->des_poses_var[0];
        numDesPos_y = cc->des_poses_var[1];
        numDesPos_z = cc->des_poses_var[2];
    }

    if (numTraj_CircTraj_Act == 0) {
        if (getSimulationTime() >= ctrlStartTime + trajStartTime) {
            if (i == 0) {
                DesPos[0] = numDesPos_x;
                EndEff_VelLinear[0] = -cc->twists_var[0];
            } else if (i == 1) {
                DesPos[1] = numDesPos_y + lineDispTrajCoef * sin(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
                EndEff_VelLinear[1] = -cc->twists_var[1] + omegaTraj * 2 * 3.14 * lineDispTrajCoef * cos(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
            } else {
                DesPos[2] = numDesPos_z;
                EndEff_VelLinear[2] = -cc->twists_var[2];
            }
        }
    } else {
        if (getSimulationTime() >= ctrlStartTime + trajStartTime) {
            if (i == 0) {
                DesPos[0] = numDesPos_x + circRadiusTraj * cos(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
                EndEff_VelLinear[0] = -cc->twists_var[0] - omegaTraj * 2 * 3.14 * circRadiusTraj * sin(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
            } else if (i == 1) {
                DesPos[1] = numDesPos_y + circRadiusTraj * sin(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
                EndEff_VelLinear[1] = -cc->twists_var[1] + omegaTraj * 2 * 3.14 * circRadiusTraj * cos(omegaTraj * 2 * 3.14 * (getSimulationTime() - (ctrlStartTime + trajStartTime)));
            } else {
                DesPos[2] = numDesPos_z;
                EndEff_VelLinear[2] = -cc->twists_var[2];
            }
        }
    }
}

void CartImpedanceControl::calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(Eigen::MatrixXf &K1, Eigen::MatrixXf &K2, Eigen::MatrixXf &Kd_Max, Eigen::VectorXf &Pmax, Eigen::VectorXf &Pmid, Eigen::VectorXf &Max_Dist_Vec, int i) {

    if (fabs(position_error[i]) > fabs(preEndEff_PosError3x1[i])) {
        Max_Dist_Vec[i] = curPos3x1[i] - desPos3x1[i];
        Kd_Max(i,i) = kConstPos(i,i) + kVarPos(i, i);
    } else{

        Max_Dist_Vec[i] = Max_Dist_Vec[i];
        Kd_Max(i,i) = Kd_Max(i,i);
    }

    // DEFINING Pmax3x1 and Pmid3x1
    Pmax[i] = desPos3x1[i] + Max_Dist_Vec[i];
    Pmid[i] = desPos3x1[i] + (Pmax[i] - desPos3x1[i])/2;

    // DEFINING K1_3x3 AND K2_3x3 w.r.t MAXIMUM STIFFNESS AT THE RELEASED POINT
    K1(i, i) = 2 * Kd_Max(i, i);
    K2(i, i) = 2 * Kd_Max(i, i);

}

void CartImpedanceControl::calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(Eigen::MatrixXf &K1_Ori, Eigen::MatrixXf &K2_Ori, Eigen::MatrixXf &Kd_Max_Ori, Eigen::VectorXf &Pmax_Ori, Eigen::VectorXf &Pmid_Ori, Eigen::VectorXf &Max_Dist_Vec_Ori, int i) {

    if ( fabs(angular_error[i]) > fabs(preEndEff_OriError3x1[i])) {
        Max_Dist_Vec_Ori[i] = curOri3x1[i] - desOri3x1[i];
        Kd_Max_Ori(i,i) = kConstOri(i,i) + kVarOri(i, i);
    } else{

        Max_Dist_Vec_Ori[i] = Max_Dist_Vec_Ori[i];
        Kd_Max_Ori(i,i) = Kd_Max_Ori(i,i);
    }

    // DEFINING Pmax3x1 and Pmid3x1
    Pmax_Ori[i] = desOri3x1[i] + Max_Dist_Vec_Ori[i];
    Pmid_Ori[i] = desOri3x1[i] + (Pmax_Ori[i] - desOri3x1[i])/2;

    // DEFINING K1_Ori_3x3 AND K2_Ori_3x3 w.r.t MAXIMUM STIFFNESS AT THE RELEASED POINT
    K1_Ori(i, i) = 2 * Kd_Max_Ori(i, i);
    K2_Ori(i, i) = 2 * Kd_Max_Ori(i, i);
}

void CartImpedanceControl::calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Eigen::MatrixXf &Kc, Eigen::MatrixXf &DampConstLinear, Eigen::VectorXf &FMidCtrl, Eigen::VectorXf &FkCtrl, int i) {

    if(getSimulationTime() >= ctrlStartTime){

        if(sgn(position_error[i],0.005f) == sgn(velLinearError[i],0.005f)){
            Kc.setZero();
            DampConstLinear(i,i) = 0;
            FMidCtrl[i] = 0;

        }else {
            Kc = K2_3x3;
            FMidCtrl = Kc * (Pmid3x1 - curPos3x1);
            FkCtrl[i] = 0;
        }
    }
}

void CartImpedanceControl::calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Eigen::MatrixXf &Kc_Ori, Eigen::MatrixXf &DampConstAngular, Eigen::VectorXf &FMidCtrl_Ori, Eigen::VectorXf &FkCtrl_Ori, int i) {

    if(getSimulationTime() >= ctrlStartTime){

        if(sgn(angular_error[i]) == sgn(velAngularError[i])){
            Kc_Ori.setZero();
            DampConstAngular(i,i) = 0;
            FMidCtrl_Ori[i] = 0;

        }else {
            Kc_Ori = K2_3x3;
            FMidCtrl_Ori = Kc_Ori * (Pmid3x1_Ori - curOri3x1);
            FkCtrl_Ori[i] = 0;
        }
    }
}

void CartImpedanceControl::curPos_on_line_projection(Eigen::VectorXf &P1, Eigen::VectorXf &P2) {
    q_projCurPos[0] = (P1(0)*pow(P2(2),2) + ((-P2(0) - P1(0))*P1(2) + (P2(0) - P1(0))*curPos3x1(2))*P2(2) + P2(0)*pow(P1(2),2) +
                        (P1(0) - P2(0))*curPos3x1(2)*P1(2) + P1(0)*pow(P2(1),2) + ((-P2(0) - P1(0))*P1(1) + (P2(0) - P1(0))*curPos3x1 (1))*P2(1) + P2(0)*pow(P1(1),2) +
                        (P1(0) - P2(0))*curPos3x1(1)*P1(1) + curPos3x1(0)*pow(P2(0),2) - 2*curPos3x1(0)*P1(0)*P2(0) + curPos3x1(0)*pow(P1(0),2))/
                                (pow(P2(2),2) - 2*P1(2)*P2(2) + pow(P1(2),2)+pow(P2(1),2) - 2*P1(1)*P2(1) + pow(P1(1),2) + pow(P2(0),2) - 2*P1(0)*P2(0) + pow(P1(0),2));

    q_projCurPos[1] = (P1(1)*pow(P2(2),2) + ((-P2(1) - P1(1))*P1(2) + (P2(1) - P1(1))*curPos3x1(2))*P2(2) + P2(1)*pow(P1(2),2) + (P1(1) - P2(1))*curPos3x1(2)*P1(2) + curPos3x1(1)*pow(P2(1),2) +
                        (-2*curPos3x1(1)*P1(1) + (curPos3x1(0) - P1(0)*P2(0)) + pow(P1(0),2) - curPos3x1(0)*P1(0))*P2(1) + curPos3x1(1)*pow(P1(1),2) + (pow(P2(0),2) + (-P1(0) - curPos3x1(0))*P2(0)
                        + curPos3x1(0)*P1(0))*P1(1))/(pow(P2(2),2) - 2*P1(2)*P2(2) + pow(P1(2),2)+pow(P2(1),2) - 2*P1(1)*P2(1) + pow(P1(1),2) + pow(P2(0),2) - 2*P1(0)*P2(0) + pow(P1(0),2));

    q_projCurPos[2] = (curPos3x1(2)*pow(P2(2),2) + (-2*curPos3x1(2)*P1(2) + (curPos3x1(1) - P1(1))*P2(1) + pow(P1(1),2) - curPos3x1(1)*P1(1) +
                        (curPos3x1(0) - P1(0))*P2(0) + pow(P1(0),2) - curPos3x1(0)*P1(0))*P2(2) + curPos3x1(2)*pow(P1(2),2) + (pow(P2(1),2) +
                                (-P1(1) - curPos3x1(1)) *P2(1) + curPos3x1(1)*P1(1) + pow(P2(0),2) + (-P1(0) - curPos3x1(0))*P2(0) + curPos3x1(0)*P1(0))*P1(2))/
                                        (pow(P2(2),2) - 2*P1(2)*P2(2) + pow(P1(2),2)+pow(P2(1),2) - 2*P1(1)*P2(1) + pow(P1(1),2) + pow(P2(0),2) - 2*P1(0)*P2(0) + pow(P1(0),2));
}

void CartImpedanceControl::compute() {

    Eigen::VectorXf error(6*cc->number_of_end_effectors);
    Eigen::VectorXf vel_error(6*cc->number_of_end_effectors);

    Eigen::MatrixXf inertia_c_inv = dc->in_inertia_var.inverse();
    lambda_c = pinv->compute(cc->in_jacobian_var * inertia_c_inv * cc->in_jacobian_var.transpose(), thresh);

    kConstPos(0,0) = Stiffness_Des[0]; kConstPos(1,1) = Stiffness_Des[1]; kConstPos(2,2) = Stiffness_Des[2];
    kConstOri(0,0) = Stiffness_Des[3]; kConstOri(1,1) = Stiffness_Des[4]; kConstOri(2,2) = Stiffness_Des[5];
    dConstLinear(0,0) = Damping_Des[0]; dConstLinear(1,1) = Damping_Des[1]; dConstLinear(2,2) = Damping_Des[2];
    dConstAngular(0,0) = Damping_Des[3]; dConstAngular(1,1) = Damping_Des[4]; dConstAngular(2,2) = Damping_Des[5];

    // NULLSPACE CTRL ELEMENTS - GAINS AND DESIRED JOINT VALUES
    float gainNullP, gainNullD;
    gainNullP = 100;
    gainNullD = 20;

    desJntPos7x1[0] = -0.154303;
    desJntPos7x1[1] = -0.300414;
    desJntPos7x1[2] = 0.222139;
    desJntPos7x1[3] = 1.21592;
    desJntPos7x1[4] = 0.0785831;
    desJntPos7x1[5] = 1.27652;
    desJntPos7x1[6] = 1.12973;

    // DELTA TIME CALCULATION
    deltaTime = time2 - time1;

    // ROBOT CONTROL ELEMENTS
    for (int i = 0; i < cc->number_of_end_effectors; i++) {

        int pose_index = i * 7;

        if(cc->poses_var.segment<4>(pose_index+3).transpose()*cc->des_poses_var.segment<4>(pose_index+3)>0){

            cc->des_poses_var.segment<4>(pose_index+3) = -cc->des_poses_var.segment<4>(pose_index+3);
        }
        Eigen::Matrix3f skew;
        CosimaUtilities::skewMatrix(cc->poses_var.segment<3>(pose_index + 4), skew);
        angular_error = cc->des_poses_var(pose_index + 3) * cc->poses_var.segment<3>(pose_index + 4) - cc->poses_var(pose_index + 3) * cc->des_poses_var.segment<3>(pose_index + 4) - skew * cc->des_poses_var.segment<3>(pose_index + 4);

        vel_error.segment<6>(i*6) = cc->des_twists_var.segment<6>(i * 6) - cc->twists_var.segment<6>(i * 6);

        for (int i = 0; i < 3; i++) {

            // TRAJECTORY ACTIVATION
            if (trajectory_activation == 1){
                trajectory_execution(desPos3x1, velLinearError, i);
            }

            // CURRENT END-EFFECTOR POSITION
            curPos3x1[i] = cc->poses_var[i];

            // POSITION AND ORIENTATION ERROR
            position_error[i] = desPos3x1[i] - curPos3x1[i];

            // COMPUTATION OF END-EFFECTOR LINEAR/ANGULAR VELOCITY ERRORS
            velLinearError[i] = cc->des_twists_var[i] - cc->twists_var[i];
            velAngularError[i] = cc->des_twists_var[i + 3] - cc->twists_var[i + 3];

            // ONLINE BOX SHAPING
            if(onlineBoxShapingSwitch == 0){

                // FMAX AND TAU_MAX CALCULATION
                calculate_F_Max (maxDisp, Fmax, i);
                calculate_Tau_Max (maxDisp_Ori, TauMax, i);

            } else
            {
                boxShaping_Online (maxDisp, i);
                calculate_F_Max (maxDisp, Fmax, i);
                calculate_Tau_Max (maxDisp_Ori, TauMax, i);

            }

            // K_SYSTEM_MAX, MAXIMUM DISPLACEMENT AND BETA CALCULATION FOR BOTH END-EEFECTOR POSITION AND ORIENTATION
            calculate_K_Max_Beta_Pos_Ori(Fmax, maxDisp, TauMax, maxDisp_Ori, KmaxSystem, beta, KmaxSystem_Ori, beta_Ori, i);

            // K_TOTAL AND K_VAR CALCULATION FOR BOTH POSITION AND OREINTATION
            get_K_Total_Pos_Ori(beta, kVarPos, kTotalDesPos, beta_Ori, kVarOri, kTotalDesOri, i);

            // APPLYING THE SYSTEM PHYSICAL LIMITATION CONDITION
            if(getSimulationTime() >= ctrlStartTime){

                if (kConstPos(i,i) + kVarPos(i,i) > KmaxSystem(i,i) && beta[i] > 0){
                    kVarPos(i,i) = Fmax[i]/fabs(position_error[i]) - kConstPos(i,i);
                }else if(kTotalDesPos(i,i) > KmaxSystem(i,i)){
                    kVarPos(i,i) = 0;
                }

                if (kConstOri(i,i) + kVarOri(i,i) > KmaxSystem_Ori(i,i) && beta_Ori[i] > 0){
                    kVarOri(i,i) = TauMax[i]/fabs(angular_error[i]) - kConstOri(i,i);
                }else if(kTotalDesOri(i,i) > KmaxSystem_Ori(i,i)){
                    kVarOri(i,i) = 0;
                }
            }

            dTotal6x6(i,i) = dConstLinear(i,i) +  dVarPos(i,i);
            dTotal6x6(i+3,i+3) = dConstAngular(i,i) + dVarOri(i,i);

            // ATTRACTIVE FORCE FOR POSITION AND ORIENTATION CTRL CALCULATION
            Fk[i] = (kConstPos(i,i) + kVarPos(i,i)) * (position_error[i]);
            Fk_Ori[i] = (kConstOri(i,i) + kVarOri(i, i)) * angular_error[i];

            // DAMPING FORCES FOR LINEAR AND ANGULAR VELOCITIES - DEFINITIONS
            Fd_Pos[i] = dConstLinear(i,i) * velLinearError[i];
            Fd_Ori[i] = dConstAngular(i,i) * velAngularError[i];

            //////////////////////// PROPOSED METHOD - END EFFECTOR POSITION ////////////////////////
            calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec(K1_3x3, K2_3x3, KdMax, Pmax3x1, Pmid3x1, maxDistVec_3x1,i);

            //////////////////////// PROPOSED METHOD - END EFFECTOR ORIENTATION ////////////////////////
            calculation_K1_K2_Kd_Max_PMax_PMid_Max_Dist_Vec_Ori(K1_Ori_3x3, K2_Ori_3x3, KdMax_Ori, Pmax3x1_Ori, Pmid3x1_Ori, maxDistVec_Ori_3x1, i);
            ////////////////////////////////////////////////////////////////////////////////////////////

            // TESTING - PROJECTION
            curPos_on_line_projection(p1,p2);

        }

        // NULLSPACE DESIRED JOINT VALUES AND CTRL GAINS INITIALIZATION
        for (int i = 0; i < 7; i++) {
            // Current Nullspace Joint Position and Velocity
            curJntPos7x1[i] = jd->in_robotstatus_var.angles[i];
            curJntVel7x1[i] = jd->in_robotstatus_var.velocities[i];

            // Matrices of Constant Stiffness and Damping gains for Nullspace Ctrl
            kConstNull (i,i) = gainNullP;
            dConstNull (i,i) = gainNullD;

        }

        //Kd_NullSpace = cc->in_jacobian_var.transpose() * kTotal6x6 * cc->in_jacobian_var;

        // NULLSPACE CTRL CALCULATION
        nullSpaceCtrl = nullSpace->calculateNullspace(cc->in_jacobian_var) * (kConstNull * (desJntPos7x1 - curJntPos7x1) + dConstNull * (desJntVel7x1 - curJntVel7x1));

        // END-EFFECTOR POSE ERROR
        for (int i = 0; i < 3; i++) {


            if (sgn(position_error[i], 0.005f) == sgn(velLinearError[i], 0.005f)) {
                posOriError_6x1[i] = position_error[i];
            }else {
                posOriError_6x1[i] = curPos3x1[i] - Pmid3x1[i];
            }

            if(sgn(angular_error[i]) == sgn(velAngularError[i])){
                posOriError_6x1[i+3] = angular_error[i];
            }else{
                posOriError_6x1[i+3] = curOri3x1[i] - Pmid3x1_Ori[i];
            }
        }

        // END-EFFECTOR VEL ERROR
        posOriVelError_6x1.head(3) = velLinearError;
        posOriVelError_6x1.tail(3) = velAngularError;

        //////////////////////// DEFINING THE CONDITIONS FOR THE PROPOSED METHOD - END EFFECTOR POSITION ////////////////////////
        calculation_Kc_DampConstLiner_FMidCtrl_FkCtrl(Kc_3x3, dConstLinear, F_midPointErrorCtrl, Fk, i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //////////////////////// DEFINING THE CONDITIONS FOR THE PROPOSED METHOD - END EFFECTOR ORIENTATION ////////////////////////
        calculation_Kc_DampConstAngular_FMidCtrl_FkCtrl_Ori(Kc_Ori_3x3, dConstAngular, F_midPointErrorCtrl_Ori, Fk_Ori, i);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // REQUIRED END-EFFECTOR CTRL FORCES
        error.segment<3>((i * 6)) = Fk + F_midPointErrorCtrl;
        error.segment<3>((i * 6) + 3) = Fk_Ori + F_midPointErrorCtrl_Ori;
        vel_error.segment<3>((i*6)) = Fd_Pos;
        vel_error.segment<3>((i*6)+3) = Fd_Ori;

    }

    // NATURAL FREQUENCY OBERVER
/*    if (getSimulationTime() >= ctrlStartTime + trajStartTime){
        //RTT::log(RTT::Error) << "KdMax.sqrt() ==> " << KdMax.sqrt()  << RTT::endlog() << RTT::endlog();
        natFreq_observer = ((kTotal6x6 * (lambda_c).inverse()).sqrt())/(2 * M_PI);
        dCritical6x6 = 2 * (kTotal6x6 * lambda_c).sqrt();
    }*/

    // POTENTIAL AND KINETIC ENERGY OF THE SYSTEM CALCULATION
    for(int i = 0; i < 3; i++) {
        if (sgn(position_error[i], 0.005f) == sgn(velLinearError[i], 0.005f)) {
            kTotal6x6(i,i) = kConstPos(i,i) + kVarPos(i,i);
        } else {
            kTotal6x6(i,i) = Kc_3x3(i,i);
        }
        if(sgn(angular_error[i]) == sgn(velAngularError[i])){
            kTotal6x6(i+3,i+3) = kConstOri(i,i) + kVarOri(i,i);
        }else{
            kTotal6x6(i+3,i+3) = Kc_Ori_3x3(i,i);
        }
    }

    potentialEnergy = 0.5 * posOriError_6x1.transpose() * kTotal6x6 * posOriError_6x1;
    kineticEnergy = 0.5 * posOriVelError_6x1.transpose() * lambda_c * posOriVelError_6x1;

    // REQUIRED CTRL JOINT TORQUES CALCULATION AND IMPLEMENTATON
    Eigen::MatrixXf JMcinv = cc->in_jacobian_var * inertia_c_inv;
    Eigen::VectorXf F = ((error + vel_error) + lambda_c * (dc->in_acceleration_var + JMcinv * (dc->in_coriolis_var + dc->in_gravity_var) - dc->in_jacobian_dot_var * jd->in_robotstatus_var.velocities));

    if (getSimulationTime() <= ctrlStartTime){
        cc->out_torques_var.torques = dc->in_gravity_var;
    } else{
        cc->out_torques_var.torques = cc->in_jacobian_var.transpose() * F  + nullSpaceCtrl;
    }

    trqSat->saturateTorqueRate(cc->out_torques_var.torques,estimated_torque_var);
    cc->out_torques_var.torques = estimated_torque_var;

    out_estimated_torque_port.write(estimated_torque_var);


    // END-EFFECTOR POSE AND VEL ERROR IN THE PREVIOUS TIME STEP
    preEndEff_PosError3x1 = position_error;
    preEndEff_VelError3x1 = velLinearError;

    preEndEff_OriError3x1 = angular_error;
    preEndEff_VelOriError3x1 = velAngularError;

    // ON SCREEN DATA PRINTING SECTION

    //<< "FT_Sensor - Torques ==> " << dc->in_forceTorqueSensor_var.torques.transpose() << RTT::endlog() << RTT::endlog();

/*    RTT::log(RTT::Error) << "Pmax3x1 = " << Pmax3x1.transpose() << RTT::endlog() << "Pmid3x1 = " << Pmid3x1.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KmaxSystem ==> " << KmaxSystem.diagonal().transpose() << RTT::endlog() << "beta ==> " << beta.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KdMax = " << (KdMax.diagonal()).transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "K1_3x3 = " << (K1_3x3.diagonal()).transpose() << RTT::endlog() << "K2_3x3 = " << (K2_3x3.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fmax ==> " << Fmax.transpose() << RTT::endlog() << "maxDisp ==> " << maxDisp.transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "dConstLinear_3x3 ==> " << (dConstLinear.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fk ==> " << Fk.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "F_midPointErrorCtrl ==> " << F_midPointErrorCtrl.transpose() << RTT::endlog() << RTT::endlog();

    RTT::log(RTT::Error) << "Pmax3x1_Ori = " << Pmax3x1_Ori.transpose() << RTT::endlog() << "Pmid3x1_Ori = " << Pmid3x1_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KmaxSystem_Ori ==> " << KmaxSystem_Ori.diagonal().transpose() << RTT::endlog() << "beta_Ori ==> " << beta_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "KdMax_Ori = " << (KdMax_Ori.diagonal()).transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "K1_Ori_3x3 = " << (K1_Ori_3x3.diagonal()).transpose() << RTT::endlog() << "K2_Ori_3x3 = " << (K2_Ori_3x3.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "TauMax ==> " << TauMax.transpose() << RTT::endlog() << "maxDisp_Ori ==> " << maxDisp_Ori.transpose() << RTT::endlog();
    RTT::log(RTT::Error) << "dConstAngular_3x3 ==> " << (dConstAngular.diagonal()).transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "Fk_Ori ==> " << Fk_Ori.transpose() << RTT::endlog() << RTT::endlog();
    RTT::log(RTT::Error) << "F_midPointErrorCtrl_Ori ==> " << F_midPointErrorCtrl_Ori.transpose() << RTT::endlog() << RTT::endlog();*/

    //RTT::log(RTT::Error) << "natFreq_observer ==> " << natFreq_observer << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Error) << "K_total ==> " << kTotal6x6 << RTT::endlog() << RTT::endlog() << "lambda_c ==> " << lambda_c << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Error) << "KmaxSystem ==> " << KmaxSystem.diagonal().transpose() << RTT::endlog() << "beta ==> " << beta.transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "maxDisp ==> " << maxDisp.transpose() << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Error) << "angular Error ==> " << angular_error.transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "desOri  ==> " << (desOri3x1).transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "curOri ==> " << (curOri3x1).transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "desired Ori ==> " << desOri3x1.transpose() << RTT::endlog() << RTT::endlog() << "current Ori ==> " << curOri3x1.transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "Fk_Ori ==> " << Fk_Ori.transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "F_midPointErrorCtrl_Ori ==> " << F_midPointErrorCtrl_Ori.transpose() << RTT::endlog() << RTT::endlog();
    //RTT::log(RTT::Error) << "K1_Ori_3x3 = " << (K1_Ori_3x3.diagonal()).transpose() << RTT::endlog() << "K2_Ori_3x3 = " << (K2_Ori_3x3.diagonal()).transpose() << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Error) << "dCritical6x6 ==>" << dCritical6x6 << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Error) << "Position Error ==> " << (desPos3x1 - curPos3x1).transpose() << RTT::endlog() << RTT::endlog() ;

/*    if (getSimulationTime() >= ctrlStartTime + trajStartTime){
        RTT::log(RTT::Error) << "KdMax.sqrt() ==> " << KdMax.sqrt()  << RTT::endlog() << RTT::endlog();
    }*/

    //RTT::log(RTT::Error) << "kConstPos ==> " << kConstPos << RTT::endlog() << RTT::endlog();

    //RTT::log(RTT::Info) << "Fmax ==> " << Fmax.transpose() << RTT::endlog() << "maxDisp ==> " << maxDisp.transpose() << RTT::endlog();
    //RTT::log(RTT::Error) << "FT_Sensor - Forces_Z ==> " << dc->in_forceTorqueSensor_var.forces[2] + 2.055 << RTT::endlog() << RTT::endlog();
    // <<  dc->in_forceTorqueSensor_var.forces[0] << " " << dc->in_forceTorqueSensor_var.forces[1] + 0.055 << " "


    // DATA PRINTING FOR DATA ANALYSIS USING MATLAB
/*    RTT::log(RTT::Error) << getSimulationTime() << "," << position_error[0] << "," << position_error[1] << "," << position_error[2] <<
                                                   "," << desPos3x1[0] << "," << desPos3x1[1]  << "," << desPos3x1[2]               <<
                                                   "," << curPos3x1[0] << "," << curPos3x1[1] << "," << curPos3x1[2]                <<
                                                   "," << angular_error[0] << "," << angular_error[1] << "," << angular_error[2]    <<
                                                   "," << desOri3x1[0] << "," << desOri3x1[1] << "," << desOri3x1[2]                <<
                                                   "," << curOri3x1[0] << "," << curOri3x1[1] << "," << curOri3x1[2]                <<
                                                   "," << kTotal6x6(0,0) << "," << kTotal6x6(1,1) << "," << kTotal6x6(2,2)          <<
                                                   "," << kTotal6x6(3,3) << "," << kTotal6x6(4,4) << "," << kTotal6x6(5,5)          <<
                                                   "," << potentialEnergy + kineticEnergy                                           <<
                                                   "," << Fk[0] << "," << Fk[1] << "," << Fk[2]                                     <<
                                                   "," << F_midPointErrorCtrl[0] << "," << F_midPointErrorCtrl[1] << "," << F_midPointErrorCtrl[2] <<
                                                   "," << F[0] << "," << F[1] << "," << F[2] << "," << F[3] << "," << F[4] << "," << F[5] << "," <<
                                                   "," << potentialEnergy << "," <<  kineticEnergy << RTT::endlog();*/
                                                   //"," <<  stabilityVerfication                                                     << RTT::endlog();

}

bool CartImpedanceControl::configureHook() {
    return cc->configureHook() && dc->configureHook() && jd->configureHook();
}

bool CartImpedanceControl::startHook() {
    time1 = getSimulationTime();
    return cc->startHook() && dc->startHook() && jd->startHook();
}

void CartImpedanceControl::updateHook() {
    time2 = getSimulationTime();

    if (!cc->getData() || !dc->getData() || !jd->getData()) {
        RTT::log(RTT::Error) << this->getName() << " Cart:" << std::to_string(cc->getData()) << "|| DC:"
                             << std::to_string(dc->getData()) << "|| JD:" << std::to_string(jd->getData())
                             << RTT::endlog();
        return;
    }

    // DESIRED END-EFFECTOR POSITION
    desPos3x1 = cc->des_poses_var.head(3);

    // PASSIVE STATE-DEPENDENT VARIABLE IMPEDANCE CONTROL
    compute();

	// WRITING THE CTRL TORQUES TO THEIR PORTS
    cc->out_torques_port.write(cc->out_torques_var);

	// CTRL TIME UPDATE
    time1 = time2;
}

void CartImpedanceControl::stopHook() {
    cc->stopHook();
    dc->stopHook();
    jd->stopHook();
}

void CartImpedanceControl::cleanupHook() {
    cc->cleanupHook();
    dc->cleanupHook();
    jd->cleanupHook();
}
ORO_CREATE_COMPONENT_LIBRARY()
ORO_LIST_COMPONENT_TYPE(CartImpedanceControl)
