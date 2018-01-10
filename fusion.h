#pragma once



struct type_ahrs {
    double roll, pitch, yaw;
    double force_x, force_y, force_z;
};

struct type_rtk  {
    double time,x, y, z;
    double speed, heading;//isNaN judge the value
    int   state_star;
    struct timeval time_stamp;
};

struct type_imu{
    double ax, ay, az;
    double gx, gy, gz;
    double mx, my, mz;
    double time;
    struct timeval time_stamp;
};

struct comp_imu {
    double ax, ay, az;
    double gx, gy, gz;
    double mx, my, mz;
    struct timeval time_stamp;
};

struct type_uwb {
    double time,x, y, z;
    int cred;//���Ŷ�
    struct timeval time_stamp;
};

struct type_indor_cal{
    double time,x, y, z;
    double heading, speed;
    int cred;
    bool is_updated;
    struct timeval time_stamp;
};

struct type_outdor_cal{
    double time,x, y, z;
    double heading, speed;
    int state_star;
    bool is_updated;
    struct timeval time_stamp;
};



class fusion {

public:
    double ax_err[60], ay_err[60], az_err[60];// ��̬�жϣ�
    //once it have it`s value can`t be change
    const double deg2rad;
    const double rad2deg;

    const double gravity;
    int num1 = 0;// ��װ���60
    struct {
        double A[2 * 2];
        double B[2 * 2];
        double B2[2 * 2];
        double I[2 * 2];
        double H[2 * 2];
        double hP[2 * 2];
        double P_prio[2 * 2];
        double A_hP[2 * 2];
        double HPR[2 * 2];
        double hX[2];
        double R[2 * 2];
        double Q[2 * 2];
        double contr_M[2];//���ƾ���
        double hX_prio[2], A_hX_prio[2], B_contr[2];
        double Pos_rtk[2], K[2 * 2], Pos_hX[2], K_Pos_hX[2];
        double IKH[2 * 2];
    } rtkkf;


    struct {
        double A[2*2];
        double B[2*2];
        double B2[2 * 2];
        double I[2*2];
        double H[2*2];
        double hP[2*2];
        double P_prio[2 * 2];
        double A_hP[2*2];
        double HPR[2 * 2];
        double hX[2];
        double R[2*2];
        double Q[2*2];
        double contr_M[2];//���ƾ���
        double hX_prio[2], A_hX_prio[2], B_contr[2];
        double Pos_uwb[2],K[2*2],Pos_hX[2], K_Pos_hX[2];
        double IKH[2 * 2];
    } uwbkf;
    struct {
        double center[3];
        double scater[9];
    } magcomp;

    struct type_ahrs ahrs,attitude;
    struct type_uwb  uwb;
    struct type_rtk  rtk;
    struct type_imu rawimu,calimu,last_imu;//
    struct type_indor_cal indor_cal;
    struct type_outdor_cal outdor_cal;

    double install_acc[2];
    int imu_count=0;
    int state_installerr;

    int cal_installerr( struct type_imu &rawimu);
    void comp_installerr(struct type_imu &rawimu,struct type_imu &calimu);
    void cal_rpy(struct type_imu &calimu, struct type_ahrs &attitude);
          //void ahrs(struct comp_imu &cal_data, double rpy[3]);
    void input_uwb(struct type_uwb &rawuwb, struct type_imu &calimu, struct type_indor_cal &indor);
    void input_imuuwb(struct type_imu &calimu,struct type_uwb &rawuwb,struct type_indor_cal &indor);
    void input_rtk(struct type_rtk &rawrtk, struct type_imu &calimu, struct type_outdor_cal &outdor);
    void input_imurtk(struct type_imu &calimu, struct type_rtk &rawrtk, struct type_outdor_cal &outdor);
    // ֻҪ�������ݴ���0.1 �Ϳ�ʼ��¼imu�����Ƶ�ֵ Ȼ�������Բ��� ��ʱֻ��ȥ��ƫ
    fusion();
};//class fusion
