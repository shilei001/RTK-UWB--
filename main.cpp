/* test main.cpp
2018.1.5
Author shilei
*/


#include <iostream>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include <fusion/fusion.h>

int main()
{

    FILE *fraw;
    FILE *fpostindor;
    FILE *fpostoutdor;

    fusion KF;// define  the real


    double install_acc[2]={0.0};

    char fppostin[1024];
    char fppostout[1024];
    char line[2048];
    char datapath[1024] = "/home/shilei/Desktop/rtkuwbimu1228.txt";
    char resultpath[1024] = "/home/shilei/Desktop/";
    fraw = fopen(datapath, "r");

    if (fraw == NULL)
    {
        printf("open 一体化终端数据 file error\n");
        return 0;
    }

    strcpy(fppostin, resultpath);
    strcat(fppostin, "indorpostdata.txt");
    fpostindor = fopen(fppostin, "w");
    if (fpostindor == NULL)
    {
        printf("write indorpostdate file error\n");
    }

    strcpy(fppostout, resultpath);
    strcat(fppostout, "outdorpostdata.txt");
    fpostoutdor = fopen(fppostout, "w");
    if (fpostoutdor == NULL)
    {
        printf("write outdorpostdate file error\n");
    }

    struct type_imu rawimu, calimu;
    struct type_ahrs ahrs;
    struct type_uwb  uwb;
    struct type_rtk  rtk;
    struct type_indor_cal  indor_cal;
    struct type_outdor_cal  outdor_cal;
    int indor_outdor=1;// indor as default


    while (!feof(fraw))//逐行开始读数据
    {
        auto i = fgets(line, sizeof(line), fraw);
        if(i==NULL){
            continue;
        }

        if (line[0] == 'i')
        {
            sscanf(line,"imu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
                   &rawimu.time,&rawimu.ax,&rawimu.ay,&rawimu.az,&rawimu.gx,&rawimu.gy,&rawimu.gz,&rawimu.mx,&rawimu.my,&rawimu.mz);

            if(KF.state_installerr==0)//we guess car in the level road.
             {
               sscanf(line,"imu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\r\n",
                   &rawimu.time,
                   &rawimu.ax,&rawimu.ay,&rawimu.az,
                   &rawimu.gx,&rawimu.gy,&rawimu.gz,
                   &rawimu.mx,&rawimu.my,&rawimu.mz);

                KF.cal_installerr(rawimu);// it will run 60times  if not wrong
            }

            else if ( indor_outdor==1 )// indoor
            {
                KF.comp_installerr(rawimu,calimu);//once it have imu this step must have
                KF.cal_rpy(calimu, ahrs );//innitialize the roll and pitch yaw calculate
                KF.input_imuuwb(calimu,uwb,indor_cal);

                fprintf(fpostindor,"indoor %12.3lf %12.3lf %4.2lf %4.2lf \n",indor_cal.x,indor_cal.y,indor_cal.z,indor_cal.cred);
            }

            else if(indor_outdor==2)// outdoor
            {
                KF.comp_installerr(rawimu,calimu);
                KF.cal_rpy(calimu, ahrs );//innitialize the roll and pitch yaw calculate
                KF.input_imurtk(calimu,rtk,outdor_cal) ;

                fprintf(fpostoutdor,"outdor %12.3lf %12.3lf %4.2lf  %d \n",outdor_cal.x,outdor_cal.y,outdor_cal.z,outdor_cal.state_star);
             }
        }

        else if (line[0] == 'u')
        {
           // printf("uwb\n");
            sscanf(line,"uwb %lf %lf %lf %lf\r\n",
                    &uwb.time,&uwb.x,&uwb.y,&uwb.z);

            KF.input_uwb(uwb,calimu,indor_cal);
            indor_outdor=1;// indoor=1 outdoor=0
             fprintf(fpostindor,"indoor %12.3lf %12.3lf %4.2lf %2.2lf \n",indor_cal.x,indor_cal.y,indor_cal.z,indor_cal.cred);
        }
        else if (line[0] == 'r')
        {
           // printf("rtk\n");
            sscanf(line,"rtk %lf %lf %lf %lf %d\r\n",
                   &rtk.time,&rtk.x,&rtk.y,&rtk.z,&rtk.state_star);

            KF.input_rtk(rtk,calimu,outdor_cal);
            indor_outdor=2;
            fprintf(fpostoutdor,"outdor %12.3lf %12.3lf %4.2lf  %d \n",outdor_cal.x,outdor_cal.y,outdor_cal.z,outdor_cal.state_star);
        }
        //fgetc(stdin);
        //system("pause");
    }

}
