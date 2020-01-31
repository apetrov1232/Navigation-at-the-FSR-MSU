#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;


int main() {
    ifstream fl,fl0;
    long long int i=0,l; long double x,y,z,a,b,c,h=0.001,e,qp1,qp2,qp3,qp4;
    fl0.open("C:\\Output_1881695243.txt");
    while (fl0>>x) i++;
    fl0.close();
    l=i/9;
    cout.precision(10);
    //qp1=0.26222;
    //qp2=-0.64984;
    //qp3=0.57530;
    //qp4=0.42188;
    fl.open("C:\\Output_1881695243.txt");
    ofstream res("C:\\quat2.txt");
    vector <long double> aksx(l),aksy(l),aksz(l),naksx(l),naksy(l),naksz(l),
    magx(l),magy(l),magz(l), nmagx(l),nmagy(l),nmagz(l),girx(l),giry(l),
    girz(l),ogirx(l),ogiry(l),ogirz(l);
    long double WW,XX,YY,ZZ,qw,qx,qy,qz,Ww,Wx,Wy,Wz,Dw,Dx,Dy,Dz,
    Norm,oXX,oYY,oZZ,B=0.165;

    vector <vector <long double> > matr(l, vector<long double> (16,0));
    for (i=0; i<l; i++){
        fl>>x>>y>>z;
        aksx[i]=x/sqrt(x*x+y*y+z*z);
        aksy[i]=y/sqrt(x*x+y*y+z*z);
        aksz[i]=z/sqrt(x*x+y*y+z*z);
        fl>>x;
        ogirx[i]=-x;
        fl>>x;
        ogiry[i]=-x;
        fl>>x;
        ogirz[i]=-x;
        fl>>x>>y>>z;
        magx[i]=(x-(-0.023323))*5.349270+(y-(-0.438206))*0.178976+(z-(0.066518))*0.066518; //îòêàëèáðîâàííûå äàííûå=ìàòðèöà êàëèáðîâêè*(äàííûå-{âåêòîð ñìåùåíèÿ})
        magz[i]=-((x-(-0.023323))*0.178976+(y-(-0.438206))*5.910940+(z-(0.066518))*(-0.012636));
        magy[i]=-((x-(-0.023323))*(-0.236383)+(y-(-0.438206))*(-0.012636)+(z-(0.066518))*5.583866);
        Norm=sqrt(magx[i]*magx[i]+magy[i]*magy[i]+magz[i]*magz[i]);
        magx[i]/=Norm;  magy[i]/=Norm;  magz[i]/=Norm;
        a=atan(-aksx[0]/aksy[0]); //óãëû ïîâîðîòîâ â ìàòðèöå ïåðåõîäà ê íîâîé ÑÊ (ïîñëåäîâàòåëüíî âîêðóã Z,X,Z), ÷òîáû g0={0,0,g}, B0={B1,0,B2}
        b=atan((aksx[0]*sin(a)-aksy[0]*cos(a))/aksz[0]);
        c=atan((cos(a)*cos(b)*magy[0]-cos(b)*magx[0]*sin(a)+magz[0]*sin(b))/(cos(a)*magx[0]+magy[0]*sin(a)));
        nmagx[i]=magz[i]*sin(b)*sin(c)+magy[i]*(cos(c)*sin(a)+cos(a)*cos(b)*sin(c))+magx[i]*(cos(a)*cos(c)-cos(b)*sin(a)*sin(c));
        nmagy[i]=cos(c)*magz[i]*sin(b)+magx[i]*(-cos(b)*cos(c)*sin(a)-cos(a)*sin(c))+magy[i]*(cos(a)*cos(b)*cos(c)-sin(a)*sin(c));
        nmagz[i]=cos(b)*magz[i]-cos(a)*magy[i]*sin(b)+magx[i]*sin(a)*sin(b);
        naksx[i]=aksz[i]*sin(b)*sin(c)+aksy[i]*(cos(c)*sin(a)+cos(a)*cos(b)*sin(c))+aksx[i]*(cos(a)*cos(c)-cos(b)*sin(a)*sin(c));
        naksy[i]=cos(c)*aksz[i]*sin(b)+aksx[i]*(-cos(b)*cos(c)*sin(a)-cos(a)*sin(c))+aksy[i]*(cos(a)*cos(b)*cos(c)-sin(a)*sin(c));
        naksz[i]=cos(b)*aksz[i]-cos(a)*aksy[i]*sin(b)+aksx[i]*sin(a)*sin(b);
        girx[i]=ogirz[i]*sin(b)*sin(c)+ogiry[i]*(cos(c)*sin(a)+cos(a)*cos(b)*sin(c))+ogirx[i]*(cos(a)*cos(c)-cos(b)*sin(a)*sin(c));
        giry[i]=cos(c)*ogirz[i]*sin(b)+ogirx[i]*(-cos(b)*cos(c)*sin(a)-cos(a)*sin(c))+ogiry[i]*(cos(a)*cos(b)*cos(c)-sin(a)*sin(c));
        girz[i]=cos(b)*ogirz[i]-cos(a)*ogiry[i]*sin(b)+ogirx[i]*sin(a)*sin(b);
    };
    WW=1; XX=0; YY=0; ZZ=0;
    oXX=XX; oYY=YY; oZZ=ZZ;
    for (i=0; i<l; i++){
        res<<WW<<'\t'<<-XX<<'\t'<<-YY<<'\t'<<-ZZ<<endl;
        Ww=0.5*(-XX*girx[i]-YY*giry[i]-ZZ*girz[i]);
        Wx=0.5*(WW*girx[i]+YY*girz[i]-ZZ*giry[i]);
        Wy=0.5*(WW*giry[i]-XX*girz[i]+ZZ*girx[i]);
        Wz=0.5*(WW*girz[i]+XX*giry[i]-YY*girx[i]);

        Dw=-2*YY*(2*(-WW*YY +XX*ZZ) - naksx[i]) + 2*XX*(2*(WW*XX + YY*ZZ) //D-ãðàäèåíò
        - naksy[i]) + (-2*ZZ*nmagx[0] +
        2*XX*nmagz[0])*(2*(XX*YY - WW*ZZ)*nmagx[0] - nmagy[i] +
        2*(WW*XX + YY*ZZ)*nmagz[0]) +
        2*YY*nmagx[0]*(2*(WW*YY + XX*ZZ)*nmagx[0] - nmagz[i] +
        2*nmagz[0]*(0.5 - XX *XX - YY *YY)) -
        2*YY*nmagz[0]*(-nmagx[i] + 2*(-WW*YY + XX*ZZ)*nmagz[0] +
        2*nmagx[0]*(0.5 - YY *YY - ZZ *ZZ));

        Dx=2*ZZ *(2*(-WW*YY + XX*ZZ) - naksx[i]) +
        2*WW*(2*(WW*XX + YY*ZZ) - naksy[i]) + (2*YY*nmagx[0] +
        2*WW*nmagz[0])*(2*(XX*YY - WW*ZZ)*nmagx[0] - nmagy[i] +
        2*(WW*XX + YY*ZZ)*nmagz[0]) -4*XX*(-naksz[i]
        + 2*(0.5 - XX *XX - YY *YY)) + (2*ZZ*nmagx[0] -
        4*XX*nmagz[0])*(2*(WW*YY + XX*ZZ)*nmagx[0] - nmagz[i] +
        2*nmagz[0]*(0.5 - XX *XX - YY *YY)) +
        2*ZZ*nmagz[0]*(-nmagx[i] + 2*(-WW*YY + XX*ZZ)*nmagz[0] +
        2*nmagx[0]*(0.5 - YY *YY - ZZ *ZZ));

        Dy=-2*WW *(2*(-WW*YY + XX*ZZ) - naksx[i]) +
        2*ZZ*(2*(WW*XX + YY*ZZ) - naksy[i]) + (2*XX*nmagx[0] +
        2*ZZ*nmagz[0])*(2*(XX*YY - WW*ZZ)*nmagx[0] - nmagy[i] +
        2*(WW*XX + YY*ZZ)*nmagz[0]) -
        4*YY*(-naksz[i] + 2*(0.5 - XX *XX - YY *YY)) + (2*WW*nmagx[0] -
        4*YY*nmagz[0])*(2*(WW*YY + XX*ZZ)*nmagx[0] - nmagz[i] +
        2*nmagz[0]*(0.5 - XX *XX - YY *YY)) + (-4*YY*nmagx[0] -
        2*WW*nmagz[0])*(-nmagx[i] + 2*(-WW*YY + XX*ZZ)*nmagz[0] +
        2*nmagx[0]*(0.5 - YY *YY - ZZ *ZZ));

        Dz=2*XX*(2*(-WW*YY + XX*ZZ) - naksx[i]) +
        2*YY*(2*(WW*XX + YY*ZZ) - naksy[i]) + (-2*WW*nmagx[0] +
        2*YY*nmagz[0])*(2*(XX*YY - WW*ZZ)*nmagx[0] - nmagy[i] +
        2*(WW*XX + YY*ZZ)*nmagz[0]) +
        2*XX*nmagx[0]*(2*(WW*YY + XX*ZZ)*nmagx[0] - nmagz[i] +
        2*nmagz[0]*(0.5 - XX *XX - YY *YY)) + (-4*ZZ*nmagx[0] +
        2*XX*nmagz[0])*(-nmagx[i] + 2*(-WW*YY + XX*ZZ)*nmagz[0] +
        2*nmagx[0]*(0.5 - YY *YY - ZZ *ZZ));

        Norm=sqrt(Dw*Dw+Dx*Dx+Dy*Dy+Dz*Dz);
        Dx/=Norm; Dz/=Norm; Dy/=Norm; Dw/=Norm;
        qw=Ww-B*Dw; qx=Wx-B*Dx; qy=Wy-B*Dy; qz=Wz-B*Dz;
        WW+=h*qw; XX+=h*qx; YY+=h*qy; ZZ+=h*qz;
        Norm=sqrt(WW*WW+YY*YY+XX*XX+ZZ*ZZ);
        WW/=Norm; XX/=Norm; YY/=Norm; ZZ/=Norm;
        oYY=YY;
        oZZ=ZZ;
    };
    res.close();
    fl.close();
    return 0;
}
