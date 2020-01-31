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



vector < vector <long double> > inverse(vector < vector <long double> > a){ //ôóíêöèÿ äëÿ íàõîæäåíèÿ îáðàòíîé ìàòðèöû
    int  n = a.size();
    vector < vector <long double> > ans(n, vector <long double> (n, 0));
    for (int i = 0; i < n; i++){
        ans[i][i] = 1.0;
    }
    for (int i = 0; i < n; i++){
        int row = i;
        long double mx = a[i][i];
        for(int k = i + 1; k < n; k++){
            if (abs(a[k][i]) > mx){
                row = k;
                mx = abs(a[k][i]);
            }
        }
        if (row != i){
            swap(a[row], a[i]);
            swap(ans[row], ans[i]);
        }
        for (int j = i+1; j < n; j++){
            long double e = a[j][i]/a[i][i];
            for (int k = 0; k < n; k++){
                a[j][k] -= e*a[i][k];
                ans[j][k] -= e*ans[i][k];
            }
        }
    }
    for (int i = n - 1; i >= 0; i--){
        for (int j = i - 1; j >= 0; j--){
            long double e = a[j][i]/a[i][i];
            for (int k = 0; k < n; k++){
                a[j][k] -= e*a[i][k];
                ans[j][k] -= e*ans[i][k];
            }
        }
        for (int j = 0; j < n; j++) {
            ans[i][j] /= a[i][i];
        }
    }
    return ans;
}




int main() {
    ifstream fl,fl0;
    long long int i=0,l; long double x,y,z,a,b,c,h=0.001,e,qp1,qp2,qp3,qp4;
    fl0.open("C:\\Output_1881695243.txt"); vector <vector <long double> > mat(3);
    mat[0].resize(3); mat[1].resize(3); mat[2].resize(3);
    while (fl0>>x) i++;
    fl0.close();
    l=i/9;
    cout.precision(10);
    //qp1=0.26222;
    //qp2=-0.64984;
    //qp3=0.57530;
    //qp4=0.42188;
    fl.open("C:\\Output_1881695243.txt");
    ofstream res("C:\\quat1.txt");
    vector <long double> aksx(l),aksy(l),aksz(l),naksx(l),naksy(l),naksz(l),
    magx(l),magy(l),magz(l), nmagx(l),nmagy(l),nmagz(l),xm1(l), xm2(l),
    xm3(l), xp1(l), xp2(l), xp3(l), q11(l), q12(l), q13(l), q21(l), q22(l),
    q23(l), q31(l), q32(l), q33(l), pm11(l), pm12(l), pm13(l), pm21(l),
    pm22(l), pm23(l), pm31(l), pm32(l), pm33(l), p11(l), p12(l), p13(l),
    p21(l), p22(l), p23(l), p31(l), p32(l), p33(l), k11(l), k12(l), k13(l),
    k21(l), k22(l), k23(l), k31(l), k32(l), k33(l), wx(l), wy(l), wz(l),
    gamma(l), beta(l),girx(l),giry(l),girz(l),ogirx(l),ogiry(l),ogirz(l),
    alpha(l),alpham(l),X(l),Y(l),Z(l);
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
    xp1[0]=0.00000000000; xp2[0]=0.00000000000; xp3[0]=1; xm1[0]=xp1[0]; xm2[0]=xp2[0]; xm3[0]=xp3[0];
    p11[0]=0.00179455; p12[0]=0; p13[0]=0; p21[0]=0; p22[0]=0.00179455;
    p23[0]=0; p31[0]=0; p32[0]=0; p33[0]=0.00179455;
    k11[0]=1; k12[0]=0; k13[0]=0; k21[0]=0; k22[0]=1;
    k23[0]=0; k31[0]=0; k32[0]=0; k33[0]=1;
    for (i=0; i<l-1; i++){
        wx[i]=0;//-h*(-0.0679264*xp2[i]+0.125511*xp3[i]); //âåêòîð øóìà ãèðîñêîïà ={0.108275463, 0.125510985, 0.067926358}
        wy[i]=0;//-h*(0.0679264*xp1[i]-0.108275*xp3[i]);  //âåêòîð w-ïîãðåøíîñòü îò ãèðîñêîïà
        wz[i]=0;//-h*(-0.125511*xp1[i]+0.108275*xp2[i]);
        xm1[i+1]=xp1[i]-h*girz[i]*xp2[i]+h*giry[i]*xp3[i];//xm-âåêòîð, êîòîðûé äîëæåí ïîëó÷èòüñÿ, èñõîäÿ èç äàííûõ ãèðîñêîïà
        xm2[i+1]=h*girz[i]*xp1[i]+xp2[i]-h*girx[i]*xp3[i];
        xm3[i+1]=-h*giry[i]*xp1[i]+h*girx[i]*xp2[i]+xp3[i];
        e=sqrt(xm1[i]*xm1[i]+xm2[i]*xm2[i]+xm3[i]*xm3[i]);
        xm1[i]=xm1[i]/e; xm2[i]=xm2[i]/e; xm3[i]=xm3[i]/e;
        q11[i]=0.0320906*h*h*(xp2[i]*xp2[i]+xp3[i]*xp3[i]); //q-ìàòðèöà êîâàðèàöèé ïðîöåññà
        q12[i]=-0.0320906*h*h*xp1[i]*xp2[i];
        q13[i]=-0.0320906*h*h*xp1[i]*xp3[i];
        q21[i]=-0.0320906*h*h*xp1[i]*xp2[i];
        q22[i]=0.0320906*h*h*(xp1[i]*xp1[i]+xp3[i]*xp3[i]);
        q23[i]=-0.0320906*h*h*xp2[i]*xp3[i];
        q31[i]=-0.0320906*h*h*xp1[i]*xp3[i];
        q32[i]=-0.0320906*h*h*xp2[i]*xp3[i];
        q33[i]=0.0320906*h*h*(xp1[i]*xp1[i]+xp2[i]*xp2[i]);

        pm11[i+1]=p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i] - h*girz[i]*(p12[i] - h*girz[i]*p22[i] + h*giry[i]*p32[i]) + h*giry[i]*(p13[i] - h*girz[i]*p23[i] + h*giry[i]*p33[i]) + q11[i]; //ìàòðèöà pm-ïðîìåæóòî÷íàÿ äëÿ âû÷èñëåíèÿ êîýôôèöèåíòà êàëìàíà
        pm12[i+1]=p12[i] - h*girz[i]*p22[i] + h*girz[i]*(p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i]) + h*giry[i]*p32[i] - h*girx[i]*(p13[i] - h*girz[i]*p23[i] + h*giry[i]*p33[i]) + q12[i];
        pm13[i+1]=p13[i] - h*girz[i]*p23[i] - h*giry[i]*(p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i]) + h*girx[i]*(p12[i] - h*girz[i]*p22[i] + h*giry[i]*p32[i]) + h*giry[i]*p33[i] + q13[i];
        pm21[i+1]=h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i] - h*girz[i]*(h*girz[i]*p12[i] + p22[i] - h*girx[i]*p32[i]) + h*giry[i]*(h*girz[i]*p13[i] + p23[i] - h*girx[i]*p33[i]) + q21[i];
        pm22[i+1]=h*girz[i]*p12[i] + p22[i] + h*girz[i]*(h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i]) - h*girx[i]*p32[i] - h*girx[i]*(h*girz[i]*p13[i] + p23[i] - h*girx[i]*p33[i]) + q22[i];
        pm23[i+1]=h*girz[i]*p13[i] + p23[i] - h*giry[i]*(h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i]) + h*girx[i]*(h*girz[i]*p12[i] + p22[i] - h*girx[i]*p32[i]) - h*girx[i]*p33[i] + q23[i];
        pm31[i+1]=-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i] - h*girz[i]*(-h*giry[i]*p12[i] + h*girx[i]*p22[i] + p32[i]) + h*giry[i]*(-h*giry[i]*p13[i] + h*girx[i]*p23[i] + p33[i]) + q31[i];
        pm32[i+1]=-h*giry[i]*p12[i] + h*girx[i]*p22[i] + h*girz[i]*(-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i]) + p32[i] - h*girx[i]*(-h*giry[i]*p13[i] + h*girx[i]*p23[i] + p33[i]) + q32[i];
        pm33[i+1]=-h*giry[i]*p13[i] + h*girx[i]*p23[i] - h*giry[i]*(-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i]) + h*girx[i]*(-h*giry[i]*p12[i] + h*girx[i]*p22[i] + p32[i]) + p33[i] + q33[i];

        mat[0][0]=pm11[i+1]+0.00179455; //0.00179455-íîðìà â êâàäðàòå âåêòîðà øóìà àêñåëåðîìåòðà
        mat[0][1]=pm12[i+1];
        mat[0][2]=pm13[i+1];
        mat[1][0]=pm21[i+1];
        mat[1][1]=pm22[i+1]+0.00179455;
        mat[1][2]=pm23[i+1];
        mat[2][0]=pm31[i+1];
        mat[2][1]=pm32[i+1];
        mat[2][2]=pm33[i+1]+0.00179455;
        mat=inverse(mat);
        k11[i+1]=pm11[i+1]*mat[0][0]+pm12[i+1]*mat[1][0]+pm13[i+1]*mat[2][0]; //ìàòðèöà êîýôôèöèåíòîâ êàëìàíà
        k12[i+1]=pm11[i+1]*mat[0][1]+pm12[i+1]*mat[1][1]+pm13[i+1]*mat[2][1];
        k13[i+1]=pm11[i+1]*mat[0][2]+pm12[i+1]*mat[1][2]+pm13[i+1]*mat[2][2];
        k21[i+1]=pm21[i+1]*mat[0][0]+pm22[i+1]*mat[1][0]+pm23[i+1]*mat[2][0];
        k22[i+1]=pm21[i+1]*mat[0][1]+pm22[i+1]*mat[1][1]+pm23[i+1]*mat[2][1];
        k23[i+1]=pm21[i+1]*mat[0][2]+pm22[i+1]*mat[1][2]+pm23[i+1]*mat[2][2];
        k31[i+1]=pm31[i+1]*mat[0][0]+pm32[i+1]*mat[1][0]+pm33[i+1]*mat[2][0];
        k32[i+1]=pm31[i+1]*mat[0][1]+pm32[i+1]*mat[1][1]+pm33[i+1]*mat[2][1];
        k33[i+1]=pm31[i+1]*mat[0][2]+pm32[i+1]*mat[1][2]+pm33[i+1]*mat[2][2];

        xp1[i+1]=k11[1 + i]*(naksx[1 + i] - xm1[1 + i]) + xm1[1 + i] + k12[1 + i]*(naksy[1 + i] - xm2[1 + i]) + k13[1 + i]*(naksz[1 + i] - xm3[1 + i]); //xp-ïðîôèëüòðîâàííûé âåêòîð
        xp2[i+1]=k21[1 + i]*(naksx[1 + i] - xm1[1 + i]) + k22[1 + i]*(naksy[1 + i] - xm2[1 + i]) + xm2[1 + i] + k23[1 + i]*(naksz[1 + i] - xm3[1 + i]);
        xp3[i+1]=k31[1 + i]*(naksx[1 + i] - xm1[1 + i]) + k32[1 + i]*(naksy[1 + i] - xm2[1 + i]) + k33[1 + i]*(naksz[1 + i] - xm3[1 + i]) + xm3[1 + i];
        e=sqrt(xp1[i+1]*xp1[i+1]+xp2[i+1]*xp2[i+1]+xp3[i+1]*xp3[i+1]);
        xp1[i+1]=xp1[i+1]/e;
        xp2[i+1]=xp2[i+1]/e;
        xp3[i+1]=xp3[i+1]/e;
        gamma[i] = atan2(xp2[i],xp3[i]); //êðåí è òàíãàæ
        beta[i] = atan2(-xp1[i],(xp2[i]/sin(gamma[i])));

        p11[i+1]=(1 - k11[1 + i])*pm11[1 + i] - k12[1 + i]*pm21[1 + i] - k13[1 + i]*pm31[1 + i];
        p12[i+1]=(1 - k11[1 + i])*pm12[1 + i] - k12[1 + i]*pm22[1 + i] - k13[1 + i]*pm32[1 + i];
        p13[i+1]=(1 - k11[1 + i])*pm13[1 + i] - k12[1 + i]*pm23[1 + i] - k13[1 + i]*pm33[1 + i];
        p21[i+1]=-k21[1 + i]*pm11[1 + i] + (1 - k22[1 + i])*pm21[1 + i] - k23[1 + i]*pm31[1 + i];
        p22[i+1]=-k21[1 + i]*pm12[1 + i] + (1 - k22[1 + i])*pm22[1 + i] - k23[1 + i]*pm32[1 + i];
        p23[i+1]=-k21[1 + i]*pm13[1 + i] + (1 - k22[1 + i])*pm23[1 + i] - k23[1 + i]*pm33[1 + i];
        p31[i+1]=-k31[1 + i]*pm11[1 + i] - k32[1 + i]*pm21[1 + i] + (1 - k33[1 + i])*pm31[1 + i];
        p32[i+1]=-k31[1 + i]*pm12[1 + i] - k32[1 + i]*pm22[1 + i] + (1 - k33[1 + i])*pm32[1 + i];
        p33[i+1]=-k31[1 + i]*pm13[1 + i] - k32[1 + i]*pm23[1 + i] + (1 - k33[1 + i])*pm33[1 + i];
    };
    gamma[l-1] = atan2(xp2[l-1],xp3[l-1]);
    beta[l-1] = atan2(-xp1[l-1],xp2[l-1]/sin(gamma[l-1]));
    for (i=0; i<l; i++){
        alpham[i]=atan2((sin(gamma[i])*nmagz[i]-cos(gamma[i])*nmagy[i]),(cos(beta[i])*nmagx[i]+nmagz[i]*cos(gamma[i])*sin(beta[i])+nmagy[i]*sin(beta[i])*sin(gamma[i])));
        X[i]=cos(alpham[i])*cos(beta[i]);
        Y[i]=cos(alpham[i])*sin(beta[i])*sin(gamma[i])-sin(alpham[i])*cos(gamma[i]);
        Z[i]=cos(alpham[i])*sin(beta[i])*cos(gamma[i])+sin(alpham[i])*sin(gamma[i]);
    };
    xp1[0]=1; xp2[0]=0.0000000000000; xp3[0]=0.000000000000;
    xm1[0]=xp1[0]; xm2[0]=xp2[0]; xm3[0]=xp3[0];
    p11[0]=0.000107364; p12[0]=0; p13[0]=0; p21[0]=0; p22[0]=0.000107364;
    p23[0]=0; p31[0]=0; p32[0]=0; p33[0]=0.000107364;
    k11[0]=1; k12[0]=0; k13[0]=0; k21[0]=0; k22[0]=1;
    k23[0]=0; k31[0]=0; k32[0]=0; k33[0]=1;
    for (i=0; i<l-1; i++){
        xm1[i+1]=xp1[i]-h*girz[i]*xp2[i]+h*giry[i]*xp3[i];//xm-âåêòîð, êîòîðûé äîëæåí ïîëó÷èòüñÿ, èñõîäÿ èç äàííûõ ãèðîñêîïà
        xm2[i+1]=h*girz[i]*xp1[i]+xp2[i]-h*girx[i]*xp3[i];
        xm3[i+1]=-h*giry[i]*xp1[i]+h*girx[i]*xp2[i]+xp3[i];
        e=sqrt(xm1[i]*xm1[i]+xm2[i]*xm2[i]+xm3[i]*xm3[i]);
        xm1[i]=xm1[i]/e; xm2[i]=xm2[i]/e; xm3[i]=xm3[i]/e;
        q11[i]=0.0320906*h*h*(xp2[i]*xp2[i]+xp3[i]*xp3[i]); //q-ìàòðèöà êîâàðèàöèé ïðîöåññà
        q12[i]=-0.0320906*h*h*xp1[i]*xp2[i];
        q13[i]=-0.0320906*h*h*xp1[i]*xp3[i];
        q21[i]=-0.0320906*h*h*xp1[i]*xp2[i];
        q22[i]=0.0320906*h*h*(xp1[i]*xp1[i]+xp3[i]*xp3[i]);
        q23[i]=-0.0320906*h*h*xp2[i]*xp3[i];
        q31[i]=-0.0320906*h*h*xp1[i]*xp3[i];
        q32[i]=-0.0320906*h*h*xp2[i]*xp3[i];
        q33[i]=0.0320906*h*h*(xp1[i]*xp1[i]+xp2[i]*xp2[i]);

        pm11[i+1]=p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i] - h*girz[i]*(p12[i] - h*girz[i]*p22[i] + h*giry[i]*p32[i]) + h*giry[i]*(p13[i] - h*girz[i]*p23[i] + h*giry[i]*p33[i]) + q11[i]; //ìàòðèöà pm-ïðîìåæóòî÷íàÿ äëÿ âû÷èñëåíèÿ êîýôôèöèåíòà êàëìàíà
        pm12[i+1]=p12[i] - h*girz[i]*p22[i] + h*girz[i]*(p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i]) + h*giry[i]*p32[i] - h*girx[i]*(p13[i] - h*girz[i]*p23[i] + h*giry[i]*p33[i]) + q12[i];
        pm13[i+1]=p13[i] - h*girz[i]*p23[i] - h*giry[i]*(p11[i] - h*girz[i]*p21[i] + h*giry[i]*p31[i]) + h*girx[i]*(p12[i] - h*girz[i]*p22[i] + h*giry[i]*p32[i]) + h*giry[i]*p33[i] + q13[i];
        pm21[i+1]=h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i] - h*girz[i]*(h*girz[i]*p12[i] + p22[i] - h*girx[i]*p32[i]) + h*giry[i]*(h*girz[i]*p13[i] + p23[i] - h*girx[i]*p33[i]) + q21[i];
        pm22[i+1]=h*girz[i]*p12[i] + p22[i] + h*girz[i]*(h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i]) - h*girx[i]*p32[i] - h*girx[i]*(h*girz[i]*p13[i] + p23[i] - h*girx[i]*p33[i]) + q22[i];
        pm23[i+1]=h*girz[i]*p13[i] + p23[i] - h*giry[i]*(h*girz[i]*p11[i] + p21[i] - h*girx[i]*p31[i]) + h*girx[i]*(h*girz[i]*p12[i] + p22[i] - h*girx[i]*p32[i]) - h*girx[i]*p33[i] + q23[i];
        pm31[i+1]=-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i] - h*girz[i]*(-h*giry[i]*p12[i] + h*girx[i]*p22[i] + p32[i]) + h*giry[i]*(-h*giry[i]*p13[i] + h*girx[i]*p23[i] + p33[i]) + q31[i];
        pm32[i+1]=-h*giry[i]*p12[i] + h*girx[i]*p22[i] + h*girz[i]*(-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i]) + p32[i] - h*girx[i]*(-h*giry[i]*p13[i] + h*girx[i]*p23[i] + p33[i]) + q32[i];
        pm33[i+1]=-h*giry[i]*p13[i] + h*girx[i]*p23[i] - h*giry[i]*(-h*giry[i]*p11[i] + h*girx[i]*p21[i] + p31[i]) + h*girx[i]*(-h*giry[i]*p12[i] + h*girx[i]*p22[i] + p32[i]) + p33[i] + q33[i];

        mat[0][0]=pm11[i+1]+0.000107364; //0.00179455-íîðìà â êâàäðàòå âåêòîðà øóìà àêñåëåðîìåòðà
        mat[0][1]=pm12[i+1];
        mat[0][2]=pm13[i+1];
        mat[1][0]=pm21[i+1];
        mat[1][1]=pm22[i+1]+0.000107364;
        mat[1][2]=pm23[i+1];
        mat[2][0]=pm31[i+1];
        mat[2][1]=pm32[i+1];
        mat[2][2]=pm33[i+1]+0.000107364;
        mat=inverse(mat);
        k11[i+1]=pm11[i+1]*mat[0][0]+pm12[i+1]*mat[1][0]+pm13[i+1]*mat[2][0]; //ìàòðèöà êîýôôèöèåíòîâ êàëìàíà
        k12[i+1]=pm11[i+1]*mat[0][1]+pm12[i+1]*mat[1][1]+pm13[i+1]*mat[2][1];
        k13[i+1]=pm11[i+1]*mat[0][2]+pm12[i+1]*mat[1][2]+pm13[i+1]*mat[2][2];
        k21[i+1]=pm21[i+1]*mat[0][0]+pm22[i+1]*mat[1][0]+pm23[i+1]*mat[2][0];
        k22[i+1]=pm21[i+1]*mat[0][1]+pm22[i+1]*mat[1][1]+pm23[i+1]*mat[2][1];
        k23[i+1]=pm21[i+1]*mat[0][2]+pm22[i+1]*mat[1][2]+pm23[i+1]*mat[2][2];
        k31[i+1]=pm31[i+1]*mat[0][0]+pm32[i+1]*mat[1][0]+pm33[i+1]*mat[2][0];
        k32[i+1]=pm31[i+1]*mat[0][1]+pm32[i+1]*mat[1][1]+pm33[i+1]*mat[2][1];
        k33[i+1]=pm31[i+1]*mat[0][2]+pm32[i+1]*mat[1][2]+pm33[i+1]*mat[2][2];

        xp1[i+1]=k11[1 + i]*(X[1 + i] - xm1[1 + i]) + xm1[1 + i] + k12[1 + i]*(Y[1 + i] - xm2[1 + i]) + k13[1 + i]*(Z[1 + i] - xm3[1 + i]); //xp-ïðîôèëüòðîâàííûé âåêòîð
        xp2[i+1]=k21[1 + i]*(X[1 + i] - xm1[1 + i]) + k22[1 + i]*(Y[1 + i] - xm2[1 + i]) + xm2[1 + i] + k23[1 + i]*(Z[1 + i] - xm3[1 + i]);
        xp3[i+1]=k31[1 + i]*(X[1 + i] - xm1[1 + i]) + k32[1 + i]*(Y[1 + i] - xm2[1 + i]) + k33[1 + i]*(Z[1 + i] - xm3[1 + i]) + xm3[1 + i];
        e=sqrt(xp1[i+1]*xp1[i+1]+xp2[i+1]*xp2[i+1]+xp3[i+1]*xp3[i+1]);
        xp1[i+1]=xp1[i+1]/e;
        xp2[i+1]=xp2[i+1]/e;
        xp3[i+1]=xp3[i+1]/e;

        p11[i+1]=(1 - k11[1 + i])*pm11[1 + i] - k12[1 + i]*pm21[1 + i] - k13[1 + i]*pm31[1 + i];
        p12[i+1]=(1 - k11[1 + i])*pm12[1 + i] - k12[1 + i]*pm22[1 + i] - k13[1 + i]*pm32[1 + i];
        p13[i+1]=(1 - k11[1 + i])*pm13[1 + i] - k12[1 + i]*pm23[1 + i] - k13[1 + i]*pm33[1 + i];
        p21[i+1]=-k21[1 + i]*pm11[1 + i] + (1 - k22[1 + i])*pm21[1 + i] - k23[1 + i]*pm31[1 + i];
        p22[i+1]=-k21[1 + i]*pm12[1 + i] + (1 - k22[1 + i])*pm22[1 + i] - k23[1 + i]*pm32[1 + i];
        p23[i+1]=-k21[1 + i]*pm13[1 + i] + (1 - k22[1 + i])*pm23[1 + i] - k23[1 + i]*pm33[1 + i];
        p31[i+1]=-k31[1 + i]*pm11[1 + i] - k32[1 + i]*pm21[1 + i] + (1 - k33[1 + i])*pm31[1 + i];
        p32[i+1]=-k31[1 + i]*pm12[1 + i] - k32[1 + i]*pm22[1 + i] + (1 - k33[1 + i])*pm32[1 + i];
        p33[i+1]=-k31[1 + i]*pm13[1 + i] - k32[1 + i]*pm23[1 + i] + (1 - k33[1 + i])*pm33[1 + i];
    };
    long double T,XX,YY,ZZ,WW,S,oXX=0,oYY=0,oZZ=0;
    res<<1<<'\t'<<0<<'\t'<<0<<'\t'<<0<<endl;
    for (i=1; i<l; i++){
        alpha[i]=atan2(-cos(gamma[i])*xp2[i]+sin(gamma[i])*xp3[i],xp1[i]/cos(beta[i]));
        //res<<gamma[i]<<" "<<beta[i]<<" "<<alpha[i]<<endl;
        matr[i][0]=cos(beta[i])*cos(alpha[i]);
        matr[i][1]=cos(beta[i])*sin(alpha[i]);
        matr[i][2]=-sin(beta[i]);
        matr[i][3]=0;
        matr[i][4]=sin(beta[i])*sin(gamma[i])*cos(alpha[i])-cos(gamma[i])*sin(alpha[i]);
        matr[i][5]=sin(beta[i])*sin(gamma[i])*sin(alpha[i]) + cos(gamma[i])*cos(alpha[i]);
        matr[i][6]=cos(beta[i])*sin(gamma[i]);
        matr[i][7]=0;
        matr[i][8]=sin(beta[i])*cos(gamma[i])*cos(alpha[i]) + sin(gamma[i])*sin(alpha[i]);
        matr[i][9]=sin(beta[i])*cos(gamma[i]) *sin(alpha[i]) - sin(gamma[i])*cos(alpha[i]);
        matr[i][10]=cos(beta[i])* cos(gamma[i]);
        matr[i][11]=0;
        matr[i][12]=0;
        matr[i][13]=0;
        matr[i][14]=0;
        matr[i][15]=1;
        T = matr[i][0] + matr[i][ 5 ] + matr[i][10] + 1;
        if (T > 0)
            {S = 0.5/sqrt(T);
            WW = 0.25/S;
            XX = (matr[i][ 9] - matr[i][ 6])*S;
            YY = (matr[i][ 2] - matr[i][ 8])*S;
            ZZ = (matr[i][ 4] - matr[i][ 1])*S;}
        else
        if (matr[i][ 0 ] >= matr[i][ 5] && matr[i][ 0 ] >= matr[i][ 10 ])
            {S = sqrt(1.0 + matr[i][ 0 ] - matr[i][ 5 ] -matr[i][ 10 ])*2;
            XX = 0.5/S;
            YY = (matr[i][ 1 ] + matr[i][ 4 ])/S;
            ZZ = (matr[i][ 2 ] + matr[i][ 8 ])/S;
            WW = (matr[i][ 6 ] + matr[i][ 9])/S;}
        else
        if (matr[i][ 5 ] >= matr[i][ 0 ] && matr[i][ 5 ] >= matr[i][ 10 ])
            {S = sqrt(1.0 + matr[i][ 5 ] - matr[i][ 0 ] - matr[i][ 10 ])*2;
            XX = (matr[i][ 1 ] + matr[i][ 4 ])/S;
            YY = 0.5/S;
            ZZ = (matr[i][ 6 ] + matr[i][ 9 ])/S;
            WW = (matr[i][ 2 ] + matr[i][ 8 ])/S;}
        else
            {S = sqrt(1.0 + matr[i][ 10 ] - matr[i][ 0 ] -matr[i][ 5 ])*2;
            XX = (matr[i][ 2 ] + matr[i][ 8 ])/S;
            YY = (matr[i][ 6 ] + matr[i][ 9 ])/S;
            ZZ = 0.5/S;
            WW = (matr[i][ 1 ] + matr[i][ 4 ])/S;};
            XX=XX/sqrt(1-WW*WW);
            YY=YY/sqrt(1-WW*WW);
            ZZ=ZZ/sqrt(1-WW*WW);
            if ( (oXX+XX)*(oXX+XX)+(oYY+YY)*(oYY+YY)+(oZZ+ZZ)*(oZZ+ZZ)<0.5) {WW=-WW; XX=-XX; YY=-YY; ZZ=-ZZ;};
            oXX=XX;
            oYY=YY;
            oZZ=ZZ;
            res<<WW<<'\t'<<-XX*sqrt(1-WW*WW)<<'\t'<<-YY*sqrt(1-WW*WW)<<'\t'<<-ZZ*sqrt(1-WW*WW)<<endl;
            //cout<<WW<<" "<<XX*sqrt(1-WW*WW)<<" "<<YY*sqrt(1-WW*WW)<<" "<<ZZ*sqrt(1-WW*WW)<<endl;
    };
    res.close();
    fl.close();
    return 0;
}
