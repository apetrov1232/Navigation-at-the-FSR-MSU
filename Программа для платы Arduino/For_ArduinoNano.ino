
#include <Wire.h>
#include <ADXL345.h>
#include <I2Cdev.h>
#include <L3G4200D.h>
#include <SPI.h>
#include <SD.h>
#include <HMC5883L.h>
#include <iarduino_Pressure_BMP.h>
iarduino_Pressure_BMP sensor;      
L3G4200D gyro;
ADXL345 adxl; 
HMC5883L mag;
int16_t avx, avy, avz,mx, my, mz;
File myFile;
unsigned long time=0,prevtime;
const int chipSelect = 4;
#define LED_PIN 13
void setup(){
  Wire.begin();
  SD.begin(chipSelect);
  myFile=SD.open("Data.txt", FILE_WRITE);
  myFile.println();
  myFile.println();
  myFile.println("_________NEWDATA_________");
  myFile.println();
  myFile.println();
  myFile.println();
  adxl.powerOn();
  gyro.initialize();
  gyro.setFullScale(2000);
  mag.initialize();
  sensor.begin();  
  pinMode(LED_PIN, OUTPUT);
}
void loop(){ 
  int x,y,z; float aX,aY,aZ,gX,gY,gZ,mX,mY,mZ,t,p; 
  gyro.getAngularVelocity(&avx, &avy, &avz);
  adxl.readAccel(&x, &y, &z); 
  mag.getHeading(&mx, &my, &mz);
  aX=(float)x; aY=(float)y; aZ=(float)z; 
  gX=(float)avx; gY=(float)avy; gZ=(float)avz;
  mX=(float)mx; mY=(float)my; mZ=(float)mz;
  sensor.read(1); t=(float)sensor.type;
  sensor.read(2); p=(float)sensor.type;
  aX=0.0382*(aX-3.5);
  aY=0.03815*(aY-4.95);
  aZ=0.039*(aZ-5.6);
  mX=0.00164*(mX-30); 
  mY=0.0015*(mY+90);
  mZ=0.0018*(mZ-115);
  gX=(gX-14.647)*0.0675;
  gY=(gY+4.498)*0.0606;
  gZ=(gZ-13.09)*0.06448;
  myFile.print(String(aX));
  myFile.print(" ");
  myFile.print(String(aY));
  myFile.print(" ");
  myFile.print(String(aZ));
  myFile.print(" ");
  myFile.print(String(gX));
  myFile.print(" ");
  myFile.print(String(gY));
  myFile.print(" ");
  myFile.print(String(gZ));
  myFile.print(" ");
  myFile.print(String(mX));
  myFile.print(" ");
  myFile.print(String(mY));
  myFile.print(" ");
  myFile.print(String(mZ));
  //myFile.print(" ");
  //myFile.print(String(t));
  //myFile.print(" ");
  //myFile.print(String(p));
  myFile.print(" ");
  prevtime=time; 
  time=millis(); 
  myFile.println(time);
  myFile.flush();
}
