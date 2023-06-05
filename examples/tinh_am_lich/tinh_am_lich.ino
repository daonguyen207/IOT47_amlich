#include "IOT47_amlich.h"

void setup() 
{
 Serial.begin(115200);
 int ngay_al,thang_al,nam_al;
 convertSolar2Lunar(19,11,2023,7,&ngay_al,&thang_al,&nam_al);//ngày dl, tháng dl, năm dl, múi giờ
 Serial.print("Âm lịch: " + String(ngay_al) + "/" + String(thang_al) + "/" + String(nam_al));
}

void loop() 
{

}
