#ifndef __IOT47_AMLICH_
#define __IOT47_AMLICH_

/*thư viện tính âm lịch
 Người viết: Đào nguyện iot47
 Email: daonguyen20798@gmail.com
 Website: iot47.com
*/

#ifndef M_PI
#define M_PI 3.1415926535898
#endif

int INT(double d) {
    return floor(d);
}
double jdFromDate(int dd, int mm, int yy) {
    double a = INT((14.0 - mm) / 12.0);
    double y = yy + 4800 - a;
    double m = mm + 12.0 * a - 3.0;
    double jd = dd + INT((153.0 * m + 2.0) / 5.0) + 365.0 * y + INT(y / 4.0) - INT(y / 100.0) + INT(y / 400.0) - 32045.0;
    if (jd < 2299161) {
        jd = dd + INT((153.0* m + 2.0)/5.0) + 365.0 * y + INT(y / 4.0) - 32083.0;
    }
    return jd;
}
void jdToDate(double jd) 
{
    double a,b,c;
    if (jd > 2299160) { // After 5/10/1582, Gregorian calendar
         a = jd + 32044;
         b = INT((4.0*a+3.0)/146097.0);
         c = a - INT((b*146097.0)/4.0);
    } else {
         b = 0;
         c = jd + 32082;
    }
    double d = INT((4.0*c+3.0)/1461.0);
    double e = c - INT((1461.0*d)/4.0);
    double m = INT((5.0*e+2.0)/153.0);
    double day = e - INT((153.0*m+2.0)/5.0) + 1;
    double month = m + 3.0 - 12.0*INT(m/10.0);
    double year = b*100.0 + d - 4800.0 + INT(m/10.0);
}
int32_t getNewMoonDay(double k, double timeZone) {
    double T = k/1236.85; // Time in Julian centuries from 1900 January 0.5
    double T2 = T * T;
    double T3 = T2 * T;
    double dr = M_PI/180;
    double Jd1 = 2415020.75933 + 29.53058868*k + 0.0001178*T2 - 0.000000155*T3;
           Jd1 = Jd1 + 0.00033*sin((166.56 + 132.87*T - 0.009173*T2)*dr); // Mean new moon
    double M = 359.2242 + 29.10535608*k - 0.0000333*T2 - 0.00000347*T3; // Sun's mean anomaly
    double Mpr = 306.0253 + 385.81691806*k + 0.0107306*T2 + 0.00001236*T3; // Moon's mean anomaly
    double F = 21.2964 + 390.67050646*k - 0.0016528*T2 - 0.00000239*T3; // Moon's argument of latitude
    double C1=(0.1734 - 0.000393*T)*sin(M*dr) + 0.0021*sin(2*dr*M);
    C1 = C1 - 0.4068*sin(Mpr*dr) + 0.0161*sin(dr*2.0*Mpr);
    C1 = C1 - 0.0004*sin(dr*3*Mpr);
    C1 = C1 + 0.0104*sin(dr*2*F) - 0.0051*sin(dr*(M+Mpr));
    C1 = C1 - 0.0074*sin(dr*(M-Mpr)) + 0.0004*sin(dr*(2.0*F+M));
    C1 = C1 - 0.0004*sin(dr*(2*F-M)) - 0.0006*sin(dr*(2.0*F+Mpr));
    C1 = C1 + 0.0010*sin(dr*(2*F-Mpr)) + 0.0005*sin(dr*(2.0*Mpr+M));
    double deltat;
    if (T < -11) {
        deltat= 0.001 + 0.000839*T + 0.0002261*T2 - 0.00000845*T3 - 0.000000081*T*T3;
    } else {
        deltat= -0.000278 + 0.000265*T + 0.000262*T2;
    };
    double JdNew = Jd1 + C1 - deltat;
    //echo "JdNew = $JdNew\n";
    return INT(JdNew + 0.5 + timeZone/24.0);
}
int32_t getSunLongitude(double jdn, double timeZone) {
    double T = (jdn - 2451545.5 - timeZone/24) / 36525; // Time in Julian centuries from 2000-01-01 12:00:00 GMT
    double T2 = T * T;
    double dr = M_PI/180.0; // degree to radian
    double M = 357.52910 + 35999.05030*T - 0.0001559*T2 - 0.00000048*T*T2; // mean anomaly, degree
    double L0 = 280.46645 + 36000.76983*T + 0.0003032*T2; // mean longitude, degree
    double DL = (1.914600 - 0.004817*T - 0.000014*T2)*sin(dr*M);
           DL = DL + (0.019993 - 0.000101*T)*sin(dr*2.0*M) + 0.000290*sin(dr*3.0*M);
    double L = L0 + DL; // true longitude, degree
    //echo "\ndr = $dr, M = $M, T = $T, DL = $DL, L = $L, L0 = $L0\n";
    L = L*dr;
    L = L - M_PI*2*(INT(L/(M_PI*2.0))); // Normalize to (0, 2*PI)
    return INT(L/M_PI*6.0);
}
int32_t getLunarMonth11(double yy, double timeZone) {
    double off = jdFromDate(31, 12, yy) - 2415021.0;
    double k = INT(off / 29.530588853);
    double nm = getNewMoonDay(k, timeZone);
    double sunLong = getSunLongitude(nm, timeZone); // sun longitude at local midnight
    if (sunLong >= 9) {
        nm = getNewMoonDay(k-1, timeZone);
    }
    return nm;
}
int32_t getLeapMonthOffset(double a11, double timeZone) {
    double k = INT((a11 - 2415021.076998695) / 29.530588853 + 0.5);
    double last = 0;
    double i = 1; // We start with the month following lunar month 11
    double arc = getSunLongitude(getNewMoonDay(k + i, timeZone), timeZone);
    do {
        last = arc;
        i = i + 1;
        arc = getSunLongitude(getNewMoonDay(k + i, timeZone), timeZone);
    } while (arc != last && i < 14);
    return i - 1;
}
void convertSolar2Lunar(double dd, double mm, double yy, double timeZone,int * ngay_al,int * thang_al, int *nam_al) {
    double lunarYear;
    double dayNumber = jdFromDate(dd, mm, yy);
    double k = INT((dayNumber - 2415021.076998695) / 29.530588853);
    double monthStart = getNewMoonDay(k+1, timeZone);
    if (monthStart > dayNumber) {
        monthStart = getNewMoonDay(k, timeZone);
    }
    double a11 = getLunarMonth11(yy, timeZone);
    double b11 = a11;
    if (a11 >= monthStart) {
        lunarYear = yy;
        a11 = getLunarMonth11(yy-1, timeZone);
    } else {
        lunarYear = yy+1;
        b11 = getLunarMonth11(yy+1, timeZone);
    }
    double lunarDay = dayNumber - monthStart + 1;
    double diff = INT((monthStart - a11)/29);
    double lunarLeap = 0;
    double lunarMonth = diff + 11;
    if (b11 - a11 > 365) {
        double leapMonthDiff = getLeapMonthOffset(a11, timeZone);
        if (diff >= leapMonthDiff) {
            lunarMonth = diff + 10;
            if (diff == leapMonthDiff) {
                lunarLeap = 1;
            }
        }
    }
    if (lunarMonth > 12) {
        lunarMonth = lunarMonth - 12;
    }
    if (lunarMonth >= 11 && diff < 4) {
        lunarYear -= 1;
    }
    *ngay_al = (int)lunarDay;
    *thang_al = (int)lunarMonth;
    *nam_al = (int)lunarYear;
}



#endif