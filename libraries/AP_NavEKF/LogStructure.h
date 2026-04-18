#pragma once

#include <AP_Logger/LogStructure.h>

// @LoggerMessage: XKY0,NKY0
// @Description: EKF Yaw Estimator States
// @Field: TimeUS: Time since system startup
// @Field: C: EKF core this data is for
// @Field: YC: GSF yaw estimate (deg)
// @Field: YCS: GSF yaw estimate 1-Sigma uncertainty (deg)
// @Field: Y0: Yaw estimate from individual EKF filter 0 (deg)
// @Field: Y1: Yaw estimate from individual EKF filter 1 (deg)
// @Field: Y2: Yaw estimate from individual EKF filter 2 (deg)
// @Field: Y3: Yaw estimate from individual EKF filter 3 (deg)
// @Field: Y4: Yaw estimate from individual EKF filter 4 (deg)
// @Field: W0: Weighting applied to yaw estimate from individual EKF filter 0
// @Field: W1: Weighting applied to yaw estimate from individual EKF filter 1
// @Field: W2: Weighting applied to yaw estimate from individual EKF filter 2
// @Field: W3: Weighting applied to yaw estimate from individual EKF filter 3
// @Field: W4: Weighting applied to yaw estimate from individual EKF filter 4
struct PACKED log_KY0 {
    LOG_PACKET_HEADER;
    uint64_t time_us;
    uint8_t core;
    float yaw_composite;
    float yaw_composite_variance;
    float yaw0;
    float yaw1;
    float yaw2;
    float yaw3;
    float yaw4;
    float wgt0;
    float wgt1;
    float wgt2;
    float wgt3;
    float wgt4;
};


// @LoggerMessage: XKY1,NKY1
// @Description: EKF Yaw Estimator Innovations
// @Field: TimeUS: Time since system startup
// @Field: C: EKF core this data is for
// @Field: IVN0: North velocity innovation from individual EKF filter 0 (m/s)
// @Field: IVN1: North velocity innovation from individual EKF filter 1 (m/s)
// @Field: IVN2: North velocity innovation from individual EKF filter 2 (m/s)
// @Field: IVN3: North velocity innovation from individual EKF filter 3 (m/s)
// @Field: IVN4: North velocity innovation from individual EKF filter 4 (m/s)
// @Field: IVE0: East velocity innovation from individual EKF filter 0 (m/s)
// @Field: IVE1: East velocity innovation from individual EKF filter 1 (m/s)
// @Field: IVE2: East velocity innovation from individual EKF filter 2 (m/s)
// @Field: IVE3: East velocity innovation from individual EKF filter 3 (m/s)
// @Field: IVE4: East velocity innovation from individual EKF filter 4 (m/s)
struct PACKED log_KY1 {
    LOG_PACKET_HEADER;
    uint64_t time_us;
    uint8_t core;
    float ivn0;
    float ivn1;
    float ivn2;
    float ivn3;
    float ivn4;
    float ive0;
    float ive1;
    float ive2;
    float ive3;
    float ive4;
};

// @LoggerMessage: XKY2,NKY2
// @Description: EKF Yaw Estimator Position Innovations
// @Field: TimeUS: Time since system startup
// @Field: C: EKF core this data is for
// @Field: IPN0: North position innovation from individual EKF filter 0 (m)
// @Field: IPN1: North position innovation from individual EKF filter 1 (m)
// @Field: IPN2: North position innovation from individual EKF filter 2 (m)
// @Field: IPN3: North position innovation from individual EKF filter 3 (m)
// @Field: IPN4: North position innovation from individual EKF filter 4 (m)
// @Field: IPE0: East position innovation from individual EKF filter 0 (m)
// @Field: IPE1: East position innovation from individual EKF filter 1 (m)
// @Field: IPE2: East position innovation from individual EKF filter 2 (m)
// @Field: IPE3: East position innovation from individual EKF filter 3 (m)
// @Field: IPE4: East position innovation from individual EKF filter 4 (m)
// @Field: PN: North component of position measurement input to fusePosData (m)
// @Field: PE: East component of position measurement input to fusePosData (m)
struct PACKED log_KY2 {
    LOG_PACKET_HEADER;
    uint64_t time_us;
    uint8_t core;
    float ipn0;
    float ipn1;
    float ipn2;
    float ipn3;
    float ipn4;
    float ipe0;
    float ipe1;
    float ipe2;
    float ipe3;
    float ipe4;
    float pin;
    float pie;
};

#define KY0_FMT "QBffffffffffff"
#define KY0_LABELS "TimeUS,C,YC,YCS,Y0,Y1,Y2,Y3,Y4,W0,W1,W2,W3,W4"
#define KY0_UNITS "s#hdhhhhh-----"
#define KY0_MULTS "F-000000000000"

#define KY1_FMT "QBffffffffff"
#define KY1_LABELS "TimeUS,C,IVN0,IVN1,IVN2,IVN3,IVN4,IVE0,IVE1,IVE2,IVE3,IVE4"
#define KY1_UNITS "s#nnnnnnnnnn"
#define KY1_MULTS "F-0000000000"

#define KY2_FMT "QBffffffffffff"
#define KY2_LABELS "TimeUS,C,IPN0,IPN1,IPN2,IPN3,IPN4,IPE0,IPE1,IPE2,IPE3,IPE4,PN,PE"
#define KY2_UNITS "s#mmmmmmmmmmmm"
#define KY2_MULTS "F-000000000000"

#define LOG_STRUCTURE_FROM_NAVEKF                                       \
    { LOG_XKY0_MSG, sizeof(log_KY0),                                    \
            "XKY0", KY0_FMT, KY0_LABELS, KY0_UNITS, KY0_MULTS, true },      \
    { LOG_XKY1_MSG, sizeof(log_KY1),                                    \
            "XKY1", KY1_FMT, KY1_LABELS, KY1_UNITS, KY1_MULTS , true },     \
    { LOG_XKY2_MSG, sizeof(log_KY2),                                    \
            "XKY2", KY2_FMT, KY2_LABELS, KY2_UNITS, KY2_MULTS , true },     \
    { LOG_NKY0_MSG, sizeof(log_KY0),                                    \
            "NKY0", KY0_FMT, KY0_LABELS, KY0_UNITS, KY0_MULTS, true },      \
    { LOG_NKY1_MSG, sizeof(log_KY1),                                    \
            "NKY1", KY1_FMT, KY1_LABELS, KY1_UNITS, KY1_MULTS , true },     \
    { LOG_NKY2_MSG, sizeof(log_KY2),                                    \
            "NKY2", KY2_FMT, KY2_LABELS, KY2_UNITS, KY2_MULTS , true },
