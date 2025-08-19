#include "AOI.h"
#include <halconcpp/HalconCpp.h>
using namespace HalconCpp;

// 深度图处理函数，用于深度图倾斜矫正和灰度转换
// 输入原始深度图和用于平面抽取的四个圆形区域，输出转换为灰度图后的深度图
HObject DepthImageHandler(HObject Imagedepth, HObject UnionedCircles) {
    HObject ho_Domain, ho_ReferencePlaneDistance;
    HObject ho_ImagedepthWithoutTilt, ho_ImagedepthCorrected, ho_ImageScaled;
    HObject ho_ImageConverted;
    HTuple hv_MRow, hv_MCol, hv_Alpha, hv_Beta, hv_Mean, hv_Pointer, hv_Type;
    HTuple hv_Width, hv_Height, hv_Area, hv_Row1, hv_Column1, hv_CosGamma;
    HTuple hv_MaxVal, hv_MinVal, hv_scaleVal;

    GetDomain(Imagedepth, &ho_Domain);
    Intersection(UnionedCircles, Imagedepth, &UnionedCircles);
    MomentsGrayPlane(UnionedCircles, Imagedepth, &hv_MRow, &hv_MCol, &hv_Alpha, &hv_Beta, &hv_Mean);
    GetImagePointer1(Imagedepth, &hv_Pointer, &hv_Type, &hv_Width, &hv_Height);
    AreaCenter(UnionedCircles, &hv_Area, &hv_Row1, &hv_Column1);
    GenImageSurfaceFirstOrder(&ho_ReferencePlaneDistance, hv_Type, hv_Alpha, hv_Beta, hv_Mean, hv_Row1, hv_Column1, hv_Width, hv_Height);
    SubImage(Imagedepth, ho_ReferencePlaneDistance, &ho_ImagedepthWithoutTilt, 1, 0);
    hv_CosGamma = 1.0 / ((((hv_Alpha * hv_Alpha) + (hv_Beta * hv_Beta)) + 1).TupleSqrt());
    ScaleImage(ho_ImagedepthWithoutTilt, &ho_ImagedepthCorrected, hv_CosGamma, 0);
    hv_MaxVal = 1;
    hv_MinVal = -1;
    hv_scaleVal = 255.0 / (hv_MaxVal - hv_MinVal);
    ScaleImage(ho_ImagedepthCorrected, &ho_ImageScaled, hv_scaleVal, hv_scaleVal * (-hv_MinVal));
    ConvertImageType(ho_ImageScaled, &ho_ImageConverted, "uint2");

    return ho_ImageConverted;
}

// 区别使用AB参数的焊点号处理函数
// 输入QVariantList形式的焊点号、该拍照点位的焊点个数，输出HTuple类型的焊点号序列，0代表使用A参数，1代表使用B参数
HTuple SParrayHandler(QVariantList sparray, int num) {
    HTuple hv_SParrayHandled;
    for (int i = 0; i < num; i++)
        hv_SParrayHandled[i] = 0;
    // 若焊点号为空则都使用A参数
    if (sparray.length() == 0) 
        return hv_SParrayHandled;
    else {
        // spdlog::error("输入的焊点号：");
        for (int j = 0; j < sparray.length(); j++) {
            int index = sparray.at(j).toInt();
            // spdlog::error(index);
            // 若焊点号超限则将序列的第一位置-1并返回
            if (index > num || index <= 0) {
                hv_SParrayHandled[0] = -1;
                return hv_SParrayHandled;
            }
            hv_SParrayHandled[index-1] = 1;
        }
        // spdlog::error("焊点号结束");
        return hv_SParrayHandled;
    }
}

// 缺焊检测函数
// 输入处理过的深度图、基准点坐标、圆形ROI、区分AB参数的焊点序列、焊点深度阈值、面积阈值，返回NG序列
HTuple CircleROIHandler(HObject ho_Imagedepth_Handled, HTuple hv_Row, HTuple hv_Column, HTuple handledsparray, float *CircleX, float *CircleY, float *CircleR, float ROIminThreshold, float KHminAreaA, float KHminAreaB) {

    HObject ho_Circle, ho_ImageReduced, ho_Region;
    HTuple hv_Value, hv_CircleNG;

    for (int i = 0; i < handledsparray.Length(); i++) {
        GenCircle(&ho_Circle, hv_Row + HTuple(CircleX[i]), hv_Column + HTuple(CircleY[i]), HTuple(CircleR[i]));
        ReduceDomain(ho_Imagedepth_Handled, ho_Circle, &ho_ImageReduced);
        Threshold(ho_ImageReduced, &ho_Region, HTuple(ROIminThreshold), 1023);
        RegionFeatures(ho_Region, "area", &hv_Value);
        // 若空焊灰度阈值内的区域面积大于空焊面积阈值则为NG
        if (handledsparray[i] == 0) {
            if (hv_Value <= HTuple(KHminAreaA)) {
                hv_CircleNG[i] = 1;
                spdlog::error("焊点{}缺焊面积{}, 阈值为A{}", i + 1, hv_Value[0].D(), KHminAreaA);
            }
            else
                hv_CircleNG[i] = 0;
        } 
        else {
            if (hv_Value <= HTuple(KHminAreaB)) {
                hv_CircleNG[i] = 1;
                spdlog::error("焊点{}缺焊面积{}, 阈值为B{}", i + 1, hv_Value[0].D(), KHminAreaB);
            } else
                hv_CircleNG[i] = 0;
        }
    }
    return hv_CircleNG;
}

// 连锡检测函数
// 输入处理过的深度图、基准点坐标、矩形ROI坐标偏移量、连锡点面积阈值，返回连锡区域
HObject RectangleROIHandler(HObject ho_Imagedepth_Handled, HTuple hv_Row, HTuple hv_Column, float *RectangleXY, float ROIminThreshold, float LXminArea, float LXmaxArea) {
    HObject ho_Rectangle, ho_ImageReduced, ho_Regions, ho_ConnectedRegions, ho_SelectedRegions;
    HTuple hv_Number, hv_AreaNG, hv_RowNG, hv_ColumnNG;
    GenRectangle1(&ho_Rectangle, hv_Row + RectangleXY[0], hv_Column + RectangleXY[1], hv_Row + RectangleXY[2], hv_Column + RectangleXY[3]);
    ReduceDomain(ho_Imagedepth_Handled, ho_Rectangle, &ho_ImageReduced);
    Threshold(ho_ImageReduced, &ho_Regions, HTuple(ROIminThreshold), 1023);
    Connection(ho_Regions, &ho_ConnectedRegions);
    SelectShape(ho_ConnectedRegions, &ho_SelectedRegions, "area", "and", HTuple(LXminArea), HTuple(LXmaxArea));
    AreaCenter(ho_SelectedRegions, &hv_AreaNG, &hv_RowNG, &hv_ColumnNG);
    for (int i = 0; i < ho_SelectedRegions.CountObj(); i++)
        spdlog::error("连锡区域{}面积{},阈值为{}到{}", i + 1, hv_AreaNG[i].D(), LXminArea, LXmaxArea);
    return ho_SelectedRegions;
}

// 缺焊连锡检测算法
AlgoOutputPtr AOI::compute(const QVariantMap& execParam){
    using namespace HalconCpp;
    auto res = new AOIOutput;
    try {
        // 从系统中获取2D图image和深度图depthimage
        BImage image = execParam.value("image").value<BImage>();
        BImage depthimage = execParam.value("depthimage").value<BImage>();
        int index = execParam.value("index").value<int>();
        // 拍照点1检测
        if (index == 1) {
            // 焊点个数F
            const int num = 7;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray1 = getParamValue<QVariantList>("SParray1");
            float ROIminThreshold1 = getParamValue<float>("ROIminThreshold1");
            float KHminAreaA1 = getParamValue<float>("KHminAreaA1");
            float KHminAreaB1 = getParamValue<float>("KHminAreaB1");
            float LXminArea1 = getParamValue<float>("LXminArea1");
            float LXmaxArea1 = getParamValue<float>("LXmaxArea1");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-1075.81, -1063.29, -1069.26, -1083.77, -1061.09, -1058.03, -1087.07};
            float CircleY[num] = {-1039.312, -597.22, -221.68, 128.98, 562.17, 946.63, 1300.62};
            float CircleR[num] = {96, 131, 127, 93, 113, 111, 90};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 100, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 12000, 90000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "ratio", "and", 0.5, 0.7);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 2000, 2600);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1000, 1600);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            float BaseX = hv_Column[0].D();
            float BaseY = hv_Row[0].D();
            res->setBaseROI(ROI::ROI(BaseX, BaseY, 10, 10));
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 1361.599).TupleConcat(hv_Row - 1364.573)).TupleConcat(hv_Row + 101.87)).TupleConcat(hv_Row - 245.28),
                (((hv_Column + 42.35).TupleConcat(hv_Column + 1096.68)).TupleConcat(hv_Column - 828.913)).TupleConcat(hv_Column + 1159.27),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray1, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold1, KHminAreaA1, KHminAreaB1);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];    
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-1215.54, -1200.796, -926.61, 1451.74};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold1, LXminArea1, LXmaxArea1);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {

                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点2检测
        else if (index == 2) {
            // 焊点个数
            const int num = 7;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray2 = getParamValue<QVariantList>("SParray2");
            float ROIminThreshold2 = getParamValue<float>("ROIminThreshold2");
            float KHminAreaA2 = getParamValue<float>("KHminAreaA2");
            float KHminAreaB2 = getParamValue<float>("KHminAreaB2");
            float LXminArea2 = getParamValue<float>("LXminArea2");
            float LXmaxArea2 = getParamValue<float>("LXmaxArea2");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-1123.3, -1117.2, -1123.9, -1128.46, -1124.8, -1129.5, -1133.2};
            float CircleY[num] = {-1773.7, -1414.4, -960.7, -600.6, -244.8, 213.3, 565.3};
            float CircleR[num] = {84, 82, 80, 80, 85, 75, 82};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 100, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 4000, 10000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "ratio", "and", 10, 20);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 2000, 3000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "height", "and", 400, 600);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 251.29).TupleConcat(hv_Row - 174.9)).TupleConcat(hv_Row - 1582.27)).TupleConcat(hv_Row - 1361.01),
                (((hv_Column + 383.12).TupleConcat(hv_Column - 1386.178)).TupleConcat(hv_Column - 854.41)).TupleConcat(hv_Column + 373.15),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray2, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold2, KHminAreaA2, KHminAreaB2);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-1267.66, -1899.524, -1023.55, 719.36};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold2, LXminArea2, LXmaxArea2);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点3检测
        else if (index == 3) {
            // 焊点个数
            const int num = 7;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray3 = getParamValue<QVariantList>("SParray3");
            float ROIminThreshold3 = getParamValue<float>("ROIminThreshold3");
            float KHminAreaA3 = getParamValue<float>("KHminAreaA3");
            float KHminAreaB3 = getParamValue<float>("KHminAreaB3");
            float LXminArea3 = getParamValue<float>("LXminArea3");
            float LXmaxArea3 = getParamValue<float>("LXmaxArea3");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {983.6, 987.5, 983.6, 982.2, 980.2, 975.8, 980.6};
            float CircleY[num] = {-1634.5, -1282.8, -920.35, -468.2, -111, 251.8, 706.5};
            float CircleR[num] = {75, 80, 80, 75, 85, 85, 85};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
             // 查找基准点
             Threshold(ho_Image, &ho_Regions0, 0, 22);
             Connection(ho_Regions0, &ho_ConnectedRegions0);
             SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 100000, 200000);
             CountObj(ho_SelectedRegions0, &hv_Number);
             if (hv_Number[0].I() != 1) {
                 spdlog::error("未找到基准点");
                 circlecount = -1;
                 rectanglecount = -1;
                 res->OK = false;
                 res->NGnum(circlecount, rectanglecount);
                 return AlgoOutputPtr(res);
             }
             AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
             // 读取深度图
             auto ho_Imagedepth = BImage2HalconImage(depthimage);
             // 处理深度图
             GenCircle(&ho_Circle,
                 (((hv_Row + 447.53).TupleConcat(hv_Row + 427.151)).TupleConcat(hv_Row + 1812.937)).TupleConcat(hv_Row + 2047.297),
                 (((hv_Column + 92.61).TupleConcat(hv_Column - 941.61)).TupleConcat(hv_Column - 1186.15)).TupleConcat(hv_Column - 126.46),
                 (((HTuple(15).Append(15)).Append(15)).Append(15)));
             Union1(ho_Circle, &ho_UnionedCircles);
             HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
             // 处理焊点号
             hv_SParrayHandled = SParrayHandler(SParray3, num);
             if (hv_SParrayHandled[0] == -1) {
                 res->OK = false;
                 spdlog::error("焊点号超限");
                 return AlgoOutputPtr(res);
             }
             // 处理圆形ROI，进行缺焊检测
             hv_CircleNG = CircleROIHandler(
                 ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold3, KHminAreaA3, KHminAreaB3);
             // 处理NG序列CircleNG[]
             for (int i = 0; i < num; i++) {
                 if (hv_CircleNG[i] == 1) {
                     res->OK = false;
                     // 生成圆形roi的NG识别框
                     circlecount++;
                     int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                     int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                     int W = CircleR[i] * 2;
                     int H = CircleR[i] * 2;
                     res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                 }
             }
             // spdlog::info(circlecount);
             // 矩形ROI坐标偏移量
             float RectangleXY[] = {858.757, -1803.181, 1113.497, 805.29};
             // 处理矩形ROI，进行连锡检测
             ho_RectangleNGRegion =
                 RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold3, LXminArea3, LXmaxArea3);
             if (ho_RectangleNGRegion.CountObj() > 0) {
                 res->OK = false;
                 SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                 for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                     rectanglecount++;
                     int X = hv_Column1[i].D();
                     int Y = hv_Row1[i].D();
                     int W = hv_Column2[i].D() - hv_Column1[i].D();
                     int H = hv_Row2[i].D() - hv_Row1[i].D();
                     res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                 }
             }
             // spdlog::info(rectanglecount);
             // 若无缺焊和连锡则为OK
             if (circlecount == 0 && rectanglecount == 0) {
                 res->OK = true;
                 return AlgoOutputPtr(res);
             } else {
                 res->NGnum(circlecount, rectanglecount);
                 return AlgoOutputPtr(res);
             }
             return AlgoOutputPtr(res);
        }
        // 拍照点4检测
        else if (index == 4) {
            // 焊点个数
            const int num = 7;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray4 = getParamValue<QVariantList>("SParray4");
            float ROIminThreshold4 = getParamValue<float>("ROIminThreshold4");
            float KHminAreaA4 = getParamValue<float>("KHminAreaA4");
            float KHminAreaB4 = getParamValue<float>("KHminAreaB4");
            float LXminArea4 = getParamValue<float>("LXminArea4");
            float LXmaxArea4 = getParamValue<float>("LXmaxArea4");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles,
                ho_RegionOpening0;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-599.855, -603.84, -604.04, -603.37, -603.87, -607.86, -606.61};
            float CircleY[num] = {-1709.217, -1252.586, -898.09, -533, -80.24, 286.02, 632.48};
            float CircleR[num] = {94.1476, 82.5617, 82.7012, 84.8149, 85, 80, 89.6064};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 160, 255);
            OpeningCircle(ho_Regions0, &ho_RegionOpening0, 3.5);
            Connection(ho_RegionOpening0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 20000, 60000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 1400, 2000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1700, 2500);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 1377.244).TupleConcat(hv_Row - 805.59)).TupleConcat(hv_Row + 260)).TupleConcat(hv_Row + 379.27),
                (((hv_Column - 1066.276).TupleConcat(hv_Column + 130.04)).TupleConcat(hv_Column - 1012.539)).TupleConcat(hv_Column + 27.98),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray4, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold4, KHminAreaA4, KHminAreaB4);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-744.04, -1847.961, -484.21, 759.94};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold4, LXminArea4, LXmaxArea4);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点5检测
        else if (index == 5) {
            // 焊点个数
            const int num = 2;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray5 = getParamValue<QVariantList>("SParray5");
            float ROIminThreshold5 = getParamValue<float>("ROIminThreshold5");
            float KHminAreaA5 = getParamValue<float>("KHminAreaA5");
            float KHminAreaB5 = getParamValue<float>("KHminAreaB5");
            float LXminArea5 = getParamValue<float>("LXminArea5");
            float LXmaxArea5 = getParamValue<float>("LXmaxArea5");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-605.35, -606.96};
            float CircleY[num] = {-200.06, 161.19};
            float CircleR[num] = {81.9562, 87.1006};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 242, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 15000, 40000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1500, 99999);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 906.903).TupleConcat(hv_Row - 912)).TupleConcat(hv_Row + 164.07)).TupleConcat(hv_Row + 335.5),
                (((hv_Column - 584.39).TupleConcat(hv_Column - 79.94)).TupleConcat(hv_Column - 554.09)).TupleConcat(hv_Column + 47.79),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray5, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold5, KHminAreaA5, KHminAreaB5);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-762.12, -376.13, -474.16, 299.41};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold5, LXminArea5, LXmaxArea5);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        } 
        // 拍照点6检测
        else if (index == 6) {
            // 焊点个数
            const int num = 7;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray6 = getParamValue<QVariantList>("SParray6");
            float ROIminThreshold6 = getParamValue<float>("ROIminThreshold6");
            float KHminAreaA6 = getParamValue<float>("KHminAreaA6");
            float KHminAreaB6 = getParamValue<float>("KHminAreaB6");
            float LXminArea6 = getParamValue<float>("LXminArea6");
            float LXmaxArea6 = getParamValue<float>("LXmaxArea6");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {408.407, 407.127, 408.687, 404.657, 412.617, 407.847, 402.537};
            float CircleY[num] = {-1680.07, -1316.703, -963.4, -511.2, -147.25, 211.06, 662.44};
            float CircleR[num] = {74, 131.359, 116.243, 95.7593, 116.904, 121.323, 87.3259};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 117, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 30000, 80000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 600, 1300);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1500, 2500);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 255.915).TupleConcat(hv_Row - 55.39)).TupleConcat(hv_Row + 756.336)).TupleConcat(hv_Row + 1158.426),
                (((hv_Column - 1716.394).TupleConcat(hv_Column + 506.57)).TupleConcat(hv_Column + 101.61)).TupleConcat(hv_Column - 327.04),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray6, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold6, KHminAreaA6, KHminAreaB6);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {256.716, -1800.712, 548.666, 793.34};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold6, LXminArea6, LXmaxArea6);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        } 
        // 拍照点7检测
        else if (index == 7) {
            // 焊点个数
            const int num = 8;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray7 = getParamValue<QVariantList>("SParray7");
            float ROIminThreshold7 = getParamValue<float>("ROIminThreshold7");
            float KHminAreaA7 = getParamValue<float>("KHminAreaA7");
            float KHminAreaB7 = getParamValue<float>("KHminAreaB7");
            float LXminArea7 = getParamValue<float>("LXminArea7");
            float LXmaxArea7 = getParamValue<float>("LXmaxArea7");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles, ho_RegionClosing0;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-1195.18, -1205.94, -1193.39, -1183.69, -1192.54, -1183.48, -1185.47, -1194.29};
            float CircleY[num] = {-1363.125, -1006.769, -564.686, -194.96, 155.09, 622.45, 879.31, 1112.97};
            float CircleR[num] = {86.9149, 86.8721, 77.8093, 71.2666, 75.1604, 117.3389, 131.528, 110.57};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 100, 255);
            ClosingCircle(ho_Regions0, &ho_RegionClosing0, 4.5);
            Connection(ho_RegionClosing0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 5000, 10000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "ratio", "and", 0, 0.2);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 2300, 2700);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 1939.092).TupleConcat(hv_Row - 1917.61)).TupleConcat(hv_Row - 257.18)).TupleConcat(hv_Row - 1052.1),
                (((hv_Column - 1031.678).TupleConcat(hv_Column + 268.78)).TupleConcat(hv_Column + 525.6)).TupleConcat(hv_Column - 829.301),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray7, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold7, KHminAreaA7, KHminAreaB7);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-1328.77, -1548.194, -1028.37, 1206.88};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold7, LXminArea7, LXmaxArea7);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点8检测
        else if (index == 8) {
            // 焊点个数
            const int num = 9;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray8 = getParamValue<QVariantList>("SParray8");
            float ROIminThreshold8 = getParamValue<float>("ROIminThreshold8");
            float KHminAreaA8 = getParamValue<float>("KHminAreaA8");
            float KHminAreaB8 = getParamValue<float>("KHminAreaB8");
            float LXminArea8 = getParamValue<float>("LXminArea8");
            float LXmaxArea8 = getParamValue<float>("LXmaxArea8");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles,
                ho_RegionClosing0;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-288.3, -269.3, -295.15, -290.97, -285.36, -284.35, -305.05, -305.05, -299.83};
            float CircleY[num] = {-1012.623, -750.689, -492.215, -113.88, 135.42, 418.97, 776.56, 1037.71, 1304.08};
            float CircleR[num] = {100, 100, 70, 100, 100, 120.211, 100, 100, 119.287};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 174, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 5000, 11000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "ratio", "and", 3, 5);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 1200, 2000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1000, 8000);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 1064.217).TupleConcat(hv_Row - 1096.974)).TupleConcat(hv_Row + 870.03)).TupleConcat(hv_Row + 793.97),
                (((hv_Column - 825.36).TupleConcat(hv_Column + 644.89)).TupleConcat(hv_Column - 923.448)).TupleConcat(hv_Column + 841.95),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray8, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold8, KHminAreaA8, KHminAreaB8);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-435.27, -1130.614, -134.58, 1422.52};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold8, LXminArea8, LXmaxArea8);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点9检测
        else if (index == 9) {
            // 焊点个数
            const int num = 2;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray9 = getParamValue<QVariantList>("SParray9");
            float ROIminThreshold9 = getParamValue<float>("ROIminThreshold9");
            float KHminAreaA9 = getParamValue<float>("KHminAreaA9");
            float KHminAreaB9 = getParamValue<float>("KHminAreaB9");
            float LXminArea9 = getParamValue<float>("LXminArea9");
            float LXmaxArea9 = getParamValue<float>("LXmaxArea9");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles,
                ho_RegionOpenings0;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {778.529, 1175.859};
            float CircleY[num] = {-783.91, -815.51};
            float CircleR[num] = {179.941, 169.318};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 200, 255);
            OpeningCircle(ho_Regions0, &ho_RegionOpenings0, 3.5);
            Connection(ho_RegionOpenings0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 30000, 110000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 2200, 2500);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 750, 950);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row + 338.722).TupleConcat(hv_Row + 657.112)).TupleConcat(hv_Row + 1483.002)).TupleConcat(hv_Row + 1023.932),
                (((hv_Column - 164.44).TupleConcat(hv_Column - 1676.357)).TupleConcat(hv_Column - 691.16)).TupleConcat(hv_Column - 1825.699),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray9, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold9, KHminAreaA9, KHminAreaB9);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {644.022, -1062.03, 1311.312, -534.01};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold9, LXminArea9, LXmaxArea9);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        } 
        // 拍照点10检测
        else if (index == 10) {
            // 焊点个数
            const int num = 4;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray10 = getParamValue<QVariantList>("SParray10");
            float ROIminThreshold10 = getParamValue<float>("ROIminThreshold10");
            float KHminAreaA10 = getParamValue<float>("KHminAreaA10");
            float KHminAreaB10 = getParamValue<float>("KHminAreaB10");
            float LXminArea10 = getParamValue<float>("LXminArea10");
            float LXmaxArea10 = getParamValue<float>("LXmaxArea10");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {362.886, 354.728, 795.548, 795.548};
            float CircleY[num] = {-451.592, -6.784, -443.43, -6.784};
            float CircleR[num] = {139.596, 132.888, 131.879, 127.58};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 186, 255);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 5000, 10000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "ratio", "and", 3, 5);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 300, 600);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1100, 1700);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row + 63.776).TupleConcat(hv_Row + 1374.511)).TupleConcat(hv_Row + 888.301)).TupleConcat(hv_Row + 1035.011),
                (((hv_Column - 251.24).TupleConcat(hv_Column + 1055.42)).TupleConcat(hv_Column - 1019.665)).TupleConcat(hv_Column - 242.24),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray10, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold10, KHminAreaA10, KHminAreaB10);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {193.514, -614.988, 947.831, 162.39};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold10, LXminArea10, LXmaxArea10);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
        // 拍照点11检测
        else if (index == 11) {
            // 焊点个数
            const int num = 10;
            // NG点个数，circlecount代表缺焊NG，rectanglecount代表连锡NG
            int circlecount = 0, rectanglecount = 0;
            // 从用户设置中获取参数
            auto SParray11 = getParamValue<QVariantList>("SParray11");
            float ROIminThreshold11 = getParamValue<float>("ROIminThreshold11");
            float KHminAreaA11 = getParamValue<float>("KHminAreaA11");
            float KHminAreaB11 = getParamValue<float>("KHminAreaB11");
            float LXminArea11 = getParamValue<float>("LXminArea11");
            float LXmaxArea11 = getParamValue<float>("LXmaxArea11");
            // Local iconic variables
            HObject ho_Regions0, ho_ConnectedRegions0, ho_SelectedRegions0, ho_RectangleNGRegion, ho_Circle, ho_UnionedCircles;
            // Local control variables
            HTuple hv_Number, hv_Area, hv_Row, hv_Column, hv_Row1, hv_Column1, hv_Row2, hv_Column2;
            HTuple hv_SParrayHandled, hv_CircleNG;
            // 圆形ROI坐标偏移量
            float CircleX[num] = {-408.13, -408.13, -408.13, -407.25, -409.22, -241.86, -242.67, -241.86, -243.83, -245.8};
            float CircleY[num] = {-404.846, -237.58, -68.34, 96.57, 261.87, -403.267, -232.34, -68.73, 96.57, 257.93};
            float CircleR[num] = {57.4031, 58.5367, 60.6432, 57.2282, 58.523, 59.7824, 60.2051, 56.376, 57.2282, 63.3923};
            // 读取原始图
            auto ho_Image = BImage2HalconImage(image);
            // 查找基准点
            Threshold(ho_Image, &ho_Regions0, 0, 20);
            Connection(ho_Regions0, &ho_ConnectedRegions0);
            SelectShape(ho_ConnectedRegions0, &ho_SelectedRegions0, "area", "and", 10000, 20000);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "row", "and", 1550, 1750);
            SelectShape(ho_SelectedRegions0, &ho_SelectedRegions0, "column", "and", 1200, 1500);
            CountObj(ho_SelectedRegions0, &hv_Number);
            if (hv_Number[0].I() != 1) {
                spdlog::error("未找到基准点");
                circlecount = -1;
                rectanglecount = -1;
                res->OK = false;
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            AreaCenter(ho_SelectedRegions0, &hv_Area, &hv_Row, &hv_Column);
            // 读取深度图
            auto ho_Imagedepth = BImage2HalconImage(depthimage);
            // 处理深度图
            GenCircle(&ho_Circle,
                (((hv_Row - 361.33).TupleConcat(hv_Row + 257.49)).TupleConcat(hv_Row + 474.01)).TupleConcat(hv_Row - 574.36),
                (((hv_Column - 1054.279).TupleConcat(hv_Column - 650.456)).TupleConcat(hv_Column + 518.63)).TupleConcat(hv_Column + 680.45),
                (((HTuple(15).Append(15)).Append(15)).Append(15)));
            Union1(ho_Circle, &ho_UnionedCircles);
            HObject ho_HandledImagedepth = DepthImageHandler(ho_Imagedepth, ho_UnionedCircles);
            // 处理焊点号
            hv_SParrayHandled = SParrayHandler(SParray11, num);
            if (hv_SParrayHandled[0] == -1) {
                res->OK = false;
                spdlog::error("焊点号超限");
                return AlgoOutputPtr(res);
            }
            // 处理圆形ROI，进行缺焊检测
            hv_CircleNG = CircleROIHandler(
                ho_HandledImagedepth, hv_Row, hv_Column, hv_SParrayHandled, CircleX, CircleY, CircleR, ROIminThreshold11, KHminAreaA11, KHminAreaB11);
            // 处理NG序列CircleNG[]
            for (int i = 0; i < num; i++) {
                if (hv_CircleNG[i] == 1) {
                    res->OK = false;
                    // 生成圆形roi的NG识别框
                    circlecount++;
                    int X = hv_Column[0].D() + CircleY[i] - CircleR[i];
                    int Y = hv_Row[0].D() + CircleX[i] - CircleR[i];
                    int W = CircleR[i] * 2;
                    int H = CircleR[i] * 2;
                    res->setCircleROI(circlecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(circlecount);
            // 矩形ROI坐标偏移量
            float RectangleXY[] = {-454.85, -443.67, -199.04, 316.7};
            // 处理矩形ROI，进行连锡检测
            ho_RectangleNGRegion =
                RectangleROIHandler(ho_HandledImagedepth, hv_Row, hv_Column, RectangleXY, ROIminThreshold11, LXminArea11, LXmaxArea11);
            if (ho_RectangleNGRegion.CountObj() > 0) {
                res->OK = false;
                SmallestRectangle1(ho_RectangleNGRegion, &hv_Row1, &hv_Column1, &hv_Row2, &hv_Column2);
                for (int i = 0; i < ho_RectangleNGRegion.CountObj(); i++) {
                    rectanglecount++;
                    int X = hv_Column1[i].D();
                    int Y = hv_Row1[i].D();
                    int W = hv_Column2[i].D() - hv_Column1[i].D();
                    int H = hv_Row2[i].D() - hv_Row1[i].D();
                    res->setRectangleROI(circlecount + rectanglecount, ROI::ROI(X, Y, W, H));
                }
            }
            // spdlog::info(rectanglecount);
            // 若无缺焊和连锡则为OK
            if (circlecount == 0 && rectanglecount == 0) {
                res->OK = true;
                return AlgoOutputPtr(res);
            } else {
                res->NGnum(circlecount, rectanglecount);
                return AlgoOutputPtr(res);
            }
            return AlgoOutputPtr(res);
        }
    } catch (HException e) {
        QString exception = QString::fromUtf8(e.ErrorMessage().ToUtf8());
        QStringList list = exception.split(":");
        spdlog::error("空焊、连锡检测异常:{}, code:{}", list.last().toStdString(), QString::number(e.ErrorCode()).toStdString());
        res->OK = false;
    } catch (...) {
        spdlog::error("空焊、连锡检测异常");
        res->OK = false;
    }
    return AlgoOutputPtr(res);
}

