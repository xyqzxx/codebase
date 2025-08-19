#pragma once
#include "algo_helper.h"
#include "algo_type.h"
#include "convertor.hpp"

#define AOI_VERSION 1
class AOI : public QObject, public Algo {
    Q_OBJECT;
    Q_INTERFACES(Algo);

public:
    Q_INVOKABLE AOI() {
        // 11个拍照点位的检测参数
        int index = 1;
        for (; index <= 11; index++)
        {
            addParam("ROIminThreshold" + QString::number(index), 160).setDesc("焊点深度阈值下限").setRange(0, 1024).setVisible(true).setSavable(true);
            addParam("SParray" + QString::number(index), QVariantList{}).setDesc("使用空焊面积下限B的焊点").setVisible(true).setSavable(true);
            addParam("KHminAreaA" + QString::number(index), 0).setDesc("空焊面积下限A").setRange(0, 500000).setVisible(true).setSavable(true);
            addParam("KHminAreaB" + QString::number(index), 4000).setDesc("空焊面积下限B").setRange(0, 500000).setVisible(true).setSavable(true);
            addParam("LXminArea" + QString::number(index), 50000).setDesc("连锡面积下限").setRange(0, 500000).setVisible(true).setSavable(true);
            addParam("LXmaxArea" + QString::number(index), 150000).setDesc("连锡面积上限").setRange(0, 500000).setVisible(true).setSavable(true);
        }
        // 对参数分组
        addGroup("拍照点位一", 1)
            .store("ROIminThreshold1")
            .store(makeGroup("空焊检测参数").store("SParray1").store("KHminAreaA1").store("KHminAreaB1"))
            .store(makeGroup("连锡检测参数").store("LXminArea1").store("LXmaxArea1"));
        addGroup("拍照点位二", 2)
            .store("ROIminThreshold2")
            .store(makeGroup("空焊检测参数").store("SParray2").store("KHminAreaA2").store("KHminAreaB2"))
            .store(makeGroup("连锡检测参数").store("LXminArea2").store("LXmaxArea2"));
        addGroup("拍照点位三", 3)
            .store("ROIminThreshold3")
            .store(makeGroup("空焊检测参数").store("SParray3").store("KHminAreaA3").store("KHminAreaB3"))
            .store(makeGroup("连锡检测参数").store("LXminArea3").store("LXmaxArea3"));
        addGroup("拍照点位四", 4)
            .store("ROIminThreshold4")
            .store(makeGroup("空焊检测参数").store("SParray4").store("KHminAreaA4").store("KHminAreaB4"))
            .store(makeGroup("连锡检测参数").store("LXminArea4").store("LXmaxArea4"));
        addGroup("拍照点位五", 5)
            .store("ROIminThreshold5")
            .store(makeGroup("空焊检测参数").store("SParray5").store("KHminAreaA5").store("KHminAreaB5"))
            .store(makeGroup("连锡检测参数").store("LXminArea5").store("LXmaxArea5"));
        addGroup("拍照点位六", 6)
            .store("ROIminThreshold6")
            .store(makeGroup("空焊检测参数").store("SParray6").store("KHminAreaA6").store("KHminAreaB6"))
            .store(makeGroup("连锡检测参数").store("LXminArea6").store("LXmaxArea6"));
        addGroup("拍照点位七", 7)
            .store("ROIminThreshold7")
            .store(makeGroup("空焊检测参数").store("SParray7").store("KHminAreaA7").store("KHminAreaB7"))
            .store(makeGroup("连锡检测参数").store("LXminArea7").store("LXmaxArea7"));
        addGroup("拍照点位八", 8)
            .store("ROIminThreshold8")
            .store(makeGroup("空焊检测参数").store("SParray8").store("KHminAreaA8").store("KHminAreaB8"))
            .store(makeGroup("连锡检测参数").store("LXminArea8").store("LXmaxArea8"));
        addGroup("拍照点位九", 9)
            .store("ROIminThreshold9")
            .store(makeGroup("空焊检测参数").store("SParray9").store("KHminAreaA9").store("KHminAreaB9"))
            .store(makeGroup("连锡检测参数").store("LXminArea9").store("LXmaxArea9"));
        addGroup("拍照点位十", 10)
            .store("ROIminThreshold10")
            .store(makeGroup("空焊检测参数").store("SParray10").store("KHminAreaA10").store("KHminAreaB10"))
            .store(makeGroup("连锡检测参数").store("LXminArea10").store("LXmaxArea10"));
        addGroup("拍照点位十一", 11)
            .store("ROIminThreshold11")
            .store(makeGroup("空焊检测参数").store("SParray11").store("KHminAreaA11").store("KHminAreaB11"))
            .store(makeGroup("连锡检测参数").store("LXminArea11").store("LXmaxArea11"));
    }
    // 在算法库管理界面显示算法描述
    virtual QStringList getAlgoDesc() const override {
        return {"缺焊连锡检测算法"};
    }
    // 在算法库管理界面显示算法版本
    virtual int getAlgoVersion() const override {
        return AOI_VERSION;
    }

    virtual AlgoOutputPtr compute(const QVariantMap &execParam) override;
};
ALGO_DECLARE_EX("AOI", AOI);

class AOIOutput : public QObject, public AlgoOutput {
    Q_OBJECT
    Q_INTERFACES(AlgoOutput)

public:
    bool OK = false;

    virtual bool isOK() override {
        return OK;
    }
    virtual bool isNG() override {
        return !OK;
    }
    // 返回NG点个数，CircleNGnum代表缺焊NG，RectangleNGnum代表连锡NG
    // 都返回-1则为未找到基准点
    void NGnum(const int circlengnum, const int rectanglengnum) {
        output.insert("CircleNGnum", QVariant::fromValue(circlengnum));
        output.insert("RectangleNGnum", QVariant::fromValue(rectanglengnum));
    }
    void setBaseROI(const ROI& roi) {
        output.insert("BaseROI", QVariant::fromValue(roi));
    }
    // 设置显示在界面上的缺焊NG框
    void setCircleROI(const int index, const ROI &roi) {
        output.insert(QString::number(index), QVariant::fromValue(roi));
    }
    // 设置显示在界面上的连锡NG框
    void setRectangleROI(const int index, const ROI &roi) {
        output.insert(QString::number(index), QVariant::fromValue(roi));
    }
    
    void setResults(const QVariantList results) {
        output.insert("results", QVariant::fromValue(results));
    }

};