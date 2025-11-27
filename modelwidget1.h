#ifndef MODELWIDGET1_H
#define MODELWIDGET1_H

#include <QWidget>
#include <QPainter>
#include <QVector>
#include <QMap>
#include <QMouseEvent>
#include <QWheelEvent>
#include <functional>
#include "modelmanager.h" // [关键] 引入 ModelManager 以使用统一数学内核

namespace Ui {
class ModelWidget1;
}

// =========================================================
// LogLogChartWidget: 绘图组件 (保持界面功能一致)
// =========================================================
class LogLogChartWidget : public QWidget
{
    Q_OBJECT

public:
    explicit LogLogChartWidget(QWidget *parent = nullptr);

    void setData(const QVector<double>& xData, const QVector<double>& yData1,
                 const QVector<double>& yData2, const QVector<double>& xData2 = QVector<double>());
    void clearData();
    void resetView();
    void autoFitData();
    void setShowOriginalData(bool show) { m_showOriginalData = show; update(); }

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;

private:
    void drawAxis(QPainter& painter, const QRect& plotRect);
    void drawData(QPainter& painter, const QRect& plotRect);
    void drawLegend(QPainter& painter, const QRect& plotRect);
    QPointF dataToPixel(double x, double y, const QRect& plotRect);

    // 原始数据
    QVector<double> m_xData;
    QVector<double> m_yData1;
    QVector<double> m_yData2;
    QVector<double> m_xData2;  // 压力导数的时间点

    QString m_title;

    // 显示范围
    double m_xMin, m_xMax, m_yMin, m_yMax;
    bool m_hasData;

    // 显示选项
    bool m_showOriginalData;

    // 鼠标交互
    bool m_isDragging;
    QPoint m_lastMousePos;
};

// =========================================================
// ModelWidget1: 界面逻辑类
// =========================================================
class ModelWidget1 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget1(QWidget *parent = nullptr);
    ~ModelWidget1();

private slots:
    void onCalculateClicked();
    void onResetParameters();
    void onExportResults();
    void onResetView();
    void onFitToData();

signals:
    void calculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private:
    // 核心计算流程 (现在只负责调用 ModelManager)
    void runCalculation();

private:
    Ui::ModelWidget1 *ui;
    LogLogChartWidget *m_chartWidget;

    // 缓存结果用于导出
    QVector<double> res_tD;
    QVector<double> res_pD;
    QVector<double> res_dpD;
    QVector<double> res_td_dpD; // 导数对应的时间点
};

#endif // MODELWIDGET1_H
