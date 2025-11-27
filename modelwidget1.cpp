#include "modelwidget1.h"
#include "ui_modelwidget1.h"
#include <cmath>
#include <algorithm>
#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QDateTime>
#include <QCoreApplication>

// =========================================================
// LogLogChartWidget 实现 (保持原有绘图逻辑不变)
// =========================================================

LogLogChartWidget::LogLogChartWidget(QWidget *parent)
    : QWidget(parent)
    , m_xMin(1e1), m_xMax(1e15)
    , m_yMin(1e-2), m_yMax(1e1)
    , m_hasData(false)
    , m_showOriginalData(true)
    , m_isDragging(false)
{
    setMinimumSize(600, 400);
    setStyleSheet("QWidget { background-color: white; border: 1px solid gray; }");
    setMouseTracking(true);
}

void LogLogChartWidget::setData(const QVector<double>& xData, const QVector<double>& yData1,
                                const QVector<double>& yData2, const QVector<double>& xData2)
{
    m_xData = xData;
    m_yData1 = yData1;
    m_yData2 = yData2;
    m_xData2 = xData2.isEmpty() ? xData : xData2;
    m_hasData = !xData.isEmpty() && !yData1.isEmpty();

    if (m_hasData) {
        autoFitData();
    }
    update();
}

void LogLogChartWidget::clearData()
{
    m_hasData = false;
    m_xData.clear(); m_yData1.clear(); m_yData2.clear(); m_xData2.clear();
    resetView();
    update();
}

void LogLogChartWidget::resetView()
{
    if (m_hasData) {
        autoFitData();
    } else {
        m_xMin = 1e1; m_xMax = 1e15;
        m_yMin = 1e-2; m_yMax = 1e1;
    }
    update();
}

void LogLogChartWidget::autoFitData()
{
    if (!m_hasData) return;

    QVector<double> allX, allY;
    for (double x : m_xData) if (x > 0 && std::isfinite(x)) allX.append(x);
    for (double x : m_xData2) if (x > 0 && std::isfinite(x)) allX.append(x);
    for (double y : m_yData1) if (y > 0 && std::isfinite(y)) allY.append(y);
    for (double y : m_yData2) if (y > 0 && std::isfinite(y)) allY.append(y);

    if (allX.isEmpty() || allY.isEmpty()) return;

    auto [xMinIt, xMaxIt] = std::minmax_element(allX.begin(), allX.end());
    auto [yMinIt, yMaxIt] = std::minmax_element(allY.begin(), allY.end());

    double xMin = *xMinIt; double xMax = *xMaxIt;
    double yMin = *yMinIt; double yMax = *yMaxIt;

    double logXMin = log10(xMin); double logXMax = log10(xMax);
    double logYMin = log10(yMin); double logYMax = log10(yMax);
    double xMargin = (logXMax - logXMin) * 0.05;
    double yMargin = (logYMax - logYMin) * 0.05;

    m_xMin = pow(10.0, logXMin - xMargin);
    m_xMax = pow(10.0, logXMax + xMargin);
    m_yMin = pow(10.0, logYMin - yMargin);
    m_yMax = pow(10.0, logYMax + yMargin);
}

void LogLogChartWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_isDragging = true; m_lastMousePos = event->pos();
        setCursor(Qt::ClosedHandCursor);
    }
}

void LogLogChartWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (m_isDragging) {
        QPoint delta = event->pos() - m_lastMousePos;
        m_lastMousePos = event->pos();
        QRect plotRect = rect().adjusted(80, 50, -50, -80);
        double logXRange = log10(m_xMax) - log10(m_xMin);
        double logYRange = log10(m_yMax) - log10(m_yMin);
        double deltaLogX = -delta.x() * logXRange / plotRect.width();
        double deltaLogY = delta.y() * logYRange / plotRect.height();
        m_xMin = pow(10.0, log10(m_xMin) + deltaLogX);
        m_xMax = pow(10.0, log10(m_xMax) + deltaLogX);
        m_yMin = pow(10.0, log10(m_yMin) + deltaLogY);
        m_yMax = pow(10.0, log10(m_yMax) + deltaLogY);
        update();
    } else {
        setCursor(Qt::OpenHandCursor);
    }
}

void LogLogChartWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_isDragging = false; setCursor(Qt::OpenHandCursor);
    }
}

void LogLogChartWidget::wheelEvent(QWheelEvent *event)
{
    const double zoomInFactor = 1.15;
    const double zoomOutFactor = 1.0 / zoomInFactor;
    double zoomFactor = (event->angleDelta().y() > 0) ? zoomInFactor : zoomOutFactor;

    QRect plotRect = rect().adjusted(80, 50, -50, -80);
    QPointF mousePos = event->position();
    double relativeX = (mousePos.x() - plotRect.left()) / plotRect.width();
    double relativeY = (plotRect.bottom() - mousePos.y()) / plotRect.height();
    relativeX = qBound(0.0, relativeX, 1.0);
    relativeY = qBound(0.0, relativeY, 1.0);

    double logXRange = log10(m_xMax) - log10(m_xMin);
    double logYRange = log10(m_yMax) - log10(m_yMin);
    double mouseLogX = log10(m_xMin) + relativeX * logXRange;
    double mouseLogY = log10(m_yMin) + relativeY * logYRange;
    double newLogXRange = logXRange / zoomFactor;
    double newLogYRange = logYRange / zoomFactor;

    m_xMin = pow(10.0, mouseLogX - relativeX * newLogXRange);
    m_xMax = pow(10.0, mouseLogX + (1.0 - relativeX) * newLogXRange);
    m_yMin = pow(10.0, mouseLogY - relativeY * newLogYRange);
    m_yMax = pow(10.0, mouseLogY + (1.0 - relativeY) * newLogYRange);
    update();
    event->accept();
}

void LogLogChartWidget::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    QRect plotRect = rect().adjusted(80, 50, -50, -80);
    painter.fillRect(rect(), Qt::white);
    painter.fillRect(plotRect, QColor(250, 250, 250));

    if (!m_title.isEmpty()) {
        painter.setFont(QFont("Arial", 14, QFont::Bold));
        painter.setPen(Qt::black);
        painter.drawText(rect().adjusted(0, 10, 0, 0), Qt::AlignHCenter | Qt::AlignTop, m_title);
    }
    drawAxis(painter, plotRect);
    if (m_hasData) {
        drawData(painter, plotRect);
        drawLegend(painter, plotRect);
    } else {
        painter.setFont(QFont("Arial", 12));
        painter.setPen(Qt::gray);
        painter.drawText(plotRect, Qt::AlignCenter, "点击'开始计算'生成压力分析图表");
    }
}

void LogLogChartWidget::drawAxis(QPainter& painter, const QRect& plotRect)
{
    painter.setPen(QPen(Qt::black, 2));
    painter.drawRect(plotRect);
    painter.setPen(QPen(Qt::lightGray, 1));
    painter.setFont(QFont("Arial", 9));

    double logXMin = log10(m_xMin);
    double logXMax = log10(m_xMax);
    int startExp = (int)floor(logXMin);
    int endExp = (int)ceil(logXMax);

    for (int exp = startExp; exp <= endExp; exp++) {
        double x = pow(10.0, exp);
        if (x >= m_xMin && x <= m_xMax) {
            QPointF point = dataToPixel(x, m_yMin, plotRect);
            if (point.x() >= plotRect.left() && point.x() <= plotRect.right()) {
                painter.drawLine(point.x(), plotRect.bottom(), point.x(), plotRect.top());
                painter.setPen(Qt::black);
                QString label = QString("1e%1").arg(exp);
                QRect labelRect(point.x() - 25, plotRect.bottom() + 5, 50, 20);
                painter.drawText(labelRect, Qt::AlignCenter, label);
                painter.setPen(Qt::lightGray);
            }
        }
    }

    double logYMin = log10(m_yMin);
    double logYMax = log10(m_yMax);
    startExp = (int)floor(logYMin);
    endExp = (int)ceil(logYMax);

    for (int exp = startExp; exp <= endExp; exp++) {
        double y = pow(10.0, exp);
        if (y >= m_yMin && y <= m_yMax) {
            QPointF point = dataToPixel(m_xMin, y, plotRect);
            if (point.y() >= plotRect.top() && point.y() <= plotRect.bottom()) {
                painter.drawLine(plotRect.left(), point.y(), plotRect.right(), point.y());
                painter.setPen(Qt::black);
                QString label = QString("1e%1").arg(exp);
                QRect labelRect(plotRect.left() - 50, point.y() - 10, 45, 20);
                painter.drawText(labelRect, Qt::AlignRight | Qt::AlignVCenter, label);
                painter.setPen(Qt::lightGray);
            }
        }
    }

    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 11, QFont::Bold));
    painter.drawText(plotRect.center().x() - 30, plotRect.bottom() + 40, "tD/CD");
    painter.save();
    painter.translate(plotRect.left() - 60, plotRect.center().y());
    painter.rotate(-90);
    painter.drawText(-50, 0, "PD & dPD");
    painter.restore();
}

void LogLogChartWidget::drawData(QPainter& painter, const QRect& plotRect)
{
    if (!m_xData.isEmpty() && !m_yData1.isEmpty()) {
        painter.setPen(QPen(Qt::red, 2));
        QVector<QPointF> points;
        for (int i = 0; i < qMin(m_xData.size(), m_yData1.size()); ++i) {
            if (m_xData[i] > 0 && m_yData1[i] > 0) {
                QPointF point = dataToPixel(m_xData[i], m_yData1[i], plotRect);
                points.append(point); // 简单处理，不裁剪，交由drawLine处理
            }
        }
        for (int i = 1; i < points.size(); ++i) {
            if(plotRect.contains(points[i-1].toPoint()) || plotRect.contains(points[i].toPoint()))
                painter.drawLine(points[i-1], points[i]);
        }
        if (m_showOriginalData) {
            painter.setPen(QPen(Qt::red, 1)); painter.setBrush(QBrush(Qt::red));
            for (const QPointF& point : points) if(plotRect.contains(point.toPoint())) painter.drawEllipse(point, 3, 3);
        }
    }

    if (!m_xData2.isEmpty() && !m_yData2.isEmpty()) {
        painter.setPen(QPen(Qt::blue, 2));
        QVector<QPointF> points;
        for (int i = 0; i < qMin(m_xData2.size(), m_yData2.size()); ++i) {
            if (m_xData2[i] > 0 && m_yData2[i] > 0) {
                QPointF point = dataToPixel(m_xData2[i], m_yData2[i], plotRect);
                points.append(point);
            }
        }
        for (int i = 1; i < points.size(); ++i) {
            if(plotRect.contains(points[i-1].toPoint()) || plotRect.contains(points[i].toPoint()))
                painter.drawLine(points[i-1], points[i]);
        }
        if (m_showOriginalData) {
            painter.setPen(QPen(Qt::blue, 1)); painter.setBrush(QBrush(Qt::blue));
            for (const QPointF& point : points) if(plotRect.contains(point.toPoint())) painter.drawEllipse(point, 3, 3);
        }
    }
}

void LogLogChartWidget::drawLegend(QPainter& painter, const QRect& plotRect)
{
    painter.setFont(QFont("Arial", 10));
    int legendX = plotRect.right() - 160;
    int legendY = plotRect.top() + 20;
    int lineSpacing = 20;

    painter.setPen(QPen(Qt::red, 2));
    painter.drawLine(legendX, legendY, legendX + 30, legendY);
    painter.setPen(Qt::black);
    painter.drawText(legendX + 35, legendY + 5, "压力");

    if (!m_yData2.isEmpty()) {
        legendY += lineSpacing;
        painter.setPen(QPen(Qt::blue, 2));
        painter.drawLine(legendX, legendY, legendX + 30, legendY);
        painter.setPen(Qt::black);
        painter.drawText(legendX + 35, legendY + 5, "压力导数");
    }
}

QPointF LogLogChartWidget::dataToPixel(double x, double y, const QRect& plotRect)
{
    x = qMax(x, 1e-20); y = qMax(y, 1e-20);
    double logX = log10(x); double logY = log10(y);
    double logXMin = log10(qMax(m_xMin, 1e-20));
    double logXMax = log10(qMax(m_xMax, 1e-20));
    double logYMin = log10(qMax(m_yMin, 1e-20));
    double logYMax = log10(qMax(m_yMax, 1e-20));

    double pixelX = plotRect.left() + (logX - logXMin) / (logXMax - logXMin) * plotRect.width();
    double pixelY = plotRect.bottom() - (logY - logYMin) / (logYMax - logYMin) * plotRect.height();
    return QPointF(pixelX, pixelY);
}

// =========================================================
// ModelWidget1 主逻辑 (重构后调用 ModelManager)
// =========================================================

ModelWidget1::ModelWidget1(QWidget *parent) :
    QWidget(parent), ui(new Ui::ModelWidget1)
{
    ui->setupUi(this);
    m_chartWidget = new LogLogChartWidget(this);
    QVBoxLayout *layout = new QVBoxLayout(ui->chartTab);
    layout->addWidget(m_chartWidget);
    layout->setContentsMargins(0,0,0,0);

    connect(ui->calculateButton, &QPushButton::clicked, this, &ModelWidget1::onCalculateClicked);
    connect(ui->resetButton, &QPushButton::clicked, this, &ModelWidget1::onResetParameters);
    connect(ui->exportButton, &QPushButton::clicked, this, &ModelWidget1::onExportResults);
    connect(ui->resetViewButton, &QPushButton::clicked, this, &ModelWidget1::onResetView);
    connect(ui->fitToDataButton, &QPushButton::clicked, this, &ModelWidget1::onFitToData);

    onResetParameters();
}

ModelWidget1::~ModelWidget1() { delete ui; }

void ModelWidget1::onResetParameters() {
    ui->omegaSpinBox->setValue(0.05);
    ui->sSpinBox->setValue(1.0);
    ui->cDSpinBox->setValue(1e-8);
    ui->lambdaSpinBox->setValue(0.1); // MATLAB: 1e-1
    ui->mfSpinBox->setValue(3);
    ui->nfSpinBox->setValue(5);
    ui->xfSpinBox->setValue(40);
    ui->yySpinBox->setValue(70);
    ui->ySpinBox->setValue(1000);
    // [关键修正] 统一设为 4
    ui->nSpinBox->setValue(4);
}

void ModelWidget1::onResetView() { m_chartWidget->resetView(); }
void ModelWidget1::onFitToData() { m_chartWidget->autoFitData(); }

void ModelWidget1::onCalculateClicked() {
    ui->calculateButton->setEnabled(false);
    ui->calculateButton->setText("计算中...");
    QCoreApplication::processEvents();

    runCalculation();

    ui->calculateButton->setEnabled(true);
    ui->calculateButton->setText("开始计算");
    ui->exportButton->setEnabled(true);
    ui->resetViewButton->setEnabled(true);
    ui->fitToDataButton->setEnabled(true);
    ui->tabWidget->setCurrentIndex(0);
}

void ModelWidget1::runCalculation() {
    // 1. 收集参数到 Map
    QMap<QString, double> params;
    params.insert("omega", ui->omegaSpinBox->value());
    params.insert("S", ui->sSpinBox->value());
    params.insert("cD", ui->cDSpinBox->value());
    params.insert("lambda", ui->lambdaSpinBox->value());
    params.insert("mf", ui->mfSpinBox->value());
    params.insert("nf", ui->nfSpinBox->value());
    params.insert("Xf", ui->xfSpinBox->value());
    params.insert("yy", ui->yySpinBox->value());
    params.insert("y", ui->ySpinBox->value());
    params.insert("N", ui->nSpinBox->value());

    // 2. [关键] 调用 ModelManager 的静态内核，确保算法唯一性
    // 注意：ModelManager::calculateModel1 是静态的，无需实例
    ModelCurveData result = ModelManager::calculateModel1(params);

    // 3. 解析结果
    res_tD = std::get<0>(result);  // tD
    res_pD = std::get<1>(result);  // PD
    QVector<double> dpVec = std::get<2>(result); // dPD (对应 td(2:end))

    // 4. 数据对齐与绘图
    // 绘图使用 tD/CD
    double cD = params["cD"];
    QVector<double> plot_tD_CD;
    for(double val : res_tD) plot_tD_CD.append(val / cD);

    // 构造导数对应的时间点 (td(2:end))
    // ModelManager 返回的 dpVec 长度比 tD 少 1
    res_dpD = dpVec;
    res_td_dpD.clear();
    if (plot_tD_CD.size() > 1) {
        // 从第2个点开始对应
        for(int i=1; i<plot_tD_CD.size(); ++i) {
            res_td_dpD.append(plot_tD_CD[i]);
        }
    }
    // 安全截断防止越界
    if (res_dpD.size() > res_td_dpD.size()) res_dpD.resize(res_td_dpD.size());

    // 5. 更新图表
    m_chartWidget->setData(plot_tD_CD, res_pD, res_dpD, res_td_dpD);

    // 6. 更新文本结果
    QString resultText = QString("计算完成\n点数: %1\n").arg(res_tD.size());
    resultText += "tD/CD\t\tPD\t\tdPD\n";
    for(int i=0; i<10 && i<res_pD.size(); ++i) {
        double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
        resultText += QString("%1\t%2\t%3\n").arg(plot_tD_CD[i],0,'e',4).arg(res_pD[i],0,'e',4).arg(dp,0,'e',4);
    }
    ui->resultTextEdit->setText(resultText);

    QMap<QString, double> rMap; rMap["PointCount"] = res_tD.size();
    emit calculationCompleted("MFHW_DualPorosity", rMap);
}

void ModelWidget1::onExportResults() {
    if (res_tD.isEmpty()) return;
    QString path = QFileDialog::getSaveFileName(this, "导出CSV", "", "CSV Files (*.csv)");
    if (path.isEmpty()) return;
    QFile f(path);
    if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&f);
        double cD = ui->cDSpinBox->value();
        out << "tD,tD/CD,PD,dPD\n";
        for (int i = 0; i < res_tD.size(); ++i) {
            double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
            out << res_tD[i] << "," << res_tD[i]/cD << "," << res_pD[i] << "," << dp << "\n";
        }
        f.close();
        QMessageBox::information(this, "导出成功", "文件已保存");
    }
}
