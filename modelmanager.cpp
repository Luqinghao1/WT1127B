#include "modelmanager.h"
#include "modelwidget1.h"
#include "modelwidget2.h"
#include "modelwidget3.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QDebug>
#include <cmath>
#include <QtMath>
#include <algorithm>

// ... [构造函数和界面初始化代码保持不变] ...

ModelManager::ModelManager(QWidget* parent)
    : QObject(parent), m_mainWidget(nullptr), m_modelTypeCombo(nullptr), m_modelStack(nullptr)
    , m_modelWidget1(nullptr), m_modelWidget2(nullptr), m_modelWidget3(nullptr), m_currentModelType(InfiniteConductive) {}

ModelManager::~ModelManager() {}

void ModelManager::initializeModels(QWidget* parentWidget) {
    if (!parentWidget) return;
    createMainWidget();
    setupModelSelection();
    m_modelStack = new QStackedWidget(m_mainWidget);
    m_modelWidget1 = new ModelWidget1(m_modelStack);
    m_modelWidget2 = new ModelWidget2(m_modelStack);
    m_modelWidget3 = new ModelWidget3(m_modelStack);
    m_modelStack->addWidget(m_modelWidget1);
    m_modelStack->addWidget(m_modelWidget2);
    m_modelStack->addWidget(m_modelWidget3);
    m_mainWidget->layout()->addWidget(m_modelStack);
    connectModelSignals();
    switchToModel(InfiniteConductive);
    if (parentWidget->layout()) parentWidget->layout()->addWidget(m_mainWidget);
    else {
        QVBoxLayout* layout = new QVBoxLayout(parentWidget);
        layout->addWidget(m_mainWidget);
        parentWidget->setLayout(layout);
    }
}

void ModelManager::createMainWidget() {
    m_mainWidget = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_mainWidget);
    mainLayout->setContentsMargins(10, 5, 10, 10);
    mainLayout->setSpacing(0);
    m_mainWidget->setLayout(mainLayout);
}

void ModelManager::setupModelSelection() {
    if (!m_mainWidget) return;
    QGroupBox* selectionGroup = new QGroupBox("模型类型选择", m_mainWidget);
    selectionGroup->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    QHBoxLayout* selectionLayout = new QHBoxLayout(selectionGroup);
    selectionLayout->setContentsMargins(9, 9, 9, 9);
    selectionLayout->setSpacing(6);
    QLabel* typeLabel = new QLabel("模型类型:", selectionGroup);
    typeLabel->setMinimumWidth(100);
    m_modelTypeCombo = new QComboBox(selectionGroup);
    m_modelTypeCombo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    m_modelTypeCombo->setMinimumWidth(200);
    m_modelTypeCombo->addItem(getModelTypeName(InfiniteConductive));
    m_modelTypeCombo->addItem(getModelTypeName(FiniteConductive));
    m_modelTypeCombo->addItem(getModelTypeName(SegmentedMultiCluster));
    m_modelTypeCombo->setStyleSheet("color: black;");
    typeLabel->setStyleSheet("color: black;");
    selectionGroup->setStyleSheet("QGroupBox { color: black; font-weight: bold; }");
    connect(m_modelTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &ModelManager::onModelTypeSelectionChanged);
    selectionLayout->addWidget(typeLabel);
    selectionLayout->addWidget(m_modelTypeCombo);
    QVBoxLayout* mainLayout = qobject_cast<QVBoxLayout*>(m_mainWidget->layout());
    if (mainLayout) { mainLayout->addWidget(selectionGroup); mainLayout->setStretchFactor(selectionGroup, 0); }
}

void ModelManager::connectModelSignals() {
    if (m_modelWidget1) connect(m_modelWidget1, &ModelWidget1::calculationCompleted, this, &ModelManager::onModel1CalculationCompleted);
    if (m_modelWidget2) connect(m_modelWidget2, &ModelWidget2::calculationCompleted, this, &ModelManager::onModel2CalculationCompleted);
    if (m_modelWidget3) connect(m_modelWidget3, &ModelWidget3::calculationCompleted, this, &ModelManager::onModel3CalculationCompleted);
}

void ModelManager::switchToModel(ModelType modelType) {
    if (!m_modelStack) return;
    if (modelType == m_currentModelType) return;
    ModelType old = m_currentModelType;
    m_currentModelType = modelType;
    m_modelStack->setCurrentIndex((int)modelType);
    if (m_modelTypeCombo) m_modelTypeCombo->setCurrentIndex((int)modelType);
    emit modelSwitched(modelType, old);
}

void ModelManager::onModelTypeSelectionChanged(int index) { switchToModel((ModelType)index); }

QString ModelManager::getModelTypeName(ModelType type) {
    switch (type) {
    case InfiniteConductive: return "无限导流双重孔隙介质页岩油藏渗流模型";
    case FiniteConductive: return "有限导流双重孔隙介质页岩油藏渗流模型";
    case SegmentedMultiCluster: return "分段多簇压裂水平井双重孔隙介质页岩油藏渗流模型";
    default: return "未知模型";
    }
}

QStringList ModelManager::getAvailableModelTypes() {
    return { getModelTypeName(InfiniteConductive), getModelTypeName(FiniteConductive), getModelTypeName(SegmentedMultiCluster) };
}

void ModelManager::onModel1CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }
void ModelManager::onModel2CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }
void ModelManager::onModel3CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }

// ==================================================================================
//  计算引擎参数定义
// ==================================================================================

QMap<QString, double> ModelManager::getDefaultParameters(ModelType type)
{
    QMap<QString, double> p;
    p.insert("k", 1.0); p.insert("S", 0.1); p.insert("cD", 1e-8);
    p.insert("N", 4.0); // 默认 Stehfest N = 4

    if (type == InfiniteConductive) {
        p.insert("omega", 0.05); p.insert("lambda", 1e-2); p.insert("mf", 3.0);
        p.insert("nf", 5.0); p.insert("Xf", 40.0); p.insert("yy", 70.0); p.insert("y", 1000.0);
    } else if (type == FiniteConductive) {
        p.insert("omega", 0.0155); p.insert("lambda", 0.083); p.insert("mf", 3.0);
        p.insert("nf", 5.0); p.insert("Xf", 193.0); p.insert("yy", 295.0);
        p.insert("y", 2758.0); p.insert("CFD", 0.9); p.insert("kpd", 0.04);
        p.insert("S", 0.81); p.insert("cD", 8.08e-8);
    } else if (type == SegmentedMultiCluster) {
        p.insert("omega1", 0.05); p.insert("omega2", 0.05); p.insert("lambda1", 1e-1);
        p.insert("lambda2", 1e-1); p.insert("mf1", 2.0); p.insert("mf2", 2.0);
        p.insert("nf", 5.0); p.insert("Xf1", 40.0); p.insert("Xf2", 40.0);
        p.insert("yy1", 70.0); p.insert("yy2", 70.0); p.insert("y", 800.0);
        p.insert("CFD1", 0.4); p.insert("CFD2", 0.4); p.insert("kpd", 0.045);
    }
    return p;
}

ModelCurveData ModelManager::calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params)
{
    switch(type) {
    case InfiniteConductive: return calculateModel1(params);
    case FiniteConductive: return calculateModel2(params);
    case SegmentedMultiCluster: return calculateModel3(params);
    default: return calculateModel1(params);
    }
}

// ==================================================================================
//  数学计算实现 (静态函数)
// ==================================================================================

ModelCurveData ModelManager::calculateModel1(const QMap<QString, double>& params)
{
    int numPoints = 100;
    int N = (int)params.value("N", 4);
    if (N % 2 != 0) N = 4;

    QVector<double> tD; tD.reserve(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        double exponent = -7.0 + 13.0 * i / (numPoints - 1);
        tD.append(pow(10.0, exponent));
    }

    QVector<double> pd(numPoints, 0.0);
    double ln2 = log(2.0);

    for (int k = 0; k < numPoints; ++k) {
        if (tD[k] <= 0) continue;
        for (int m = 1; m <= N; ++m) {
            double s = m * ln2 / tD[k];
            double L = flaplace1(s, params);
            double coeff = stefestCoefficient(m, N);
            pd[k] += coeff * ln2 * L / tD[k];
        }
    }

    QVector<double> dpVec, tdpVec;
    double cD = params.value("cD", 1e-8);
    computePressureDerivative(tD, pd, cD, dpVec, tdpVec);
    return std::make_tuple(tD, pd, dpVec);
}

ModelCurveData ModelManager::calculateModel2(const QMap<QString, double>& params)
{
    int numPoints = 100;
    int N = (int)params.value("N", 4);
    double kpd = params.value("kpd", 0.04);

    QVector<double> tD;
    for (int i = 0; i < numPoints; ++i) {
        double exponent = -11.0 + 16.0 * i / (numPoints - 1);
        tD.append(pow(10.0, exponent));
    }

    QVector<double> pd(numPoints, 0.0);
    double ln2 = log(2.0);

    for (int k = 0; k < numPoints; ++k) {
        for (int m = 1; m <= N; ++m) {
            double s = m * ln2 / tD[k];
            double L = flaplace2(s, params);
            pd[k] += stefestCoefficient(m, N) * ln2 * L / tD[k];
        }
    }

    QVector<double> pd1(numPoints);
    for(int i=0; i<numPoints; ++i) pd1[i] = -1.0 / kpd * log(1.0 - kpd * pd[i]);

    QVector<double> dpVec, tdpVec;
    double cD = params.value("cD", 1e-8);
    computePressureDerivative(tD, pd1, cD, dpVec, tdpVec);
    return std::make_tuple(tD, pd1, dpVec);
}

ModelCurveData ModelManager::calculateModel3(const QMap<QString, double>& params) { return calculateModel1(params); }

void ModelManager::computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                             double cD, QVector<double>& dpd, QVector<double>& td_dpd)
{
    dpd.clear(); td_dpd.clear();
    if (tD.size() < 2) return;

    QVector<double> plot_tD_CD;
    for(double val : tD) plot_tD_CD.append(val / cD);

    for (int k = 0; k < tD.size() - 1; ++k) {
        double t_next = plot_tD_CD[k+1];
        double t_curr = plot_tD_CD[k];
        double p_next = pd[k+1];
        double p_curr = pd[k];

        if (std::abs(t_next - t_curr) > 1e-15) {
            double val = t_next * (p_next - p_curr) / (t_next - t_curr);
            dpd.append(val);
            td_dpd.append(t_next);
        }
    }
}

// ==================================================================================
//  数学基础 (静态实现)
// ==================================================================================

static const double GL15_X[] = { 0.0, 0.2011940939974345, 0.3941513470775634, 0.5709721726085388, 0.7244177313601701, 0.8482065834104272, 0.9372985251687639, 0.9879925180204854 };
static const double GL15_W[] = { 0.2025782419255613, 0.1984314853271116, 0.1861610000155622, 0.1662692058169939, 0.1395706779049514, 0.1071592204671719, 0.0703660474881081, 0.0307532419961173 };

double ModelManager::gauss15(std::function<double(double)> f, double a, double b) {
    double halfLen = 0.5 * (b - a);
    double center = 0.5 * (a + b);
    double sum = GL15_W[0] * f(center);
    for (int i = 1; i < 8; ++i) {
        double dx = halfLen * GL15_X[i];
        sum += GL15_W[i] * (f(center - dx) + f(center + dx));
    }
    return sum * halfLen;
}

double ModelManager::adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth) {
    double c = (a + b) / 2.0;
    double v1 = gauss15(f, a, b);
    double v2 = gauss15(f, a, c) + gauss15(f, c, b);
    if (depth >= maxDepth) return v2;
    if (std::abs(v1 - v2) < 1e-9 * std::abs(v2) + eps) return v2; // 稍微放宽精度以提升拟合速度
    return adaptiveGauss(f, a, c, eps/2, depth+1, maxDepth) + adaptiveGauss(f, c, b, eps/2, depth+1, maxDepth);
}

double ModelManager::integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double a, double b)
{
    double sqrt_fz = sqrt(fz);
    double dist_y_sq = (YDkv - yDij) * (YDkv - yDij);

    auto func = [=](double xwD) -> double {
        double dist_x = XDkv - xwD;
        double arg = sqrt(dist_x * dist_x + dist_y_sq) * sqrt_fz;
        return besselK0(arg);
    };

    // [优化] 远场判断：如果距离较远，直接使用固定高斯积分，避免昂贵的自适应递归
    double segmentLen = b - a;
    double distToCenter = sqrt(pow(XDkv - (a+b)/2.0, 2) + dist_y_sq);
    // 阈值设为 2.0 倍线段长，此时 K0 函数变化足够平缓
    if (distToCenter > 2.0 * segmentLen) {
        return gauss15(func, a, (a+b)/2.0) + gauss15(func, (a+b)/2.0, b);
    }

    // [优化] 奇异性处理：只有在真正重合时才深度递归
    if (dist_y_sq < 1e-16) {
        if (XDkv > a + 1e-9 && XDkv < b - 1e-9) {
            // 奇异点在区间内，必须拆分，限制深度 8 足够拟合使用
            return adaptiveGauss(func, a, XDkv, 1e-6, 0, 8) + adaptiveGauss(func, XDkv, b, 1e-6, 0, 8);
        }
    }
    // 近场非奇异，中等精度
    return adaptiveGauss(func, a, b, 1e-6, 0, 8);
}

double ModelManager::besselK0(double x)
{
    if (x <= 1e-15) return 50.0;
#if __cplusplus >= 201703L
    return std::cyl_bessel_k(0, x);
#else
    // 快速近似实现
    if (x <= 2.0) {
        double y = x * x / 4.0;
        return (-log(x / 2.0) * std::cyl_bessel_i(0, x)) +
               (0.42278420 + y * (0.23069756 + y * (0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.00000740)))));
    } else {
        double y = 2.0 / x;
        return (exp(-x) / sqrt(x)) * (1.25331414 + y * (-0.07832358 + y * (0.02189568 + y * (-0.01062446 + y * (0.00587872 + y * (-0.00251540 + y * 0.00053208))))));
    }
#endif
}

QVector<double> ModelManager::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b)
{
    if(A.isEmpty() || b.isEmpty() || A.size() != b.size()) return QVector<double>();
    int n = b.size();
    QVector<QVector<double>> M = A;
    for(int i=0; i<n; ++i) M[i].append(b[i]);

    for (int k = 0; k < n - 1; ++k) {
        if (std::abs(M[k][k]) < 1e-12) continue;
        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j < n + 1; ++j) M[i][j] -= factor * M[k][j];
        }
    }
    QVector<double> x(n);
    if (std::abs(M[n-1][n-1]) > 1e-12) x[n-1] = M[n-1][n] / M[n-1][n-1];
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += M[i][j] * x[j];
        if (std::abs(M[i][i]) > 1e-12) x[i] = (M[i][n] - sum) / M[i][i];
    }
    return x;
}

double ModelManager::stefestCoefficient(int i, int N)
{
    double sum = 0.0;
    int start = (i + 1) / 2;
    int end = std::min(i, N / 2);
    for (int k = start; k <= end; ++k) {
        double num = pow(k, N/2.0) * factorial(2*k);
        double den = factorial(N/2-k) * factorial(k) * factorial(k-1) * factorial(i-k) * factorial(2*k-i);
        if(den > 0) sum += num/den;
    }
    return ((N/2+i)%2==0?1:-1)*sum;
}

double ModelManager::factorial(int n) {
    if (n<=1) return 1.0;
    static QVector<double> c={1.0,1.0};
    if (n < c.size()) return c[n];
    for(int i=c.size();i<=n;++i) c.append(c.last()*i);
    return c[n];
}

double ModelManager::flaplace1(double z, const QMap<QString, double>& p) {
    int mf = (int)p.value("mf", 3);
    int nf = (int)p.value("nf", 5);
    double omega = p.value("omega", 0.05);
    double lambda = p.value("lambda", 1e-2);
    double Xf = p.value("Xf", 40.0);
    double yy = p.value("yy", 70.0);
    double y = p.value("y", 1000.0);
    double S = p.value("S", 0.1);
    double cD = p.value("cD", 1e-8);

    if (mf <= 0 || nf <= 0) return 0.0;

    int sz = mf * 2 * nf;
    QVector<QVector<double>> E(sz + 1, QVector<double>(sz + 1, 0.0));
    QVector<double> F(sz + 1, 0.0);

    for(int i=1; i<=mf; ++i) {
        for(int j=1; j<=2*nf; ++j) {
            for(int k=1; k<=mf; ++k) {
                for(int v=1; v<=2*nf; ++v) {
                    E[(i-1)*2*nf+j-1][(k-1)*2*nf+v-1] = e_function(z,i,j,k,v,mf,nf,omega,lambda,Xf,yy,y);
                }
            }
            E[(i-1)*2*nf+j-1][sz] = -1.0;
            E[sz][(i-1)*2*nf+j-1] = f_function(j,nf,Xf,y);
        }
    }
    F[sz] = 1.0/z;

    QVector<double> res = solveLinearSystem(E, F);
    if (res.size() <= sz) return 0.0;

    double pwd = z * res[sz] + S;
    return pwd / (z * (1.0 + z * cD * pwd));
}

double ModelManager::flaplace2(double z, const QMap<QString, double>& p) {
    int mf = (int)p.value("mf", 3);
    int nf = (int)p.value("nf", 5);
    double omega = p.value("omega", 0.0155);
    double lambda = p.value("lambda", 0.083);
    double Xf = p.value("Xf", 193.0);
    double yy = p.value("yy", 295.0);
    double y = p.value("y", 2758.0);
    double S = p.value("S", 0.81);
    double cD = p.value("cD", 8.08e-8);
    double CFD = p.value("CFD", 0.9);

    if(mf<=0 || nf<=0) return 0.0;

    int sz = mf * nf;
    double deltaL = Xf / (nf * y);

    QVector<QVector<double>> I(sz+mf+1, QVector<double>(sz+mf+1));
    QVector<double> F(sz+mf+1, 0);

    for(int i=1; i<=mf; ++i) {
        for(int j=1; j<=nf; ++j) {
            int r = (i-1)*nf + j - 1;
            for(int k=1; k<=mf; ++k) {
                for(int v=1; v<=nf; ++v) {
                    double val = integralBesselK0(
                        ((2*v-2*nf-1)/(2.0*nf)*Xf)/y, (yy+(y-2*yy)/(mf-1)*(k-1))/y, (yy+(y-2*yy)/(mf-1)*(i-1))/y,
                        z*(omega*z*(1-omega)+lambda)/(lambda+(1-omega)*z), ((j-nf-1)/(double)nf*Xf)/y, ((j-nf)/(double)nf*Xf)/y
                        );
                    if(i==k) val -= (j==v ? M_PI/(CFD*8)*deltaL*deltaL : M_PI/CFD*std::abs(j-v)*deltaL*deltaL);
                    I[r][(k-1)*nf + v - 1] = val;
                }
            }
            I[r][sz] = (2*M_PI/CFD) * (((j-nf)/(double)nf*Xf/y) - ((j-nf-1)/(double)nf*Xf/y));
            I[r][sz+mf] = -1.0;
        }
    }
    for(int i=1; i<=mf; ++i) {
        for(int j=(i-1)*nf; j<i*nf; ++j) I[sz+i-1][j] = deltaL;
        I[sz+i-1][sz+i-1] = -1.0;
    }
    for(int i=sz; i<sz+mf; ++i) I[sz+mf][i] = 1.0;

    F[sz+mf] = 1.0/z;
    QVector<double> res = solveLinearSystem(I, F);
    if(res.size() <= sz+mf) return 0.0;
    double pd1 = res[sz+mf];
    return (z * pd1 + S) / (z + z*z*cD * (z * pd1 + S));
}

double ModelManager::e_function(double z, int i, int j, int k, int v, int mf, int nf, double omega, double lambda, double Xf, double yy, double y)
{
    double fz = (z * (omega * z * (1 - omega) + lambda)) / (lambda + (1 - omega) * z);
    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;
    double Xkv = (2*v - 2*nf - 1) / (double)(2*nf) * Xf;
    double XDkv = Xkv / y;
    double yij = yy + (y - 2*yy) / (mf - 1) * (i - 1);
    double yDij = yij / y;
    double Ykv = yy + (y - 2*yy) / (mf - 1) * (k - 1);
    double YDkv = Ykv / y;
    return integralBesselK0(XDkv, YDkv, yDij, fz, xDij, xDij1);
}

double ModelManager::f_function(int j, int nf, double Xf, double y)
{
    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;
    return xDij1 - xDij;
}
