#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <QObject>
#include <QWidget>
#include <QComboBox>
#include <QStackedWidget>
#include <QMap>
#include <tuple>
#include <QVector>
#include <functional>

class ModelWidget1;
class ModelWidget2;
class ModelWidget3;

typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

class ModelManager : public QObject
{
    Q_OBJECT
public:
    enum ModelType {
        InfiniteConductive = 0,
        FiniteConductive,
        SegmentedMultiCluster
    };
    Q_ENUM(ModelType)

    explicit ModelManager(QWidget* parent = nullptr);
    ~ModelManager();

    void initializeModels(QWidget* parentWidget);
    ModelType getCurrentModelType() const { return m_currentModelType; }

    static QStringList getAvailableModelTypes();
    static QString getModelTypeName(ModelType type);

    QMap<QString, double> getDefaultParameters(ModelType type);

    // [核心修改] 计算接口，供 FittingWidget 和 ModelWidget1 共同使用
    ModelCurveData calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params);

signals:
    void modelSwitched(ModelType newType, ModelType oldType);
    void calculationCompleted(const QString& analysisType, const QMap<QString, double>& results);

private slots:
    void onModelTypeSelectionChanged(int index);
    void onModel1CalculationCompleted(const QString& type, const QMap<QString, double>& results);
    void onModel2CalculationCompleted(const QString& type, const QMap<QString, double>& results);
    void onModel3CalculationCompleted(const QString& type, const QMap<QString, double>& results);

private:
    QWidget* m_mainWidget;
    QComboBox* m_modelTypeCombo;
    QStackedWidget* m_modelStack;
    ModelWidget1* m_modelWidget1;
    ModelWidget2* m_modelWidget2;
    ModelWidget3* m_modelWidget3;
    ModelType m_currentModelType;

    void createMainWidget();
    void setupModelSelection();
    void connectModelSignals();
    void switchToModel(ModelType modelType);

    // --- 静态数学内核 (供内部和外部通用调用) ---
public:
    // 计算 Model 1
    static ModelCurveData calculateModel1(const QMap<QString, double>& params);
    static ModelCurveData calculateModel2(const QMap<QString, double>& params);
    static ModelCurveData calculateModel3(const QMap<QString, double>& params);

private:
    static void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                          double cD, QVector<double>& dpd, QVector<double>& td_dpd);

    static double stefestCoefficient(int i, int N);
    static double factorial(int n);

    // Laplace 空间函数
    static double flaplace1(double z, const QMap<QString, double>& p);
    static double flaplace2(double z, const QMap<QString, double>& p);

    static double e_function(double z, int i, int j, int k, int v, int mf, int nf, double omega, double lambda, double Xf, double yy, double y);
    static double f_function(int j, int nf, double Xf, double y);

    // 积分与 Bessel
    static double integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double a, double b);
    static double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);
    static double gauss15(std::function<double(double)> f, double a, double b);
    static double besselK0(double x);

    static QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
};

#endif // MODELMANAGER_H
