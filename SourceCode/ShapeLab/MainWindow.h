#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSignalMapper>
#include <QStandardItemModel>
#include "../GLKLib/GLKLib.h"
#include "../QMeshLib/PolygenMesh.h"
#include <omp.h>
#include <QTimer>
#include <QLabel>

#define PI		3.141592654
#define DEGREE_TO_ROTATE(x)		0.0174532922222*x
#define ROTATE_TO_DEGREE(x)		57.295780490443*x

using namespace std;

class DeformTet;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = 0);
    ~MainWindow();

public slots:
    //void doTimerGcodeMoving();

private:
    Ui::MainWindow* ui;
    GLKLib* pGLK;
    GLKObList polygenMeshList;

private:
    void createActions();
    void createTreeView();
    PolygenMesh* getSelectedPolygenMesh();

    QSignalMapper* signalMapper;
    QStandardItemModel* treeModel;

private:
    PolygenMesh* _buildPolygenMesh(mesh_type type, std::string name);
    PolygenMesh* _detectPolygenMesh(mesh_type type);
    std::vector<int> _setLayer_idxSet(std::string model_name);

protected:
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);

private slots:
    void open();
    void save();
    void saveSelection();
    void readSelection();

    void signalNavigation(int flag);
    void shiftToOrigin();
    void updateTree();
    void mouseMoveEvent(QMouseEvent* event);
    void on_pushButton_clearAll_clicked();
    void on_treeView_clicked(const QModelIndex& index);

    /*This is Strip Path Printing*/
    void offset_Layer();
    void input_Layer_fStress();
    void directionField_Opt();
    void stripGeneration();
    void toolpathGeneration_fromStripe();
    void postProcessing_4_stripPath();

    /*This is for Display*/
    void update_Display();
    void update_EdgeDisplay();
    void show_All();
    void show_AllEdge();
};

#endif // MAINWINDOW_H
