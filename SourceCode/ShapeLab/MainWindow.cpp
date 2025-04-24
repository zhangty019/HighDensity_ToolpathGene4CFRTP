#include "stdafx.h"

#include "MainWindow.h"
#include "ui_MainWindow.h"
#include <QFileDialog>
#include <QtDebug>
#include <QDesktopWidget>
#include <QCoreApplication>
#include <QMimeData>
#include <QTreeView>
#include <QThread>
#include <QTimer>
#include <QDateTime>
#include <QMessageBox>
#include <QScreen>
#include <QStyleFactory>
#include <fstream>

#include "../GLKLib/GLKCameraTool.h"
#include "../GLKLib/InteractiveTool.h"
#include "../GLKLib/GLKMatrixLib.h"
#include "../GLKLib/GLKGeometry.h"
#include "../QMeshLib/QMeshPatch.h"
#include "../QMeshLib/QMeshTetra.h"
#include "../QMeshLib/QMeshFace.h"
#include "../QMeshLib/QMeshEdge.h"
#include "../QMeshLib/QMeshNode.h"

#include "fileIO.h"
#include "directionFieldOpt.h"
#include "meshOperation.h"
#include "DDGMeshProcessing.h"
#include "toolpathGeneration_stripe.h"


MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QApplication::setStyle(QStyleFactory::create("Fusion"));

    signalMapper = new QSignalMapper(this);
    addToolBar(ui->toolBar);
    addToolBar(ui->navigationToolBar);
    addToolBar(ui->selectionToolBar);

    createTreeView();
    createActions();

    pGLK = new GLKLib();
    ui->horizontalLayout->addWidget(pGLK);
    ui->horizontalLayout->setMargin(0);
    pGLK->setFocus();

    pGLK->clear_tools();
    pGLK->set_tool(new GLKCameraTool(pGLK, ORBITPAN));

    //connect timer with timer function
    //connect(&Gcode_timer, SIGNAL(timeout()), this, SLOT(doTimerGcodeMoving()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::createActions()
{
    // file IO
    connect(ui->actionOpen, SIGNAL(triggered(bool)), this, SLOT(open()));
    connect(ui->actionSave, SIGNAL(triggered(bool)), this, SLOT(save()));
    connect(ui->actionSaveSelection, SIGNAL(triggered(bool)), this, SLOT(saveSelection()));
    connect(ui->actionReadSelection, SIGNAL(triggered(bool)), this, SLOT(readSelection()));

    // navigation
    connect(ui->actionFront, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionBack, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionTop, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionBottom, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionLeft, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionRight, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionIsometric, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionZoom_In, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionZoom_Out, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionZoom_All, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionZoom_Window, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    signalMapper->setMapping(ui->actionFront, 0);
    signalMapper->setMapping(ui->actionBack, 1);
    signalMapper->setMapping(ui->actionTop, 2);
    signalMapper->setMapping(ui->actionBottom, 3);
    signalMapper->setMapping(ui->actionLeft, 4);
    signalMapper->setMapping(ui->actionRight, 5);
    signalMapper->setMapping(ui->actionIsometric, 6);
    signalMapper->setMapping(ui->actionZoom_In, 7);
    signalMapper->setMapping(ui->actionZoom_Out, 8);
    signalMapper->setMapping(ui->actionZoom_All, 9);
    signalMapper->setMapping(ui->actionZoom_Window, 10);

    // view
    connect(ui->actionShade, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionMesh, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionNode, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionProfile, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionFaceNormal, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionNodeNormal, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    signalMapper->setMapping(ui->actionShade, 20);
    signalMapper->setMapping(ui->actionMesh, 21);
    signalMapper->setMapping(ui->actionNode, 22);
    signalMapper->setMapping(ui->actionProfile, 23);
    signalMapper->setMapping(ui->actionFaceNormal, 24);
    signalMapper->setMapping(ui->actionNodeNormal, 25);
    ui->actionShade->setChecked(true);

    connect(ui->actionShifttoOrigin, SIGNAL(triggered(bool)), this, SLOT(shiftToOrigin()));

    // select
    connect(ui->actionSelectNode, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionSelectEdge, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionSelectFace, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionSelectFix, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));
    connect(ui->actionSelectHandle, SIGNAL(triggered(bool)), signalMapper, SLOT(map()));

    signalMapper->setMapping(ui->actionSelectNode, 30);
    signalMapper->setMapping(ui->actionSelectEdge, 31);
    signalMapper->setMapping(ui->actionSelectFace, 32);
    signalMapper->setMapping(ui->actionSelectFix, 33);
    signalMapper->setMapping(ui->actionSelectHandle, 34);


    connect(signalMapper, SIGNAL(mapped(int)), this, SLOT(signalNavigation(int)));

    //Button for display
    connect(ui->pushButton_ShowAll, SIGNAL(released()), this, SLOT(show_All()));
    connect(ui->spinBox_ShowIndex, SIGNAL(valueChanged(int)), this, SLOT(update_Display()));
    connect(ui->pushButton_showAllEdge, SIGNAL(released()), this, SLOT(show_AllEdge()));
    connect(ui->spinBox_showEdgeIndex, SIGNAL(valueChanged(int)), this, SLOT(update_EdgeDisplay()));

    //Button for stripPath generation
    connect(ui->pushButton_readLayerAndfStress, SIGNAL(released()), this, SLOT(input_Layer_fStress()));
    connect(ui->pushButton_DirectionFieldOpt, SIGNAL(released()), this, SLOT(directionField_Opt()));
    connect(ui->pushButton_stripGeneration, SIGNAL(released()), this, SLOT(stripGeneration()));
    connect(ui->pushButton_toolpathGeneration_strip, SIGNAL(released()), this, SLOT(toolpathGeneration_fromStripe()));
    connect(ui->pushButton_offsetLayers, SIGNAL(released()), this, SLOT(offset_Layer()));
    connect(ui->pushButton_postProcessing4path, SIGNAL(released()), this, SLOT(postProcessing_4_stripPath()));
    
}

void MainWindow::open()
{
    QString filenameStr = QFileDialog::getOpenFileName(this, tr("Open File,"), "..", tr(""));
    QFileInfo fileInfo(filenameStr);
    QString fileSuffix = fileInfo.suffix();
    QByteArray filenameArray = filenameStr.toLatin1();
    char* filename = filenameArray.data();

    // set polygen name
    std::string strFilename(filename);
    std::size_t foundStart = strFilename.find_last_of("/");
    std::size_t foundEnd = strFilename.find_last_of(".");
    std::string modelName;
    modelName = strFilename.substr(0, foundEnd);
    modelName = modelName.substr(foundStart + 1);

    if (QString::compare(fileSuffix, "obj") == 0) {
        PolygenMesh* polygenMesh = new PolygenMesh(UNDEFINED);
        polygenMesh->ImportOBJFile(filename, modelName);
        polygenMesh->BuildGLList(polygenMesh->m_bVertexNormalShading);
        pGLK->AddDisplayObj(polygenMesh, true);
        polygenMeshList.AddTail(polygenMesh);
    }

    else if (QString::compare(fileSuffix, "tet") == 0) {
        PolygenMesh* polygenMesh = new PolygenMesh(TET_MODEL);
        std::cout << filename << std::endl;
        std::cout << modelName << std::endl;
        polygenMesh->ImportTETFile(filename, modelName);
        polygenMesh->BuildGLList(polygenMesh->m_bVertexNormalShading);
        pGLK->AddDisplayObj(polygenMesh, true);
        polygenMeshList.AddTail(polygenMesh);
    }

    updateTree();

    shiftToOrigin();
    pGLK->refresh(true);
}

void MainWindow::save()
{
    PolygenMesh* polygenMesh = getSelectedPolygenMesh();
    if (!polygenMesh)
        polygenMesh = (PolygenMesh*)polygenMeshList.GetHead();
    if (!polygenMesh)
        return;
    QString filenameStr = QFileDialog::getSaveFileName(this, tr("OBJ File Export,"), "..", tr("OBJ(*.obj)"));
    QFileInfo fileInfo(filenameStr);
    QString fileSuffix = fileInfo.suffix();

    if (QString::compare(fileSuffix, "obj") == 0) {
        QFile exportFile(filenameStr);
        if (exportFile.open(QFile::WriteOnly | QFile::Truncate)) {
            QTextStream out(&exportFile);
            for (GLKPOSITION posMesh = polygenMesh->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
                QMeshPatch* patch = (QMeshPatch*)polygenMesh->GetMeshList().GetNext(posMesh);
                for (GLKPOSITION posNode = patch->GetNodeList().GetHeadPosition(); posNode != nullptr;) {
                    QMeshNode* node = (QMeshNode*)patch->GetNodeList().GetNext(posNode);
                    double xx, yy, zz;
                    node->GetCoord3D(xx, yy, zz);
                    float r, g, b;
                    node->GetColor(r, g, b);
                    out << "v " << xx << " " << yy << " " << zz << " " << node->value1 << endl;
                }
                for (GLKPOSITION posFace = patch->GetFaceList().GetHeadPosition(); posFace != nullptr;) {
                    QMeshFace* face = (QMeshFace*)patch->GetFaceList().GetNext(posFace);
                    out << "f " << face->GetNodeRecordPtr(0)->GetIndexNo() << " " << face->GetNodeRecordPtr(1)->GetIndexNo() << " " << face->GetNodeRecordPtr(2)->GetIndexNo() << endl;
                }
            }
        }
        exportFile.close();
    }
}

void MainWindow::saveSelection()
{
    //printf("%s exported\n", Model->ModelName);

    PolygenMesh* polygenMesh = getSelectedPolygenMesh();
    if (!polygenMesh)
        polygenMesh = (PolygenMesh*)polygenMeshList.GetHead();
    QMeshPatch* patch = (QMeshPatch*)polygenMesh->GetMeshList().GetHead();

    std::string filename = polygenMesh->getModelName();
    const char* c = filename.c_str();
    char* cstr = new char[filename.length() + 1];
    strcpy(cstr, filename.c_str());

    const char* split = ".";
    char* p = strtok(cstr, split);

    char output_filename[256];
    strcpy(output_filename, "..\\selection_file\\");
    strcat(output_filename, cstr);
    char filetype[64];
    strcpy(filetype, ".txt");
    strcat(output_filename, filetype);

    ofstream nodeSelection(output_filename);
    if (!nodeSelection)
        cerr << "Sorry!We were unable to build the file NodeSelect!\n";
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* CheckNode = (QMeshNode*)patch->GetNodeList().GetNext(Pos);
        nodeSelection << CheckNode->GetIndexNo() << ":";
        //for the selection of fixing part
        if (CheckNode->isFixed == true) nodeSelection << "1:";
        else nodeSelection << "0:";
        //for the selection of hard part
        if (CheckNode->isHandle == true) nodeSelection << "1:" << endl;
        else nodeSelection << "0:" << endl;
    }

    nodeSelection.close();
    printf("Finish output selection \n");
}

void MainWindow::readSelection()
{
    PolygenMesh* polygenMesh = getSelectedPolygenMesh();
    if (!polygenMesh)
        polygenMesh = (PolygenMesh*)polygenMeshList.GetHead();
    QMeshPatch* patch = (QMeshPatch*)polygenMesh->GetMeshList().GetHead();

    std::string filename = polygenMesh->getModelName();
    const char* c = filename.c_str();

    char* cstr = new char[filename.length() + 1];
    strcpy(cstr, filename.c_str());

    const char* split = ".";
    char* p = strtok(cstr, split);

    char input_filename[256];
    strcpy(input_filename, "..\\selection_file\\");
    strcat(input_filename, cstr);
    char filetype[64];
    strcpy(filetype, ".txt");
    strcat(input_filename, filetype);

    ifstream nodeSelect(input_filename);
    if (!nodeSelect)
        cerr << "Sorry!We were unable to open the file!\n";
    vector<int> NodeIndex(patch->GetNodeNumber()), checkNodeFixed(patch->GetNodeNumber()), checkNodeHandle(patch->GetNodeNumber());
    //string line;
    int LineIndex1 = 0;
    string sss;
    while (getline(nodeSelect, sss)) {
        const char* c = sss.c_str();
        sscanf(c, "%d:%d:%d", &NodeIndex[LineIndex1], &checkNodeFixed[LineIndex1], &checkNodeHandle[LineIndex1]);
        LineIndex1++;
    }

    nodeSelect.close();
    for (GLKPOSITION Pos = patch->GetNodeList().GetHeadPosition(); Pos;) {
        QMeshNode* CheckNode = (QMeshNode*)patch->GetNodeList().GetNext(Pos);
        if (checkNodeFixed[CheckNode->GetIndexNo() - 1] == 1) CheckNode->isFixed = true;
        if (checkNodeHandle[CheckNode->GetIndexNo() - 1] == 1) CheckNode->isHandle = true;
    }

    for (GLKPOSITION Pos = patch->GetFaceList().GetHeadPosition(); Pos != NULL;)
    {
        QMeshFace* Face = (QMeshFace*)patch->GetFaceList().GetNext(Pos);
        if (Face->GetNodeRecordPtr(0)->isHandle == true &&
            Face->GetNodeRecordPtr(1)->isHandle == true &&
            Face->GetNodeRecordPtr(2)->isHandle == true)
            Face->isHandleDraw = true;
        else Face->isHandleDraw = false;

        if (Face->GetNodeRecordPtr(0)->isFixed == true &&
            Face->GetNodeRecordPtr(1)->isFixed == true &&
            Face->GetNodeRecordPtr(2)->isFixed == true)
            Face->isFixedDraw = true;
        else Face->isFixedDraw = false;
    }
    printf("Finish input selection \n");
    pGLK->refresh(true);

}

void MainWindow::mouseMoveEvent(QMouseEvent* event)
{
    //QMouseEvent *e = (QMouseEvent*)event;
    //QPoint pos = e->pos();
    //cout << "Mouse position updated" << endl;
    //double wx, wy, wz;
    //pGLK->screen_to_wcl(100.0, 100.0, wx, wy, wz);
    //ui->CorrdinateMouse->setText(QString("X = %1").arg(wx));

    //QString text;
    //text = QString("%1 X %2").arg(event->pos().x()).arg(event->pos().y());
    ///** Update the info text */
    //ui->statusBar->showMessage(text);
}

void MainWindow::signalNavigation(int flag)
{
    if (flag <= 10)
        pGLK->setNavigation(flag);
    if (flag >= 20 && flag <= 25) {
        pGLK->setViewModel(flag - 20);
        switch (flag) {
        case 20:
            ui->actionShade->setChecked(pGLK->getViewModel(0));
            break;
        case 21:
            ui->actionMesh->setChecked(pGLK->getViewModel(1));
            break;
        case 22:
            ui->actionNode->setChecked(pGLK->getViewModel(2));
            break;
        case 23:
            ui->actionProfile->setChecked(pGLK->getViewModel(3));
            break;
        case 24:
            ui->actionFaceNormal->setChecked(pGLK->getViewModel(4));
            break;
        case 25:
            ui->actionNodeNormal->setChecked(pGLK->getViewModel(5));
            break;
        }
    }
    //  if (flag==30 || flag==31 || flag==32 || flag == 33 || flag == 34){
    //      InteractiveTool *tool;
    //      switch (flag) {
    //      case 30:
    //          tool = new InteractiveTool(pGLK, &polygenMeshList, (GLKMouseTool*)pGLK->GetCurrentTool(), NODE, ui->boxDeselect->isChecked());
    //          break;
    //      case 31:
    //          tool = new InteractiveTool(pGLK, &polygenMeshList, (GLKMouseTool*)pGLK->GetCurrentTool(), EDGE, ui->boxDeselect->isChecked());
    //          break;
    //      case 32:
    //          tool = new InteractiveTool(pGLK, &polygenMeshList, (GLKMouseTool*)pGLK->GetCurrentTool(), FACE, ui->boxDeselect->isChecked());
    //          break;
          //case 33:
          //	tool = new InteractiveTool(pGLK, &polygenMeshList, (GLKMouseTool*)pGLK->GetCurrentTool(), FIX, ui->boxDeselect->isChecked());
          //	break;
          //case 34:
          //	tool = new InteractiveTool(pGLK, &polygenMeshList, (GLKMouseTool*)pGLK->GetCurrentTool(), NHANDLE, ui->boxDeselect->isChecked());
          //	break;
    //      }
    //      pGLK->set_tool(tool);
    //  }
}

void MainWindow::shiftToOrigin()
{

}

void MainWindow::dragEnterEvent(QDragEnterEvent* event)
{
    if (event->mimeData()->hasUrls())
        event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* event)
{
    QString filenameStr;
    foreach(const QUrl & url, event->mimeData()->urls())
        filenameStr = url.toLocalFile();
    QByteArray filenameArray = filenameStr.toLatin1();
    char* filename = filenameArray.data();

    PolygenMesh* polygenMesh = new PolygenMesh(UNDEFINED);

    // set polygen name
    std::string strFilename(filename);
    std::size_t foundStart = strFilename.find_last_of("/");
    std::size_t foundEnd = strFilename.find_last_of(".");
    std::string modelName;
    modelName = strFilename.substr(0, foundEnd);
    modelName = modelName.substr(foundStart + 1);
    int i = 0;
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygen = (PolygenMesh*)polygenMeshList.GetNext(pos);
        std::string name = (polygen->getModelName()).substr(0, (polygen->getModelName()).find(' '));
        if (name == modelName)
            i++;
    }
    if (i > 0)
        modelName += " " + std::to_string(i);

    QFileInfo fileInfo(filenameStr);
    QString fileSuffix = fileInfo.suffix();
    if (QString::compare(fileSuffix, "obj") == 0) {
        polygenMesh->ImportOBJFile(filename, modelName);
    }
    else if (QString::compare(fileSuffix, "tet") == 0) {
        polygenMesh->ImportTETFile(filename, modelName);
        polygenMesh->meshType = TET_MODEL;
    }
    polygenMesh->m_bVertexNormalShading = false;
    polygenMesh->BuildGLList(polygenMesh->m_bVertexNormalShading);
    pGLK->AddDisplayObj(polygenMesh, true);
    polygenMeshList.AddTail(polygenMesh);

    updateTree();
}

void MainWindow::createTreeView()
{
    treeModel = new QStandardItemModel();
    ui->treeView->setModel(treeModel);
    ui->treeView->setHeaderHidden(true);
    ui->treeView->setContextMenuPolicy(Qt::CustomContextMenu);
    ui->treeView->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->treeView->expandAll();
}

void MainWindow::updateTree()
{
    treeModel->clear();
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);
        QString modelName = QString::fromStdString(polygenMesh->getModelName());
        QStandardItem* modelListItem = new QStandardItem(modelName);
        modelListItem->setCheckable(true);
        modelListItem->setCheckState(Qt::Checked);
        treeModel->appendRow(modelListItem);
    }
    pGLK->refresh(true);
}

PolygenMesh* MainWindow::getSelectedPolygenMesh()
{
    if (!treeModel->hasChildren())
        return nullptr;
    QModelIndex index = ui->treeView->currentIndex();
    QString selectedModelName = index.data(Qt::DisplayRole).toString();
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);
        QString modelName = QString::fromStdString(polygenMesh->getModelName());
        if (QString::compare(selectedModelName, modelName) == 0)
            return polygenMesh;
    }
    return nullptr;
}

void MainWindow::on_pushButton_clearAll_clicked()
{
    int i = 0;
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr; i++) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);
        QMeshPatch* patch = (QMeshPatch*)polygenMesh->GetMeshList().GetHead();
        if (i < 2)
            continue;
        for (GLKPOSITION pos2 = patch->GetFaceList().GetHeadPosition(); pos2 != nullptr;) {
            QMeshFace* face = (QMeshFace*)patch->GetFaceList().GetNext(pos2);
            face->m_nIdentifiedPatchIndex = 0;
        }
    }
    pGLK->refresh(true);
}

void MainWindow::on_treeView_clicked(const QModelIndex& index)
{
    ui->treeView->currentIndex();
    QStandardItem* modelListItem = treeModel->itemFromIndex(index);
    ui->treeView->setCurrentIndex(index);
    PolygenMesh* polygenMesh = getSelectedPolygenMesh();
    if (modelListItem->checkState() == Qt::Checked)
        polygenMesh->bShow = true;
    else
        polygenMesh->bShow = false;
    pGLK->refresh(true);
}

PolygenMesh* MainWindow::_buildPolygenMesh(mesh_type type, std::string name) {

    PolygenMesh* newMesh = new PolygenMesh(type);
    newMesh->setModelName(name);
    newMesh->BuildGLList(newMesh->m_bVertexNormalShading);
    pGLK->AddDisplayObj(newMesh, true);
    polygenMeshList.AddTail(newMesh);
    updateTree();
    return newMesh;

}

PolygenMesh* MainWindow::_detectPolygenMesh(mesh_type type) {

    PolygenMesh* detectedMesh = NULL;
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* thispolygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);
        if (thispolygenMesh->meshType == type) {
            detectedMesh = thispolygenMesh; break;
        }
    }
    return detectedMesh;
}

void MainWindow::update_Display() {

    bool single = ui->checkBox_EachSwitch->isChecked();
    int currentLayerIndex = ui->spinBox_ShowIndex->value();

    int selectedField = ui->comboBox_selectShow->currentIndex();

    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);
        if (polygenMesh->meshType != CURVED_LAYER
            && polygenMesh->meshType != TOOL_PATH) continue;

        for (GLKPOSITION posMesh = polygenMesh->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
            QMeshPatch* Patch = (QMeshPatch*)polygenMesh->GetMeshList().GetNext(posMesh);

            Patch->drawThisPatch = false;

            if (single) {
                if (Patch->GetIndexNo() == currentLayerIndex)
                    Patch->drawThisPatch = true;
            }
            else {
                if (Patch->GetIndexNo() <= currentLayerIndex)
                    Patch->drawThisPatch = true;
            }

            //update the field to show
            if (selectedField == 0) {
                Patch->draw_stripPattern = false;
                Patch->draw_Field_onFace = false;
                Patch->draw_directionField_onFace = false;
            }
            else if (selectedField == 1) {
                Patch->draw_stripPattern = false;
                Patch->draw_Field_onFace = true;
                Patch->draw_directionField_onFace = true;
            }
            else if (selectedField == 2) {
                Patch->draw_stripPattern = false;
                Patch->draw_Field_onFace = true;
                Patch->draw_directionField_onFace = false;
            }
            else if (selectedField == 3) {
                Patch->draw_stripPattern = true;
            }
            else {

            }
        }
    }
    pGLK->refresh(true);
}

void MainWindow::update_EdgeDisplay() {

    int currentEdgeIndex = ui->spinBox_showEdgeIndex->value();

    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);

        if (polygenMesh->meshType != TOOL_PATH) continue;

        for (GLKPOSITION posMesh = polygenMesh->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
            QMeshPatch* Patch = (QMeshPatch*)polygenMesh->GetMeshList().GetNext(posMesh);

            for (GLKPOSITION Pos = Patch->GetEdgeList().GetHeadPosition(); Pos;) {
                QMeshEdge* Edge = (QMeshEdge*)Patch->GetEdgeList().GetNext(Pos);

                Edge->drawThisEdge = false;
                if (Edge->GetIndexNo() < currentEdgeIndex)
                    Edge->drawThisEdge = true;
            }
        }
    }
    pGLK->refresh(true);
}


void MainWindow::show_All() {
    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);

        if (polygenMesh->meshType != CURVED_LAYER
            && polygenMesh->meshType != TOOL_PATH) continue;

        for (GLKPOSITION posMesh = polygenMesh->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
            QMeshPatch* Patch = (QMeshPatch*)polygenMesh->GetMeshList().GetNext(posMesh);
            Patch->drawThisPatch = true;
        }
    }
    pGLK->refresh(true);
}

void MainWindow::show_AllEdge() {

    for (GLKPOSITION pos = polygenMeshList.GetHeadPosition(); pos != nullptr;) {
        PolygenMesh* polygenMesh = (PolygenMesh*)polygenMeshList.GetNext(pos);

        if (polygenMesh->meshType != TOOL_PATH) continue;

        for (GLKPOSITION posMesh = polygenMesh->GetMeshList().GetHeadPosition(); posMesh != nullptr;) {
            QMeshPatch* Patch = (QMeshPatch*)polygenMesh->GetMeshList().GetNext(posMesh);

            for (GLKPOSITION Pos = Patch->GetEdgeList().GetHeadPosition(); Pos;) {
                QMeshEdge* Edge = (QMeshEdge*)Patch->GetEdgeList().GetNext(Pos);

                Edge->drawThisEdge = true;
            }
        }
    }
    pGLK->refresh(true);
}

void MainWindow::offset_Layer() {

    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();
    PolygenMesh* isoLayerSet = this->_buildPolygenMesh(CURVED_LAYER, model_name + "_Layers");

    fileIO* IO_operator = new fileIO();
    std::string path = "../DataSet/" + model_name + "/CurvedLayer";
    IO_operator->read_layer_files(isoLayerSet, path);

    PolygenMesh* offsetLayerSet = this->_buildPolygenMesh(CURVED_LAYER, "offset_Layers");

    meshOperation* meshOperator = new meshOperation();
    meshOperator->initial(isoLayerSet);
    meshOperator->offsetMesh(isoLayerSet, offsetLayerSet);
    meshOperator->initial(offsetLayerSet);
    delete meshOperator;

    path = "../DataSet/" + model_name + "/CurvedLayer_Offset";
    IO_operator->outputIsoSurface(offsetLayerSet, path);
    delete IO_operator;

    pGLK->refresh(true);
    pGLK->Zoom_All_in_View();
    std::cout << "--> Finish Layer Offset and output into: " << path << std::endl;
}

void MainWindow::input_Layer_fStress() {

    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();
    PolygenMesh* isoLayerSet = this->_buildPolygenMesh(CURVED_LAYER, model_name + "_Layers");

    fileIO* IO_operator = new fileIO();
    std::string path = "../DataSet/" + model_name + "/CurvedLayer_Offset";
    IO_operator->read_layer_files(isoLayerSet, path);
    std::cout << "The index of face start from 1!" << std::endl;

    meshOperation* meshOperator = new meshOperation();
    meshOperator->setBoundaryNodes(isoLayerSet);
    delete meshOperator;

    path = "../DataSet/" + model_name + "/StressField_onFace";
    IO_operator->input_Projeted_stressField(isoLayerSet, path);
    delete IO_operator;

    updateTree();
    pGLK->Zoom_All_in_View();

    std::cout << "--> Finish input the layers and attached stress field on the face" << endl;
}

void MainWindow::directionField_Opt() {

    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();

    PolygenMesh* isoLayerSet = this->_detectPolygenMesh(CURVED_LAYER);

    directionFieldOpt* dField_opt = new directionFieldOpt();
    dField_opt->initialize(isoLayerSet, model_name);
    dField_opt->run();
    delete dField_opt;

    
    fileIO* IO_operator = new fileIO();
    std::string path = "../DataSet/" + model_name + "/DirectionField_onNode";
    IO_operator->outputDirectionField_onNodes(isoLayerSet, path);
    delete IO_operator;

    pGLK->refresh(true);

    std::cout << "--> Finish direction field optimization." << endl;
}

void MainWindow::stripGeneration() {

    //preparation
    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();
    fileIO* IO_operator = new fileIO();

    PolygenMesh* isoLayerSet = this->_detectPolygenMesh(CURVED_LAYER);
    if (isoLayerSet) {
        isoLayerSet->ClearAll();
        isoLayerSet->setModelName(isoLayerSet->getModelName() + std::string("New"));
    }
    else
        isoLayerSet = this->_buildPolygenMesh(CURVED_LAYER, model_name + "_LayersNew");

    QImage texture = QImage(QString("../DataSet/textures/texture2.png"));
    isoLayerSet->SetTexture(texture.width(), texture.height(), texture.bits(), texture.byteCount());

    //get the nameSet of the curved layers
    std::string layerPath = "../DataSet/" + model_name + "/CurvedLayer_Offset";
    std::vector<std::string> layer_name_set;
    IO_operator->get_layer_files_name_set(layerPath, layer_name_set);

    std::vector<int> layer_idxSet = this->_setLayer_idxSet(model_name);

    //get the nameSet of the directionField_onNode files
    std::string dFieldPath = "../DataSet/" + model_name + "/DirectionField_onNode";
    std::vector<std::string> dField_name_set;
    IO_operator->get_dField_files_name_set(dFieldPath, dField_name_set);

    //process the each layer in the loop
    for (size_t _id = 0; _id < layer_name_set.size(); _id++) {

        if (std::find(layer_idxSet.begin(), layer_idxSet.end(), _id) == layer_idxSet.end()) continue;

        if (!IO_operator->_find_matched_dField(layer_name_set[_id], dField_name_set)) continue;

        std::string matched_dField_name =
            layer_name_set[_id].substr(0, layer_name_set[_id].find_last_of('.')) + ".txt";
        std::string one_dFieldPath = dFieldPath + "/" + matched_dField_name;
        std::vector<std::vector<double>> dField_matrix;
        IO_operator->read_dField_onNode(one_dFieldPath, dField_matrix);

        DDGMeshProcessing* ddgmp = new DDGMeshProcessing();
        std::string one_layerPath = layerPath + "/" + layer_name_set[_id];

        ddgmp->setMesh(one_layerPath);
        ddgmp->setInitField(dField_matrix);
        ddgmp->meshProcessing();
        //from DDG::Mesh to QMeshPatch
        QMeshPatch* one_newPatch = ddgmp->toQMeshPatch();
        delete ddgmp;
        one_newPatch->SetIndexNo(isoLayerSet->GetMeshList().GetCount()); //id-> 0
        isoLayerSet->GetMeshList().AddTail(one_newPatch);
        one_newPatch->patchName = layer_name_set[_id].data();
        one_newPatch->draw_stripPattern = true;
    }

    updateTree();
    pGLK->Zoom_All_in_View();
    ui->comboBox_selectShow->setCurrentIndex(3);
    std::cout << "--> Finish strip pattern generation." << endl;
}

std::vector<int> MainWindow::_setLayer_idxSet(std::string model_name) {

    std::vector<int> layerRange = { 0 };

    if (model_name == "TshapeBracketNew") layerRange = { 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 
        ,29 ,30 ,31 ,32 ,33 ,34 ,35 ,36 ,37, 38, 39 }; //19,39
    if (model_name == "bladeSmall") layerRange = { 0, 1, 2, 3 };
    if (model_name == "ncc2") layerRange = { 9 ,12 ,15 ,18 ,21 ,24 ,27 ,30 ,33 ,36 ,39 ,44 ,
        47 ,50 ,54 ,56 ,59 ,62 ,69 ,71 ,75 };
    if (model_name == "AB2d") layerRange = { 0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,
     16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,26 ,27 ,28 ,29};
    if (model_name == "BA3d") layerRange = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
    if (model_name == "GEBracketModified") layerRange = { 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
        75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98,
        99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113 };
    if (model_name == "bridge") layerRange = { 49 ,50 ,51 ,52 ,53 ,54 ,55 ,56 ,57 ,58 ,59 ,60 ,61 ,62 ,63 ,64 ,65 ,
        66 ,67 ,68 ,69 ,70 ,71 ,72 ,73 };
    if (model_name == "ncc3_half") layerRange = { 30 ,31 ,32 ,33 ,34 ,35 };
    if (model_name == "shelf_ty") layerRange = { 10, 11, 12, 13, 14, 15, 16, 141, 142, 143, 144, 145, 146, 147 };

    std::cout << "Please check _setLayer_idxSet to decide which layer is picked for strip gene." << std::endl;

    return layerRange;
}

void MainWindow::toolpathGeneration_fromStripe() {

    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();
    PolygenMesh* isoLayerSet = this->_detectPolygenMesh(CURVED_LAYER);
    if (!isoLayerSet) {
        std::cout << "There is no isoLayer, please check!" << endl;
        return;
    }

    PolygenMesh* toolpathSet = this->_buildPolygenMesh(TOOL_PATH, "strip_ToolPath");
    toolpathGeneration_stripe* ToolPathComp_ = new toolpathGeneration_stripe(isoLayerSet, toolpathSet);

    ToolPathComp_->generate_toolPath_marchingSquare();
    delete ToolPathComp_;

    fileIO* IO_operator = new fileIO();
    std::string path = "../DataSet/" + model_name + "/StripPath_Raw";
    IO_operator->output_rawStripToolpath(toolpathSet, path);
    delete IO_operator;

    pGLK->refresh(true);
    pGLK->Zoom_All_in_View();
    std::cout << "--> Finish generating toolpath from stripe and output into: " << path  << std::endl;
}

void MainWindow::postProcessing_4_stripPath(){

    std::string model_name = (ui->lineEdit_SorceDataDir->text()).toStdString();

    PolygenMesh* isoLayerSet = this->_detectPolygenMesh(CURVED_LAYER);
    if (isoLayerSet == NULL) {
        isoLayerSet = _buildPolygenMesh(CURVED_LAYER, model_name + "_LayersNew");
    }
    else {
        isoLayerSet->ClearAll();
        std::cout << "\nThere is already existing a isoLayers PolygenMesh, it has been reconstructed!" << std::endl;
    }

    fileIO* IO_operator = new fileIO();
    std::string path = "../DataSet/" + model_name + "/CurvedLayer_Offset";
    IO_operator->read_layer_files(isoLayerSet, path);
    std::cout << "--> Finish offset layers from: " << path << std::endl;

    PolygenMesh* toolpathSet = this->_detectPolygenMesh(TOOL_PATH);
    if (toolpathSet == NULL) {
        toolpathSet = this->_buildPolygenMesh(TOOL_PATH, "strip_ToolPath");
    }
    else {
        toolpathSet->ClearAll();
        toolpathSet->setModelName("strip_ToolPath");
        std::cout << "\nThere is already existing a toolpathSet PolygenMesh, it has been reconstructed!" << std::endl;
    }

    path = "../DataSet/" + model_name + "/StripPath_Raw";
    std::cout << "Please ensure stripPathes name are the same as the layers in " << path << std::endl;
    IO_operator->input_raw_stripPath(isoLayerSet, toolpathSet, path);

    toolpathGeneration_stripe* ToolPathComp_ = new toolpathGeneration_stripe(isoLayerSet, toolpathSet);
    ToolPathComp_->postProcessing_stripPath();
    delete ToolPathComp_;

    path = "../DataSet/" + model_name + "/StripPath";
    std::cout << "Please ensure stripPathes name are the same as the layers in " << path << std::endl;
    IO_operator->output_StripToolpath(toolpathSet, path);
    delete IO_operator;

    ui->checkBox_EachSwitch->setChecked(true);
    pGLK->refresh(true);
    pGLK->Zoom_All_in_View();
    std::cout << "--> Finish CCF Strip Covering Toolpath Generation.\n" << std::endl;
}