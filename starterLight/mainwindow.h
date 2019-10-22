#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDebug>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);

    void boiteEnglobante(MyMesh *_mesh);

    void barycentre(MyMesh * _mesh);

    void resetAllColorsAndThickness(MyMesh* _mesh);

    void calcul_nb_truc(MyMesh *_mesh);

    MyMesh::Normal calcul_normale_une_face(MyMesh *_mesh, int faceId);
    void calcul_normale_sommet(MyMesh *_mesh);

    float aire_face(MyMesh* _mesh, unsigned int faceID);
    void faceObjet(MyMesh* _mesh);


private slots:
    void on_pushButton_chargement_clicked();

    void on_pushButton_generer_clicked();

private:

    MyMesh mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
