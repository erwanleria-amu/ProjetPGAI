#include "mainwindow.h"
#include "ui_mainwindow.h"


/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh, true);

    calcul_nb_truc(&mesh);
    boiteEnglobante(&mesh);
    barycentre(&mesh);

    histoAngleDiedre(&mesh);

    calcul_normale_sommet(&mesh);
    aire_Objet(&mesh);

    K_Curv(&mesh);

    displayMesh(&mesh, true);


}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }

    //====================== NOS MODIFS=====================
    /*EdgeHandle vh = _mesh->edge_handle(2);
     _mesh->set_color(vh,MyMesh::Color(0, 0, 100));
     _mesh->data(vh).thickness = 120;
     qDebug("nb faces : %d", _mesh->n_faces()); // nombre de faces
     qDebug("nb sommets : %d", _mesh->n_vertices()); // nombre de sommets
     //faces sans voisins

#if 0
     for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
     {
        if(curFace)
     }
#endif*/
}

void MainWindow:: calcul_nb_truc(MyMesh *_mesh) {
   qDebug("YOOOOO");
    qDebug("nb sommets : %d",(int) _mesh->n_vertices()); // nombre de sommets
   qDebug("nb faces : %d",(int) _mesh->n_faces());

   qDebug("nb aretes : %d",(int) _mesh->n_edges());
}

// ==== a continuer
MyMesh::Normal MainWindow:: calcul_normale_une_face(MyMesh *_mesh,
                                          int faceId){

    MyMesh::Normal  normalFh;
    FaceHandle fh = _mesh->face_handle(faceId);
    normalFh=_mesh->normal(fh);
    return normalFh;
}

void MainWindow:: calcul_normale_sommet(MyMesh *_mesh){
    // pour toutes les faces de tous les sommets
    //=> calculer leurs normales => faire moyenne des normales
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)//tous les sommets
    {   int nb_faces=0;
        MyMesh::Normal  normalFh;

        for (MyMesh::VertexFaceIter curFace = _mesh->vf_iter(*curVert); curFace.is_valid(); curFace ++)// ttes les faces
        {
            normalFh += calcul_normale_une_face(_mesh,curFace->idx());
            nb_faces++;
        }

        /*qDebug("pour le sommet %d : ",curVert->idx());
        for(int i=0; i<normalFh.size();i++){
            qDebug("x: %d\n y: %d\n z: %d\n",normalFh[0],normalFh[1],normalFh[2]);
        }
        */
    }
}
// ========
//calcul aire d'une face
float MainWindow::aire_face(MyMesh* _mesh, unsigned int faceID)
{

      VertexHandle vertexA;
      VertexHandle vertexB;
      VertexHandle vertexC;

      FaceHandle curFace = _mesh->face_handle(faceID); // face actuelle selon parametre

      // recuperer demi-arete de notre face actuelle
      MyMesh::HalfedgeHandle curHalfEdge = _mesh->halfedge_handle(curFace);
      //recuperer un sommet
      vertexA = _mesh->to_vertex_handle(curHalfEdge);

      //demi arete suivante
      curHalfEdge = _mesh->next_halfedge_handle(curHalfEdge);
      //sommet suivant
      vertexB = _mesh->to_vertex_handle(curHalfEdge);

      //demi arete suivante
      curHalfEdge = _mesh->next_halfedge_handle(curHalfEdge);
      //sommet suivant
      vertexC = _mesh->to_vertex_handle(curHalfEdge);

      //faire calcul aire
      return (((_mesh->point(vertexB)-_mesh->point(vertexA))%(_mesh->point(vertexC)-_mesh->point(vertexA))).norm())/2;

}

void MainWindow::aire_Objet(MyMesh* _mesh){
    float aireObjet=0;
    float aireFace=0;
    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        aireFace=aire_face(_mesh,curFace->idx());
        qDebug("aireface %f",aireFace);
        aireObjet+=aireFace;
        //qDebug("aire de la face %d: %d\n",curFace->idx(),aireFace);
    }
    qDebug("aire totale de l'objet:%f",aireObjet);
}

void MainWindow::boiteEnglobante(MyMesh *_mesh) {

   //Initialisation avec le premier vertex
    float xMin = _mesh->point(VertexHandle(0))[0], xMax = _mesh->point(VertexHandle(0))[0],
            yMin = _mesh->point(VertexHandle(0))[1], yMax = _mesh->point(VertexHandle(0))[1],
            zMin = _mesh->point(VertexHandle(0))[2], zMax = _mesh->point(VertexHandle(0))[2];

    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end() ; v++) {

        MyMesh::Point p = _mesh->point(*v);

        if(xMin > p[0])
            xMin = p[0];

        if(xMax < p[0])
            xMax = p[0];

        if(yMin > p[1])
            yMin = p[1];

        if(yMax < p[1])
            yMax = p[1];

        if(zMin < p[2])
            zMin = p[2];

        if(zMax > p[2])
            zMax = p[2];

    }

    qDebug("\n\nParametres de la boite englobante \n");
    qDebug() << "xMin : " << xMin << ", xMax : " << xMax << endl;
    qDebug() << "yMin : " << yMin << ", yMax : " << yMax << endl;
    qDebug() << "zMin : " << zMin << ", zMax : " << zMax << endl;

    qDebug()<< "debut boite englobante";
    float x,y,z;
    float x_max,x_min;
    float y_max,y_min;
    float z_max,z_min;

    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        VertexHandle vh = *curVert;
         x = _mesh->point(vh)[0];
         y = _mesh->point(vh)[1];
         z = _mesh->point(vh)[2];

         x_min=x;
         x_max=x;
         y_min=y, y_max=y;
         z_min=z, z_max=z;

        if(x<x_min) x_min=x; if(x>x_max) x_max=x;
        if(y<y_min) y_min=y; if(y>y_max) y_max=y;
        if(z<z_min) z_min=z; if(z>z_max) z_max=z;
    }
    qDebug()<< "Boite englobante: ";
    qDebug()<< "point maximum" << x_max << "," << y_max <<","<<z_max ;
    qDebug()<< "point minimum" << x_min << "," << y_min <<","<<z_min ;
}

void MainWindow::barycentre(MyMesh *_mesh) {

    double bary_accX = 0;
    double bary_accY = 0;
    double bary_accZ = 0;

    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end() ; v++) {

        MyMesh::Point p = _mesh->point(*v);

        bary_accX += p[0];
        bary_accY += p[1];
        bary_accZ += p[2];
    }

    bary_accX = int(bary_accX/_mesh->n_vertices());
    bary_accY = int(bary_accY/_mesh->n_vertices());
    bary_accZ = int(bary_accZ/_mesh->n_vertices());

    qDebug("\n\nBarycentre \n");
    qDebug() << "(X, Y, Z) = " << bary_accX << ", " << bary_accY << ", " << bary_accZ << endl;

}


float MainWindow::angleDiedre(MyMesh* _mesh, int faceID0,  int faceID1)
{

    FaceHandle face0(faceID0);
    FaceHandle face1(faceID1);

    MyMesh::Point n1 = _mesh->calc_face_normal(face0);
    MyMesh::Point n2 = _mesh->calc_face_normal(face1);

    n1.normalize();
    n2.normalize();

    int coef = 1;

    if(asin(n1|n2) < 0)
    {
        coef = -1;
    }

    //qDebug() << "angle " << acos(n1|n2)*180/M_PI << endl;
    return acos(n1|n2)*180/M_PI; //angle en degré
}

void MainWindow::histoAngleDiedre(MyMesh *_mesh)
{
    static unsigned long accAngle[36] = {0};

    for(MyMesh::FaceIter f_it = _mesh->faces_begin() ; f_it != _mesh->faces_end() ; f_it++)
    {
        for(MyMesh::FaceFaceCWIter ffcw_it = _mesh->ff_cwbegin(*f_it) ; ffcw_it.is_valid() ; ffcw_it++)
        {
            /*
            float d = angleDiedre(_mesh, f_it->idx(), ffcw_it->idx());
            float idx = d/10;
            int i = int(idx);
                accAngle[i] += 1;
                */
        }
    }


    qDebug() << "Histogramme d'angles dièdres\n" ;
    for(int angle = 0 ; angle < 36 ; angle++)
    {
        qDebug() << "Angles inférieurs à " << (angle+1)*10 << " degrés :" << accAngle[angle] << endl;
    }
}


float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    FaceHandle face0(faceID);

    VertexHandle vex0(vertexID);

    QVector<MyMesh::Point> points;

    for(MyMesh::FaceVertexIter fv_it = _mesh->fv_begin(face0) ; fv_it.is_valid() ; fv_it++)
    {
        if(_mesh->point(fv_it) != _mesh->point(vex0))
            points.append(_mesh->point(*fv_it));
    }

    MyMesh::Point p0 = _mesh->point(vex0);

    MyMesh::Point u = points[0] - p0;
    MyMesh::Point v = points[1] - p0;
    v.normalize();
    u.normalize();

    return acos(u|v);
}

double MainWindow::aireBarycentre(MyMesh* _mesh, VertexHandle v)
{
    double r = 0;
    for(MyMesh::VertexFaceIter vf_it = _mesh->vf_begin(v) ; vf_it.is_valid() ; vf_it++)
    {
        r += aire_face(_mesh, vf_it->idx()) ;
    }
    return (1/3)*r;
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for(MyMesh::VertexIter v_it = _mesh->vertices_begin() ; v_it != _mesh->vertices_end() ; v_it++)
    {
         double k = (1/(aireBarycentre(_mesh, v_it)));
         double c = 2*M_PI;
         double sum = 0;
         for(MyMesh::VertexFaceIter vf_it = _mesh->vf_begin(v_it) ; vf_it.is_valid() ; vf_it++)
         {
             sum += angleEE(_mesh, v_it->idx(), vf_it->idx());
         }

         double total = k * (c - sum);
         _mesh->data(v_it).value = total;
    }
    qDebug() << "End of function" << endl;
}



//======================FIN  NOS MODIFS=====================


// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


