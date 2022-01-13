
#include <matrices.h>
#include <primitives.h>

#include <QGLViewer/simple_viewer.h>
#include "raytracing.h"

//Légende :
//AC : à compléter
//DC : déjà complet (ne pas toucher sauf en cas d'extrême urgence)
//TC : théoriquement complet (mais modifications possibles en fonction de votre projet)

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	SimpleViewer::init_gl();
	SimpleViewer viewer(NOIR,10);
    viewer.setGeometry(10, 10, 800, 800);

	RTracer rt(&viewer);

	// GL init
	viewer.f_init = [&] ()
	{
		Node::prim.gl_init();
		//AC

        // Initialiser avec la matrice de transformation 
        auto mat = glm::mat4(1);

        // Initialiser l'objet BVH du ray tracer
        BVH* bvh = new BVH(mat);
        rt.bvh = bvh;

        // Incliner le volet de la plaque pour'il correspond à la scène
        mat = glm::rotate(mat, -60.0f*3.1428f/180.0f, glm::vec3(1.0f, 0.0f, 0.0f));
        mat = glm::scale(mat, Vec3(0.8f, 0.8f, 0.8f));

        // Faire une copie de la scène de la plaque, tous les objets suivants seront transformés par rapport à la plate
        auto transf = mat;

        // Ajouter la plaque à la scène
        mat = glm::scale(mat, Vec3(10.0f, 10.0f, 0.25f));
        rt.add_cylinder_bvh(rt.bvh, mat, Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7, 0.5);

        // Ajouter l'éponge de milieu à la plaque
        mat = transf;
        mat = glm::scale(mat, Vec3(7.0f, 7.0f, 7.0f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        mat = glm::rotate(mat, 60.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(1.0f,0.0f,0.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 2);

        // Ajouter les 6 éponges distribuées uniformément sur les marques de la plaque
        mat = transf;
        mat = glm::rotate(mat, 25.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        transf = mat;
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(1.0f,0.5f,0.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        mat = transf;
        mat = glm::rotate(mat, 60.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(0.0f,1.0f,1.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        mat = transf;
        mat = glm::rotate(mat, 120.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(1.0f,0.0f,1.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        mat = transf;
        mat = glm::rotate(mat, 180.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(220/255.0f,20/255.0f,60/255.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        mat = transf;
        mat = glm::rotate(mat, 240.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(1.0f,0.0f,0.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        mat = transf;
        mat = glm::rotate(mat, 300.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, -8.0f, 0.0f));
        mat = glm::scale(mat, Vec3(3.2f, 3.2f, 3.2f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(0.0f,1.0f,0.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        // Ajouter l'éponge de haut à celle de milieu
        mat = transf;
        mat = glm::rotate(mat, 60.0f*3.1428f/180.0f, glm::vec3(0.0f, 0.0f, 1.0f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 7.0f));
        mat = glm::scale(mat, Vec3(2.0f, 2.0f, 2.0f));
        mat = glm::translate(mat, Vec3(0.0f, 0.0f, 0.5f));
        rt.add_menger_sponge_bvh(rt.bvh, mat,
                                 Vec3(0.0f,1.0f,0.0f), 0.4f, 0.7f,  // Cube Prameters
                                 Vec3(170/255.0f,169/255.0f,173/255.0f), 0.7f, 0.5f, // Sphere PArameters
                                 1);

        // Initialiser les bounding boxes
        Vec3 minV, maxV;

        // Initialiser le sommet max avec la valeur de moins inifi
        maxV.x = -std::numeric_limits<float>::max();
        maxV.y = -std::numeric_limits<float>::max();
        maxV.z = -std::numeric_limits<float>::max();

        // Initialiser le sommet min avec la valeur de moins inifi
        minV.x = std::numeric_limits<float>::max();
        minV.y = std::numeric_limits<float>::max();
        minV.z = std::numeric_limits<float>::max();

        rt.bvh->getBoundingBox(&minV, &maxV);
	};

	viewer.f_draw = [&] ()
	{
		Node::prim.set_matrices(viewer.getCurrentModelViewMatrix(), viewer.getCurrentProjectionMatrix());
        draw_prim_bvh (rt.bvh); //segfault s'il y a pas d'éléments ajoutés au bvh
        //TC: le bvh est censé englober toute la scène.
	};


	viewer.f_keyPress = [&] (int key, Qt::KeyboardModifiers /*mod*/)
	{
		//TC : entre 0 et 4, le nombre de rebonds pour votre rendu
		switch(key)
		{
			case Qt::Key_0:
				rt.CalcImage(0)->show();
				break;
			case Qt::Key_1:
				rt.CalcImage(1)->show();
				break;
			case Qt::Key_2:
				rt.CalcImage(2)->show();
				break;
			case Qt::Key_3:
				rt.CalcImage(3)->show();
				break;
			case Qt::Key_4:
				rt.CalcImage(4)->show();
				break;
			default:
				break;
		}
	};


	viewer.show();
	return a.exec();
	//TC
}
