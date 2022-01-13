#ifndef __RAYTRACING_H__
#define __RAYTRACING_H__

#include "QGLViewer/simple_viewer.h"
#include "matrices.h"
#include "primitives.h"
#include <thread>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/quaternion.hpp>

#include <QApplication>
#include <QLabel>

const Vec3 ROUGE = {1, 0, 0};      // 0
const Vec3 VERT = {0, 1, 0};       // 1
const Vec3 BLEU = {0, 0, 1};       // 2
const Vec3 JAUNE = {1, 1, 0};      // 3
const Vec3 CYAN = {0, 1, 1};       // 4
const Vec3 MAGENTA = {1, 0, 1};    // 5
const Vec3 BLANC = {1, 1, 1};      // 6
const Vec3 GRIS = {0.5, 0.5, 0.5}; // 7
const Vec3 NOIR = {0, 0, 0};       // 8
const Vec3 ROUGE2 = {0.5, 0, 0};   // 9
const Vec3 VERT2 = {0, 0.5, 0};    // 10
const Vec3 ORANGE = {1, 0.5, 0};   // 11

class Node;
struct Inter {
  float a_min; // La distance de l'intersection depuis l'origine du rayon
  Node *node;  // Les nodes, s'il y en a
  Vec3 Pos; // Le point d'impact du rayon = o + dt (o: origine du rayon, d: la direction, t: la distance de l'intersection depuis l'origine du rayon)
  Vec3 Dir; // La direction du rayon
};

class Node {
public:
  static const float Epsilon;
  Mat4 transfo;     // matrice de transformation.
  Mat4 inv_transfo;
  Vec3 Col;         // couleur.
  float spec;       // indice de spécularité.
  float transp;     // indice de transparence.
  Vec3 maxVertex;   // Le bord maximum
  Vec3 minVertex;   // Le bord minimum


  Node(const Mat4 &m, const Vec3 &c, float sp, float tr);

  inline virtual ~Node() {}
  virtual void draw_gl() const = 0; // la fonction permettant de dessiner la
                                    // primitive dans la vue temps réel.
  virtual bool intersecte(ray r, Inter *I) = 0;  // test d'intersection de la primitive
                                          // (modifier contenu de I).
  virtual Vec3 normal(const Vec3 &P) = 0; // la normale à la primitive
  virtual void getBoundingBox(Vec3* min, Vec3* max) = 0; // Initialiser le bounding box de la primitive
  virtual bool intersectionBB(ray r) = 0; // Vérifier si le rayon et le bounding box de la primitive s'intersecte
  static Primitives prim;
};

class Cube : public Node {
public:
  inline Cube(const Mat4 &_transfo, const Vec3 &_color, float _spec,
              float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(ray r, Inter *I);
  Vec3 normal(const Vec3 &P);
  void getBoundingBox(Vec3* minV, Vec3* maxV);
  bool intersectionBB(ray r);
};

class Sphere : public Node {
public:
  inline Sphere(const Mat4 &_transfo, const Vec3 &_color, float _spec,
                float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(ray r, Inter *I);
  Vec3 normal(const Vec3 &P);
  void getBoundingBox(Vec3* minV, Vec3* maxV);
  bool intersectionBB(ray r);
};

class Cylinder : public Node {
public:
  inline Cylinder(const Mat4 &_transfo, const Vec3 &_color, float _spec,
                  float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(ray r, Inter *I);
  Vec3 normal(const Vec3 &P);
  void getBoundingBox(Vec3* minV, Vec3* maxV);
  bool intersectionBB(ray r);
};

class BVH {
public:
  Mat4 transfo;
  Mat4 inv_transfo;
  std::vector<BVH *> children;
  std::vector<Node *> nodes;

  Vec3 maxVertex; // Le bord maximum
  Vec3 minVertex; // Le bord minimum

  inline BVH(const Mat4 &m) : transfo(m), inv_transfo(glm::inverse(m)) {}

  inline BVH *add_child(const Mat4 &m) {
    BVH *b = new BVH(m);
    children.push_back(b);
    return b;
  }

  float max(float a, float b);
  float min(float a, float b);

  // Vérifie si le rayon et le bounding box du BVH s'intersecte
  bool intersectionBB(ray r);

  inline void add_node(Node *n) { nodes.push_back(n); }

  // Le chemin du BVH pour calculer l'intersection entre les rayons primaires et les primitive de la scène
  void closestIntersection(ray r, Inter *I);

  // La route du BVH pour calculer les shadows, radius qui a comme point d'origine du point de l'intersection d'un rayon primaire de la primitive et la direction vers la source lumière
  void intersecteShadow(ray r, float &sha);

  // Initialiser le bounding box du BVH, où le sommet minimum est le bord minimum et le sommet maximum est la bord maximum
  void getBoundingBox(Vec3* minV, Vec3* maxV);
};

class RTracer {
public:
  BVH *bvh;
  SimpleViewer *sv;
  std::vector<Node *> primitives;
  int depth;
  Vec3 posLum;  // Position de la source de lumière
  static const int nbthr = 3;

  RTracer(SimpleViewer *s);

  void add_cube_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                    float tr);
  void add_sphere_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                      float tr);
  void add_cylinder_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                        float tr);

  void add_menger_sponge_bvh(BVH *b, const Mat4 &m, const Vec3 &colorCube, float specCube,
                        float trCube, const Vec3 &colorSphere, float specSphere, float trSphere, int Depth);

  QLabel *CalcImage(int depth); // calcule l'image par lancer de rayons. Renvoie
                                // de quoi l'afficher.

  // v   à utiliser si vous voulez tester un rendu de la scène avant
  // d'implémenter un bvh.
  Vec3 ColorRay(const Vec3 &P, const Vec3 &V, int depth);
  Vec3 ColorRayBVH(const Vec3 &Origin, const Vec3 &Dir, int rec);

  float max(float a, float b);
};

void draw_prim_bvh(BVH *b);

#endif
