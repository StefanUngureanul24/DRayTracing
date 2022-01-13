#include "raytracing.h"

Primitives Node::prim;

float clamp01(float v) {
  if (v < 0.0f)
    return 0.0f;
  if (v > 1.0f)
    return 1.0f;
  return v;
}

Vec3 Color(const Vec3 &c) {
  return Vec3(clamp01(c.x), clamp01(c.y), clamp01(c.z));
}

/*
 *  /////////////////////////////////////////
 *  /////////////// PRIMITIVES //////////////
 *  /////////////////////////////////////////
 */
Node::Node(const Mat4 &m, const Vec3 &c, float sp, float tr)
    : transfo(m), inv_transfo(glm::inverse(m)), Col(c), spec(sp), transp(tr) {}

void Cube::draw_gl() const {
  Node::prim.draw_cube(this->transfo, this->Col);
  // DC
}
void Sphere::draw_gl() const {
  Node::prim.draw_sphere(this->transfo, this->Col);
  // DC
}
void Cylinder::draw_gl() const {
  Node::prim.draw_cylinder(this->transfo, this->Col);
  // DC
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * I : L'informations sur le point d'intersection
 * La fonction renvoie true s'il on a bien une intersection
*/
bool Cube::intersecte(ray r, Inter *I) {

    // Si le bounding box et le rayon ne s'intersecte pas
    if(!intersectionBB(r))
    {
        return false;
    }

    // Calculer l'intersection avec le rayon avec la primitive
    Vec4 intersec = this->prim.intersect_cube(transfo, r);

    // Si le rayon intersecte la primitive, on ajout un offset pour la stabilité
    if (intersec.w>-I->node->Epsilon)
    {
        // Mettre à jour la structure I si le node est plus proche que celui contenu dans I
        if (intersec.w < I->a_min)
        {
            I->Dir = r.dir;
            I->node = this;
            I->a_min = intersec.w;
            I->Pos = Vec3(intersec);
        }
        return true;
    }
    return false;
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * renvoie la normale du cube
*/
Vec3 Cube::normal(const Vec3 &P) {
    return this->prim.normal_cube(transfo, P);
}

/*
 * minV : Le bord minimum du bounding box.
 * maxV : Le bord maximum du bounding box.
 */
void Cube::getBoundingBox(Vec3* minV, Vec3* maxV)
{
    // Tableau où on stocke les sommets du cube
    Vec3 vertex[8];
    for (int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++)
        {
            for (int k=0; k<2; k++)
            {
                vertex[i*4+j*2+k].x = i*2-1;
                vertex[i*4+j*2+k].y = j*2-1;
                vertex[i*4+j*2+k].z = k*2-1;
            }
        }
    }

    // Initialiser le sommet max avec moins l'infini
    maxVertex.x = -std::numeric_limits<float>::max();
    maxVertex.y = -std::numeric_limits<float>::max();
    maxVertex.z = -std::numeric_limits<float>::max();

    // Initialiser le sommet min avec l'infini
    minVertex.x = std::numeric_limits<float>::max();
    minVertex.y = std::numeric_limits<float>::max();
    minVertex.z = std::numeric_limits<float>::max();

    // Parcourir les sommets et mettre à jour le bounding box
    for (int i=0; i<8;i++)
    {
        // Transformer vertex[i]
        Vec4 T = transfo*Vec4(vertex[i], 1.0f);

        // Mettre à jour maxVertex
        maxVertex.x = prim.max(maxVertex.x, T.x);
        maxVertex.y = prim.max(maxVertex.y, T.y);
        maxVertex.z = prim.max(maxVertex.z, T.z);

        // Mettre à jour minVertex
        minVertex.x = prim.min(minVertex.x, T.x);
        minVertex.y = prim.min(minVertex.y, T.y);
        minVertex.z = prim.min(minVertex.z, T.z);
    }

    // Renvoie le sommet maximum du bounding box
    maxV->x = prim.max(maxV->x, maxVertex.x);
    maxV->y = prim.max(maxV->y, maxVertex.y);
    maxV->z = prim.max(maxV->z, maxVertex.z);

    // Renvoie le sommet minimum du bounding box
    minV->x = prim.min(minV->x, minVertex.x);
    minV->y = prim.min(minV->y, minVertex.y);
    minV->z = prim.min(minV->z, minVertex.z);
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * La fonction renvoie true s'il on a bien une intersection 
*/
bool Cube::intersectionBB(ray r)
{
    /*
        Ce test utilise l'algorithme suivant:
        Majercik, Alexander, Cyril Crassin, Peter Shirley, and Morgan McGuire. "A ray-box 
        intersection algorithm and efficient dynamic voxel rendering." Journal of Computer Graphics Techniques Vol 7, no. 3 (2018).
    */

    // Renvoie true si le rayon provient de l'interieur du box
    if( (r.x0.x > minVertex.x) && (r.x0.x < maxVertex.x) )
    {
        if( (r.x0.y > minVertex.y) && (r.x0.y < maxVertex.y) )
        {
            if( (r.x0.z > minVertex.z) && (r.x0.z < maxVertex.z) )
            {
                return true;
            }
        }
    }

    // Calculer l'intersection du rayon avec chaque 6 côté du box
    float t1 = (minVertex.x - r.x0.x)*r.n_inv.x;
    float t2 = (maxVertex.x - r.x0.x)*r.n_inv.x;
    float t3 = (minVertex.y - r.x0.y)*r.n_inv.y;
    float t4 = (maxVertex.y - r.x0.y)*r.n_inv.y;
    float t5 = (minVertex.z - r.x0.z)*r.n_inv.z;
    float t6 = (maxVertex.z - r.x0.z)*r.n_inv.z;

    // Calculer la composante max de tmin
    float tmin = prim.max(prim.max(prim.min(t1, t2), prim.min(t3, t4)), prim.min(t5, t6));

    // Calculer la composante min de tmax
    float tmax = prim.min(prim.min(prim.max(t1, t2), prim.max(t3, t4)), prim.max(t5, t6));

    // Si le rayon intersecte le box, mais le box se situe en arrière
    if (tmax < 0)
    {
        return false;
    }

    // S'il n'y a pas d'intersection
    if (tmin > tmax)
    {
        return false;
    }

    return true;
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * I : L'informations sur le point d'intersection
 * La fonction renvoie true s'il on a bien une intersection
*/
bool Sphere::intersecte(ray r, Inter *I)
{
    // La cas où le rayon n'intersecte pas le bounding box
    if(!intersectionBB(r))
    {
        return false;
    }

    // Calculer l'intersection du rayon avec la primitive
    Vec4 intersec = this->prim.intersect_sphere(transfo, r);

    // Si le rayon intersecte la primitive, on ajout un petit offset pour la stabilité
    if (intersec.w>-I->node->Epsilon)
    {
        // Mettre à jour la structure I si le node est plus proche que celui contenu dans I
        if (intersec.w < I->a_min)
        {
            I->Dir = r.dir;
            I->node = this;
            I->a_min = intersec.w;
            I->Pos = Vec3(intersec);
        }
        return true;
    }
    return false;
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Sphere::normal(const Vec3 &P) {
    return this->prim.normal_sphere(transfo, P);
}

/*
 * minV : Le sommet minimum du bounding box.
 * maxV : Le sommet maximum du bounding box.
 */
void Sphere::getBoundingBox(Vec3* minV, Vec3* maxV)
{
    // Tableau où on stocke les sommets du cube
    Vec3 vertex[8];
    for (int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++)
        {
            for (int k=0; k<2; k++)
            {
                vertex[i*4+j*2+k].x = i*2-1;
                vertex[i*4+j*2+k].y = j*2-1;
                vertex[i*4+j*2+k].z = k*2-1;
            }
        }
    }

    // Initialiser le sommet max avec moins l'infini
    maxVertex.x = -std::numeric_limits<float>::max();
    maxVertex.y = -std::numeric_limits<float>::max();
    maxVertex.z = -std::numeric_limits<float>::max();

    // Initialiser le sommet min avec l'infini
    minVertex.x = std::numeric_limits<float>::max();
    minVertex.y = std::numeric_limits<float>::max();
    minVertex.z = std::numeric_limits<float>::max();

    // Parcourir les sommets et mettre à jour le bounding box
    for (int i=0; i<8;i++)
    {
        // Transformer vertex[i]
        Vec4 T = transfo*Vec4(vertex[i], 1.0f);

        // Mettre à jour maxVertex
        maxVertex.x = prim.max(maxVertex.x, T.x);
        maxVertex.y = prim.max(maxVertex.y, T.y);
        maxVertex.z = prim.max(maxVertex.z, T.z);

        // Mettre à jour minVertex
        minVertex.x = prim.min(minVertex.x, T.x);
        minVertex.y = prim.min(minVertex.y, T.y);
        minVertex.z = prim.min(minVertex.z, T.z);
    }

    // Renvoie le sommet maximum du bounding box
    maxV->x = prim.max(maxV->x, maxVertex.x);
    maxV->y = prim.max(maxV->y, maxVertex.y);
    maxV->z = prim.max(maxV->z, maxVertex.z);

    // Renvoie le sommet minimum du bounding box
    minV->x = prim.min(minV->x, minVertex.x);
    minV->y = prim.min(minV->y, minVertex.y);
    minV->z = prim.min(minV->z, minVertex.z);
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * La fonction renvoie true s'il on a bien une intersection 
*/
bool Sphere::intersectionBB(ray r)
{
    /*
        Ce teste utilise l'algorithme suivant:
        Majercik, Alexander, Cyril Crassin, Peter Shirley, and Morgan McGuire. "A ray-box 
        intersection algorithm and efficient dynamic voxel rendering." Journal of Computer Graphics Techniques Vol 7, no. 3 (2018).
    */

    // Renvoie true si le rayon provient de l'interieur du box
    if( (r.x0.x > minVertex.x) && (r.x0.x < maxVertex.x) )
    {
        if( (r.x0.y > minVertex.y) && (r.x0.y < maxVertex.y) )
        {
            if( (r.x0.z > minVertex.z) && (r.x0.z < maxVertex.z) )
            {
                return true;
            }
        }
    }

    // Calculer l'intersection entre le rayon et les 6 côtés du box
    float t1 = (minVertex.x - r.x0.x)*r.n_inv.x;
    float t2 = (maxVertex.x - r.x0.x)*r.n_inv.x;
    float t3 = (minVertex.y - r.x0.y)*r.n_inv.y;
    float t4 = (maxVertex.y - r.x0.y)*r.n_inv.y;
    float t5 = (minVertex.z - r.x0.z)*r.n_inv.z;
    float t6 = (maxVertex.z - r.x0.z)*r.n_inv.z;

    // Calculer la composante max du tmin
    float tmin = prim.max(prim.max(prim.min(t1, t2), prim.min(t3, t4)), prim.min(t5, t6));

    // Calculer la composante min du tmax
    float tmax = prim.min(prim.min(prim.max(t1, t2), prim.max(t3, t4)), prim.max(t5, t6));

    // Si le rayon intersecte le bounding box (qui est en arrière)
    if (tmax < 0)
    {
        return false;
    }

    // Dans le cas contraire
    if (tmin > tmax)
    {
        return false;
    }

    return true;
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * I : L'informations sur le point d'intersection
 * La fonction renvoie true s'il on a bien une intersection
*/
bool Cylinder::intersecte(ray r, Inter *I) {

    // renvoie le rayon s'il n'y a pas d'intersection
    if(!intersectionBB(r))
    {
        return false;
    }

    // Calculer l'intersection du rayon avec la primitive
    Vec4 intersec = this->prim.intersect_cylinder(transfo, r);

    // S'il y a une intersection, on ajout un petit offset pour la stabilité
    if (intersec.w>-I->node->Epsilon)
    {
        // Mettre à jour I si le node est plus proche que celui contenu dans I
        if (intersec.w < I->a_min)
        {
            I->Dir = r.dir;
            I->node = this;
            I->a_min = intersec.w;
            I->Pos = Vec3(intersec);
        }
        return true;
    }
    return false;
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Cylinder::normal(const Vec3 &P) {
    return this->prim.normal_cylinder(transfo, P);
}

/*
 * minV : Le bord minimum du bounding box.
 * maxV : Le bord maximum du bounding box.
 */
void Cylinder::getBoundingBox(Vec3* minV, Vec3* maxV)
{
    // Stocker tous les sommets du cube dans le tableau vertex
    Vec3 vertex[8];
    for (int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++)
        {
            for (int k=0; k<2; k++)
            {
                vertex[i*4+j*2+k].x = i*2-1;
                vertex[i*4+j*2+k].y = j*2-1;
                vertex[i*4+j*2+k].z = k*2-1;
            }
        }
    }

    // Initialiser maxVertex avec moins l'infini
    maxVertex.x = -std::numeric_limits<float>::max();
    maxVertex.y = -std::numeric_limits<float>::max();
    maxVertex.z = -std::numeric_limits<float>::max();

    // Initialiser minVertex avec l'infini
    minVertex.x = std::numeric_limits<float>::max();
    minVertex.y = std::numeric_limits<float>::max();
    minVertex.z = std::numeric_limits<float>::max();

    // Parcourir tous les sommets du cube et mettre à jour le bounding box
    for (int i=0; i<8;i++)
    {
        // Transformer vertex[i]
        Vec4 T = transfo*Vec4(vertex[i], 1.0f);

        // Mettre à jour maxVertex
        maxVertex.x = prim.max(maxVertex.x, T.x);
        maxVertex.y = prim.max(maxVertex.y, T.y);
        maxVertex.z = prim.max(maxVertex.z, T.z);

        // Mettre à jour minVertex
        minVertex.x = prim.min(minVertex.x, T.x);
        minVertex.y = prim.min(minVertex.y, T.y);
        minVertex.z = prim.min(minVertex.z, T.z);
    }

    // Renvoie le sommet max du bounding box
    maxV->x = prim.max(maxV->x, maxVertex.x);
    maxV->y = prim.max(maxV->y, maxVertex.y);
    maxV->z = prim.max(maxV->z, maxVertex.z);

    // Renvoie le sommet min du bounding box
    minV->x = prim.min(minV->x, minVertex.x);
    minV->y = prim.min(minV->y, minVertex.y);
    minV->z = prim.min(minV->z, minVertex.z);
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * La fonction renvoie true s'il on a bien une intersection 
*/
bool Cylinder::intersectionBB(ray r)
{
    /*
        Ce test utilise l'algo dans
        Majercik, Alexander, Cyril Crassin, Peter Shirley, and Morgan McGuire. "A ray-box
        intersection algorithm and efficient dynamic voxel rendering." Journal of Computer Graphics
    */
    
    // Renvoie true s'il le rayon provient de l'intérieur du box
    if( (r.x0.x > minVertex.x) && (r.x0.x < maxVertex.x) )
    {
        if( (r.x0.y > minVertex.y) && (r.x0.y < maxVertex.y) )
        {
            if( (r.x0.z > minVertex.z) && (r.x0.z < maxVertex.z) )
            {
                return true;
            }
        }
    }

    // Calculer l'intersection du rayon avec chaque 6 côtes du box
    float t1 = (minVertex.x - r.x0.x)*r.n_inv.x;
    float t2 = (maxVertex.x - r.x0.x)*r.n_inv.x;
    float t3 = (minVertex.y - r.x0.y)*r.n_inv.y;
    float t4 = (maxVertex.y - r.x0.y)*r.n_inv.y;
    float t5 = (minVertex.z - r.x0.z)*r.n_inv.z;
    float t6 = (maxVertex.z - r.x0.z)*r.n_inv.z;

    // Calculer la composante max de tmin
    float tmin = prim.max(prim.max(prim.min(t1, t2), prim.min(t3, t4)), prim.min(t5, t6));

    // Calculer la composante min de tmax
    float tmax = prim.min(prim.min(prim.max(t1, t2), prim.max(t3, t4)), prim.max(t5, t6));

    // Si le rayon intersecte avec le box qui est en arrière
    if (tmax < 0)
    {
        return false;
    }

    // Dans le cas contraire
    if (tmin > tmax)
    {
        return false;
    }

    return true;
}

/*
 *  /////////////////////////////////////////
 *  ////////////////// BVH //////////////////
 *  /////////////////////////////////////////
 */

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Intersection point information
 */
void BVH::closestIntersection(ray r, Inter *I) {
  // AC
    /*
        On va parcourir le BVH, en testant l'intersection avec Origin et Dir des paramètres 
        et on cherche l'intersection la plus proche de l'origine du rayon 
        et on stocke l'information de l'intersection dans I
    */

    // Si le rayon n'intersecte pas le box
    if(intersectionBB(r) == false)
    {
        return;
    }

    // On vérifie l'intersection dans tous les nodes de BVH
    for (int i=0; i<nodes.size();i++)
    {
        nodes[i]->intersecte(r, I);
    }

    // On fait la même chose pour les enfants de BVH
    for (int i=0; i<children.size();i++)
    {
        children[i]->closestIntersection(r, I);
    }
}

/*
 * Origin: l'origine du rayon
 * Dir : la direction du rayon
 * sha : le coefficient de l'ombre / shadow coefficient (entre 0.0 et 1.0)
 * retour true s'il y a bien une intersection
*/
void BVH::intersecteShadow(ray r, float &sha) {
  // AC
  
    /* 
        On va parcourir le BVH, en testant l'intersection avec Origin et Dir des paramètres 
        et on cherche l'intersection avec les nodes intermediaires du BVH et après les Nodes
        Ici on cherche juste s'il y a d'intersection ou pas
    */

    /* 
        S'il y a une intersection alors le point d'origine est le shadow car il y a un autre objet
        entre lui et la source de la lumière
    */
    sha = 0;

    // Calculer la contribution de l'ombre (shadow contribution) de tous les nodes de BVH
    for (int i=0; i<nodes.size();i++)
    {
        Inter It;
        It.a_min = std::numeric_limits<float>::max();
        It.node = nullptr;

        // Vérifier si le rayon intersecte le node
        if (nodes[i]->intersecte(r, &It))
        {
            // Mettre un treshold pour la stabilité
            if (It.a_min > Node::Epsilon)
            {
                // Mettre à jour le coéfficient de l'ombre (shadow co-efficient) en ajoutant la valeur de transparence du node
                sha = sha + nodes[i]->transp;

                if (sha>1.0f)
                {
                    sha = 1.0f;
                    return;
                }
            }
        }
    }

    // Calculer le shadow contribution des enfant de BVH
    for (int i=0; i<children.size();i++)
    {
        // Vérifier si le rayon intersecte l'enfant
        children[i]->intersecteShadow(r, sha);

        
        if (sha>1.0f)
        {
            sha = 1.0f;
            return;
        }
    }

    if (sha>1.0f)
    {
        sha = 1.0f;
    }

    /* 
        Si l'objet est transparent, on va juste ajouter le coefficient de transparence et on continue à tester
        tant qu'on a pas arrivé à la valeur 1.0
        On met à jour le shadow coefficient la variable sha a pris comme paramètre 
        s'il y en a d'intersections, en utilisant le coefficient de transparence de l'objet en question
    */
}

/*
 * minV : Le bord minimum du bounding box.
 * maxV : Le bord maximum du bounding box.
*/
void BVH::getBoundingBox(Vec3* minV, Vec3* maxV)
{
    // Initialiser le sommet max avec moins l'infini
    maxVertex.x = -std::numeric_limits<float>::max();
    maxVertex.y = -std::numeric_limits<float>::max();
    maxVertex.z = -std::numeric_limits<float>::max();

    // Initialiser le sommet min avec l'infini
    minVertex.x = std::numeric_limits<float>::max();
    minVertex.y = std::numeric_limits<float>::max();
    minVertex.z = std::numeric_limits<float>::max();

    // Calculer les bounding boxes des nodes dans BVH
    for(int i=0; i<nodes.size();i++)
    {
        nodes[i]->getBoundingBox(&minVertex, &maxVertex);
    }

    // Calculer les bounding boxes des nodes des enfants de BVH
    for(int i=0; i<children.size();i++)
    {
        children[i]->getBoundingBox(&minVertex, &maxVertex);
    }

    // Renvoie le sommet max du bounding box
    maxV->x = max(maxV->x, maxVertex.x);
    maxV->y = max(maxV->y, maxVertex.y);
    maxV->z = max(maxV->z, maxVertex.z);

    // Renvoie le sommet min du bounding box
    minV->x = min(minV->x, minVertex.x);
    minV->y = min(minV->y, minVertex.y);
    minV->z = min(minV->z, minVertex.z);
}

/*
 * a : Première nombre à comparer
 * b : Deuxième nombre à comparer
 * renvoie true si a > b
*/
float BVH::max(float a, float b)
{
    if (a>b)
    {
        return a;
    }

    return b;
}

/*
 * a : Première nombre à comparer
 * b : Deuxième nombre à comparer
 * renvoie true si a < b
 */
float BVH::min(float a, float b)
{
    if (a<b)
    {
        return a;
    }

    return b;
}

/*
 * r : L'objet rayon qui contient l'origine, la direction ainsi que l'inverse
 * La fonction renvoie true s'il on a bien une intersection 
*/
bool BVH::intersectionBB(ray r) {

    /*
        Ce test utilise l'algo dans
        Majercik, Alexander, Cyril Crassin, Peter Shirley, and Morgan McGuire. "A ray-box
        intersection algorithm and efficient dynamic voxel rendering." Journal of Computer Graphics
    */

    // Renvoie true s'il le rayon provient de l'intérieur du box
    if( (r.x0.x > minVertex.x) && (r.x0.x < maxVertex.x) )
    {
        if( (r.x0.y > minVertex.y) && (r.x0.y < maxVertex.y) )
        {
            if( (r.x0.z > minVertex.z) && (r.x0.z < maxVertex.z) )
            {
                return true;
            }
        }
    }

    // Calculer l'intersection du rayon avec chaque 6 côtes du box
    float t1 = (minVertex.x - r.x0.x)*r.n_inv.x;
    float t2 = (maxVertex.x - r.x0.x)*r.n_inv.x;
    float t3 = (minVertex.y - r.x0.y)*r.n_inv.y;
    float t4 = (maxVertex.y - r.x0.y)*r.n_inv.y;
    float t5 = (minVertex.z - r.x0.z)*r.n_inv.z;
    float t6 = (maxVertex.z - r.x0.z)*r.n_inv.z;

    // Calculer la composante max du tmin
    float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));

    // Calculer la composante min du tmax
    float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    // Si le rayon intersecte avec le box qui est en arrière
    if (tmax < 0)
    {
        return false;
    }

    // Dans le cas contraire
    if (tmin > tmax)
    {
        return false;
    }

    return true;
}

const float Node::Epsilon = 0.001f;

/*
 *  /////////////////////////////////////////
 *  /////////////// RAY TRACER //////////////
 *  /////////////////////////////////////////
 */

RTracer::RTracer(SimpleViewer *s) : sv(s), depth(0) {}

void RTracer::add_cube_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                           float tr) {
  b->nodes.push_back(new Cube(m, color, spec, tr));
  // DC
}
void RTracer::add_sphere_bvh(BVH *b, const Mat4 &m, const Vec3 &color,
                             float spec, float tr) {

  b->nodes.push_back(new Sphere(m, color, spec, tr));
  // DC
}
void RTracer::add_cylinder_bvh(BVH *b, const Mat4 &m, const Vec3 &color,
                               float spec, float tr) {
  b->nodes.push_back(new Cylinder(m, color, spec, tr));
  // DC
}

/*
 * b           : L'objet BVH auquel on ajout l'Éponge de Menger 
 * m           : La matrice de transformation de l'éponge
 * colorCube   : La couleur des cubes de l'éponge
 * specCube    : La spécularité des cubes dans l'éponge
 * trCube      : La transparence des cubes dans l'éponge
 * colorSphere : La couleur des sphères dans l'éponge
 * specSphere  : La spécularité des sphères dans l'éponge 
 * trSphere    : La transparence des sphères dans l'éponge
 * Depth       : Nombre de récursions dans l'éponge
 */
void RTracer::add_menger_sponge_bvh(BVH *b, const Mat4 &m, const Vec3 &colorCube, float specCube,
                      float trCube, const Vec3 &colorSphere, float specSphere, float trSphere, int Depth)
{
    Mat4 mat = glm::translate(m, Vec3(-1/3.0f,-1/3.0f, -1/3.0f));
    Mat4 transX = glm::scale(mat, Vec3(1/3.0f, 1/3.0f, 1/3.0f));
    for(int i=0; i<3; i++)
    {
        Mat4 transY = transX;
        transX = glm::translate(transX, Vec3(1.0f, 0.0f, 0.0f));
        for(int j=0; j<3; j++)
        {
            Mat4 transZ = transY;
            transY = glm::translate(transY, Vec3(0.0f, 1.0f, 0.0f));
            for(int k=0; k<3; k++)
            {
                Mat4 trans = transZ;
                transZ = glm::translate(transZ, Vec3(0.0f, 0.0f, 1.0f));

                int score = (i==1) + (j==1) + (k==1);

                if (score == 3)
                {
                    trans = glm::scale(trans, Vec3(0.4f, 0.4f, 0.4f));
                    add_sphere_bvh(bvh, trans, colorSphere, specSphere, trSphere);
                }
                else if (score == 2)
                {
                    trans = glm::scale(trans, Vec3(0.3f, 0.3f, 0.3f));
                    add_sphere_bvh(bvh, trans, colorSphere, specSphere, trSphere);
                }
                else
                {
                    if (Depth == 0)
                    {
                        trans = glm::scale(trans, Vec3(0.45f, 0.45f, 0.45f));
                        add_cube_bvh(bvh, trans, colorCube, specCube, trCube);
                    }
                    else
                    {
                        trans = glm::scale(trans, Vec3(0.9f, 0.9f, 0.9f));
                        BVH* bvhChild = new BVH(trans);
                        add_menger_sponge_bvh(bvhChild, trans, colorCube, specCube, trCube, colorSphere, specSphere, trSphere, Depth-1);
                        b->children.push_back(bvhChild);
                    }
                }
            }
        }
    }
}

/*
 * a : Premièr nombre à comparer
 * b : Deuxième nombre à comparer
 * renvoie true si a > b
 */
float RTracer::max(float a, float b)
{
    if (a>b)
    {
        return a;
    }

    return b;
}

/*
 * Origin : L'origine du rayon
 * Dir    : La direction du rayon
 * rec    : Nombre de rebonds (rebounds)
 * renvoie la couleur de pixels
 */
Vec3 RTracer::ColorRayBVH(const Vec3 &Origin, const Vec3 &Dir, int rec) {

  Inter I;
  I.a_min = std::numeric_limits<float>::max();
  I.node = nullptr;

  /*
    L'objet rayon est calculé, en calculant aussi l'inverse de la direction
    utilisée dans le calcul d'intersections de box du rayon 
  */
  ray r = ray(Origin, Dir);

  // Trouver le noeud ayant l'intersection la plus proche avec le rayon
  bvh->closestIntersection(r, &I);

  // Si le rayon a une intersection avec un noeud
  if (I.node != nullptr)
  {
    // Calculer et normaliser le vecteur directeur
    const Vec3 &D = glm::normalize(I.Dir);

    // Position
    Vec3 Po = I.Pos;

    // Light direction
    // AC : calculer Lu, la direction de la lumière
    Vec3 Lu = glm::normalize(this->posLum - I.Pos);

    // Normal
    Vec3 No = I.node->normal(Po);

    // Shadow
    float sha = 0.0f;
    // AC : calculer le coefficient d'ombre
    ray rb = ray(Po, Lu);
    bvh->intersecteShadow(rb, sha);


    // Lambert (diffuse) BRDF
    float lambert = 0.0f;
    // AC : calculer la diffusion de l'objet avec -> (1 - sha) * col ∗ coeff ∗ max(0, No · Lu)
    lambert = (1-sha)*0.8*max(0.0, glm::dot(No, Lu)) + 0.2;

    // Ajouter un petit offset au point de l'intersection pour garder la stabilité
    Po = Po +I.node->Epsilon*No;

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // The rest is already completed, it is about other behaviors of light and objects: specularity, transparency and reflections
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // View direction for specular BRDF
    Vec3 inci = glm::normalize(D);
    Vec3 R = inci - 2.0f * glm::dot(No, inci) * No;
    Vec3 RR = glm::normalize(R);

    if (rec > 0) {
      Vec3 col2(0);
      if (I.node->spec > 0)
      {
          col2 = ColorRayBVH(Po, RR, rec - 1);
      }

      Vec3 col3(0);
      if (I.node->transp > 0)
        col3 = ColorRayBVH(Po, D, rec - 1);

      Vec3 final = (1.0f - I.node->transp) *
                       ((1.0f - I.node->spec) * I.node->Col * lambert +
                        I.node->spec * col2) +
                   I.node->transp * col3;

      // Specular BRDF
      final += (1.0f - sha) * I.node->spec *
               Vec3(std::pow(std::max(0.0f, glm::dot(RR, Lu)), 40.0f));

      return Color(final);
    }

    Vec3 final = I.node->Col * lambert;

    // Specular BRDF
    final += (1.0f - sha) * I.node->spec *
             Vec3(std::pow(std::max(0.0f, glm::dot(RR, Lu)), 40.0f));

    return Color(final);
  }

  return Vec3(0, 0, 0);
}

QLabel *RTracer::CalcImage(int depth) {
  std::cout << "Ray tracing calculation with " << depth << " rebounds - on "
            << this->nbthr << " threads ..." << std::endl;

  FILE* f = fopen("log.txt", "w");
  fclose(f);


  QImage im(this->sv->width(), this->sv->height(), QImage::Format_RGB32);
  Vec4 SHADER_POSLUM(300, 1000, 1000, 1);
  this->posLum =
      Vec3(glm::inverse(this->sv->getCurrentModelViewMatrix()) * SHADER_POSLUM);
  this->depth = depth;

  std::cout<<"Light at ("<<this->posLum.x<<","<<this->posLum.y<<","<<this->posLum.z<<")\n";


  int x, y;
  auto start = std::chrono::system_clock::now();
#pragma omp parallel for private(x) schedule(dynamic) num_threads(this->nbthr)
  for (y = 0; y < this->sv->height(); ++y) {
    auto start1 = std::chrono::system_clock::now();

    for (x = 0; x < this->sv->width(); ++x) {
      qglviewer::Vec Pq = this->sv->camera()->unprojectedCoordinatesOf(
          qglviewer::Vec(x, y, 0.0));
      qglviewer::Vec Qq = this->sv->camera()->unprojectedCoordinatesOf(
          qglviewer::Vec(x, y, 1.0));
      Vec3 P(Pq[0], Pq[1], Pq[2]);
      Vec3 D(Qq[0] - Pq[0], Qq[1] - Pq[1], Qq[2] - Pq[2]);
      Vec3 C = this->ColorRayBVH(P, D, this->depth);
      // ^ ou ColorRay(...) : launches a ray of origin P and direction D,
      // bouncing this-> depth times.
      im.setPixel(
          x, y,
          QColor(int(C.r * 255.0f), int(C.g * 255.0f), int(C.b * 255.0f))
              .rgb());
    }

    auto end1 = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end1 - start1;
    FILE* f = fopen("log.txt", "a");
    fprintf(f, "%f\n", elapsed.count());
    fclose(f);

  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;

  QLabel *label_img = new QLabel(this->sv);
  label_img->setWindowFlags(Qt::Window);
  label_img->setPixmap(QPixmap::fromImage(im));
  return label_img;
  // DC
}

void draw_prim_bvh(BVH *b) {
  for (const auto *n : b->nodes)
    n->draw_gl();
  for (auto *c : b->children)
    draw_prim_bvh(c);
  // DC
}
