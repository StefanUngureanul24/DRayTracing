#include <cmath>
#include "primitives.h"


void Primitives::set_matrices(const Mat4& view, const Mat4& projection)
{
	viewMatrix = view;
	projectionMatrix = projection;
}

void Primitives::add_cylinder(int sides, float radius, std::vector<int>& indices)
{
	int beg = m_points.size();
	float a = M_PI/4.0f;
	float b = 2.0f*M_PI/sides;
	for (int i=0;i<sides;++i)
	{
		m_points.push_back(Vec3(radius*cos(a),radius*sin(a),-1.0f));
		m_points.push_back(Vec3(radius*cos(a),radius*sin(a),1.0f));
		a+=b;
	}

	int cb = m_points.size();
	m_points.push_back(Vec3(0.0,0.0,-1.0f));
	int ch = m_points.size();
	m_points.push_back(Vec3(0.0,0.0, 1.0f));

	int sides2 = 2*sides;

	for (int i=0;i<sides;++i)
	{
        // le côté
        indices.push_back(beg+2*i);
        indices.push_back(beg+(2*i+2)%sides2);
        indices.push_back(beg+(2*i+1)%sides2);

        // le côté
        indices.push_back(beg+(2*i+3)%sides2);
        indices.push_back(beg+(2*i+1)%sides2);
        indices.push_back(beg+(2*i+2)%sides2);

        // l'arrière
        indices.push_back(beg+(2*i+2)%sides2);
        indices.push_back(beg+2*i);
        indices.push_back(cb);

        // la face
        indices.push_back(beg+(2*i+1)%sides2);
        indices.push_back(beg+(2*i+3)%sides2);
        indices.push_back(ch);
	}
}

void Primitives::add_cone(int sides, float radius, std::vector<int>& indices)
{
	int beg = m_points.size();
	float a = 0.0f;
	float b = 2.0f*M_PI/sides;
	for (int i=0;i<sides;++i)
	{
		m_points.push_back(Vec3(radius*cos(a),radius*sin(a),-0.5));
		a+=b;
	}

	int cb = m_points.size();
	m_points.push_back(Vec3(0.0,0.0,-0.5));
	int ch = m_points.size();
	m_points.push_back(Vec3(0.0,0.0, 0.5));

	for (int i=0;i<sides;++i)
	{
		indices.push_back(beg+i);
		indices.push_back(beg+(i+1)%sides);
		indices.push_back(ch);

		indices.push_back(beg+(i+1)%sides);
		indices.push_back(beg+i);
		indices.push_back(cb);
	}
}

void Primitives::add_sphere(int sides, float radius, std::vector<int>& indices)
{
    int beg = m_points.size();

    int nbPara = sides;
    int nbMeri = sides;

    float a1 = (180.0/(nbPara+1))*M_PI/180.0;
    float a2 = (360.0/nbMeri)*M_PI/180.0;

    // les paralleles
    for (int i= 0; i< nbPara; ++i)
    {
        float angle = -M_PI/2.0 + a1*(i+1);
        float z = radius*sin(angle);
        float rad = radius*cos(angle);

        for  (int j=0; j< nbMeri; ++j)
        {
            m_points.push_back(Vec3(rad*cos(a2*j), rad*sin(a2*j),z));
        }
    }
    // les poles
    m_points.push_back(Vec3(0.0,0.0,-radius));
    m_points.push_back(Vec3(0.0,0.0, radius));

    // triangles
    for (int i= 0; i< (nbPara-1); ++i)
    {
        for  (int j=0; j< nbMeri; ++j)
        {
            indices.push_back(beg+nbMeri*i+j);
            indices.push_back(beg+nbMeri*i+(j+1)%nbMeri);
            indices.push_back(beg+nbMeri*(i+1)+(j+1)%nbMeri);
            indices.push_back(beg+nbMeri*((i+1))+(j+1)%nbMeri);
            indices.push_back(beg+nbMeri*((i+1))+j);
            indices.push_back(beg+nbMeri*i+j);
        }
    }
    // poles
    for  (int j=0; j< nbMeri; ++j)
    {
        indices.push_back(beg+nbMeri*nbPara);
        indices.push_back(beg+(j+1)%nbMeri);
        indices.push_back(beg+j);
    }
    for  (int j=0; j< nbMeri; ++j)
    {
        indices.push_back(beg+nbMeri*nbPara+1);
        indices.push_back(beg+nbMeri*(nbPara-1)+j);
        indices.push_back(beg+nbMeri*(nbPara-1)+(j+1)%nbMeri);
    }
}

/*
 * p : Le point à tester
 * vertex1 : Le premièr sommet du triangle
 * vertex2 : Le deuxième sommet du triangle
 * vertex3 : Le troisième sommet du triangle
 * renvoyer true si p se trouve dans le triangle formé par ces 3 sommets
*/
bool Primitives::inTriangle(const Vec3 &p, const Vec3 &vertex1, const Vec3 &vertex2, const Vec3 &vertex3)
{
    /* 
        Le point est à l'intérieur du triangle si les points qui se croisent ont la même direction rélative au plan
        Soit les points sont tous à l'intérieur, soit tous à l'exterieur
    */
    
    // Deplacer le triangle pour que le point devient l'origine du triangle
    Vec3 a = vertex1 - p;
    Vec3 b = vertex2 - p;
    Vec3 c = vertex3 - p;

    /*
        Calcul de normes des vecteurs des triangles:
        * u = normal of PBC
        * v = normal of PCA
        * w = normal of PAB
    */
    Vec3 u = glm::cross(b, c);
    Vec3 v = glm::cross(c, a);
    Vec3 w = glm::cross(a, b);

    /* 
        Tester si les normes ont la même direction:
        * true si oui
        * false si non
    */
    if (glm::dot(u, v) < 0.00f) {
        return false;
    }
    if (glm::dot(u, w) <0.00f) {
        return false;
    }

    return true;
}

/*
 * transfo : la matrice de transformation de la primitive du cube
 * p : le point où la norme doit être calculée
 * renvoie la norme de la primitive
*/
Vec3 Primitives::normal_cube(const Mat4 &transfo, const Vec3& p)
{
    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_cube.size(); i=i+3)
    {
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cube[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cube[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cube[i+2]], 1.0f));

        // Vérifier si le point est sur le triangle
        if (inTriangle(p, a, b, c))
        {
            Vec3 V = glm::normalize(b - a);
            Vec3 W = glm::normalize(c - a);

            // Calculer la norme            
            Vec3 N = glm::normalize(glm::cross(V,W));
            return N;
        }
    }

    return Vec3(0.0f, 0.0f, 0.0f);
}

/*
 * transfo : la matrice de transformation de la primitive de la sphère
 * p : le point où la norme doit être calculée
 * renvoie la norme de la primitive
*/
Vec3 Primitives::normal_sphere(const Mat4 &transfo, const Vec3& p)
{
    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_sphere.size(); i=i+3)
    {
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_sphere[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_sphere[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_sphere[i+2]], 1.0f));

        // Vérifier si le point est sur le triangle
        if (inTriangle(p, a, b, c))
        {
            Vec3 V = glm::normalize(b - a);
            Vec3 W = glm::normalize(c - a);

            // Calculer la norme
            Vec3 N = glm::normalize(glm::cross(V,W));
            return N;
        }
    }

    return Vec3(0.0f, 0.0f, 0.0f);
}

/*
 * transfo : la matrice de transformation de la primitive du cylindre
 * p : le point où la norme doit être calculée
 * renvoie la norme de la primitive
*/
Vec3 Primitives::normal_cylinder(const Mat4 &transfo, const Vec3& p)
{
    // iterate through all triangles in the primitive
    for(int i=0; i<m_indices_cylinder.size(); i=i+3)
    {
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i+2]], 1.0f));

        // Vérifier si le point est sur le triangle
        if (inTriangle(p, a, b, c))
        {
            Vec3 V = b - a;
            Vec3 W = c - a;

            // Calculer la norme
            Vec3 N = glm::normalize(glm::cross(V,W));
            return N;
        }
    }

    return Vec3(0.0f, 0.0f, 0.0f);
}

/*
 * transfo : la matrice de transformation de la primitive du cône
 * p : le point où la norme doit être calculée
 * renvoie la norme de la primitive
*/
Vec3 Primitives::normal_cone(const Mat4 &transfo, const Vec3& p)
{
    // iterate through all triangles in the primitive
    for(int i=0; i<m_indices_cone.size(); i=i+3)
    {
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cone[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cone[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cone[i+2]], 1.0f));

        // Vérifier si le point est sur le triangle
        if (inTriangle(p, a, b, c))
        {
            Vec3 V = b - a;
            Vec3 W = c - a;

            // Calculer la norme
            Vec3 N = glm::normalize(glm::cross(V,W));
            return N;
        }
    }

    return Vec3(0.0f, 0.0f, 0.0f);
}

/*
 * a : Le premièr nombre à comparer
 * b : Le deuxième nombre à comparer
 * renvoie true si a > b, false sinon
*/
float Primitives::max(float a, float b)
{
    if (a>b)
    {
        return a;
    }

    return b;
}

/*
 * a : Le premièr nombre à comparer
 * b : Le deuxième nombre à comparer
 * renvoie true si a < b, false sinon
 */
float Primitives::min(float a, float b)
{
    if (a<=b)
    {
        return a;
    }

    return b;
}

ray::ray(const Vec3 &Origin, const Vec3 &Dir)
{
    this->x0 = Origin;
    this->n_inv.x = 1.0f/Dir.x;
    this->n_inv.y = 1.0f/Dir.y;
    this->n_inv.z = 1.0f/Dir.z;
    this->dir = glm::normalize(Dir);
}

/*                                  
 * r : L'objet rayon - il contient l'origine du rayon, la direction, ainsi que l'inverse
 * vertex1 : Le première sommet du triangle 
 * vertex2 : Le deuxième sommet du triangle
 * vertex3 : Le troisième sommet du triangle
 * La fonction renvoie le point d'intersection concaténé avec la distance d'intersection
*/
Vec4 Primitives::intersectTriangle(ray r, const Vec3 &vertex1, const Vec3 &vertex2, const Vec3 &vertex3)
{
    /*
        Vérifier si le rayon et le triangle s'intersecte en utilisant l'algo de Tomas Moller et Ben Trumbore
    */

    // Calcul des vecteurs avec les bords du triangle
    Vec3 edge1 = vertex2-vertex1;
    Vec3 edge2 = vertex3-vertex1;

    // Calcul du déterminant
    Vec3 directionCrossEdge2 = glm::cross(r.dir, edge2);
    float determinant = glm::dot(edge1, directionCrossEdge2);

    // Si le rayon et le triangle sont parallèles alors il n'y a pas de collision
    if (determinant > -0.0000001 && determinant < 0.0000001)
    {
        return Vec4(0.0f, 0.0f, 0.0f, -1.0f);
    }

    float inverseDeterminant = 1.0f / determinant;

    // Calcul du paramètre U du point de l'intersection 
    Vec3 distanceVector = r.x0 - vertex1;

    float triangleU = glm::dot(distanceVector, directionCrossEdge2);
    triangleU *= inverseDeterminant;

    // En s'assurant que ça se passe à l'interieur du triangle
    if (triangleU < 0 || triangleU > 1)
    {
        return Vec4(0.0f, 0.0f, 0.0f, -1.0f);
    }

    // Calcul du paramètre V du point de l'intersection
    Vec3 distanceCrossEdge1 = glm::cross(distanceVector, edge1);

    float triangleV = glm::dot(r.dir, distanceCrossEdge1);
    triangleV *= inverseDeterminant;

    // En s'assurant que ça se passe à l'interieur du triangle
    if (triangleV < 0 || triangleU + triangleV > 1)
    {
        return Vec4(0.0f, 0.0f, 0.0f, -1.0f);
    }

    // Calcul de la distance entre le rayon et le triangle
    float rayDistance;
    rayDistance = glm::dot(edge2, distanceCrossEdge1);
    rayDistance *= inverseDeterminant;

    // Vérifier si le triangle est derrière l'origine du rayon
    if (rayDistance < -1.0)
    {
        return Vec4(0.0f, 0.0f, 0.0f, -1.0f);
    }

    // Calcul du point de l'intersection à partir des coordonnées barycentriques
    Vec3 Point = (1-triangleU-triangleV)*vertex1 + triangleU*vertex2 + triangleV*vertex3;

    // Ça renvoie le point d'intersection concaténé avec la distance de l'intersection
    return Vec4(Point.x, Point.y, Point.z, rayDistance);
}

/*
 * transfo: La matrice de transformation de la primitive du cube
 * r: Le rayon à tester
 * La fonction renvoie le point de l'intersection concaténé avec la distance de l'intersection 
*/
Vec4 Primitives::intersect_cube(const Mat4 &transfo, ray r)
{
    Vec4 dist;

    // La distance de l'intersection est initialisé à l'infini
    dist.w = std::numeric_limits<float>::max();
    int R = -1;

    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_cube.size(); i=i+3)
    {
        // Transformer les sommets du triangle
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cube[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cube[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cube[i+2]], 1.0f));

        // Tester s'il y a de l'intersection avec le rayon
        Vec4 res = intersectTriangle(r, a, b, c);

        // On ajout un petit offset pour stabiliser
        if (res.w>-0.001)
        {
            // Si le triangle est plus proche que le point d'intersection, on le met à jour
            if(res.w<dist.w)
            {
                dist = res;
                R = 1;
            }
        }
    }

    // S'il n'y a pas d'intersection on renvoie la distance négative
    if(R<0)
    {
        dist.w = -1;
    }

    return dist;
}

/*
 * transfo: La matrice de transformation de la primitive de la sphère
 * r: Le rayon à tester
 * La fonction renvoie le point de l'intersection concaténé avec la distance de l'intersection 
*/
Vec4 Primitives::intersect_sphere(const Mat4 &transfo, ray r)
{
    Vec4 dist;

    // La distance de l'intersection est initialisé à l'infini
    dist.w = std::numeric_limits<float>::max();
    int R = -1;

    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_sphere.size(); i=i+3)
    {
        // Transformer les sommets du triangle
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_sphere[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_sphere[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_sphere[i+2]], 1.0f));

        // Tester s'il y a de l'intersection avec le rayon
        Vec4 res = intersectTriangle(r, a, b, c);

        // On ajout un petit offset pour stabiliser
        if (res.w>-0.001)
        {
            // Si le triangle est plus proche que le point d'intersection, on le met à jour
            if(res.w<dist.w)
            {
                dist = res;
                R = 1;
            }
        }
    }

    // S'il n'y a pas d'intersection on renvoie la distance négative
    if(R<0)
    {
        dist.w = -1;
    }

    return dist;
}

/*
 * transfo: La matrice de transformation de la primitive du cylindre
 * r: Le rayon à tester
 * La fonction renvoie le point de l'intersection concaténé avec la distance de l'intersection 
*/
Vec4 Primitives::intersect_cylinder(const Mat4 &transfo, ray r)
{
    Vec4 dist;

    // La distance de l'intersection est initialisé à l'infini
    dist.w = std::numeric_limits<float>::max();
    int R = -1;

    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_cylinder.size(); i=i+3)
    {
        // Transformer les sommets du triangle
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cylinder[i+2]], 1.0f));

        // Tester s'il y a de l'intersection avec le rayon
        Vec4 res = intersectTriangle(r, a, b, c);

        // On ajout un petit offset pour stabiliser
        if (res.w>-0.001)
        {
            // Si le triangle est plus proche que le point d'intersection, on le met à jour
            if(res.w<dist.w)
            {
                dist = res;
                R = 1;
            }
        }
    }

    // S'il n'y a pas d'intersection on renvoie la distance négative
    if(R<0)
    {
        dist.w = -1;
    }

    return dist;
}

/*
 * transfo: La matrice de transformation de la primitive du cône
 * r: Le rayon à tester
 * La fonction renvoie le point de l'intersection concaténé avec la distance de l'intersection 
*/
Vec4 Primitives::intersect_cone(const Mat4 &transfo, ray r)
{
    Vec4 dist;

    // La distance de l'intersection est initialisé à l'infini
    dist.w = std::numeric_limits<float>::max();
    int R = -1;

    // Parcourir les triangles de la primitive
    for(int i=0; i<m_indices_cone.size(); i=i+3)
    {
        // Transformer les sommets du triangle
        Vec3 a = Vec3(transfo*Vec4(m_points[m_indices_cone[i]], 1.0f));
        Vec3 b = Vec3(transfo*Vec4(m_points[m_indices_cone[i+1]], 1.0f));
        Vec3 c = Vec3(transfo*Vec4(m_points[m_indices_cone[i+2]], 1.0f));

        // Tester s'il y a de l'intersection avec le rayon
        Vec4 res = intersectTriangle(r, a, b, c);

        // On ajout un petit offset pour stabiliser
        if (res.w>-0.001)
        {
            // Si le triangle est plus proche que le point d'intersection, on le met à jour
            if(res.w<dist.w)
            {
                dist = res;
                R = 1;
            }
        }
    }

    // S'il n'y a pas d'intersection on renvoie la distance négative
    if(R<0)
    {
        dist.w = -1;
    }

    return dist;
}

void Primitives::gl_init()
{
	m_shader_flat = new ShaderProgramFlat();

	//VBO
	glGenBuffers(1, &m_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	const std::vector<Vec3>& points = getPoints();
	glBufferData(GL_ARRAY_BUFFER, 3 * points.size() * sizeof(GLfloat), &(points.front()[0]), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//VAO
	glGenVertexArrays(1, &m_vao);
	glBindVertexArray(m_vao);
	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	glEnableVertexAttribArray(m_shader_flat->idOfVertexAttribute);
	glVertexAttribPointer(m_shader_flat->idOfVertexAttribute, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glBindVertexArray(0);

	//EBO indices
	glGenBuffers(1, &m_ebo_cube);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_cube);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,m_indices_cube.size() * sizeof(int), &(m_indices_cube.front()), GL_STATIC_DRAW);

	glGenBuffers(1, &m_ebo_cone);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_cone);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,m_indices_cone.size() * sizeof(int), &(m_indices_cone.front()), GL_STATIC_DRAW);

	glGenBuffers(1, &m_ebo_cylinder);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_cylinder);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,m_indices_cylinder.size() * sizeof(int), &(m_indices_cylinder.front()), GL_STATIC_DRAW);

	glGenBuffers(1, &m_ebo_sphere);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo_sphere);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,m_indices_sphere.size() * sizeof(int), &(m_indices_sphere.front()), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

void Primitives::draw_cube(const Mat4& transfo, const Vec3& color) const
{
	m_shader_flat->startUseProgram();

	m_shader_flat->sendViewMatrix(viewMatrix*transfo);
	m_shader_flat->sendProjectionMatrix(projectionMatrix);

	glUniform3fv(m_shader_flat->idOfColorUniform, 1, glm::value_ptr(color));
	glUniform3fv(m_shader_flat->idOfBColorUniform, 1, glm::value_ptr(color));

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,m_ebo_cube);
	glDrawElements(GL_TRIANGLES, m_indices_cube.size(),GL_UNSIGNED_INT,0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	glBindVertexArray(0);

	m_shader_flat->stopUseProgram();
}

void Primitives::draw_cylinder(const Mat4& transfo, const Vec3& color) const
{
	m_shader_flat->startUseProgram();

	m_shader_flat->sendViewMatrix(viewMatrix*transfo);
	m_shader_flat->sendProjectionMatrix(projectionMatrix);

	glUniform3fv(m_shader_flat->idOfColorUniform, 1, glm::value_ptr(color));
	glUniform3fv(m_shader_flat->idOfBColorUniform, 1, glm::value_ptr(color));

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,m_ebo_cylinder);
	glDrawElements(GL_TRIANGLES, m_indices_cylinder.size(),GL_UNSIGNED_INT,0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	glBindVertexArray(0);

	m_shader_flat->stopUseProgram();
}

void Primitives::draw_cone(const Mat4& transfo, const Vec3& color) const
{
	m_shader_flat->startUseProgram();

	m_shader_flat->sendViewMatrix(viewMatrix*transfo);
	m_shader_flat->sendProjectionMatrix(projectionMatrix);

	glUniform3fv(m_shader_flat->idOfColorUniform, 1, glm::value_ptr(color));
	glUniform3fv(m_shader_flat->idOfBColorUniform, 1, glm::value_ptr(color));

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,m_ebo_cone);
	glDrawElements(GL_TRIANGLES, m_indices_cone.size(),GL_UNSIGNED_INT,0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	glBindVertexArray(0);

	m_shader_flat->stopUseProgram();
}

void Primitives::draw_sphere(const Mat4& transfo, const Vec3& color) const
{
	m_shader_flat->startUseProgram();

	m_shader_flat->sendViewMatrix(viewMatrix*transfo);
	m_shader_flat->sendProjectionMatrix(projectionMatrix);

	glUniform3fv(m_shader_flat->idOfColorUniform, 1, glm::value_ptr(color));
	glUniform3fv(m_shader_flat->idOfBColorUniform, 1, glm::value_ptr(color));

	glBindVertexArray(m_vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,m_ebo_sphere);
	glDrawElements(GL_TRIANGLES, m_indices_sphere.size(),GL_UNSIGNED_INT,0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
	glBindVertexArray(0);

	m_shader_flat->stopUseProgram();
}


Primitives::Primitives()
{
    add_cylinder(64, 1.0f, m_indices_cylinder);
	add_cone(32,1.0f,m_indices_cone);
	add_cylinder(4,1.414,m_indices_cube);
	add_sphere(32,1.0f,m_indices_sphere);
}
