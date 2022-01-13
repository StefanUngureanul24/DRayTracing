#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <vector>
#include <OGLRender/shaderprogramflat.h>
#include <matrices.h>

class ray {
    public:

    ray(const Vec3 &Origin, const Vec3 &Dir);
    Vec3 x0;     // L'origine du rayon
    Vec3 n_inv;  // L'inverse de la direction du rayon
    Vec3 dir;    // La direction du rayon
};

class Primitives
{
	Mat4 viewMatrix;
	Mat4 projectionMatrix;

	std::vector<Vec3> m_points;
	std::vector<int> m_indices_cube;
	std::vector<int> m_indices_sphere;
	std::vector<int> m_indices_cylinder;
	std::vector<int> m_indices_cone;

	void add_cylinder(int sides, float radius, std::vector<int>& indices);
	void add_cone(int sides, float radius, std::vector<int>& indices);
	void add_sphere(int sides, float radius, std::vector<int>& indices);

	ShaderProgramFlat* m_shader_flat;
	GLuint m_vao;
	GLuint m_vbo;

	GLuint m_ebo_cube;
	GLuint m_ebo_cylinder;
	GLuint m_ebo_cone;
	GLuint m_ebo_sphere;

    Vec4 intersectTriangle(ray r, const Vec3 &vertex1, const Vec3 &vertex2, const Vec3 &vertex3);
    bool inTriangle(const Vec3 &p, const Vec3 &vertex1, const Vec3 &vertex2, const Vec3 &vertex3);

public:
	Primitives();

    float max(float a, float b);
    float min(float a, float b);

	void gl_init();

	void set_matrices(const Mat4& view, const Mat4& projection);

	void draw_cube(const Mat4& transfo, const Vec3& color) const;

	void draw_cone(const Mat4& transfo, const Vec3& color) const;

	void draw_sphere(const Mat4& transfo, const Vec3& color) const;

	void draw_cylinder(const Mat4& transfo, const Vec3& color) const;

    Vec3 normal_cube(const Mat4 &transfo, const Vec3& p);
    Vec4 intersect_cube(const Mat4& transfo, ray r);

    Vec3 normal_cylinder(const Mat4 &transfo, const Vec3& P);
    Vec4 intersect_cylinder(const Mat4 &transfo, ray r);

    Vec3 normal_cone(const Mat4 &transfo, const Vec3& P);
    Vec4 intersect_cone(const Mat4 &transfo, ray r);

    Vec3 normal_sphere(const Mat4 &transfo, const Vec3& P);
    Vec4 intersect_sphere(const Mat4& transfo, ray r);

	inline const std::vector<Vec3>& getPoints() const { return m_points;}

	inline const std::vector<int>& getCubeIndices() const {return m_indices_cube;}

	inline const std::vector<int>& getConeIndices() const {return m_indices_cone;}

	inline const std::vector<int>& getCylinderIndices() const {return m_indices_cylinder;}

};

#endif // PRIMITIVES_H
