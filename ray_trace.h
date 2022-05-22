#pragma once
#include "lineTriangle3DIntersection.h"
#define M_PI 3.141592653589793
#define CUR_INFINITY 1e8
#define MAX_RAY_DEPTH 1
struct my_draw_color
{
	float r;
	float g;
	float b;
	my_draw_color operator +(const my_draw_color& input_color)
	{
		return my_draw_color{ r + input_color.r,g + input_color.g,b + input_color.b };
	}
	my_draw_color operator *(const float& cons)
	{
		return my_draw_color{ r * cons,g * cons,b * cons };
	}
};
//对给定的投影平面进行离散点采样(直接在投影面上采点会无法显示）
//——本函数要求与视点与世界坐标系下某条轴重合，且投影平面与某个坐标平面平行，上述坐标轴垂直于上述坐标平面
void samplepoint_sonprojectionplan(float left_x, float right_x, float bottom_y, float up_y, float
	nearplane_distance, float eye_z, std::map<my_3D_point_coord*, my_draw_color*>&
	render_vertices, unsigned& width, unsigned& height, int gap)
{
	float x_delt, y_delt;
	if (gap == 0) {
		x_delt = (right_x - left_x) / (ceil(right_x) - floor(left_x));
		y_delt = (up_y - bottom_y) / (ceil(up_y) - floor(bottom_y));
	}
	else {
		x_delt = gap;
		y_delt = gap;
	}
	float z_val = eye_z - nearplane_distance - 1;
	width = 0;
	height = 0;
	bool counted = false;
	for (float x_iter = left_x; x_iter <= right_x; x_iter += x_delt)
	{
		width++;
		for (float y_iter = bottom_y; y_iter <= up_y; y_iter += y_delt)
		{
			my_3D_point_coord* tempPoint_ptr = new my_3D_point_coord(x_iter, y_iter, z_val);
			my_draw_color* tempColor_ptr = new my_draw_color{ 0,0,0 };
			render_vertices.insert(pair<my_3D_point_coord*, my_draw_color*>(tempPoint_ptr, tempColor_ptr));
			if (counted == false) height++;
		}
		counted = true;
	}
}
//Blinn-Phong Reflection Model
my_draw_color calculate_direct_light_on_one_vertex_usingBPRM(my_3D_point_coord
	input_vertex, my_3Dvector vertex_normal,
	my_3D_point_coord eye_position, my_3D_point_coord light_position,
	float light_rgb_ambient[], float material_ambient_rgb_reflection[],
	float light_rgb_diffuse_specular[], float material_diffuse_rgb_reflection[],
	float material_specular_rgb_reflection[], float ns)
{
	//环境光反射能量
	float rval = light_rgb_ambient[0] * material_ambient_rgb_reflection[0];
	float gval = light_rgb_ambient[1] * material_ambient_rgb_reflection[1];
	float bval = light_rgb_ambient[2] * material_ambient_rgb_reflection[2];
	//漫反射光能量——需保证非零增长角度在0-90
	my_3Dvector pointTolight_vector(input_vertex, light_position);
	pointTolight_vector.normalized();
	float costheta = vertex_normal.dot(pointTolight_vector);
	costheta = (costheta > 0.0) ? costheta : 0.0;
	rval += light_rgb_diffuse_specular[0] * material_diffuse_rgb_reflection[0] * costheta;
	gval += light_rgb_diffuse_specular[1] * material_diffuse_rgb_reflection[1] * costheta;
	bval += light_rgb_diffuse_specular[2] * material_diffuse_rgb_reflection[2] * costheta;
	//镜面高光能量
	my_3Dvector pointToView_vector(input_vertex, eye_position);
	pointToView_vector.normalized();
	my_3Dvector half_vector = pointTolight_vector + pointToView_vector;
	half_vector = half_vector * 0.5;
	float costheta2 = vertex_normal.dot(half_vector);
	costheta2 = (costheta2 > 0.0) ? costheta2 : 0.0;
	float special_coeff = powf(costheta2, ns);
	rval += light_rgb_diffuse_specular[0] * material_specular_rgb_reflection[0] * special_coeff;
	gval += light_rgb_diffuse_specular[1] * material_specular_rgb_reflection[1] * special_coeff;
	bval += light_rgb_diffuse_specular[2] * material_specular_rgb_reflection[2] * special_coeff;
	my_draw_color output_color = { rval,gval,bval };
	return output_color;
}
/***************************************
* 计算折射光线方向的函数
* inpuray_dir为入射光向量 nhit为交点法向 refracted_dir为输出的折射光线方向
* ni_over_nt为所离开物体与被进入物体之间的反射率比值 用斯奈尔定律计算
******************************************/
bool get_refract_dir_my(const my_3Dvector& inpuray_dir, my_3Dvector& nhit, float ni_over_nt,
	my_3Dvector& refracted_dir)
{
	my_3Dvector uv(inpuray_dir.dx, inpuray_dir.dy, inpuray_dir.dz);
	uv.normalized();
	float dt = uv.dot(nhit);
	float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
	if (discriminant > 0)
	{
		refracted_dir = (uv - nhit * dt) * ni_over_nt - nhit * sqrt(discriminant);
		return true;
	}
	else
		return false;
}
//光线跟踪算法、、raydir需要是一个单位向量
my_draw_color one_ray_trace_my(my_3D_point_coord rayorig, my_3Dvector raydir, const
	std::vector<my_triangle_3DModel>& all_models, const int& depth,
	const my_3D_point_coord& eye_position, const my_3D_point_coord& light_position, float
	light_rgb_ambient[], float light_rgb_diffuse_specular[])//用于交点直接颜色计算
{
	//超过递归层数，返回000
	if (depth > MAX_RAY_DEPTH)
		return my_draw_color{ 0,0,0 };
	//计算模型（光线的交点）——三角网格遍历与光线求交——一定要注意起点在三角面上的问题
	float nearest_t = INFINITY; //射线参数方程中的t
	my_3Dvector nearestTrangleNormal;
	const my_triangle_3DModel* nearestModel = NULL;
	my_3D_line curRay(rayorig.x, rayorig.y, rayorig.z, raydir.dx, raydir.dy, raydir.dz);
	for (unsigned model_index = 0; model_index < all_models.size(); model_index++)
	{
		for (unsigned int tri_index = 0; tri_index < all_models[model_index].faceSets.size();
			tri_index++)
		{
			//从面实例中取出三个顶点
			int firstPointIndex =
				all_models[model_index].faceSets[tri_index].first_point_index;//取出顶点索引
			int secondPointIndex =
				all_models[model_index].faceSets[tri_index].second_point_index;
			int thirdPointIndex =
				all_models[model_index].faceSets[tri_index].third_point_index;
			my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex];//第一个顶点
			my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex];
			//第二个顶点
			my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //第三个顶点
			my_3D_triangle curTriangle = { p1, p2, p3, true };
			IntersectionBetweenLineAndTriangle newIntTest(curRay, curTriangle);
			//不仅有交点，还要求不能是出发点newIntTest.GetLineParameter() > 1e-3 
			if (newIntTest.Find() && newIntTest.GetLineParameter() > 0.002f &&
				(newIntTest.GetLineParameter() < nearest_t))
			{
				nearestModel = &all_models[model_index];
				nearestTrangleNormal = newIntTest.GetHitPointNormal();
				nearest_t = newIntTest.GetLineParameter();
			}
		}
	}
	if (!nearestModel)
		return depth == 0 ? my_draw_color{ 1,1,1 } : my_draw_color{ 0,0,0 };
	my_3Dvector added_valVec = raydir * nearest_t;
	my_3D_point_coord phit = rayorig.add(added_valVec.dx, added_valVec.dy, added_valVec.dz); // 获得交点
	my_3Dvector nhit = nearestTrangleNormal; // 获得交点处的法向
	nhit.normalized();
	//用Blinn-Phong Reflection Model计算交点颜色 ——只要有击中就算一下是否能被光源覆盖到，不管是直射、反射还是折射
	float ambient_rgb_reflection[3] = {
	nearestModel->material_ambient_rgb_reflection[0],
	nearestModel->material_ambient_rgb_reflection[1],
	nearestModel->material_ambient_rgb_reflection[2] };
	float diffuse_rgb_reflection[3] = {
	nearestModel->material_diffuse_rgb_reflection[0],
	nearestModel->material_diffuse_rgb_reflection[1],
	nearestModel->material_diffuse_rgb_reflection[2] };
	float specular_rgb_reflection[3] = {
	nearestModel->material_specular_rgb_reflection[0] ,
	nearestModel->material_specular_rgb_reflection[1] ,
	nearestModel->material_specular_rgb_reflection[2] };
	my_draw_color surface_directColor =
		calculate_direct_light_on_one_vertex_usingBPRM(phit, nhit, eye_position, light_position,
			light_rgb_ambient, ambient_rgb_reflection,
			light_rgb_diffuse_specular, diffuse_rgb_reflection,
			specular_rgb_reflection, nearestModel->ns);
	//若是内部点，说明光线在模型内部走动，则交到的法向反向，目前此种情况只考虑折射光继续，反射光终止
	bool inside = false;
	if (raydir.dot(nhit) > 0)
	{
		inside = true;
		nhit = nhit * -1;
	}
	//计算反射光贡献的颜色
	my_draw_color reflectionColor = { 0,0,0 };
	if (nearestModel->reflection > 0 && inside == false)
	{
		my_3Dvector refldir = raydir - nhit * 2 * raydir.dot(nhit);//计算反射光方向
		refldir.normalized();
		reflectionColor = one_ray_trace_my(phit, refldir, all_models, depth + 1, eye_position,
			light_position, light_rgb_ambient, light_rgb_diffuse_specular); //+ nhit*bias
		reflectionColor = reflectionColor * 0.5;//0.5 反射系数（能量减一半）
	}
	//计算折射光贡献的颜色
	my_draw_color refractionColor = { 0,0,0 };
	if (nearestModel->transparency > 0)
	{
		my_3Dvector refrdir = raydir;
		float ni_over_nt = inside ? nearestModel->transparency / 1.00029 : 1.00029 /
			nearestModel->transparency; //1.00029 为空气的折射率
		bool refrected = get_refract_dir_my(raydir, nhit, ni_over_nt, refrdir); //计算折射光贡献的颜色
		if (refrected)
		{
			refrdir.normalized();
			refractionColor = one_ray_trace_my(phit, refrdir, all_models, depth + 1,
				eye_position, light_position, light_rgb_ambient, light_rgb_diffuse_specular); //- nhit*bias
			refractionColor = refractionColor * 0.5; //0.5折射系数 能量减一半
		}
	}
	return surface_directColor + reflectionColor + refractionColor;
}