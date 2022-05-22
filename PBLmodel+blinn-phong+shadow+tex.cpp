#define _CRT_SECURE_NO_WARNINGS
#define GLUT_DISABLE_ATEXIT_HACK
#include "ObjLoader.h"
#include "lineTriangle3DIntersection.h"
#include "ray_trace.h"
#include "GLUT.H"
#include <math.h>

int texenable = 0;
GLuint texGround0;
GLuint texGround1;
GLuint texGround3;
GLuint texGround5;
#define BMP_Header_Length 54									//图像数据在内存块中的偏移量

// 判断整数是不是2的整数次幂
int power_of_two(int n) {
	if (n <= 0)
		return 0;
	return (n & (n - 1)) == 0;
}

// 读取一个BMP文件作为纹理如果失败，返回0，如果成功，返回纹理编号
GLuint load_texture(const char* file_name) {
	GLint width, height, total_bytes;
	GLubyte* pixels = 0;
	GLuint last_texture_ID = 0, texture_ID = 0;

	// 如果失败
	FILE* pFile = fopen(file_name, "rb");
	if (pFile == 0)
		return 0;

	// 读取图的宽度和高度
	fseek(pFile, 0x0012, SEEK_SET);
	fread(&width, 4, 1, pFile);
	fread(&height, 4, 1, pFile);
	fseek(pFile, BMP_Header_Length, SEEK_SET);

	// 计算每行像素所占字节数，并根据此数据计算总像素字节数
	{
		GLint line_bytes = width * 3;
		while (line_bytes % 4 != 0)
			++line_bytes;
		total_bytes = line_bytes * height;
	}

	// 根据总像素字节数分配内存
	pixels = (GLubyte*)malloc(total_bytes);
	if (pixels == 0) {
		fclose(pFile);
		return 0;
	}

	// 读取像素数据
	if (fread(pixels, total_bytes, 1, pFile) <= 0) {
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// 如果图的宽度和高度不是2的整数次方，缩放
	// 若图像宽高超过了OpenGL规定的最大值，缩放
	{
		GLint max;
		glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max);
		if (!power_of_two(width)
			|| !power_of_two(height)
			|| width > max
			|| height > max) {
			const GLint new_width = 256;
			const GLint new_height = 256; // 规定缩放后新的大小为边长的正方形
			GLint new_line_bytes, new_total_bytes;
			GLubyte* new_pixels = 0;

			// 计算每行需要的字节数和总字节数
			new_line_bytes = new_width * 3;
			while (new_line_bytes % 4 != 0)
				++new_line_bytes;
			new_total_bytes = new_line_bytes * new_height;

			// 分配内存
			new_pixels = (GLubyte*)malloc(new_total_bytes);
			if (new_pixels == 0) {
				free(pixels);
				fclose(pFile);
				return 0;
			}

			// 进行像素缩放
			gluScaleImage(GL_RGB, width, height, GL_UNSIGNED_BYTE, pixels, new_width, new_height, GL_UNSIGNED_BYTE, new_pixels);

			// 释放原来的像素数据，把pixels指向新的像素数据，并重新设置width和height
			free(pixels);
			pixels = new_pixels;
			width = new_width;
			height = new_height;
		}
	}

	// 分配一个新的纹理编号
	glGenTextures(1, &texture_ID);
	if (texture_ID == 0) {
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// 绑定新的纹理，载入纹理并设置纹理参数
	// 在绑定前，先获得原来绑定的纹理编号，以便在最后进行恢复
	GLint lastTextureID = last_texture_ID;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &lastTextureID);
	glBindTexture(GL_TEXTURE_2D, texture_ID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);
	glBindTexture(GL_TEXTURE_2D, lastTextureID);  //恢复之前的纹理绑定
	free(pixels);
	return texture_ID;
}



#define PI 3.14159265
#define nearplane_width 400										//视景体宽度
#define nearplane_height 400									//视景体高度
int nearplane_distance = 300;									//视景体近平面与视点距离
int farplane_distance = nearplane_distance + 1800;				//视景体远平面与视点距离

float theta = 5;												//视点旋转角度(角度制)
float delta = 10;												//位移量
float ratio = 0.5;												//投影面大小比例
float diffuse = 0.8;											//漫反射系数

my_3Dvector eyesight;											//朝向面的视线
float rotateoffset = 0;											//旋转偏移
my_3D_point_coord cameraoffset(0, 0, 0);						//摄像机平移偏移
my_3D_point_coord lightoffset(0, 0, 0);							//灯1平移偏移
my_3D_point_coord light1offset(0, 0, 0);						//灯2平移偏移
my_3D_point_coord original(0, 0, 0);							//原点位置
my_3D_point_coord eye_default_position(0, 0, 500);				//视点默认位置
my_3D_point_coord light_default_position(100.0, 50.0, 0);		//灯1默认位置
my_3D_point_coord light1_default_position(-100.0, 0, 10.0);		//灯2默认位置
my_3D_point_coord eye_position = { eye_default_position.x,
								  eye_default_position.y,
								  eye_default_position.z };		//视点位置
my_3D_point_coord light_position(light_default_position.x,
								light_default_position.y,
								light_default_position.z);		//光源1位置
my_3D_point_coord light1_position(light1_default_position.x,
								light1_default_position.y,
								light1_default_position.z);		//光源2位置

float light_rgb_ambient[] = { 0.9, 0.9, 0.9 };
float light1_rgb_ambient[] = { 0.5, 0.05, 0.05 };
float light_rgb_diffuse_specular[] = { diffuse, diffuse, diffuse };
float light1_rgb_diffuse_specular[] = { diffuse, diffuse, diffuse };

int selected = 1;												//已选中物体						//	1 2 3
float renderoffset = 0;											//纵向微调 — 解决光追蜡笔效果		//	4
int sighttype = 0;												//视线是否跟随物体					//	5
int drawaxis = 0;												//是否画模型中轴和原点				//	6

int selectedrgb = 0;											//选中的光源rgb — 0r 1g 2b 3亮度	//	, . /  [\]
int zbuffer = 0;												//自带zbuffer						//	z
bool open_light = true;											//灯1开关							//	o
bool open_light1 = false;										//灯2开关							//	o
bool rendered = false;											//光线追踪							//	r
int showshadow = 0;												//自写zbuffer是否显示阴影			//	l
int skipmodel[7] = {1, 1, 1, 1, 1, 1, 1};						//开关模型							//	b n m - g h j - t

std::vector< my_triangle_3DModel> all_models;					//场景中所有模型
std::map<my_3D_point_coord*, my_draw_color*> render_vertices;	//最终需要绘制的点以及采样点

unsigned image_w, image_h;										//光追投影面尺寸
int gap = 0;													//光追采样密度
void printManual();




#define FAREST 9999
int my_zbuffer = 0;												//自写zbuffer						//	x
float Zbuffer[nearplane_width + 1][nearplane_width + 1];
my_draw_color Framebuffer[nearplane_width + 1][nearplane_width + 1];
int shadowbuffer[nearplane_width + 1][nearplane_width + 1];

// 判断点是否在投影面内
int inprojection(my_3D_point_coord p, int width, int height) {
	if (p.x >= -width / 2 && p.x <= width / 2 && p.y >= -height / 2 && p.y <= height / 2)
		return 1;
	return 0;
}

// 判断点是否在投影面上的三角形内
int inface(my_3D_point_coord p, my_3D_point_coord a, my_3D_point_coord b, my_3D_point_coord c) {
	bool flag = false;
	//扫描线法，每当和边相交，就判断求反
	if ((((a.y <= p.y) && (p.y < b.y)) ||
		((b.y <= p.y) && (p.y < a.y))) &&
		(p.x < (a.x - b.x) * (p.y - b.y) / (a.y - b.y) + b.x)) {
		flag = !flag;
	}
	if ((((a.y <= p.y) && (p.y < c.y)) ||
		((c.y <= p.y) && (p.y < a.y))) &&
		(p.x < (a.x - c.x) * (p.y - c.y) / (a.y - c.y) + c.x)) {
		flag = !flag;
	}
	if ((((c.y <= p.y) && (p.y < b.y)) ||
		((b.y <= p.y) && (p.y < c.y))) &&
		(p.x < (c.x - b.x) * (p.y - b.y) / (c.y - b.y) + b.x)) {
		flag = !flag;
	}
	return flag;
}

// 判断点是否在三维空间中的三角面内
int inface3D(my_3D_point_coord p, my_3D_point_coord a, my_3D_point_coord b, my_3D_point_coord c) {
	//面积法
	//三个小三角形面积 1/2*a*b*sin<a,b>
	my_3Dvector A, B, C, M, N, P, n1, n2, n3;
	A = my_3Dvector(p, a);
	B = my_3Dvector(p, b);
	C = my_3Dvector(p, c);
	n1 = A; n1.normalized();
	n2 = B; n2.normalized();
	n3 = C; n3.normalized();
	float S1, S2, S3, angle1, angle2, angle3;
	angle1 = acosf(n1.dot(n2));
	angle2 = acosf(n1.dot(n3));
	angle3 = acosf(n2.dot(n3));
	S1 = 0.5 * A.len * B.len * sinf(angle1);
	S2 = 0.5 * A.len * C.len * sinf(angle2);
	S3 = 0.5 * B.len * C.len * sinf(angle3);

	//一个大三角形面积
	my_3Dvector Q,R;
	M = my_3Dvector(a, b);
	N = my_3Dvector(a, c);
	Q = M; Q.normalized();
	R = N; R.normalized();
	float S, angle;
	angle = acosf(Q.dot(R));
	S = 0.5 * M.len * N.len * sinf(angle);

	//若近似相等 硬阴影部分
	if (fabsf(S1 + S2 + S3 - S) < 0.015)
		return 1;
	////若近似程度降低 软边缘部分
	//else if (fabsf(S1 + S2 + S3 - S) >= 0.5 && fabsf(S1 + S2 + S3 - S) < 100)
	//	return 2;
	return 0;
	/*//叉乘法
	my_3Dvector N1, N2, N3;
	my_3Dvector A, B;

	N1 = my_3Dvector(a, p); 
	N2 = my_3Dvector(a, b); 
	N3 = my_3Dvector(a, c); 
	A = N1.cross(N2);
	B = N1.cross(N3);
	if (A.dot(B) > 0)
		return 0;

	N1 = my_3Dvector(b, p);
	N2 = my_3Dvector(b, a);
	N3 = my_3Dvector(b, c);
	A = N1.cross(N2);
	B = N1.cross(N3);
	if (A.dot(B) > 0)
		return 0;

	N1 = my_3Dvector(c, p);
	N2 = my_3Dvector(c, a);
	N3 = my_3Dvector(c, b);
	A = N1.cross(N2);
	B = N1.cross(N3);
	if (A.dot(B) > 0)
		return 0;
	return 1;*/
}

// 判断点是否和光源之间有遮挡
void findBarrier(my_3D_point_coord start, my_3D_point_coord end, my_3D_point_coord proj, 
	int model_index, int face_index, my_draw_color formal_color, int lightcount) {
	int nobarrier = 1;
	for (int i = 0; i < all_models.size(); i++) {
		if (skipmodel[i] == 0) continue;
		for (int j = 0; j < all_models[i].faceSets.size(); j++) {
			if (i == model_index && j == face_index) continue;//若遍历到start点所在模型，跳过

			int firstPointIndex = all_models[i].faceSets[j].first_point_index;
			int secondPointIndex = all_models[i].faceSets[j].second_point_index;
			int thirdPointIndex = all_models[i].faceSets[j].third_point_index;
			my_3D_point_coord p4 = all_models[i].pointSets[firstPointIndex];
			my_3D_point_coord p5 = all_models[i].pointSets[secondPointIndex];
			my_3D_point_coord p6 = all_models[i].pointSets[thirdPointIndex];

			//求 灯光与start点连线之间 与该平面的交点坐标
			my_3Dvector N;
			my_3D_point_coord res;
			N.dx = all_models[i].faceSets[j].n.dx;
			N.dy = all_models[i].faceSets[j].n.dy;
			N.dz = all_models[i].faceSets[j].n.dz;
			my_3Dvector A(start, p4);
			my_3Dvector B(start, end);
			long double w = fabsf(N.dot(A)) / fabsf(N.dot(B));
			res.x = B.dx * w + start.x;
			res.y = B.dy * w + start.y;
			res.z = B.dz * w + start.z;

			my_3Dvector a1(start, res); // 判断start和end在面两侧
			my_3Dvector a2(end, res);

			//第一个灯
			if (lightcount == 1) {
				//若 交点在三角面内 并 投影点没有施加过阴影
				if (inface3D(res, p4, p5, p6) && a1.dot(a2) < 0
					&& shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] == 0) {
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].r -= 0.4;
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].g -= 0.4;
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].b -= 0.4;
					shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = 1;
				}
			}
			//第二个灯 判断第二个灯和点连线之间有没有遮挡物
			else if (lightcount == 2) {
				if (inface3D(res, p4, p5, p6) && a1.dot(a2) < 0) {
					nobarrier = 0;
					break;
				}
			}
		}
		if (!nobarrier && lightcount == 2) break; //若到第二个灯有遮挡物 跳出循环
	}
	//若到第二个灯没有遮挡物 该点颜色置为无阴影颜色
	if (nobarrier && lightcount == 2) {
		if (shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] == 1) {
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].r = formal_color.r;
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].g = formal_color.g;
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].b = formal_color.b;
			shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = 0;
		}
	}
}

// 计算投影点的真实坐标
void finddepth_BPRM_shadow(my_3D_point_coord eye, my_3D_point_coord proj, int model_index, int face_index) {
	//求真实坐标res
	my_3Dvector N;
	my_3D_point_coord a, res;
	a = all_models[model_index].pointSets[all_models[model_index].faceSets[face_index].third_point_index];
	N.dx = all_models[model_index].faceSets[face_index].n.dx;
	N.dy = all_models[model_index].faceSets[face_index].n.dy;
	N.dz = all_models[model_index].faceSets[face_index].n.dz;
	my_3Dvector A(eye, a);
	my_3Dvector B(eye, proj);
	long double w = fabsf(N.dot(A)) / fabsf(N.dot(B));
	res.x = B.dx * w + eye.x;
	res.y = B.dy * w + eye.y;
	res.z = B.dz * w + eye.z;
	double pointdepth = sqrtf(res.x * res.x + res.y * res.y + (res.z - eye.z) * (res.z - eye.z));

	//若 res比zbuffer中更靠近视点
	if (Zbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] > pointdepth) {
		Zbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = (float)pointdepth;

		N.normalized();
		my_draw_color pcolor, colortemp;
		pcolor.r = pcolor.g = pcolor.b = 0;

		//从面实例中取出三个顶点
		int firstPointIndex = all_models[model_index].faceSets[face_index].first_point_index;//取出顶点索引
		int secondPointIndex = all_models[model_index].faceSets[face_index].second_point_index;
		int thirdPointIndex = all_models[model_index].faceSets[face_index].third_point_index;

		my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex];//第一个顶点
		my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //第二个顶点
		my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //第三个顶点

		////glClearColor(1.f, 1.f, 1.f, 0.f);
		////glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		////设置纹理映射方式（混合模式）
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		////实施纹理坐标和几何坐标映射
		//if (model_index == 0)
		//	glBindTexture(GL_TEXTURE_2D, texGround0);
		//else if (model_index == 1)
		//	glBindTexture(GL_TEXTURE_2D, texGround1);
		//else if (model_index == 3)
		//	glBindTexture(GL_TEXTURE_2D, texGround3);
		//else if (model_index == 5)
		//	glBindTexture(GL_TEXTURE_2D, texGround5);
		////GL_BUFFER
		////glPointSize(0.001);
		////glBegin(GL_POINTS);//开始绘制
		//glColor3f(1, 1, 1);
		////glColor3f(p1color.r, p1color.g, p1color.b);
		//glTexCoord2f(0.0f, 0.0f);
		//glVertex3f(p1.x, p1.y, p1.z);
		//glColor3f(1, 1, 1);
		////glColor3f(p2color.r, p2color.g, p2color.b);
		//glTexCoord2f(0.0f, 1.0f);
		//glVertex3f(p2.x, p2.y, p2.z);
		//glColor3f(1, 1, 1);
		////glColor3f(p3color.r, p3color.g, p3color.b);
		//glTexCoord2f(1.0f, 0.0f);
		//glVertex3f(p3.x, p3.y, p3.z);
		////glEnd();

		if (open_light) {
			colortemp.r = colortemp.g = colortemp.b = 0;
			colortemp = calculate_direct_light_on_one_vertex_usingBPRM(res,
				N, eye_position, light_position, light_rgb_ambient,
				all_models[model_index].material_ambient_rgb_reflection,
				light_rgb_diffuse_specular,
				all_models[model_index].material_diffuse_rgb_reflection,
				all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
			pcolor.r += colortemp.r; pcolor.g += colortemp.g; pcolor.b += colortemp.b;
		}
		if (open_light1) {
			colortemp.r = colortemp.g = colortemp.b = 0;
			colortemp = calculate_direct_light_on_one_vertex_usingBPRM(res,
				N, eye_position, light1_position, light1_rgb_ambient,
				all_models[model_index].material_ambient_rgb_reflection,
				light_rgb_diffuse_specular,
				all_models[model_index].material_diffuse_rgb_reflection,
				all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
			pcolor.r += colortemp.r; pcolor.g += colortemp.g; pcolor.b += colortemp.b;
		}

		Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].r = pcolor.r;
		Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].g = pcolor.g;
		Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].b = pcolor.b;

		if (showshadow) {
			//查找该点和光源之间是否有遮挡物 若有 该点亮度降低
			if(open_light && !open_light1)
				findBarrier(res, light_position, proj, model_index, face_index, pcolor, 1);
			else if(open_light1 && !open_light)
				findBarrier(res, light1_position, proj, model_index, face_index, pcolor, 1);
			else if (open_light && open_light1) {
				findBarrier(res, light_position, proj, model_index, face_index, pcolor, 1);
				findBarrier(res, light1_position, proj, model_index, face_index, pcolor, 2);
			}
		}	
	}
}

// 自写zbuffer
void My_zbuffer() {
	//初始化zbuffer和framebuffer和shadowbuffer
	for (int i = 0; i <= nearplane_width; i++) {
		for (int j = 0; j <= nearplane_height; j++) {
			Zbuffer[i][j] = FAREST;
			Framebuffer[i][j].r = Framebuffer[i][j].g = Framebuffer[i][j].b = 0;
			shadowbuffer[i][j] = 0;
		}
	}

	//填充zbuffer和framebuffer
	for (int model_index = 0; model_index < all_models.size(); model_index++) {
		if (skipmodel[model_index] == 0) continue;
		for (int face_index = 0; face_index < all_models[model_index].faceSets.size(); face_index++) {
			//先背面剔除
			eyesight = my_3Dvector(eye_position,
				all_models[model_index].pointSets[all_models[model_index].faceSets[face_index].first_point_index]);
			if (all_models[model_index].faceSets[face_index].n.dot(eyesight) > 0) continue;

			//从面实例中取出三个顶点
			int firstPointIndex = all_models[model_index].faceSets[face_index].first_point_index; //取出顶点索引
			int secondPointIndex = all_models[model_index].faceSets[face_index].second_point_index;
			int thirdPointIndex = all_models[model_index].faceSets[face_index].third_point_index;

			int firstNormalIndex = all_models[model_index].faceSets[face_index].first_point_normal_index; //取出法向量索引
			int secondNormalIndex = all_models[model_index].faceSets[face_index].second_point_normal_index;
			int thirdNormalIndex = all_models[model_index].faceSets[face_index].third_point_normal_index;

			my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex]; //第一个顶点
			my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //第二个顶点
			my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //第三个顶点

			my_3Dvector p1Normal = all_models[model_index].pointNormalSets[firstNormalIndex];//第一个顶点法向
			my_3Dvector p2Normal = all_models[model_index].pointNormalSets[secondNormalIndex];//第二个顶点法向
			my_3Dvector p3Normal = all_models[model_index].pointNormalSets[thirdNormalIndex];//第三个顶点法向

			my_draw_color p1color, p2color, p3color, colortemp;
			p1color.r = p1color.g = p1color.b = 0;
			p2color.r = p2color.g = p2color.b = 0;
			p3color.r = p3color.g = p3color.b = 0;
			if (open_light) {
				//用Blinn-Phong Reflection Model计算每一点到光源位置的能量
				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p1,
					p1Normal, eye_position, light_position, light_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p1color.r += colortemp.r; p1color.g += colortemp.g; p1color.b += colortemp.b;

				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p2,
					p2Normal, eye_position, light_position, light_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p2color.r += colortemp.r; p2color.g += colortemp.g; p2color.b += colortemp.b;

				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p3,
					p3Normal, eye_position, light_position, light_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p3color.r += colortemp.r; p3color.g += colortemp.g; p3color.b += colortemp.b;
			}
			if (open_light1) {
				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p1,
					p1Normal, eye_position, light1_position, light1_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light1_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p1color.r += colortemp.r; p1color.g += colortemp.g; p1color.b += colortemp.b;

				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p2,
					p2Normal, eye_position, light1_position, light1_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light1_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p2color.r += colortemp.r; p2color.g += colortemp.g; p2color.b += colortemp.b;

				colortemp.r = colortemp.g = colortemp.b = 0;
				colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p3,
					p3Normal, eye_position, light1_position, light1_rgb_ambient,
					all_models[model_index].material_ambient_rgb_reflection,
					light1_rgb_diffuse_specular,
					all_models[model_index].material_diffuse_rgb_reflection,
					all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
				p3color.r += colortemp.r; p3color.g += colortemp.g; p3color.b += colortemp.b;
			}

			my_3D_point_coord p1proj, p2proj, p3proj;

			p1proj.x = (floor)(300.0 / (500.0 - p1.z) * p1.x + 0.5);
			p1proj.y = (floor)(300.0 / (500.0 - p1.z) * p1.y + 0.5);
			p1proj.z = 200;

			p2proj.x = (floor)(300.0 / (500.0 - p2.z) * p2.x + 0.5);
			p2proj.y = (floor)(300.0 / (500.0 - p2.z) * p2.y + 0.5);
			p2proj.z = 200;

			p3proj.x = (floor)(300.0 / (500.0 - p3.z) * p3.x + 0.5);
			p3proj.y = (floor)(300.0 / (500.0 - p3.z) * p3.y + 0.5);
			p3proj.z = 200;

			if (inprojection(p1proj, nearplane_width, nearplane_height) //若 面的顶点在投影面内 && zbuffer中在对应位置小于新值
				&& Zbuffer[(int)(p1proj.x) + nearplane_width / 2 + 1][(int)(p1proj.y) + nearplane_height / 2 + 1]
				> my_3Dvector(eye_position, p1).len) {
				Zbuffer[(int)(p1proj.x) + nearplane_width / 2 + 1]
					[(int)(p1proj.y) + nearplane_height / 2 + 1]
				= my_3Dvector(eye_position, p1).len;
				Framebuffer[(int)(p1proj.x) + nearplane_width / 2 + 1][(int)(p1proj.y) + nearplane_height / 2 + 1].r = p1color.r;
				Framebuffer[(int)(p1proj.x) + nearplane_width / 2 + 1][(int)(p1proj.y) + nearplane_height / 2 + 1].g = p1color.g;
				Framebuffer[(int)(p1proj.x) + nearplane_width / 2 + 1][(int)(p1proj.y) + nearplane_height / 2 + 1].b = p1color.b;
			}

			if (inprojection(p2proj, nearplane_width, nearplane_height)
				&& Zbuffer[(int)(p2proj.x) + nearplane_width / 2 + 1][(int)(p2proj.y) + nearplane_height / 2 + 1]
				> my_3Dvector(eye_position, p2).len) {
				Zbuffer[(int)(p2proj.x) + nearplane_width / 2 + 1]
					[(int)(p2proj.y) + nearplane_height / 2 + 1]
				= my_3Dvector(eye_position, p2).len;
				Framebuffer[(int)(p2proj.x) + nearplane_width / 2 + 1][(int)(p2proj.y) + nearplane_height / 2 + 1].r = p2color.r;
				Framebuffer[(int)(p2proj.x) + nearplane_width / 2 + 1][(int)(p2proj.y) + nearplane_height / 2 + 1].g = p2color.g;
				Framebuffer[(int)(p2proj.x) + nearplane_width / 2 + 1][(int)(p2proj.y) + nearplane_height / 2 + 1].b = p2color.b;
			}

			if (inprojection(p3proj, nearplane_width, nearplane_height)
				&& Zbuffer[(int)(p3proj.x) + nearplane_width / 2 + 1][(int)(p3proj.y) + nearplane_height / 2 + 1]
				> my_3Dvector(eye_position, p3).len) {
				Zbuffer[(int)(p3proj.x) + nearplane_width / 2 + 1]
					[(int)(p3proj.y) + nearplane_height / 2 + 1]
				= my_3Dvector(eye_position, p3).len;
				Framebuffer[(int)(p3proj.x) + nearplane_width / 2 + 1][(int)(p3proj.y) + nearplane_height / 2 + 1].r = p3color.r;
				Framebuffer[(int)(p3proj.x) + nearplane_width / 2 + 1][(int)(p3proj.y) + nearplane_height / 2 + 1].g = p3color.g;
				Framebuffer[(int)(p3proj.x) + nearplane_width / 2 + 1][(int)(p3proj.y) + nearplane_height / 2 + 1].b = p3color.b;
			}

			for (int i = 1; i <= nearplane_width; i++) //遍历投影面每一个像素
				for (int j = 1; j <= nearplane_height; j++)
					//若 遍历像素在 现在处理中的面的范围内
					if (inface(my_3D_point_coord(i - nearplane_width / 2 - 1, j - nearplane_height / 2 - 1, 200), p1proj, p2proj, p3proj))
						//求 该点深度
						finddepth_BPRM_shadow(eye_position, my_3D_point_coord(i - nearplane_width / 2 - 1, j - nearplane_height / 2 - 1, 200),
							model_index, face_index);
			printf("第 %d / %d 个模型，第 %d / %d 个面zbuffer完成\t(面总数未考虑背面剔除)\n",
				model_index + 1, all_models.size(), face_index + 1, all_models[model_index].faceSets.size());
		}
	}
	printManual();
}




// 打印说明书
void printManual() {
	cout << "\n\t\t\t\b\bManual" << endl;
	cout << "\t/--------------------------------|" << endl;
	cout << "\t|[1] 操作摄像机\t\t\t |" << endl;
	cout << "\t|[2] 操作灯1\t\t\t |" << endl;
	cout << "\t|[3] 操作灯2\t\t\t |" << endl;
	cout << "\t|[4] 修复光追蜡笔效果\t\t |" << endl;
	cout << "\t|[5] 视线跟随物体\t\t |" << endl;
	cout << "\t|[6] 显示旋转轴\t\t\t |" << endl;
	cout << "\t|[`] 显示新说明书\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[w][a][s][d] 控制平移\t\t |" << endl;
	cout << "\t|[q][e] 控制旋转\t\t |" << endl;
	cout << "\t|[i][k] 控制升降\t\t |" << endl;
	cout << "\t|[z] 开关自带zbuffer\t\t |" << endl;
	cout << "\t|[x] 开关自写zbuffer\t\t |" << endl;
	cout << "\t|[o] 开关灯\t\t\t |" << endl;
	cout << "\t|[p] 物体复位\t\t\t |" << endl;
	cout << "\t|[r] 开关光追\t\t\t |" << endl;
	cout << "\t|[l] 开关阴影显示\t\t |" << endl;
	cout << "\t|[y] 开关纹理\t\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[t][ ][ ]\t\t\t |" << endl;
	cout << "\t|[g][h][j] 开关模型显示\t\t |" << endl;
	cout << "\t|[b][n][m]\t\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[,][.][/] 切换调整rgb\t\t |" << endl;
	cout << "\t|[\\] 切换调整灯亮度\t\t |" << endl;
	cout << "\t|[[][]] 调整rgb\t\t\t |" << endl;
	cout << "\t|--------------------------------/\n" << endl;
}

// 打印光追进度
int printprogress(long double i) {
	printf("已完成 [ ");
	if (100 * i / (image_w * image_h) == 100.0) {
		printf("* * * * * * * * * * ] 100.00 %%");
		return 1;
	}
	for (int u = 0; u < (int)(i * 10 / (image_w * image_h)) % 10; u++) {
		printf("* ");
	}
	for (int v = ((int)(i * 10 / (image_w * image_h)) % 10); v < 10; v++) {
		printf("- ");
	}
	printf("] %.2f %%\n", 100 * i / (image_w * image_h));
	return 1;
}

// 模型平移
void translate(float delta, int dir, int light) {
	/*平移摄像机*/
	if (light == 0) {
		//平移模型
		for (int model_index = 0; model_index < all_models.size(); model_index++) {
			for (int i = 0; i < all_models[model_index].pointSets.size(); i++) {
				if (dir == 1)/*123456前后左右上下*/ {
					all_models[model_index].pointSets[i].z -= delta;
				}
				else if (dir == 2) {
					all_models[model_index].pointSets[i].z += delta;
				}
				else if (dir == 3) {
					all_models[model_index].pointSets[i].x -= delta;
				}
				else if (dir == 4) {
					all_models[model_index].pointSets[i].x += delta;
				}
				else if (dir == 5) {
					all_models[model_index].pointSets[i].y += delta;
				}
				else if (dir == 6) {
					all_models[model_index].pointSets[i].y -= delta;
				}
			}
		}
		//平移灯和原点
		if (dir == 1) /*123456前后左右上下*/ {
			light_position.z -= delta;
			light1_position.z -= delta;
			original.z -= delta;
		}
		else if (dir == 2) {
			light_position.z += delta;
			light1_position.z += delta;
			original.z += delta;
		}
		else if (dir == 3) {
			light_position.x -= delta;
			light1_position.x -= delta;
			original.x -= delta;
		}
		else if (dir == 4) {
			light_position.x += delta;
			light1_position.x += delta;
			original.x += delta;
		}
		else if (dir == 5) {
			light_position.y += delta;
			light1_position.y += delta;
			original.y += delta;
		}
		else if (dir == 6) {
			light_position.y -= delta;
			light1_position.y -= delta;
			original.y -= delta;
		}
	}
	/* 平移灯1  */
	else if (light == 1) {
		if (dir == 1) /*123456前后左右上下*/
			light_position.z -= delta;
		else if (dir == 2)
			light_position.z += delta;
		else if (dir == 3)
			light_position.x -= delta;
		else if (dir == 4)
			light_position.x += delta;
		else if (dir == 5)
			light_position.y += delta;
		else if (dir == 6)
			light_position.y -= delta;
	}
	/* 平移灯2  */
	else if (light == 2) {
		if (dir == 1) /*123456前后左右上下*/
			light1_position.z -= delta;
		else if (dir == 2)
			light1_position.z += delta;
		else if (dir == 3)
			light1_position.x -= delta;
		else if (dir == 4)
			light1_position.x += delta;
		else if (dir == 5)
			light1_position.y += delta;
		else if (dir == 6)
			light1_position.y -= delta;
	}
}

// 矩阵相乘
my_3D_point_coord matrix_multiply_vector(float matrix[][4], my_3D_point_coord input_v) {
	my_3D_point_coord translated_v;
	translated_v.x = matrix[0][0] * input_v.x + matrix[0][1] * input_v.y + matrix[0][2] * input_v.z + matrix[0][3] * 1;
	translated_v.y = matrix[1][0] * input_v.x + matrix[1][1] * input_v.y + matrix[1][2] * input_v.z + matrix[1][3] * 1;
	translated_v.z = matrix[2][0] * input_v.x + matrix[2][1] * input_v.y + matrix[2][2] * input_v.z + matrix[2][3] * 1;
	return translated_v;
}

// 模型旋转
void rotate(float angle, int target) {
	float rotate_matrix[4][4];
	memset(rotate_matrix, 0, sizeof(float) * 16);
	rotate_matrix[1][1] = rotate_matrix[3][3] = 1;
	rotate_matrix[0][0] = rotate_matrix[2][2] = cosf(angle / 180 * PI);
	rotate_matrix[0][2] = -sinf(angle / 180 * PI);
	rotate_matrix[2][0] = sinf(angle / 180 * PI);

	my_3D_point_coord input_v;
	if (target == 0)
		for (int model_index = 0; model_index < all_models.size(); model_index++) {
			for (int i = 0; i < all_models[model_index].pointSets.size(); i++) {
				input_v.x = all_models[model_index].pointSets[i].x - original.x;
				input_v.y = all_models[model_index].pointSets[i].y - original.y;
				input_v.z = all_models[model_index].pointSets[i].z - original.z;
				input_v = matrix_multiply_vector(rotate_matrix, input_v);
				all_models[model_index].pointSets[i].x = input_v.x + original.x;
				all_models[model_index].pointSets[i].y = input_v.y + original.y;
				all_models[model_index].pointSets[i].z = input_v.z + original.z;
			}
		}
	if (target == 1 || target == 0) {
		input_v.x = light_position.x - original.x;
		input_v.y = light_position.y - original.y;
		input_v.z = light_position.z - original.z;
		input_v = matrix_multiply_vector(rotate_matrix, input_v);
		light_position.x = input_v.x + original.x;
		light_position.y = input_v.y + original.y;
		light_position.z = input_v.z + original.z;
	}

	if (target == 2 || target == 0) {
		input_v.x = light1_position.x - original.x;
		input_v.y = light1_position.y - original.y;
		input_v.z = light1_position.z - original.z;
		input_v = matrix_multiply_vector(rotate_matrix, input_v);
		light1_position.x = input_v.x + original.x;
		light1_position.y = input_v.y + original.y;
		light1_position.z = input_v.z + original.z;
	}
}

// 初始化加载模型
void init(void) {
	//物体对光线的反射率
	float material_ambient_rgb_reflection[] = { 0.2, 0.2, 0.2 };
	float material_specular_rgb_reflection[] = { 0.2, 0.2, 0.2 };
	float ns = 40; //聚光指数
	//加载模型
	ObjLoader objLoader1 = ObjLoader("model\\book1.obj");
	float obj1_material_diffuse_rgb_reflection[] = { 1, 0.87, 0.87 };
	float obj1_specular_rgb_reflection[] = { 0.2,0.2,0.2 };
	objLoader1.my_3DModel.modify_color_configuration(0, 0,
		material_ambient_rgb_reflection, obj1_material_diffuse_rgb_reflection,
		obj1_specular_rgb_reflection, ns);
	all_models.push_back(objLoader1.my_3DModel);

	ObjLoader objLoader2 = ObjLoader("model\\bookface.obj");
	float obj2_material_diffuse_rgb_reflection[] = { 1, 0.99, 0.99 };
	float obj2_specular_rgb_reflection[] = { 0.2,0.2,0.2 };
	objLoader2.my_3DModel.modify_color_configuration(128, 0,
		material_ambient_rgb_reflection, obj2_material_diffuse_rgb_reflection,
		obj2_specular_rgb_reflection, ns);
	all_models.push_back(objLoader2.my_3DModel);

	ObjLoader objLoader3 = ObjLoader("model\\house.obj");
	float obj3_material_diffuse_rgb_reflection[] = { 0.99, 0.99, 0.996 };
	float obj3_specular_rgb_reflection[] = { 0.2,0.2,0.2 };
	objLoader3.my_3DModel.modify_color_configuration(0, 0,
		material_ambient_rgb_reflection, obj3_material_diffuse_rgb_reflection,
		obj3_specular_rgb_reflection, ns);
	all_models.push_back(objLoader3.my_3DModel);

	ObjLoader objLoader4 = ObjLoader("model\\lamp.obj");
	float obj4_material_diffuse_rgb_reflection[] = { 0.2, 0.2, 0.5 };
	float obj4_specular_rgb_reflection[] = { 0.2,0.2,0.2 };
	objLoader4.my_3DModel.modify_color_configuration(0, 0,
		material_ambient_rgb_reflection, obj4_material_diffuse_rgb_reflection,
		obj4_specular_rgb_reflection, ns);
	all_models.push_back(objLoader4.my_3DModel);

	ObjLoader objLoader5 = ObjLoader("model\\pen.obj");
	float obj5_material_diffuse_rgb_reflection[] = { 0.6, 0.6, 0 };
	objLoader5.my_3DModel.modify_color_configuration(9, 0,
		material_ambient_rgb_reflection, obj5_material_diffuse_rgb_reflection,
		material_specular_rgb_reflection, ns);
	all_models.push_back(objLoader5.my_3DModel);

	ObjLoader objLoader6 = ObjLoader("model\\table1.obj");
	float obj6_material_diffuse_rgb_reflection[] = { 0.9, 0.7, 0.3 };
	objLoader6.my_3DModel.modify_color_configuration(0, 0.00000001,
		material_ambient_rgb_reflection, obj6_material_diffuse_rgb_reflection,
		material_specular_rgb_reflection, ns);
	all_models.push_back(objLoader6.my_3DModel);

	ObjLoader objLoader7 = ObjLoader("model\\tableleg.obj");
	float obj7_material_diffuse_rgb_reflection[] = { 0.8, 0.8, 0.8 };
	float obj7_specular_rgb_reflection[] = { 0.2,0.2,0.2 };
	objLoader7.my_3DModel.modify_color_configuration(9, 0,
		material_ambient_rgb_reflection, obj7_material_diffuse_rgb_reflection,
		obj7_specular_rgb_reflection, ns);
	all_models.push_back(objLoader7.my_3DModel);
}

// 绘制内容
void display(void) {
	//求各面法向
	my_3Dvector N1, N2, N3;
	for (int j = 0; j < all_models.size(); j++)
		for (int i = 0; i < all_models[j].faceSets.size(); i++) {
			N1 = my_3Dvector(all_models[j].pointSets[all_models[j].faceSets[i].first_point_index],
				all_models[j].pointSets[all_models[j].faceSets[i].second_point_index]);
			N2 = my_3Dvector(all_models[j].pointSets[all_models[j].faceSets[i].first_point_index],
				all_models[j].pointSets[all_models[j].faceSets[i].third_point_index]);
			N3 = N1.cross(N2);
			all_models[j].faceSets[i].n.dx = N3.dx;
			all_models[j].faceSets[i].n.dy = N3.dy;
			all_models[j].faceSets[i].n.dz = N3.dz;
		}

	glClearColor(1.f, 1.f, 1.f, 0.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// 自写zbuffer + shadow绘制
	if (my_zbuffer) {
		My_zbuffer();
		//glClearColor(1.f, 1.f, 1.f, 0.f);
		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		for (int i = 1; i <= nearplane_width; i++) {
			for (int j = 1; j <= nearplane_height; j++) {
				glPointSize(1);
				glBegin(GL_POINTS);
				glColor3f(Framebuffer[i][j].r, Framebuffer[i][j].g, Framebuffer[i][j].b);
				glVertex3f(i - nearplane_width / 2 - 1, j - nearplane_height / 2 - 1, 200);
				glEnd();
			}
		}
	}

	// 直接光绘制
	else if (!rendered && !my_zbuffer) {
		glShadeModel(GL_SMOOTH);
		for (unsigned int model_index = 0; model_index < all_models.size(); model_index++)
		{
			if (skipmodel[model_index] == 0) continue;
			for (unsigned int i = 0; i < all_models[model_index].faceSets.size(); i++)
			{
				//先背面剔除
				eyesight = my_3Dvector(eye_position,
					all_models[model_index].pointSets[all_models[model_index].faceSets[i].first_point_index]);
				if (all_models[model_index].faceSets[i].n.dot(eyesight) > 0)
					continue;

				//从面实例中取出三个顶点
				int firstPointIndex = all_models[model_index].faceSets[i].first_point_index;//取出顶点索引
				int secondPointIndex = all_models[model_index].faceSets[i].second_point_index;
				int thirdPointIndex = all_models[model_index].faceSets[i].third_point_index;

				int firstNormalIndex = all_models[model_index].faceSets[i].first_point_normal_index; //取出法向量索引
				int secondNormalIndex = all_models[model_index].faceSets[i].second_point_normal_index;
				int thirdNormalIndex = all_models[model_index].faceSets[i].third_point_normal_index;

				my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex];//第一个顶点
				my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //第二个顶点
				my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //第三个顶点

				my_3Dvector p1Normal = all_models[model_index].pointNormalSets[firstNormalIndex];//第一个顶点法向
				my_3Dvector p2Normal = all_models[model_index].pointNormalSets[secondNormalIndex];//第二个顶点法向
				my_3Dvector p3Normal = all_models[model_index].pointNormalSets[thirdNormalIndex];//第三个顶点法向

				my_draw_color p1color, p2color, p3color, colortemp;
				p1color.r = p1color.g = p1color.b = 0;
				p2color.r = p2color.g = p2color.b = 0;
				p3color.r = p3color.g = p3color.b = 0;
				if (open_light)
				{
					glPointSize(3);
					glBegin(GL_POINTS);
					glColor3f(0, 0, 1);
					glVertex3f(light_position.x, light_position.y, light_position.z);
					glEnd();

					//用Blinn-Phong Reflection Model计算每一点到光源位置的能量
					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p1,
						p1Normal, eye_position, light_position, light_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p1color.r += colortemp.r; p1color.g += colortemp.g; p1color.b += colortemp.b;

					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p2,
						p2Normal, eye_position, light_position, light_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p2color.r += colortemp.r; p2color.g += colortemp.g; p2color.b += colortemp.b;

					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p3,
						p3Normal, eye_position, light_position, light_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p3color.r += colortemp.r; p3color.g += colortemp.g; p3color.b += colortemp.b;
				}
				if (open_light1) {
					glPointSize(3);
					glBegin(GL_POINTS);
					glColor3f(1, 0, 0);
					glVertex3f(light1_position.x, light1_position.y, light1_position.z);
					glEnd();

					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p1,
						p1Normal, eye_position, light1_position, light1_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light1_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p1color.r += colortemp.r; p1color.g += colortemp.g; p1color.b += colortemp.b;

					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p2,
						p2Normal, eye_position, light1_position, light1_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light1_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p2color.r += colortemp.r; p2color.g += colortemp.g; p2color.b += colortemp.b;

					colortemp.r = colortemp.g = colortemp.b = 0;
					colortemp = calculate_direct_light_on_one_vertex_usingBPRM(p3,
						p3Normal, eye_position, light1_position, light1_rgb_ambient,
						all_models[model_index].material_ambient_rgb_reflection,
						light1_rgb_diffuse_specular,
						all_models[model_index].material_diffuse_rgb_reflection,
						all_models[model_index].material_specular_rgb_reflection, all_models[model_index].ns);
					p3color.r += colortemp.r; p3color.g += colortemp.g; p3color.b += colortemp.b;
				}
				// 不加纹理的模型
				if (model_index != 5 && model_index != 3 && model_index != 1 && model_index != 0) {
					glBegin(GL_TRIANGLES);//开始绘制
					glColor3f(p1color.r, p1color.g, p1color.b);
					glVertex3f(p1.x, p1.y, p1.z);
					glColor3f(p2color.r, p2color.g, p2color.b);
					glVertex3f(p2.x, p2.y, p2.z);
					glColor3f(p3color.r, p3color.g, p3color.b);
					glVertex3f(p3.x, p3.y, p3.z);
					glEnd();
				}
				// 加纹理的模型
				else{
					//设置纹理映射方式（混合模式）
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
					
					//实施纹理坐标和几何坐标映射
					if (model_index == 0)
						glBindTexture(GL_TEXTURE_2D, texGround0);
					else if (model_index == 1)
						glBindTexture(GL_TEXTURE_2D, texGround1);
					else if (model_index == 3)
						glBindTexture(GL_TEXTURE_2D, texGround3);
					else if (model_index == 5)
						glBindTexture(GL_TEXTURE_2D, texGround5);

					glBegin(GL_TRIANGLES);//开始绘制
					glColor3f(p1color.r, p1color.g, p1color.b);
					glTexCoord2f(0.0f, 0.0f);
					glVertex3f(p1.x, p1.y, p1.z);
					glColor3f(p2color.r, p2color.g, p2color.b);
					glTexCoord2f(0.0f, 1.0f);
					glVertex3f(p2.x, p2.y, p2.z);
					glColor3f(p3color.r, p3color.g, p3color.b);
					glTexCoord2f(1.0f, 0.0f);
					glVertex3f(p3.x, p3.y, p3.z);
					glEnd();
				}
			}
		}
	}

	// 光线跟踪绘制
	else if (rendered && !my_zbuffer) {
		glBegin(GL_POINTS);//开始绘制
		for (std::map<my_3D_point_coord*, my_draw_color*>::iterator piter = render_vertices.begin(); piter != render_vertices.end(); piter++) {
			glColor3f(piter->second->r, piter->second->g, piter->second->b);
			glVertex3f(piter->first->x, piter->first->y, piter->first->z);
		}
		glEnd();
	}

	if (drawaxis) {
		//绘制中轴
		glPointSize(5);
		glBegin(GL_POINTS);
		glColor4f(1, 1, 0, 0.2);
		for (int i = -300; i <= 300; i++)
			glVertex3f(original.x, original.y + i, original.z);
		glEnd();

		//绘制original
		glPointSize(8);
		glBegin(GL_POINTS);
		glColor4f(1, 0, 1, 0.2);
		glVertex3f(original.x, original.y, original.z);
		glEnd();
	}

	glutSwapBuffers();
}

// 投影方式、modelview方式设置、对投影面采样
void reshape(int w, int h) {
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//if (w <= h) {//设置平行投影视景体
	//	glOrtho(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
	//		(GLfloat)nearplane_height / (GLfloat)nearplane_width, ratio * nearplane_height *
	//		(GLfloat)nearplane_height / (GLfloat)nearplane_width,
	//		nearplane_distance, farplane_distance); //相对于视点
	//}
	//else {//设置平行投影视景体
	//	glOrtho(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
	//		(GLfloat)nearplane_width / (GLfloat)nearplane_height, ratio * nearplane_height *
	//		(GLfloat)nearplane_width / (GLfloat)nearplane_height,
	//		nearplane_distance, farplane_distance);
	//}

	if (w <= h) {//设置透视视景体
		glFrustum(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
			(GLfloat)nearplane_height / (GLfloat)nearplane_width, ratio * nearplane_height *
			(GLfloat)nearplane_height / (GLfloat)nearplane_width,
			nearplane_distance, farplane_distance); //相对于视点
	}
	else {//设置透视视景体
		glFrustum(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
			(GLfloat)nearplane_width / (GLfloat)nearplane_height, ratio * nearplane_height *
			(GLfloat)nearplane_width / (GLfloat)nearplane_height,
			nearplane_distance, farplane_distance);
	}
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if (sighttype == 1)
		gluLookAt(eye_position.x, eye_position.y, eye_position.z, original.x, original.y, original.z, 0, 1, 0);
	else if (sighttype == 0)
		gluLookAt(eye_position.x, eye_position.y, eye_position.z, 0, 0, 0, 0, 1, 0);
}

// 键盘交互事件
void keyboard(unsigned char key, int x, int y) {
	switch (key)
	{
	case '1':			/*    摄像机  */ {
		selected = 1;
		printf("已选中 [摄像机]\n");
		break;
	}
	case '2':			/*     灯1    */ {
		selected = 2;
		printf("已选中 [灯1]\n");
		break;
	}
	case '3':			/*     灯2    */ {
		selected = 3;
		printf("已选中 [灯2]\n");
		break;
	}
	case '4':			/* 渲染后偏移 */ {
		if (renderoffset == 0) {
			renderoffset = 9.9;
			printf("光追微调 [ON]\n");
		}
		else {
			renderoffset = 0;
			printf("光追微调 [OFF]\n");
		}
		break;
	}
	case '5':			/*  视点朝向  */ {
		if (sighttype == 1) {
			sighttype = 0;
			printf("视线锁定 [OFF]\n");
		}
		else {
			sighttype = 1;
			printf("视线锁定 [ON]\n");
		}
		reshape(nearplane_width, nearplane_height);
		glutPostRedisplay();
		break;
	}
	case '6':			/*  中轴绘制  */ {
		if (drawaxis)
			drawaxis = 0;
		else
			drawaxis = 1;
		glutPostRedisplay();
		break;
	}
	case '`':			/*  显示指南  */ {
		system("cls");
		printManual();
		break;
	}
	case 'q': case 'Q': /* 摄像机旋转 */ {
		if (selected == 1) {
			rotate(-theta, 0);//逆时针旋转
			rotateoffset -= theta;
			glutPostRedisplay();
		}
		break;
	}
	case 'e': case 'E': /* 摄像机旋转 */ {
		if (selected == 1) {
			rotate(theta, 0);//顺时针旋转
			rotateoffset += theta;
			glutPostRedisplay();
		}
		break;
	}
	case 'w': case 'W': /*    前进    */ {
		if (selected == 1) {//摄像机
			translate(delta, 2, 0);
			cameraoffset.z += delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 1, 1);
			lightoffset.z -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 1, 2);
			light1offset.z -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 's': case 'S': /*    后退    */ {
		if (selected == 1) {//摄像机
			translate(delta, 1, 0);
			cameraoffset.z -= delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 2, 1);
			lightoffset.z += 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 2, 2);
			light1offset.z += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'a': case 'A': /*    向左    */ {
		if (selected == 1) {//摄像机
			translate(delta, 4, 0);
			cameraoffset.x += delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 3, 1);
			lightoffset.x -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 3, 2);
			light1offset.x -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'd': case 'D': /*    向右    */ {
		if (selected == 1) {//摄像机
			translate(delta, 3, 0);
			cameraoffset.x -= delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 4, 1);
			lightoffset.x += 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 4, 2);
			light1offset.x += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'i': case 'I': /*    向上    */ {
		if (selected == 1) {//摄像机
			translate(delta - renderoffset, 6, 0);
			cameraoffset.y -= (delta - renderoffset);
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 5, 1);
			lightoffset.y += 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 5, 2);
			light1offset.y += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'k': case 'K': /*    向下    */ {
		if (selected == 1) {//摄像机
			translate(delta - renderoffset, 5, 0);
			cameraoffset.y += (delta - renderoffset);
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//灯1
			translate(10, 6, 1);
			lightoffset.y -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//灯2
			translate(10, 6, 2);
			light1offset.y -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'p': case 'P': /*    复位    */ {
		if (selected == 1) {
			rotate(-rotateoffset, 0);
			translate(cameraoffset.x, 3, 0);
			translate(cameraoffset.y, 6, 0);
			translate(cameraoffset.z, 1, 0);
			cameraoffset.x = cameraoffset.y = cameraoffset.z = 0;
			rotateoffset = 0;
			original.x = original.y = original.z = 0;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {
			light_position.x = light_default_position.x + original.x;
			light_position.y = light_default_position.y + original.y;
			light_position.z = light_default_position.z + original.z;
			rotate(rotateoffset, 1);
			lightoffset.x = lightoffset.y = lightoffset.z = 0;
		}
		else if (selected == 3 && open_light1 == true) {
			light1_position.x = light1_default_position.x + original.x;
			light1_position.y = light1_default_position.y + original.y;
			light1_position.z = light1_default_position.z + original.z;
			rotate(rotateoffset, 2);
			light1offset.x = light1offset.y = light1offset.z = 0;
		}
		glutPostRedisplay();
		break;
	}
	case 'o': case 'O': /*   开关灯   */ {
		if (selected == 2) {
			if (open_light)
				open_light = false;
			else
				open_light = true;
		}
		else if (selected == 3) {
			if (open_light1)
				open_light1 = false;
			else
				open_light1 = true;
		}
		glutPostRedisplay();
		break;
	}
	case 'z': case 'Z': /*   zbuffer  */ {
		zbuffer += 1;
		if (zbuffer) {
			glEnable(GL_DEPTH_TEST); //打开深度缓冲测试
			glDepthFunc(GL_LESS); //判断遮挡关系时，离视点近的物体遮挡离视点远的物体	
			zbuffer -= 2;
		}
		else if (!zbuffer) {
			glDisable(GL_DEPTH_TEST); //关闭深度缓冲测试
		}
		glutPostRedisplay();
		break;
	}
	case 'x': case 'X': /* My_zbuffer */ {
		if (!my_zbuffer) {
			my_zbuffer = 1;
		}
		else if (my_zbuffer) {
			my_zbuffer = 0;
		}
		glutPostRedisplay();
		break;
	}
	case 'l': case 'L': /*  阴影显示  */ {
		if (!showshadow) {
			showshadow = 1;
			printf("阴影显示 [ON]\n");
		}
		else if (showshadow) {
			showshadow = 0;
			printf("阴影显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'r': case 'R': /*  光线追踪  */ {
		//依据当前视景体设置对投影平面进行顶点采样
		if (rendered == true) {
			rendered = false;
			render_vertices.clear();
			break;
		}
		samplepoint_sonprojectionplan(-ratio * nearplane_width, ratio * nearplane_width,
			-ratio * nearplane_height * (GLfloat)nearplane_height / (GLfloat)nearplane_width,
			ratio * nearplane_height * (GLfloat)nearplane_height / (GLfloat)nearplane_width,
			nearplane_distance, eye_position.z, render_vertices, image_w, image_h, gap);
		long double i = 0;
		for (std::map<my_3D_point_coord*, my_draw_color*>::iterator piter = render_vertices.begin(); piter != render_vertices.end(); piter++) {
			i++;
			my_3Dvector raydir(eye_position, *(piter->first));
			raydir.normalized();
			my_draw_color newColor2; newColor2.r = newColor2.g = newColor2.b = 0;
			my_draw_color newColor3; newColor3.r = newColor3.g = newColor3.b = 0;
			if (open_light)
				newColor2 = one_ray_trace_my(eye_position, raydir, all_models,
					0, eye_position, light_position, light_rgb_ambient, light_rgb_diffuse_specular);
			if (open_light1)
				newColor3 = one_ray_trace_my(eye_position, raydir, all_models,
					0, eye_position, light1_position, light1_rgb_ambient, light1_rgb_diffuse_specular);
			my_draw_color drawcolor;
			drawcolor.r = newColor2.r + newColor3.r;
			drawcolor.g = newColor2.g + newColor3.g;
			drawcolor.b = newColor2.b + newColor3.b;
			*(piter->second) = drawcolor;
			printprogress(i);

		}
		printManual();
		rendered = true;
		glutPostRedisplay();
		break;
	}
	case 'b': case 'B': /* 开关模型1  */ {
		if (skipmodel[0] != 1) {
			skipmodel[0] = 1;
			printf("模型1 [书皮] 显示 [ON]\n");
		}
		else if (skipmodel[0] == 1) {
			skipmodel[0] = 0;
			printf("模型1 [书皮] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'n': case 'N': /* 开关模型2  */ {
		if (skipmodel[1] != 1) {
			skipmodel[1] = 1;
			printf("模型2 [书页] 显示 [ON]\n");
		}
		else if (skipmodel[1] == 1) {
			skipmodel[1] = 0;
			printf("模型2 [书页] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'm': case 'M': /* 开关模型3  */ {
		if (skipmodel[2] != 1) {
			skipmodel[2] = 1;
			printf("模型3 [房间] 显示 [ON]\n");
		}
		else if (skipmodel[2] == 1) {
			skipmodel[2] = 0;
			printf("模型3 [房间] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'g': case 'G': /* 开关模型4  */ {
		if (skipmodel[3] != 1) {
			skipmodel[3] = 1;
			printf("模型4 [台灯] 显示 [ON]\n");
		}
		else if (skipmodel[3] == 1) {
			skipmodel[3] = 0;
			printf("模型4 [台灯] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'h': case 'H': /* 开关模型5  */ {
		if (skipmodel[4] != 1) {
			skipmodel[4] = 1;
			printf("模型5 [铅笔] 显示 [ON]\n");
		}
		else if (skipmodel[4] == 1) {
			skipmodel[4] = 0;
			printf("模型5 [铅笔] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'j': case 'J': /* 开关模型6  */ {
		if (skipmodel[5] != 1) {
			skipmodel[5] = 1;
			printf("模型6 [桌面] 显示 [ON]\n");
		}
		else if (skipmodel[5] == 1) {
			skipmodel[5] = 0;
			printf("模型6 [桌面] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 't': case 'T': /* 开关模型7  */ {
		if (skipmodel[6] != 1) {
			skipmodel[6] = 1;
			printf("模型7 [桌腿] 显示 [ON]\n");
		}
		else if (skipmodel[6] == 1) {
			skipmodel[6] = 0;
			printf("模型7 [桌腿] 显示 [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'y': case 'Y': /*  开关纹理  */ {
		if (!texenable) {
			texenable = 1;
			glEnable(GL_TEXTURE_2D);
			texGround0 = load_texture("tex\\book.bmp");  //加载纹理
			texGround1 = load_texture("tex\\bookface.bmp");
			texGround3 = load_texture("tex\\lamp.bmp");
			texGround5 = load_texture("tex\\table.bmp");
		}
		else if (texenable) {
			texenable = 0;
			glDisable(GL_TEXTURE_2D);
			texGround0 = NULL;
			texGround1 = NULL;
			texGround3 = NULL;
			texGround5 = NULL;
		}
		glutPostRedisplay();
		break;
	}
	case ',':		    /* 灯光选中r  */ {
		selectedrgb = 0;
		printf("光照已选中 [R]\n");
		break;
	}
	case '.':		    /* 灯光选中g  */ {
		selectedrgb = 1;
		printf("光照已选中 [G]\n");
		break;
	}
	case '/':		    /* 灯光选中b  */ {
		selectedrgb = 2;
		printf("光照已选中 [B]\n");
		break;
	}
	case '\\':		    /*灯光选中亮度*/ {
		selectedrgb = 3;
		printf("光照已选中 [亮度]\n");
		break;
	}
	case '[':		    /*上调灯rgb值 */ {
		if (selected == 2) {
			if (light_rgb_ambient[selectedrgb] < 1 && open_light && selectedrgb != 3)
				light_rgb_ambient[selectedrgb] += 0.05;
			else if (selectedrgb == 3 && light_rgb_ambient[0] < 1 && light_rgb_ambient[1] < 1 && light_rgb_ambient[2] < 1) {
				light_rgb_ambient[0] += 0.05;
				light_rgb_ambient[1] += 0.05;
				light_rgb_ambient[2] += 0.05;
			}
		}
		else if (selected == 3) {
			if (light1_rgb_ambient[selectedrgb] < 1 && open_light1 && selectedrgb != 3)
				light1_rgb_ambient[selectedrgb] += 0.05;
			else if (selectedrgb == 3 && light1_rgb_ambient[0] < 1 && light1_rgb_ambient[1] < 1 && light1_rgb_ambient[2] < 1) {
				light1_rgb_ambient[0] += 0.05;
				light1_rgb_ambient[1] += 0.05;
				light1_rgb_ambient[2] += 0.05;
			}
		}
		glutPostRedisplay();
		break;
	}
	case ']':		    /*下调灯rgb值 */ {
		if (selected == 2) {
			if (light_rgb_ambient[selectedrgb] > 0 && open_light && selectedrgb != 3)
				light_rgb_ambient[selectedrgb] -= 0.05;
			else if (selectedrgb == 3 && light_rgb_ambient[0] > 0 && light_rgb_ambient[1] > 0 && light_rgb_ambient[2] > 0) {
				light_rgb_ambient[0] -= 0.05;
				light_rgb_ambient[1] -= 0.05;
				light_rgb_ambient[2] -= 0.05;
			}
		}
		else if (selected == 3) {
			if (light1_rgb_ambient[selectedrgb] > 0 && open_light1 && selectedrgb != 3)
				light1_rgb_ambient[selectedrgb] -= 0.05;
			else if (selectedrgb == 3 && light1_rgb_ambient[0] > 0 && light1_rgb_ambient[1] > 0 && light1_rgb_ambient[2] > 0) {
				light1_rgb_ambient[0] -= 0.05;
				light1_rgb_ambient[1] -= 0.05;
				light1_rgb_ambient[2] -= 0.05;
			}
		}
		glutPostRedisplay();
		break;
	}
	case 27:			/*  退出程序  */ {
		exit(0);
		break;
	}
	}
}

// 鼠标交互事件
void mouse(int button, int state, int x, int y) {
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		if (state == GLUT_DOWN) {

		}
		break;
	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN) {

		}
		break;
	default:
		break;
	}
}

// 主调函数
int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(nearplane_width, nearplane_height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("四男生组 - 期末大作业");
	init();
	printManual();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();
	return 0;
}