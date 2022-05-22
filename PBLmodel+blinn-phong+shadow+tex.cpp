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
#define BMP_Header_Length 54									//ͼ���������ڴ���е�ƫ����

// �ж������ǲ���2����������
int power_of_two(int n) {
	if (n <= 0)
		return 0;
	return (n & (n - 1)) == 0;
}

// ��ȡһ��BMP�ļ���Ϊ�������ʧ�ܣ�����0������ɹ�������������
GLuint load_texture(const char* file_name) {
	GLint width, height, total_bytes;
	GLubyte* pixels = 0;
	GLuint last_texture_ID = 0, texture_ID = 0;

	// ���ʧ��
	FILE* pFile = fopen(file_name, "rb");
	if (pFile == 0)
		return 0;

	// ��ȡͼ�Ŀ�Ⱥ͸߶�
	fseek(pFile, 0x0012, SEEK_SET);
	fread(&width, 4, 1, pFile);
	fread(&height, 4, 1, pFile);
	fseek(pFile, BMP_Header_Length, SEEK_SET);

	// ����ÿ��������ռ�ֽ����������ݴ����ݼ����������ֽ���
	{
		GLint line_bytes = width * 3;
		while (line_bytes % 4 != 0)
			++line_bytes;
		total_bytes = line_bytes * height;
	}

	// �����������ֽ��������ڴ�
	pixels = (GLubyte*)malloc(total_bytes);
	if (pixels == 0) {
		fclose(pFile);
		return 0;
	}

	// ��ȡ��������
	if (fread(pixels, total_bytes, 1, pFile) <= 0) {
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// ���ͼ�Ŀ�Ⱥ͸߶Ȳ���2�������η�������
	// ��ͼ���߳�����OpenGL�涨�����ֵ������
	{
		GLint max;
		glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max);
		if (!power_of_two(width)
			|| !power_of_two(height)
			|| width > max
			|| height > max) {
			const GLint new_width = 256;
			const GLint new_height = 256; // �涨���ź��µĴ�СΪ�߳���������
			GLint new_line_bytes, new_total_bytes;
			GLubyte* new_pixels = 0;

			// ����ÿ����Ҫ���ֽ��������ֽ���
			new_line_bytes = new_width * 3;
			while (new_line_bytes % 4 != 0)
				++new_line_bytes;
			new_total_bytes = new_line_bytes * new_height;

			// �����ڴ�
			new_pixels = (GLubyte*)malloc(new_total_bytes);
			if (new_pixels == 0) {
				free(pixels);
				fclose(pFile);
				return 0;
			}

			// ������������
			gluScaleImage(GL_RGB, width, height, GL_UNSIGNED_BYTE, pixels, new_width, new_height, GL_UNSIGNED_BYTE, new_pixels);

			// �ͷ�ԭ�����������ݣ���pixelsָ���µ��������ݣ�����������width��height
			free(pixels);
			pixels = new_pixels;
			width = new_width;
			height = new_height;
		}
	}

	// ����һ���µ�������
	glGenTextures(1, &texture_ID);
	if (texture_ID == 0) {
		free(pixels);
		fclose(pFile);
		return 0;
	}

	// ���µ������������������������
	// �ڰ�ǰ���Ȼ��ԭ���󶨵������ţ��Ա��������лָ�
	GLint lastTextureID = last_texture_ID;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &lastTextureID);
	glBindTexture(GL_TEXTURE_2D, texture_ID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);
	glBindTexture(GL_TEXTURE_2D, lastTextureID);  //�ָ�֮ǰ�������
	free(pixels);
	return texture_ID;
}



#define PI 3.14159265
#define nearplane_width 400										//�Ӿ�����
#define nearplane_height 400									//�Ӿ���߶�
int nearplane_distance = 300;									//�Ӿ����ƽ�����ӵ����
int farplane_distance = nearplane_distance + 1800;				//�Ӿ���Զƽ�����ӵ����

float theta = 5;												//�ӵ���ת�Ƕ�(�Ƕ���)
float delta = 10;												//λ����
float ratio = 0.5;												//ͶӰ���С����
float diffuse = 0.8;											//������ϵ��

my_3Dvector eyesight;											//�����������
float rotateoffset = 0;											//��תƫ��
my_3D_point_coord cameraoffset(0, 0, 0);						//�����ƽ��ƫ��
my_3D_point_coord lightoffset(0, 0, 0);							//��1ƽ��ƫ��
my_3D_point_coord light1offset(0, 0, 0);						//��2ƽ��ƫ��
my_3D_point_coord original(0, 0, 0);							//ԭ��λ��
my_3D_point_coord eye_default_position(0, 0, 500);				//�ӵ�Ĭ��λ��
my_3D_point_coord light_default_position(100.0, 50.0, 0);		//��1Ĭ��λ��
my_3D_point_coord light1_default_position(-100.0, 0, 10.0);		//��2Ĭ��λ��
my_3D_point_coord eye_position = { eye_default_position.x,
								  eye_default_position.y,
								  eye_default_position.z };		//�ӵ�λ��
my_3D_point_coord light_position(light_default_position.x,
								light_default_position.y,
								light_default_position.z);		//��Դ1λ��
my_3D_point_coord light1_position(light1_default_position.x,
								light1_default_position.y,
								light1_default_position.z);		//��Դ2λ��

float light_rgb_ambient[] = { 0.9, 0.9, 0.9 };
float light1_rgb_ambient[] = { 0.5, 0.05, 0.05 };
float light_rgb_diffuse_specular[] = { diffuse, diffuse, diffuse };
float light1_rgb_diffuse_specular[] = { diffuse, diffuse, diffuse };

int selected = 1;												//��ѡ������						//	1 2 3
float renderoffset = 0;											//����΢�� �� �����׷����Ч��		//	4
int sighttype = 0;												//�����Ƿ��������					//	5
int drawaxis = 0;												//�Ƿ�ģ�������ԭ��				//	6

int selectedrgb = 0;											//ѡ�еĹ�Դrgb �� 0r 1g 2b 3����	//	, . /  [\]
int zbuffer = 0;												//�Դ�zbuffer						//	z
bool open_light = true;											//��1����							//	o
bool open_light1 = false;										//��2����							//	o
bool rendered = false;											//����׷��							//	r
int showshadow = 0;												//��дzbuffer�Ƿ���ʾ��Ӱ			//	l
int skipmodel[7] = {1, 1, 1, 1, 1, 1, 1};						//����ģ��							//	b n m - g h j - t

std::vector< my_triangle_3DModel> all_models;					//����������ģ��
std::map<my_3D_point_coord*, my_draw_color*> render_vertices;	//������Ҫ���Ƶĵ��Լ�������

unsigned image_w, image_h;										//��׷ͶӰ��ߴ�
int gap = 0;													//��׷�����ܶ�
void printManual();




#define FAREST 9999
int my_zbuffer = 0;												//��дzbuffer						//	x
float Zbuffer[nearplane_width + 1][nearplane_width + 1];
my_draw_color Framebuffer[nearplane_width + 1][nearplane_width + 1];
int shadowbuffer[nearplane_width + 1][nearplane_width + 1];

// �жϵ��Ƿ���ͶӰ����
int inprojection(my_3D_point_coord p, int width, int height) {
	if (p.x >= -width / 2 && p.x <= width / 2 && p.y >= -height / 2 && p.y <= height / 2)
		return 1;
	return 0;
}

// �жϵ��Ƿ���ͶӰ���ϵ���������
int inface(my_3D_point_coord p, my_3D_point_coord a, my_3D_point_coord b, my_3D_point_coord c) {
	bool flag = false;
	//ɨ���߷���ÿ���ͱ��ཻ�����ж���
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

// �жϵ��Ƿ�����ά�ռ��е���������
int inface3D(my_3D_point_coord p, my_3D_point_coord a, my_3D_point_coord b, my_3D_point_coord c) {
	//�����
	//����С��������� 1/2*a*b*sin<a,b>
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

	//һ�������������
	my_3Dvector Q,R;
	M = my_3Dvector(a, b);
	N = my_3Dvector(a, c);
	Q = M; Q.normalized();
	R = N; R.normalized();
	float S, angle;
	angle = acosf(Q.dot(R));
	S = 0.5 * M.len * N.len * sinf(angle);

	//��������� Ӳ��Ӱ����
	if (fabsf(S1 + S2 + S3 - S) < 0.015)
		return 1;
	////�����Ƴ̶Ƚ��� ���Ե����
	//else if (fabsf(S1 + S2 + S3 - S) >= 0.5 && fabsf(S1 + S2 + S3 - S) < 100)
	//	return 2;
	return 0;
	/*//��˷�
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

// �жϵ��Ƿ�͹�Դ֮�����ڵ�
void findBarrier(my_3D_point_coord start, my_3D_point_coord end, my_3D_point_coord proj, 
	int model_index, int face_index, my_draw_color formal_color, int lightcount) {
	int nobarrier = 1;
	for (int i = 0; i < all_models.size(); i++) {
		if (skipmodel[i] == 0) continue;
		for (int j = 0; j < all_models[i].faceSets.size(); j++) {
			if (i == model_index && j == face_index) continue;//��������start������ģ�ͣ�����

			int firstPointIndex = all_models[i].faceSets[j].first_point_index;
			int secondPointIndex = all_models[i].faceSets[j].second_point_index;
			int thirdPointIndex = all_models[i].faceSets[j].third_point_index;
			my_3D_point_coord p4 = all_models[i].pointSets[firstPointIndex];
			my_3D_point_coord p5 = all_models[i].pointSets[secondPointIndex];
			my_3D_point_coord p6 = all_models[i].pointSets[thirdPointIndex];

			//�� �ƹ���start������֮�� ���ƽ��Ľ�������
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

			my_3Dvector a1(start, res); // �ж�start��end��������
			my_3Dvector a2(end, res);

			//��һ����
			if (lightcount == 1) {
				//�� �������������� �� ͶӰ��û��ʩ�ӹ���Ӱ
				if (inface3D(res, p4, p5, p6) && a1.dot(a2) < 0
					&& shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] == 0) {
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].r -= 0.4;
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].g -= 0.4;
					Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].b -= 0.4;
					shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = 1;
				}
			}
			//�ڶ����� �жϵڶ����ƺ͵�����֮����û���ڵ���
			else if (lightcount == 2) {
				if (inface3D(res, p4, p5, p6) && a1.dot(a2) < 0) {
					nobarrier = 0;
					break;
				}
			}
		}
		if (!nobarrier && lightcount == 2) break; //�����ڶ��������ڵ��� ����ѭ��
	}
	//�����ڶ�����û���ڵ��� �õ���ɫ��Ϊ����Ӱ��ɫ
	if (nobarrier && lightcount == 2) {
		if (shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] == 1) {
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].r = formal_color.r;
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].g = formal_color.g;
			Framebuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1].b = formal_color.b;
			shadowbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = 0;
		}
	}
}

// ����ͶӰ�����ʵ����
void finddepth_BPRM_shadow(my_3D_point_coord eye, my_3D_point_coord proj, int model_index, int face_index) {
	//����ʵ����res
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

	//�� res��zbuffer�и������ӵ�
	if (Zbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] > pointdepth) {
		Zbuffer[(int)proj.x + nearplane_width / 2 + 1][(int)proj.y + nearplane_width / 2 + 1] = (float)pointdepth;

		N.normalized();
		my_draw_color pcolor, colortemp;
		pcolor.r = pcolor.g = pcolor.b = 0;

		//����ʵ����ȡ����������
		int firstPointIndex = all_models[model_index].faceSets[face_index].first_point_index;//ȡ����������
		int secondPointIndex = all_models[model_index].faceSets[face_index].second_point_index;
		int thirdPointIndex = all_models[model_index].faceSets[face_index].third_point_index;

		my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex];//��һ������
		my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //�ڶ�������
		my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //����������

		////glClearColor(1.f, 1.f, 1.f, 0.f);
		////glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		////��������ӳ�䷽ʽ�����ģʽ��
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		////ʵʩ��������ͼ�������ӳ��
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
		////glBegin(GL_POINTS);//��ʼ����
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
			//���Ҹõ�͹�Դ֮���Ƿ����ڵ��� ���� �õ����Ƚ���
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

// ��дzbuffer
void My_zbuffer() {
	//��ʼ��zbuffer��framebuffer��shadowbuffer
	for (int i = 0; i <= nearplane_width; i++) {
		for (int j = 0; j <= nearplane_height; j++) {
			Zbuffer[i][j] = FAREST;
			Framebuffer[i][j].r = Framebuffer[i][j].g = Framebuffer[i][j].b = 0;
			shadowbuffer[i][j] = 0;
		}
	}

	//���zbuffer��framebuffer
	for (int model_index = 0; model_index < all_models.size(); model_index++) {
		if (skipmodel[model_index] == 0) continue;
		for (int face_index = 0; face_index < all_models[model_index].faceSets.size(); face_index++) {
			//�ȱ����޳�
			eyesight = my_3Dvector(eye_position,
				all_models[model_index].pointSets[all_models[model_index].faceSets[face_index].first_point_index]);
			if (all_models[model_index].faceSets[face_index].n.dot(eyesight) > 0) continue;

			//����ʵ����ȡ����������
			int firstPointIndex = all_models[model_index].faceSets[face_index].first_point_index; //ȡ����������
			int secondPointIndex = all_models[model_index].faceSets[face_index].second_point_index;
			int thirdPointIndex = all_models[model_index].faceSets[face_index].third_point_index;

			int firstNormalIndex = all_models[model_index].faceSets[face_index].first_point_normal_index; //ȡ������������
			int secondNormalIndex = all_models[model_index].faceSets[face_index].second_point_normal_index;
			int thirdNormalIndex = all_models[model_index].faceSets[face_index].third_point_normal_index;

			my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex]; //��һ������
			my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //�ڶ�������
			my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //����������

			my_3Dvector p1Normal = all_models[model_index].pointNormalSets[firstNormalIndex];//��һ�����㷨��
			my_3Dvector p2Normal = all_models[model_index].pointNormalSets[secondNormalIndex];//�ڶ������㷨��
			my_3Dvector p3Normal = all_models[model_index].pointNormalSets[thirdNormalIndex];//���������㷨��

			my_draw_color p1color, p2color, p3color, colortemp;
			p1color.r = p1color.g = p1color.b = 0;
			p2color.r = p2color.g = p2color.b = 0;
			p3color.r = p3color.g = p3color.b = 0;
			if (open_light) {
				//��Blinn-Phong Reflection Model����ÿһ�㵽��Դλ�õ�����
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

			if (inprojection(p1proj, nearplane_width, nearplane_height) //�� ��Ķ�����ͶӰ���� && zbuffer���ڶ�Ӧλ��С����ֵ
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

			for (int i = 1; i <= nearplane_width; i++) //����ͶӰ��ÿһ������
				for (int j = 1; j <= nearplane_height; j++)
					//�� ���������� ���ڴ����е���ķ�Χ��
					if (inface(my_3D_point_coord(i - nearplane_width / 2 - 1, j - nearplane_height / 2 - 1, 200), p1proj, p2proj, p3proj))
						//�� �õ����
						finddepth_BPRM_shadow(eye_position, my_3D_point_coord(i - nearplane_width / 2 - 1, j - nearplane_height / 2 - 1, 200),
							model_index, face_index);
			printf("�� %d / %d ��ģ�ͣ��� %d / %d ����zbuffer���\t(������δ���Ǳ����޳�)\n",
				model_index + 1, all_models.size(), face_index + 1, all_models[model_index].faceSets.size());
		}
	}
	printManual();
}




// ��ӡ˵����
void printManual() {
	cout << "\n\t\t\t\b\bManual" << endl;
	cout << "\t/--------------------------------|" << endl;
	cout << "\t|[1] ���������\t\t\t |" << endl;
	cout << "\t|[2] ������1\t\t\t |" << endl;
	cout << "\t|[3] ������2\t\t\t |" << endl;
	cout << "\t|[4] �޸���׷����Ч��\t\t |" << endl;
	cout << "\t|[5] ���߸�������\t\t |" << endl;
	cout << "\t|[6] ��ʾ��ת��\t\t\t |" << endl;
	cout << "\t|[`] ��ʾ��˵����\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[w][a][s][d] ����ƽ��\t\t |" << endl;
	cout << "\t|[q][e] ������ת\t\t |" << endl;
	cout << "\t|[i][k] ��������\t\t |" << endl;
	cout << "\t|[z] �����Դ�zbuffer\t\t |" << endl;
	cout << "\t|[x] ������дzbuffer\t\t |" << endl;
	cout << "\t|[o] ���ص�\t\t\t |" << endl;
	cout << "\t|[p] ���帴λ\t\t\t |" << endl;
	cout << "\t|[r] ���ع�׷\t\t\t |" << endl;
	cout << "\t|[l] ������Ӱ��ʾ\t\t |" << endl;
	cout << "\t|[y] ��������\t\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[t][ ][ ]\t\t\t |" << endl;
	cout << "\t|[g][h][j] ����ģ����ʾ\t\t |" << endl;
	cout << "\t|[b][n][m]\t\t\t |" << endl;
	cout << "\t|--------------------------------|" << endl;
	cout << "\t|[,][.][/] �л�����rgb\t\t |" << endl;
	cout << "\t|[\\] �л�����������\t\t |" << endl;
	cout << "\t|[[][]] ����rgb\t\t\t |" << endl;
	cout << "\t|--------------------------------/\n" << endl;
}

// ��ӡ��׷����
int printprogress(long double i) {
	printf("����� [ ");
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

// ģ��ƽ��
void translate(float delta, int dir, int light) {
	/*ƽ�������*/
	if (light == 0) {
		//ƽ��ģ��
		for (int model_index = 0; model_index < all_models.size(); model_index++) {
			for (int i = 0; i < all_models[model_index].pointSets.size(); i++) {
				if (dir == 1)/*123456ǰ����������*/ {
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
		//ƽ�Ƶƺ�ԭ��
		if (dir == 1) /*123456ǰ����������*/ {
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
	/* ƽ�Ƶ�1  */
	else if (light == 1) {
		if (dir == 1) /*123456ǰ����������*/
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
	/* ƽ�Ƶ�2  */
	else if (light == 2) {
		if (dir == 1) /*123456ǰ����������*/
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

// �������
my_3D_point_coord matrix_multiply_vector(float matrix[][4], my_3D_point_coord input_v) {
	my_3D_point_coord translated_v;
	translated_v.x = matrix[0][0] * input_v.x + matrix[0][1] * input_v.y + matrix[0][2] * input_v.z + matrix[0][3] * 1;
	translated_v.y = matrix[1][0] * input_v.x + matrix[1][1] * input_v.y + matrix[1][2] * input_v.z + matrix[1][3] * 1;
	translated_v.z = matrix[2][0] * input_v.x + matrix[2][1] * input_v.y + matrix[2][2] * input_v.z + matrix[2][3] * 1;
	return translated_v;
}

// ģ����ת
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

// ��ʼ������ģ��
void init(void) {
	//����Թ��ߵķ�����
	float material_ambient_rgb_reflection[] = { 0.2, 0.2, 0.2 };
	float material_specular_rgb_reflection[] = { 0.2, 0.2, 0.2 };
	float ns = 40; //�۹�ָ��
	//����ģ��
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

// ��������
void display(void) {
	//����淨��
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

	// ��дzbuffer + shadow����
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

	// ֱ�ӹ����
	else if (!rendered && !my_zbuffer) {
		glShadeModel(GL_SMOOTH);
		for (unsigned int model_index = 0; model_index < all_models.size(); model_index++)
		{
			if (skipmodel[model_index] == 0) continue;
			for (unsigned int i = 0; i < all_models[model_index].faceSets.size(); i++)
			{
				//�ȱ����޳�
				eyesight = my_3Dvector(eye_position,
					all_models[model_index].pointSets[all_models[model_index].faceSets[i].first_point_index]);
				if (all_models[model_index].faceSets[i].n.dot(eyesight) > 0)
					continue;

				//����ʵ����ȡ����������
				int firstPointIndex = all_models[model_index].faceSets[i].first_point_index;//ȡ����������
				int secondPointIndex = all_models[model_index].faceSets[i].second_point_index;
				int thirdPointIndex = all_models[model_index].faceSets[i].third_point_index;

				int firstNormalIndex = all_models[model_index].faceSets[i].first_point_normal_index; //ȡ������������
				int secondNormalIndex = all_models[model_index].faceSets[i].second_point_normal_index;
				int thirdNormalIndex = all_models[model_index].faceSets[i].third_point_normal_index;

				my_3D_point_coord p1 = all_models[model_index].pointSets[firstPointIndex];//��һ������
				my_3D_point_coord p2 = all_models[model_index].pointSets[secondPointIndex]; //�ڶ�������
				my_3D_point_coord p3 = all_models[model_index].pointSets[thirdPointIndex]; //����������

				my_3Dvector p1Normal = all_models[model_index].pointNormalSets[firstNormalIndex];//��һ�����㷨��
				my_3Dvector p2Normal = all_models[model_index].pointNormalSets[secondNormalIndex];//�ڶ������㷨��
				my_3Dvector p3Normal = all_models[model_index].pointNormalSets[thirdNormalIndex];//���������㷨��

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

					//��Blinn-Phong Reflection Model����ÿһ�㵽��Դλ�õ�����
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
				// ���������ģ��
				if (model_index != 5 && model_index != 3 && model_index != 1 && model_index != 0) {
					glBegin(GL_TRIANGLES);//��ʼ����
					glColor3f(p1color.r, p1color.g, p1color.b);
					glVertex3f(p1.x, p1.y, p1.z);
					glColor3f(p2color.r, p2color.g, p2color.b);
					glVertex3f(p2.x, p2.y, p2.z);
					glColor3f(p3color.r, p3color.g, p3color.b);
					glVertex3f(p3.x, p3.y, p3.z);
					glEnd();
				}
				// �������ģ��
				else{
					//��������ӳ�䷽ʽ�����ģʽ��
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
					glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
					
					//ʵʩ��������ͼ�������ӳ��
					if (model_index == 0)
						glBindTexture(GL_TEXTURE_2D, texGround0);
					else if (model_index == 1)
						glBindTexture(GL_TEXTURE_2D, texGround1);
					else if (model_index == 3)
						glBindTexture(GL_TEXTURE_2D, texGround3);
					else if (model_index == 5)
						glBindTexture(GL_TEXTURE_2D, texGround5);

					glBegin(GL_TRIANGLES);//��ʼ����
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

	// ���߸��ٻ���
	else if (rendered && !my_zbuffer) {
		glBegin(GL_POINTS);//��ʼ����
		for (std::map<my_3D_point_coord*, my_draw_color*>::iterator piter = render_vertices.begin(); piter != render_vertices.end(); piter++) {
			glColor3f(piter->second->r, piter->second->g, piter->second->b);
			glVertex3f(piter->first->x, piter->first->y, piter->first->z);
		}
		glEnd();
	}

	if (drawaxis) {
		//��������
		glPointSize(5);
		glBegin(GL_POINTS);
		glColor4f(1, 1, 0, 0.2);
		for (int i = -300; i <= 300; i++)
			glVertex3f(original.x, original.y + i, original.z);
		glEnd();

		//����original
		glPointSize(8);
		glBegin(GL_POINTS);
		glColor4f(1, 0, 1, 0.2);
		glVertex3f(original.x, original.y, original.z);
		glEnd();
	}

	glutSwapBuffers();
}

// ͶӰ��ʽ��modelview��ʽ���á���ͶӰ�����
void reshape(int w, int h) {
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//if (w <= h) {//����ƽ��ͶӰ�Ӿ���
	//	glOrtho(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
	//		(GLfloat)nearplane_height / (GLfloat)nearplane_width, ratio * nearplane_height *
	//		(GLfloat)nearplane_height / (GLfloat)nearplane_width,
	//		nearplane_distance, farplane_distance); //������ӵ�
	//}
	//else {//����ƽ��ͶӰ�Ӿ���
	//	glOrtho(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
	//		(GLfloat)nearplane_width / (GLfloat)nearplane_height, ratio * nearplane_height *
	//		(GLfloat)nearplane_width / (GLfloat)nearplane_height,
	//		nearplane_distance, farplane_distance);
	//}

	if (w <= h) {//����͸���Ӿ���
		glFrustum(-ratio * nearplane_width, ratio * nearplane_width, -ratio * nearplane_height *
			(GLfloat)nearplane_height / (GLfloat)nearplane_width, ratio * nearplane_height *
			(GLfloat)nearplane_height / (GLfloat)nearplane_width,
			nearplane_distance, farplane_distance); //������ӵ�
	}
	else {//����͸���Ӿ���
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

// ���̽����¼�
void keyboard(unsigned char key, int x, int y) {
	switch (key)
	{
	case '1':			/*    �����  */ {
		selected = 1;
		printf("��ѡ�� [�����]\n");
		break;
	}
	case '2':			/*     ��1    */ {
		selected = 2;
		printf("��ѡ�� [��1]\n");
		break;
	}
	case '3':			/*     ��2    */ {
		selected = 3;
		printf("��ѡ�� [��2]\n");
		break;
	}
	case '4':			/* ��Ⱦ��ƫ�� */ {
		if (renderoffset == 0) {
			renderoffset = 9.9;
			printf("��׷΢�� [ON]\n");
		}
		else {
			renderoffset = 0;
			printf("��׷΢�� [OFF]\n");
		}
		break;
	}
	case '5':			/*  �ӵ㳯��  */ {
		if (sighttype == 1) {
			sighttype = 0;
			printf("�������� [OFF]\n");
		}
		else {
			sighttype = 1;
			printf("�������� [ON]\n");
		}
		reshape(nearplane_width, nearplane_height);
		glutPostRedisplay();
		break;
	}
	case '6':			/*  �������  */ {
		if (drawaxis)
			drawaxis = 0;
		else
			drawaxis = 1;
		glutPostRedisplay();
		break;
	}
	case '`':			/*  ��ʾָ��  */ {
		system("cls");
		printManual();
		break;
	}
	case 'q': case 'Q': /* �������ת */ {
		if (selected == 1) {
			rotate(-theta, 0);//��ʱ����ת
			rotateoffset -= theta;
			glutPostRedisplay();
		}
		break;
	}
	case 'e': case 'E': /* �������ת */ {
		if (selected == 1) {
			rotate(theta, 0);//˳ʱ����ת
			rotateoffset += theta;
			glutPostRedisplay();
		}
		break;
	}
	case 'w': case 'W': /*    ǰ��    */ {
		if (selected == 1) {//�����
			translate(delta, 2, 0);
			cameraoffset.z += delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 1, 1);
			lightoffset.z -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 1, 2);
			light1offset.z -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 's': case 'S': /*    ����    */ {
		if (selected == 1) {//�����
			translate(delta, 1, 0);
			cameraoffset.z -= delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 2, 1);
			lightoffset.z += 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 2, 2);
			light1offset.z += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'a': case 'A': /*    ����    */ {
		if (selected == 1) {//�����
			translate(delta, 4, 0);
			cameraoffset.x += delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 3, 1);
			lightoffset.x -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 3, 2);
			light1offset.x -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'd': case 'D': /*    ����    */ {
		if (selected == 1) {//�����
			translate(delta, 3, 0);
			cameraoffset.x -= delta;
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 4, 1);
			lightoffset.x += 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 4, 2);
			light1offset.x += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'i': case 'I': /*    ����    */ {
		if (selected == 1) {//�����
			translate(delta - renderoffset, 6, 0);
			cameraoffset.y -= (delta - renderoffset);
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 5, 1);
			lightoffset.y += 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 5, 2);
			light1offset.y += 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'k': case 'K': /*    ����    */ {
		if (selected == 1) {//�����
			translate(delta - renderoffset, 5, 0);
			cameraoffset.y += (delta - renderoffset);
			reshape(nearplane_width, nearplane_height);
		}
		else if (selected == 2 && open_light == true) {//��1
			translate(10, 6, 1);
			lightoffset.y -= 10;
		}
		else if (selected == 3 && open_light1 == true) {//��2
			translate(10, 6, 2);
			light1offset.y -= 10;
		}
		glutPostRedisplay();
		break;
	}
	case 'p': case 'P': /*    ��λ    */ {
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
	case 'o': case 'O': /*   ���ص�   */ {
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
			glEnable(GL_DEPTH_TEST); //����Ȼ������
			glDepthFunc(GL_LESS); //�ж��ڵ���ϵʱ�����ӵ���������ڵ����ӵ�Զ������	
			zbuffer -= 2;
		}
		else if (!zbuffer) {
			glDisable(GL_DEPTH_TEST); //�ر���Ȼ������
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
	case 'l': case 'L': /*  ��Ӱ��ʾ  */ {
		if (!showshadow) {
			showshadow = 1;
			printf("��Ӱ��ʾ [ON]\n");
		}
		else if (showshadow) {
			showshadow = 0;
			printf("��Ӱ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'r': case 'R': /*  ����׷��  */ {
		//���ݵ�ǰ�Ӿ������ö�ͶӰƽ����ж������
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
	case 'b': case 'B': /* ����ģ��1  */ {
		if (skipmodel[0] != 1) {
			skipmodel[0] = 1;
			printf("ģ��1 [��Ƥ] ��ʾ [ON]\n");
		}
		else if (skipmodel[0] == 1) {
			skipmodel[0] = 0;
			printf("ģ��1 [��Ƥ] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'n': case 'N': /* ����ģ��2  */ {
		if (skipmodel[1] != 1) {
			skipmodel[1] = 1;
			printf("ģ��2 [��ҳ] ��ʾ [ON]\n");
		}
		else if (skipmodel[1] == 1) {
			skipmodel[1] = 0;
			printf("ģ��2 [��ҳ] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'm': case 'M': /* ����ģ��3  */ {
		if (skipmodel[2] != 1) {
			skipmodel[2] = 1;
			printf("ģ��3 [����] ��ʾ [ON]\n");
		}
		else if (skipmodel[2] == 1) {
			skipmodel[2] = 0;
			printf("ģ��3 [����] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'g': case 'G': /* ����ģ��4  */ {
		if (skipmodel[3] != 1) {
			skipmodel[3] = 1;
			printf("ģ��4 [̨��] ��ʾ [ON]\n");
		}
		else if (skipmodel[3] == 1) {
			skipmodel[3] = 0;
			printf("ģ��4 [̨��] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'h': case 'H': /* ����ģ��5  */ {
		if (skipmodel[4] != 1) {
			skipmodel[4] = 1;
			printf("ģ��5 [Ǧ��] ��ʾ [ON]\n");
		}
		else if (skipmodel[4] == 1) {
			skipmodel[4] = 0;
			printf("ģ��5 [Ǧ��] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'j': case 'J': /* ����ģ��6  */ {
		if (skipmodel[5] != 1) {
			skipmodel[5] = 1;
			printf("ģ��6 [����] ��ʾ [ON]\n");
		}
		else if (skipmodel[5] == 1) {
			skipmodel[5] = 0;
			printf("ģ��6 [����] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 't': case 'T': /* ����ģ��7  */ {
		if (skipmodel[6] != 1) {
			skipmodel[6] = 1;
			printf("ģ��7 [����] ��ʾ [ON]\n");
		}
		else if (skipmodel[6] == 1) {
			skipmodel[6] = 0;
			printf("ģ��7 [����] ��ʾ [OFF]\n");
		}
		glutPostRedisplay();
		break;
	}
	case 'y': case 'Y': /*  ��������  */ {
		if (!texenable) {
			texenable = 1;
			glEnable(GL_TEXTURE_2D);
			texGround0 = load_texture("tex\\book.bmp");  //��������
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
	case ',':		    /* �ƹ�ѡ��r  */ {
		selectedrgb = 0;
		printf("������ѡ�� [R]\n");
		break;
	}
	case '.':		    /* �ƹ�ѡ��g  */ {
		selectedrgb = 1;
		printf("������ѡ�� [G]\n");
		break;
	}
	case '/':		    /* �ƹ�ѡ��b  */ {
		selectedrgb = 2;
		printf("������ѡ�� [B]\n");
		break;
	}
	case '\\':		    /*�ƹ�ѡ������*/ {
		selectedrgb = 3;
		printf("������ѡ�� [����]\n");
		break;
	}
	case '[':		    /*�ϵ���rgbֵ */ {
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
	case ']':		    /*�µ���rgbֵ */ {
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
	case 27:			/*  �˳�����  */ {
		exit(0);
		break;
	}
	}
}

// ��꽻���¼�
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

// ��������
int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(nearplane_width, nearplane_height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("�������� - ��ĩ����ҵ");
	init();
	printManual();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();
	return 0;
}