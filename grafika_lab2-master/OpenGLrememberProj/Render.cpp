#include "Render.h"

#include <sstream>
#include <iostream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"
#include <math.h>
#include <vector>
#include <array>
std::array<double, 3> normal(double*, double*, double*);
void shep();
void povbezie();
double fact(int);
void vorota();

double a = sqrt(0.5);   //���������� �� � ������� �����������
double b = sqrt(0.5);	//���������� �� � ������� �����������
double t = 0;			//���� ���������
double alfa = atan2(a, b); //���� ����������� � ��������
int flag = 0;
int flagh = 0;
int flagp = 0;
int flago = 0;

bool textureMode = true;
bool lightMode = true;

//����� ��� ��������� ������
class CustomCamera : public Camera
{
public:
	//��������� ������
	double camDist;
	//���� �������� ������
	double fi1, fi2;

	
	//������� ������ �� ���������
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//������� ������� ������, ������ �� ����� ��������, ���������� �������
	void SetUpCamera()
	{
		//�������� �� ������� ������ ������
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//������� ��������� ������
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //������� ������ ������


//����� ��� ��������� �����
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//��������� ������� �����
		pos = Vector3(1, 1, 3);
	}

	
	//������ ����� � ����� ��� ���������� �����, ���������� �������
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//����� �� ��������� ����� �� ����������
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//������ ���������
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// ��������� ��������� �����
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// �������������� ����������� �����
		// ������� ��������� (���������� ����)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// ��������� ������������ �����
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// ��������� ���������� ������������ �����
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //������� �������� �����




//������ ���������� ����
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//������ ���� ������ ��� ������� ����� ������ ����
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//������� ���� �� ���������, � ����� ��� ����
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
	if (key == VK_LEFT) //�����
	{
		a = a*cos( -3.14 / 180)+b* sin(  -3.14 / 180);
		b = -a*sin( -3.14 / 180) + b* cos(-3.14 / 180); 
	}
	if (key == VK_UP) //�����
	{
		t+=0.5;
	}
	if (key == VK_RIGHT) //������
	{
		a = a * cos( 3.14 / 180) + b * sin( 3.14 / 180);
		b = -a * sin( 3.14 / 180) + b * cos( 3.14 / 180);
	}
	if (key == VK_DOWN) //����
	{
		if (t > 0)
		{
			t-=0.5;
		}
		else
			t = 0;
	}
	if (key == 'H')
	{
		flagh=1;
	}
	if (key == 'P')
	{
		flagp = 1;
	}
	if (key == 'O')
	{
		flago = 1;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	if (key == 'H')
	{
		flagh = 0;
	}
	if (key == 'P')
	{
		flagp = 0;
	}
	if (key == 'O')
	{
		flago = 0;
	}

}


GLuint texId;
GLuint tex2Id;
GLuint tex3Id;


//����������� ����� ������ ��������
void initRender(OpenGL *ogl)
{
	
	
	
	
	//��������� �������

	//4 ����� �� �������� �������
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//��������� ������ ��������� �������
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//�������� ��������
	glEnable(GL_TEXTURE_2D);
	



	//������ ����������� ���������  (R G B)
	RGBTRIPLE* texarray3;

	//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
	char* texCharArray3;
	int texW3, texH3;
	OpenGL::LoadBMP("watertexture.bmp", &texW3, &texH3, &texarray3);
	OpenGL::RGBtoChar(texarray3, texW3, texH3, &texCharArray3);



	//���������� �� ��� ��������
	glGenTextures(1, &tex3Id);
	//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
	glBindTexture(GL_TEXTURE_2D, tex3Id);

	//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW3, texH3, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray3);

	//�������� ������
	free(texCharArray3);
	free(texarray3);

	//������� ����
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);










	//������ ����������� ���������  (R G B)
	RGBTRIPLE* texarray2;

	//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
	char* texCharArray2;
	int texW2, texH2;
	OpenGL::LoadBMP("coraltexture.bmp", &texW2, &texH2, &texarray2);
	OpenGL::RGBtoChar(texarray2, texW2, texH2, &texCharArray2);



	//���������� �� ��� ��������
	glGenTextures(1, &tex2Id);
	//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
	glBindTexture(GL_TEXTURE_2D, tex2Id);

	//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW2, texH2, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray2);

	//�������� ������
	free(texCharArray2);
	free(texarray2);

	//������� ����
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);





	//������ ����������� ���������  (R G B)
	RGBTRIPLE *texarray;

	//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
	char *texCharArray;
	int texW, texH;
	OpenGL::LoadBMP("vorotatexture.bmp", &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	
	
	//���������� �� ��� ��������
	glGenTextures(1, &texId);
	//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
	glBindTexture(GL_TEXTURE_2D, texId);

	//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//�������� ������
	free(texCharArray);
	free(texarray);

	//������� ����
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);



	





	//������ � ���� ����������� � "������"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// ������������ �������� : �� ����� ����� ����� 1
	glEnable(GL_NORMALIZE);

	// ���������� ������������� ��� �����
	glEnable(GL_LINE_SMOOTH); 


	//   ������ ��������� ���������
	//  �������� GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  ������� � ���������� �������� ���������(�� ���������), 
	//                1 - ������� � ���������� �������������� ������� ��������       
	//                �������������� ������� � ���������� ��������� ����������.    
	//  �������� GL_LIGHT_MODEL_AMBIENT - ������ ������� ���������, 
	//                �� ��������� �� ���������
	// �� ��������� (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}




double shep1[3] ;  //����� ��� ���������� �������
double shep2[3] ;
double shep3[3] ;
double shep4[3] ;

double shepcentr[3] = { 0,0,0 };//����� �������





void Render(OpenGL *ogl)
{
	

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//��������������
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//��������� ���������
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//�������
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//��������
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//����������
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//������ �����
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//���� ���� �������, ��� ����������� (����������� ���������)
	glShadeModel(GL_SMOOTH);
	//===================================
	//������� ���  


	
	double sea1[3] = { 20,20,0 };
	double sea2[3] = { -20,20,0 };
	double sea3[3] = { -20,-20,0 };
	double sea4[3] = { 20,-20,0 };

	
	
	povbezie();
	flag = 1;
	

	
	
	
	shep();

	vorota();
	
	//"����"
	
	glBindTexture(GL_TEXTURE_2D, tex3Id);

	glBegin(GL_QUADS);

	glNormal3d(0, 0, 1);
	glColor3d(0.7, 0.2, 0);
	glTexCoord2d(1, 1);
	glVertex3dv(sea1);
	glTexCoord2d(0, 1);
	glVertex3dv(sea2);
	glTexCoord2d(0, 0);
	glVertex3dv(sea3);
	glTexCoord2d(1, 0);
	glVertex3dv(sea4);

	glEnd;
	
	
	
	
	
	


   //��������� ������ ������

	
	glMatrixMode(GL_PROJECTION);	//������ �������� ������� ��������. 
	                                //(���� ��������� ��������, ����� �� ������������.)
	glPushMatrix();   //��������� ������� ������� ������������� (������� ��������� ������������� ��������) � ���� 				    
	glLoadIdentity();	  //��������� ��������� �������
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //������� ����� ������������� ��������

	glMatrixMode(GL_MODELVIEW);		//������������� �� �����-��� �������
	glPushMatrix();			  //��������� ������� ������� � ���� (��������� ������, ����������)
	glLoadIdentity();		  //���������� �� � ������

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //������� ����� ��������� ��� ������� ������ � �������� ������.
	rec.setSize(300, 150);
	rec.setPosition(10, ogl->getHeight() - 150 - 10);


	std::stringstream ss;
	ss << "T - ���/���� �������" << std::endl;
	ss << "L - ���/���� ���������" << std::endl;
	ss << "F - ���� �� ������" << std::endl;
	ss << "G - ������� ���� �� �����������" << std::endl;
	ss << "G+��� ������� ���� �� ���������" << std::endl;
	ss << "�����. �����: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "�����. ������: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "��������� ������: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //��������������� ������� �������� � �����-��� �������� �� �����.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();


	

}



std::array<double, 3> normal(double* T1, double* T2, double* T3)
{
	double ax, ay, az, bx, by, bz, xt1, xt2, xt3, yt1, yt2, yt3, zt1, zt2, zt3, nx, ny, nz, r;
	std::array<double, 3> N;
	xt1 = T1[0];
	yt1 = T1[1];
	zt1 = T1[2];
	xt2 = T2[0];
	yt2 = T2[1];
	zt2 = T2[2];
	xt3 = T3[0];
	yt3 = T3[1];
	zt3 = T3[2];
	ax = xt3 - xt1;
	ay = yt3 - yt1;
	az = zt3 - zt1;
	bx = xt2 - xt1;
	by = yt2 - yt1;
	bz = zt2 - zt1;

	nx = ay * bz - az * by;
	ny = ax * bz - az * bx;
	nz = ax * by - ay * bx;

	r = sqrt(nx * nx + ny * ny + nz * nz);
	nx = nx / r;
	ny = ny / r;
	nz = nz / r;

	N[0] = nx;
	N[1] = -ny;
	N[2] = nz;
	return N;
}




void shep()
{
	
	std::array<double, 3>N;

	shepcentr[0] += a *  t / 10;
	shepcentr[1] += b *  t /10;
	
	
	glBegin(GL_QUADS);
	
	shep1[0] = shepcentr[0] +b/2;
	shep1[1] = shepcentr[1] -a/2;
	shep1[2] = 1;

	shep2[0] = shepcentr[0] -b/2;
	shep2[1] = shepcentr[1] +a/2;
	shep2[2] = 1;

	shep3[0] = shepcentr[0] +2*a+b/2;
	shep3[1] = shepcentr[1] +2*b-a/2;
	shep3[2] = 1;

	shep4[0] = shepcentr[0] +2*a-b/2;
	shep4[1] = shepcentr[1] +2*b+a/2;
	shep4[2] = 1;


	N = normal(shep1, shep2, shep3);


	glNormal3dv(N._Elems);
	glColor3d(1, 0.2, 1);
	glVertex3dv(shep2);
	glVertex3dv(shep1);
	glVertex3dv(shep3);
	glVertex3dv(shep4);

	shep1[0] = shepcentr[0] ;
	shep1[1] = shepcentr[1] ;
	shep1[2] = 0;

	shep2[0] = shepcentr[0] -b/2;
	shep2[1] = shepcentr[1] +a/2;
	shep2[2] = 1;

	shep3[0] = shepcentr[0] +2*a;
	shep3[1] = shepcentr[1] +2*b;
	shep3[2] = 0;

	shep4[0] = shepcentr[0] +2*a-b/2;
	shep4[1] = shepcentr[1] +2*b+a/2;
	shep4[2] = 1;

	N = normal(shep1, shep3, shep2);


	glNormal3dv(N._Elems);
	glColor3d(1, 0.2, 1);
	glVertex3dv(shep2);
	glVertex3dv(shep1);
	glVertex3dv(shep3);
	glVertex3dv(shep4);



	shep1[0] = shepcentr[0] + b / 2;
	shep1[1] = shepcentr[1] - a / 2;
	shep1[2] = 1;

	shep2[0] = shepcentr[0] ;
	shep2[1] = shepcentr[1] ;
	shep2[2] = 0;

	shep3[0] = shepcentr[0] + 2 * a + b / 2;
	shep3[1] = shepcentr[1] + 2 * b - a / 2;
	shep3[2] = 1;

	shep4[0] = shepcentr[0] +2*a;
	shep4[1] = shepcentr[1] +2*b;
	shep4[2] = 0;

	N = normal(shep1, shep3, shep2);


	glNormal3dv(N._Elems);
	glColor3d(1, 0.2, 1);
	glVertex3dv(shep2);
	glVertex3dv(shep1);
	glVertex3dv(shep3);
	glVertex3dv(shep4);




	glEnd;
}

double ep[400][3];
double e[100][3];


void povbezie()
{
	//������������ ����������� �����
	
	double epuse[] = { 0,0,0 };
	double u = 0;
	double v = 0;
	double x, y;
	int i, j, k, l;
	
	if (flag == 0 || flago==1)
	{
		for (i = 0; i < 100; i++)
		{
			e[i][0] = (-20 + (i % 10) * 4.4);
			e[i][1] = (-20 + ((i/10)%10) * 4.4);
			e[i][2] = (-25) + rand() % 30;
		}
	}
	if (flagh == 1)
	{
		glPointSize(5);
		glColor3d(0, 1, 0);
		glBegin(GL_POINTS);
		for (i = 0; i < 100; i++)
		{
			glVertex3dv(e[i]);

		}

		glEnd();

		glColor3d(0, 1, 0);
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < 100; i++)
		{
			if (i % 10 == 0)
			{
				glEnd();
				glBegin(GL_LINE_STRIP);
			}
			glVertex3dv(e[i]);
		}

		glEnd();

		glColor3d(0, 1, 0);
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < 10; i++)
		{
			for (j = 0; j < 10; j++)
			{
				if ((j + i * 10) % 10 == 0)
				{
					glEnd();
					glBegin(GL_LINE_STRIP);
				}
				glVertex3dv(e[j * 10 + i]);
			}
		}

		glEnd();
	}
	double x1, x2, x3,x4,x5;
	
	for (k = 0; k < 20; k++)
	{
		u = 0.053 * k ;
		for (l = 0; l < 20; l++)
		{
			v = 0.053 * l ;
			for (i = 0; i < 10; i++)
			{
				for (j = 0; j < 10; j++)
				{
					x1 = fact(9);
					x2 = fact(9 - i);
					x3 = fact(9 - j);
					x4 = fact(i);
					x5 = fact(j);

					epuse[0] += ((x1 /(x2*x4)) * pow(u, i) * pow(1 - u, 9 - i)) * ((x1 / (x3 * x5)) * pow(v, j) * pow(1 - v, 9 - j)) * e[i * 10 + j][0];
					epuse[1] += ((x1 / (x2 * x4)) * pow(u, i) * pow(1 - u, 9 - i)) * ((x1 / (x3 * x5)) * pow(v, j) * pow(1 - v, 9 - j)) * e[i * 10 + j][1];
					epuse[2] += ((x1 / (x2 * x4)) * pow(u, i) * pow(1 - u, 9 - i)) * ((x1 / (x3 * x5)) * pow(v, j) * pow(1 - v, 9 - j)) * e[i * 10 + j][2];



				}

			}
			ep[k * 20 + l][0] = epuse[0];
			ep[k * 20 + l][1] = epuse[1];
			ep[k * 20 + l][2] = epuse[2];

			epuse[0] = 0;
			epuse[1] = 0;
			epuse[2] = 0;

		}
	}


	

	


	glPointSize(5);
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (k = 0; k < 400; k++)
	{
		
		glVertex3dv(ep[k]);
	}

	glEnd();

	glColor3d(1, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (k = 0; k < 400; k++)
	{
		if (k % 20 == 0)
		{
			glEnd();
			glBegin(GL_LINE_STRIP);
		}
		glVertex3dv(ep[k]);
	}

	glEnd();

	glColor3d(1, 0, 0);
	glBegin(GL_LINE_STRIP);
	for (k = 0; k < 20; k++)
	{
		for (l = 0; l < 20; l++)
		{
			if ((l + k * 20) % 20 == 0)
			{
				glEnd();
				glBegin(GL_LINE_STRIP);
			}
			glVertex3dv(ep[l * 20 + k]);
		}
	}

	glEnd();


	glBindTexture(GL_TEXTURE_2D, tex2Id);
	glColor3d(1, 1, 0);
	glBegin(GL_QUAD_STRIP);
	for (k = 1; k < 20; k++)
	{
		for (l = 0; l < 20; l++)
		{
			if ((l + k * 20) % 20 == 0)
			{
				glEnd();
				glBegin(GL_QUAD_STRIP);
				
			}
			
			glTexCoord2d((k-1+0.0001)/20, (l+0.0001)/20);
			glVertex3dv(ep[(k-1) * 20 + l]);
			glTexCoord2d((k + 0.01)/20, (l + 0.01)/20);
			glVertex3dv(ep[k * 20 + l]);
		}
	}

	glEnd();



}


double fact(int n)
{
	int b = n;
	int i = 0;
	int r = 1;
	if (b==0)
	{
		return 1;
	}
	else
	{
		for (i = b; i >= 1; i--)
		{
			r *= i;
		}
		return r;
	}
}

double vorotamass[5][3] = { {-10,-10,0}, {-10,10,0}, {10,10,0} ,{10,-10,0},{17,17,0} };
double vorx, vory;
double schet=0;
double flagschet1 = 0;
double flagschet2 = 0;

void vorota()
{
	double tower1[3] = { 1, 1, 0 };
	double tower2[3] = { -1, 1, 0 };
	double tower3[3] = { -1, -1, 0 };
	double tower4[3] = { 1, -1, 0 };

	double vorota1[3] ;
	double vorota2[3];
	double vorota3[3] ;

	


	if (flagp == 1)
	{
		for (int i = 0; i < 5; i++)
		{
			vorx = -20 + rand() % 40;
			vory = -20 + rand() % 40;
			vorotamass[i][0] = vorx;
			vorotamass[i][1] = vory;
			vorotamass[i][2] = 0;
		}
	}

	
	flagschet1++;

	for (int i = 0; i < 5; i++)
	{
		if (shepcentr[0] >= vorotamass[i][0] && shepcentr[0]  <= (vorotamass[i][0]+2) && (shepcentr[1]- vorotamass[i][1]<=0.1)&& (shepcentr[1] - vorotamass[i][1] >= -0.1))
		{
			flagschet2++;
		}
		
	}
	if (flagschet1 == flagschet2 && flagschet1==1)
	{
		schet++;
	}
	else if (flagschet1 == flagschet2+1)
	{
		flagschet1 = 0;
		flagschet2 = 0;
	}
	
	for (int i = 0; i <= schet; i++)
	{
		tower1[2] = i;
		tower2[2] = i;
		tower3[2] = i;
		tower4[2] = i;
		
		
		glBegin(GL_QUADS);

		glNormal3d(0, 0, 1);
		glColor3d(0.7, 0.2, 0);
		glVertex3dv(tower1);
		glVertex3dv(tower2);
		glVertex3dv(tower3);
		glVertex3dv(tower4);

		glEnd;
	}


	for (int i = 0; i < 5; i++)
	{
		vorota1[0] = vorotamass[i][0];
		vorota1[1] = vorotamass[i][1]+0.1;
		vorota1[2] = vorotamass[i][2]+3;

		vorota2[0] = vorotamass[i][0]+2;
		vorota2[1] = vorotamass[i][1]+0.1;
		vorota2[2] = vorotamass[i][2] + 3;

		vorota3[0] = vorotamass[i][0]+2;
		vorota3[1] = vorotamass[i][1]+0.1;
		vorota3[2] = vorotamass[i][2] ;

		vorotamass[i][1] += 0.1;

		glBindTexture(GL_TEXTURE_2D, texId);

	
		glBegin(GL_QUADS);

		glNormal3d(0, 1, 0);
		glColor3d(0.2, 0.2, 0);

		glTexCoord2d(0, 0);
		glVertex3dv(vorotamass[i]);
		glTexCoord2d(0, 1);
		glVertex3dv(vorota1);
		glTexCoord2d(1, 1);
		glVertex3dv(vorota2);
		glTexCoord2d(1, 0);
		glVertex3dv(vorota3);



		vorota1[0] = vorotamass[i][0];
		vorota1[1] = vorotamass[i][1] - 0.1;
		vorota1[2] = vorotamass[i][2] + 3;

		vorota2[0] = vorotamass[i][0] + 2;
		vorota2[1] = vorotamass[i][1] - 0.1;
		vorota2[2] = vorotamass[i][2] + 3;

		vorota3[0] = vorotamass[i][0] + 2;
		vorota3[1] = vorotamass[i][1] - 0.1;
		vorota3[2] = vorotamass[i][2];

		vorotamass[i][1] -= 0.1;

		glNormal3d(0, -1, 0);
		glColor3d(0.2, 0.2, 0);

		glTexCoord2d(0, 0);
		glVertex3dv(vorotamass[i]);
		glTexCoord2d(0, 1);
		glVertex3dv(vorota1);
		glTexCoord2d(1, 1);
		glVertex3dv(vorota2);
		glTexCoord2d(1, 0);
		glVertex3dv(vorota3);




		glEnd;
	}
}