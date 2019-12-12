#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <string.h>

#include <GL/glut.h>


#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))

class point
{
public:
	double x,y,z;
    point(double a,double b,double c)
    {
        x=a;
        y=b;
        z=c;
    }
    point(){    }

    void normalize()
    {
        double val=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        x=x/val;
        y=y/val;
        z=z/val;
    }
    
    point cross(point a, point b)
    {
        point v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    double dot(point a, point b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    double val()
    {
        return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    }

    
};

class color
{
public:
    double r;
    double g;
    double b;
    color(double x,double y,double z)
    {
        r=x;
        g=y;
        b=z;
    }
    color(){ }
};


class Sphere{
public:
    double radius;
    point center;
    color myColor;
    double ambient,diffuse,specular,reflection,shininess;

    Sphere(){ }

    Sphere(point a,double r,color cl)
    {
        center=a;
        radius=r;
        myColor=cl;
    }

    point getNormal(point P)
    {
        point v(P.x-center.x,P.y-center.y,P.z-center.z);
        v.normalize();
        return v;
    }

    void drawSphere()
    {
        point points[100][100];
        int i,j;
        int stacks=20;
        int slices=30;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            glColor3f(myColor.r,myColor.g,myColor.b);
            //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x+center.x,points[i][j].y+center.y,points[i][j].z+center.z);
                    glVertex3f(points[i][j+1].x+center.x,points[i][j+1].y+center.y,points[i][j+1].z+center.z);
                    glVertex3f(points[i+1][j+1].x+center.x,points[i+1][j+1].y+center.y,points[i+1][j+1].z+center.z);
                    glVertex3f(points[i+1][j].x+center.x,points[i+1][j].y+center.y,points[i+1][j].z+center.z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x+center.x,points[i][j].y+center.y,-points[i][j].z+center.z);
                    glVertex3f(points[i][j+1].x+center.x,points[i][j+1].y+center.y,-points[i][j+1].z+center.z);
                    glVertex3f(points[i+1][j+1].x+center.x,points[i+1][j+1].y+center.y,-points[i+1][j+1].z+center.z);
                    glVertex3f(points[i+1][j].x+center.x,points[i+1][j].y+center.y,-points[i+1][j].z+center.z);
                }glEnd();
            }
        }
    }

};

class PyramidBase
{
public:
    point p1,p2,p3,p4;
    color myColor;
    double width;
    double ambient,diffuse,specular,reflection,shininess;

    PyramidBase(point a,point b,point c,point d,color cl)
    {
        p1=a;p2=b;p3=c;p4=d;
        myColor=cl;
    }

    point getNormal()
    {
        point v(0,0,-1);
        return v;
    }

    void drawRectangle()
    {
        glColor3f(myColor.r,myColor.g,myColor.b);
        glBegin(GL_QUADS);{
            glVertex3f(p1.x,p1.y,p1.z);
            glVertex3f(p2.x,p2.y,p2.z);
            glVertex3f(p3.x,p3.y,p3.z);
            glVertex3f(p4.x,p4.y,p4.z);
        }glEnd();
    }
};

class PyramidTriangle
{
public:
    point p1,p2,p3;
    color col;
    double ambient,diffuse,specular,reflection,shininess;

    PyramidTriangle(point a,point b,point c,color cl)
    {
        p1=a;p2=b;p3=c;
        col=cl;
    }

    point getNormal()
    {
        point vec1(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
        point vec2(p3.x-p2.x,p3.y-p2.y,p3.z-p2.z);
        point v=vec1.cross(vec1,vec2);
        v.normalize();
        if(v.z<0){
            v.z=-v.z;
            v.x=-v.x;
            v.y=-v.y;
        } 
        
        return v;
    }

    void drawTriangle()
    {
        glColor3f(col.r,col.g,col.b);
        glBegin(GL_TRIANGLES);{
            glVertex3f(p1.x,p1.y,p1.z);
            glVertex3f(p2.x,p2.y,p2.z);
            glVertex3f(p3.x,p3.y,p3.z);
        }glEnd();
    }
};

class NormalLight
{
public:
    point position;
    double falloff;
    NormalLight(){ }
    NormalLight(point p, double d)
    {
        position=p;
        falloff=d;
    }
};

class SpotLight
{
public:
    point position,lookat;
    double falloff,cutoff;
    SpotLight(){ }
    SpotLight(point p, double d,point look,double cut)
    {
        position=p;
        falloff=d;
        lookat=look;
        cutoff=cut;
    }
};

class CheckerBoard
{
public:
    double cellWidth;
    double ambient,diffuse,reflection,specular,shininess;
    CheckerBoard(){ }
    CheckerBoard(double x)
    {
        cellWidth=x;
    }

    point getNormal()
    {
        point v(0,0,1);
        return v;
    }



    void drawCheckerBoard()
    {
        int num_Cells=ceil(10000/cellWidth);
        point initialPoint(-10000/2,-10000/2,0);

        for (int i=0;i<num_Cells;i++)
        {
            for(int j=0;j<num_Cells;j++)
            {
                if((i+j)%2==0) glColor3f(1,1,1);
                else glColor3f(0,0,0);
                glBegin(GL_QUADS);{
                    glVertex3f(initialPoint.x+cellWidth*i,initialPoint.y+cellWidth*j,initialPoint.z);
                    glVertex3f(initialPoint.x+cellWidth*i,initialPoint.y+cellWidth*(j+1),initialPoint.z);
                    glVertex3f(initialPoint.x+cellWidth*(i+1),initialPoint.y+cellWidth*(j+1),initialPoint.z);
                    glVertex3f(initialPoint.x+cellWidth*(i+1),initialPoint.y+cellWidth*j,initialPoint.z);
                }glEnd();
            }
            
        }
    }
};



bitmap_image b_img ("texture.bmp");

int texture;
color **textureBuffer;
int Texheight, Texwidth;

ofstream outp;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double T;
int globalCount,GPP;

double nearDist,farDist,fovY,aspectRatio,cellWidth,scene_Y,scene_X;
int numRecursion,numPixels,numObjects;
point l_Captured,r_Captured,u_Captured;
point camPos_Captured;

vector<NormalLight> lightSources;
vector<SpotLight> spotLights;
vector<Sphere> sphereObjs;
vector<PyramidBase> pyramidBases;
vector<PyramidTriangle> pyramidTriangles;

point bullets[100];
point pos,u,l,r;
point Bu,Bl,Br,su,sl,sr;
point endGun,bullet;
Sphere S,P;
CheckerBoard checkerboard;




double determinant(point p,point q,point r)
{
	return p.x*(q.y*r.z-q.z*r.y)-q.x*(p.y*r.z-p.z*r.y)+r.x*(p.y*q.z-p.z*q.y);
}

bool illuminatesPoint(point P,point S)
{
	point toSource(S.x-P.x,S.y-P.y,S.z-P.z);
	double tmin=sqrt(pow(toSource.x,2)+pow(toSource.y,2)+pow(toSource.z,2));
	toSource.normalize();
	point BufferPoint(P.x+0.001*toSource.x,P.y+0.001*toSource.y,P.z+0.001*toSource.z);
	point Rd=toSource;

	for(int i=0;i<sphereObjs.size();i++)
	{
		point center=sphereObjs.at(i).center;
		point Ro(BufferPoint.x-center.x,BufferPoint.y-center.y,BufferPoint.z-center.z);
		double a=1;
		double b=2*(Rd.x*Ro.x+Rd.y*Ro.y+Rd.z*Ro.z);
		double c=(Ro.x*Ro.x+Ro.y*Ro.y+Ro.z*Ro.z)-pow(sphereObjs.at(i).radius,2);
		double d= sqrt(pow(b,2)- 4*a*c);
		double t1=(-b+d)/(2*a);
		double t2=(-b-d)/(2*a);
		
		if((pow(b,2)- 4*a*c)<0) continue;

		if((t1>0 && t1<tmin)|| (t2>0 && t2<tmin))
		{
			return false;
		}
		
	}
	
	point Ro=BufferPoint;
	
	for(int i=0;i<pyramidBases.size();i++)
	{
		double Z=pyramidBases.at(i).p1.z;
		double t=(Z-Ro.z)/Rd.z;
		
		point intersection(Ro.x+t*Rd.x,Ro.y+t*Rd.y,Ro.z+t*Rd.z);
		if(intersection.x>=pyramidBases.at(i).p1.x && intersection.x<=pyramidBases.at(i).p1.x+pyramidBases.at(i).width)
		{
			if(intersection.y>=pyramidBases.at(i).p1.y && intersection.y<=pyramidBases.at(i).p1.y+pyramidBases.at(i).width)
			{
				if(t>0 && t<tmin )
				{
					return false;
				}
			}
		}	
		
	}

	for(int i=0;i<pyramidTriangles.size();i++)
	{
		PyramidTriangle temp=pyramidTriangles.at(i);

		point a1(temp.p1.x-temp.p2.x,temp.p1.y-temp.p2.y,temp.p1.z-temp.p2.z);
		point a2(temp.p1.x-temp.p3.x,temp.p1.y-temp.p3.y,temp.p1.z-temp.p3.z);
		point a3(Rd.x,Rd.y,Rd.z);
		point a4(temp.p1.x-Ro.x,temp.p1.y-Ro.y,temp.p1.z-Ro.z);

		double beta=determinant(a4,a2,a3)/determinant(a1,a2,a3);
		double gamma=determinant(a1,a4,a3)/determinant(a1,a2,a3);
		double t=determinant(a1,a2,a4)/determinant(a1,a2,a3);

		double alpha=1-beta-gamma;
		if(beta>=0 && beta<=1 && gamma>=0 && gamma<=1 && alpha>=0 && alpha<=1 && (alpha+beta+gamma)<=1)
		{
			if(t<tmin && t>0 )
			{
				return false;
			}
		}
	}

	//checkerboard
	double t=-Ro.z/Rd.z;
	
	if(t>0 && t<tmin)
	{
		return false;
	}

	return true; 
}


color addTexture(point p)
{
	point lowest;
	lowest.x=floor(p.x/cellWidth)*cellWidth;
	lowest.y=floor(p.y/cellWidth)*cellWidth;
	double diff_x=p.x-lowest.x;
	double diff_y=p.y-lowest.y;
	int x_cord=floor((Texwidth/cellWidth)*diff_x);
	int y_cord=floor((Texheight/cellWidth)*diff_y);
	color pix=textureBuffer[x_cord][y_cord];
	return pix;
}
int global=0;

double farClip(point p)
{
	double dist=p.val();
	double x=(farDist/nearDist)*dist;
	return x-dist;
}


double distance(point a,point b)
{
	double s=pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2);
	return sqrt(s);
}

color generate_pixels(point BufferPoint,point Rd,int depth)
{
	if(depth==0){
		color c(0,0,0);
		return c;
	}

	double maxdist=farClip(Rd);
	Rd.normalize();
	double tmin=100000;
	color pixel(0,0,0);
	point P;
	point normal;
	bool doesIntersect=false;
	bool isSphere=false;
	bool checker=false;
	double ambient,diffuse,specular,reflection,shininess;
	ambient=0;
	 
	for(int i=0;i<sphereObjs.size();i++)
	{
		point center=sphereObjs.at(i).center;
		point Ro(BufferPoint.x-center.x,BufferPoint.y-center.y,BufferPoint.z-center.z);
		double a=1;
		double b=2*(Rd.x*Ro.x+Rd.y*Ro.y+Rd.z*Ro.z);
		double c=(Ro.x*Ro.x+Ro.y*Ro.y+Ro.z*Ro.z)-pow(sphereObjs.at(i).radius,2);
		double d= sqrt(pow(b,2)- 4*a*c);
		double t1=(-b+d)/(2*a);
		double t2=(-b-d)/(2*a);

		if((pow(b,2)- 4*a*c)<0) continue;
		double temp_T=-1;

		if(t1<0 && t2<0) continue;
		if(t1>0 && t2<0) temp_T=t1;
		if(t1<0 && t2>0) temp_T=t2;
		if(t1>0 && t2>0) temp_T=fmin(t1,t2);

		if(temp_T>0 && temp_T<tmin && temp_T<=farDist)
		{
			tmin=temp_T;
			doesIntersect=true;
			point intersection(BufferPoint.x+tmin*Rd.x,BufferPoint.y+tmin*Rd.y,BufferPoint.z+tmin*Rd.z);
			P=intersection;
			normal=sphereObjs.at(i).getNormal(P);
			pixel=sphereObjs.at(i).myColor;
			ambient=sphereObjs.at(i).ambient;
			diffuse=sphereObjs.at(i).diffuse;
			specular=sphereObjs.at(i).specular;
			reflection=sphereObjs.at(i).reflection;
			shininess=sphereObjs.at(i).shininess;
			isSphere=true;
			
		}

		
		else continue;
	}
	
	point Ro=BufferPoint;
	
	for(int i=0;i<pyramidBases.size();i++)
	{
		double Z=pyramidBases.at(i).p1.z;
		double t=(Z-Ro.z)/Rd.z;
		
		point intersection(Ro.x+t*Rd.x,Ro.y+t*Rd.y,Ro.z+t*Rd.z);
		if(intersection.x>=pyramidBases.at(i).p1.x && intersection.x<=pyramidBases.at(i).p1.x+pyramidBases.at(i).width)
		{
			if(intersection.y>=pyramidBases.at(i).p1.y && intersection.y<=pyramidBases.at(i).p1.y+pyramidBases.at(i).width)
			{
				if(t>0 && t<tmin && t<=maxdist)
				{
					doesIntersect=true;
					tmin=t;
					P=intersection;
					normal=pyramidBases.at(i).getNormal();
					pixel=pyramidBases.at(i).myColor;
					ambient=pyramidBases.at(i).ambient;
					diffuse=pyramidBases.at(i).diffuse;
					specular=pyramidBases.at(i).specular;
					reflection=pyramidBases.at(i).reflection;
					shininess=pyramidBases.at(i).shininess;
				}
			}
		}	
		
	}

	for(int i=0;i<pyramidTriangles.size();i++)
	{
		PyramidTriangle temp=pyramidTriangles.at(i);

		point a1(temp.p1.x-temp.p2.x,temp.p1.y-temp.p2.y,temp.p1.z-temp.p2.z);
		point a2(temp.p1.x-temp.p3.x,temp.p1.y-temp.p3.y,temp.p1.z-temp.p3.z);
		point a3(Rd.x,Rd.y,Rd.z);
		point a4(temp.p1.x-Ro.x,temp.p1.y-Ro.y,temp.p1.z-Ro.z);

		double beta=determinant(a4,a2,a3)/determinant(a1,a2,a3);
		double gamma=determinant(a1,a4,a3)/determinant(a1,a2,a3);
		double t=determinant(a1,a2,a4)/determinant(a1,a2,a3);

		double alpha=1-beta-gamma;
		if(beta>=0 && beta<=1 && gamma>=0 && gamma<=1 && alpha>=0 && alpha<=1 && (alpha+beta+gamma)<=1)
		{
			if(t<tmin && t>0 && t<=maxdist)
			{
				tmin=t;
				doesIntersect=true;
				point intersection(Ro.x+t*Rd.x,Ro.y+t*Rd.y,Ro.z+t*Rd.z);
				P=intersection;
				normal=temp.getNormal();
				pixel=temp.col;
				ambient=pyramidTriangles.at(i).ambient;
				diffuse=pyramidTriangles.at(i).diffuse;
				specular=pyramidTriangles.at(i).specular;
				reflection=pyramidTriangles.at(i).reflection;
				shininess=pyramidTriangles.at(i).shininess;
			}
		}
	}

	//checkerboard
	Ro=BufferPoint;
	Rd.normalize();
	double t=-Ro.z/Rd.z;
	point w(Ro.x+t*Rd.x,Ro.y+t*Rd.y,Ro.z+t*Rd.z);
	
	if(t>0 && t<tmin && t<=maxdist)
	{
		tmin=t;
		doesIntersect=true;
		point intersection(Ro.x+t*Rd.x,Ro.y+t*Rd.y,Ro.z+t*Rd.z);
		P=intersection;
		normal=checkerboard.getNormal();
		checker=true;
		ambient=checkerboard.ambient;
		diffuse=checkerboard.diffuse;
		reflection=checkerboard.reflection;
		specular=1;
		shininess=1;
		
		
		if(texture==1){
			pixel=addTexture(P);
			
		}
		else{
			int x=floor((Ro.x+Rd.x*t)/checkerboard.cellWidth);
			int y=floor((Ro.y+Rd.y*t)/checkerboard.cellWidth);
			
			if((x+y)%2==0){
				pixel.r=1;
				pixel.g=1;
				pixel.b=1;
			} 
			else{
				pixel.r=0;
				pixel.g=0;
				pixel.b=0;
			}
		}
		
	}

	pixel.r=max(pixel.r,0.0);
	pixel.g=max(pixel.g,0.0);
	pixel.b=max(pixel.b,0.0);
	
	
	if(doesIntersect==true){
		double lambert=0;
		double phong=0;
		point NNormal=normal;
	
		point V(P.x-camPos_Captured.x, P.y-camPos_Captured.y, P.z-camPos_Captured.z);
		V.normalize();
		for(int i=0;i<lightSources.size();i++)
		{
			bool res=illuminatesPoint(P,lightSources.at(i).position);
			if(res==false) continue;
			else{
				NormalLight S=lightSources.at(i);
				globalCount++;
				point toSource(S.position.x-P.x,S.position.y-P.y,S.position.z-P.z);
				
				double distance=sqrt(toSource.x*toSource.x+toSource.y*toSource.y+toSource.z*toSource.z);
				
				toSource.normalize();
				normal.normalize();
				double scaling_factor =exp(-distance*distance*S.falloff);
				
				lambert += max((toSource.dot(toSource,normal)*scaling_factor),0.0);


				double tp=V.dot(V,normal)*2;
				point new_ref(V.x-normal.x*tp,V.y-normal.y*tp,V.z-normal.z*tp);
				new_ref.normalize();
				
				toSource.normalize();
				if(checker==false)
				{
					double some=new_ref.dot(new_ref, toSource);
					if(some>0) phong += max(pow(some, shininess) * scaling_factor, 0.0);
				}
				
			}
		}

		//Spotlight
		for(int i=0;i<spotLights.size();i++)
		{
			SpotLight S=spotLights.at(i);
			bool res=illuminatesPoint(P,spotLights.at(i).position);
			if(res==false) continue;
			else{
				
				SpotLight S=spotLights.at(i);
				point toSource(S.position.x-P.x,S.position.y-P.y,S.position.z-P.z);
				double distance=sqrt(toSource.x*toSource.x+toSource.y*toSource.y+toSource.z*toSource.z);
				
				toSource.normalize();
				point v2(spotLights.at(i).lookat.x-spotLights.at(i).position.x,spotLights.at(i).lookat.y-spotLights.at(i).position.y,spotLights.at(i).lookat.z-spotLights.at(i).position.z);
				v2.normalize();
				point vv(-toSource.x,-toSource.y,-toSource.z);
				double angle = acos(vv.dot(vv,v2));
				if(angle>S.cutoff*(pi/180)) continue;

				normal.normalize();
				double scaling_factor = exp(-distance*distance*S.falloff);
				
				lambert += max((toSource.dot(toSource,normal)*scaling_factor),0.0);
				toSource.x=-toSource.x;  toSource.y=-toSource.y;  toSource.z=-toSource.z;
				double tmp=(toSource.dot(toSource,normal))*2;
				point tmp1(normal.x*tmp,normal.y*tmp,normal.z*tmp);
				point R_cap(tmp1.x-toSource.x,tmp1.y-toSource.y,tmp1.z-toSource.z);
				R_cap.normalize();
				
				if(checker==false)
				{
					phong+= max((pow(R_cap.dot(R_cap,V), shininess)*scaling_factor),0.0);
				}
				
			}
		}

		color pixel_illum(diffuse*lambert *pixel.r+specular*phong*pixel.r, diffuse*lambert*pixel.g+specular*phong*pixel.g, diffuse* lambert*pixel.b+specular*phong*pixel.b);

		pixel_illum.r=pixel.r*ambient+pixel_illum.r;
		pixel_illum.g=pixel.g*ambient+pixel_illum.g;
		pixel_illum.b=pixel.b*ambient+pixel_illum.b;

		pixel_illum.r=max(pixel_illum.r,0.0);
		pixel_illum.g=max(pixel_illum.g,0.0);
		pixel_illum.b=max(pixel_illum.b,0.0);

		
		NNormal.normalize();

		double some=2*(Rd.dot(Rd,NNormal));
		NNormal.x=NNormal.x*some;  NNormal.y=NNormal.y*some;  NNormal.z=NNormal.z*some;
		point reflected_ray(Rd.x-NNormal.x,Rd.y-NNormal.y,Rd.z-NNormal.z);
		point ray=reflected_ray;
		reflected_ray.normalize();

		point nextPoint(P.x+0.001*reflected_ray.x,P.y+0.001*reflected_ray.y,P.z+0.001*reflected_ray.z);
		
		color c=generate_pixels(nextPoint,ray,depth-1);

		pixel_illum.r=pixel_illum.r+c.r*reflection;
		pixel_illum.g=pixel_illum.g+c.g*reflection;
		pixel_illum.b=pixel_illum.b+c.b*reflection;

		return pixel_illum;
	}
	
	else return pixel;
}

void generateImage(vector<vector<color>> pixelBuffer)
{
	bitmap_image image(numPixels, numPixels);
    for (int x = 0; x < numPixels; x++) {
        for (int y = 0; y < numPixels; y++) {
            image.set_pixel(x, y, min(pixelBuffer[y][x].r,1.0)*255, min(pixelBuffer[y][x].g,1.0)*255, min(pixelBuffer[y][x].b,1.0)*255);
        }
    }
    image.save_image("out.bmp");
}

void generateRays(vector<vector<point>> pointBuffer)
{
	vector<vector<color>> pixelBuffer;
	for(int i=0;i<numPixels;i++)
	{
		vector<color> pixels;
		for(int j=0;j<numPixels;j++)
		{
			point Rd(pointBuffer.at(i).at(j).x-camPos_Captured.x, pointBuffer.at(i).at(j).y-camPos_Captured.y, pointBuffer.at(i).at(j).z-camPos_Captured.z);
			
			pixels.push_back(generate_pixels(pointBuffer.at(i).at(j),Rd,3)); 
		}
		pixelBuffer.push_back(pixels);
	}
	
	generateImage(pixelBuffer);
	cout<<"DONE :"<<endl;
}


void generatePoints()
{
	cout << std::fixed;
    cout << std::setprecision(8);
	point midPoint;

	vector<vector<point>> pointBuffer;

	midPoint.x=camPos_Captured.x+l_Captured.x*nearDist;
	midPoint.y=camPos_Captured.y+l_Captured.y*nearDist;
	midPoint.z=camPos_Captured.z+l_Captured.z*nearDist;

	scene_Y=2*nearDist*tan((pi/180.0)*(fovY/2));
	scene_X=2*nearDist*tan((pi/180.0)*((fovY*aspectRatio)/2));

	double pixelWidth=scene_X/numPixels;
	double pixelHeight=scene_Y/numPixels;
	
	

	point midUp(midPoint.x+u_Captured.x*(scene_Y/2),midPoint.y+u_Captured.y*(scene_Y/2),midPoint.z+u_Captured.z*(scene_Y/2));
	point initialP(midUp.x-r_Captured.x*(scene_X/2),midUp.y-r_Captured.y*(scene_X/2),midUp.z-r_Captured.z*(scene_X/2));

	initialP.x=initialP.x-u_Captured.x*0.5*pixelHeight+r_Captured.x*0.5*pixelWidth;
	initialP.y=initialP.y-u_Captured.y*0.5*pixelHeight+r_Captured.y*0.5*pixelWidth;
	initialP.z=initialP.z-u_Captured.z*0.5*pixelHeight+r_Captured.z*0.5*pixelWidth;

	for(int i=0;i<numPixels;i++)
	{
		vector<point> temp;
		for(int j=0;j<numPixels;j++)
		{
			point p(initialP.x+r_Captured.x*j*pixelWidth-u_Captured.x*i*pixelHeight, initialP.y+r_Captured.y*j*pixelWidth-u_Captured.y*i*pixelHeight,  initialP.z+r_Captured.z*j*pixelWidth-u_Captured.z*i*pixelHeight);		
			temp.push_back(p);
		}
		pointBuffer.push_back(temp);
		temp.clear();		
	}
	
	generateRays(pointBuffer);
}


void drawAxes()
{
	
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES);{
		glVertex3f( 500,0,0);
		glVertex3f(-500,0,0);

		glVertex3f(0,-500,0);
		glVertex3f(0, 500,0);

		glVertex3f(0,0, 500);
		glVertex3f(0,0,-500);
	}glEnd();
}

void drawSphere(int radius, point center)
{
	point points[100][100];
	int i,j;
	int stacks=20;
	int slices=30;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
				//upper hemisphere
				glVertex3f(points[i][j].x+center.x,points[i][j].y+center.y,points[i][j].z+center.z);
				glVertex3f(points[i][j+1].x+center.x,points[i][j+1].y+center.y,points[i][j+1].z+center.z);
				glVertex3f(points[i+1][j+1].x+center.x,points[i+1][j+1].y+center.y,points[i+1][j+1].z+center.z);
				glVertex3f(points[i+1][j].x+center.x,points[i+1][j].y+center.y,points[i+1][j].z+center.z);
				//lower hemisphere
				glVertex3f(points[i][j].x+center.x,points[i][j].y+center.y,-points[i][j].z+center.z);
				glVertex3f(points[i][j+1].x+center.x,points[i][j+1].y+center.y,-points[i][j+1].z+center.z);
				glVertex3f(points[i+1][j+1].x+center.x,points[i+1][j+1].y+center.y,-points[i+1][j+1].z+center.z);
				glVertex3f(points[i+1][j].x+center.x,points[i+1][j].y+center.y,-points[i+1][j].z+center.z);
			}glEnd();
		}
	}
}


void rotate_left_right(double angle){
	point temp = l;
	l = {l.x * cos(angle) + r.x * sin(angle), l.y * cos(angle) + r.y * sin(angle),
							l.z * cos(angle) + r.z * sin(angle)};
	
	r = {r.x * cos(angle) - temp.x * sin(angle), r.y * cos(angle) - temp.y * sin(angle),
													r.z * cos(angle) - temp.z * sin(angle)};
}


void rotate_up_down(double angle){
	point temp = l;
	l = {l.x * cos(angle) - u.x * sin(angle), l.y * cos(angle) - u.y * sin(angle),
							l.z * cos(angle) - u.z * sin(angle)};
	u = {u.x * cos(angle) + temp.x * sin(angle), u.y * cos(angle) + temp.y * sin(angle),
													u.z * cos(angle) + temp.z * sin(angle)};
}


void tilt(double angle){
	point temp = r;
	r = {r.x * cos(angle) + u.x * sin(angle), r.y * cos(angle) + u.y * sin(angle),
							r.z * cos(angle) + u.z * sin(angle)};
	u = {u.x * cos(angle) - temp.x * sin(angle), u.y * cos(angle) - temp.y * sin(angle),
													u.z * cos(angle) - temp.z * sin(angle)};
}


void keyboardListener(unsigned char key, int x,int y){	
	switch(key){
		case '1':
			rotate_left_right(-0.01);
			break;
		case '2':
			rotate_left_right(0.01);
			break;
		case '3':
			rotate_up_down(0.01);
			break;
		case '4':
			rotate_up_down(-0.01);
			break;
		case '5':
			tilt(0.01);
			break;
		case '6':
			tilt(-0.01);
			break;
		case '0':
			l_Captured=l;
			l_Captured.normalize();

			r_Captured=r;
			r_Captured.normalize();
			
			u_Captured=u;
			u_Captured.normalize();
			
			camPos_Captured=pos;
			
			generatePoints();
			break;
		case ' ':
			texture=1-texture;
			break;

		default:
			break;
	}
}



void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
			pos.x=pos.x-l.x*2;
			pos.y=pos.y-l.y*2;
			pos.z=pos.z-l.z*2;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
			pos.x=pos.x+l.x*2;
			pos.y=pos.y+l.y*2;
			pos.z=pos.z+l.z*2;
			break;

		case GLUT_KEY_RIGHT:
			pos.x=pos.x+r.x*2;
			pos.y=pos.y+r.y*2;
			pos.z=pos.z+r.z*2;
			//cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			pos.x=pos.x-r.x*2;
			pos.y=pos.y-r.y*2;
			pos.z=pos.z-r.z*2;
			//cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			pos.x=pos.x+u.x*2;
			pos.y=pos.y+u.y*2;
			pos.z=pos.z+u.z*2;			
			break;
		case GLUT_KEY_PAGE_DOWN:
			pos.x=pos.x-u.x*2;
			pos.y=pos.y-u.y*2;
			pos.z=pos.z-u.z*2;
			break;
		

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/

	
	checkerboard.drawCheckerBoard();
	for(int i=0;i<sphereObjs.size();i++)
	{
		sphereObjs.at(i).drawSphere();
	}

	for(int i=0;i<pyramidBases.size();i++)
	{
		pyramidBases.at(i).drawRectangle();
	}

	for(int i=0;i<pyramidTriangles.size();i++)
	{
		pyramidTriangles.at(i).drawTriangle();
	}

	glColor3f(1,1,1);
	for(int i=0;i<lightSources.size();i++)
	{
		point tmp=lightSources.at(i).position;
		drawSphere(5,tmp);
	}
	for(int i=0;i<spotLights.size();i++)
	{
		glColor3f(1,0.8,0);
		point tmp=spotLights.at(i).position;
		drawSphere(5,tmp);
	}
	
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	globalCount=0;
	GPP=0;
	

	texture=0;
	pos.x=100.0; // 70.0 100.0;
	pos.y=100;
	pos.z=100;

	u.x=0;
	u.y=0;
	u.z=1;

	l.x=-1/sqrt(2);
	l.y=-1/sqrt(2);
	l.z=0;

	r.x=-1/sqrt(2);
	r.y=1/sqrt(2);
	r.z=0;
	
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fovY, aspectRatio, nearDist,	farDist);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}


int main(int argc, char **argv){
	//take input
	ifstream input;
	
    input.open ("description.txt");

	input>>nearDist>>farDist>>fovY>>aspectRatio;
	
	input>>numRecursion>>numPixels;

	input>>cellWidth;
	CheckerBoard tmp(cellWidth);
	checkerboard=tmp;
	input>>checkerboard.ambient>>checkerboard.diffuse>>checkerboard.reflection;

	input>>numObjects;
	string objType;
	
	for(int i=0;i<numObjects;i++)
	{
		input>>objType;
		if(objType=="sphere")
		{
			point p;
			input>>p.x>>p.y>>p.z;
			double r;
			input>>r;
			color c;
			input>>c.r>>c.g>>c.b;
			Sphere tmp(p,r,c);
			input>>tmp.ambient>>tmp.diffuse>>tmp.specular>>tmp.reflection>>tmp.shininess;
			sphereObjs.push_back(tmp);
		}
		else if(objType=="pyramid")
		{
			point p1;
			double width,height;
			input>>p1.x>>p1.y>>p1.z>>width>>height;
			point p2(p1.x+width,p1.y,p1.z);
			point p3(p1.x+width,p1.y+width,p1.z);
			point p4(p1.x,p1.y+width,p1.z);

			point mid(p1.x+width/2,p1.y+width/2,p1.z);
			point peak(mid.x,mid.y,mid.z+height);

			color c;
			input>>c.r>>c.g>>c.b;
			PyramidBase pyr(p1,p2,p3,p4,c);
			pyr.width=width;
			input>>pyr.ambient>>pyr.diffuse>>pyr.specular>>pyr.reflection>>pyr.shininess;
			pyramidBases.push_back(pyr);

			PyramidTriangle A(p1,p2,peak,c);
			A.ambient=pyr.ambient; A.diffuse=pyr.diffuse; A.specular=pyr.specular; 
			A.reflection=pyr.reflection; A.shininess=pyr.shininess;
			pyramidTriangles.push_back(A);

			PyramidTriangle B(p2,p3,peak,c);
			B.ambient=pyr.ambient; B.diffuse=pyr.diffuse; B.specular=pyr.specular; 
			B.reflection=pyr.reflection; B.shininess=pyr.shininess;
			pyramidTriangles.push_back(B);

			PyramidTriangle C(p3,p4,peak,c);
			C.ambient=pyr.ambient; C.diffuse=pyr.diffuse; C.specular=pyr.specular; 
			C.reflection=pyr.reflection; C.shininess=pyr.shininess;
			pyramidTriangles.push_back(C);

			PyramidTriangle D(p4,p1,peak,c);
			D.ambient=pyr.ambient; D.diffuse=pyr.diffuse; D.specular=pyr.specular; 
			D.reflection=pyr.reflection; D.shininess=pyr.shininess;
			pyramidTriangles.push_back(D);
		}
	}
	

	int numLights;
	input>>numLights;
	for(int i=0;i<numLights;i++)
	{
		point source;
		double falloff;
		input>>source.x>>source.y>>source.z>>falloff;
		NormalLight lightSrc(source,falloff);
		lightSources.push_back(lightSrc);
	}

	input>>numLights;
	for(int i=0;i<numLights;i++)
	{
		point source,look;
		double falloff,cut;
		input>>source.x>>source.y>>source.z>>falloff>>look.x>>look.y>>look.z>>cut;
		SpotLight lightSrc(source,falloff,look,cut);
		spotLights.push_back(lightSrc);
	}

	outp.open("ooo.txt");

	Texheight = b_img.height();
	Texwidth = b_img.width();
	textureBuffer = new color* [Texwidth];
	for (int i = 0; i < Texwidth; i++) {
		textureBuffer[i] = new color [Texheight];
		for (int j = 0; j < Texheight; j++) {
			unsigned char r, g, b;
			b_img.get_pixel(i, j, r, g, b);
			color c(r/255.0, g/255.0, b/255.0);
			textureBuffer[i][j] = c;
		}
	}

	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
