//cs371 Computer Graphics - Fall 2015
//hw5.cpp
//Author: Gordon Griesel
//Date:   Summer 2013 to present
//XWindows program for graphics output
//
//This is a simple orthographic ray tracer.
//
//modified by: Daniel Stieber
//date: 10/16/2015
//purpose: Trace rays through a triange mesh of my initials (DS).
//
//
//
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/Xutil.h>
#include "log.h"
#include "ppm.h"

#define rnd() (float)rand() / (float)RAND_MAX

//XWindows globals
static Display *dpy;
static Window win;
static GC gc;
//
//Variable types...
typedef double Flt;
typedef Flt Vec[3];
typedef int iVec[3];
//
//Structures
typedef struct t_ray {
	//Ray origin and direction
	Vec o, d;
} Ray;

typedef struct t_hit {
	//t     = distance to hit point
	//p     = hit point
	//norm  = normal of surface hit
	//color = color of surface hit
	Flt t;
	Vec p, norm, color;
} Hit;

enum {
	TYPE_NONE=0,
	TYPE_RING,
	TYPE_SPHERE,
	TYPE_TRIANGLE,
	TYPE_CYLINDER,
	TYPE_CONE
};

typedef struct t_object {
	int type;
	Vec center;
	Vec norm;
    Vec tri[3];
	Flt radius;
	Flt iRadius;
    Vec color;
} Object;
#define MAXOBJECTS 10000
Object object[MAXOBJECTS];
int nobjects=0;
//
#define XRES 640.0
#define YRES 480.0
static int xres=640, yres=480;
static int done=0;
static int mode=0;
int perspective=0;
//
//Function prototypes
void initXWindows(void);
void cleanupXWindows(void);
void init(void);
void checkResize(XEvent *e);
void checkMouse(XEvent *e);
void checkKeys(XEvent *e);
void render();
//vector functions...
void vecZero(Vec v);
void vecMake(Flt a, Flt b, Flt c, Vec v);
void vecCopy(Vec source, Vec dest);
void vecSub(Vec v0, Vec v1, Vec dest);
void vecNormalize(Vec v);
void vecCrossProduct(Vec v0, Vec v1, Vec dest);
void getTriangleNormal(Vec tri[3], Vec norm);
void RayPoint(Ray *ray, Flt t, Vec h);
void scale(float x, Vec vec);
void move(Vec offset, Vec vec);
int rayTriangleIntersect(Object *o, Ray *ray, Hit *hit);
int getBarycentric2(Vec a, Vec b, Vec c, Vec p, Flt *u, Flt *v);
int pointInTriangle(Vec a, Vec b, Vec c, Vec p, Flt *u, Flt *v);
Flt getTriangleArea(Vec t0, Vec t1, Vec t2);
Flt vecDotProduct(Vec v0, Vec v1);
Flt vecLength(Vec v);

int main(void)
{
	initXWindows();
	//init();
	while(!done) {
		while(XPending(dpy)) {
			XEvent e;
			XNextEvent(dpy, &e);
			checkResize(&e);
			checkMouse(&e);
			checkKeys(&e);
			//render();
		}
	}
	cleanupXWindows();
	return 0;
}

void cleanupXWindows(void)
{
	XDestroyWindow(dpy, win);
	XCloseDisplay(dpy);
	logClose();
}

void initXWindows(void)
{
	int scr;

	logOpen();
	//Log("init_xwindows()...\n");
	if(!(dpy=XOpenDisplay(NULL))) {
		fprintf(stderr, "ERROR: could not open display\n");
		exit(1);
	}
	scr = DefaultScreen(dpy);
	win = XCreateSimpleWindow(dpy, RootWindow(dpy, scr), 1, 1,
							xres, yres, 0, 0x00ffffff, 0x00000000);
	XStoreName(dpy, win, "CS371. Press R to render");
	gc = XCreateGC(dpy, win, 0, NULL);
	XMapWindow(dpy, win);
	XSelectInput(dpy, win, ExposureMask | StructureNotifyMask |
							PointerMotionMask | ButtonPressMask |
							ButtonReleaseMask | KeyPressMask);
	srand((unsigned)time(NULL));
}


void takeScreenshot(const char *filename, int reset)
{
    static int picnum = 0;
    int x,y;
    XWindowAttributes gwa;
    XGetWindowAttributes(dpy, win, &gwa);
    int width = gwa.width;
    int height = gwa.height;
    if (reset)
        picnum = 0;
    XImage *image = XGetImage(dpy, win,
                    0, 0 , width, height, AllPlanes, ZPixmap);
    Ppmimage *ppmImage = ppm6CreateImage(width, height);
    unsigned char *ptr;
    unsigned long red_mask   = image->red_mask;
    unsigned long green_mask = image->green_mask;
    unsigned long blue_mask  = image->blue_mask;
    unsigned long pixel;
    for (x=0; x<width; x++) {
        for (y=0; y<height; y++) {
            pixel = XGetPixel(image,x,y);
            ptr = (unsigned char *)ppmImage->data + ((x + width * y) * 3);
            *(ptr+0) = (pixel & red_mask) >> 16;
            *(ptr+1) = (pixel & green_mask) >> 8;
            *(ptr+2) = (pixel & blue_mask);
        }
    }
    //save image
    //generate a random filename...
    char ts[256] = "";
    strcpy(ts, filename);
    if (ts[0] == '\0') {
        sprintf(ts,"./sshot%i.ppm", picnum);
        picnum++;
    }
    ppm6SaveImage(ppmImage, ts);
    ppm6CleanupImage(ppmImage);
    XFree(image);
}


char *findSpace(char *ptr) {
    while (*ptr && !isspace(*ptr))
        ptr++;
    return ptr;
}

void initCSUB(void)
{
    Object *o;
    nobjects=0;
    char ts[200], *ptr;
    int nverts=0, nfaces=0;
    Vec *vert = new Vec[3000];
    iVec *face = new iVec[3000];
    FILE *fin = fopen("ship.txt", "r");
    if (!fin) printf("file not open\n");
    if (fin) {
        while(1) {
            if (feof(fin)) break;
            fgets(ts,100,fin);
            ptr = ts;
            if (*ptr == 'v') {
                ptr = findSpace(ptr)+1;
                vert[nverts][0] = atof(ptr) * 200.0;
                ptr = findSpace(ptr)+1;
                vert[nverts][1] = atof(ptr) * 200.0;
                ptr = findSpace(ptr)+1;
                vert[nverts][2] = atof(ptr) * 200.0;
                nverts++;
                continue;
            }
            if (*ptr == 'f') {
                ptr = findSpace(ptr)+1;
                face[nfaces][0] = atoi(ptr);
                ptr = findSpace(ptr)+1;
                face[nfaces][1] = atoi(ptr);
                ptr = findSpace(ptr)+1;
                face[nfaces][2] = atoi(ptr);
                nfaces++;
            }
        }
        fclose(fin);
    }
    Vec mv;
    for (int i=0; i<nfaces; i++) {
        o = &object[nobjects];
        o->type = TYPE_TRIANGLE;
        int f;
        f=face[i][0]-1; vecMake(vert[f][0],vert[f][1],vert[f][2],o->tri[1]);
        f=face[i][1]-1; vecMake(vert[f][0],vert[f][1],vert[f][2],o->tri[0]);
        f=face[i][2]-1; vecMake(vert[f][0],vert[f][1],vert[f][2],o->tri[2]);
        vecMake((float)i/nfaces + .2,(float)i/nfaces + .2,(float)i/nfaces + .2, o->color);
        //vecMake(0.9,0.9,0.9, o->color);
        vecMake(-20.0, -300.0, 0.0, mv);
        for (int j=0; j<3; j++) {
            move(mv, o->tri[j]);
            scale(0.4, o->tri[j]);
        }
        getTriangleNormal(o->tri, o->norm);
        //change the color based on normal
        //if (o->norm[2] > 0.5 || o->norm[2] < -0.5) {
        //    vecMake(.1,.1,.1, o->color);
        //}
        nobjects++;
    }
    delete vert;
    delete face;
}

void init(void)
{
	Object *o;
	nobjects=0;
    int nTriangles = 24;
    Vec mv;
	vecMake(-6.5,-3.5,0, mv);

    Vec colors[10] = {
        {0.9,0.1,0.3},{0.6,0.1,0.9},
        {0.8,0.2,0.5},{0.3,0.2,0.7},
        {0.6,0.2,0.4},{0.9,0.2,0.7},
        {0.8,0.1,0.8},{0.8,0.7,1.0},
        {0.1,0.3,0.8},{0.5,0.7,0.9}
    };

    Vec triVerts[72] = {
        {0,7,0},{4,7,0},{1,6,0},
        {1,6,0},{4,7,0},{4,6,0},
        {4,6,0},{4,7,0},{5,5,0},
        {4,7,0},{6,5,0},{5,5,0},
        {5,5,0},{6,5,0},{5,2,0},
        {5,2,0},{6,5,0},{6,2,0},
        {4,1,0},{5,2,0},{6,2,0},
        {4,0,0},{4,1,0},{6,2,0},
        {1,1,0},{4,1,0},{4,0,0},
        {0,0,0},{1,1,0},{4,0,0},
        {1,1,0},{0,0,0},{0,7,0}, //11
        {1,1,0},{0,7,0},{1,6,0}, //12
        {8,0,0},{8,1,0},{10,1,0},
        {8,0,0},{10,1,0},{11,0,0},
        {10,1,0},{13,2,0},{11,0,0},
        {10,1,0},{12,3,0},{13,2,0},
        {13,2,0},{12,3,0},{13,4,0},
        {9,4,0},{13,4,0},{12,3,0},
        {8,3,0},{9,4,0},{12,3,0},
        {8,3,0},{8,5,0},{9,4,0},
        {8,5,0},{10,7,0},{9,4,0},
        {9,4,0},{10,7,0},{11,6,0},
        {10,7,0},{13,7,0},{11,6,0},
        {11,6,0},{13,7,0},{13,6,0}
    };

    for (int i=0; i<nTriangles; i++) {
        o = &object[nobjects];
	    o->type = TYPE_TRIANGLE;
        for (int j=0; j<3; j++) {
            vecCopy(triVerts[i*3+j], o->tri[j]);       
            move(mv, o->tri[j]);
            scale(30.0, o->tri[j]);
        }
    getTriangleNormal(o->tri, o->norm);
	vecMake(colors[i%10][0],colors[i%10][1],
            colors[i%10][2], o->color);
	nobjects++;
    }
}

#define VecCross(a,b,c) \
(c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; \
(c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; \
(c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]

void vecCrossProduct(Vec v0, Vec v1, Vec dest)
{
	dest[0] = v0[1]*v1[2] - v1[1]*v0[2];
	dest[1] = v0[2]*v1[0] - v1[2]*v0[0];
	dest[2] = v0[0]*v1[1] - v1[0]*v0[1];
}

void getTriangleNormal(Vec tri[3], Vec norm)
{
	Vec v0,v1;
	vecSub(tri[1], tri[0], v0);
	vecSub(tri[2], tri[0], v1);
	vecCrossProduct(v0, v1, norm);
	//VecCross(v0, v1, norm);
	vecNormalize(norm);
}

Flt vecDotProduct(Vec v0, Vec v1)
{
	return (v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]);
}

void vecZero(Vec v)
{
	v[0] = v[1] = v[2] = 0.0;
}

void vecMake(Flt a, Flt b, Flt c, Vec v)
{
	v[0] = a;
	v[1] = b;
	v[2] = c;
}

void vecCopy(Vec source, Vec dest)
{
	dest[0] = source[0];
	dest[1] = source[1];
	dest[2] = source[2];
}

void move(Vec offset, Vec vec)
{
    vec[0]+=offset[0];
    vec[1]+=offset[1];
    vec[2]+=offset[2];
}

void scale(float x, Vec vec)
{
    vec[0] = vec[0] * x;
    vec[1] = vec[1] * x;
    vec[2] = vec[2] * x;
}

Flt vecLength(Vec vec)
{
	return sqrt(vecDotProduct(vec, vec));
}

void vecNormalize(Vec v)
{
    Flt length = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if (length==0.0)
        return;
    length = 1.0 / sqrt(length);
    v[0] *= length;
    v[1] *= length;
    v[2] *= length;
}

void vecSub(Vec v0, Vec v1, Vec dest)
{
	dest[0] = v0[0] - v1[0];
	dest[1] = v0[1] - v1[1];
	dest[2] = v0[2] - v1[2];
}

Flt getTriangleArea(Vec t0, Vec t1, Vec t2)
{
	// a+b+c = perimeter
	// s = (a+b+c)/2
	// area = sqrt(s * (s-a * s-b * s-c) )
	//
	Flt a,b,c,s,x;
	Vec v0,v1,v2;
	vecSub(t1, t0, v0);
	vecSub(t2, t1, v1);
	vecSub(t2, t0, v2);
	a = vecLength(v0);
	b = vecLength(v1);
	c = vecLength(v2);
	if (a+b < c || a+c < b || c+b < a)
		return 0.0;
	s = (a+b+c) / 2.0;
	x = s * (s-a) * (s-b) * (s-c);
	return sqrt(x);
}

int getBarycentric2(Vec a, Vec b, Vec c, Vec p, Flt *u, Flt *v)
{
	//From Christer Ericson's Real-Time Collision Detection
	Vec v0,v1,v2;
	vecSub(b, a, v0);
	vecSub(c, a, v1);
	vecSub(p, a, v2);
	Flt d00 = vecDotProduct(v0, v0);
	Flt d01 = vecDotProduct(v0, v1);
	Flt d11 = vecDotProduct(v1, v1);
	Flt d20 = vecDotProduct(v2, v0);
	Flt d21 = vecDotProduct(v2, v1);
	Flt denom = d00 * d11 - d01 * d01;
	if (denom == 0.0) return 0;
	*u = (d11 * d20 - d01 * d21) / denom;
	*v = (d00 * d21 - d01 * d20) / denom;
	return (*u >= 0.0 && *v >= 0.0 && *u + *v <= 1.0);
}

int pointInTriangle(Vec a, Vec b, Vec c, Vec p, Flt *u, Flt *v)
{
	//source:
	//http://blogs.msdn.com/b/rezanour/archive/2011/08/07/
	//
	Vec cross0,cross1,cross2;
	Vec ba,ca,pa;
	vecSub(b,a,ba);
	vecSub(c,a,ca);
	vecSub(p,a,pa);
    //This is a half-space test
	vecCrossProduct(ca,pa,cross1);
	vecCrossProduct(ca,ba,cross0);
	if (vecDotProduct(cross0, cross1) < 0.0)
		return 0;
	//This is a half-space test
	vecCrossProduct(ba,pa,cross2);
	vecCrossProduct(ba,ca,cross0);
	if (vecDotProduct(cross0, cross2) < 0.0)
		return 0;
	//Point is within 2 half-spaces
	//Get area proportions
	//Area is actually length/2
	Flt areaABC = vecLength(cross0);
	Flt areaV = vecLength(cross1);
	Flt areaU = vecLength(cross2);
	*u = areaU / areaABC;
	*v = areaV / areaABC;
	return (*u >= 0.0 && *v >= 0.0 && *u + *v <= 1.0);
}


void checkResize(XEvent *e)
{
	//Respond to ConfigureNotify.
	//Event sent by the server when window is resized.
	if (e->type != ConfigureNotify)
		return;
	XConfigureEvent xce = e->xconfigure;
	xres = xce.width;
	yres = xce.height;
}

void checkMouse(XEvent *e)
{
	if (e->type == ButtonPress) {
		Log("ButtonPress\n");
	}
}

void checkKeys(XEvent *e)
{
	if (e->type == KeyPress) {
		int key = XLookupKeysym(&e->xkey, 0);
		if (key == XK_a) {
            takeScreenshot("", 0);
            return;
        }
        if (key == XK_c) {
            initCSUB();
            render();
            mode = 1;
            return;
        }
        if (key == XK_r) {
			init();
            render();
			mode = 1;
			return;
		}
		if (key == XK_p) {
			perspective=1;
			render();
			perspective=0;
			mode = 1;
			return;
		}
		if (key == XK_Escape) {
			if (mode) {
				mode = 0;
				XClearWindow(dpy, win);
				return;
			}
			done=1;
			return;
		}
	}
}

int rayPlaneIntersect(Vec center, Vec normal, Ray *ray, Hit *hit)
{
    Vec dest;
    Flt numer, denom;
    vecSub(center, ray->o, dest);
    numer = vecDotProduct(dest, normal);
    denom = vecDotProduct(ray->d, normal);
    if (denom > 1e-6) {
        hit->t = numer/denom;
        hit->p[0] = ray->o[0] + hit->t*ray->d[0];
        hit->p[1] = ray->o[1] + hit->t*ray->d[1];
        hit->p[2] = ray->o[2] + hit->t*ray->d[2];
        return 1;
    }
    return 0;
}

int rayRingIntersect(Object *o, Ray *ray, Hit *hit)
{
	if (rayPlaneIntersect(o->center, o->norm, ray, hit)) {
		Flt d0,d1,d2,dist;
		d0 = o->center[0] - hit->p[0];
        d1 = o->center[1] - hit->p[1];
        d2 = o->center[2] - hit->p[2];
        dist = d0*d0+d1*d1+d2*d2;
		if (dist <= o->radius * o->radius && 
                dist >= o->iRadius * o->iRadius) {
			return 1;
		}
	}
	return 0;
}

void RayPoint(Ray *ray, Flt t, Vec h)
{
	h[0] = ray->o[0] + ray->d[0] * t;
	h[1] = ray->o[1] + ray->d[1] * t;
	h[2] = ray->o[2] + ray->d[2] * t;
}

int rayTriangleIntersect(Object *o, Ray *ray, Hit *hit)
{
	//Does the ray intersect the triangle's plane?
	if (rayPlaneIntersect(o->tri[1], o->norm, ray, hit)) {
		//yes
		Flt u,v;
		//Flt w;
		if (pointInTriangle(o->tri[0], o->tri[1], o->tri[2], hit->p, &u, &v)) {
			//w = 1.0 - u - v;
			return 1;
		}
	}
	return 0;
}

void trace(Ray *ray, Vec rgb)
{
	//Trace a ray through the scene.
	int i;
	Hit hit, closehit;
	Object *o;
	int h = -1;
	closehit.t = 9e9;
	for (i=0; i<nobjects; i++) {
		o = &object[i];
		switch(o->type) {
			case TYPE_RING:
				if (rayRingIntersect(o, ray, &hit)) {
					if (hit.t < closehit.t) {
						closehit.t = hit.t;
						vecCopy(hit.p, closehit.p);
						vecCopy(o->color, closehit.color);
						h=i;
					}
				}
				break;
			case TYPE_SPHERE:
				break;
			case TYPE_TRIANGLE:
                if (rayTriangleIntersect(o, ray, &hit)) {
					if (hit.t < closehit.t) {
						closehit.t = hit.t;
						vecCopy(hit.p, closehit.p);
						vecCopy(o->color, closehit.color);
						vecCopy(o->norm, closehit.norm);
						h=i;
					}
				}
				break;
			case TYPE_CYLINDER:
				break;
			case TYPE_CONE:
				break;
		}
	}
	if (h < 0) {
		return;
	}
	//The ray hit an object
	//Set the color of the pixel to the color of the object
	rgb[0] = closehit.color[0];
	rgb[1] = closehit.color[1];
	rgb[2] = closehit.color[2];
	return;
}

int rgb_to_int(Vec rgb)
{
	//Convert rgb[3] into an integer
	const float range = 255.999;
	int i;
	int cref = 0;
	for (i=0; i<3; i++) {
		//Don't let a color get too bright.
		if (rgb[i] > 1.0)
			rgb[i] = 1.0;
		//Shift left 8 bits
		cref <<= 8;
		//Put next color component in low-order byte
		cref += (int)(rgb[i]*range);
	}
	return cref;
}

void drawPixel(int x, int y, Vec rgb)
{
	//Draw color rgb to screen pixel x, y
	XSetForeground(dpy, gc, rgb_to_int(rgb));
	XDrawPoint(dpy, win, gc, x, yres-y);
}

void shoot(Vec eye, Vec pixel, Ray *ray)
{
	//Make a ray from eye through a screen pixel
	vecCopy(eye, ray->o);
	vecSub(pixel, eye, ray->d);
	vecNormalize(ray->d);
}

void render()
{
	int i,j;
	//Assume our scene was setup for a 640x480 screen.
	Flt xStart = -XRES / 2.0;
	Flt yStart = -YRES / 2.0;
	Flt xStep = XRES / (Flt)xres;
	Flt yStep = YRES / (Flt)yres;
	Vec eye, pixel, rgb;
	Ray ray;
	printf("xStep: %lf %lf\n",xStep,yStep);

	//Orthographic projection.
	//Set eye in front of each pixel.
	//Make a ray straight through each screen pixel.
	pixel[1] = yStart;
	for (i=0; i<yres; i++) {
		pixel[0] = xStart;
		for (j=0; j<xres; j++) {
			vecCopy(pixel, eye);
			if (perspective) {
				//Place code here
               Vec zero = {0, 0, 0};
               vecCopy(zero, eye);
			}
			//This moves the eye back from the screen.
			eye[2] = 1000.0;
			shoot(eye, pixel, &ray);
			vecZero(rgb);
			trace(&ray, rgb);
			drawPixel(j, i, rgb);
			pixel[0] += xStep;
		}
		pixel[1] += yStep;
	}
}


