/************************************************************************
	 File:        3DUtils.cpp

	 Author:
				  Michael Gleicher, gleicher@cs.wisc.edu
	 Modifier
				  Yu-Chi Lai, yu-chi@cs.wisc.edu

	 Comment:     some useful routines for writing 3D interactive
					   programs written for CS638 -
						Michael Gleicher, November 1999
				  re-written and expanded, October 2005

						+ Routines to draw objects
							- drawCube
							- drawFloor
						+	Drop shadow code
						+	MousePole code
						+ Quick and Dirty Quaternions

	  Note:        The utilities in this file are meant to be
						an example for	CS559 students to refer to. please
						follow the course policy on
						using example code!

	 Platform:    Visio Studio.Net 2003/2005

*************************************************************************/
#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <windows.h>
#include <GL/gl.h>
#include <FL/Fl.h>
#include <GL/glu.h>

#include "3DUtils.H"

#include <vector>
using std::vector;

void drawCuboid(float x, float y, float z, float sx, float sy, float sz, float rx, float ry, float rz)
//===============================================================================
{
	float unit = 1;

	glPushMatrix();
	glTranslatef(x, y, z);
	glScalef(sx, sy, sz);
	float theta1 = -radiansToDegrees(atan2(rz, rx));
	glRotatef(theta1, 0, 1, 0);
	float theta2 = -radiansToDegrees(acos(ry));
	glRotatef(theta2, 0, 0, 1);

	glBegin(GL_QUADS);
	glNormal3f(0, 0, 1);
	glVertex3f(unit, unit, unit);
	glVertex3f(-unit, unit, unit);
	glVertex3f(-unit, -unit, unit);
	glVertex3f(unit, -unit, unit);

	glNormal3f(0, 0, -1);
	glVertex3f(unit, unit, -unit);
	glVertex3f(unit, -unit, -unit);
	glVertex3f(-unit, -unit, -unit);
	glVertex3f(-unit, unit, -unit);

	// no top - it will be the point

	glNormal3f(0, -1, 0);
	glVertex3f(unit, -unit, unit);
	glVertex3f(-unit, -unit, unit);
	glVertex3f(-unit, -unit, -unit);
	glVertex3f(unit, -unit, -unit);

	glNormal3f(1, 0, 0);
	glVertex3f(unit, unit, unit);
	glVertex3f(unit, -unit, unit);
	glVertex3f(unit, -unit, -unit);
	glVertex3f(unit, unit, -unit);

	glNormal3f(-1, 0, 0);
	glVertex3f(-unit, unit, unit);
	glVertex3f(-unit, unit, -unit);
	glVertex3f(-unit, -unit, -unit);
	glVertex3f(-unit, -unit, unit);
	glEnd();

	//glBegin(GL_TRIANGLE_FAN);
	//glNormal3f(0, 1.0f, 0);
	//glVertex3f(0, 3.0f * unit, 0);
	//glNormal3f(1.0f, 0.0f, 1.0f);
	//glVertex3f(unit, unit, unit);
	//glNormal3f(-1.0f, 0.0f, 1.0f);
	//glVertex3f(-unit, unit, unit);
	//glNormal3f(-1.0f, 0.0f, -1.0f);
	//glVertex3f(-unit, unit, -unit);
	//glNormal3f(1.0f, 0.0f, -1.0f);
	//glVertex3f(unit, unit, -unit);
	//glNormal3f(1.0f, 0.0f, 1.0f);
	//glVertex3f(unit, unit, unit);
	//glEnd();
	glPopMatrix();


}


//*************************************************************************
//
// A utility function - draw a little cube at the origin (use a 
// transform to put it in the appropriate place)
// note: we pass the size of the cube (rather than using a scale
// transform) since we want our normals to stay unit length if possible
//
// Note: 
//  1. This isn't necesarily the fastest way since I recompute each
//       vertex.
//  2. Notice that I don't keep all my polygons with the same orientations!
//===============================================================================
void drawCube(float x, float y, float z, float l)
//===============================================================================
{
	glPushMatrix();
	glTranslated(x, y, z);
	glScalef(l, l, l);
	glRotated(45, 0, 0, 1);
	glBegin(GL_QUADS);
	glNormal3d(0, 0, 1);
	glVertex3d(0.5, 0.5, 0.5);
	glVertex3d(-0.5, 0.5, 0.5);
	glVertex3d(-0.5, -0.5, 0.5);
	glVertex3d(0.5, -0.5, 0.5);

	glNormal3d(0, 0, -1);
	glVertex3d(0.5, 0.5, -0.5);
	glVertex3d(0.5, -0.5, -0.5);
	glVertex3d(-0.5, -0.5, -0.5);
	glVertex3d(-0.5, 0.5, -0.5);

	glNormal3d(0, 1, 0);
	glVertex3d(0.5, 0.5, 0.5);
	glVertex3d(0.5, 0.5, -0.5);
	glVertex3d(-0.5, 0.5, -0.5);
	glVertex3d(-0.5, 0.5, 0.5);

	glNormal3d(0, -1, 0);
	glVertex3d(0.5, -0.5, 0.5);
	glVertex3d(-0.5, -0.5, 0.5);
	glVertex3d(-0.5, -0.5, -0.5);
	glVertex3d(0.5, -0.5, -0.5);

	glNormal3d(1, 0, 0);
	glVertex3d(0.5, 0.5, 0.5);
	glVertex3d(0.5, -0.5, 0.5);
	glVertex3d(0.5, -0.5, -0.5);
	glVertex3d(0.5, 0.5, -0.5);

	glNormal3d(-1, 0, 0);
	glVertex3d(-0.5, 0.5, 0.5);
	glVertex3d(-0.5, 0.5, -0.5);
	glVertex3d(-0.5, -0.5, -0.5);
	glVertex3d(-0.5, -0.5, 0.5);
	glEnd();
	glPopMatrix();
}

//*************************************************************************
//
// the two colors for the floor for the check board
//
//*************************************************************************
float GREEN_COLOR[3] = { 0.815, 0.976, 0.545 };
//float floorColor1[3] = { .7f, .7f, .7f }; // Light color
float floorColor1[3] = { 0.815, 0.976, 0.545 }; // Green color
//float floorColor1[3] = { 0.396, 0.11, 0.467 }; // Purple color


//float floorColor2[3] = { .3f, .3f, .3f }; // Dark color
float floorColor2[3] = { 0.815, 0.976, 0.545 }; // Dark color

float* vertices = 0;
float* colors = 0;
int* indices = 0;
float prevGroundScale = 0;

int getVerticesCount(int width, int height) {
	return width * height * 3;
}

int getIndicesCount(int width, int height) {
	return (width * height) + (width - 1) * (height - 2);
}

float* unscaleVertices(int width, int height, double scaleXZ, double scaleY) {
	int i = 0;

	for (int row = 0; row < height; row++) {
		for (int col = 0; col < width; col++) {
			vertices[i++] /= scaleXZ;
			vertices[i++] /= scaleY;
			vertices[i++] /= scaleXZ;
		}
	}

	return vertices;
}
float* getVertices(int scaledWidth, int scaledHeight, double scaleXZ, double scaleY, int waveLength) {

	/*if (vertices && scaleXZ == prevGroundScale) {
		return vertices;
	}*/
	prevGroundScale = scaleXZ;
	delete[] indices;
	indices = 0;

	const int width = scaledWidth / scaleXZ;
	const int height = scaledHeight / scaleXZ;
	const int translateOffsetW = scaledWidth / 2;
	const int translateOffsetH = scaledHeight / 2;
	const int floorOffset = 10;
	delete[] vertices;
	delete[] colors;
	vertices = new float[getVerticesCount(scaledWidth, scaledHeight)];
	colors = new float[getVerticesCount(scaledWidth, scaledHeight)];
	int i = 0;


	for (int row = 0; row < scaledHeight; row++) {
		for (int col = 0; col < scaledWidth; col++) {
			const float realRow = row / scaleXZ;
			const float realCol = col / scaleXZ;
			const float unscaledA = 1.0 / 2 * sin(realRow / waveLength) + 1.0 / 2 * sin(realCol / waveLength);
			const float normalizedA = (unscaledA + 1) / 2; // 0 ~ 1
			const float max = 0.9;
			const float min = 0.4;
			const float a = (normalizedA) * (max - min) + min; // 0 ~ 1
			//const float a = 1-(normalizedA) * (max - min) + min; // 0 ~ 1
			colors[i] = 0.3;
			//colors[i] = a * floorColor1[0];
			vertices[i++] = (float)col - translateOffsetH;
			//colors[i] = a * floorColor1[1];
			colors[i] = min(1, a);
			vertices[i++] = 5 * sin(realRow / waveLength) + 5 * sin(realCol / waveLength) - floorOffset;
			//colors[i] = a * floorColor1[2];
			colors[i] = 0.3;
			vertices[i++] = (float)row - translateOffsetW;
		}
	}

	unscaleVertices(scaledWidth, scaledHeight, scaleXZ, scaleY);

	return vertices;
}

int* getIndices(int width, int height) {
	//if (indices) return indices;
	const int iSize = getIndicesCount(width, height);

	indices = new int[iSize];
	int i = 0;

	for (int row = 0; row < height - 1; row++) {
		if ((row & 1) == 0) { // even rows
			for (int col = 0; col < width; col++) {
				indices[i++] = col + row * width;
				indices[i++] = col + (row + 1) * width;
			}
		}
		else { // odd rows
			for (int col = width - 1; col > 0; col--) {
				indices[i++] = col + (row + 1) * width;
				indices[i++] = col - 1 + +row * width;
			}
		}
	}
	/*if ((mHeight & 1) && mHeight > 2) {
		mpIndices[i++] = (mHeight - 1) * mWidth;
	}*/

	return indices;
}

void renderTerrain(float size, float scale, int waveLength) {
	glColor4fv(floorColor1);
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	const float scaledSize = size * scale;
	glVertexPointer(3, GL_FLOAT, 0, getVertices(scaledSize, scaledSize, scale, 1, waveLength));
	glColorPointer(3, GL_FLOAT, 0, colors);

	glDrawElements(GL_TRIANGLE_STRIP, getIndicesCount(scaledSize, scaledSize), GL_UNSIGNED_INT, getIndices(scaledSize, scaledSize));
	glDisableClientState(GL_VERTEX_ARRAY);
}

//*************************************************************************
//
// Draw the check board floor without texturing it
//===============================================================================
void drawFloor(float size, int scale, int waveLength)
//===============================================================================
{
	renderTerrain(size, scale, waveLength);
	return;

}

void drawGridFloor(float size, int nSquares) {
	// parameters:
	float maxX = size / 2, maxY = size / 2;
	float minX = -size / 2, minY = -size / 2;

	int x, y, v[3], i;
	float xp, yp, xd, yd;
	v[2] = 0;
	xd = (maxX - minX) / ((float)nSquares);
	yd = (maxY - minY) / ((float)nSquares);
	glBegin(GL_QUADS);
	for (x = 0, xp = minX; x < nSquares; x++, xp += xd) {
		for (y = 0, yp = minY, i = x; y < nSquares; y++, i++, yp += yd) {
			glColor4fv(i % 2 == 1 ? floorColor1 : floorColor2);
			glNormal3f(0, 1, 0);
			glVertex3d(xp, 0, yp);
			glVertex3d(xp, 0, yp + yd);
			glVertex3d(xp + xd, 0, yp + yd);
			glVertex3d(xp + xd, 0, yp);

		} // end of for j
	}// end of for i
	glEnd();
}

//*************************************************************************
//
// Shadows
// * Before drawing the floor, setup the stencil buffer and enable depth testing
// * A trick: because we only want to draw the shadows
//            on the ground plane, when we draw the ground plane
//            we'll set the stencil buffer to 1. then we'll only
//            draw shadows where the stencil buffer is one, so we 
//            don't have shadows floating in space
//===============================================================================
void setupFloor(void)
//===============================================================================
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_ALWAYS, 0x1, 0x1);
	glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
	glStencilMask(0x1);		// only deal with the 1st bit
}

//*************************************************************************
//
// * now draw the objects normally - be sure to have them write the stencil
//   to show the floor isn't there anymore
//===============================================================================
void setupObjects(void)
//===============================================================================
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_ALWAYS, 0x0, 0x0);
	glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
	glStencilMask(0x1);		// only deal with the 1st bit
}


//*************************************************************************
//
// These are cheap "hack" shadows - basically just squish the objects onto the floor.
//
// * We'll use a projection matrix to do the squishing (basically set the 
//   Y component to zero)
// * We'll also need to turn off lighting (since we want the objects to be black)
// * To make things look nice, we'll draw the shadows as transparent (so
//   we can see the groundplane)
// * To avoid Z-Fighting (since the shadows are on the ground), we'll turn
//   off the Z-Buffer - although, this means we'll see the shadows through
//   the floor
// * Finally, we'll use the stencil buffer to only draw where the ground
//   has already been drawn - this way we won't have shadows floating in
//   space!
//===============================================================================
void setupShadows(void)
//===============================================================================
{
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_STENCIL_TEST);
	glStencilFunc(GL_EQUAL, 0x1, 0x1);
	glStencilOp(GL_KEEP, GL_ZERO, GL_ZERO);
	glStencilMask(0x1);		// only deal with the 1st bit

	glPushMatrix();
	// a matrix that squishes things onto the floor
	float sm[16] = { 1,0,0,0, 0,0,0,0, 0,0,1,0, 0,0,0,1 };
	glMultMatrixf(sm);
	// draw in transparent black (to dim the floor)
	glColor4f(0, 0, 0, .5);
}

//*************************************************************************
//
// * Warning - this puts things back to a "normal" state, not
//   necessarily where they were before setupShadows
//===============================================================================
void unsetupShadows(void)
//===============================================================================
{
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_STENCIL_TEST);
	glDisable(GL_BLEND);
}


//*************************************************************************
//
// * Convert from mouse coordinates to world coordinates
//   in a sense, this mimics the old fashioned gl's mapw
//   the only advantage is that we don't have to pick the near clipping
//   plane, so we can be a little more well-balanced
//   this code mimics page 147 of the OpenGL book
//===============================================================================
int getMouseLine(double& x1, double& y1, double& z1,
	double& x2, double& y2, double& z2)
	//===============================================================================
{
	int x = Fl::event_x();
	int iy = Fl::event_y();

	double mat1[16], mat2[16];		// we have to deal with the projection matrices
	int viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, mat1);
	glGetDoublev(GL_PROJECTION_MATRIX, mat2);

	int y = viewport[3] - iy; // originally had an extra -1?

	int i1 = gluUnProject((double)x, (double)y, .25, mat1, mat2, viewport, &x1, &y1, &z1);
	int i2 = gluUnProject((double)x, (double)y, .75, mat1, mat2, viewport, &x2, &y2, &z2);

	return i1 && i2;
}


//*************************************************************************
//
// More of an explanation in the header file.
// The beginings of a mousepole handler, given a pair of points on the
// mouse ray, and the point in question, find a new point. By default,
// it assumes the plane parallel to the ground, unless the elevator
// button is pushed, in which case it is the YZ plane
//===============================================================================
inline double ABS(double x) { return (x < 0) ? -x : x; };
//===============================================================================

//*************************************************************************
//
// * When you have a mouse line, you want to pick a point where you think the user 
//   wants to drag to
// * Idea: given a plane parallel to the floor, pick a point on that
//         plane (where the line intersects it)
// * Problem: what to do when the plane is parallel to the line?
// * Problem: how to make things go vertically?
// * Answer:
//   1. Have an "elevator mode" where we use a plane perpindicular to the floor
//   2. r1,r2 are two points on the line
//      a. l is the initial position of the point - we need this to know where
//         the plane is. 
//      b. r is the resulting position. it will share 1 of its coordinates
//         with l, but will be on the line
//===============================================================================
void mousePoleGo(double r1x, double r1y, double r1z,
	double r2x, double r2y, double r2z,
	double lx, double ly, double lz,
	double& rx, double& ry, double& rz,
	bool elevator)
	//===============================================================================
{
	rx = lx; ry = ly; rz = lz;
	if (elevator || (ABS(r1y - r2y) < .01)) {
		if (ABS(r1z - r2z) > ABS(r1x - r2x)) {
			double zd = r1z - r2z;
			if (ABS(zd) > .01) {
				double zp = (lz - r1z) / zd;
				rx = r1x + (r1x - r2x) * zp;
				ry = r1y + (r1y - r2y) * zp;
			}
		}
		else {
			double xd = r1x - r2x;
			if (ABS(xd) > .01) {
				double xp = (lx - r1x) / xd;
				rz = r1z + (r1z - r2z) * xp;
				ry = r1y + (r1y - r2y) * xp;
			}
		}
	}
	else {
		double yd = r1y - r2y;
		// we have already made sure that the elevator is not singular
		double yp = (ly - r1y) / yd;
		rx = r1x + (r1x - r2x) * yp;
		rz = r1z + (r1z - r2z) * yp;
	}
}

//*******************************************************************************
//
// *
//
//*******************************************************************************
static const float rtdf = static_cast<float>(180.0 / M_PI);

//*******************************************************************************
//
// *
//===============================================================================
float radiansToDegrees(const float x)
//===============================================================================
{
	return rtdf * x;
}

//*******************************************************************************
//
// * save the lighting state on the stack so it can be restored
//
//*******************************************************************************
struct LightState {
	bool lighting;
	bool smooth;

	LightState(bool l, bool s) : lighting(l), smooth(s) {};
};
static vector<LightState> lightStateStack;

//*******************************************************************************
//
// *
//===============================================================================
void setLighting(const LightOnOff lighting, const LightOnOff smoothi)
//===============================================================================
{
	int smooth;
	bool lights = glIsEnabled(GL_LIGHTING) > 0;
	glGetIntegerv(GL_SHADE_MODEL, &smooth);
	lightStateStack.push_back(LightState(lights, (smooth == GL_SMOOTH)));

	if (lighting != keep) {
		if (lighting == on) glEnable(GL_LIGHTING);
		else glDisable(GL_LIGHTING);
	}
	if (smoothi != keep) {
		if (smoothi == on) glShadeModel(GL_SMOOTH);
		else glShadeModel(GL_FLAT);
	}

}

//*******************************************************************************
//
// *
//===============================================================================
void restoreLighting()
//===============================================================================
{
	if (!lightStateStack.empty()) {
		if ((lightStateStack.end() - 1)->lighting)
			glEnable(GL_LIGHTING);
		else
			glDisable(GL_LIGHTING);

		if ((lightStateStack.end() - 1)->smooth)
			glShadeModel(GL_SMOOTH);
		else
			glShadeModel(GL_FLAT);

	}
}