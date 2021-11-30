/************************************************************************
	 File:        TrainView.cpp

	 Author:
				  Michael Gleicher, gleicher@cs.wisc.edu

	 Modifier
				  Yu-Chi Lai, yu-chi@cs.wisc.edu

	 Comment:
						The TrainView is the window that actually shows the
						train. Its a
						GL display canvas (Fl_Gl_Window).  It is held within
						a TrainWindow
						that is the outer window with all the widgets.
						The TrainView needs
						to be aware of the window - since it might need to
						check the widgets to see how to draw

	  Note:        we need to have pointers to this, but maybe not know
						about it (beware circular references)

	 Platform:    Visio Studio.Net 2003/2005

*************************************************************************/

#include <iostream>
#include <Fl/fl.h>
#include <vector>
#include <algorithm>
#include <tuple>

// we will need OpenGL, and OpenGL needs windows.h
#include <windows.h>
//#include "GL/gl.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include "GL/glu.h"

#include "TrainView.H"
#include "TrainWindow.H"
#include "Utilities/3DUtils.H"


#ifdef EXAMPLE_SOLUTION
#	include "TrainExample/TrainExample.H"
#endif


//************************************************************************
//
// * Constructor to set up the GL window
//========================================================================
TrainView::
TrainView(int x, int y, int w, int h, const char* l)
	: Fl_Gl_Window(x, y, w, h, l)
	//========================================================================
{
	mode(FL_RGB | FL_ALPHA | FL_DOUBLE | FL_STENCIL);

	resetArcball();
}

//************************************************************************
//
// * Reset the camera to look at the world
//========================================================================
void TrainView::
resetArcball()
//========================================================================
{
	// Set up the camera to look at the world
	// these parameters might seem magical, and they kindof are
	// a little trial and error goes a long way
	arcball.setup(this, 40, 250, .2f, .4f, 0);
}

//************************************************************************
//
// * FlTk Event handler for the window
//########################################################################
// TODO: 
//       if you want to make the train respond to other events 
//       (like key presses), you might want to hack this.
//########################################################################
//========================================================================
int TrainView::handle(int event)
{
	// see if the ArcBall will handle the event - if it does, 
	// then we're done
	// note: the arcball only gets the event if we're in world view
	if (tw->worldCam->value())
		if (arcball.handle(event))
			return 1;

	// remember what button was used
	static int last_push;

	switch (event) {
		// Mouse button being pushed event
	case FL_PUSH:
		last_push = Fl::event_button();
		// if the left button be pushed is left mouse button
		if (last_push == FL_LEFT_MOUSE) {
			doPick();
			damage(1);
			return 1;
		};
		break;

		// Mouse button release event
	case FL_RELEASE: // button release
		damage(1);
		last_push = 0;
		return 1;

		// Mouse button drag event
	case FL_DRAG:

		// Compute the new control point position
		if ((last_push == FL_LEFT_MOUSE) && (selectedCube >= 0)) {
			ControlPoint* cp = &m_pTrack->points[selectedCube];

			double r1x, r1y, r1z, r2x, r2y, r2z;
			getMouseLine(r1x, r1y, r1z, r2x, r2y, r2z);

			double rx, ry, rz;
			mousePoleGo(r1x, r1y, r1z, r2x, r2y, r2z,
				static_cast<double>(cp->pos.x),
				static_cast<double>(cp->pos.y),
				static_cast<double>(cp->pos.z),
				rx, ry, rz,
				(Fl::event_state() & FL_CTRL) != 0);

			cp->pos.x = (float)rx;
			cp->pos.y = (float)ry;
			cp->pos.z = (float)rz;
			damage(1);
		}
		break;

		// in order to get keyboard events, we need to accept focus
	case FL_FOCUS:
		return 1;

		// every time the mouse enters this window, aggressively take focus
	case FL_ENTER:
		focus(this);
		break;

	case FL_KEYBOARD:
		int k = Fl::event_key();
		int ks = Fl::event_state();
		if (k == 'p') {
			// Print out the selected control point information
			if (selectedCube >= 0)
				printf("Selected(%d) (%g %g %g) (%g %g %g)\n",
					selectedCube,
					m_pTrack->points[selectedCube].pos.x,
					m_pTrack->points[selectedCube].pos.y,
					m_pTrack->points[selectedCube].pos.z,
					m_pTrack->points[selectedCube].orient.x,
					m_pTrack->points[selectedCube].orient.y,
					m_pTrack->points[selectedCube].orient.z);
			else
				printf("Nothing Selected\n");

			return 1;
		};
		break;
	}

	return Fl_Gl_Window::handle(event);
}

//************************************************************************
//
// * this is the code that actually draws the window
//   it puts a lot of the work into other routines to simplify things
//========================================================================
void TrainView::draw()
{

	//*********************************************************************
	//
	// * Set up basic opengl informaiton
	//
	//**********************************************************************
	//initialized glad
	if (gladLoadGL())
	{
		//initiailize VAO, VBO, Shader...
	}
	else
		throw std::runtime_error("Could not initialize GLAD!");

	// Set up the view port
	glViewport(0, 0, w(), h());

	// clear the window, be sure to clear the Z-Buffer too
	glClearColor(0, 0, .3f, 0);		// background should be blue

	// we need to clear out the stencil buffer since we'll use
	// it for shadows
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glEnable(GL_DEPTH);

	// Blayne prefers GL_DIFFUSE
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

	// prepare for projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	setProjection();		// put the code to set up matrices here

	//######################################################################
	// TODO: 
	// you might want to set the lighting up differently. if you do, 
	// we need to set up the lights AFTER setting up the projection
	//######################################################################
	// enable the lighting
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// top view only needs one light
	if (tw->topCam->value()) {
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHT2);
	}
	else {
		glEnable(GL_LIGHT1);
		glEnable(GL_LIGHT2);
	}

	//*********************************************************************
	//
	// * set the light parameters
	//
	//**********************************************************************
	GLfloat lightPosition1[] = { 0,1,1,0 }; // {50, 200.0, 50, 1.0};
	GLfloat lightPosition2[] = { 1, 0, 0, 0 };
	GLfloat lightPosition3[] = { 0, -1, 0, 0 };
	GLfloat yellowLight[] = { 0.5f, 0.5f, .1f, 1.0 };
	GLfloat whiteLight[] = { 1.0f, 1.0f, 1.0f, 1.0 };
	GLfloat blueLight[] = { .1f,.1f,.3f,1.0 };
	GLfloat grayLight[] = { .3f, .3f, .3f, 1.0 };

	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteLight);
	glLightfv(GL_LIGHT0, GL_AMBIENT, grayLight);

	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, yellowLight);

	glLightfv(GL_LIGHT2, GL_POSITION, lightPosition3);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, blueLight);



	//*********************************************************************
	// now draw the ground plane
	//*********************************************************************
	// set to opengl fixed pipeline(use opengl 1.x draw function)
	glUseProgram(0);

	setupFloor();
	glDisable(GL_LIGHTING);
	drawFloor(200, 10);


	//*********************************************************************
	// now draw the object and we need to do it twice
	// once for real, and then once for shadows
	//*********************************************************************
	glEnable(GL_LIGHTING);
	setupObjects();

	drawStuff();

	// this time drawing is for shadows (except for top view)
	if (!tw->topCam->value()) {
		setupShadows();
		drawStuff(true);
		unsetupShadows();
	}
}

//************************************************************************
//
// * This sets up both the Projection and the ModelView matrices
//   HOWEVER: it doesn't clear the projection first (the caller handles
//   that) - its important for picking
//========================================================================
void TrainView::
setProjection()
//========================================================================
{
	// Compute the aspect ratio (we'll need it)
	float aspect = static_cast<float>(w()) / static_cast<float>(h());

	// Check whether we use the world camp
	if (tw->worldCam->value())
		arcball.setProjection(false);
	// Or we use the top cam
	else if (tw->topCam->value()) {
		float wi, he;
		if (aspect >= 1) {
			wi = 110;
			he = wi / aspect;
		}
		else {
			he = 110;
			wi = he * aspect;
		}

		// Set up the top camera drop mode to be orthogonal and set
		// up proper projection matrix
		glMatrixMode(GL_PROJECTION);
		glOrtho(-wi, wi, -he, he, 200, -200);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(-90, 1, 0, 0);
	}
	// Or do the train view or other view here
	//####################################################################
	// TODO: 
	// put code for train view projection here!	
	//####################################################################
	else {
#ifdef EXAMPLE_SOLUTION
		trainCamView(this, aspect);
#endif
	}
}



vector<double> normalized(vector<double> arr) {
	double max = *std::max_element(arr.begin(), arr.end());
	double min = *std::min_element(arr.begin(), arr.end());
	std::transform(arr.begin(), arr.end(), arr.begin(), [max, min](double num) {
		return  (num - min) / (max - min);
		});
	return arr;
}
// t: 0~1
// w: 0~w
vector<vector<double>> calCubicBSplineN(
	int n,
	int k,
	vector<double> pt,
	// Total number of points
	int w
	// scaleOfValue = 100,
) {
	// A knot vector (t0 , t1 , ... , tn+k )
	const int sizeOfT = n + k + 1;
	vector<vector<double>>  N(w, std::vector<double>(sizeOfT, 0.0));

	// t: 0~1
	vector<double> t = normalized(pt);

	const double dtPerW = (t[sizeOfT - 1] - t[0]) / w;
	int i1 = 0;
	double _t = t[0];
	for (int _w = 0; _w < w; _w++) {
		while (_t >= t[i1]) i1++;
		const int i = i1 - 1;
		N[_w][i] = 1;
		// k bottom-up
		for (int _k = 2; _k <= k; _k++) {
			//  basis functions calculation
			int jb = i - _k + 1;
			if (jb < 0) jb = 0;
			for (int ii = jb; ii <= i; ii++)
			{

				if (ii + _k >= sizeOfT) {
					continue;
				}
				// ºu°Ê
				N[_w][ii] =
					(N[_w][ii] * (_t - t[ii])) / (t[ii + _k - 1] - t[ii]) +
					(N[_w][ii + 1] * (t[ii + _k] - _t)) / (t[ii + _k] - t[ii + 1]);
			}
		}
		_t += dtPerW;
	}

	return N;
};

std::tuple <vector<double>, vector<double>> calCubicBSplineP(
	int n,
	int k,
	vector<double> ptX,
	vector<double> ptY,
	vector<double> pt,
	// Total number of points
	int w,
	// scaleOfValue = 100,
	vector<vector<double>>  N
) {
	// A knot vector (t0 , t1 , ... , tn+k )
	vector<double> Px = ptX;
	vector<double> Py = ptY;
	vector<double> ti = normalized(pt);

	const int startW = floor((ti[k - 1] - ti[0]) / w) + 1;
	int currentSplineStartW = startW;

	vector<double> pointXsInW(w);
	vector<double> pointYsInW(w);

	double sX = 0;
	double sY = 0;
	for (int i = 0; i < n + 1; i++) {
		sX += Px[i] * N[currentSplineStartW][i];
		sY += Py[i] * N[currentSplineStartW][i];
	}
	pointXsInW[currentSplineStartW] = sX;
	pointYsInW[currentSplineStartW] = sY;

	// for spline
	for (int curSpline = k - 1; curSpline < n + 1; curSpline++) {
		const double currentSplineEndT = ti[curSpline + 1] - ti[0];
		const int currentEndW = floor(currentSplineEndT * w);

		// draw points in spline
		for (int _w = currentSplineStartW; _w < currentEndW; _w++) {
			sX = sY = 0;
			for (int i = 0; i < n + 1; i++) {
				sX += Px[i] * N[_w][i];
				sY += Py[i] * N[_w][i];
			}
			pointXsInW[_w] = sX;
			pointYsInW[_w] = sY;
		}
		currentSplineStartW = currentEndW;
	}

	return make_tuple(pointXsInW, pointYsInW);
}


void drawSpline(
	int n,
	int k,
	vector<double> ptX,
	vector<double> ptY,
	vector<double> pt,
	// Total number of points
	int w,
	int h,
	int totalNumOfPoints,
	// scaleOfValue = 100,
	vector<double>  pointXsInW,
	vector<double> pointYsInW
) {
	// A knot vector (t0 , t1 , ... , tn+k )
	vector<double> Px = ptX;
	vector<double> Py = ptY;
	vector<double> ti = normalized(pt);
	const int h1 = h - 1;


	// Draw linear control point lines
	//glBegin(GL_LINE_STRIP);
	//for (int i = 1; i < n + 1; i++) {

	//	glVertex3f(Px[i] * w, 1, Py[i] * h);

	//}
	//glEnd();




	// Draw splines
	glLineWidth(5);
	const int startW = floor((ti[k - 1] - ti[0]) * totalNumOfPoints) + 1;
	int currentSplineStartW = startW;
	// for spline
	for (int curSpline = k - 1; curSpline < n + 1; curSpline++) {
		glBegin(GL_LINE_STRIP);
		const double currentEndT = (ti[curSpline + 1] - ti[0]);

		const int currentSplineEndW = floor(currentEndT * totalNumOfPoints);
		std::cout << "currentSplineStartW " << currentSplineStartW << " " << currentSplineEndW << std::endl;
		glColor3f(double(curSpline % 10) / 10, 1 - double(curSpline % 10) / 10, double(curSpline % 4) / 5);
		/*	ctx.strokeStyle = iColor[(curSpline - k + 1) % 7];
			ctx.beginPath();
			ctx.moveTo(
				pointXsInW[currentSplineStartW] * w,
				h - pointYsInW[currentSplineStartW] * h,
				);*/
				// draw points in spline
		for (int _w = currentSplineStartW; _w < currentSplineEndW; _w++) {
			glVertex3f(pointXsInW[_w] * w, 1, pointYsInW[_w] * h);
		}

		//ctx.stroke();
		currentSplineStartW = currentSplineEndW;
		glEnd();
	}

	// draw control points
	/*ctx.lineWidth = d;
	ctx.strokeStyle = '#0000f0';
	for (let i = 0; i < n + 1; i++)
		ctx.strokeRect(Px[i] * w - d, h - Py[i] * h - d, d2, d2);*/

	for (int i = 0; i < Px.size(); i++) {
		drawCube(Px[i], 5, Py[i], 5);

	}
}

void TrainView::drawTrack(bool doingShadows)
{

	vector<ControlPoint> controlPoints = m_pTrack->points;

	if (false) {

		// Create and draw the curves for every 4 points
		// Use GL_LINE_STRIP instead to fill the gaps
		glBegin(GL_LINE_LOOP);
		// Draw the curve
		for (int i = 0; i < controlPoints.size(); i++) {
			ControlPoint point = controlPoints[i];
			glVertex3f(point.pos.x, point.pos.y, point.pos.z);
		}
		glEnd();
	}


	if (controlPoints.size() < 4) {
		return;
	}
	int n = controlPoints.size() + 2;
	int k = 4;
	/*vector<double> ptX = { 0.2, 0.8, 0.8, 0.2, 0.8, 0.8 };
	vector<double> ptY = { 0.9, 0.9, 0.2, 0.9, 0.9, 0.2 };
	vector<double> pt = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };*/
	vector<double> ptX;
	vector<double> ptY;
	for (int i = 0; i < controlPoints.size(); i++)
	{
		ControlPoint point = controlPoints[i];
		ptX.push_back(point.pos.x);
		ptY.push_back(point.pos.z);
	}
	// To make a C k-2 continuous closed loop you need only, that the last k - 1 control points repeat the first k - 1 ones, i.e. [P0 , P1 , P2 , P0 , P1 , P2 ] for n = 5, k = 4 (in this example the last 3 points are displaced a bit to make them visible).
	for (int i = 0; i < k - 1; i++)
	{
		ControlPoint point = controlPoints[i];
		ptX.push_back(point.pos.x);
		ptY.push_back(point.pos.z);

	}

	vector<double> pt;
	int sizeOfT = n + k + 1;
	for (int i = 0; i < sizeOfT; i++)
	{
		pt.push_back(i);

	}
	int w = 1;
	int h = 1;
	int totalNumOfPoints = 400;
	auto N = calCubicBSplineN(n, k, pt, totalNumOfPoints);
	/*for (int _w = 0; _w < 10; _w++)
	{
		for (int i = 0; i < n + k + 1; i++)
		{
			std::cout << N[_w][i] << " ";

		}
		std::cout << std::endl;
	}
	std::cout << std::endl;*/
	vector<double> pointXsInW, pointYsInW;
	tie(pointXsInW, pointYsInW) = calCubicBSplineP(n, k, ptX, ptY, pt, totalNumOfPoints, N);
	/*for (int i = 0; i < pointXsInW.size(); i++)
	{
		std::cout << pointXsInW[i] << " ";
	}
	std::cout << std::endl;*/

	// drawFun(N);
	drawSpline(n, k, ptX, ptY, pt, w, h, totalNumOfPoints, pointXsInW, pointYsInW);
	for (int i = 0; i < ptX.size(); i++)
	{
		std::cout << ptX[i] << " " << ptY[i] << std::endl;
	}
	std::cout << std::endl;

	vector<double> ti = normalized(pt);

	double trainT = m_pTrack->trainU;
	double numOfTrackPoints = double((n + 1) - (k - 1)) * totalNumOfPoints / (sizeOfT - 1);
	const int startW = floor((ti[k - 1] - ti[0]) * totalNumOfPoints);
	int trainW = startW + (trainT * numOfTrackPoints);
	//std::cout << "numOfTrackPoints" << numOfTrackPoints << std::endl;
	//std::cout << "trainW" << trainW << std::endl;
	double trainX = pointXsInW[trainW % totalNumOfPoints];
	double trainY = pointYsInW[trainW % totalNumOfPoints];
	drawCube(trainX, +5, trainY, 5);
	std::cout << m_pTrack->trainU << std::endl;

}

void TrainView::drawTrain(bool doingShadows)
{
	vector<ControlPoint> controlPoints = m_pTrack->points;
	drawCube(controlPoints[0].pos.x, controlPoints[0].pos.y + 5, controlPoints[0].pos.z, 5);
	m_pTrack->trainU;
	//glBegin(GL_LINE_LOOP);
	//// Draw the curve
	//for (int i = 0; i < controlPoints.size(); i++) {
	//	ControlPoint point = controlPoints[i];
	//	glVertex3f(point.pos.x, point.pos.y, point.pos.z);
	//}
	//glEnd();

}

//************************************************************************
//
// * this draws all of the stuff in the world
//
//	NOTE: if you're drawing shadows, DO NOT set colors (otherwise, you get 
//       colored shadows). this gets called twice per draw 
//       -- once for the objects, once for the shadows
//########################################################################
// TODO: 
// if you have other objects in the world, make sure to draw them
//########################################################################
//========================================================================
void TrainView::drawStuff(bool doingShadows)
{
	// Draw the control points
	// don't draw the control points if you're driving 
	// (otherwise you get sea-sick as you drive through them)
	if (!tw->trainCam->value()) {
		for (size_t i = 0; i < m_pTrack->points.size(); ++i) {
			if (!doingShadows) {
				if (((int)i) != selectedCube)
					glColor3ub(240, 60, 60);
				else
					glColor3ub(240, 240, 30);
			}
			m_pTrack->points[i].draw();
		}
	}
	// draw the track
	//####################################################################
	// TODO: 
	// call your own track drawing code
	//####################################################################

	this->drawTrack(doingShadows);

	// draw the train
	//####################################################################
	// TODO: 
	//	call your own train drawing code
	//####################################################################
	// don't draw the train if you're looking out the front window
	if (!tw->trainCam->value()) {
		this->drawTrain(doingShadows);

	}
}


// 
//************************************************************************
//
// * this tries to see which control point is under the mouse
//	  (for when the mouse is clicked)
//		it uses OpenGL picking - which is always a trick
//########################################################################
// TODO: 
//		if you want to pick things other than control points, or you
//		changed how control points are drawn, you might need to change this
//########################################################################
//========================================================================
void TrainView::
doPick()
//========================================================================
{
	// since we'll need to do some GL stuff so we make this window as 
	// active window
	make_current();

	// where is the mouse?
	int mx = Fl::event_x();
	int my = Fl::event_y();

	// get the viewport - most reliable way to turn mouse coords into GL coords
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// Set up the pick matrix on the stack - remember, FlTk is
	// upside down!
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPickMatrix((double)mx, (double)(viewport[3] - my),
		5, 5, viewport);

	// now set up the projection
	setProjection();

	// now draw the objects - but really only see what we hit
	GLuint buf[100];
	glSelectBuffer(100, buf);
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	// draw the cubes, loading the names as we go
	for (size_t i = 0; i < m_pTrack->points.size(); ++i) {
		glLoadName((GLuint)(i + 1));
		m_pTrack->points[i].draw();
	}

	// go back to drawing mode, and see how picking did
	int hits = glRenderMode(GL_RENDER);
	if (hits) {
		// warning; this just grabs the first object hit - if there
		// are multiple objects, you really want to pick the closest
		// one - see the OpenGL manual 
		// remember: we load names that are one more than the index
		selectedCube = buf[3] - 1;
	}
	else // nothing hit, nothing selected
		selectedCube = -1;

	printf("Selected Cube %d\n", selectedCube);
}