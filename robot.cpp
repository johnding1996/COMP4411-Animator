// The sample model.  You should build a file
// very similar to this for when you make your model.
#include "modelerview.h"
#include "modelerapp.h"
#include "modelerdraw.h"
#include "camera.h"
#include <FL/gl.h>
#include <FL/glut.h>
#include <FL/Fl.H>
#include <math.h>
#include "vec.h"
#include "mat.h"
#include "texture.h"
#include "modelerglobals.h"
#include "metaball.h"
#include "ParticleSystem.h"

#define PI 3.14159265

// Animation Period (30Hz)
#define ANIMATION_PERIOD (30*3)


// To make a SampleModel, we inherit off of ModelerView
class SampleModel : public ModelerView
{
public:
	SampleModel(int x, int y, int w, int h, char *label) :
		ModelerView(x, y, w, h, label)
	{
	}
	virtual void draw();
private:
	int animationTime = 0;

	int animationHeadAngle;
	int animationUpperArmAngle;
	int animationForeArmAngle;
	int animationThighAngle;
	int animationShankAngle;
	double animationFireStrength;

	void drawHead();
	void drawHat();
	void drawFace();
	void drawNeck();

	void drawChest();
	void drawBelly();

	void drawRightShoulder();
	void drawRightUpperArm();
	void drawRightForeArm();
	void drawRightHand();

	void drawLeftShoulder();
	void drawLeftUpperArm();
	void drawLeftForeArm();
	void drawLeftHand();

	void drawRightThigh();
	void drawRightShank();
	void drawRightFoot();

	void drawLeftThigh();
	void drawLeftShank();
	void drawLeftFoot();

	void drawBackPack();

	void drawSkybox();

	void drawParticles();


	Mat4d CameraMatrix;
};

// We need to make a creator function, mostly because of
// nasty API stuff that we'd rather stay away from.
ModelerView* createSampleModel(int x, int y, int w, int h, char *label)
{
	return new SampleModel(x, y, w, h, label);
}

// We are going to override (is that the right word?) the draw()
// method of ModelerView to draw out SampleModel
void SampleModel::draw()
{
	ModelerView::draw();
	CameraMatrix = getModelViewMatrix();
	animationTime = 0;
	animationHeadAngle = 0;
	animationForeArmAngle = 0;
	animationUpperArmAngle = 0;
	animationThighAngle = 0;
	animationShankAngle = 0;
	animationFireStrength = 1.0;

	// draw the sample modelS
	setAmbientColor((float)(GETVAL(LIGHTR)*GETVAL(INTENSITY)) / 255,
		(float)(GETVAL(LIGHTG)*GETVAL(INTENSITY)) / 255,
		(float)(GETVAL(LIGHTB)*GETVAL(INTENSITY)) / 255);

	glPushMatrix(); // push identity
	glTranslated(GETVAL(XPOS), GETVAL(YPOS), GETVAL(ZPOS)); // values set by the sliders

	setDiffuseColor(1.0f, .863f, .694f);

	drawHead();
	drawHat();
	drawFace();
	drawNeck();

	drawChest();
	drawBelly();

	drawRightShoulder();
	drawRightUpperArm();
	drawRightForeArm();
	drawRightHand();

	drawLeftShoulder();
	drawLeftUpperArm();
	drawLeftForeArm();
	drawLeftHand();

	drawRightThigh();
	drawRightShank();
	drawRightFoot();

	drawLeftThigh();
	drawLeftShank();
	drawLeftFoot();

	drawBackPack();


	if (GETVAL(METABALLSKIN)) {
		setDiffuseColor(1.0f, .863f, .694f);
		MetaBallSystem mb((float)GETVAL(METABALLSIZE), (float)GETVAL(METABALLTHRESHOULD));
		mb.grid();
		mb.makeMetaballs();
		mb.evaluationField();
		mb.draw();
	}

	if (GETVAL(SKYBOX)) drawSkybox();

	glPopMatrix();
}
void SampleModel::drawHead() {
	glPushMatrix();
	glTranslated(0, CHEST_HALFHEIGHT + 0.1, 0);
	glRotated(animationHeadAngle, 0, 1, 0);
	glRotated(-90, 1, 0, 0);
	drawCylinder(HEAD_HALFHEIGHT * 2, 0.5, 0.4);
	glRotated(90, 1, 0, 0);
	glTranslated(0, HEAD_HALFHEIGHT, 0);
}
void SampleModel::drawHat() {
	glPushMatrix();
	glTranslated(0, HEAD_HALFHEIGHT, 0);
	glRotated(-90, 1, 0, 0);
	setDiffuseColor(0.6f, 0.2f, 0.2f);
	drawCylinder(0.4, 0.8, 0);
	if (GETVAL(TORUSHALO)) {
		setDiffuseColor(0.9f, 0.9f, 0.1f);
		glRotated(90, 1, 0, 0);
		glTranslated(0, 0.5, 0);
		drawTorus(0.5, 0.08);
	}
	glPopMatrix();
}
void SampleModel::drawFace() {
	glPushMatrix();
	// eyes
	setDiffuseColor(0, 0, 0);
	glTranslated(0.2, 0.3, 0.4);
	drawSphere(0.05);
	glTranslated(-0.4, 0, 0);
	drawSphere(0.05);
	// nose
	setDiffuseColor(1.0f, .8f, .6f);
	glTranslated(0.2, -0.3, 0);
	drawSphere(0.1);
	//mouth
	setDiffuseColor(1.0f, 0.0f, 0.1f);
	glTranslated(-0.2, -0.5, 0.08);
	glRotated(GETVAL(MOUTHANGLE), 1.0, 0.0, 0.0);
	glTranslated(0, -0.05, 0);
	drawBox(0.4, 0.1, 0);
	glPopMatrix();

	// reset color
	setDiffuseColor(1.0f, .863f, .694f);
	glPopMatrix();
}
void SampleModel::drawNeck() {
	glPushMatrix();
	glTranslated(0, CHEST_HALFHEIGHT, 0);
	glScaled(0.6, 1.0, 0.6);
	glTranslated(-0.5, -0.5, -0.5);
	drawBox(1, 1, 1);
	glPopMatrix();
}
void SampleModel::drawChest() {
	glPushMatrix();
	glTranslated(0, -CHEST_HALFHEIGHT, 0);
	glRotated(-90, 1, 0, 0);
	glScaled(1, GETVAL(FAT), 1);
	drawCylinder(CHEST_HALFHEIGHT * 2, WAIST_RADIUS, SHOULDER_RADIUS); // center at (0, 0, 0)
	glPopMatrix();
}
void SampleModel::drawBelly() {
	glPushMatrix();
	glTranslated(0, -CHEST_HALFHEIGHT, 0); // move down
	glRotated(90, 1.0, 0.0, 0.0);
	glScaled(1, GETVAL(FAT), 1);
	drawCylinder(BELLY_HALFHEIGHT, WAIST_RADIUS, HIP_RADIUS);
	glPopMatrix();
}
void SampleModel::drawRightShoulder() {
	glPushMatrix();
	glTranslated(SHOULDER_RADIUS, 0.6, 0);
	glRotated(20, 0.0, 0.0, 1.0);
	glRotated(90, 0.0, 1.0, 0.0);
	glTranslated(0, 0, -0.2);
	drawCylinder(0.4, 0.2, 0.2);
	glPopMatrix();

	glPushMatrix();
	glTranslated(SHOULDER_RADIUS + 0.2, 0.7, 0);
	drawSphere(0.3);
	glPopMatrix();
}
void SampleModel::drawRightUpperArm() {
	glPushMatrix();
	glTranslated(SHOULDER_RADIUS + 0.2, 0.8, 0);
	glRotated(GETVAL(RIGHTARMZ), 0, 0, 1.0);
	glRotated(GETVAL(RIGHTARMY), 0, -1.0, 0);
	glRotated(GETVAL(RIGHTARMX), -1.0, 0, 0);

	glRotated(20, 0, 0, 1);
	glTranslated(-UPPERARM_RADIUS / 2, 0, -UPPERARM_RADIUS / 2);
	drawBox(UPPERARM_RADIUS, -UPPERARM_LENGTH, UPPERARM_RADIUS);

	glTranslated(UPPERARM_RADIUS / 2, -UPPERARM_LENGTH, UPPERARM_RADIUS / 2);
	glRotated(-20, 0, 0, 1);
}
void SampleModel::drawRightForeArm() {
	drawSphere(0.2);
	glRotated(GETVAL(RIGHTELBOWX), -1.0, 0, 0);
	glRotated(GETVAL(RIGHTELBOWY), 0, 1.0, 0);

	glPushMatrix();
	glTranslated(0, -0.7, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.6);
	drawCylinder(FOREARM_LENGTH, FOREARM_RADIUS, FOREARM_RADIUS / 2);
	glPopMatrix();

}
void SampleModel::drawRightHand() {
	glTranslated(0, -1.4, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.1);
	glRotated(GETVAL(RIGHTHANDX), -1.0, 0, 0);
	drawSphere(0.1);


	if (GETVAL(DETAILEDHAND)) {
		glScaled(1, 2, 1);
		glRotated(30, 0, 1, 0);
		drawCylinder(HAND_LENGTH*0.7, 0.05, 0.03);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH*0.8, 0.05, 0.02);
	}
	else {
		glScaled(2, 1, 1);
		drawCylinder(HAND_LENGTH*0.8, 0.05, 0.1);
	}
	//Particle System
	drawParticles();
	glPopMatrix();
}
void SampleModel::drawLeftShoulder() {
	glPushMatrix();
	glTranslated(-SHOULDER_RADIUS, 0.6, 0);
	glRotated(20, 0.0, 0.0, 1.0);
	glRotated(90, 0.0, 1.0, 0.0);
	glTranslated(0, 0, -0.2);
	drawCylinder(0.4, 0.2, 0.2);
	glPopMatrix();

	glPushMatrix();
	glTranslated(-SHOULDER_RADIUS - 0.2, 0.7, 0);
	drawSphere(0.3);
	glPopMatrix();
}
void SampleModel::drawLeftUpperArm() {
	glPushMatrix();
	glTranslated(-SHOULDER_RADIUS - 0.2, 0.8, 0);
	glRotated(GETVAL(LEFTARMZ), 0, 0, 1.0);
	glRotated(GETVAL(LEFTARMY), 0, -1.0, 0);
	glRotated(GETVAL(LEFTARMX), -1.0, 0, 0);

	glRotated(-20, 0, 0, 1);
	glTranslated(-UPPERARM_RADIUS / 2, 0, -UPPERARM_RADIUS / 2);
	drawBox(UPPERARM_RADIUS, -UPPERARM_LENGTH, UPPERARM_RADIUS);
	glTranslated(UPPERARM_RADIUS / 2, -UPPERARM_LENGTH, UPPERARM_RADIUS / 2);
	glRotated(+20, 0, 0, 1);
}
void SampleModel::drawLeftForeArm() {
	drawSphere(0.2);
	glRotated(GETVAL(LEFTELBOWX), -1.0, 0, 0);
	glRotated(GETVAL(LEFTELBOWY), 0, 1.0, 0);

	glPushMatrix();
	glTranslated(0, -0.7, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.6);
	drawCylinder(FOREARM_LENGTH, FOREARM_RADIUS, FOREARM_RADIUS / 2);
	glPopMatrix();
}
void SampleModel::drawLeftHand() {
	glTranslated(0, -1.4, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.1);
	glRotated(GETVAL(LEFTHANDX), -1.0, 0, 0);
	drawSphere(0.1);
	if (GETVAL(DETAILEDHAND)) {
		glScaled(1, 2, 1);
		glRotated(30, 0, 1, 0);
		drawCylinder(HAND_LENGTH*0.7, 0.05, 0.03);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH, 0.05, 0.02);
		glRotated(-15, 0, 1, 0);
		drawCylinder(HAND_LENGTH*0.8, 0.05, 0.02);
	}
	else {
		glScaled(2, 1, 1);
		drawCylinder(HAND_LENGTH*0.8, 0.05, 0.1);
	}
	//Particle System
	drawParticles();
	glPopMatrix();
}
void SampleModel::drawRightThigh() {
	glPushMatrix();
	glTranslated(0.9, -CHEST_HALFHEIGHT - BELLY_HALFHEIGHT - 0.1, 0);
	glRotated(GETVAL(RIGHTKNEE), -1.0, 0.0, 0.0);
	glRotated(GETVAL(RIGHTLEGX), -1.0, 0, 0);
	glRotated(GETVAL(RIGHTLEGZ), 0.0, 0, 1.0);
	drawSphere(0.3);

	glPushMatrix();
	glTranslated(0, -0.7, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.6);
	drawCylinder(THIGH_LENGTH, THIGH_RADIUS, THIGH_RADIUS*0.75);
	glTranslated(0, 0, THIGH_LENGTH + 0.1);
	drawSphere(THIGH_RADIUS*0.75);
	glPopMatrix();
}
void SampleModel::drawRightShank() {
	glRotated(GETVAL(RIGHTKNEE), 1.0, 0.0, 0.0);
	glTranslated(0, -2.1, 0);
	glPushMatrix();
	glRotated(GETVAL(RIGHTKNEE) * 2, 1.0, 0.0, 0.0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.8);
	drawCylinder(SHANK_LENGTH, SHANK_RADIUS, SHANK_RADIUS*0.75);
	glPopMatrix();
}
void SampleModel::drawRightFoot() {

	glTranslated(0.15, -0.4, 0.4);
	glTranslated(0, 0, -GETVAL(RIGHTKNEE) / 180 * 2 * SHANK_LENGTH);
	glRotated(20, 0.0, 1.0, 0.0);
	glTranslated(0, 0, -0.5);
	glScaled(1, 0.9, 1);
	drawCylinder(0.9, 0.25, 0.1);
	glPopMatrix();
}
void SampleModel::drawLeftThigh() {
	glPushMatrix();
	glTranslated(-0.9, -CHEST_HALFHEIGHT - BELLY_HALFHEIGHT - 0.1, 0);
	glRotated(GETVAL(LEFTKNEE), -1.0, 0.0, 0.0);
	glRotated(GETVAL(LEFTLEGX), -1.0, 0, 0);
	glRotated(GETVAL(LEFTLEGZ), 0.0, 0, 1.0);
	drawSphere(0.3);

	glPushMatrix();
	glTranslated(0, -0.7, 0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.6);
	drawCylinder(THIGH_LENGTH, THIGH_RADIUS, THIGH_RADIUS*0.75);
	glTranslated(0, 0, THIGH_LENGTH + 0.1);
	drawSphere(THIGH_RADIUS*0.75);
	glPopMatrix();
}
void SampleModel::drawLeftShank() {

	glRotated(GETVAL(LEFTKNEE), 1.0, 0.0, 0.0);
	glTranslated(0, -2.1, 0);
	glPushMatrix();
	glRotated(GETVAL(LEFTKNEE) * 2, 1.0, 0.0, 0.0);
	glRotated(90, 1.0, 0.0, 0.0);
	glTranslated(0, 0, -0.8);
	drawCylinder(SHANK_LENGTH, SHANK_RADIUS, SHANK_RADIUS*0.75);
	glPopMatrix();
}
void SampleModel::drawLeftFoot() {

	glTranslated(-0.15, -0.4, 0.4);
	glTranslated(0, 0, -GETVAL(LEFTKNEE) / 180 * 2 * SHANK_LENGTH);
	glRotated(-20, 0.0, 1.0, 0.0);
	glTranslated(0, 0, -0.5);
	glScaled(1, 0.9, 1);
	drawCylinder(0.9, 0.25, 0.1);
	glPopMatrix();
}
void SampleModel::drawBackPack() {
	glPushMatrix();
	setDiffuseColor(.60f, .20f, .10f);
	glTranslated(0, 0, -GETVAL(FAT)*SHOULDER_RADIUS);
	glRotated(180, 0, 1, 0);
	glScaled(0.7, 1, GETVAL(FAT));

	// back
	drawTriangle(0, 1.2, 1.0,
		0, 0.4, 1.3,
		1.1, 0.9, 0.7);
	drawTriangle(0, 1.2, 1.0,
		0, 0.4, 1.3,
		-1.1, 0.9, 0.7);
	drawTriangle(0, 0.4, 1.3,
		1.1, 0.9, 0.7,
		1.1, 0.1, 0.7);
	drawTriangle(0, 0.4, 1.3,
		-1.1, 0.9, 0.7,
		-1.1, 0.1, 0.7);
	drawTriangle(0, 0.4, 1.3,
		1.1, 0.1, 0.7,
		0, -1.0, 1.3);
	drawTriangle(0, 0.4, 1.3,
		-1.1, 0.1, 0.7,
		0, -1.0, 1.3);
	drawTriangle(1.1, 0.1, 0.7,
		0, -1.0, 1.3,
		1.1, -1.5, 0.7);
	drawTriangle(-1.1, 0.1, 0.7,
		0, -1.0, 1.3,
		-1.1, -1.5, 0.7);
	drawTriangle(0, -1.0, 1.3,
		0, -1.9, 1.0,
		1.1, -1.5, 0.7);
	drawTriangle(0, -1.0, 1.3,
		0, -1.9, 1.0,
		-1.1, -1.5, 0.7);

	// side
	drawTriangle(0, 1.2, 1.0,
		1.1, 0.9, 0.7,
		1.3, 1.0, -0.7);
	drawTriangle(1.3, 1.0, -0.7,
		0.0, 1.2, -0.7,
		0, 1.2, 1.0);
	drawTriangle(0, 1.2, 1.0,
		-1.1, 0.9, 0.7,
		0.0, 1.2, -0.7);
	drawTriangle(-1.3, 1.0, -0.7,
		0.0, 1.2, -0.7,
		-1.1, 0.9, 0.7);
	drawTriangle(-1.3, 1.0, -0.7,
		-1.1, 0.9, 0.7,
		-1.3, -1.3, -0.7);
	drawTriangle(-1.1, 0.9, 0.7,
		-1.1, -1.5, 0.7,
		-1.3, -1.3, -0.7);
	drawTriangle(1.3, -1.5, -0.7,
		0, -1.9, 1.0,
		1.1, -1.5, 0.7);
	drawTriangle(1.3, -1.5, -0.7,
		0.0, -1.6, -0.7,
		0, -1.9, 1.0);
	drawTriangle(-1.3, -1.3, -0.7,
		0.0, -1.6, -0.7,
		0, -1.9, 1.0);
	drawTriangle(0, -1.9, 1.0,
		-1.1, -1.5, 0.7,
		-1.3, -1.3, -0.7);
	drawTriangle(1.1, 0.9, 0.7,
		1.1, -1.5, 0.7,
		1.3, -1.5, -0.7);
	drawTriangle(1.3, 1.0, -0.7,
		1.1, 0.9, 0.7,
		1.3, -1.5, -0.7);

	// front
	drawTriangle(-1.3, 1.0, -0.7,
		-1.3, -1.3, -0.7,
		0.0, 1.2, -0.7);
	drawTriangle(1.3, 1.0, -0.7,
		1.3, -1.3, -0.7,
		0.0, 1.2, -0.7);
	drawTriangle(0, 1.2, -0.7,
		-1.3, -1.3, -0.7,
		0.0, -1.5, -0.7);
	drawTriangle(0, 1.2, -0.7,
		1.3, -1.3, -0.7,
		0.0, -1.5, -0.7);

	if (GETVAL(JETBACKPACK)) {

		glPushMatrix();
		setDiffuseColor(.35f, .35f, .35f);
		glTranslated(-0.9, -1.5, 0.6);
		glRotated(90, 1, 0, 0);
		drawCylinder(0.4, 0.3, 0.6);
		glTranslated(0, 0, 0.4);
		setDiffuseColor(0.9f, 0.9f, 0.0f);
		drawCylinder(1.0*animationFireStrength, 0.2*(2 - animationFireStrength), 0);
		glPopMatrix();

		glPushMatrix();
		setDiffuseColor(.35f, .35f, .35f);
		glTranslated(0.9, -1.5, 0.6);
		glRotated(90, 1, 0, 0);
		drawCylinder(0.4, 0.3, 0.6);
		glTranslated(0, 0, 0.4);
		setDiffuseColor(0.9f, 0.9f, 0.0f);
		drawCylinder(1.0*animationFireStrength, 0.2*(2 - animationFireStrength), 0);
		glPopMatrix();

	}

	glPopMatrix();
}

// Draw Particles
void SampleModel::drawParticles()
{
	Mat4d WorldMatrix = CameraMatrix.inverse() * getModelViewMatrix();
	Vec4d pos = WorldMatrix * Vec4d(0, 0, 0.2, 1);
	ParticleSystem *ps = ModelerApplication::Instance()->GetParticleSystem();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	ps->SpawnParticles(Vec3f((float)pos[0], (float)pos[1], (float)pos[2]), (int)GETVAL(PARTICLE_NUM), m_camera);
	glPopAttrib();
}
void SampleModel::drawSkybox()
{
	// Store the current matrix
	glPushMatrix();
	// Reset and transform the matrix.

	// Enable/Disable features
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_BLEND);
	glColor3f(1, 1, 1);

	// Just in case we set all vertices to white.
	setAmbientColor(0.4f, 0.4f, 0.4f);
	Vec3f cameraPosition = m_camera->getPosition();
	glTranslated(cameraPosition[0], cameraPosition[1], cameraPosition[2]);
	glScaled(80,80,80);
	Texture skyboxTexture;
	// Render the front quad
	GLuint textureID = skyboxTexture.loadBMP_custom("Image/cliffFront.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(1, 0); glVertex3f(0.5f, -0.5f, -0.5f);
	glTexCoord2f(0, 0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glTexCoord2f(0, 1); glVertex3f(-0.5f, 0.5f, -0.5f);
	glTexCoord2f(1, 1); glVertex3f(0.5f, 0.5f, -0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);


	// Render the left quad
	textureID = skyboxTexture.loadBMP_custom("Image/cliffLeft.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(1, 0); glVertex3f(0.5f, -0.5f, 0.5f);
	glTexCoord2f(0, 0); glVertex3f(0.5f, -0.5f, -0.5f);
	glTexCoord2f(0, 1); glVertex3f(0.5f, 0.5f, -0.5f);
	glTexCoord2f(1, 1); glVertex3f(0.5f, 0.5f, 0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);


	// Render the back quad
	textureID = skyboxTexture.loadBMP_custom("Image/cliffBack.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(1, 0); glVertex3f(-0.5f, -0.5f, 0.5f);
	glTexCoord2f(0, 0); glVertex3f(0.5f, -0.5f, 0.5f);
	glTexCoord2f(0, 1); glVertex3f(0.5f, 0.5f, 0.5f);
	glTexCoord2f(1, 1); glVertex3f(-0.5f, 0.5f, 0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);

	// Render the right quad
	textureID = skyboxTexture.loadBMP_custom("Image/cliffRight.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(1, 0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glTexCoord2f(0, 0); glVertex3f(-0.5f, -0.5f, 0.5f);
	glTexCoord2f(0, 1); glVertex3f(-0.5f, 0.5f, 0.5f);
	glTexCoord2f(1, 1); glVertex3f(-0.5f, 0.5f, -0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);

	// Render the top quad
	textureID = skyboxTexture.loadBMP_custom("Image/cliffTop.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(0, 1); glVertex3f(-0.5f, 0.5f, 0.5f);
	glTexCoord2f(0, 0); glVertex3f(0.5f, 0.5f, 0.5f);
	glTexCoord2f(1, 0); glVertex3f(0.5f, 0.5f, -0.5f);
	glTexCoord2f(1, 1); glVertex3f(-0.5f, 0.5f, -0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);

	// Render the bottom quad
	textureID = skyboxTexture.loadBMP_custom("Image/cliffBottom.bmp");
	glBegin(GL_QUADS);
	glTexCoord2f(0, 0); glVertex3f(-0.5f, -0.5f, -0.5f);
	glTexCoord2f(0, 1); glVertex3f(-0.5f, -0.5f, 0.5f);
	glTexCoord2f(1, 1); glVertex3f(0.5f, -0.5f, 0.5f);
	glTexCoord2f(1, 0); glVertex3f(0.5f, -0.5f, -0.5f);
	glEnd();
	glDeleteTextures(1, &textureID);

	// Clear the texture memory
	glFinish();
	// Restore enable bits and matrix
	glPopAttrib();
	glPopMatrix();
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);

}

int main()
{
	// Initialize the controls
	// Constructor is ModelerControl(name, minimumvalue, maximumvalue, stepsize, defaultvalue)
	ModelerControl controls[NUMCONTROLS];
	controls[LIGHTR] = ModelerControl("Light Source Red", 0, 255, 1, 255);
	controls[LIGHTG] = ModelerControl("Light Source Green", 0, 255, 1, 255);
	controls[LIGHTB] = ModelerControl("Light Source Blue", 0, 255, 1, 255);
	controls[INTENSITY] = ModelerControl("Light Intensity", 0, 3, 0.1f, 1);
	controls[XPOS] = ModelerControl("X Position", -10, 10, 0.1f, 0);
	controls[YPOS] = ModelerControl("Y Position", -10, 10, 0.1f, 0);
	controls[ZPOS] = ModelerControl("Z Position", -10, 10, 0.1f, 0);
	controls[RIGHTARMX] = ModelerControl("Right Upper Arm Angle X", -30, 180, 1, 0);
	controls[RIGHTARMY] = ModelerControl("Right Upper Arm Angle Y", -30, 60, 1, 0);
	controls[RIGHTARMZ] = ModelerControl("Right Upper Arm Angle Z", 0, 120, 1, 0);
	controls[RIGHTELBOWX] = ModelerControl("Right Fore Arm Angle X", -10, 120, 1, 0);
	controls[RIGHTELBOWY] = ModelerControl("Right Fore Arm Angle Y", -90, 90, 1, 0);
	controls[RIGHTHANDX] = ModelerControl("Right Hand Angle X", 60, -60, 1, 0);
	controls[LEFTARMX] = ModelerControl("Left Upper Arm Angle X", -30, 180, 1, 0);
	controls[LEFTARMY] = ModelerControl("Left Upper Arm Angle Y", -60, 30, 1, 0);
	controls[LEFTARMZ] = ModelerControl("Left Upper Arm Angle Z", -120, 0, 1, 0);
	controls[LEFTELBOWX] = ModelerControl("Left Fore Arm Angle X", -10, 120, 1, 0);
	controls[LEFTELBOWY] = ModelerControl("Left Fore Arm Angle Y", -90, 90, 1, 0);
	controls[LEFTHANDX] = ModelerControl("Left Hand Angle X", 60, -60, 1, 0);
	controls[RIGHTLEGX] = ModelerControl("Right Leg Angle X", -60, 60, 1, 0);
	controls[RIGHTLEGZ] = ModelerControl("Right Leg Angle Z", -30, 60, 1, 0);
	controls[RIGHTKNEE] = ModelerControl("Right Knee Angle X", 0, 10, 1, 0);
	controls[LEFTLEGX] = ModelerControl("Left Leg Angle X", -60, 60, 1, 0);
	controls[LEFTLEGZ] = ModelerControl("Left Leg Angle Z", -60, 30, 1, 0);
	controls[LEFTKNEE] = ModelerControl("Left Knee Angle X", 0, 10, 1, 0);
	controls[FAT] = ModelerControl("Plumpness", 0.3f, 1.0f, 0.05f, 0.6f);
	controls[MOUTHANGLE] = ModelerControl("Mouth Angle (Happy or Sad)", -40, 40, 1, -40);
	controls[DETAILEDHAND] = ModelerControl("Enable Detailed Hands", 0, 1, 1, 1);
	controls[JETBACKPACK] = ModelerControl("Enable Jet Backpack", 0, 1, 1, 0);
	controls[TORUSHALO] = ModelerControl("Enable Torus Halo", 0, 1, 1, 0);
	controls[METABALLSKIN] = ModelerControl("Enable Metaball Skin", 0, 1, 1, 0);
	controls[METABALLSIZE] = ModelerControl("Metaball Radius", 0.0f, 0.8f, 0.005f, 0.2f);
	controls[METABALLTHRESHOULD] = ModelerControl("Metaball Threshould", 0.0f, 1.0f, 0.01f, 0.50f);
	controls[SKYBOX] = ModelerControl("Enable Skybox Background", 0, 1, 1, 1);
	controls[PARTICLE_NUM] = ModelerControl("Number of particles", 0, 50, 1, 20);

	// You should create a ParticleSystem object ps here and then
	// call ModelerApplication::Instance()->SetParticleSystem(ps)
	// to hook it up to the animator interface.
	ParticleSystem *ps = new ParticleSystem(3, 0.5);
	ModelerApplication::Instance()->SetParticleSystem(ps);


	ModelerApplication::Instance()->Init(&createSampleModel, controls, NUMCONTROLS);
	return ModelerApplication::Instance()->Run();
}
