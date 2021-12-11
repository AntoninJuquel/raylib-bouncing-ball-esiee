/*******************************************************************************************
*
*   raylib [core] example - Basic window
*
*   Welcome to raylib!
*
*   To test examples, just press F6 and execute raylib_compile_execute script
*   Note that compiled executable is placed in the same folder as .c file
*
*   You can find all basic examples on C:\raylib\raylib\examples folder or
*   raylib official webpage: www.raylib.com
*
*   Enjoy using raylib. :)
*
*   This example has been created using raylib 1.0 (www.raylib.com)
*   raylib is licensed under an unmodified zlib/libpng license (View raylib.h for details)
*
*   Copyright (c) 2014 Ramon Santamaria (@raysan5)
*
********************************************************************************************/

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

#pragma region Structs
struct Cylindrical {
	float rho;
	float theta;
	float y;

	inline Cylindrical operator+(Cylindrical a) {
		return { a.rho + rho,a.theta + theta,a.y + y };
	}
};

struct Spherical {
	float rho;
	float theta;
	float phi;

	inline Spherical operator+(Spherical a) {
		return { a.rho + rho,a.theta + theta,a.phi + phi };
	}

};

struct Segment
{
	Vector3 p1, p2;
};

struct Sphere
{
	Quaternion rotation;
	Vector3 position;
	float radius;
};

struct Plane
{
	Quaternion rotation;
	Vector3 position;
	Vector2 scale;
};

struct Cylinder
{
	Quaternion rotation;
	Vector3 position;
	Vector2 scale;
};

struct Disk
{
	Quaternion rotation;
	Vector3 position;
	float radius;
};

struct Referential
{
	Vector3 origin;
	Vector3 i, j, k;
};

struct Capsule
{
	Vector3 position;
	Quaternion rotation;
	Vector2 scale;

	Sphere sph1;
	Sphere sph2;
	Cylinder cyl;
};

struct Box
{
	Vector3 position;
	Vector3 scale;
	Quaternion rotation;

	Plane quads[6] = {};
};

struct Ball
{
	Vector3 position;
	Vector3 velocity;
	float radius;
	float bounciness;
};
#pragma endregion

#pragma region Conversions
Cylindrical CartesianToCylindrical(Vector3 cart)
{
	Cylindrical cyl;
	cyl.rho = sqrtf(cart.x * cart.x + cart.z * cart.z);
	cyl.y = cart.y;

	if (cyl.rho < EPSILON)cyl.theta = 0;
	else
	{
		cyl.theta = atan2f(cart.x, cart.z);
		if (cyl.theta < 0)cyl.theta += PI * 2;
	}
	return cyl;
}

Vector3 CylindricalToCartesian(Cylindrical cyl)
{
	return Vector3{ cyl.rho * sinf(cyl.theta),cyl.y,cyl.rho * cosf(cyl.theta) };
}

Vector3 SphericalToCartesian(Spherical sph)
{
	return Vector3{ sph.rho * sinf(sph.phi) * sinf(sph.theta),
	sph.rho * cosf(sph.phi),
	sph.rho * sinf(sph.phi) * cosf(sph.theta) };
}
#pragma endregion

#pragma region Drawers
void MyDrawSphere(Sphere sphere, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(sphere.position.x, sphere.position.y, sphere.position.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//

	rlScalef(sphere.radius, sphere.radius, sphere.radius);


	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ 1,theta + deltaTheta,phi + deltaPhi });


			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = tmpBottomLeft;
			tmpBottomLeft = bottomRight;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}

void MyDrawSphereWires(Sphere sphere, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 4;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(sphere.position.x, sphere.position.y, sphere.position.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(sphere.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//

	rlScalef(sphere.radius, sphere.radius, sphere.radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });
			Vector3 topRight = vertexBufferTheta[j + 1];

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = bottomLeft;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}

void MyDrawSpherePortion(Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {

}

void MyDrawSphereWiresPortion(Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {

}


void MyDrawQuad(Plane plane, Color color) {

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(plane.position.x, plane.position.y, plane.position.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(plane.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Front face
	float width = plane.scale.x;
	float length = plane.scale.y;

	// by default facing up 

	rlVertex3f(+width / 2, 0, -length / 2);  // Top Left
	rlVertex3f(-width / 2, 0, -length / 2);  // Bottom Left
	rlVertex3f(-width / 2, 0, +length / 2);  // Bottom Right

	rlVertex3f(+width / 2, 0, +length / 2);  // Top Right
	rlVertex3f(+width / 2, 0, -length / 2);  // Top Left
	rlVertex3f(-width / 2, 0, +length / 2);  // Bottom Right

	rlEnd();
	rlPopMatrix();
}
void MyDrawQuadWire(Plane plane, Color color) {

	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	float width = plane.scale.x;
	float height = plane.scale.y;
	float length = height;


	if (rlCheckBufferLimit(36)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(plane.position.x, plane.position.y, plane.position.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(plane.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// facing up by default

	rlVertex3f(x - width / 2, y, z - length / 2);  // Bottom Left
	rlVertex3f(x + width / 2, y, z - length / 2);  // Bottom Right

	// Left Line
	rlVertex3f(x + width / 2, y, z - length / 2);  // Bottom Right
	rlVertex3f(x + width / 2, y, z + length / 2);  // Top Right

	// Top Line
	rlVertex3f(x + width / 2, y, z + length / 2);  // Top Right
	rlVertex3f(x - width / 2, y, z + length / 2);  // Top Left

	// Right Line
	rlVertex3f(x - width / 2, y, z + length / 2);  // Top Left
	rlVertex3f(x - width / 2, y, z - length / 2);  // Bottom Left

	//Diagonal
	rlVertex3f(x - width / 2, y, z + length / 2);  // Top Left
	rlVertex3f(x + width / 2, y, z - length / 2);  // Bottom Right

	rlVertex3f(x + width / 2, y, z + length / 2);  // Top Right
	rlVertex3f(x - width / 2, y, z - length / 2);  // Bottom Left

	rlEnd();
	rlPopMatrix();
}

void MyDrawDisk(Disk disk, int nSegmentsTheta, Color color) {
	int sides = 100;
	int numVertex = sides * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();
	Vector3 position = disk.position;
	rlTranslatef(position.x, position.y, position.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(disk.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	for (int i = 0; i < 360; i += 360 / sides)
	{
		rlVertex3f(0, 0, 0);
		rlVertex3f(sinf(DEG2RAD * i) * disk.radius, 0, cosf(DEG2RAD * i) * disk.radius);
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * disk.radius, 0, cosf(DEG2RAD * (i + 360 / sides)) * disk.radius);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinder(Cylinder cyl, int nSegmentsTheta, Color color) {
	int sides = 100;
	int numVertex = sides * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();
	Vector3 position = cyl.position;
	rlTranslatef(position.x, position.y, position.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cyl.rotation, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	float radius = cyl.scale.x;
	float height = cyl.scale.y;

	Vector3 top = { 0,height * .5f,0 };
	Vector3 bottom = { 0,-height * .5f,0 };


	// Draw Body -------------------------------------------------------------------------------------
	for (int i = 0; i < 360; i += 360 / sides)
	{

		rlVertex3f(sinf(DEG2RAD * i) * radius, bottom.y, cosf(DEG2RAD * i) * radius); //Bottom Left
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * radius, bottom.y, cosf(DEG2RAD * (i + 360 / sides)) * radius); //Bottom Right
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * radius, top.y, cosf(DEG2RAD * (i + 360 / sides)) * radius); //Top Right
		rlVertex3f(sinf(DEG2RAD * i) * radius, top.y, cosf(DEG2RAD * i) * radius); //Top Left
		rlVertex3f(sinf(DEG2RAD * i) * radius, bottom.y, cosf(DEG2RAD * i) * radius); //Bottom Left
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * radius, top.y, cosf(DEG2RAD * (i + 360 / sides)) * radius); //Top Right
	}
	Disk dT = Disk{};
	dT.position = { 0,height * 0.5f,0 };
	dT.radius = radius;
	dT.rotation = cyl.rotation;
	//
	// TODO racalculer la rotation pour face la base
	//MyDrawDisk(dT, nSegmentsTheta, color);

	// Draw Cap --------------------------------------------------------------------------------------
	for (int i = 0; i < 360; i += 360 / sides)
	{
		rlVertex3f(0, top.y, 0);
		rlVertex3f(sinf(DEG2RAD * i) * radius, top.y, cosf(DEG2RAD * i) * radius);
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * radius, top.y, cosf(DEG2RAD * (i + 360 / sides)) * radius);
	}
	// Draw Base -----------------------------------------------------------------------------------------
	for (int i = 0; i < 360; i += 360 / sides)
	{
		rlVertex3f(0, bottom.y, 0);
		rlVertex3f(sinf(DEG2RAD * (i + 360 / sides)) * radius, bottom.y, cosf(DEG2RAD * (i + 360 / sides)) * radius);
		rlVertex3f(sinf(DEG2RAD * i) * radius, bottom.y, cosf(DEG2RAD * i) * radius);
	}
	rlEnd();
	rlPopMatrix();
}
void MyDrawCylinderWires(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {

}

void MyDrawCapsule(Capsule capsule, float nSegmentsTheta, Color color) {
	Cylinder cyl = capsule.cyl;
	MyDrawCylinder(cyl, nSegmentsTheta, color);
	Sphere s1 = capsule.sph1;
	MyDrawSphere(s1, nSegmentsTheta, nSegmentsTheta, color);
	Sphere s2 = capsule.sph2;
	MyDrawSphere(s2, nSegmentsTheta, nSegmentsTheta, color);
}
void MyDrawBox(Box box, Color color) {
	for each (Plane quad in box.quads)
	{
		MyDrawQuad(quad, color);
	}
}
void MyDrawDiskWires(Quaternion q, Vector3 position, float radius, int nSegmentsTheta, Color color) {

}

void MyDrawDiskPortion(Quaternion q, Vector3 position, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color) {

}

void MyDrawDiskWiresPortion(Quaternion q, Vector3 position, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color) {

}

void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {

}
void MyDrawCylinderWiresPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {

}
#pragma endregion

#pragma region Intersections
bool InterSegmentSphere(Segment seg, Sphere s, Vector3& interPt, Vector3& interNormal) {

	double a, b, c, mu1, mu2;
	double bb4ac;
	Vector3 dp = Vector3Subtract(seg.p2, seg.p1);

	a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
	b = 2 * (dp.x * (seg.p1.x - s.position.x) + dp.y * (seg.p1.y - s.position.y) + dp.z * (seg.p1.z - s.position.z));
	c = s.position.x * s.position.x + s.position.y * s.position.y + s.position.z * s.position.z;
	c += seg.p1.x * seg.p1.x + seg.p1.y * seg.p1.y + seg.p1.z * seg.p1.z;
	c -= 2 * (s.position.x * seg.p1.x + s.position.y * seg.p1.y + s.position.z * seg.p1.z);
	c -= s.radius * s.radius;
	bb4ac = b * b - 4 * a * c;

	if (abs(a) < EPSILON || bb4ac < 0) {
		mu1 = 0;
		mu2 = 0;
		return false;
	}

	mu1 = (-b + sqrt(bb4ac)) / (2 * a);
	mu2 = (-b - sqrt(bb4ac)) / (2 * a);

	interPt = Vector3Add(seg.p1, Vector3Scale(dp, mu1));
	interNormal = Vector3Normalize(Vector3Subtract(interPt, s.position));

	return true;
}

bool InterSegmentPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal) {

	Vector3 diff = Vector3Subtract(seg.p1, plane.position);
	Vector3 lineVector = Vector3Normalize(Vector3Subtract(seg.p2, seg.p1));
	Vector3 planeNormal = Vector3RotateByQuaternion({ 0,1,0 }, plane.rotation);
	Vector3 planePoint = plane.position;

	interPt = Vector3Add(Vector3Add(diff, planePoint), Vector3Scale(lineVector, -Vector3DotProduct(diff, planeNormal) / Vector3DotProduct(lineVector, planeNormal)));
	interNormal = planeNormal;

	return Vector3Distance(seg.p1, interPt) <= Vector3Distance(seg.p1, seg.p2);
}

bool InterSegmentQuad(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal) {
	if (InterSegmentPlane(seg, plane, interPt, interNormal)) {
		//Determiner C A B les sommets du quad
		//Vérifier que 0 < AinterPt . AC < AC . AC et que 0 < AinterPt . AB < AB . AB

		Vector3 A = { plane.position.x - plane.scale.x * .5f, plane.position.y,  plane.position.z - plane.scale.y * .5f };
		Vector3 B = { plane.position.x - plane.scale.x * .5f, plane.position.y,  plane.position.z + plane.scale.y * .5f };
		Vector3 C = { plane.position.x + plane.scale.x * .5f, plane.position.y,  plane.position.z - plane.scale.y * .5f };

		A = Vector3RotateByQuaternion(A, plane.rotation);
		B = Vector3RotateByQuaternion(B, plane.rotation);
		C = Vector3RotateByQuaternion(C, plane.rotation);

		Vector3 AinterPt = Vector3Subtract(interPt, A);
		Vector3 AC = Vector3Subtract(C, A);
		Vector3 AB = Vector3Subtract(B, A);

		float dot1 = Vector3DotProduct(AinterPt, AC);
		float dot2 = Vector3DotProduct(AinterPt, AB);

		return  0 < dot1 && dot1 < Vector3DotProduct(AC, AC) && 0 < dot2 && dot2 < Vector3DotProduct(AB, AB);
	}

	return false;
}

bool InterSegmentBox(Segment seg, Box box, Vector3& interPt, Vector3& interNormal) {
	for each (Plane plane in box.quads) {
		if (InterSegmentQuad(seg, plane, interPt, interNormal)) 
			return true;
	}

	return false;
}

bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	return false;
}

bool InterSegmentFiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	return false;
}
bool InterSegmentCapsule(Segment seg, Capsule capsule, Vector3& interPt, Vector3& interNormal) {

	Sphere s1 = capsule.sph1;
	Sphere s2 = capsule.sph2;
	Cylinder cyl = capsule.cyl;

	return InterSegmentFiniteCylinder(seg, cyl, interPt, interNormal) || InterSegmentSphere(seg, s1, interPt, interNormal) || InterSegmentSphere(seg, s2, interPt, interNormal);
}

bool InterSegmentDisk(Segment seg, Disk disk, Vector3& interPt, Vector3& interNormal) {
	return false;
}
#pragma endregion

#pragma region Methods
void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 10,PI / 4.f,PI / 4.f };
	static Spherical sphSpeed = { 10,.4f,.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime : 0,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime : 0 };

	Spherical newSphPos = sphPos + sphDelta;
	newSphPos = { Clamp(newSphPos.rho,rhoMin,rhoMax),
	newSphPos.theta,
	Clamp(newSphPos.phi,PI / 100.f,.99f * PI) };

	sphPos = newSphPos;

	camera->position = SphericalToCartesian(sphPos);
}

void UpdateBall(Ball* ball, float deltaTime, std::vector<Plane> planes, Box box) {
	ball->velocity.y += -9.81 * deltaTime;


	Vector3 nextPosition = Vector3Add(ball->position, Vector3Scale(ball->velocity, deltaTime));

	Segment collisionSeg = {};
	collisionSeg.p1 = ball->position;
	collisionSeg.p2 = nextPosition;

	for each (Plane plane in planes)
	{
		Vector3 interPt;
		Vector3 interNormal;
		if (InterSegmentQuad(collisionSeg, plane, interPt, interNormal)) {
			ball->velocity = Vector3Add(ball->velocity, Vector3Scale(interNormal, ball->bounciness));
			nextPosition = Vector3Add(ball->position, Vector3Scale(ball->velocity, deltaTime));
		}
		collisionSeg.p1 = ball->position;
		collisionSeg.p2 = nextPosition;
	}

	Vector3 interPt;
	Vector3 interNormal;
	
	if (InterSegmentBox(collisionSeg, box, interPt, interNormal)) {
		ball->velocity = Vector3Add(ball->velocity, Vector3Scale(interNormal, ball->bounciness));
		nextPosition = Vector3Add(ball->position, Vector3Scale(ball->velocity, deltaTime));
	}
	collisionSeg.p1 = ball->position;
	collisionSeg.p2 = nextPosition;

	ball->velocity = Vector3Scale(Vector3Normalize(ball->velocity), 5);
	ball->position = nextPosition;
}

void CreateCapsule(Capsule* capsule) {
	Cylinder cyl = {};
	cyl.position = capsule->position;
	cyl.scale = capsule->scale;
	cyl.rotation = capsule->rotation;

	Vector3 XAxis = Vector3RotateByQuaternion({ 0,1,0 }, capsule->rotation);

	Sphere s1 = {};
	s1.position = Vector3Add(capsule->position, Vector3Scale(XAxis, capsule->scale.y * 0.5f));
	s1.radius = capsule->scale.x;

	Sphere s2 = {};
	s2.position = Vector3Add(capsule->position, Vector3Scale(XAxis, -capsule->scale.y * 0.5f));
	s2.radius = capsule->scale.x;

	capsule->cyl = cyl;
	capsule->sph1 = s1;
	capsule->sph2 = s2;
}

void CreateBox(Box* box) {
	Plane planes[6] = {};

	Vector3 XAxis = Vector3RotateByQuaternion({ 1,0,0 }, box->rotation);
	Vector3 YAxis = Vector3RotateByQuaternion({ 0,1,0 }, box->rotation);
	Vector3 ZAxis = Vector3RotateByQuaternion({ 0,0,1 }, box->rotation);

	planes[0].position = Vector3Scale(ZAxis, box->scale.z * .5f);
	planes[0].scale = { box->scale.x, box->scale.y };
	planes[0].rotation = QuaternionFromAxisAngle({ 1,0,0 }, 90 * DEG2RAD);

	planes[1].position = Vector3Scale(ZAxis, -box->scale.z * .5f);
	planes[1].scale = { box->scale.x, box->scale.y };
	planes[1].rotation = QuaternionFromAxisAngle({ 1,0,0 }, -90 * DEG2RAD);

	planes[2].position = Vector3Scale(YAxis, box->scale.y * 0.5f);
	planes[2].scale = { box->scale.x, box->scale.z };
	planes[2].rotation = QuaternionIdentity();

	planes[3].position = Vector3Scale(YAxis, -box->scale.y * 0.5f);
	planes[3].scale = { box->scale.x, box->scale.z };
	planes[3].rotation = QuaternionFromAxisAngle({ 1,0,0 }, 180 * DEG2RAD);

	planes[4].position = Vector3Scale(XAxis, box->scale.x * .5f);
	planes[4].scale = { box->scale.z, box->scale.y };
	planes[4].rotation = QuaternionMultiply(QuaternionFromAxisAngle({ 0,1,0 }, -90 * DEG2RAD), QuaternionFromAxisAngle({ 1,0,0 }, -90 * DEG2RAD)); //QuaternionMultiply(box->rotation, QuaternionFromAxisAngle({ 0,0,1 }, 90 * DEG2RAD));

	planes[5].position = Vector3Scale(XAxis, -box->scale.x * .5f);
	planes[5].scale = { box->scale.z, box->scale.y };
	planes[5].rotation = QuaternionMultiply(QuaternionFromAxisAngle({ 0,1,0 }, -90 * DEG2RAD), QuaternionFromAxisAngle({ 1,0,0 }, 90 * DEG2RAD));

	for (size_t i = 0; i < 6; i++)
	{
		box->quads[i].position = Vector3Add(box->position, planes[i].position);
		box->quads[i].scale = planes[i].scale;
		box->quads[i].rotation = QuaternionMultiply(box->rotation, planes[i].rotation);
	}
}
#pragma endregion

int main(int argc, char* argv[])
{
	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = 0.7f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "Bouncy Sphere");

	SetTargetFPS(60);

	//CAMERA
	Vector3 cameraPos = { 8.0f, 15.0f, 14.0f };
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0.0f, 0.0f, 0.0f };
	camera.up = { 0.0f, 1.0f, 0.0f };
	camera.fovy = 45.0f;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//BALL
	Ball ball = {};
	ball.position = { 0,2,0 };
	ball.radius = .1f;
	ball.velocity = { 0,0,0 };
	ball.bounciness = 10;

	//Planes
	std::vector<Plane> planes;

	Plane plane = {};
	plane.position = { 0,0,0 };
	plane.scale = { 5, 5 };
	plane.rotation = QuaternionIdentity();

	Plane plane1 = {};
	plane1.position = { 2.5f,2.5f,0 };
	plane1.scale = { 5, 5 };
	plane1.rotation = QuaternionFromAxisAngle({ 0,0,1 }, 90 * DEG2RAD);

	Plane plane2 = {};
	plane2.position = { -2.5f,2.5f,0 };
	plane2.scale = { 5, 5 };
	plane2.rotation = QuaternionFromAxisAngle({ 0,0,1 }, -90 * DEG2RAD);

	Plane plane3 = {};
	plane3.position = { 0,2.5f,2.5f };
	plane3.scale = { 5, 5 };
	plane3.rotation = QuaternionFromAxisAngle({ 1,0,0 }, -90 * DEG2RAD);

	Plane plane4 = {};
	plane4.position = { 0,2.5f,-2.5f };
	plane4.scale = { 5, 5 };
	plane4.rotation = QuaternionFromAxisAngle({ 1,0,0 }, 90 * DEG2RAD);

	Plane plane5 = {};
	plane5.position = { 0,5,0 };
	plane5.scale = { 5, 5 };
	plane5.rotation = QuaternionFromAxisAngle({ 1,0,0 }, 180 * DEG2RAD);

	planes.push_back(plane);
	planes.push_back(plane1);
	planes.push_back(plane2);
	planes.push_back(plane3);
	planes.push_back(plane4);
	planes.push_back(plane5);

	Box box = {};
	box.position = { 0, 1, 0 };
	box.rotation = QuaternionFromAxisAngle({ 1,1,1 }, 35 * DEG2RAD);
	box.scale = { 3,.5f,3 };
	CreateBox(&box);
	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();
		Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,3,-4 }), time);
		Quaternion qOrient2 = QuaternionFromAxisAngle(Vector3Normalize({ 1,3,-4 }), time);

		MyUpdateOrbitalCamera(&camera, deltaTime);
		UpdateBall(&ball, deltaTime, planes, box);

		Sphere sphere = Sphere{};
		sphere.position = ball.position;
		sphere.rotation = { 0,0,0,1 };
		sphere.radius = ball.radius;

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			for each (Plane plane in planes)
			{
				MyDrawQuad(plane, BLUE);
			}

			MyDrawSphere(sphere, 10, 10, GREEN);

			MyDrawBox(box, RED);

			//3D REFERENTIAL
			DrawGrid(20, 1.0f);        // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);
		}
		EndMode3D();

		EndDrawing();
		//----------------------------------------------------------------------------------
	}

	// De-Initialization
	//--------------------------------------------------------------------------------------  
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}

