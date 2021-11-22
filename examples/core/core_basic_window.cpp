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
	Vector3 center;
	float radius;
};

struct Plane
{
	Vector3 center;
	Vector2 size;
};

struct Cylinder
{
	Vector3 pt1;
	Vector3 pt2;
	float radius;
};

struct Disk
{

};

struct Referential
{
	Vector3 origin;
	Vector3 i, j, k;
};

struct Capsule
{

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
	rlTranslatef(sphere.center.x, sphere.center.y, sphere.center.z);

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
	rlTranslatef(sphere.center.x, sphere.center.y, sphere.center.z);

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

void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {

}

void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {

}


void MyDrawQuad(Quaternion q, Vector3 center, Vector2 size, Color color) {

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(center.x, center.y, center.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	// Front face
	float width = size.x;
	float height = size.y;
	float length = height;
	rlVertex3f(center.x - width / 2, center.y - height / 2, center.z + length / 2);  // Bottom Left
	rlVertex3f(center.x + width / 2, center.y - height / 2, center.z + length / 2);  // Bottom Right
	rlVertex3f(center.x - width / 2, center.y + height / 2, center.z + length / 2);  // Top Left
	rlVertex3f(center.x + width / 2, center.y + height / 2, center.z + length / 2);  // Top Right
	rlVertex3f(center.x - width / 2, center.y + height / 2, center.z + length / 2);  // Top Left
	rlVertex3f(center.x + width / 2, center.y - height / 2, center.z + length / 2);  // Bottom Right

	rlEnd();
	rlPopMatrix();
}
void MyDrawQuadWire(Vector3 center, Vector2 size, Color color) {

}

void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {

}
void MyDrawCylinderWires(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {

}

void MyDrawDisk(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {

}
void MyDrawDiskWires(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {

}

void MyDrawDiskPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color) {

}

void MyDrawDiskWiresPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color) {

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
	b = 2 * (dp.x * (seg.p1.x - s.center.x) + dp.y * (seg.p1.y - s.center.y) + dp.z * (seg.p1.z - s.center.z));
	c = s.center.x * s.center.x + s.center.y * s.center.y + s.center.z * s.center.z;
	c += seg.p1.x * seg.p1.x + seg.p1.y * seg.p1.y + seg.p1.z * seg.p1.z;
	c -= 2 * (s.center.x * seg.p1.x + s.center.y * seg.p1.y + s.center.z * seg.p1.z);
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
	interNormal = Vector3Normalize(Vector3Subtract(interPt, s.center));

	return true;
}

bool InterSegmentPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal) {
	return false;
}

bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	return false;
}

bool InterSegmentFiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	return false;
}
bool InterSegmentCapsule(Segment seg, Capsule capsule, Vector3& interPt, Vector3& interNormal) {
	return false;
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

	//TEST CONVERSION CARTESIAN->CYLINDRICAL
	Vector3 pos = { 1,1,1 };
	Cylindrical cyl = CartesianToCylindrical(pos);
	printf("cyl = (%f,%f,%f) ", cyl.rho, cyl.theta, cyl.y);
	cyl = cyl + cyl;
	printf("cyl = (%f,%f,%f) ", cyl.rho, cyl.theta, cyl.y);

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


		Segment seg = Segment{ };
		seg.p1 = { 5, 5, 0 };
		seg.p2 = { -5, -2, 0 };

		Sphere sphere = Sphere{ };
		sphere.rotation = qOrient;
		sphere.center = Vector3{ 0 };
		sphere.radius = 2;

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			//
			MyDrawSphere(sphere, 40, 20, BLUE);
			MyDrawSphereWires(sphere, 40, 20, WHITE);

			DrawLine3D(seg.p1, seg.p2, RED);
			DrawSphere(seg.p1, .1f, RED);
			DrawSphere(seg.p2, .1f, RED);

			Vector3 intersection = { 0,0,0 };
			Vector3 normal = { 0,0,0 };
			InterSegmentSphere(seg, sphere, intersection, normal);
			DrawLine3D(intersection, Vector3Add(intersection, normal), RED);

			DrawSphere(intersection, .25f, GREEN);


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

