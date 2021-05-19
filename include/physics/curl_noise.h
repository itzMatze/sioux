#include "util.h"

void potential_occluder(
	float p[],        // center of occluder
	float radius[],   // radii of ellipsoid occluder
	float phi[],      // potential so far
	float x[], // point to evaluate
	float result[]) // resulting potential
{
	float local_x[3] = { x[0], x[1], x[2] };
	for (int k = 0; k < 3; k++) local_x[k] -= p[k];
	for (int k = 0; k < 3; k++) local_x[k] /= radius[k];
	float n[3] = { x[0], x[1], x[2] };
	normalise(n);
	float dist = length(local_x);
	float alpha = smoothstep(1.0f, 1.5f, dist);
	float dot = dot(n, phi);
	// ramp down normal component of potential field
	for (int k = 0; k < 3; k++) result[k] = (1.0f - alpha) * n[k] * dot + phi[k] * alpha;
	// clamp potential inside of occluder, to prevent moving inside it
	// this is just for visualisation, because if lines spawn inside occluder, they will move which looks wrong
	if (dist < 1.0f)
	{
		for (int k = 0; k < 3; k++) result[k] = 0.0f;
	}
}

void potential_vortex(
	float R,      // radius of influence
	float x_c[],     // center point of vortex
	float omega_c[], // angular velocity of vortex
	float x[], // point to evaluate
	float result[]) // resulting potential
{
	
	float dist[3] = { 0.0f, 0.0f, 0.0f };
	for (int k = 0; k < 3; k++) dist[k] = x[k] - x_c[k];
	// rotate with speed s
	float s = CLAMP((1.0f - length(dist) / R), 0.0f, 1.0f) * (R * R - dot(dist, dist)) / 2.0;
	for (int k = 0; k < 3; k++) result[k] = omega_c[k] * s;
}

void potential_vortex_ring(
	float R,       // radius of vortex around ring
	float r,       // radius of ring itself
	float c[],        // center of ring
	float n[],        // normal of ring
	float x[],      // point to evaluate
	float result[]) // resulting potential
{
	float dist[3] = { 0.0f, 0.0f, 0.0f };
	for (int k = 0; k < 3; k++) dist[k] = x[k] - c[k];
	float dot = dot(dist, n);
	float n_dot[3] = { n[0], n[1], n[2] };
	for (int k = 0; k < 3; k++) n_dot[k] *= dot;
	for (int k = 0; k < 3; k++) dist[k] -= n_dot[k];
	// closest point on circular curve
	float x_c[3] = { 0.0f, 0.0f, 0.0f };
	normalise(dist);
	for (int k = 0; k < 3; k++) x_c[k] = dist[k] * r;
	float omega_c[3] = { 0.0f, 0.0f, 0.0f };
	cross(n, dist, omega_c);
	for (int k = 0; k < 3; k++) x_c[k] += c[k];
	potential_vortex(R, x_c, omega_c, x, result);
}

void potential_field(float x[], float potential[], float center[], float radius)
{
	float radius_function = -abs(radius - 5.95f) + 5.95f;
	float av[] = { 0.0f, -0.5f, 0.0f }; // angular velocity
	// rotation of downwash
	float phi[3] = { 0.0f, 0.0f, 0.0f };
	potential_vortex(
		5.95f,           // radius
		center,            // center
		av, // angular velocity
		x,
		phi);
	//for (int k = 0; k < 3; k++) phi[k] = 0.0f; // removing downwash rotation for debugging
	float axis[3] = { 0.0f, 1.0f, 0.0f };
	// vortex ring of main rotor
	float temporal_potential[3] = { 0.0f, 0.0f, 0.0f };
	potential_vortex_ring(
		radius_function,        // radius of vortex around ring
		5.95f,         // radius of ring itself
		center, // center of ring
		axis, // normal of ring
		x,
		temporal_potential);
	for (int k = 0; k < 3; k++) phi[k] += temporal_potential[k];
	// second vortex ring, starts after first vortex ring covers the whole main rotor
	if (radius > 5.95f)
	{
		axis[1] = -1.0f;
		potential_vortex_ring(
			5.95f - radius_function,        // max: 11.9/4 radius of vortex around ring
			5.95f - radius_function,        // max: 11.9/4 radius of ring itself
			center, // center of ring
			axis, // normal of ring
			x,
			temporal_potential);
		for (int k = 0; k < 3; k++) phi[k] += temporal_potential[k];
	}

	// fuselage obstacle
	float com[3] = { 0.0f, 0.0f, 0.0f }; // center of mass
	float occluder_radius[3] = { 1.0f, 1.6f, 4.5f }; // dimensions of helicopter
	potential_occluder(com, occluder_radius, phi, x, temporal_potential);
	for (int k = 0; k < 3; k++) phi[k] = temporal_potential[k];
	for (int k = 0; k < 3; k++) potential[k] = phi[k];
}

// calculate partial differential equations
void potential_deriv(
	float center[],
	float radius,
	float x[],
	float dpdx[],
	float dpdy[],
	float dpdz[])
{
	float eps = 1e-4f;
	float c[3] = { 0.0f, 0.0f, 0.0f };
	float point[3] = { 0.0f, 0.0f, 0.0f };
	float delta[3] = { eps, 0.0f, 0.0f };
	float potential[3] = { 0.0f, 0.0f, 0.0f };
	potential_field(x, c, center, radius);
	for (int k = 0; k < 3; k++) point[k] = x[k] + delta[k];
	potential_field(point, potential, center, radius);
	for (int k = 0; k < 3; k++) dpdx[k] = c[k] - potential[k];
	for (int k = 0; k < 3; k++) dpdx[k] /= eps;
	delta[1] = eps;
	delta[0] = 0.0f;
	for (int k = 0; k < 3; k++) point[k] = x[k] + delta[k];
	potential_field(point, potential, center, radius);
	for (int k = 0; k < 3; k++) dpdy[k] = c[k] - potential[k];
	for (int k = 0; k < 3; k++) dpdy[k] /= eps;
	delta[2] = eps;
	delta[1] = 0.0f;
	for (int k = 0; k < 3; k++) point[k] = x[k] + delta[k];
	potential_field(point, potential, center, radius);
	for (int k = 0; k < 3; k++) dpdz[k] = c[k] - potential[k];
	for (int k = 0; k < 3; k++) dpdz[k] /= eps;
}

// compute divergence free noise by using curl (grad x)
// occluder is at position (0, 0, 0), so move center of vortex ring up to main rotor
void velocity_field(float x[], float center[], float radius, float result[])
{
	float dpdx[3] = { 0.0f, 0.0f, 0.0f };
	float dpdy[3] = { 0.0f, 0.0f, 0.0f };
	float dpdz[3] = { 0.0f, 0.0f, 0.0f };
	potential_deriv(center, radius, x, dpdx, dpdy, dpdz);
	result[0] = dpdy[2] - dpdz[1];
	result[1] = dpdz[0] - dpdx[2];
	result[2] = dpdx[1] - dpdy[0];
}
