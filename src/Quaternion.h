// Copyright (C) 2019 Martin Weigel <mail@MartinWeigel.com>
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

/**
 * @file    Quaternion.h
 * @brief   A basic quaternion library written in C
 * @date    2019-11-28
 */

#ifndef _EMBEDDEDQUATERNION_H_
#define _EMBEDDEDQUATERNION_H_

#include <math.h>
#include <stdbool.h>

#define QUATERNION_EPS (1e-4)

/*
a complex number of the form w + xi + yj + zk, where w, x, y, z are real numbers
and i, j, k are imaginary units that satisfy certain conditions.
*/
typedef struct Quaternion {
  double w;  /* real part */
  double x;  // x
  double y;  // y
  double z;  // z
} Quaternion;

// Sets the given values to the output quaternion.
void Quaternion_set(Quaternion* output, const double w, const double x,
                    const double y, const double z);

// Sets quaternion to its identity (1 + 0i + 0j + 0k).
void Quaternion_setIdentity(Quaternion* q);

// Copies one quaternion to another.
void Quaternion_copy(Quaternion* output, const Quaternion* q);

// check if all quaternion values are equal (using QUATERNION_EPS).
bool Quaternion_equal(const Quaternion* q1, const Quaternion* q2);

/*
The conjugate of a quaternion number is a quaternion with the same magnitudes
but with the sign of the imaginary parts changed
*/
void Quaternion_conjugate(Quaternion* dest, const Quaternion* source);

/* Quaternion Scalar Multiplication */
void Quaternion_scale(Quaternion* dest, const Quaternion* source,
                      const double s);

/* Quaternion Add */
void Quaternion_add(Quaternion* dest, const Quaternion* q1,
                    const Quaternion* q2);

/**
 * Calculates the squared norm of a given quaternion:
 * squared norm = w*w + v1*v1 + v2*v2 + v3*v3
 */
double Quaternion_squaredNorm(const Quaternion* q);

/**
 * Calculates the norm of a given quaternion:
 * norm = sqrt(w*w + v1*v1 + v2*v2 + v3*v3)
 */
double Quaternion_norm(const Quaternion* q);

/**
 * Normalizes the quaternion.
 */
void Quaternion_normalize(Quaternion* unit_quat, const Quaternion* source);

/**
 * convert the quaternion to a unit quaternion
 */
void Quaternion_normalize_self(Quaternion* output);

/**
 * compute the reciprocal (inverse) of a nonzero quaternion.
 */
void Quaternion_reciprocal(Quaternion* quat, const Quaternion* source);

/**
 * Multiplies two quaternions: output = q1 * q2
 * @param q1
 *      The rotation to apply on q2.
 * @param q2
 *      The orientation to be rotated.
 */
void Quaternion_multiply(Quaternion* output, const Quaternion* q1,
                         const Quaternion* q2);

/**
 * Set the unit quaternion to the equivalent of axis-angle rotation.
 * @param axis
 *      The axis (x,y,z) of the rotation (should be a unit vector).
 * @param angle
 *      Rotation angle in radians (around the axis).
 */
void Quaternion_fromAxisAngle(Quaternion* unit_quat, const double axis[3],
                              const double angle);

/**
 * convert a quaternion to axis angle (should be a unit quaternion)
 * @param output
 *      A 3D vector of the quaternion rotation axis.
 * @return
 *      The rotation angle in radians.
 */
double Quaternion_toAxisAngle(double output[3], const Quaternion* q);

/**
 * Set the quaternion to the equivalent a rotation around the X-axis.
 * @param angle
 *      Rotation angle in radians (around the x-axis).
 */
void Quaternion_fromXRotation(Quaternion* output, const double angle);

/**
 * Set the quaternion to the equivalent a rotation around the Y-axis.
 * @param angle
 *      Rotation angle in radians (around the y-axis).
 */
void Quaternion_fromYRotation(Quaternion* output, const double angle);

/**
 * Set the quaternion to the equivalent a rotation around the Z-axis.
 * @param angle
 *      Rotation angle in radians (around the z-axis).
 */
void Quaternion_fromZRotation(Quaternion* output, const double angle);

/**
 * Set the quaternion to the equivalent of euler angles.
 * @param eulerZYX
 *      Euler angles in ZYX, but stored in array as [x'', y', z].
 *      roll (x), pitch (Y), yaw (z)
 */
void Quaternion_fromEulerZYX(Quaternion* output, const double eulerZYX[3]);

/*
 * Calculates the euler angles of a unit quaternion.
 * @param output
 *      Euler angles in ZYX, but stored in array as [x'', y', z].
 */
void Quaternion_toEulerZYX(double output[3], const Quaternion* q);

/*
 * Applies quaternion rotation to a given vector.
 * Vout = q * Vin * conj(q)
 */
void Quaternion_rotate(double output[3], const Quaternion* rotation_quat,
                       const double source_v[3]);

/*
 *  Conversion rotation Matrix (DCM) to a unit Quaternion
 */
void Quaternion_rotm2quat(Quaternion* output, const double rotm[3 * 3]);

/*
Convert the quaternion to a 3x3 rotation matrix (DCM). The quaternion is
required to be normalized, otherwise the result is undefined.
*/
void Quaternion_quat2rotm(double rotm[3 * 3], const Quaternion* unit_quat);

/**
 * Interpolates between two quaternions.
 * @param t
 *      Interpolation between the two quaternions [0, 1].
 *      0 is equal with q1, 1 is equal with q2,
 *      0.5 is the middle between q1 and q2.
 */
void Quaternion_slerp(Quaternion* output, const Quaternion* q1,
                      const Quaternion* q2, const double t);

#endif /* _EMBEDDEDQUATERNION_H_ */