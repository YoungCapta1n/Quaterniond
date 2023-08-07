#include "Quaternion.h"

void Quaternion_set(Quaternion* output, const double w, const double x,
                    const double y, const double z) {
  output->w = w;
  output->x = x;
  output->y = y;
  output->z = z;
}  // Quaternion_set

void Quaternion_setIdentity(Quaternion* q) {
  Quaternion_set(q, 1, 0, 0, 0);
}  // Quaternion_setIdentity

void Quaternion_copy(Quaternion* dest, const Quaternion* source) {
  Quaternion_set(dest, source->w, source->x, source->y, source->z);
}  // Quaternion_copy

bool Quaternion_equal(const Quaternion* q1, const Quaternion* q2) {
  bool equalW = fabs(q1->w - q2->w) <= QUATERNION_EPS;
  bool equalX = fabs(q1->x - q2->x) <= QUATERNION_EPS;
  bool equalY = fabs(q1->y - q2->y) <= QUATERNION_EPS;
  bool equalZ = fabs(q1->z - q2->z) <= QUATERNION_EPS;
  return equalW && equalX && equalY && equalZ;
}  // Quaternion_equal

void Quaternion_conjugate(Quaternion* dest, const Quaternion* source) {
  Quaternion_set(dest, source->w, -source->x, -source->y, -source->z);
}  // Quaternion_conjugate

void Quaternion_scale(Quaternion* dest, const Quaternion* source,
                      const double s) {
  Quaternion_set(dest, s * source->w, s * source->x, s * source->y,
                 s * source->z);
}  // Quaternion_scale

void Quaternion_add(Quaternion* dest, const Quaternion* q1,
                    const Quaternion* q2) {
  Quaternion_set(dest, q1->w + q2->w, q1->x + q2->x, q1->y + q2->y,
                 q1->z + q2->z);
}  // Quaternion_add

double Quaternion_squaredNorm(const Quaternion* q) {
  return q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z;
}  // Quaternion_squaredNorm

double Quaternion_norm(const Quaternion* q) {
  return sqrt(Quaternion_squaredNorm(q));
}  // Quaternion_norm

void Quaternion_normalize(Quaternion* unit_quat, const Quaternion* q) {
  double _len = 1.0 / Quaternion_norm(q);
  Quaternion_set(unit_quat, _len * q->w, _len * q->x, _len * q->y, _len * q->z);
}  // Quaternion_normalize

void Quaternion_normalize_self(Quaternion* output) {
  double _len = 1.0 / Quaternion_norm(output);
  output->w *= _len;
  output->x *= _len;
  output->y *= _len;
  output->z *= _len;
}  // Quaternion_normalize_self

void Quaternion_reciprocal(Quaternion* reci_quat, const Quaternion* q) {
  double _len = 1.0 / Quaternion_squaredNorm(q);
  Quaternion_set(reci_quat, _len * q->w, -_len * q->x, -_len * q->y,
                 -_len * q->z);
}  // Quaternion_reciprocal

void Quaternion_multiply(Quaternion* output, const Quaternion* q1,
                         const Quaternion* q2) {
  /*
  Formula from
  http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/index.htm
           a*e - b*f - c*g - d*h
      + i (b*e + a*f + c*h- d*g)
      + j (a*g - b*h + c*e + d*f)
      + k (a*h + b*g - c*f + d*e)
  */
  output->w = q1->w * q2->w - q1->x * q2->x - q1->y * q2->y - q1->z * q2->z;
  output->x = q1->x * q2->w + q1->w * q2->x + q1->y * q2->z - q1->z * q2->y;
  output->y = q1->w * q2->y - q1->x * q2->z + q1->y * q2->w + q1->z * q2->x;
  output->z = q1->w * q2->z + q1->x * q2->y - q1->y * q2->x + q1->z * q2->w;
}  // Quaternion_multiply

void Quaternion_fromAxisAngle(Quaternion* unit_quat, const double axis[3],
                              const double angle) {
  // Formula from
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/

  double cvalue = cos(0.5 * angle);
  double svalue = sin(0.5 * angle);

  unit_quat->w = cvalue;
  unit_quat->x = svalue * axis[0];
  unit_quat->y = svalue * axis[1];
  unit_quat->z = svalue * axis[2];
}  // Quaternion_fromAxisAngle

double Quaternion_toAxisAngle(double output[3], const Quaternion* q) {
  // Formula from
  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

  Quaternion _q;
  Quaternion_copy(&_q, q);

  if (_q.w > 1) {
    // if w>1 acos and sqrt will produce errors,
    // this cant happen if quaternion is normalised
    Quaternion_normalize_self(&_q);
  }

  double angle = 2.0 * acos(_q.w);
  double divider = sqrt(1.0 - _q.w * _q.w);
  if (divider > QUATERNION_EPS) {
    // Calculate the axis
    output[0] = _q.x / divider;
    output[1] = _q.y / divider;
    output[2] = _q.z / divider;
  } else {  // test to avoid divide by zero, s is always positive due to sqrt
    // if s close to zero then direction of axis not important

    // if it is important that axis is normalised then replace with x=1; y=z=0;
    output[0] = 1;
    output[1] = 0;
    output[2] = 0;
    /*
    output[0] = _q.x;
    output[1] = _q.y;
    output[2] = _q.z;
    */
  }
  return angle;
}  // Quaternion_toAxisAngle

void Quaternion_fromXRotation(Quaternion* output, const double angle) {
  double axis[3] = {1.0, 0, 0};
  Quaternion_fromAxisAngle(output, axis, angle);
}  // Quaternion_fromXRotation

void Quaternion_fromYRotation(Quaternion* output, const double angle) {
  double axis[3] = {0, 1.0, 0};
  Quaternion_fromAxisAngle(output, axis, angle);
}  // Quaternion_fromYRotation

void Quaternion_fromZRotation(Quaternion* output, const double angle) {
  double axis[3] = {0, 0, 1.0};
  Quaternion_fromAxisAngle(output, axis, angle);
}  // Quaternion_fromZRotation

void Quaternion_fromEulerZYX(Quaternion* output, const double eulerZYX[3]) {
  // Based on
  // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

  // yaw
  double cy = cos(eulerZYX[2] * 0.5);
  double sy = sin(eulerZYX[2] * 0.5);

  // roll
  double cr = cos(eulerZYX[0] * 0.5);
  double sr = sin(eulerZYX[0] * 0.5);

  // pitch
  double cp = cos(eulerZYX[1] * 0.5);
  double sp = sin(eulerZYX[1] * 0.5);

  output->w = cy * cr * cp + sy * sr * sp;
  output->x = cy * sr * cp - sy * cr * sp;
  output->y = cy * cr * sp + sy * sr * cp;
  output->z = sy * cr * cp - cy * sr * sp;
}  // Quaternion_fromEulerZYX

void Quaternion_toEulerZYX(double output[3], const Quaternion* q) {
  // assert q be a unit quaternion

  // Roll (x-axis rotation)
  double sinr_cosp = 2.0 * (q->w * q->x + q->y * q->z);
  double cosr_cosp = 1.0 - 2.0 * (q->x * q->x + q->y * q->y);
  output[0] = atan2(sinr_cosp, cosr_cosp);

  // Pitch (y-axis rotation)
  double sinp = 2.0 * (q->w * q->y - q->z * q->x);
  if (fabs(sinp) >= 1)
    output[1] = copysign(0.5 * M_PI, sinp);  // use 90 degrees if out of range
  else
    output[1] = asin(sinp);

  // Yaw (z-axis rotation)
  double siny_cosp = 2.0 * (q->w * q->z + q->x * q->y);
  double cosy_cosp = 1.0 - 2.0 * (q->y * q->y + q->z * q->z);
  output[2] = atan2(siny_cosp, cosy_cosp);
}  // Quaternion_toEulerZYX

void Quaternion_rotate(double output[3], const Quaternion* rotation_quat,
                       const double source_v[3]) {
  double ww = rotation_quat->w * rotation_quat->w;
  double xx = rotation_quat->x * rotation_quat->x;
  double yy = rotation_quat->y * rotation_quat->y;
  double zz = rotation_quat->z * rotation_quat->z;
  double wx = rotation_quat->w * rotation_quat->x;
  double wy = rotation_quat->w * rotation_quat->y;
  double wz = rotation_quat->w * rotation_quat->z;
  double xy = rotation_quat->x * rotation_quat->y;
  double xz = rotation_quat->x * rotation_quat->z;
  double yz = rotation_quat->y * rotation_quat->z;

  // Formula from
  // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm
  // p2.x = w*w*p1.x + 2*y*w*p1.z - 2*z*w*p1.y + x*x*p1.x + 2*y*x*p1.y +
  // 2*z*x*p1.z - z*z*p1.x - y*y*p1.x; p2.y = 2*x*y*p1.x + y*y*p1.y + 2*z*y*p1.z
  // + 2*w*z*p1.x - z*z*p1.y + w*w*p1.y - 2*x*w*p1.z - x*x*p1.y; p2.z =
  // 2*x*z*p1.x + 2*y*z*p1.y + z*z*p1.z - 2*w*y*p1.x - y*y*p1.z + 2*w*x*p1.y -
  // x*x*p1.z + w*w*p1.z;
  double vx = source_v[0];
  double vy = source_v[1];
  double vz = source_v[2];

  output[0] = ww * vx + 2.0 * wy * vz - 2.0 * wz * vy + xx * vx +
              2.0 * xy * vy + 2.0 * xz * vz - zz * vx - yy * vx;
  output[1] = 2.0 * xy * vx + yy * vy + 2.0 * yz * vz + 2.0 * wz * vx -
              zz * vy + ww * vy - 2.0 * wx * vz - xx * vy;
  output[2] = 2.0 * xz * vx + 2.0 * yz * vy + zz * vz - 2.0 * wy * vx -
              yy * vz + 2.0 * wx * vy - xx * vz + ww * vz;

}  // Quaternion_rotate

void Quaternion_rotm2quat(Quaternion* output, const double rotm[3 * 3]) {
  // Formula from (Angel's Method)
  // https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

  double m00 = rotm[0];
  double m01 = rotm[1];
  double m02 = rotm[2];
  double m10 = rotm[3];
  double m11 = rotm[4];
  double m12 = rotm[5];
  double m20 = rotm[6];
  double m21 = rotm[7];
  double m22 = rotm[8];

  double trace = m00 + m11 + m22;
  if (trace > 0) {
    double s = 0.5 / sqrt(trace + 1.0);
    output->w = 0.25 / s;
    output->x = (m21 - m12) * s;
    output->y = (m02 - m20) * s;
    output->z = (m10 - m01) * s;
  } else {
    if (m00 > m11 && m00 > m22) {
      double s = 2.0 * sqrt(1.0 + m00 - m11 - m22);
      output->w = (m21 - m12) / s;
      output->x = 0.25 * s;
      output->y = (m01 + m10) / s;
      output->z = (m02 + m20) / s;
    } else if (m11 > m22) {
      double s = 2.0 * sqrt(1.0 + m11 - m00 - m22);
      output->w = (m02 - m20) / s;
      output->x = (m01 + m10) / s;
      output->y = 0.25 * s;
      output->z = (m12 + m21) / s;
    } else {
      double s = 2.0 * sqrt(1.0 + m22 - m00 - m11);
      output->w = (m10 - m01) / s;
      output->x = (m02 + m20) / s;
      output->y = (m12 + m21) / s;
      output->z = 0.25 * s;
    }
  }  // end if

  Quaternion_normalize_self(output);

}  // Quaternion_rotm2quat

void Quaternion_quat2rotm(double rotm[3 * 3], const Quaternion* unit_quat) {
  // https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm

  /*
    Quaternion _quat;
    Quaternion_normalize(&_quat, quat);

    double xx = _quat.x * _quat.x;
    double xy = _quat.x * _quat.y;
    double xz = _quat.x * _quat.z;
    double xw = _quat.x * _quat.w;

    double yy = _quat.y * _quat.y;
    double yz = _quat.y * _quat.z;
    double yw = _quat.y * _quat.w;

    double zz = _quat.z * _quat.z;
    double zw = _quat.z * _quat.w;
  */
  double xx = unit_quat->x * unit_quat->x;
  double xy = unit_quat->x * unit_quat->y;
  double xz = unit_quat->x * unit_quat->z;
  double xw = unit_quat->x * unit_quat->w;

  double yy = unit_quat->y * unit_quat->y;
  double yz = unit_quat->y * unit_quat->z;
  double yw = unit_quat->y * unit_quat->w;

  double zz = unit_quat->z * unit_quat->z;
  double zw = unit_quat->z * unit_quat->w;

  rotm[0] = 1.0 - 2.0 * (yy + zz);
  rotm[1] = 2.0 * (xy - zw);
  rotm[2] = 2.0 * (xz + yw);

  rotm[3] = 2.0 * (xy + zw);
  rotm[4] = 1.0 - 2.0 * (xx + zz);
  rotm[5] = 2.0 * (yz - xw);

  rotm[6] = 2.0 * (xz - yw);
  rotm[7] = 2.0 * (yz + xw);
  rotm[8] = 1.0 - 2.0 * (xx + yy);

}  // Quaternion_quat2rotm

void Quaternion_slerp(Quaternion* output, const Quaternion* q1,
                      const Quaternion* q2, const double t) {
  // Based on
  // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
  double cosHalfTheta =
      q1->w * q2->w + q1->x * q2->x + q1->y * q2->y + q1->z * q2->z;

  // if q1=q2 or qa=-q2 then theta = 0 and we can return qa
  if (fabs(cosHalfTheta) >= 1.0) {
    Quaternion_copy(output, q1);
    return;
  }

  double halfTheta = acos(cosHalfTheta);
  double sinHalfTheta = sqrt(1.0 - cosHalfTheta * cosHalfTheta);
  // If theta = 180 degrees then result is not fully defined
  // We could rotate around any axis normal to q1 or q2
  if (fabs(sinHalfTheta) < QUATERNION_EPS) {
    output->w = (q1->w + q2->w) * 0.5;
    output->x = (q1->x + q2->x) * 0.5;
    output->y = (q1->y + q2->y) * 0.5;
    output->z = (q1->z + q2->z) * 0.5;
    return;
  }
  // Default quaternion calculation
  double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
  double ratioB = sin(t * halfTheta) / sinHalfTheta;
  output->w = (q1->w * ratioA + q2->w * ratioB);
  output->x = (q1->x * ratioA + q2->x * ratioB);
  output->y = (q1->y * ratioA + q2->y * ratioB);
  output->z = (q1->z * ratioA + q2->z * ratioB);
  return;
}  // Quaternion_slerp