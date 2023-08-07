/*
***********************************************************************
* Quaternion_test.c: Test Quaternion
* This header file can be read by C compilers
*
*  by Hu.ZH(CrossOcean.ai)
***********************************************************************
*/
#include "../src/Quaternion.h"

#include <stdio.h>
#include <stdlib.h>

#include "minunit.h"

#define TO_RAD(x) (x / 180.0 * M_PI)

MU_TEST(testQuaternion_set) {
  Quaternion q;
  Quaternion_set(&q, 5.1, 4.2, 3.3, 2.4);
  mu_assert_double_eq(q.w, 5.1);
  mu_assert_double_eq(q.x, 4.2);
  mu_assert_double_eq(q.y, 3.3);
  mu_assert_double_eq(q.z, 2.4);

  Quaternion q1 = {
      -0.2,  // w
      1.2,   // x
      3.9,   // y
      7.3    // z
  };
  mu_assert_double_eq(q1.w, -0.2);
  mu_assert_double_eq(q1.x, 1.2);
  mu_assert_double_eq(q1.y, 3.9);
  mu_assert_double_eq(q1.z, 7.3);
}

MU_TEST(testQuaternion_setIdentity) {
  Quaternion q;
  Quaternion_set(&q, 5.1, 4.2, 3.3, 2.4);
  Quaternion_setIdentity(&q);
  mu_assert_double_eq(q.w, 1.0);
  mu_assert_double_eq(q.x, 0.0);
  mu_assert_double_eq(q.y, 0.0);
  mu_assert_double_eq(q.z, 0.0);
}

MU_TEST(testQuaternion_copy) {
  Quaternion q, r;
  Quaternion_set(&q, 5.1, 4.2, 3.3, 2.4);
  Quaternion_copy(&r, &q);
  mu_assert_double_eq(r.w, q.w);
  mu_assert_double_eq(r.x, q.x);
  mu_assert_double_eq(r.y, q.y);
  mu_assert_double_eq(r.z, q.z);
}

MU_TEST(testQuaternion_conjugate) {
  Quaternion q, c;
  Quaternion_set(&q, 5.1, 4.2, 3.3, 2.4);
  Quaternion_conjugate(&c, &q);
  mu_assert_double_eq(c.w, q.w);
  mu_assert_double_eq(c.x, -q.x);
  mu_assert_double_eq(c.y, -q.y);
  mu_assert_double_eq(c.z, -q.z);
}

MU_TEST(testQuaternion_equal) {
  Quaternion q1, q2, q3, q4, q5;
  Quaternion_set(&q1, 5.1, 4.2, 3.3, 2.4);
  Quaternion_set(&q2, 5.1, 4.2, 3.3, 2.4);
  Quaternion_set(&q3, 5.1, 4.2 + (QUATERNION_EPS / 2.0), 3.3, 2.4);
  Quaternion_set(&q4, 9.1, 4.2, 3.3, 2.4);
  Quaternion_set(&q5, 5.1, 4.2, 9.3, 2.4);
  mu_assert(Quaternion_equal(&q1, &q1), "Quaternion_equal with same reference");
  mu_assert(Quaternion_equal(&q1, &q2), "Quaternion_equal with same values");
  mu_assert(Quaternion_equal(&q1, &q3), "Quaternion_equal with similar values");
  mu_assert(!Quaternion_equal(&q1, &q4), "Quaternion_equal with different w");
  mu_assert(!Quaternion_equal(&q1, &q5), "Quaternion_equal with different v");
}

MU_TEST(testQuaternion_addscale) {
  Quaternion q, c;
  Quaternion_set(&q, 5.1, -4.2, -3.3, 2.4);
  Quaternion_scale(&c, &q, 2.0);
  mu_assert_double_eq(c.w, 10.2);
  mu_assert_double_eq(c.x, -8.4);
  mu_assert_double_eq(c.y, -6.6);
  mu_assert_double_eq(c.z, 4.8);

  Quaternion sq;
  Quaternion_add(&sq, &q, &c);

  mu_assert_double_eq(sq.w, 15.3);
  mu_assert_double_eq(sq.x, -12.6);
  mu_assert_double_eq(sq.y, -9.9);
  mu_assert_double_eq(sq.z, 7.2);
}

MU_TEST(testQuaternion_norm) {
  Quaternion q, c;
  Quaternion_set(&q, 0.5 * sqrt(2), 0.5 * sqrt(2), 0, 0);
  Quaternion_setIdentity(&c);
  mu_assert_double_eq(Quaternion_norm(&q), 1);
  mu_assert_double_eq(Quaternion_norm(&c), 1);

  Quaternion_set(&q, 4, 5, 6, 7);
  mu_assert_double_eq(Quaternion_norm(&q), sqrt(126));
}

MU_TEST(testQuaternion_normalize) {
  Quaternion q, c;
  Quaternion_set(&q, -4, 5, -6, 7);
  Quaternion_normalize(&c, &q);

  mu_assert_double_eq(c.w, -4 / sqrt(126));
  mu_assert_double_eq(c.x, 5 / sqrt(126));
  mu_assert_double_eq(c.y, -6 / sqrt(126));
  mu_assert_double_eq(c.z, 7 / sqrt(126));

  Quaternion_normalize_self(&q);
  mu_assert_double_eq(q.w, -4 / sqrt(126));
  mu_assert_double_eq(q.x, 5 / sqrt(126));
  mu_assert_double_eq(q.y, -6 / sqrt(126));
  mu_assert_double_eq(q.z, 7 / sqrt(126));
}

MU_TEST(testQuaternion_multiply) {
  // Examples from https://www.mathworks.com/help/aerotbx/ug/quatmultiply.html
  Quaternion q1, q2, result, real;
  Quaternion_set(&q1, 1.0, 0.0, 1.0, 0.0);
  Quaternion_set(&q2, 1.0, 0.5, 0.5, 0.75);
  Quaternion_set(&real, 0.5, 1.25, 1.5, 0.25);
  Quaternion_multiply(&result, &q1, &q2);
  mu_assert(Quaternion_equal(&result, &real), "Quaternion_multiply example 1");

  Quaternion_set(&q1, 1.0, 0.0, 1.0, 0.0);
  Quaternion_set(&real, 0, 0, 2, 0);
  Quaternion_multiply(&result, &q1, &q1);
  mu_assert(Quaternion_equal(&result, &real), "Quaternion_multiply example 2");

  Quaternion_set(&q1, 1.0, 0.0, 1.0, 0.0);
  Quaternion_set(&q2, 2.0, 1.0, 0.1, 0.1);
  Quaternion_set(&real, 1.9, 1.1, 2.1, -0.9);
  Quaternion_multiply(&result, &q1, &q2);
  mu_assert(Quaternion_equal(&result, &real), "Quaternion_multiply example 3");
}

MU_TEST(testQuaternion_fromAxisAngle) {
  Quaternion q;
  double vector[3] = {1, 0, 0};
  double angle = TO_RAD(90.0);
  Quaternion_fromAxisAngle(&q, vector, angle);
  mu_assert_double_eq(q.w, 0.5 * sqrt(2));
  mu_assert_double_eq(q.x, 0.5 * sqrt(2));
  mu_assert_double_eq(q.y, 0);
  mu_assert_double_eq(q.z, 0);
  mu_assert_double_eq(Quaternion_norm(&q), 1.0);
}

MU_TEST(testQuaternion_fromXRotation) {
  Quaternion q;
  double vector[3] = {1, 0, 0};
  double angle = TO_RAD(90.0);
  Quaternion_fromAxisAngle(&q, vector, angle);
  mu_assert_double_eq(q.w, 0.5 * sqrt(2));
  mu_assert_double_eq(q.x, 0.5 * sqrt(2));
  mu_assert_double_eq(q.y, 0);
  mu_assert_double_eq(q.z, 0);

  Quaternion c;
  Quaternion_fromXRotation(&c, angle);
  mu_assert_double_eq(q.w, c.w);
  mu_assert_double_eq(q.x, c.x);
  mu_assert_double_eq(q.y, c.y);
  mu_assert_double_eq(q.z, c.z);
}

MU_TEST(testQuaternion_fromYRotation) {
  Quaternion q;
  double vector[3] = {0, 1, 0};
  double angle = TO_RAD(90.0);
  Quaternion_fromAxisAngle(&q, vector, angle);
  mu_assert_double_eq(q.w, 0.5 * sqrt(2));
  mu_assert_double_eq(q.x, 0);
  mu_assert_double_eq(q.y, 0.5 * sqrt(2));
  mu_assert_double_eq(q.z, 0);
}

MU_TEST(testQuaternion_fromZRotation) {
  Quaternion q;
  double vector[3] = {0, 0, 1};
  double angle = TO_RAD(90.0);
  Quaternion_fromAxisAngle(&q, vector, angle);
  mu_assert_double_eq(q.w, 0.5 * sqrt(2));
  mu_assert_double_eq(q.x, 0);
  mu_assert_double_eq(q.y, 0);
  mu_assert_double_eq(q.z, 0.5 * sqrt(2));

  Quaternion c;
  Quaternion_fromZRotation(&c, -0.5 * M_PI);
  mu_assert_double_eq(c.w, 0.5 * sqrt(2));
  mu_assert_double_eq(c.x, 0);
  mu_assert_double_eq(c.y, 0);
  mu_assert_double_eq(c.z, -0.5 * sqrt(2));
}

MU_TEST(testQuaternion_toAxisAngle) {
  double v90[3];
  Quaternion rot90;
  Quaternion_set(&rot90, 1.5 * sqrt(2), 1.5 * sqrt(2), 0, 0);
  double a90 = Quaternion_toAxisAngle(v90, &rot90);
  mu_assert_double_eq(a90, TO_RAD(90.0));
  mu_assert_double_eq(v90[0], 1);
  mu_assert_double_eq(v90[1], 0);
  mu_assert_double_eq(v90[2], 0);

  double v0[3];
  Quaternion rot0;
  Quaternion_setIdentity(&rot0);
  double a0 = Quaternion_toAxisAngle(v0, &rot0);
  mu_assert_double_eq(a0, 0);

  double arbAxis[3] = {0.802, 0.267, -0.534};
  double vArb[3];
  Quaternion rotArb;
  Quaternion_fromAxisAngle(&rotArb, arbAxis, 1.11);
  double aArb = Quaternion_toAxisAngle(vArb, &rotArb);
  mu_assert_double_eq(aArb, 1.11);
  mu_assert_double_eq(vArb[0], arbAxis[0]);
  mu_assert_double_eq(vArb[1], arbAxis[1]);
  mu_assert_double_eq(vArb[2], arbAxis[2]);
}

MU_TEST(testQuaternion_fromEulerZYX) {
  // Check test with https://www.andre-gaschler.com/rotationconverter/
  double e0[3] = {TO_RAD(90.0), 0, 0};
  Quaternion q, real;
  Quaternion_set(&real, 0.7071, 0.7071, 0, 0);
  Quaternion_fromEulerZYX(&q, e0);
  mu_assert(Quaternion_equal(&q, &real),
            "Quaternion_fromEulerZYX for x-rotation only");

  double e1[3] = {TO_RAD(90.0), 0, TO_RAD(180.0)};
  Quaternion_set(&real, 0, 0, 0.7071, 0.7071);
  Quaternion_fromEulerZYX(&q, e1);
  mu_assert(Quaternion_equal(&q, &real), "Quaternion_fromEulerZYX example 1");

  double e2[3] = {TO_RAD(90.0), 0, TO_RAD(90.0)};
  Quaternion_set(&real, 0.5, 0.5, 0.5, 0.5);
  Quaternion_fromEulerZYX(&q, e2);
  mu_assert(Quaternion_equal(&q, &real), "Quaternion_fromEulerZYX example 2");

  double e3[3] = {TO_RAD(165.0), TO_RAD(63.0), TO_RAD(122.0)};
  Quaternion_set(&real, 0.5070333, 0.3501829, 0.7724199, -0.1538071);
  Quaternion_fromEulerZYX(&q, e3);
  mu_assert(Quaternion_equal(&q, &real), "Quaternion_fromEulerZYX example 3");
}

MU_TEST(testQuaternion_toEulerZYX) {
  // Check test with https://www.andre-gaschler.com/rotationconverter/
  double euler[3];
  Quaternion q;
  Quaternion_set(&q, 0.5 * sqrt(2), 0.5 * sqrt(2), 0, 0);
  Quaternion_toEulerZYX(euler, &q);
  mu_assert_double_eq(euler[0], TO_RAD(90.0));
  mu_assert_double_eq(euler[1], 0);
  mu_assert_double_eq(euler[2], 0);

  Quaternion_set(&q, 0, 0, 0.5 * sqrt(2), 0.5 * sqrt(2));
  Quaternion_toEulerZYX(euler, &q);
  mu_assert_double_eq(euler[0], TO_RAD(90.0));
  mu_assert_double_eq(euler[1], 0);
  mu_assert_double_eq(euler[2], TO_RAD(180.0));

  Quaternion_set(&q, 0.5, 0.5, 0.5, 0.5);
  Quaternion_toEulerZYX(euler, &q);
  mu_assert_double_eq(euler[0], TO_RAD(90.0));
  mu_assert_double_eq(euler[1], 0);
  mu_assert_double_eq(euler[2], TO_RAD(90.0));

  Quaternion_set(&q, 0.5070333, 0.3501829, 0.7724199, -0.1538071);
  Quaternion_toEulerZYX(euler, &q);
  mu_assert_double_eq(euler[0], 2.879793116146);
  mu_assert_double_eq(euler[1], 1.09955727504);
  mu_assert_double_eq(euler[2], 2.129301506661);
}

MU_TEST(testQuaternion_rotate) {
  double result[3];
  double v1[3] = {5.1, 6.8, -5.3};
  Quaternion identity;
  Quaternion_setIdentity(&identity);
  Quaternion_rotate(result, &identity, v1);
  mu_assert_double_eq(result[0], v1[0]);
  mu_assert_double_eq(result[1], v1[1]);
  mu_assert_double_eq(result[2], v1[2]);

  // Example 2 from http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf
  Quaternion q;
  double v2[3] = {1, 0, 0};
  Quaternion_set(&q, 0.5, 0.5, 0.5, 0.5);
  Quaternion_rotate(result, &q, v2);
  mu_assert_double_eq(result[0], 0);
  mu_assert_double_eq(result[1], 1);
  mu_assert_double_eq(result[2], 0);

  // Example from
  // https://gamedev.stackexchange.com/questions/119725/why-does-transforming-a-vector3-by-a-quaternion-result-in-reversed-z
  double v3[3] = {1, 0, 0};
  Quaternion_set(&q, 0.6532815, -0.270598, 0.270598, 0.6532815);
  Quaternion_rotate(result, &q, v3);
  mu_assert_double_eq(result[0], 0);
  mu_assert_double_eq(result[1], 0.7071068812765);
  mu_assert_double_eq(result[2], -0.707106669348);

  // Example from
  // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm
  double v4[3] = {1, 0, 0};
  Quaternion_set(&q, 0.5 * sqrt(2), 0, 0, 0.5 * sqrt(2));
  Quaternion_rotate(result, &q, v4);
  mu_assert_double_eq(result[0], 0);
  mu_assert_double_eq(result[1], 1);
  mu_assert_double_eq(result[2], 0);
}

MU_TEST(testQuaternion_slerp) {
  Quaternion q1, q2, result;
  Quaternion_set(&q1, 0.6532815, -0.270598, 0.270598, 0.6532815);
  Quaternion_set(&q2, 0.5, 0.5, 0.5, 0.5);

  Quaternion_slerp(&result, &q1, &q2, 0);
  mu_assert(Quaternion_equal(&result, &q1), "Quaternion_slerp with t=0");

  Quaternion_slerp(&result, &q1, &q2, 1);
  mu_assert(Quaternion_equal(&result, &q2), "Quaternion_slerp with t=1");

  Quaternion_slerp(&result, &q1, &q2, 0.62);
  mu_assert_double_eq(result.w, 0.6119266025696755);
  mu_assert_double_eq(result.x, 0.22069444274723088);
  mu_assert_double_eq(result.y, 0.4498729015909088);
  mu_assert_double_eq(result.z, 0.6119266025696755);
}

MU_TEST(testQuaternion_rotm2quat) {
  // Example from
  // https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/examples/index.htm

  {
    double rotm[3 * 3] = {1, 0, 0,   //
                          0, 0, -1,  //
                          0, 1, 0};

    Quaternion q1;
    Quaternion_rotm2quat(&q1, rotm);

    mu_assert_double_eq(q1.w, 0.5 * sqrt(2));
    mu_assert_double_eq(q1.x, 0.5 * sqrt(2));
    mu_assert_double_eq(q1.y, 0.0);
    mu_assert_double_eq(q1.z, 0.0);
  }

  {
    double rotm[3 * 3] = {0,  0, 1,  //
                          0,  1, 0,  //
                          -1, 0, 0};
    Quaternion q1;
    Quaternion_rotm2quat(&q1, rotm);

    mu_assert_double_eq(q1.w, 0.5 * sqrt(2));
    mu_assert_double_eq(q1.x, 0.0);
    mu_assert_double_eq(q1.y, 0.5 * sqrt(2));
    mu_assert_double_eq(q1.z, 0.0);
  }

  {
    double rotm[3 * 3] = {-1, 0, 0,  //
                          0,  1, 0,  //
                          0,  0, -1};
    Quaternion q1;
    Quaternion_rotm2quat(&q1, rotm);

    mu_assert_double_eq(q1.w, 0.0);
    mu_assert_double_eq(q1.x, 0.0);
    mu_assert_double_eq(q1.y, 1.0);
    mu_assert_double_eq(q1.z, 0.0);
  }

  {
    double rotm[3 * 3] = {0, 0,  -1,  //
                          1, 0,  0,   //
                          0, -1, 0};
    Quaternion q1;
    Quaternion_rotm2quat(&q1, rotm);

    mu_assert_double_eq(q1.w, 0.5);
    mu_assert_double_eq(q1.x, -0.5);
    mu_assert_double_eq(q1.y, -0.5);
    mu_assert_double_eq(q1.z, 0.5);
  }
}

MU_TEST(testQuaternion_quat2rotm) {
  {
    double rotm[3 * 3];
    Quaternion q = {
        1.0,  // w
        0.0,  // x
        0.0,  // y
        0.0   // z
    };
    Quaternion_quat2rotm(rotm, &q);

    double correct_rotm[3 * 3] = {
        1.0, 0.0, 0.0,  //
        0.0, 1.0, 0.0,  //
        0.0, 0.0, 1.0   //
    };
    for (int i = 0; i != 9; ++i) mu_assert_double_eq(rotm[i], correct_rotm[i]);
  }

  {
    double rotm[3 * 3];
    Quaternion q = {
        0.5 * sqrt(2),  // w
        0.0,            // x
        0.5 * sqrt(2),  // y
        0.0             // z
    };
    Quaternion_quat2rotm(rotm, &q);

    double correct_rotm[3 * 3] = {
        0.0,  0.0, 1.0,  //
        0.0,  1.0, 0.0,  //
        -1.0, 0.0, 0.0   //
    };
    for (int i = 0; i != 9; ++i) mu_assert_double_eq(rotm[i], correct_rotm[i]);
  }
}

void test_setup(void) { /* Nothing */
}

void test_teardown(void) { /* Nothing */
}

MU_TEST_SUITE(test_suite) {
  MU_SUITE_CONFIGURE(&test_setup, &test_teardown);

  MU_RUN_TEST(testQuaternion_set);
  MU_RUN_TEST(testQuaternion_setIdentity);
  MU_RUN_TEST(testQuaternion_copy);
  MU_RUN_TEST(testQuaternion_conjugate);
  MU_RUN_TEST(testQuaternion_equal);
  MU_RUN_TEST(testQuaternion_addscale);
  MU_RUN_TEST(testQuaternion_norm);
  MU_RUN_TEST(testQuaternion_normalize);
  MU_RUN_TEST(testQuaternion_multiply);
  MU_RUN_TEST(testQuaternion_fromAxisAngle);
  MU_RUN_TEST(testQuaternion_fromXRotation);
  MU_RUN_TEST(testQuaternion_fromYRotation);
  MU_RUN_TEST(testQuaternion_fromZRotation);
  MU_RUN_TEST(testQuaternion_toAxisAngle);
  MU_RUN_TEST(testQuaternion_fromEulerZYX);
  MU_RUN_TEST(testQuaternion_toEulerZYX);
  MU_RUN_TEST(testQuaternion_rotate);
  MU_RUN_TEST(testQuaternion_slerp);
  MU_RUN_TEST(testQuaternion_rotm2quat);
  MU_RUN_TEST(testQuaternion_quat2rotm);
}

int main(void) {
  MU_RUN_SUITE(test_suite);
  MU_REPORT();

  return EXIT_SUCCESS;
}