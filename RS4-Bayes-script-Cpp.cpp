
/*
 * Preamble: 
 * 
 * Title: C++ functions for Hidden-state Bayesian simulation of online motor responses to sonification feedback.
 * 
 */

#include <Rcpp.h>
using namespace Rcpp;

/*
 * Sonification constants
 */

const int base_wave_freq = 152380; // highest available with 9 bit resolution is 152380
const int pwm_resolution = (int) log2( 80000000 / base_wave_freq );
const int amp_max = pow(2,pwm_resolution) - 1;
const int pitchinitial = 440;
const int PitchValueMax = 6000;
const int PitchValueMin = 220;

// [[Rcpp::export]]
int get_amp_max() {
  return amp_max;
}
// [[Rcpp::export]]
int get_PitchValueMax() {
  return PitchValueMax;
}

/*
 * Helper Functions
 */

// Vector-Valued Functions

// Function #1: For computing the distance between two 4D vectors (quaternions)
float vdst4 ( NumericVector thisVector1, NumericVector thisVector2 ) {
  float x = thisVector1[0] - thisVector2[0];
  float y = thisVector1[1] - thisVector2[1];
  float z = thisVector1[2] - thisVector2[2];
  float r = thisVector1[3] - thisVector2[3];
  float distance = sqrtf ( x*x + y*y + z*z + r*r );
  return distance;
}

// Function #2: For computing the squared distance between two 4D vectors (quaternions)
float vdstSq4 ( NumericVector thisVector1, NumericVector thisVector2 ) {
  float x = thisVector1[0] - thisVector2[0];
  float y = thisVector1[1] - thisVector2[1];
  float z = thisVector1[2] - thisVector2[2];
  float r = thisVector1[3] - thisVector2[3];
  float distance = x*x + y*y + z*z + r*r;
  return distance;
}

// Function #3: For subtracting two 3D vectors 
NumericVector vsubtract ( NumericVector v1, NumericVector v2 ) {
  int input_length = v1.size();
  int input_length2 = v2.size();
  if ( input_length != input_length2 ) stop("Error: Vectors must be of the same length.");
  NumericVector output(input_length);
  output[0] = v1[0] - v2[0];
  output[1] = v1[1] - v2[1];
  output[2] = v1[2] - v2[2];
  if ( input_length == 4 ) output[3] = v1[3] - v2[3];
  return output;
}

// Function #4: For adding two vectors
NumericVector vadd ( NumericVector v1, NumericVector v2 ) {
  int input_length = v1.size();
  int input_length2 = v2.size();
  if ( input_length != input_length2 ) stop("Error: Vectors must be of the same length.");
  NumericVector output(input_length);
  output[0] = v1[0] + v2[0];
  output[1] = v1[1] + v2[1];
  output[2] = v1[2] + v2[2];
  if ( input_length == 4 ) output[3] = v1[3] + v2[3];
  return output;
}

// Function #5: For multiplying a vector by a scalar
NumericVector vscaler_mult ( NumericVector thisVector1, float thisfloat ) {
  int input_length = thisVector1.size();
  NumericVector output(input_length);
  output[0] = thisVector1[0] * thisfloat;
  output[1] = thisVector1[1] * thisfloat;
  output[2] = thisVector1[2] * thisfloat;
  if ( input_length == 4 ) output[3] = thisVector1[3] * thisfloat;
  return output;
}

// Function #6: For computing the cross product of two vectors
NumericVector crossp ( NumericVector thisVector1, NumericVector thisVector2 ) {
  // For 4D quaternions, returns cross product of the real x-y-z components only. 
  int input_length = thisVector1.size();
  int input_length2 = thisVector2.size();
  if ( input_length != input_length2 ) stop("Error: Vectors must be of the same length.");
  NumericVector crossproduct(input_length);
  crossproduct[0] = thisVector1[1] * thisVector2[2] - thisVector1[2] * thisVector2[1];
  crossproduct[1] = thisVector1[2] * thisVector2[0] - thisVector1[0] * thisVector2[2];
  crossproduct[2] = thisVector1[0] * thisVector2[1] - thisVector1[1] * thisVector2[0];
  if ( input_length == 4 ) crossproduct[3] = 0;
  return crossproduct;
}

// Function #7: For computing the dot product of two vectors
float dotp ( NumericVector thisVector1, NumericVector thisVector2 ) {
  int input_length = thisVector1.size();
  int input_length2 = thisVector2.size();
  if ( input_length != input_length2 ) stop("Error: Vectors must be of the same length.");
  float dotproduct = thisVector1[0] * thisVector2[0] + thisVector1[1] * thisVector2[1] + thisVector1[2] * thisVector2[2];
  // Note: Always only takes the dot product of the first 3 components of the vectors.
  return dotproduct;
}

// Function #8: For computing the magnitude of a vector
float vmag ( NumericVector thisVector ) {
  int input_length = thisVector.size();
  float MagV;
  if (input_length == 3) MagV = sqrtf ( thisVector[0] * thisVector[0] + thisVector[1] * thisVector[1] + thisVector[2] * thisVector[2] );
  else MagV = sqrtf ( thisVector[0] * thisVector[0] + thisVector[1] * thisVector[1] + thisVector[2] * thisVector[2] + thisVector[3] * thisVector[3] );
  return MagV;
}

// Function #9: For dividing a vector by a scalar
NumericVector vdivide ( NumericVector v, float r ) {
  int input_length = v.size();
  NumericVector w(input_length);
  w[0] = v[0] /r;
  w[1] = v[1] /r;
  w[2] = v[2] /r;
  if (input_length == 4) w[3] = v[3] /r;
  return w;
}

// Function #10: For normalizing a vector
NumericVector normalizeV ( NumericVector v ) {
  int input_length = v.size();
  NumericVector w(input_length);
  w = vdivide(v,vmag(v));
  return w; 
}

// Quaternion Functions

// Function #0: Function for converting to Euler angles in 3-2-1 sequence
NumericVector ToEulerAngles( NumericVector q) {
  
  NumericVector euler(3);
  q = normalizeV(q);
  
  // roll (x-axis rotation)
  double sinr_cosp = 2 * (q[3] * q[0] + q[1] * q[2]);
  double cosr_cosp = 1 - 2 * (q[0] * q[0] + q[1] * q[1]);
  euler[0] = std::atan2(sinr_cosp, cosr_cosp);
  
  // pitch (y-axis rotation)
  double sinp = std::sqrt(1 + 2 * (q[3] * q[1] - q[0] * q[2]));
  double cosp = std::sqrt(1 - 2 * (q[3] * q[1] - q[0] * q[2]));
  euler[1] = 2 * std::atan2(sinp, cosp) - M_PI / 2;
  
  // yaw (z-axis rotation)
  double siny_cosp = 2 * (q[3] * q[2] + q[0] * q[1]);
  double cosy_cosp = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
  euler[2] = std::atan2(siny_cosp, cosy_cosp);

  return euler;
  
}

// Function #1: For forming a quaternion from an angle theta and a rotation axis v,
//  ... this quaterion rotates space by angle theta (radians) about axis v
NumericVector formquat ( float theta, NumericVector v ) {
  NumericVector q(4);
  float cos2theta = cosf(theta/2.0);
  float sin2theta = sinf(theta/2.0);
  q[3] = cos2theta;
  q[0] = sin2theta * v[0];
  q[1] = sin2theta * v[1];
  q[2] = sin2theta * v[2];
  q = normalizeV(q);
  return q; 
}

// Function #2: For multiplying two quaternions
NumericVector quat_mult ( NumericVector q1, NumericVector q2 ) {
  NumericVector q(4);
  NumericVector qv(4); 
  q1 = normalizeV(q1);
  q2 = normalizeV(q2);
  qv = vadd( vscaler_mult(q2,q1[3]) , vadd( vscaler_mult(q1,q2[3]) , crossp(q1,q2) ) );
  q[3] = q1[3] * q2[3] - dotp(q1,q2); 
  q[0] = qv[0];
  q[1] = qv[1];
  q[2] = qv[2];
  q = normalizeV(q);
  return q;
}

// Function #3: For formating a quat from Euler rotations
NumericVector rot_quat ( NumericVector r, NumericVector a1, NumericVector a2, NumericVector a3, NumericVector u ) {
  
  // This code was written awhile back and assumes r is a set of Euler angles
  r = ToEulerAngles(r);
  
  // Shift coordinates so that u is the origin
  NumericVector a1shift = vsubtract(a1,u); 
  NumericVector a2shift = vsubtract(a2,u); 
  NumericVector a3shift = vsubtract(a3,u);
  
  // Form new quat ... 
  NumericVector qx = formquat(r[0],a1shift);
  NumericVector qy = formquat(r[1],a2shift);
  NumericVector qz = formquat(r[2],a3shift);
  NumericVector q = quat_mult( qx, quat_mult( qy, qz) );
  
  q = normalizeV(q); // we must ensure this is normalized, although it should be.
  return q;
  
}

// Function #4: For rotating a vector by a quaternion
NumericVector qvq ( NumericVector thisQuat, NumericVector thisPosition ) {
  
  // Will rotate the vector thisPosition by the quaternion thisQuat
  // Fewer computations than the other definition
  // Also, note that we're assuming our quaternions are unit quaternions (they should be, or should be very close). 
  
  thisQuat = normalizeV(thisQuat);
  
  int input_length = thisPosition.size();
  int input_length2 = thisQuat.size();
  if ( input_length2 != 4 ) stop("Error: qvq requires 4D quat.");
  if ( input_length != 4 ) {
    NumericVector new_thisPosition(4);
    new_thisPosition[0] = thisPosition[0];
    new_thisPosition[1] = thisPosition[1];
    new_thisPosition[2] = thisPosition[2];
    new_thisPosition[3] = 0;
    thisPosition = new_thisPosition;
  }
  
  thisQuat = normalizeV(thisQuat);
  
  float qv_r = thisQuat[3] * thisPosition[3] - thisQuat[0] * thisPosition[0] - thisQuat[1] * thisPosition[1] - thisQuat[2] * thisPosition[2];
  float qv_x = thisQuat[3] * thisPosition[0] + thisQuat[0] * thisPosition[3] + thisQuat[1] * thisPosition[2] - thisQuat[2] * thisPosition[1]; 
  float qv_y = thisQuat[3] * thisPosition[1] - thisQuat[0] * thisPosition[2] + thisQuat[1] * thisPosition[3] + thisQuat[2] * thisPosition[0];
  float qv_z = thisQuat[3] * thisPosition[2] + thisQuat[0] * thisPosition[1] - thisQuat[1] * thisPosition[0] + thisQuat[2] * thisPosition[3]; 

  // This is the same operation as above, expect the sign is flipped everywhere we have a non-real quat component (since the second q in qvq is the inverse of q). 
  float qvq_r = qv_r * thisQuat[3] + qv_x * thisQuat[0] + qv_y * thisQuat[1] + qv_z * thisQuat[2]; 
  float qvq_x = qv_r * thisQuat[0] * -1.0 + qv_x * thisQuat[3] - qv_y * thisQuat[2] + qv_z * thisQuat[1];
  float qvq_y = qv_r * thisQuat[1] * -1.0 + qv_x * thisQuat[2] + qv_y * thisQuat[3] - qv_z * thisQuat[0];
  float qvq_z = qv_r * thisQuat[2] * -1.0 - qv_x * thisQuat[1] + qv_y * thisQuat[0] + qv_z * thisQuat[3];

  NumericVector thisFinalPosition(input_length);
  thisFinalPosition[0] = qvq_x;
  thisFinalPosition[1] = qvq_y;
  thisFinalPosition[2] = qvq_z;
  if ( input_length==4 ) thisFinalPosition[3] = qvq_r; // physically meaningful vectors have r = 0

  return thisFinalPosition;
  
}

// Function #5: For rotating a vector by a quaternion
NumericVector quatrot ( NumericVector q, NumericVector v, NumericVector p) {
  q = normalizeV(q);
  NumericVector vshift = vsubtract(v,p);
  NumericVector rotatedposition = qvq(q,vshift);
  rotatedposition = vadd(rotatedposition,p);
  return rotatedposition;
}

// Function #6: Compute quaternion conjugate
NumericVector quat_conj ( NumericVector q ) {
  q = normalizeV(q); 
  NumericVector qc(4);
  qc[3] = q[3];
  qc[0] = -1 * q[0];
  qc[1] = -1 * q[1];
  qc[2] = -1 * q[2];
  qc = normalizeV(qc);
  return(qc); 
}

// Function #7: Find quat difference
// [[Rcpp::export]]
NumericVector quat_diff ( NumericVector q_old, NumericVector q_new ) {
  NumericVector qd(4);
  qd = quat_mult(quat_conj(q_old),q_new);
  qd = normalizeV(qd);
  return(qd); 
}

// Function #8: Jitter quaternion trajecotry
// [[Rcpp::export]]
NumericMatrix jitter_quat_trajectory ( NumericMatrix q1_traj, NumericMatrix q2_traj, NumericVector jitter1, NumericVector jitter2 ) {
  
  int num_samples = q1_traj.nrow();
  NumericMatrix q_traj_jittered(num_samples,8);
  
  float deg1 = jitter1[0] * 5.0;
  NumericVector jitter_vec1(3);
  jitter_vec1[0] = jitter1[1];
  jitter_vec1[1] = jitter1[2];
  jitter_vec1[2] = jitter1[3];
  jitter_vec1 = normalizeV(jitter_vec1);
  NumericVector jitter_quat1 = formquat(deg1*(3.1416/180.0), jitter_vec1);
  
  float deg2 = jitter2[0] * 5.0;
  NumericVector jitter_vec2(3);
  jitter_vec2[0] = jitter2[1];
  jitter_vec2[1] = jitter2[2];
  jitter_vec2[2] = jitter2[3];
  jitter_vec2 = normalizeV(jitter_vec2);
  NumericVector jitter_quat2 = formquat(deg2*(3.1416/180.0), jitter_vec2);
  
  for ( int s = 0; s < num_samples; ++s ) {
    NumericVector q1_jittered = normalizeV(qvq(jitter_quat1,q1_traj(s,_)));
    NumericVector q2_jittered = normalizeV(qvq(jitter_quat2,q2_traj(s,_)));
    q_traj_jittered(s,0) = q1_jittered[0];
    q_traj_jittered(s,1) = q1_jittered[1];
    q_traj_jittered(s,2) = q1_jittered[2];
    q_traj_jittered(s,3) = q1_jittered[3];
    q_traj_jittered(s,4) = q2_jittered[0];
    q_traj_jittered(s,5) = q2_jittered[1];
    q_traj_jittered(s,6) = q2_jittered[2];
    q_traj_jittered(s,7) = q2_jittered[3];
  }
  
  return q_traj_jittered;

}

// Misc: 

// Function #1: For keeping a value constrained (similar to Arduino library constrain)
int constrain(int value, int minVal, int maxVal) {
  int output = value;
  if (output < minVal) output = minVal;
  if (output > maxVal) output = maxVal;
  return output;
}

/*
 * Sonification simulation
 * 
 * The below function uses the code from the sonification system to simulate the sonification response to
 *    a reach. It includes the real-time time-warping algorithm.
 */

// [[Rcpp::export]]
NumericMatrix simulate_sonification(
    NumericMatrix Sensor1_Quat, 
    NumericMatrix Sensor2_Quat,
    NumericMatrix model_q1,
    NumericMatrix model_q2
) {
  
  float to_end_start = vdst4(model_q1(0, _),model_q1(model_q1.nrow()-1, _)) + vdst4(model_q2(0, _),model_q2(model_q2.nrow()-1, _));
  float a_amp = to_end_start / ( -1.0 * logf( 1.0 / (float)amp_max ) );
  float a_pitch = to_end_start / ( -1.0 * logf( (float)pitchinitial / (float)PitchValueMax ) );
  int MOTsamplecount = 0;
  int MOTsamplecount_max = model_q1.nrow() - 1;
  int num_samples = Sensor1_Quat.nrow();
    
  NumericMatrix output(num_samples, 3);
  
  for (int s = 0; s < num_samples; ++s) {
    
    float error1 = vdst4(Sensor1_Quat(s, _),model_q1(MOTsamplecount, _)); // max error = 2
    float error2 = vdst4(Sensor2_Quat(s, _),model_q2(MOTsamplecount, _)); // max error = 2
    float error_total = error1 + error2; // max of 4
    float to_end = vdst4(Sensor1_Quat(s, _),model_q1(model_q1.nrow()-1, _)) + vdst4(Sensor2_Quat(s, _),model_q2(model_q2.nrow()-1, _));
    
    int amp_temp = (int) ( (float)amp_max * expf( -1.0 * ( error_total / a_amp ) ) );
    int pitch_temp = (int) ( (float)PitchValueMax * expf( -1.0 * ( to_end / a_pitch ) ) ); 
    
    int amp_current = constrain(amp_temp,0,amp_max);
    int pitch_current = constrain(pitch_temp,PitchValueMin,PitchValueMax);
    
    output(s, 0) = amp_current;
    output(s, 1) = pitch_current;
    output(s, 2) = MOTsamplecount;
    
    // determine what MOTsample we should be at, based on current QUAT vectors
    int MOTadvance = 1; // used in real-time time-warping algorithm
    if ( MOTsamplecount < 100 ) { // first 100 MOTsamplecount are a special case; just advance once (initiation too noisey)
      MOTsamplecount++;
    }
    else {
      float distn1 = vdstSq4(Sensor1_Quat(s, _),model_q1(MOTsamplecount-1, _)) + vdstSq4(Sensor2_Quat(s, _),model_q2(MOTsamplecount-1, _));
      float dist0 = vdstSq4(Sensor1_Quat(s, _),model_q1(MOTsamplecount, _)) + vdstSq4(Sensor2_Quat(s, _),model_q2(MOTsamplecount, _));
      float dist1 = vdstSq4(Sensor1_Quat(s, _),model_q1(MOTsamplecount+1, _)) + vdstSq4(Sensor2_Quat(s, _),model_q2(MOTsamplecount+1, _));
      if ( MOTadvance == 0 ) {
        if ( dist1 < dist0 ) {
          MOTadvance = 2;
        } else {
          MOTadvance = 1;
        }
      } else {
        if ( distn1 < dist0 ) {
          MOTadvance = 0;
        } else if ( dist1 < dist0 ) {
          MOTadvance = 2;
        } else {
          MOTadvance = 1;
        }
      }
      MOTsamplecount = MOTsamplecount + MOTadvance; 
    }
    if (MOTsamplecount >= MOTsamplecount_max) MOTsamplecount = MOTsamplecount_max-1; 
    
  }
  
  return output;
  
}

/*
 * Two-pivot upper-extremity reach simulation.
 * 
 * This function will take a pair of quaternion trajectories and turn them into a simulated reach trajectory
 *    in x-y-z Cartesian space. 
 *    
 * This function is not for final use in the full simulation; it is merely for testing and getting a MWE 
 *    up and running. The final version will include a more biomechanically accurate model from OpenSim.
 */

// [[Rcpp::export]]
List Initiate_two_pivot_system () {
  
  /*
   * This is a motion model for a two-pivot system. 
   *  - Assume pivot point p1 is fixed. 
   *  - Attached to p1 is a rigid rod with termination at point p2. 
   *  - p2 is a second pivot point with another rigid rod attached, terminating at point t. 
   *  - On the first rod between p1 and p2 is a 3-axis gyro, located at position s1, with gyro axes (unit) pointing to a1x, a1y, a1z.
   *  - Similarly, between p2 and t is another gyro, located at s2, with axes pointing to a2x, a2y, and a2z. 
   *  - We will initiate the system in the "initial position" of the experiment (seated, arm across the body resting on table). 
   * These point coordinates would all actually be empirical values measured values from the participant's arm; 
   *    However, for a rought-and-ready simulation to test our code, we can make them up.
   */
  
  // Set coordinates so that x is pointed forward, y is pointed up, and z is pointed right. 
  
  // Body segment lengths and orientations
  float upper_arm_length = 30.0; // Assume upper arm is 30cm long.
  float forearm_length = 28.0; // Assume forearm is 28cm long.
  float upper_arm_angle_off_vertical = 15.0; // Assume upper arm is 30 degrees off vertical.
  
  // Coordinate system 
  NumericVector x_unit = NumericVector::create(1.0, 0.0, 0.0);
  NumericVector y_unit = NumericVector::create(0.0, 1.0, 0.0);
  NumericVector z_unit = NumericVector::create(0.0, 0.0, 1.0);
  NumericVector origin = NumericVector::create(0.0, 0.0, 0.0);
  
  // Place body segments into coordinate system 
  NumericVector upperarm = NumericVector::create(0.0, upper_arm_length, 0.0); 
  NumericVector upperarm_straight = NumericVector::create(0.0, upper_arm_length, 0.0); 
  NumericVector shoulder_quat = formquat(upper_arm_angle_off_vertical*(3.1416/180.0), z_unit); // Assume shoulder is pointing straight up
  upperarm = qvq(shoulder_quat,upperarm);
  NumericVector forearm = NumericVector::create(0.0, 0.0, -1*forearm_length); // Assume forearm is pointing straight left.
  
  // Set up two-pivot system 
  List two_pivot_system = List::create(
    Named("p1") = upperarm, // place shoulder at upper end of upper arm
    Named("p2") = origin, // place elbow at origin
    Named("t") = forearm, // place wrist at end of forearm
    Named("s1") = vscaler_mult(forearm,0.67), // place sensor 1 2/3 of the way down the forearm
    Named("a1x") = vadd(vscaler_mult(forearm,0.67),vscaler_mult(x_unit,-1.0)), 
    Named("a1y") = vadd(vscaler_mult(forearm,0.67),vscaler_mult(z_unit,-1.0)), 
    Named("a1z") = vadd(vscaler_mult(forearm,0.67),y_unit), 
    Named("s2") = vscaler_mult(upperarm,0.5), // place sensor 2 1/2 of the way down the upper arm
    Named("a2x") = qvq(shoulder_quat,vadd(vscaler_mult(upperarm_straight,0.5),z_unit)), 
    Named("a2y") = qvq(shoulder_quat,vadd(vscaler_mult(upperarm_straight,0.5),y_unit)), 
    Named("a2z") = qvq(shoulder_quat,vadd(vscaler_mult(upperarm_straight,0.5),x_unit))
   );
  
  return two_pivot_system;
  
}

// [[Rcpp::export]]
List compute_two_pivot_next_step ( List two_pivot_system, NumericVector rotd1, NumericVector rotd2) {
  
  /*
   * This is a motion model for a two-pivot system. 
   * Assume pivot point p1 is fixed. 
   * Attached to p1 is a rigid rod with termination at point p2. 
   * p2 is a second pivot point with another rigid rod attached, terminating at point t. 
   * On the first rod between p1 and p2 is a 3-axis gyro, located at position s1, with gyro axes (unit) pointing to a1x, a1y, a1z.
   * Similarly, between p2 and t is another gyro, located at s2, with axes pointing to a2x, a2y, and a2z. 
   * If rotd1 gives the gyro readings (in radians) of the sensor at s1 and rotd2 gives the readings of the sensor at s2, 
   * then the below code computes the new positions of all points for a given sensor reading. 
   */
  
  // Extract the points from the two-pivot system object
  NumericVector p1 = two_pivot_system["p1"];
  NumericVector p2 = two_pivot_system["p2"];
  NumericVector t = two_pivot_system["t"];
  NumericVector s1 = two_pivot_system["s1"];
  NumericVector a1x = two_pivot_system["a1x"];
  NumericVector a1y = two_pivot_system["a1y"];
  NumericVector a1z = two_pivot_system["a1z"];
  NumericVector s2 = two_pivot_system["s2"];
  NumericVector a2x = two_pivot_system["a2x"];
  NumericVector a2y = two_pivot_system["a2y"];
  NumericVector a2z = two_pivot_system["a2z"];
  
  // Rotate points on limb segment 1 as if pivoting about s1
  NumericVector q1 = rot_quat(rotd1,a1x,a1y,a1z,s1);
  NumericVector p1_new = quatrot(q1,p1,s1); // p1 (shoulder joint) not actually rotating, but need to compute this to find translation for system
  NumericVector p2_new = quatrot(q1,p2,s1);
  NumericVector s1_new = s1; // quatrot(q1,s1,s1); // unecessary, as s1 is translated into origin
  NumericVector a1x_new = quatrot(q1,a1x,s1);
  NumericVector a1y_new = quatrot(q1,a1y,s1);
  NumericVector a1z_new = quatrot(q1,a1z,s1);

  // This would have moved p1, which is actually stationary (shoulder joint in reach), so, translate the L1 system back to p1
  NumericVector trans1 = vsubtract(p1,p1_new);
  p2_new = vadd(p2_new, trans1); // because p1 is stationary, we want to add the translation, i.e. p1_new + trans1 = p1
  s1_new = vadd(s1_new, trans1);
  a1x_new = vadd(a1x_new, trans1);
  a1y_new = vadd(a1y_new, trans1);
  a1z_new = vadd(a1z_new, trans1);

  // Rotate points on limb segment 2 as if pivoting about s2
  NumericVector q2 = rot_quat(rotd2,a2x,a2y,a2z,s2);
  NumericVector p2o_new = quatrot(q2,p2,s2); // The origin of limb segment 2; the first time (when computing p2_new) we treated p2 as the termination of LS1, now treating it as origin of LS2
  NumericVector t_new = quatrot(q2,t,s2);
  NumericVector s2_new = s2; // quatrot(q2,s2,s2); // unecessary, as s2 is translated into origin
  NumericVector a2x_new = quatrot(q2,a2x,s2);
  NumericVector a2y_new = quatrot(q2,a2y,s2);
  NumericVector a2z_new = quatrot(q2,a2z,s2);

  // Must reattach the origin of limb segment 2 to the termination of limb segment 1
  NumericVector trans2 = vsubtract(p2_new, p2o_new);
  t_new = vadd(t_new, trans2);
  s2_new = vadd(s2_new, trans2);
  a2x_new = vadd(a2x_new, trans2);
  a2y_new = vadd(a2y_new, trans2);
  a2z_new = vadd(a2z_new, trans2);
  
  // Finally, update the positions
  List updated_two_pivot_system = List::create(
    Named("p1") = p1, // Remember, p1 should *not* be updated to p1 new!
    Named("p2") = p2_new, 
    Named("t") = t_new, 
    Named("s1") = s1_new, 
    Named("a1x") = a1x_new, 
    Named("a1y") = a1y_new, 
    Named("a1z") = a1z_new, 
    Named("s2") = s2_new, 
    Named("a2x") = a2x_new, 
    Named("a2y") = a2y_new, 
    Named("a2z") = a2z_new
  );
  
  return updated_two_pivot_system;
  
}