uniform vec2 iResolution;
uniform int iFrame;
#define PI 3.14159265
#define dt 0.5
#define R iResolution.xy
// useful functions
#define GS(x) exp(-dot(x,x))
#define GS0(x) exp(-length(x))
#define CI(x) smoothstep(1.0, 0.9, length(x))
#define Dir(ang) vec2(cos(ang), sin(ang))
#define Rot(ang) mat2(cos(ang), sin(ang), -sin(ang), cos(ang))
#define loop(i,x) for(int i = 0; i < x; i++)
#define range(i,a,b) for(int i = a; i <= b; i++)

// squishy solid
uniform float relax; //0.000
uniform float distribution_size; // 1.0
// Lamé parameters for stress-strain relationship
uniform float elastic_lambda; // 3.2
uniform float elastic_mu; // 4.2
uniform float incompressible_viscosity; // 1.0

#define max_vel 1.0
#define rad 2
#define diffusion 0.4
#define gravity 0.0001
#define rest_rho 1.0
#define eos_pressure 0.5

// viscous fluid
//#define relax 0.05
//#define distribution_size 0.98
//// Lamé parameters for stress-strain relationship
//#define elastic_lambda 0.2
//#define elastic_mu 0.1
//#define incompressible_viscosity 0.05

// estimation str
#define difd 2.0
// target density
#define trho 0.
// density target strenght
#define rhoe 0.0

//pressure equation of state
float pressure(float rho)
{
    return 2.2*(rho - 1.2); //gas
   // return eos_pressure*(pow(rho/rest_rho,4.) - 1.0); //Tait EOS (water)
}

void InitialConditions(inout float m, inout vec2 v, vec2 P)
{
    vec2 dx = P - R*vec2(0.3, 0.5);

    float d = smoothstep(R.y*0.5, R.y*0.49, length(dx));  
    m = 0.0*d;

    v = d*0.3*normalize(vec2(dx.y,-dx.x));
}

//KERNEL FUNCTIONS
// Linear kernel
float k0(vec2 dx) {
    vec2 k = max(1.0 - abs(dx), 0.);
    return k.x*k.y;
}
// Linear kernel with the center of mass
vec3 K0(vec2 dx) {
    vec2 k = max(1.0 - abs(dx), 0.);
    return vec3(dx*0.5, k.x*k.y);
}
// Quadratic kernel
float k1(vec2 dx) 
{
    vec2 f = max(1.5 - abs(dx), 0.0);
    vec2 k = min(max(0.75 - dx*dx, 0.5), 0.5*f*f);
    return k.x*k.y;
}
//// Box size enstimator (FROM FLUID SIMULATION: https://www.shadertoy.com/view/7lV3DD )
//vec2 destimator(vec2 dx){
//    return diffusion*dt+clamp(1.0 - 2.0*abs(dx), 0.001, 1.0);
//}
// estimating the in-cell distribution size
vec2 destimator(vec2 dx, float M){
    //size estimate by in-cell location
    vec2 ds = distribution_size*clamp(1.0 - difd*abs(dx), 0.001, 1.0);
    return ds + 0.3*max(M/(ds.x*ds.y) - 1.1, 0.)*dt;
}
// box overlap with center of mass
vec3 overlap(vec2 dx, vec2 box) {
    vec2 min0 = max(dx - box*0.5, -0.5); 
    vec2 max0 = min(dx + box*0.5, 0.5); 
    vec2 size = max(max0 - min0, 0.); 
    return vec3(0.5*(max0 + min0), size.x*size.y/(box.x*box.y));
}
float deformation_energy(mat2 D){
    D = transpose(D)*D;
    return 2.*(D[0][0]*D[0][0] + D[1][1]*D[1][1] - 2.0);
}
// MD force
float MF(vec2 dx, vec2 dv){
    return incompressible_viscosity*dot(dx,dv)*GS(0.8*dx);
}
float Ha(vec2 x){
    return ((x.x >= 0.)?1.:0.)*((x.y >= 0.)?1.:0.);
}
float Hb(vec2 x){
    return ((x.x > 0.)?1.:0.)*((x.y > 0.)?1.:0.);
}
float sdBox( in vec2 p, in vec2 b ){
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}
float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }
vec2 opRepLim(in vec2 p, in vec2 c, in vec2 l){
    return p-c*clamp(round(p/c),-l,l);
}
float border(vec2 p) {
  float bound = -sdBox( p - R*0.5            , R*vec2(0.5, 0.495)); // <--- USED TO BE .48,.48 
  float box   =  sdBox((p - R*vec2(0.5, 0.6)), R*vec2(0.05, 0.01));
  float drain = -sdBox( p - R*vec2(0.5, 0.7) , R*vec2(0.0 , 0.0 ));
  return bound;
}
#define h 1.
vec3 bN(vec2 p) {
  vec3 dx = vec3(-h,0,h);
  vec4 idx = vec4(-1./h, 0., 1./h, 0.25);
  vec3 r = idx.zyw*border(p + dx.zy)
         + idx.xyw*border(p + dx.xy)
         + idx.yzw*border(p + dx.yz)
         + idx.yxw*border(p + dx.yx);
  return vec3(normalize(r.xy), r.z + 1e-4);
}
vec2 normalize2(vec2 a) { return all(equal(a, vec2(0.)))?vec2(0.):normalize(a); }
mat2 strain(mat2 D) {
  float J = abs(determinant(D)) + 0.001;

  // MPM course, page 46
  float volume = J;

  // useful matrices for Neo-Hookean model
  mat2 F_T = transpose(D);
  mat2 F_inv_T = inverse(F_T);
  mat2 F_minus_F_inv_T = D - F_inv_T;

  // MPM course equation 48
  mat2 P_term_0 = elastic_mu * (F_minus_F_inv_T);
  mat2 P_term_1 = elastic_lambda * log(J) * F_inv_T;
  mat2 P = P_term_0 + P_term_1;

  // equation 38, MPM course
  mat2 stress = P*F_T;

  return volume * stress;
}

void Particle2Grid(in vec2 pos, in sampler2D[3] iChannel2, in sampler2D scene, out vec2 X, out vec2 V, out float M, out vec2 C, out mat2 D) {
  // Integrate over all updated neighbor distributions that fall inside of this pixel
  // This makes the tracking conservative
  range(i, -1, 1) range(j, -1, 1) {
      vec2 di = vec2(i,j);
      vec2 tpos = pos + di;
      vec4   XV =      texelFetch(iChannel2[0], ivec2(mod(tpos,R)), 0);
      vec4   MC =      texelFetch(iChannel2[1], ivec2(mod(tpos,R)), 0);
      mat2   D0 = mat2(texelFetch(iChannel2[2], ivec2(mod(tpos,R)), 0));

      vec2  X0 = XV.xy + tpos;
      vec2  V0 = XV.zw;
      float M0 = MC.x;
      
      
      // particle distribution size
      vec2 box = destimator(X0 - tpos , MC.x);
      //*
      X0 += V0*dt; //integrate position
      
      vec4 aabbX = vec4(max(pos - 0.5, X0 - box*0.5), min(pos + 0.5, X0 + box*0.5)); //overlap aabb
      vec2 center = 0.5*(aabbX.xy + aabbX.zw); //center of mass
      vec2 size = max(aabbX.zw - aabbX.xy, 0.); //only positive
      
      
      // the deposited mass into this cell
      vec3 m = MC.x*vec3(center, 1.0)*size.x*size.y/(box.x*box.y);

      X += m.xy;   // Add weighted by mass
      V += V0*m.z;
      C += m.z*MC.yz;
      M += m.z;    // Add mass
      D += D0*m.z; // Add deformation grad weighted by mass

      /*/
      
      //vec2 box = destimator(X0);   // Estimate the shape of the distribution
      X0      += di + V0*dt;        // Update particle position
      vec3 o   = overlap(X0, box);  // Find cell contribution
      X       += M0*o.xy*o.z;       // Update distribution
      M       += M0*o.z;            // Add mass
      D       += D0*o.z;            // Add deformation grad weighted by mass
      C       += o.z*MC.yz;
      float w  = k1(X0);            // Find grid node contribution
      V       += (V0 + D0*X0)*w*M0; // Distribute momentum onto grid

      //rho     += M0*w;
      //*/
  }
  
  // normalization
  if(M != 0.) {
      X /= M;
      V /= M;
      C /= M;
      D /= M;
  } else { D = mat2(1.0); }

  // initial condition
  if(iFrame < 1) {
    // Set the initial position and velocity
    X = pos;
    V = vec2(0.);

    // Initialize Each Pixel based on the original Texture (TODO: MOVE THIS TO THE TEXTURE CONSTRUCTOR)
    vec4 sceneT = texelFetch(scene, ivec2(pos), 0);
    M = float(!(sceneT.rgb == vec3(1.0, 0.0, 1.0) || sceneT.rgb == vec3(128.0/255.0, 0.0, 128.0/255.0) || sceneT.rgb == vec3(0.0)));
    C = pos/R;//mod(pos/R, 1.);
    D = mat2(1.0);
  }
  
  X = X - pos;
}


void Grid2Particle(in vec2 pos, in sampler2D[3] iChannel0, in vec3 iMouse, out vec2 X, out vec2 V, out float M, out vec2 C, out mat2 D) {
  vec4 XV = texelFetch(iChannel0[0], ivec2(mod(pos,R)), 0);
  vec4 MC = texelFetch(iChannel0[1], ivec2(mod(pos,R)), 0);

  D = mat2(texelFetch(iChannel0[2], ivec2(mod(pos,R)), 0));
  X = XV.xy + pos;
  V = XV.zw;
  M = clamp(MC.x, 0., 2.0);
  C = MC.yz;

  // not vacuum
  if (M > 0.0) {

      // Compute the velocity gradient matrix
      float a = 0.01;
      mat2 B = mat2(0.);
      float rho = 0.;

      // Compute the force
      vec2 F = vec2(0.);
      float b = 0.;

      mat2 local_strain = strain(D);

      range(i, -2,2) range(j, -2, 2) {
        vec2 tpos = pos + vec2(i,j);
        vec4 XV0 = texelFetch(iChannel0[0], ivec2(mod(tpos,R)), 0);
        vec4 MC0 = texelFetch(iChannel0[1], ivec2(mod(tpos,R)), 0);

        vec2  X0 = XV0.xy + tpos;
        vec2  V0 = XV0.zw;
        float M0 = MC0.x;
        vec2  dx = X0 - X;
        vec2  dv = V0 - V;
        mat2  D0 = mat2(texelFetch(iChannel0[2], ivec2(mod(tpos,R)), 0));

        if(!(i == 0 && j == 0)) {
          float weight = GS(0.8*dx);
          //F += M0*strain((D0*M + D*M0)/(M+M0))*dx*weight;
          mat2 strain0 = (M0*strain(D0) + M*local_strain)/(M0+M) + mat2(2.4*dot(dx,dv));
          F += M0*strain0*dx*weight;
          b += weight;
        }
        if(M>0.01){
          float weight = clamp(M0, 0.0, 1.0)*GS(0.8*dx);
          rho += M0*weight;
          B += mat2(dv*dx.x,dv*dx.y)*weight;
          a += weight;
        }
      }
      B   /= a;
      rho /= a;
      F    = clamp(F/b, -0.4,0.4);

      if(iMouse.z > 0.){
          vec2 dx= pos - iMouse.xy;
          F += 0.01*normalize2(dx)*GS(dx/80.);
      }
      
      // Gravity
      F += gravity*vec2(0,-1);
      
      // integrate velocity
      V += F*dt;
      //X +=  0.*F*dt;
      
      // Apply Border effects
      vec3 BORD = bN(X);
      V += 0.1*smoothstep(0., 5., -BORD.z)*BORD.xy;
      V *= 1. - 0.5*smoothstep(-30., 0., -pos.y);
      
      //velocity limit
      float v = length(V);
      V /= (v > 1.)?1.*v:1.;

      float drho = rho - 1.0;
      B -= 0.007*mat2(drho)*abs(drho);
     
      //integrate deformation gradient
      D += 1.*dt*B*D;
     
      //smoothing
      
      float r = relax + 0.05*smoothstep(-30., 0., -pos.y);
      D = D*(1. - r) + mat2(1.)*r;
      
      //clamp the gradient to not go insane
      D = mat2(clamp(vec4(D - mat2(1.)), -5.0, 5.0)) + mat2(1.);
  }
  
  //save
  X = X - pos;
  M = MC.x;
}
