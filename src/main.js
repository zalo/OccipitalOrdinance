import * as THREE from 'three';
import { GUI } from '../node_modules/three/examples/jsm/libs/lil-gui.module.min.js';
import { MultiTargetGPUComputationRenderer } from './MultiTargetGPUComputationRenderer.js';

class OccipitalOrdinance {
  constructor() {
    //let debug = new Debug();

    this.container = document.createElement( 'div' );
    document.body.appendChild( this.container );

    this.scene = new THREE.Scene();
    this.scene.name = 'Scene';

    this.camera = new THREE.PerspectiveCamera( 35, window.innerWidth / window.innerHeight, 1.0, 2000.0 );
    this.camera.position.set( 0, 0, 100 );
    this.scene.add(this.camera);

    this.time   = 0;
    this.width  = 800;
    this.height = 450;

    this.renderer = new THREE.WebGLRenderer({ antialias: true, transparent: true, alpha: true });
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setSize(window.innerWidth, window.innerHeight);
    this.renderer.setAnimationLoop(this.render.bind(this));
    this.renderer.setClearAlpha(0.0);
    this.renderer.setClearColor(new THREE.Color(1, 1, 1), 0.0);
    this.renderer.domElement.style.position = 'fixed';
    this.renderer.domElement.style.zIndex   = '-1000';
    this.container.appendChild(this.renderer.domElement);
    this.container.style.touchAction = "none";
    this.container.style.userSelect  = "none";

    this.raycaster = new THREE.Raycaster();
    this.pointer   = new THREE.Vector3();

    new THREE.TextureLoader().load('./assets/Back.bmp', (backTexture) => {
      new THREE.TextureLoader().load('./assets/Test-BG.png', (bgTexture) => {
        new THREE.TextureLoader().load('./assets/Test-FG.png', (texture) => {
          this.testTexture               = texture;
          this.testTexture.minFilter     = THREE.NearestFilter;
          this.testTexture.magFilter     = THREE.NearestFilter;
          this.testTextureBG             = bgTexture;
          this.testTextureBG.minFilter   = THREE.NearestFilter;
          this.testTextureBG.magFilter   = THREE.NearestFilter;
          this.testTextureBack           = backTexture;
          this.testTextureBack.minFilter = THREE.NearestFilter;
          this.testTextureBack.magFilter = THREE.NearestFilter;

          this.width  = this.testTexture.source.data.width;
          this.height = this.testTexture.source.data.height;

          this.uniforms = {
            scene       : { value: this.testTexture },
            sceneBG     : { value: this.testTextureBG },
            sceneBack   : { value: this.testTextureBack },
            outputResolutionScale : { value: 1.0 },
            iterations  : { value: 32 },
            iFrame      : { value: 0 },
            iMouse      : { value: new THREE.Vector3() },
            iChannel0   : { value: null },
            //iChannel1   : { value: null },
            iResolution : { value: new THREE.Vector2(this.width, this.height) }
          }

          this.gui = new GUI()
          this.gui.add(this.uniforms.iterations, 'value', 1,  64).name('Iterations Per Frame');
          this.gui.add(this.uniforms.outputResolutionScale, 'value', 0.1,  1.0).name('Resolution Scale')
            .onChange(() => { this.renderer.setPixelRatio(window.devicePixelRatio * this.uniforms.outputResolutionScale.value); });
          this.gui.open();

          this.createReintegrationSystem();
        });
      });
    });

    window.addEventListener('resize', this.resize.bind(this));
    window.addEventListener('pointermove', this.onPointerMove.bind(this));
    window.addEventListener('pointerdown', ()=>{ this.pointer.z = 1.0; });
    window.addEventListener('pointerup'  , ()=>{ this.pointer.z = 0.0; });
    this.resize();

    this.lastTime = this.time;
  }

  onPointerMove( event ) {
    this.pointer.x =   ( event.clientX / window.innerWidth  ) * 2 - 1;
    this.pointer.y = - ( event.clientY / window.innerHeight ) * 2 + 1;
  }

  createReintegrationSystem() {
    if (this.reintegrationComputation) {
      this.scene.remove(this.labelMesh);
      this.scene.remove(this.reintegrationMesh);
      this.reintegrationComputation.dispose();
      this.labelMaterial.dispose();
      this.reprojectionMaterial.dispose();
    }

    this.reintegrationComputation = new MultiTargetGPUComputationRenderer(this.width, this.height, this.renderer);
    this.bufferA = this.reintegrationComputation.addVariable("iChannel0", undefined, 3);
    this.bufferC = this.reintegrationComputation.addVariable("iChannel2", undefined, 2);
    this.bufferD = this.reintegrationComputation.addVariable("iChannel3");

    this.commonFunctions = `
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
      #define velocity_averaging 0.
      
      // squishy solid
      #define relax 0.000
      #define distribution_size 1.0
      // Lamé parameters for stress-strain relationship
      #define elastic_lambda 3.2
      #define elastic_mu 4.2
      #define incompressible_viscosity 1.0
      // viscous fluid
      /*
      #define relax 0.05
      #define distribution_size 0.98
      // Lamé parameters for stress-strain relationship
      #define elastic_lambda 0.2
      #define elastic_mu 0.1
      #define incompressible_viscousity 0.05
      */ 

      // estimation str
      #define difd 2.0
      // target density
      #define trho 0.
      // density target strenght
      #define rhoe 0.0
      
      // estimating the in-cell distribution size
      vec2 destimator(vec2 dx, float M){
          //size estimate by in-cell location
          vec2 ds = distribution_size*clamp(1.0 - difd*abs(dx), 0.001, 1.0);
          return ds + 0.3*max(M/(ds.x*ds.y) - 1.1, 0.)*dt;
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
      }`;

    this.bufferAPass = this.reintegrationComputation.addPass(this.bufferA, [this.bufferC, this.bufferD], `
      uniform sampler2D scene;
      uniform float outputResolutionScale;
      `+this.commonFunctions+`

      layout(location = 0) out vec4 XVOut;
      layout(location = 1) out vec4 MCOut;
      layout(location = 2) out vec4 DOut;

      // Reintegration tracking
      void main(){
          vec2 pos = gl_FragCoord.xy;// / resolution.xy
          ivec2 p  = ivec2(pos);
          
          vec2  X = vec2(0);
          vec2  V = vec2(0);
          float M = 0.;
          vec2  C = vec2(0.);

          // deformation gradient
          mat2 D = mat2(0);

          // basically integral over all updated neighbor distributions
          // that fall inside of this pixel
          // this makes the tracking conservative
          range(i, -1, 1) range(j, -1, 1) {
              vec2 tpos = pos + vec2(i,j);
              vec4 XV = texelFetch(iChannel2[0], ivec2(mod(tpos,R)), 0);
              vec4 MC = texelFetch(iChannel2[1], ivec2(mod(tpos,R)), 0);
              mat2 D0 = mat2(texelFetch(iChannel3, ivec2(mod(tpos,R)), 0));

              vec2 X0 = XV.xy + tpos;
              vec2 V0 = XV.zw;
              
              // particle distribution size
              vec2 K = destimator(X0 - tpos , MC.x);
             
              X0 += V0*dt; //integrate position
              
              vec4 aabbX = vec4(max(pos - 0.5, X0 - K*0.5), min(pos + 0.5, X0 + K*0.5)); //overlap aabb
              vec2 center = 0.5*(aabbX.xy + aabbX.zw); //center of mass
              vec2 size = max(aabbX.zw - aabbX.xy, 0.); //only positive
              
              // the deposited mass into this cell
              vec3 m = MC.x*vec3(center, 1.0)*size.x*size.y/(K.x*K.y);
              
              // add weighted by mass
              X += m.xy;
              V += V0*m.z;
              C += m.z*MC.yz;
              //add mass
              M += m.z;

              // add deformation grad weighted by mass
              D += D0*m.z;
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
              X = pos;
              V = vec2(0.);

              vec4 sceneT = texelFetch(scene, ivec2(p), 0);
              M = float(!(sceneT.rgb == vec3(1.0, 0.0, 1.0) || sceneT.rgb == vec3(128.0/255.0, 0.0, 128.0/255.0) || sceneT.rgb == vec3(0.0)));
              
              C = mod(pos/R, 1.);

              D = mat2(1.0);
          }
          
          X = X - pos;
          XVOut = clamp(vec4(X, V), -1.0, 1.0);
          MCOut = vec4(M, clamp(C, 0.0, 1.0), 1.0);
          DOut  = vec4(D);
      }`);
    Object.assign(this.bufferAPass.material.uniforms, this.uniforms);
    //this.bufferAPass.material.uniforms["map"] = { value: this.testTexture };
    this.bufferAPass.material.uniformsNeedUpdate = true;
    this.bufferAPass.material.needsUpdate = true;

    this.bufferCPass = this.reintegrationComputation.addPass(this.bufferC, [this.bufferA], `
      uniform float outputResolutionScale;
      uniform vec3 iMouse;
      `+this.commonFunctions+`
      //velocity update

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
      
      vec2 normalize2(vec2 a) {
          return all(equal(a, vec2(0.)))?vec2(0.):normalize(a);
      }
      
      
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
      
      layout(location = 0) out vec4 XVOut;
      layout(location = 1) out vec4 MCOut;

      void main() {
          vec2 pos = gl_FragCoord.xy;// / resolution.xy
          vec2 uv = pos/R;
          ivec2 p = ivec2(pos);

          vec4 XV = texelFetch(iChannel0[0], ivec2(mod(pos,R)), 0);
          vec4 MC = texelFetch(iChannel0[1], ivec2(mod(pos,R)), 0);

          mat2  D = mat2(texelFetch(iChannel0[2], ivec2(mod(pos,R)), 0));
          vec2  X = XV.xy + pos;
          vec2  V = XV.zw;
          float M = clamp(MC.x, 0., 2.0);
          vec2  C = MC.yz;
          // not vacuum
          if (M > 0.0) {
              // Compute the force
            
              vec2 F = vec2(0.);
              float b = 0.;
         
              mat2 local_strain = strain(D);
              if(M > 0.0)
              {
                  range(i, -2,2) range(j, -2, 2)
                  {
                      if(!(i == 0 && j == 0))
                      {
                          vec2 tpos = pos + vec2(i,j);
                          vec4 XV0 = texelFetch(iChannel0[0], ivec2(mod(tpos,R)), 0);
                          vec4 MC0 = texelFetch(iChannel0[1], ivec2(mod(tpos,R)), 0);

                          vec2  X0 = 0.*XV0.xy + tpos; // INVESTIGATE WHY THIS IS MULTIPLIED BY 0.0
                          vec2  V0 = XV0.zw;
                          float M0 = MC0.x;
                          vec2  dx = X0 - X;
                          vec2  dv = V0 - V;
                          mat2  D0 = mat2(texelFetch(iChannel0[2], ivec2(mod(tpos,R)), 0));
                          float weight = GS(0.8*dx);
                         
                          //F += M0*strain((D0*M + D*M0)/(M+M0))*dx*weight;
                          mat2 strain0 = (M0*strain(D0) + M*local_strain)/(M0+M) + mat2(2.4*dot(dx,dv));
                          F += M0*strain0*dx*weight;
                         
                          b += weight;
                      }
                  }
             
                  F /= b;
                  F = clamp(F, -0.4,0.4);
              }
              if(iMouse.z > 0.)
              {
                  vec2 dx= pos - iMouse.xy;
                  F += 0.01*normalize2(dx)*GS(dx/80.);
              }
              
              //gravity
              F += 0.0001*vec2(0,-1);
              
              //integrate velocity
              V += F*dt;
              //X +=  0.*F*dt;
              
              vec3 BORD = bN(X);
              V += 0.1*smoothstep(0., 5., -BORD.z)*BORD.xy;
              V *= 1. - 0.5*smoothstep(-30., 0., -pos.y);
              
              //velocity limit
              float v = length(V);
              V /= (v > 1.)?1.*v:1.;
          }
          
          //save
          X = X - pos;
          XVOut = clamp(vec4(X, V), -1.0, 1.0);
          MCOut = vec4(MC.x, clamp(C, 0.0, 1.0), 1.0);
      }`);
    Object.assign(this.bufferCPass.material.uniforms, this.uniforms);
    this.bufferCPass.material.uniformsNeedUpdate = true;
    this.bufferCPass.material.needsUpdate = true;

    this.bufferDPass = this.reintegrationComputation.addPass(this.bufferD, [this.bufferA], `
      out highp vec4 pc_fragColor;
      uniform float outputResolutionScale;
      `+this.commonFunctions+`
      vec2 normalize2(vec2 a) { return all(equal(a, vec2(0.)))?vec2(0.):normalize(a); }
      
      void main() {
          vec2 pos = gl_FragCoord.xy;// / resolution.xy
          ivec2 p = ivec2(pos);

          vec4 XV = texelFetch(iChannel0[0], ivec2(mod(pos,R)), 0);
          vec4 MC = texelFetch(iChannel0[1], ivec2(mod(pos,R)), 0);

          mat2  D = mat2(texelFetch(iChannel0[2], ivec2(mod(pos,R)), 0));
          vec2  X = XV.xy + pos;
          vec2  V = XV.zw;
          float M = MC.x;
          
          if(M > 0.01) //not vacuum
          {
              //Compute the velocity gradient matrix
              mat2 B = mat2(0.);
              float a = 0.01;
              float rho = 0.;
              range(i, -2, 2) range(j, -2, 2)
              {
                  vec2 tpos = pos + vec2(i,j);
                  vec4 XV0 = texelFetch(iChannel0[0], ivec2(mod(tpos,R)), 0);
                  vec4 MC0 = texelFetch(iChannel0[1], ivec2(mod(tpos,R)), 0);

                  vec2  X0 = XV0.xy + tpos;
                  vec2  V0 = XV0.zw;
                  float M0 = MC0.x;
                  vec2 dx = X0 - X;
                  vec2 dv = V0 - V;
                  vec2 dsize = clamp(destimator(X0 - tpos, M0), 0.3, 1.0);
                  float weight = clamp(M0, 0.0, 1.0)*GS(0.8*dx);
                  rho += M0*weight;
                  B += mat2(dv*dx.x,dv*dx.y)*weight;
                  a += weight;
              }
              B /= a;
              rho /= a;
            
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
          pc_fragColor = vec4(D);
      }`);
    Object.assign(this.bufferDPass.material.uniforms, this.uniforms);
    //this.bufferDPass.material.uniforms["map"] = { value: this.testTexture };
    this.bufferDPass.material.uniformsNeedUpdate = true;
    this.bufferDPass.material.needsUpdate = true;


    const error = this.reintegrationComputation.init();
    if ( error !== null ) { console.error( error ); }

    ////console.log(this.isovistComputation.getCurrentRenderTarget(this.isovist));
    //this.labelMaterial = new THREE.MeshBasicMaterial( 
    //  { map: this.reintegrationComputation.getCurrentRenderTarget(this.bufferB).texture, side: THREE.DoubleSide });
    //this.labelPlane = new THREE.PlaneGeometry(50, 50);
    //this.labelMesh = new THREE.Mesh(this.labelPlane, this.labelMaterial);
    //this.labelMesh.position.set(25, 0, 0);
    //this.labelMesh.scale   .set(1, this.height/this.width, 1);
    //this.scene.add(this.labelMesh);

    //this.uniforms["isovist"] = { value: null };
    //this.uniforms["map"]     = { value: this.testTexture };

    this.reprojectionMaterial = new THREE.ShaderMaterial( {
      side: THREE.FrontSide,
      //dithering: true,
      //transparent: true,
      uniforms: this.uniforms,
      vertexShader: `
        varying vec2 vUv;
        void main() {
            vUv = uv;
            gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );
            //gl_Position = vec4( ( uv - 0.5 ) * 2.0, 0.0, 1.0 );
        }`,
      fragmentShader: `
        uniform sampler2D[3] iChannel0;
        uniform sampler2D scene, sceneBG, sceneBack;
        uniform float outputResolutionScale;
        varying vec2 vUv;

        `+this.commonFunctions+`
        vec3 hsv2rgb( in vec3 c ) {
          vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
          rgb = rgb*rgb*(3.0-2.0*rgb); // cubic smoothing	
          return c.z * mix( vec3(1.0), rgb, c.y);
        }
        
        #define radius 0.25
        #define zoom 0.2
        
        void main() {
            vec2 pos = vec2(vUv.x * iResolution.x, vUv.y * iResolution.y);
            float rho = 0.001;
            vec2 c = vec2(0.);
            float De = 0.;
            vec2 vel = vec2(0., 0.);
            vec2 grad = vec2(0.);
        
            float rho2 = 0.;
            //compute the smoothed density and velocity
            range(i, -2, 2) range(j, -2, 2)
            {
                vec2 tpos = floor(pos) + vec2(i,j);
                vec4 XV = texelFetch(iChannel0[0], ivec2(mod(tpos,R)), 0);
                vec4 MC = texelFetch(iChannel0[1], ivec2(mod(tpos,R)), 0);

                vec2  X0 = XV.xy + tpos;
                vec2  V0 = XV.zw;
                float M0 = MC.x;
                vec2  dx = X0 - pos;
                vec2 dx0 = X0 - tpos;
                mat2  D0 = mat2(texelFetch(iChannel0[2], ivec2(mod(tpos,R)), 0));
                
                float K = GS(dx/radius)/(radius*radius);
                rho  += M0*K;
                grad += normalize(dx)*K;
                c    += M0*K*MC.yz;
                De   += M0*K*abs(deformation_energy(D0));
                vel  += M0*K*V0;
                vec2 dsize = destimator(dx0,  M0);
                float bsdf = sdBox(pos - X0,0.5*dsize);
                //float bsdf = length(pos - X0) - 0.5*length(destimator(dx0));
                rho2 += M0*smoothstep(0.1, -0.1, bsdf)/(dsize.x*dsize.y);
            }
        
          grad /= rho; 
          c    /= rho;
          vel  /= rho;
          De   /= rho;

          float d = smoothstep(0.3,0.7,mix(rho, rho2,1.0));
          gl_FragColor.rgb = vec3(0.0,0.5,1.0)*De*De + 0.04*rho2 + texture(scene, c).xyz; //

          vec3 background = texture(sceneBG, vUv).rgb;
          vec3 skyBackground = texture(sceneBack, vUv).rgb;
          background = (background == vec3(1.0, 0.0, 1.0)) ? skyBackground : background;

          gl_FragColor.rgb = mix(background, gl_FragColor.rgb, d);
          gl_FragColor.a = 1.0;
        }`
    });
  
    this.reintegrationMesh = new THREE.Mesh( new THREE.PlaneGeometry( 80, 80 ), this.reprojectionMaterial );
    this.reintegrationMesh.position.set(0, 0, 0);
    this.reintegrationMesh.scale   .set(1, this.height/this.width, 1);
    this.scene.add(this.reintegrationMesh);
  }

  resize() {
    this.camera.aspect = window.innerWidth / window.innerHeight;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize( window.innerWidth, window.innerHeight);
    this.lastTime=this.time;
  }

  render(timeMS) {
    this.time = timeMS;
    if (this.time == 0) { this.lastTime = this.time; }

    if (this.reintegrationComputation) {
      this.raycaster.setFromCamera( this.pointer, this.camera );
      let intersects = this.raycaster.intersectObjects( [this.reintegrationMesh] );
      if(intersects.length > 0) { 
        this.uniforms["iMouse"].value = new THREE.Vector3(
          intersects[0].uv.x * this.width, 
          intersects[0].uv.y * this.height, this.pointer.z);
      }

      for(let i = 0; i < this.uniforms.iterations.value; i++){
        this.reintegrationComputation.compute();
        this.uniforms["iFrame"].value += 1;
      }

      this.uniforms["iChannel0"].value = this.reintegrationComputation.getCurrentRenderTarget(this.bufferA).texture;
    }

    this.renderer.render(this.scene, this.camera);
  }
}

window.OccipitalOrdinance = new OccipitalOrdinance();
