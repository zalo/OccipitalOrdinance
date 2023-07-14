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

    fetch('./assets/ReintegrationTracking.glsl')
      .then(data => data.text())
      .then(shaderText => {
      new THREE.TextureLoader().load('./assets/Back.bmp', (backTexture) => {
        new THREE.TextureLoader().load('./assets/Test-BG.png', (bgTexture) => {
          new THREE.TextureLoader().load('./assets/Test-FG.png', (texture) => {
            this.commonFunctions           = shaderText;
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
    });

    window.addEventListener('resize'     , this.resize.bind(this));
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
    this.bufferC = this.reintegrationComputation.addVariable("iChannel2", undefined, 3);

    this.bufferAPass = this.reintegrationComputation.addPass(this.bufferA, [this.bufferC], `
      uniform sampler2D scene;
      uniform float outputResolutionScale;
      `+this.commonFunctions+`

      layout(location = 0) out vec4 XVOut;
      layout(location = 1) out vec4 MCOut;
      layout(location = 2) out vec4 DOut;

      // Reintegration tracking
      void main(){
          vec2 pos = gl_FragCoord.xy;
          
          vec2   X = vec2(0);  // Position (relative to pixel)
          vec2   V = vec2(0);  // Velocity
          float  M = 0.;       // Mass
          vec2   C = vec2(0.); // Original UV Position
          mat2   D = mat2(0);  // Deformation Gradient

          Grid2Particle(pos, iChannel2, scene, X, V, M, C, D);

          XVOut = clamp(vec4(X, V), -1.0, 1.0);
          MCOut = vec4(M, clamp(C, 0.0, 1.0), 1.0);
           DOut = vec4(D);
      }`);
    Object.assign(this.bufferAPass.material.uniforms, this.uniforms);
    this.bufferAPass.material.uniformsNeedUpdate = true;
    this.bufferAPass.material.needsUpdate = true;

    this.bufferCPass = this.reintegrationComputation.addPass(this.bufferC, [this.bufferA], `
      uniform float outputResolutionScale;
      uniform vec3 iMouse;
      `+this.commonFunctions+`

      //velocity update

      layout(location = 0) out vec4 XVOut;
      layout(location = 1) out vec4 MCOut;
      layout(location = 2) out vec4 DOut;

      void main() {
          vec2 pos = gl_FragCoord.xy;// / resolution.xy
          ivec2 p = ivec2(pos);

          vec2   X = vec2(0);  // Position (relative to pixel)
          vec2   V = vec2(0);  // Velocity
          float  M = 0.;       // Mass
          vec2   C = vec2(0.); // Original UV Position
          mat2   D = mat2(0);  // Deformation Gradient

          Particle2Grid(pos, iChannel0, iMouse, X, V, M, C, D);
          
          XVOut = clamp(vec4(X, V), -1.0, 1.0);
          MCOut = vec4(M, clamp(C, 0.0, 1.0), 1.0);
           DOut = vec4(D);
      }`);
    Object.assign(this.bufferCPass.material.uniforms, this.uniforms);
    this.bufferCPass.material.uniformsNeedUpdate = true;
    this.bufferCPass.material.needsUpdate = true;

    const error = this.reintegrationComputation.init();
    if ( error !== null ) { console.error( error ); }

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
            // compute the smoothed density and velocity
            range(i, -2, 2) range(j, -2, 2) {
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

          vec3 background    = texture(sceneBG  , vUv).rgb;
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
