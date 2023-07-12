# [OccipitalOrdinance](https://zalo.github.io/OccipitalOrdinance/)

<p align="left">
  <a href="https://github.com/zalo/OccipitalOrdinance/deployments/activity_log?environment=github-pages">
      <img src="https://img.shields.io/github/deployments/zalo/OccipitalOrdinance/github-pages?label=Github%20Pages%20Deployment" title="Github Pages Deployment"></a>
  <a href="https://github.com/zalo/OccipitalOrdinance/commits/main">
      <img src="https://img.shields.io/github/last-commit/zalo/OccipitalOrdinance" title="Last Commit Date"></a>
  <!--<a href="https://github.com/zalo/OccipitalOrdinance/blob/main/LICENSE">
      <img src="https://img.shields.io/github/license/zalo/OccipitalOrdinance" title="License: Apache V2"></a> -->
</p>

A tribute to Cortex Command; re-rendering one of the levels using [Neohookean MLS-MPM by Mykhailo Moroz](https://www.shadertoy.com/view/wdGcRK).

 # Building

This can either be run without building (in Chrome/Edge/Opera since raw three.js examples need [Import Maps](https://caniuse.com/import-maps)), or built with:
```
npm install
npm run build
```
If building manually, make sure to edit the index .html to point from `"./src/main.js"` to `"./build/main.js"`.

 # Dependencies
 - [three.js](https://github.com/mrdoob/three.js/) (3D Rendering Engine)
 - [esbuild](https://github.com/evanw/esbuild/) (Bundler)
