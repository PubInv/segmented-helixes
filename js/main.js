// Copyright (C) 2019 by
//   Robert L. Read <read.robert@gmail.com>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// var tm = UGLY_GLOBAL_SINCE_I_CANT_GET_MY_MODULE_INTO_THE_BROWSER;
// var OPERATION = "normal"; // "normal" or "helices"

var WINDOW_HEIGHT_FACTOR = 0.68;

var TET_DISTANCE = 0.5;

var NUM_PRISMS = 10;
var NUM_SEGMENTS = (2 * NUM_PRISMS) + 1;

// Detects webgl
if (!Detector.webgl) {
  Detector.addGetWebGLMessage();
  document.getElementById('threecontainer').innerHTML = "";
}

const PLATONIC_SOLIDS = ["TETRAHEDRON","CUBE","OCTAHEDRON","DODECAHEDRON","ICOSAHEDRON"];


// Here I attempt to create an abstract prism object.

const PRISM_FACE_RATIO_LENGTH = 1/2;


var PHI_SPRITE;

var bncG;

function renderPrismInstance(p_i,tau) {
  // now that we have the points, we can
  // construct the objects....
  var objects = [];
  let colors = [d3.color("DarkRed"), d3.color("DarkOrange"), d3.color("Blue")];

  let L = p_i.p.L;
  let PRISM_FACE_LENGTH = L * PRISM_FACE_RATIO_LENGTH;
  let SR = PRISM_FACE_LENGTH/20;

  let TB = p_i.tb;
  let LB = p_i.lb;
  let RB = p_i.rb;

  let TC = p_i.tc;
  let LC = p_i.lc;
  let RC = p_i.rc;
  objects.push(createSphere(SR,TB,colors[0].hex()));
  objects.push(createSphere(SR,LB,colors[1].hex()));
  objects.push(createSphere(SR,RB,colors[2].hex()));

  objects.push(createSphere(SR,TC,colors[0].hex()));
  objects.push(createSphere(SR,LC,colors[1].hex()));
  objects.push(createSphere(SR,RC,colors[2].hex()));

  objects.push(createSphere(SR,p_i.b,colors[2].hex()));
  objects.push(createSphere(SR,p_i.c,colors[2].hex()));

  // These are the joint axes....
  objects.push(createSphere(SR,new THREE.Vector3(0,0,-L/2),colors[0].hex()));
  objects.push(createSphere(SR,new THREE.Vector3(0,0,L/2),colors[0].hex()));


  var scolors = [d3.color("DarkRed"), d3.color("DarkOrange"), d3.color("Indigo")];
  var smats = [new THREE.Color(0x8B0000),
               new THREE.Color(0xFF8C00),
               new THREE.Color(0x000082)];

  let w = SR/4;
  objects.push(create_actuator_pure(TB,LB,SR,SR/2,memo_color_mat(smats[0])));
  objects.push(create_actuator_pure(LB,RB,SR,SR/2,memo_color_mat(smats[1])));
  objects.push(create_actuator_pure(RB,TB,SR,SR/2,memo_color_mat(smats[2])));

  objects.push(create_actuator_pure(TC,LC,SR,SR/2,memo_color_mat(smats[0])));
  objects.push(create_actuator_pure(LC,RC,SR,SR/2,memo_color_mat(smats[1])));
  objects.push(create_actuator_pure(RC,TC,SR,SR/2,memo_color_mat(smats[2])));

  objects.push(create_actuator_pure(TB,TC,SR,SR/2,memo_color_mat(smats[0])));
  objects.push(create_actuator_pure(LB,LC,SR,SR/2,memo_color_mat(smats[1])));
  objects.push(create_actuator_pure(RB,RC,SR,SR/2,memo_color_mat(smats[2])));

  // Need to return the points here, and the up vector, not just the meshes...
  // That is an "instance of an abstract prism".
  // POSITION SUPERSTRUCTURE
  // Not having computed the transform for this,
  // But rather having computed all the points in the
  // world space, we can postion SUP by moving it
  // to the mid point of nu.b and nu.c, "looking at" C,
  // and then making sure the local Y axis points at
  // the up vectors by rotating around BC.

  var transG;

  // The rejection of a onto b is perpendicular to b.
  function rejection(a,b) {
    const bu = b.clone().normalize();
    const a1 = a.dot(b);
    const a1v = bu.multiplyScalar(a1);
    const a2v = a.clone().sub(a1v);
    console.assert(near(a2v.dot(b),0));
    return a2v;
  }

  function positionSuperStructure(nu,tau) {
    // First we point axis in bc direction....
    var bc = nu.c.clone().sub(nu.b);
    let bcn = bc.clone().normalize();
    bncG = bcn.clone();

    const m = new THREE.Vector3().addVectors(nu.b,nu.c);
    m.multiplyScalar(1/2);
    const tm = new THREE.Vector3().addVectors(nu.tb,nu.tc);
    tm.multiplyScalar(1/2);
    let axis = nu.c.clone().sub(nu.b);

    // This is actually not quite right, as it is not
    // perpendicular to the BC axis.
    const prism_v = tm.clone().sub(m);
    const prism_up = rejection(prism_v,axis);
    //    console.log("V,UP",prism_v,prism_up);

    const q0 = new THREE.Quaternion();


    q0.setFromUnitVectors(new THREE.Vector3(0,0,1),bcn);

    //    nu.sup.applyQuaternion(q);
    // finally we move the object to the midpoint...

    var transm = new THREE.Matrix4();
    transm.makeTranslation(m.x,m.y,m.z);
    var rot0 = new THREE.Matrix4();
    rot0.makeRotationFromQuaternion(q0);
    var trans = transm.clone();
    trans.multiply(rot0);
    nu.sup.applyMatrix(trans);
    nu.sup.updateMatrix();
    // Now that we are in the correct postion, we want to
    // rotate about the B-C axis to make the prism up match the
    // super_structure up


    { // The entire matching of the up vector must be rethought...
      // We need to create a rotation along the BC axis, to
      // match faces, but I am not sure how to do this, and
      // I should have it from the other prism work.
      // Then we make up direction correct...
      //    const super_up = nu.sup.localToWorld(nu.sup.up);

      // I guess instead of using up, I must attach one to the object
      // or I must compute it from the centroid.
      var rot1 = new THREE.Matrix4().identity();

      var cent = findCentroid(nu.sup.children[0].geometry);

//      console.log("PRISM_UP",prism_up);
      if (vnear(cent, new THREE.Vector3(0,0,0))) {
        console.log("WARNING! Centroid near zero!");
        // In this case (at least one of them), we are
        // In the "straight line" case and the rotation is just a Z rotation
        // We have to use a different up vector here!!
      } else {
        cent.normalize();
        cent.applyMatrix4(rot0);
        cent.normalize();
        prism_up.normalize();
        console.assert(near(axis.dot(cent),0));
        if (!near(axis.dot(cent),0)) {
          debugger;
        }
        const q1 = new THREE.Quaternion();
        q1.setFromUnitVectors(cent,prism_up);

        let theta = -prism_up.angleTo(cent);
//        console.log("COMPUTED THETA",theta * 180/Math.PI,q1);

        // TODO: In our linear cases, try computing this
        // a different way and see if it comes out.
        // Test that this really does a 60 rotation.
        rot1.makeRotationFromQuaternion(q1);
        let x = new THREE.Vector3(1,0,0);
        let y = x.clone();
        x.applyMatrix4(rot1);
//        console.log("PRE, POST",y,x);
      }

      // Now that we have a rotation axis, we will have
      // to translate a point on the axis to the origin, then
      // invert that transform to get back. We'll use m.
      var transm_i = new THREE.Matrix4();
      transm_i.getInverse(transm);
      var up_v_trans = new THREE.Matrix4().identity();
      up_v_trans.multiply(transm);
      up_v_trans.multiply(rot1);
      up_v_trans.multiply(transm_i);

      nu.sup.applyMatrix(up_v_trans);
      nu.sup.updateMatrix();
      const super_up2 = findCentroid(nu.sup.children[0].geometry);
      super_up2.applyMatrix4(up_v_trans);
    }
  }

  if (p_i.p.superstructure_prototype) {
    p_i.sup =
      p_i.p.superstructure_prototype.GdeepCloneMaterials();
    positionSuperStructure(p_i,tau);
    objects.push(p_i.sup);
  }
  return objects;
}

var INITIAL_NORM_POINT_Y = -0.7;
var INITIAL_NORM_POINT_X = -0.62;

var WORLD_HEIGHT = 2.0;
var GTRANS = new THREE.Matrix4().makeTranslation(0,WORLD_HEIGHT,0);
var GLOBAL_P0 = new AbstractPrism(
  1,
  new THREE.Vector3(INITIAL_NORM_POINT_X,INITIAL_NORM_POINT_Y,-1),
  new THREE.Vector3(-INITIAL_NORM_POINT_X,INITIAL_NORM_POINT_Y,1));

function testCreatePrism() {
  var p_i = CreatePrism(GLOBAL_P0,PRISM_FACE_RATIO_LENGTH);

  // We shall place this upward, for the purpose of
  // making it easier to see...
  var TP = renderPrismInstance(p_i,0);
  console.log(TP);
  TP.forEach(o => { am.scene.add(o); });
}

// In order to do some labeling, it is valuble to return
// some of these prisms; I will return a triple:
// [origin,positives,negatives], which respects the
// symmetry of how they are created and makes it easy
// to pick them out.
function createAdjoinedPrisms(p_i,tau,num) {
  var TP = renderPrismInstance(p_i,tau);
  TP.forEach(o => { am.scene.add(o); });
  // DANGER
//  p_i.sup = null;
  var orign = TP;
  var positives = [];
  var negatives = [];
  var cur = p_i;
  for(let i = 0; i < num; i++) {
    var p_c = adjoinPrism(cur,tau,true)[0];
    var TP = renderPrismInstance(p_c,tau);
    TP.forEach(o => { am.scene.add(o); });
    positives.push(TP);
    cur = p_c;
  }
  var cur = p_i;
  for(let i = 0; i < num; i++) {
    var p_c = adjoinPrism(cur,tau,false,false)[0];
    var TP = renderPrismInstance(p_c,tau);
    TP.forEach(o => { am.scene.add(o); });
    negatives.push(TP);
    cur = p_c;
  }
  return [orign,positives,negatives];
}

function addShadowedLight(scene, x, y, z, color, intensity) {
  var directionalLight = new THREE.DirectionalLight(color, intensity);
  directionalLight.position.set(x, y, z);
  scene.add(directionalLight);
  directionalLight.castShadow = true;
  var d = 1;
  directionalLight.shadow.camera.left = -d;
  directionalLight.shadow.camera.right = d;
  directionalLight.shadow.camera.top = d;
  directionalLight.shadow.camera.bottom = -d;
  directionalLight.shadow.camera.near = 1;
  directionalLight.shadow.camera.far = 4;
  directionalLight.shadow.mapSize.width = 1024;
  directionalLight.shadow.mapSize.height = 1024;
  directionalLight.shadow.bias = -0.005;
}
function createParalellepiped(sx, sy, sz, pos, quat, material) {
  var pp = new THREE.Mesh(new THREE.BoxGeometry(sx, sy, sz, 1, 1, 1), material);
  pp.castShadow = false;;
  pp.receiveShadow = true;
  pp.position.set(pos.x, pos.y, pos.z);
  return pp;

}
// Not sure how to use the quaternion here,
function createSphere(r, pos, color) {
  //    var cmat = memo_color_mat(tcolor);
  var tcolor = new THREE.Color(color);
  var cmat = new THREE.MeshPhongMaterial({ color: tcolor });
  var ball = new THREE.Mesh(new THREE.SphereGeometry(r, 18, 16), cmat);
  ball.position.set(pos.x, pos.y, pos.z);
  ball.castShadow = false;;
  ball.receiveShadow = true;

  return ball;
}

function get_member_color(gui, len) {
  if (len < am.MIN_EDGE_LENGTH)
    return d3.color("black");
  else if (len > am.MAX_EDGE_LENGTH)
    return d3.color("black");
  else {
    var p = (len - am.MIN_EDGE_LENGTH) / (am.MAX_EDGE_LENGTH - am.MIN_EDGE_LENGTH);
    return d3.rgb(gui.color_scale(len));
  }
}

function create_actuator(b_a, b_z, pos, cmat) {
  var len = b_z.distanceTo(b_a) + -am.JOINT_RADIUS;
  var quat = new THREE.Quaternion();

  var pos = new THREE.Vector3(b_z.x, b_z.y, b_z.z);
  pos.add(b_a);
  pos.divideScalar(2);

  var mesh = createParalellepiped(
    am.INITIAL_EDGE_WIDTH,
    am.INITIAL_EDGE_WIDTH,
    len,
    pos,
    quat,
    cmat);

  mesh.lookAt(b_z);

  mesh.castShadow = false;;
  mesh.receiveShadow = true;
  am.scene.add(mesh);
  mesh.structureKind = "member";
  mesh.name = b_a.name + " " + b_z.name;
  return mesh;
}

function create_actuator_pure(b_a, b_z,jr,w, cmat) {
  var len = b_z.distanceTo(b_a) + -jr;
  var quat = new THREE.Quaternion();

  var pos = new THREE.Vector3(b_z.x, b_z.y, b_z.z);
  pos.add(b_a);
  pos.divideScalar(2);

  var mesh = createParalellepiped(
    w,
    w,
    len,
    pos,
    quat,
    cmat);

  mesh.lookAt(b_z);

  mesh.castShadow = false;;
  mesh.receiveShadow = true;
  //    am.scene.add(mesh);
  mesh.structureKind = "member";
  mesh.name = b_a.name + " " + b_z.name;
  return mesh;
}

function memo_color_mat(tcolor) {
  var string = tcolor.getHexString();
  if (!(string in am.color_material_palette)) {
    var cmat = new THREE.MeshPhongMaterial({ color: tcolor });
    am.color_material_palette[string] = cmat;
  }
  return am.color_material_palette[string]
}

var scolors = [d3.color("DarkRed"), d3.color("DarkOrange"), d3.color("Indigo")];
var smats = [new THREE.Color(0x8B0000),
             new THREE.Color(0xFF8C00),
             new THREE.Color(0x000082)];

function create_vertex_mesh(pos, c) {
  var mesh = createSphere(am.JOINT_RADIUS/2, pos, c.hex());
  mesh.castShadow = false;
  mesh.receiveShadow = false;
  am.scene.add(mesh);
  return mesh;
}

function cto3(c) {
  return new THREE.Color(c.hex());
}

function get_random_int(max) {
  return Math.floor(Math.random() * Math.floor(max));
}

function get_direction(n, v, i) {
  if (n < 10)
    return get_random_int(3);
  else return -1;
}

function get_vertex(n, v, i, pa, pb, pc, s, l, m) {
  var valid = { v: true };
  var l0 = pa.distanceTo(pb);
  var l1 = pc.distanceTo(pa);
  var l2 = pb.distanceTo(pc);
  var ad = m ? s[0]-l0 : (n % 2 == 0) ? l[0] : l[3];
  var bd = m ? s[1]-l1 : (n % 2 == 0) ? l[1] : l[4];
  var cd = m ? s[2]-l2 : (n % 2 == 0) ? l[2] : l[5];
  var pd = find_fourth_point_given_three_points_and_three_distances(
    CHIRALITY_CCW,
    pa, pb, pc,
    ad, bd, cd,
    valid);
  return pd;
}
var colors = [d3.color("DarkRed"), d3.color("DarkOrange"), d3.color("Indigo"), d3.color("purple"), d3.color("black")];
function get_colors(n, v, i) {
  return [d3.color("DarkRed"), d3.color("DarkOrange"), d3.color("Indigo"), d3.color("purple")];
}

var AM = function () {
  this.container,
  this.stats;
  this.camera;
  this.controls;
  this.scene;
  this.sceneOrtho;
  this.renderer;
  this.textureLoader;
  this.clock = new THREE.Clock();
  this.clickRequest = false;
  this.mouseCoords = new THREE.Vector2();
  this.raycaster = new THREE.Raycaster();
  this.ballMaterial = new THREE.MeshPhongMaterial({ color: 0x202020 });
  this.pos = new THREE.Vector3();
  this.quat = new THREE.Quaternion();


  this.BT_CONSTRAINT_STOP_CFM = 3;
  this.BT_CONSTRAINT_STOP_ERP = 1
  this.myCFMvalue = 0.0;
  this.myERPvalue = 0.8;

  this.jointBody = null;

  this.playgroundDimensions = {
    w: 10,
    d: 10,
    h: 3
  };
  this.GROUND_WIDTH = 1.0;

  this.gravity_on = true;
  this.margin = 0.05;

  this.armMovement = 0;

  //    this.window_height_factor = 1/4.0;
  this.window_height_factor = WINDOW_HEIGHT_FACTOR;
  // Sadly, this seems to do nothing!
  this.CAMERA_RADIUS_FACTOR = 1;

  this.grid_scene = null;
  // Used in manipulation of objects
  this.gplane = false;


  this.INITIAL_EDGE_LENGTH = TET_DISTANCE;
  this.INITIAL_EDGE_WIDTH = this.INITIAL_EDGE_LENGTH / 40;
  this.INITIAL_HEIGHT = 3 * this.INITIAL_EDGE_LENGTH / 2;
  this.NUMBER_OF_TETRAHEDRA = 70;
  //       this.NUMBER_OF_TETRAHEDRA = 5;


  this.JOINT_RADIUS = 0.09 * this.INITIAL_EDGE_LENGTH; // This is the current turret joint ball.

  this.LENGTH_FACTOR = 20;

  // Helices look like this...
  // {
  // 	helix_joints: [],
  // 	helix_members: []
  // }
  this.helices = [];



  this.meshes = [];
  this.bodies = [];


  // This is sometimes useful for debugging.
  //    this.jointGeo = new THREE.BoxGeometry( this.JOINT_RADIUS*2,this.JOINT_RADIUS*2,this.JOINT_RADIUS*2);
  this.jointGeo = new THREE.SphereGeometry(this.JOINT_RADIUS, 32, 32);
  this.jointMaterial = new THREE.MeshPhongMaterial({ color: 0x0000ff });

  this.floorTexture = new THREE.ImageUtils.loadTexture("images/logo-white-background.png");

  this.MIN_EDGE_LENGTH = this.INITIAL_EDGE_LENGTH / 2;
  this.MAX_EDGE_LENGTH = this.INITIAL_EDGE_LENGTH * 2;
  this.color_scale = d3.scale.quantile().domain([this.MIN_EDGE_LENGTH, this.MAX_EDGE_LENGTH])
    .range(['violet', 'indigo', '#8A2BE2', 'blue', 'green', 'yellow', '#FFD700', 'orange', '#FF4500']);
  this.color_material_palette = {};

  this.GROUND_PLANE_MESH;
  this.GROUND_BODY;

  this.latestLookAt = new THREE.Vector3(0, 0, 0);

  this.helix_params = [];

  // a final adjustment
  this.INITIAL_EDGE_WIDTH *= 4;
  this.JOINT_RADIUS *= 3;

}

AM.prototype.clear_non_floor_body_mesh_pairs = function () {
  this.meshes = [];
  this.bodies = [];
  this.meshes.push(am.GROUND_PLANE_MESH);
  this.bodies.push(am.GROUND_BODY);
}

var am = new AM();


var bulbLight, bulbMat, ambientLight, object, loader, stats;
var ballMat, cubeMat, floorMat;
// ref for lumens: http://www.power-sure.com/lumens.htm
var bulbLuminousPowers = {
  "110000 lm (1000W)": 110000,
  "3500 lm (300W)": 3500,
  "1700 lm (100W)": 1700,
  "800 lm (60W)": 800,
  "400 lm (40W)": 400,
  "180 lm (25W)": 180,
  "20 lm (4W)": 20,
  "Off": 0
};
// ref for solar irradiances: https://en.wikipedia.org/wiki/Lux
var hemiLuminousIrradiances = {
  "0.0001 lx (Moonless Night)": 0.0001,
  "0.002 lx (Night Airglow)": 0.002,
  "0.5 lx (Full Moon)": 0.5,
  "3.4 lx (City Twilight)": 3.4,
  "50 lx (Living Room)": 50,
  "100 lx (Very Overcast)": 100,
  "350 lx (Office Room)": 350,
  "400 lx (Sunrise/Sunset)": 400,
  "1000 lx (Overcast)": 1000,
  "18000 lx (Daylight)": 18000,
  "50000 lx (Direct Sun)": 50000
};
var params = {
  shadows: true,
  exposure: 0.68,
  bulbPower: Object.keys(bulbLuminousPowers)[4],
  hemiIrradiance: Object.keys(hemiLuminousIrradiances)[0]
};


function initGraphics() {

  am.container = document.getElementById('threecontainer');

  var PERSPECTIVE_NEAR = 0.3;

  am.camera = new THREE.PerspectiveCamera(60, window.innerWidth / (window.innerHeight * am.window_height_factor), PERSPECTIVE_NEAR, 2000);

  //   am.camera.aspect = window.innerWidth / (window.innerHeight * am.window_height_factor);

  var origin = new THREE.Vector3(0, 0, 0);
  am.camera.lookAt(origin);

  //    am.camera.quaternion.setFromAxisAngle(new THREE.Vector3(0,1,0), (Math.PI/2));

  am.scene = new THREE.Scene();
  am.scene.fog = new THREE.Fog(0x000000, 500, 10000);

  am.camera.position.x = -0.25;
  am.camera.position.y = 2.5;
  am.camera.position.z = 4;

  am.controls = new THREE.OrbitControls(am.camera, am.container);
  am.controls.target.set(0, 0, 0);

  am.renderer = new THREE.WebGLRenderer({ antialias: true });
  am.renderer.setClearColor(0xffffff);
  am.renderer.autoClearColor = true;

  am.renderer.setPixelRatio(window.devicePixelRatio);
  am.renderer.setSize(window.innerWidth, window.innerHeight * am.window_height_factor);
  am.SCREEN_WIDTH = am.renderer.getSize().width;
  am.SCREEN_HEIGHT = am.renderer.getSize().height;
  am.camera.radius = (am.SCREEN_WIDTH + am.SCREEN_HEIGHT) / this.CAMERA_RADIUS_FACTOR;


  am.cameraOrtho = new THREE.OrthographicCamera(0, am.SCREEN_WIDTH, am.SCREEN_HEIGHT, 0, - 10, 10);

  hemiLight = new THREE.HemisphereLight(0xffffff, 0xffffff, 1);
  am.scene.add(hemiLight);

  var directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
  directionalLight.position = new THREE.Vector3(100, 5, 0);
  am.scene.add(directionalLight);

  var ambientLight = new THREE.AmbientLight(0x404040);

  am.grid_scene = new THREE.Scene();
  am.grid_scene.fog = new THREE.Fog(0x000000, 500, 10000);

  // GROUND
  var groundGeo = new THREE.PlaneBufferGeometry(10000, 10000);
  var groundMat;
  groundMat = new THREE.MeshPhongMaterial({ color: 0x777777, specular: 0x050505 });

  var ground = new THREE.Mesh(groundGeo, groundMat);
  ground.name = "GROUND";
  ground.rotation.x = -Math.PI / 2;
  ground.position.y = 0;
  am.scene.add(ground);

  ground.receiveShadow = true;


  // HACK:  These diemensions are probably not right here!
  gridInit(am.grid_scene, am.playgroundDimensions);

  am.container.innerHTML = "";

  am.container.appendChild(am.renderer.domElement);

  am.sceneOrtho = new THREE.Scene();

  window.addEventListener('resize', onWindowResize, false);

}

AM.prototype.push_body_mesh_pair = function (body, mesh) {
  this.meshes.push(mesh);
  this.bodies.push(body);
}
AM.prototype.remove_body_mesh_pair = function (body, mesh) {
  for (var i = this.meshes.length - 1; i >= 0; i--) {
    if (this.meshes[i].name === mesh.name) {
      this.meshes.splice(i, 1);
      this.bodies.splice(i, 1);
    }
  }
  //    delete mesh["ammo_obj"];
  for (var i = this.rigidBodies.length - 1; i >= 0; i--) {
    if (this.rigidBodies[i].name === body.name) {
      this.rigidBodies.splice(i, 1);
    }
  }
}


function onWindowResize() {
  am.camera.aspect = window.innerWidth / (window.innerHeight * am.window_height_factor);
  am.renderer.setSize(window.innerWidth, window.innerHeight * am.window_height_factor);

  am.camera.updateProjectionMatrix();
  am.SCREEN_WIDTH = am.renderer.getSize().width;
  am.SCREEN_HEIGHT = am.renderer.getSize().height;
  am.camera.radius = (am.SCREEN_WIDTH + am.SCREEN_HEIGHT) / this.CAMERA_RADIUS_FACTOR;

  am.cameraOrtho = new THREE.OrthographicCamera(0, am.SCREEN_WIDTH, am.SCREEN_HEIGHT, 0, - 10, 10);
}

function animate() {
  requestAnimationFrame(animate);
  render();
}

var sprite_controls = new function () {
  this.size = 50;
  this.sprite = 0;
  this.transparent = true;
  this.opacity = 0.6;
  this.colorize = 0xffffff;
  this.textcolor = "yellow";
  this.rotateSystem = true;

  this.clear = function (x, y) {
    am.sceneOrtho.children.forEach(function (child) {
      if (child instanceof THREE.Sprite) am.sceneOrtho.remove(child);
    })
  };

  this.draw_and_create = function (sprite, x, y, message) {
    var fontsize = 128;
    var ctx, texture,
        spriteMaterial,
        canvas = document.createElement('canvas');
    ctx = canvas.getContext('2d');
    ctx.font = fontsize + "px Arial";

    // setting canvas width/height before ctx draw, else canvas is empty
    canvas.width = ctx.measureText(message).width;
    canvas.height = fontsize * 1; // fontsize * 1.5

    // after setting the canvas width/height we have to re-set font to apply!?! looks like ctx reset
    ctx.font = fontsize + "px Arial";
    ctx.fillStyle = this.textcolor;
    ctx.fillText(message, 0, fontsize);

    texture = new THREE.Texture(canvas);
    texture.minFilter = THREE.LinearFilter; // NearestFilter;
    texture.needsUpdate = true;

    var spriteMaterial = new THREE.SpriteMaterial({
      opacity: this.opacity,
      color: this.colorize,
      transparent: this.transparent,
      map: texture
    });

    spriteMaterial.scaleByViewport = true;
    spriteMaterial.blending = THREE.AdditiveBlending;

    if (!sprite) {
      sprite = new THREE.Sprite(spriteMaterial);
    }

    sprite.scale.set(this.size, this.size, this.size);
    sprite.position.set(x, y, 0);

    am.sceneOrtho.add(sprite);
    return sprite;
  };
};

function render() {
  var deltaTime = am.clock.getDelta();

//  sprite_controls.clear();
  am.controls.update(deltaTime);


  // note this....
  //    am.renderer.autoClear = true;
  am.renderer.render(am.scene, am.camera);
  am.renderer.render(am.grid_scene, am.camera);
  am.renderer.autoClear = false;
  am.renderer.render(am.sceneOrtho, am.cameraOrtho);
}

function initiation_stuff() {
  // Initialize Three.js
  if (!Detector.webgl) Detector.addGetWebGLMessage();
}


function init() {
  initGraphics();
}


initiation_stuff();

init();
animate();



// Find the normal to a triangle in 3space: https://stackoverflow.com/questions/19350792/calculate-normal-of-a-single-triangle-in-3d-space
// arguments THREE.js Vector3's
function normal(a, b, c) {
  var U = b.sub(a);
  var V = c.sub(a);
  return U.cross(V);
}

function clearAm() {
  am.clear_non_floor_body_mesh_pairs();
  for (var i = am.scene.children.length - 1; i >= 0; i--) {
    var obj = am.scene.children[i];
    if (obj.type == "Mesh" && obj.name != "GROUND") {
      am.scene.remove(obj);
    }
    if ((obj.name == "HELIX") || (obj.name == "AXIS") || (obj.name == "SUPERSTRUCTURE")) {
      am.scene.remove(obj);
    }
    if ((obj.type == "PROTRACTOR_LINE") || (obj.type == "PROTRACTOR_SPHERE")) {
      am.scene.remove(obj);
    }

  }
  am.helices = [];
  am.helix_params = [];
}

function getLegalTauValues(solid) {
  var taus = [];
  switch(solid) {
  case "TETRAHEDRON":
    taus = [-(Math.PI * 2) / 3, 0, (Math.PI * 2) / 3];
    break;
  case "CUBE":
    taus = [-(Math.PI * 2) / 4, 0, (Math.PI * 2) / 4];
    break;
  case "OCTAHEDRON":
    taus = [-(Math.PI * 2) / 3, 0, (Math.PI * 2) / 3];
    break;
  case "DODECAHEDRON":
    taus = [- 2 * (Math.PI * 2) / 5, -(Math.PI * 2) / 5, 0, (Math.PI * 2) / 5, 2 * (Math.PI * 2) / 5];
    break;
  case "ICOSAHEDRON":
    // DANGER!!! This works, but I have no explanation...it is undoubtedly
    // dependent on the face that we choose. At present, I only allow
    // one face to be chosen, but this means that the legal values are probably
    // dependent on that face.
    const fudge_factor = 75 * Math.PI / 180;
    taus = [(-(Math.PI * 2) / 3) + fudge_factor, fudge_factor, ((Math.PI * 2) / 3) + fudge_factor];
    break;
  }
  return taus;
}

function getLegalTauValuesX(solid,face) {
  var taus = [];
  const n = numFaces(solid);
  [numrots,base,delta] = twistBaseIncrement(solid,face);
  const values = [numrots];
  for(var i = 0; i < numrots; i++) {
    const degrees = base + ((i - (Math.floor(numrots/2))))*delta;
    values[i] = degrees * Math.PI / 180;
  }
  return values;
}

function updateLegalTauValues(solid,face) {
  var taus = getLegalTauValuesX(solid,face).map(r => r*180/Math.PI);
  if (taus.length == 3) {
    $("#radio-t0-l").html("N/A");
    $("#radio-t1-l").html(format_num(taus[0],0));
    $("#radio-t2-l").html(format_num(taus[1],0));
    $("#radio-t3-l").html(format_num(taus[2],0));
    $("#radio-t4-l").html("N/A");

    $("#radio-t4-l").hide();
    $("#radio-t0-l").hide();
  } else if (taus.length == 4) {
    $("#radio-t0-l").html(format_num(taus[0],0));
    $("#radio-t1-l").html(format_num(taus[1],0));
    $("#radio-t2-l").html(format_num(taus[2],0));
    $("#radio-t3-l").html(format_num(taus[3],0));
    $("#radio-t4-l").html("N/A");
    $("#radio-t4-l").hide();
  } else {
    $("#radio-t0-l").html(format_num(taus[0],0));
    $("#radio-t1-l").html(format_num(taus[1],0));
    $("#radio-t2-l").html(format_num(taus[2],0));
    $("#radio-t3-l").html(format_num(taus[3],0));
    $("#radio-t4-l").html(format_num(taus[4],0));

    $("#radio-t4-l").show();
    $("#radio-t0-l").show();
  }
  return taus;
}


function getSelectedTaus(taus) {
  var modex = $("input:radio[name=radio-2]:checked").attr('id');
  if (!modex) {
    modex = "0";
  }
  let n = modex.substr(-1);
  var v;
  if (taus.length == 3) {
    if (n == 0)
      n++;
    if (n == 4)
      n--;
    n = n - 1;
  }
  return taus[n];
}
function setSelectedTau(tau_radians) {

  var f = format_num(tau_radians * 180/ Math.PI,0);
  for(var i = 0; i < 5; i++) {
    var id = "#radio-t"+i+"-l";
    var html = $(id).html();
    if (f == html) {
      $(id).prop("checked", true).change();
    }
  }

}

function onComputeDelix() {
  const SOLID = getPlatonicSolidInput();
  const tau_v = TAU_d * Math.PI / 180;
  RenderSegmentedHelix(SOLID,tau_v);
 }

function main() {
  $( "input[type='radio']" ).checkboxradio();
  $('fieldset[id="superstructure"] input').change(function () {
    var modex = $("input:radio[name=radio-2]:checked")
    modex.removeAttr('checked');
    modex.button( "refresh" )

    var SOLID = getPlatonicSolidInput();
    if (SOLID) {
      const face = $( "#face-spinner" ).spinner( "value" );
      var taus = updateLegalTauValues(SOLID,face);
      tau_v = getSelectedTaus(taus);
      tau_v = tau_v * Math.PI / 180;
    }
    if (SOLID == null) {
      $("#tau-fieldset").hide();
      $("#face-spinner-span").hide();
      $("#platonic_parameters").hide();
    } else {
      $("#tau-fieldset").show();
      $("#face-spinner-span").show();
      $("#platonic_parameters").show();
    }
  });

  $('fieldset[id="tau-fieldset"] input').change(function () {
    var SOLID = getPlatonicSolidInput();
    const face = $( "#face-spinner" ).spinner( "value" );
    var taus = updateLegalTauValues(SOLID,face);
    var tau_v = getSelectedTaus(taus);
    tau_v = tau_v * Math.PI / 180;
    setTauValue(tau_v);
  });

  $('fieldset input').change(function () {

    onComputeDelix();
  });
}



// Render a Helix of radius r, with theta, v is the vector
// The helix is parallel to the vector v.
// The helix is centered on the y axis, and the two points
// at n = -1, n = 0, or centered on the z axis.
var gmat = new THREE.LineBasicMaterial({color: "green"});
function RenderHelix(l,r,d,theta,v,phi,wh,MAX_POINTS) {
  // One way to effect this is to compute a z-axis aligned helix,
  // Then rotate it parallel to v, then translate it on the
  // z axis so that the certain points on the on the z-axis.
  // In fact the rotation is purely about the y-axis.

  var init_y  = r * Math.cos(0.5*theta);
  var trans = new THREE.Matrix4().makeTranslation(0,wh - init_y,0);
  var points3D = new THREE.Geometry();
  // We'll tack on some extra segments to make it look better.
  let POINTS =  (2 + Math.floor(MAX_POINTS / 2)) * 2;
//  console.log("d,theta,phi",d,theta,phi);
  const travel = d;
  const sidedness = (2 * Math.PI / theta);
  const pitch = sidedness * travel;
//  console.log("TIGHTNESS: ",pitch/r);
  for (var i=0; i < POINTS; i++) {
    var n = i - (POINTS/2) + 0.5;
    var y = r * Math.cos(n*theta);
    // Not entirely sure why this is negated...
    var x = r * Math.sin(n*theta);
    x = -x;
    var z = n * d;
    // We will apply the global translation here...
    var p = new THREE.Vector3(x,y,z);
    p.applyMatrix4(trans);
    points3D.vertices.push(p);
  }
  var line2 = new THREE.Line(points3D, gmat);
  line2.rotation.y = phi;
  line2.name = "HELIX";
  am.scene.add(line2);
}



function set_outputs(radius,theta,travel,phi,ScaleToUnitEdgeLength,tightness_stats) {
  $( "#radius_output" ).val( format_num(radius,4) );
  $( "#theta_output" ).val( format_num(theta * 180 / Math.PI,4));
  $( "#travel_output" ).val( format_num(travel,4) );
  $( "#phi_output" ).val( format_num(phi * 180 / Math.PI,4) );

  const [[max_tightness,max_tightness_tau],
         [min_tightness,min_tightness_tau]] = tightness_stats;

  if (theta != 0) {
    const sidedness = (2 * Math.PI / theta);
    const p = sidedness * travel;
    $( "#pitch_output" ).val( format_num(p,4) );
    $( "#sidedness_output" ).val( format_num(sidedness,4) );
    const b = travel / theta;
    const torsion = travel / (radius**2 + b**2);
    const tightness = travel / radius;
    const curvature = Math.abs(radius) / (radius**2 + b**2);
    $( "#tightness_output" ).val( format_num(tightness,4) );
    $( "#torsion_output" ).val( format_num(torsion,4) );
    $( "#curvature_output" ).val( format_num(curvature,4) );

    $( "#max_tightness_output" ).val( format_num(max_tightness,4) );
    $( "#max_tightness_tau_output" ).val( format_num(max_tightness_tau * 180 / Math.PI,4) );

    $( "#min_tightness_output" ).val( format_num(min_tightness,4) );
    $( "#min_tightness_tau_output" ).val( format_num(min_tightness_tau * 180 / Math.PI,4) );

    // If we are in the word of solids, we want
    // to express these interms of unit length.
    const s = ScaleToUnitEdgeLength;
    $( "#radius_output_plato" ).val( format_num(radius * s,4) );
    $( "#travel_output_plato" ).val( format_num(travel * s,4) );
    $( "#pitch_output_plato" ).val( format_num(p * s,4) );
    $( "#torsion_output_plato" ).val( format_num(torsion * s,4) );
    $( "#tightness_output_plato" ).val( format_num(tightness * s,4) );
    $( "#curvature_output_plato" ).val( format_num(curvature * s,4) );


  } else {
    $( "#pitch_output" ).val( "NA" );
    $( "#sidedness_output" ).val( "NA" );
    $( "#torsion_output" ).val( "NA" );
    $( "#curvature_output" ).val( "NA" );
  }

}

// Here I will attempt to do several things:
// First, to compute the 4 points corresponding to rho and omega.
// Secondly, I will compute the intrinsic parameters as I have done
// in Mathematica.
// However, the point here is to render something.
// I suppose at first I can render lines.

function format_num(num,digits) {
  return parseFloat(Math.round(num * 10**digits) / 10**digits).toFixed(digits);
}



// I'm treating a label spreat as an object having postion p,
// color c, and text t.
var A_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "A"};
var B_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "B"};
var C_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "C"};
var D_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "D"};
var PHI_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "D"};
var TAU_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "D"};
var THETA_SPRITE = { p: new THREE.Vector3(0,0,0),
                 c: "green",
                 t: "D"};
var LABEL_SPRITES = [
  A_SPRITE,
  B_SPRITE,
  C_SPRITE,
  D_SPRITE,
  PHI_SPRITE,
  TAU_SPRITE,
  THETA_SPRITE
                    ];
// Create a visual protractor betwen points A, B, C in 3space
// This should really use an ellipse curver to make a fine
// protractor. However, I will just use a straightline instead
// for now.
// obj is a sprite object to attach the label two
  function lineBetwixt(A,B,color) {
    var BApoints = new THREE.Geometry();
    BApoints.vertices.push(B.clone());
    BApoints.vertices.push(A.clone());
    var BAline = new THREE.Line(BApoints, new THREE.LineBasicMaterial({color: color,linewidth: 10}));
    am.scene.add(BAline);
    BAline.type = "PROTRACTOR_LINE";
    return BAline;
  }
  function cSphere(size,p,color) {
    var mesh = createSphere(size, p, color);
    mesh.castShadow = false;
    mesh.receiveShadow = false;
    mesh.debugObject = true;
    mesh.type = "PROTRACTOR_SPHERE";
    am.scene.add(mesh);
  }

// Create an arc centered on B, touch A and C.
// At the moment, I will insist that BA = BC.
function createArc(color,A,B,C) {
  const radius = B.distanceTo(A);
  const BA = A.clone().sub(B);
  const BC = C.clone().sub(B);
  const M = vMidPoint(A,C);
  const BM = M.clone().sub(B); // Vector from B to Midpoint M.

  const angle = BA.angleTo(BC);

  const normal = BA.cross(BC);
  normal.normalize();
  // we compute the normal so we have one vector to rotate


  const Mc = BM.clone();
  const Z = new THREE.Vector3(0,0,1);
  const X = new THREE.Vector3(1,0,0);
  let qz = new THREE.Quaternion();
  qz.setFromUnitVectors(normal,Z);
  // q is the inverse quaternion to qz..
  let q = new THREE.Quaternion();
  q.setFromUnitVectors(Z,normal);

  // qz will now rotate the "normal" into Z.
  Mc.applyQuaternion(qz);
  // Mc is now rotated such that the normal is at Z...
  // that puts Mc somewhere in the Z = 0 XY plane.

  // We we want to compute the angle around the X axis...
  const rotAboutNorm = Math.atan2(Mc.y,Mc.x);

  // now we create the actual ellipse...
  // The idea is to split the difference between the angles
  // by taking the midpoint.
  var curve = new THREE.EllipseCurve(
    0,  0,            // ax, aY
    radius, radius,           // xRadius, yRadius
    rotAboutNorm-angle/2,  rotAboutNorm + angle/2,  // aStartAngle, aEndAngle
    false,            // aClockwise
    0                 // aRotation
  );

  var points = curve.getPoints( 50 );
  var geometry = new THREE.BufferGeometry().setFromPoints( points );

  var material = new THREE.LineBasicMaterial( { color : color } );

  // Create the final object to add to the scene
  var ellipse = new THREE.Line( geometry, material );

  ellipse.quaternion = q.clone();

  // This is the inverse quaternion; it seems to be required to
  // set it this qy.
  ellipse.quaternion.setFromUnitVectors(Z,normal);

  ellipse.type = "PROTRACTOR_LINE";
  ellipse.position.copy(B);

  am.scene.add(ellipse);
}
function vMidPoint(A,B) {
    return new THREE.Vector3((A.x + B.x)/2,(A.y + B.y)/2,(A.z + B.z)/2);
}
function createProtractor(obj,prefix,color,A,B,C) {

  const size = am.JOINT_RADIUS/5;

  function cSphere(size,p,color) {
    var mesh = createSphere(size, p, color);
    mesh.castShadow = false;
    mesh.receiveShadow = false;
    mesh.debugObject = true;
    mesh.type = "PROTRACTOR_SPHERE";
    am.scene.add(mesh);
  }

  cSphere(size,A,color);
  cSphere(size,B,color);
  cSphere(size,C,color);
  lineBetwixt(B,A,color);
  lineBetwixt(B,C,color);


  // now compute intermediate points....
  // how should this work if one is very small?
  const BtoA = A.clone().sub(B);
  const BtoC = C.clone().sub(B);
  const b_to_a = BtoA.length();
  const b_to_c = BtoC.length();
  var minLength = Math.min(b_to_a,b_to_c);
  var avgMidLength = (b_to_a + b_to_c)/4;
  var lengthToDraw = Math.min(minLength,avgMidLength);
  // Now that we have the length, we'll move along our vectors
  // to create new points...
  BtoA.clampLength(lengthToDraw,lengthToDraw);
  BtoC.clampLength(lengthToDraw,lengthToDraw);
  const Ap = B.clone().add(BtoA);
  const Cp = B.clone().add(BtoC);
  cSphere(size/4.0,Ap,color);
  cSphere(size/4.0,Cp,color);

  // TODO: This sometimes fails based on order
  createArc(color,Ap,B,Cp);
//  createArc(color,Cp,B,Ap);

  const PpCp_mid = vMidPoint(Ap,Cp);

  obj.p = PpCp_mid.clone();

  // There is a problem here that this is unsigned...
  const angle_rads = BtoA.angleTo(BtoC);
  obj.t = prefix + format_num((angle_rads * 180 / Math.PI),1) + " deg";
  obj.c = color;
}


const LABEL_SPRITE_FONT_SIZE = 20;
function renderSprite(obj) {
  const FifteenSpaces = "               ";

  if (obj.s) {
    am.grid_scene.remove(obj.s);
  }

  obj.s = makeTextSprite(FifteenSpaces + obj.t,
                              {fontsize: LABEL_SPRITE_FONT_SIZE},
                              obj.c );

  obj.s.position.set(obj.p.x,obj.p.y,obj.p.z);
  am.grid_scene.add(obj.s);
}

function renderSprites() {
  LABEL_SPRITES.forEach(s => renderSprite(s));
}

function setPlatonicSolidInput(mode) {
  var id;
  switch(mode) {
  case "TETRAHEDRON":
    id = "radio-tet";
    break;
  case "CUBE":
    id = "radio-cub";
    break;
  case "OCTAHEDRON":
    id = "radio-oct";
    break;
  case "DODECAHEDRON":
    id = "radio-dod";
    break;
  case "ICOSAHEDRON":
    id = "radio-ico";
    break;
  }
  $("#"+id).prop("checked", true).trigger("click").change();
}

function getPlatonicSolidInput() {
  var mode = $(":radio:checked").attr('id');
  var SOLID;
  switch (mode) {
    case "radio-arb":
    SOLID = null;
    break;
    case "radio-tet":
    SOLID = "TETRAHEDRON";
    break;
    case "radio-cub":
    SOLID = "CUBE";
    break;
    case "radio-oct":
    SOLID = "OCTAHEDRON";
    break;
    case "radio-dod":
    SOLID = "DODECAHEDRON";
    break;
    case "radio-ico":
    SOLID = "ICOSAHEDRON";
    break;
  }
  return SOLID;
}

function computeInternal(L0,B,C,tau_v,p_i,Nb,Nc) {

  // Apparently, at present I need the
  // prism to be properly rotated before my call
  // call to ajoinPrism can produce the right rotations
  // for the computeThetaAxis that way to work...a catch22
  // This is probably not right.

  let Arot = AfromLtauNbNc(L0,tau_v,Nb,Nc,false);
  let A = Arot[0];
  // I am not sure whey this is negated...
  //  var rotation = -Arot[1];

  // The rotation here is the amount we rotated
  // the initial prism to get it into a "balanced" position.
  // Balanced means the midpoint of the project of the normals
  // is pointing straight down. So we roation by this distance
  // to form rt, that must be applied to the prism  below.
  var rotation = Arot[1];

  var rt = new THREE.Matrix4();

  rt.makeRotationAxis(new THREE.Vector3(0,0,1),rotation);

  // Why wouldn't this be -A.x?
  let D = new THREE.Vector3(-A.x,A.y,-A.z);
  resK = PointAxis(L0,A);

   // NOTE: the rotations here is translated into world
  // coordinates (+2 y upward.) This makes everything
  // terribly confusing. I need to rework this with clarity.

  [p_b,rotations] = adjoinPrism(p_i,tau_v,true,false);
  //  console.log("rotations, pre",rotations);

  // I have to use this point instead of B because my
  // rotations matrix is computed in world coordinates!!!
//  let Btrans = B.clone().applyMatrix4(GTRANS);
  let Btrans = B.clone();
  resM = computeThetaAxisFromMatrix4(L0,rotations,Btrans);
  let psi = Math.PI - Math.atan2(A.x,A.y);
  if (psi > Math.PI) psi = 2*Math.PI - psi;
  let ratio = psi/tau_v;
//  console.log(ratio,psi * 180/Math.PI,tau_v * 180/Math.PI);
  return [resK,resM,p_i,rt];
}

function RenderSegmentedHelix(solid,tau_v) {
  clearAm();

  var wfgui = document.getElementById('wireframe');
  var bcgui = document.getElementById('blendcolor');
  var wf = wfgui ? wfgui.checked : true;
  var bc = bcgui ? bcgui.checked : false;


  var resK;
  var resM;
  var r;
  var theta;
  var d;
  var phi;
  let L0 = 1;

  let B = new THREE.Vector3(0,0,-L0/2);
  let C = new THREE.Vector3(0,0,L0/2);

  var obj;
  var Bn_solid;
  var Cn_solid;
  var Nb, Nc;

  var  ScaleToUnitEdgeLength;
  if (solid) {
    const face = $( "#face-spinner" ).spinner( "value" );

    [obj,Nb,Nc,ScaleToUnitEdgeLength] =
      createZAlignedPlatonic(solid,face,B,C);

    console.assert(Nb.length() > 0.01);
    console.assert(Nc.length() > 0.01);

//    taus = getLegalTauValues(solid);
//    tau_v = taus[0];
    // now we should set the legal taus radio set..
    // now we also want to set Nb and Nc....
    setTauValue(tau_v);
    setNBValues(Nb);
    setNCValues(Nc);

  } else {
    var Nb = new THREE.Vector3(NORMAL_B_X,
                               NORMAL_B_Y,
                               NORMAL_B_Z
                              );
    var Nc = new THREE.Vector3(NORMAL_C_X,
                               NORMAL_C_Y,
                               NORMAL_C_Z
                              );
    Nb.normalize();
    Nc.normalize();
    obj = null;
  }


   GLOBAL_P0 = new AbstractPrism(
    L0,
    Nb,
    Nc,
    obj);


  var p_temp = CreatePrism(GLOBAL_P0,PRISM_FACE_RATIO_LENGTH);

  [resK,resM,p_i,rt] = computeInternal(L0,B,C,tau_v,p_temp,Nb,Nc);

  const [[max_tightness,max_tightness_tau],
          [min_tightness,min_tightness_tau]]
        = bruteForceTauComputation(Nb.clone(),Nc.clone());

  console.log("normals",Nb,Nc);

  console.log("MAX_TIGHTNESS:",max_tightness_tau * 180 / Math.PI,max_tightness);
  console.log("MIN_TIGHTNESS:",min_tightness_tau * 180 / Math.PI,min_tightness);

  // console.log("computed_tau,tau",computed_tau*180/Math.PI, tau_v*180/Math.PI);

//  console.log("SHOULD MATCH");
//  console.log(resK);
//  console.log(resM);

  // We shall place this upward, for the purpose of
  // making it easier to see...

  // Take this out, and input an instance!
  applyMatrix4ToPrism(p_i,rt);

  // GTRANS sets us into position in the world
  applyMatrix4ToPrism(p_i,GTRANS);

  // if (p_i.sup) {
  //    p_i.sup.applyMatrix(GTRANS);
  // }
  p_i.p.Nb.applyMatrix4(rt);
  p_i.p.Nc.applyMatrix4(rt);

  let Cp = resM[6];
  if (Cp) {
    var Cpa = new THREE.Vector3(Cp.x,Cp.y,Cp.z);
    cSphere(am.JOINT_RADIUS/4,Cpa,"red");
    var B_ = resM[7];
    cSphere(am.JOINT_RADIUS/3,B_,"green");
  }
  if (false) {
    console.assert(near(resK[0],resM[0]));
    console.assert(near(resK[1],resM[1]));
    if (!near(resK[1],resM[1])) {
      console.log("Kahn, Matrix",resK[1] * 180 / Math.PI, resM[1] * 180 / Math.PI);
//      debugger;
    }
    console.assert(near(resK[2],resM[2]));
    console.assert(near(resK[3],resM[3]));
    console.assert(near(resK[4],resM[4]));
    if (!(near(resK[4],resM[4]))) {
      console.log("Kahn Phi, Matrix Phi", resK[4] * 180 / Math.PI, resM[4] * 180 / Math.PI);
//      debugger;
    }
    console.assert(vnear(resK[5],resM[5]));
    if (!vnear(resK[5],resM[5])) {
      console.log("AXIS KAHN, MATRIX", resK[5],resM[5]);
    }
  }

  // TODO: Figure out why I can't use this tomorrow!
  var USE_MATRIX = false;
  var res;
  if (USE_MATRIX)
    res = resM;
  else
    res = resK;

  r = res[0];
  theta = res[1];
  d = res[2];
  phi = res[4];

  var prisms = createAdjoinedPrisms(p_i,tau_v,(NUM_SEGMENTS -1)/2);

  B = prisms[0][6].position;

  // vector pointing from B to Ba
  var Ba = res[6];
  // console.log("BA ",Ba);
  // Ba may not be defined!
  if (Ba) {
    if (!USE_MATRIX) {
      Ba.applyMatrix4(GTRANS);
    }

    // We'll put a Ball at Ba ...
    cSphere(am.JOINT_RADIUS/5,new THREE.Vector3(Ba.x,Ba.y,Ba.z),"red");
  }

  // Ba and Ca need to be on the axis, that is an assertion.
  C = prisms[0][7].position;

  B_SPRITE.p = B.clone();
  C_SPRITE.p = C.clone();

  D_SPRITE.p = prisms[1][0][7].position.clone();
  A_SPRITE.p = prisms[2][0][6].position.clone();


  if (Ba) {
    // Since we've set the first prism up symmetrically, Ca
    // mirrors Ba...
    var Ca = new THREE.Vector3(-Ba.x,Ba.y,-Ba.z);
//    var Ca = new THREE.Vector3(Ba.x,Ba.y,Ba.z);
    cSphere(am.JOINT_RADIUS/5,Ca,"blue");
  }
  // r = res[0];
  // theta = res[1];
  // d = res[2];
  // phi = res[4];
  const tightness_stats = [[max_tightness,max_tightness_tau],
                           [min_tightness,min_tightness_tau]];
  set_outputs(r,theta,d,phi,ScaleToUnitEdgeLength,tightness_stats);




  RenderHelix(L0,r,d,theta,new THREE.Vector3(0,0,1),phi,
              WORLD_HEIGHT,NUM_SEGMENTS);

  // These are the "axes" markers...
  create_vertex_mesh(new THREE.Vector3(0,0,0),d3.color("white"));
  create_vertex_mesh(new THREE.Vector3(1,0,0),d3.color("red"));
  create_vertex_mesh(new THREE.Vector3(0,1,0),d3.color("green"));
  create_vertex_mesh(new THREE.Vector3(0,0,1),d3.color("blue"));

  // now we would like to draw the axis of the helix...
  // we have the vector H from the PointAxi algorithm.
  // We have to find one point on the helix---
  // We know the helix intersects the y axis,
  // and we have the radius, which is the distance
  // to the joints.
  // First, let me just draw one at the origin
  // in the correct direction
  var points3D = new THREE.Geometry();

  H = res[5].clone();
  // we compute y via Pythagoras from the a line
  // from the y-axis to a joint--- yd is always slightly
  // less than radius because it is the distance to
  // the midpoint of a segment.
  let Qsq = r**2 + (d/2)**2 - (L0/2)**2;
  if (near(Qsq,0,1e-4)) {
    Qsq = 0;
  }
  let yd = Math.sqrt(Qsq);
  // I unfortunately have some kind of sign error here...
  const FACTOR = 6;
  H.multiplyScalar(FACTOR);
  H.setY(WORLD_HEIGHT - yd);
  // WTF?
  H.setX(H.x);
  var Hn = H.clone(); // H "negative", not H normal!

  Hn.multiplyScalar(-1);
  Hn.setY(H.y);


  points3D.vertices.push(H);
  points3D.vertices.push(Hn);

  var axis_line = new THREE.Line(points3D, new THREE.LineBasicMaterial({color: "green",linewidth: 10}));
  axis_line.name = "AXIS";
  am.scene.add(axis_line);

  var hex = 0x008000;

  let Hnmlzd = H.clone().sub(Hn).normalize();

  var arrowHelper = new THREE.ArrowHelper( Hnmlzd, Hn, 2*FACTOR, hex,0.5,0.15 );
  arrowHelper.name = "AXIS";
  am.scene.add(arrowHelper);


  // Technically, we could draw the Theta protractor but we would
  // have to extend the lines
  // TODO: My Theta protract is wrong when tau = 180 degrees.
  if (Ba) {
    // Now we will attempt to render the B-BA line...
    lineBetwixt(B,Ba,"red");
    lineBetwixt(C,Ca,"green");

    // Now, in order to be able todraw the theta
    // protractor, we will translate Ca to Cpara in the -H
    // direction to place it in a circle at Ba, then
    // add a protractor between them.p
    var Cpara = C.clone();
    var Hdir = Hnmlzd.clone().clampLength(d,d);
    Cpara.sub(Hdir);

    cSphere(am.JOINT_RADIUS/5,Cpara,"green");
    // a nice greenline parallel to the helix axis should help..

    lineBetwixt(C,Cpara,"back");
    // here I attempt to create the visually important
    // theta protractor
    // A test here is that the cross of the vectors
    // should be equal axis vector...
    {
      // This is just a test....
      const BBa = new THREE.Vector3().subVectors(B,Ba);
      const CparaBa = new THREE.Vector3().subVectors(Cpara,Ba);
      const perpToProtractor = new THREE.Vector3().crossVectors(BBa,CparaBa);
      const crossPerpVH = new THREE.Vector3().crossVectors(perpToProtractor,Hdir);
      console.assert(vnear(crossPerpVH,new THREE.Vector3()));
    }

    {
      createProtractor(THETA_SPRITE,"theta = ","black",B,Ba,Cpara);
    }
  }

  let O = new THREE.Vector3(0,0,0);
  {
    let Z = new THREE.Vector3(0,0,1);
    let Hyplane = new THREE.Vector3(H.x,0,H.z);
    Hyplane.clampLength(1,1);
    createProtractor(PHI_SPRITE,"phi = ","purple",Z,O,Hyplane);
  }

  // here I attempt to create the visually important
  // tau protractor
  {
    // tau comes from the prisms, the center of the
    // joint is just B
    // The other two points are corresponding points
    // on at the joint face. We'll use the TOP elements.
    let Bface = prisms[0][3].position;
    let Cface = prisms[1][0][0].position;
    createProtractor(TAU_SPRITE,"tau = ","green",Bface,C,Cface);
  }

  renderSprites();
}




function addDebugSphere(am,pos,color) {
  if (!color) {
    color = "yellow";
  }
  var mesh = createSphere(am.JOINT_RADIUS/5, pos, color);
  mesh.castShadow = false;
  mesh.receiveShadow = false;
  mesh.debugObject = true;
  am.scene.add(mesh);
}

var ANGLE_RHO_d = 0;

$(function() {
  $( "#angle_rho_slider" ).slider({
    range: "max",
    min: 0,
    max: 180,
    value: ANGLE_RHO_d,
    step: 0.001,
    slide: function( event, ui ) {
      $( "#angle_rho" ).val( ui.value );
      ANGLE_RHO_d = ui.value;
      console.log(ANGLE_RHO_d);
      onComputeDelix();
    }
  });
  $( "#angle_rho" ).val( $( "#angle_rho_slider" ).slider( "value" ) );
});

var ANGLE_OMEGA_d = 0;

$(function() {
  $( "#angle_omega_slider" ).slider({
    range: "max",
    min: 0,
    max: 180,
    value: ANGLE_OMEGA_d,
    step: 0.001,
    slide: function( event, ui ) {
      $( "#angle_omega" ).val( ui.value );
      ANGLE_OMEGA_d = ui.value;
      console.log(ANGLE_OMEGA_d);
      onComputeDelix();
    }
  });
  $( "#angle_omega" ).val( $( "#angle_omega_slider" ).slider( "value" ) );
});


var NORMAL_B_X = -0.5;
var NORMAL_B_Y = -0.5;
var NORMAL_B_Z = -1;
var NORMAL_C_X = 0;
var NORMAL_C_Y = 0;
var NORMAL_C_Z = 1;
var TAU_d = -10;
var TAU = 0;

{

  function setTauValue(tau) {
    TAU_d = format_num(tau * 180 / Math.PI,3);
    $( "#tau_slider" ).slider( "value", TAU_d );
    $( "#tau_txt" ).val( TAU_d );
    $( "#tau_d" ).val( TAU_d );
  }

  $(function() {
    $( "#tau_slider" ).slider({
      range: "max",
      min: -180,
      max: 180,
      value: TAU_d,
      step: 0.01,
      slide: function( event, ui ) {
	$( "#tau_txt" ).val( ui.value );
	$( "#tau_d" ).val( ui.value );
	TAU_d = ui.value;
        console.log(TAU_d);
        onComputeDelix();
      }
    });
    $( "#tau_d" ).val( $( "#tau_slider" ).slider( "value" ) );
  });

  function setNBValues(Nb) {
    NORMAL_B_X =  Nb.x;
    NORMAL_B_Y =  Nb.y;
    NORMAL_B_Z =  Nb.z;
    $( "#normal_b_x_slider" ).slider( "value", NORMAL_B_X );
    $( "#b_x" ).val( format_num(NORMAL_B_X,3) );
    $( "#normal_b_x" ).val( format_num(NORMAL_B_X,3) );

    $( "#normal_b_y_slider" ).slider( "value", NORMAL_B_Y );
    $( "#b_y" ).val( format_num(NORMAL_B_Y,3) );
    $( "#normal_b_y" ).val( format_num(NORMAL_B_Y,3) );

    $( "#normal_b_z_slider" ).slider( "value", NORMAL_B_Z );
    $( "#b_z" ).val( format_num(NORMAL_B_Z,3) );
    $( "#normal_b_z" ).val( format_num(NORMAL_B_Z,3) );
  }
  function setNCValues(Nc) {
    NORMAL_C_X =  Nc.x;
    NORMAL_C_Y =  Nc.y;
    NORMAL_C_Z =  Nc.z;
    $( "#normal_c_x_slider" ).slider( "value", NORMAL_C_X );
    $( "#c_x" ).val( format_num(NORMAL_C_X,3) );
    $( "#normal_c_x" ).val( format_num(NORMAL_C_X,3) );

    $( "#normal_c_y_slider" ).slider( "value", NORMAL_C_Y );
    $( "#c_y" ).val( format_num(NORMAL_C_Y,3) );
    $( "#normal_c_y" ).val( format_num(NORMAL_C_Y,3) );

    $( "#normal_c_z_slider" ).slider( "value", NORMAL_C_Z );
    $( "#c_z" ).val( format_num(NORMAL_C_Z,3) );
    $( "#normal_c_z" ).val( format_num(NORMAL_C_Z,3) );
  }
  $(function() {
    $( "#normal_b_x_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_B_X,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#b_x" ).val( ui.value );
	$( "#normal_b_x" ).val( ui.value );
	NORMAL_B_X = ui.value;
        console.log(NORMAL_B_X);
        onComputeDelix();
      }
    });
    $( "#b_x" ).val( $( "#normal_b_slider" ).slider( "value" ) );
  });

  $(function() {
    $( "#normal_b_y_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_B_Y,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#b_y" ).val( ui.value );
	$( "#normal_b_y" ).val( ui.value );
	NORMAL_B_Y = ui.value;
        console.log(NORMAL_B_Y);
        onComputeDelix();
      }
    });
    $( "#b_y" ).val( $( "#normal_b_slider" ).slider( "value" ) );
  });

  $(function() {
    $( "#normal_b_z_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_B_Z,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#b_z" ).val( ui.value );
	$( "#normal_b_z" ).val( ui.value );
	NORMAL_B_Z = ui.value;
        console.log(NORMAL_B_Z);
        onComputeDelix();
      }
    });
    $( "#b_z" ).val( $( "#normal_b_slider" ).slider( "value" ) );
  });
}

{
  $(function() {
    $( "#normal_c_x_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_C_X,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#c_x" ).val( ui.value );
	$( "#normal_c_x" ).val( ui.value );
	NORMAL_C_X = ui.value;
        console.log(NORMAL_C_X);
        onComputeDelix();
      }
    });
    $( "#c_x" ).val( $( "#normal_c_slider" ).slider( "value" ) );
  });

  $(function() {
    $( "#normal_c_y_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_C_Y,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#c_y" ).val( ui.value );
	$( "#normal_c_y" ).val( ui.value );
	NORMAL_C_Y = ui.value;
        console.log(NORMAL_C_Y);
        onComputeDelix();
      }
    });
    $( "#c_y" ).val( $( "#normal_c_slider" ).slider( "value" ) );
  });

  $(function() {
    $( "#normal_c_z_slider" ).slider({
      range: "max",
      min: -1,
      max: 1,
      value: NORMAL_C_Z,
      step: 0.001,
      slide: function( event, ui ) {
	$( "#c_z" ).val( ui.value );
	$( "#normal_c_z" ).val( ui.value );
	NORMAL_C_Z = ui.value;
        console.log(NORMAL_C_Z);
        onComputeDelix();
      }
    });
    $( "#c_z" ).val( $( "#normal_c_slider" ).slider( "value" ) );
  });
}


const INIT_RHO = 10;
const INIT_OMEGA = 55;

function setup_input_molecule(slider,ro,txt,x,set)
{
  $( slider ).slider( "value",x );
  $( ro ).val( x );
  $( txt ).val( "" );

  $( txt ).keypress(function(event) {
    if (event.which == 13) {
      // Does this change the value or the parameter?
      x = event.currentTarget.value;
      set(x);
      $( slider ).slider( "value",x );
      $( ro ).val( x );
      onComputeDelix();
    }
  });
}

function numFaces(solid) {
  switch(solid) {
  case "TETRAHEDRON":
    return 4;
  case "CUBE":
    return 6;
  case "OCTAHEDRON":
    return 8;
  case "DODECAHEDRON":
    return 12;
  case "ICOSAHEDRON":
    return 20;
  }
}

function twistBaseIncrement(solid,face) {
  switch(solid) {
  case "TETRAHEDRON":
    return [3,0,120];
  case "CUBE":
    return [4,0,90];
  case "OCTAHEDRON": {
    const map_face = [];
    map_face[1] = 0;
    map_face[2] = 0;
    map_face[3] = 0;
    map_face[4] = 0;
    map_face[5] = 60;
    map_face[6] = 0;
    map_face[7] = 0;
    return [3,map_face[face],120];
  }
  case "DODECAHEDRON": {
    const map_face = [];
    map_face[1] = 0;
    map_face[2] = 0;
    map_face[3] = 0;
    map_face[4] = 0;
    map_face[5] = 0;
    map_face[6] = 0;
    map_face[7] = 0;
    map_face[8] = 36; // opposite face
    map_face[9] = 0;
    map_face[10] = 0;
    map_face[11] = 0;
    return [5,map_face[face],72];
  }
  case "ICOSAHEDRON":
    // Thewe were developed by hand;
    // the represent the angles for
    // exact face-to-face match (with +- 120 degrees)
    const map_face = [];
    map_face[1] = 120;
    map_face[2] = 120;
    map_face[3] = 120;
    map_face[4] = 120;
    map_face[5] = 120;
    map_face[6] = 120;
    map_face[7] = 120;
    map_face[8] = 45;
    map_face[9] = 45;
    map_face[10] = 45;
    map_face[11] = 45;
    map_face[12] = 120;
    map_face[13] = 60; // opposite face
    map_face[14] = 120;
    map_face[15] = 120;
    map_face[16] = 120;
    map_face[17] = 45;
    map_face[18] = 120;
    map_face[19] = 45;
    return [3,map_face[face],120];
  }
}

function renderhelixrow(solid,face,tau) {
  setPlatonicSolidInput(solid);
  updateLegalTauValues(solid,face);
  setSelectedTau(tau);
  $( "#face-spinner" ).spinner( "value",face );
  RenderSegmentedHelix(solid,tau);
}


function hideClassNum(solid,cnum,visibility) {
  var table = document.getElementById("platonichelices");
  var cnt = 0;
  for (var i = 0, row; row = table.rows[i]; i++) {
    //iterate through rows
    //rows would be accessed using the "row" variable assigned in the for loop
    var class_num = row.getAttribute("class_num");
    var solid_name = row.getAttribute("solid_name");

    if ((solid == solid_name) && (class_num == cnum)) {
      if (cnt != 0) {
        row.style.display = visibility;
      }
      cnt++;
    }
  }
}

function registerHelix(name,number,solid_num,analogs,face,tau,radius,theta,travel,helix_angle,class_num) {
  var table = document.getElementById("platonichelices");

  var row = table.insertRow(-1);

  row.id = number;
  row.setAttribute("class_num", class_num);
  row.setAttribute("solid_name", name);

  var cnt = 0;
  var button_c = row.insertCell(cnt++);
  var totalnum_c = row.insertCell(cnt++);
  var solidnum_c = row.insertCell(cnt++);
  var name_c = row.insertCell(cnt++);
  var analogs_c = row.insertCell(cnt++);
  var face_c = row.insertCell(cnt++);
  var tau_c = row.insertCell(cnt++);
  var radius_c = row.insertCell(cnt++);
  var theta_c = row.insertCell(cnt++);
  var travel_c = row.insertCell(cnt++);
  var angle_c = row.insertCell(cnt++);
  var class_c = row.insertCell(cnt++);

  button_c.innerHTML = "<button onclick='renderhelixrow(\""+
    name+"\","+face+","+
    tau+
  ")'>Draw</button>";
  totalnum_c.innerHTML = number;
  solidnum_c.innerHTML = solid_num;
  name_c.innerHTML = name;

  analogs_c.innerHTML = analogs;
  face_c.innerHTML = face;
  tau_c.innerHTML = format_num(tau * 180/Math.PI,0);
  radius_c.innerHTML = format_num(radius,3);
  theta_c.innerHTML = format_num(theta * 180 / Math.PI,3);
  travel_c.innerHTML = format_num(travel,3);
  angle_c.innerHTML = format_num(helix_angle  * 180 / Math.PI,3);
  var button1 = "<button class='collapse' onclick='hideClassNum(\""+name+"\","+class_num+",\"none\")'>Collapse "+class_num+"</button>";
  var button2 = "<button onclick='hideClassNum(\""+name+"\","+class_num+",\"table-row\")'>Expand "+class_num+"</button>";
  class_c.innerHTML = button1 + button2;
}

function countClass(measures,solid,cnum) {
  var cnt = 0;
  for(var i = 0; i < measures.length; i++) {
    if ((measures[i][8] == cnum) && (measures[i][0] == solid))  {
      cnt++;
    }
  }
  return cnt;
}

function addMeasures(measures,s,cnt,f,tau,r,theta,d,phi) {
  var class_num = -1;
  var max_class = -1;
  for(var i = 0; i < measures.length; i++) {
    if ((s == measures[i][0]) &&
        (near(r,measures[i][4])) &&
        (near(theta,measures[i][5])) &&
        (near(Math.abs(d),Math.abs(measures[i][6])))
       ) {
      class_num = measures[i][8];
    }
    if ((measures[i][8] > max_class)) max_class = measures[i][8];
  }
  measures.push([s,cnt,f,tau,r,theta,d,phi,(class_num >= 0) ? class_num : max_class + 1 ]);
}
// Our goal here is to create a list of the classe, which
// elements belongs to which, how many are in which, and wether they contain both chiralities
// (I think they will have to.
function analyzeClasses(mAndC) {
  let analysis = [];
  const len = mAndC.length;
  for(var i = 0; i < len; i++) {
    const c = mAndC[i][8];
    const s = Math.sign(mAndC[i][6]);
    if (analysis[c] != null) {
      const v = analysis[c];
      const h = (v[2] != s) ? "B" : s;
      analysis[c] = [v[0]+1,[...v[1],i],h];
    } else {
      analysis[c] = [1,[i],s];
    }
  }
  return analysis;
}
function populatePlatonicHelixTable() {
  var cnt = 0;
  var rownum = 0;
  for(var s of PLATONIC_SOLIDS) {
    const n = numFaces(s);
    let measuresAndClasses = [];
    var scnt = 0;
    // We don't operated on the zero face, which is technically a weakness on my part....
    for(var f = 1; f < n; f++ ) {
      [numrots,base,delta] = twistBaseIncrement(s,f);

      let L0 = 1;
      let B = new THREE.Vector3(0,0,-L0/2);
      let C = new THREE.Vector3(0,0,L0/2);

      [obj,Nb,Nc,ScaleToUnitEdgeLength] = createZAlignedPlatonic(s,f,B,C);

      // now we loop over num_rotations...
      for(var j = 0; j < numrots; j++) {
        const tau = (base + ((j - (Math.floor(numrots/2))))*delta) * Math.PI / 180;
        // unsure about this.....
        var p_temp = CreatePrism(GLOBAL_P0,PRISM_FACE_RATIO_LENGTH);
        [resK,resM,p_i,rt] = computeInternal(L0,B,C,tau,p_temp,Nb,Nc);
        var USE_MATRIX = false;
        var res = (USE_MATRIX) ? resM : resK;
// Note that here we apply the scaling factor...
        const r = res[0] * ScaleToUnitEdgeLength;
        const theta = res[1];
        const d = res[2] * ScaleToUnitEdgeLength;
        const phi = res[4];
        addMeasures(measuresAndClasses,s,scnt,f,tau,r,theta,d,phi);
        scnt++
      }
    }
    var analysis = analyzeClasses(measuresAndClasses);
    analysis.sort((a,b) =>
                  {
                    // We will look at exexmplary members of the classes
                    [sa,cnta,fa,taua,ra,thetaa,da,phia,cna] = measuresAndClasses[a[1][0]];
                    [sb,cntb,fb,taub,rb,thetab,db,phib,cnb] = measuresAndClasses[b[1][0]];

                    return (ra == rb) ? fa - fb : ra - rb;
                  }
                 );

    for (var a of analysis) {
      for (var idx of a[1]) {
        [s,cnt,f,tau,r,theta,d,phi,cn] = measuresAndClasses[idx];
        var analogs = countClass(measuresAndClasses,s,cn);
        registerHelix(s,rownum,idx,analogs,f,tau,r,theta,d,phi,cn);
        rownum++;
      }
    }

    console.log(s);
    console.log(analysis);
  }
}

function createZAlignedPlatonic(solid,face,B,C) {
  let Bf = 0;
  // We disallow using face 0. Altought this is mathematically
  // possible, it is the ultimate self-intersection and not
  // terribly interesting, and major pain to try to handle
  // algorithmically
  face = Math.max(1,face);
  face = Math.min(numFaces(solid)-1,face);
  switch(solid) {
  case "TETRAHEDRON":
    face = Math.min(3,face);
    return createZAlignedPlatonicAux(solid,B,C,Bf,face);
  case "CUBE":
    face = Math.min(5,face);
    return createZAlignedPlatonicAux(solid,B,C,Bf,face);
  case "OCTAHEDRON":
    face = Math.min(7,face);
    return createZAlignedPlatonicAux(solid,B,C,Bf,face);
  case "DODECAHEDRON":
    face = Math.min(11,face);
    return createZAlignedPlatonicAux(solid,B,C,Bf,face);
  case "ICOSAHEDRON":
    face = Math.min(19,face);
    return createZAlignedPlatonicAux(solid,B,C,Bf,face);
  }
}

// This is subtle and needs to be renamed.
// Using the "centroid" as the "up" vector works well
// except in some configurations when it is zero (that is,
// the face are perfectly opposed and parallel.
// In that case we need to choose the center of a DIFFERENT face.

function findCentroid(geo) {
  var sum = new THREE.Vector3(0,0,0);
  var n = 0;
  geo.vertices.forEach(v => { sum.add(v); n++; });
  var cent = sum.multiplyScalar(1/n);
  if (vnear(cent,new THREE.Vector3(0,0,0))) {
    return new THREE.Vector3(0,1,0);
  } else {
    return cent;
  }
}

function findFaceCentroid(geo,face) {
  var sum = new THREE.Vector3(0,0,0);
  var n = 0;
  sum.add(geo.vertices[face.a]);
  sum.add(geo.vertices[face.b]);
  sum.add(geo.vertices[face.c]);
  var cent = sum.multiplyScalar(1/3);
  return cent;
}
function averageVectors(vs) {
  var sum = new THREE.Vector3(0,0,0);
  const n = vs.length;
  for(var i = 0; i < n; i++) {
    sum.add(vs[i]);
  }
  var avg = sum.multiplyScalar(1/n);
  return avg;

}

// This routine needs to return not only the object
// for rending, but also the normal vectors!
// And legal tau values!
function createZAlignedPlatonicAux(solid, B, C, Bf, Cf) {
  // Each of the platonic solids actually requires you to also
  // make a choice as to the face...
  var geo;
  switch(solid) {
  case "TETRAHEDRON":
    geo = new THREE.TetrahedronGeometry(1,0);
    break;
  case "CUBE":
    geo = new THREE.CubeGeometry(1,0);
    break;
  case "OCTAHEDRON":
    geo = new THREE.OctahedronGeometry(1,0);
    break;
  case "DODECAHEDRON":
    geo = new THREE.DodecahedronGeometry(1,0);
    break;
  case "ICOSAHEDRON":
    geo = new THREE.IcosahedronGeometry(1,0);
    break;
  }
  // This probably doesn't work when the faces aren't triangles...
  // we need to choose the vertices more generally.
  // In particualr for a cube and dodecahedron this is likely to be wrong.

  let Bc = new THREE.Vector3();
  let Cc = new THREE.Vector3();

  // The are LOCAL face normals....
  var Bnl;
  var Cnl;
  // There are a fixed set of legal "tau" values for platonic helices

  if ((solid != "CUBE") && (solid != "DODECAHEDRON")) {
    Bc.add(geo.vertices[geo.faces[Bf].a]);
    Bc.add(geo.vertices[geo.faces[Bf].b]);
    Bc.add(geo.vertices[geo.faces[Bf].c]);

    Bc.multiplyScalar(1/3);
    Bnl = geo.faces[Bf].normal.clone();

    Cc.add(geo.vertices[geo.faces[Cf].a]);
    Cc.add(geo.vertices[geo.faces[Cf].b]);
    Cc.add(geo.vertices[geo.faces[Cf].c]);
    Cc.multiplyScalar(1/3);
    Cnl = geo.faces[Cf].normal.clone();
  } else if (solid == "CUBE") {
    // here we know we are the midpoint of opposite vertices...
    // the concept of "faces" for a cube and a dodecahedron
    // does not match the geometry, which uses multiple triangles.
    // So it is easier to the mapping our selves...
    // Although technically for a cube you could choose
    // four faces, the are all symmetric, so we will just
    // wire one in and not worry about it.
    // This is technically an error, and does not support the
    // fifth face (which is a very boring helix, having parallel
    // normals, but so be it. I will capture this by aborting if
    // you don't choose 0,1 for the faces.

    const vFF = []; // vertices crossing this face
    // How do I figure out how to fill this in?
    vFF[0] = [0,3];
    vFF[1] = [4,7];
    vFF[2] = [0,4];
    vFF[3] = [3,7];
    vFF[4] = [2,5];
    vFF[5] = [1,6];

    const geoF = []; // map our face into triangle normals
    geoF[0] = 0;
    geoF[1] = 3;
    geoF[2] = 5;
    geoF[3] = 7;
    geoF[4] = 9;
    geoF[5] = 11;

    Bc.add(geo.vertices[vFF[Bf][0]]);
    Bc.add(geo.vertices[vFF[Bf][1]]);
    Bc.multiplyScalar(1/2);

    Cc.add(geo.vertices[vFF[Cf][0]]);
    Cc.add(geo.vertices[vFF[Cf][1]]);
    Cc.multiplyScalar(1/2);

    Bnl = geo.faces[geoF[Bf]].normal.clone();
    Cnl = geo.faces[geoF[Cf]].normal.clone();
  } else if (solid == "DODECAHEDRON") {
    // Now, sadly, we really have to just know the vertices of
    // the first face...this numbering comes from THREE and is
    // effectively arbitrary as far as we are concerned.
    const vFF = []; // vertices crossing this face
    // How do I figure out how to fill this in?
    vFF[0] = [0,1,2,3,4];
    vFF[1] = [1,3,5,6,7];
    vFF[2] = [6,7,8,9,10];
    vFF[3] = [9,11,12,13,10];
    vFF[4] = [11,12,14,15,16];
    vFF[5] = [3,4,7,10,13];
    vFF[6] = [2,4,12,13,16];
    vFF[7] = [0,2,15,16,17];
    vFF[8] = [8,9,11,14,18];
    vFF[9] = [0,1,5,17,19];
    vFF[10] = [5,6,8,18,19];
    vFF[11] = [14,15,17,18,19];

    const geoF = []; // map our face into triangle normals
    geoF[0] = 0;
    geoF[1] = 3;
    geoF[2] = 6;
    geoF[3] = 9;
    geoF[4] = 12;
    geoF[5] = 15;
    geoF[6] = 18;
    geoF[7] = 21;
    geoF[8] = 24;
    geoF[9] = 27;
    geoF[10] = 30;
    geoF[11] = 33;

    Bc.add(geo.vertices[vFF[Bf][0]]);
    Bc.add(geo.vertices[vFF[Bf][1]]);
    Bc.add(geo.vertices[vFF[Bf][2]]);
    Bc.add(geo.vertices[vFF[Bf][3]]);
    Bc.add(geo.vertices[vFF[Bf][4]]);
    Bc.multiplyScalar(1/5);

    Cc.add(geo.vertices[vFF[Cf][0]]);
    Cc.add(geo.vertices[vFF[Cf][1]]);
    Cc.add(geo.vertices[vFF[Cf][2]]);
    Cc.add(geo.vertices[vFF[Cf][3]]);
    Cc.add(geo.vertices[vFF[Cf][4]]);
    Cc.multiplyScalar(1/5);

    Bnl = geo.faces[geoF[Bf]].normal.clone();
    Cnl = geo.faces[geoF[Cf]].normal.clone();
  }
  console.assert(Bnl.length() > 0.1);
  console.assert(Cnl.length() > 0.1);

  let O = new THREE.Vector3();

  // The problem is when the faces are the same, this
  // becomes zero...I think this
  var d = Cc.clone();
  d.sub(Bc);
  var olen = d.length();
  let bclen = C.clone().sub(B).length();
  let scale_m = new THREE.Matrix4().identity();
  // If the olen is zero, which means they are same face,
  // we cannot determine scale in this way, so we will default to 1.
  let s = (olen == 0) ? 1 : bclen/olen;
  // we could do this via world transformations, but
  // we are attempting to build a locally aligned object...
  scale_m.makeScale(s,s,s);
  geo.scale(s,s,s);
  const SCALE_TO_UNIT_EDGE = 1/s;

  Bc.applyMatrix4(scale_m);
  Cc.applyMatrix4(scale_m);


  d = Cc.clone();
  d.sub(Bc);

  let F = Cc.clone();
  let T = Bc.clone();
  const Z = new THREE.Vector3(0,0,1);
  Z.normalize();
  d.normalize();

  let rotation = new THREE.Matrix4().identity();
   let q = new THREE.Quaternion();
  q.setFromUnitVectors(d,Z);
//   console.log(Z,d);
  rotation.makeRotationFromQuaternion(q);
  // Note: This is a local transformation!
  geo.applyMatrix(rotation);
  Bc.applyMatrix4(rotation);
  Cc.applyMatrix4(rotation);
  console.assert(Bnl.length() > 0.1);
  console.assert(Cnl.length() > 0.1);

  Bnl.applyMatrix4(rotation);
  Cnl.applyMatrix4(rotation);
  O.applyMatrix4(rotation);
  console.assert(Bnl.length() > 0.1);
  console.assert(Cnl.length() > 0.1);


  let trans = B.clone().sub(Bc);
  let trans_m = new THREE.Matrix4().makeTranslation(trans.x,trans.y,trans.z);

  // Note: This is a local tranformation
  geo.applyMatrix(trans_m);
//  arrowHelper.applyMatrix(trans_m);
  Bc.applyMatrix4(trans_m);
  Cc.applyMatrix4(trans_m);
//  Bnl.applyMatrix4(trans_m);
//  Cnl.applyMatrix4(trans_m);
  O.applyMatrix4(trans_m);
//  console.log("Bc,Cc",Bc,Cc);

    // now we may have Bc and Cc in the right location.
  // We may have to rotate about the Z (Bc-Cc) axis to
  // make sure the up vector is correct. But how to even define?
  // ideally would use the centroid of the object, but Bc-Cc might
  // pass through them!!
  // However, we can at least try to compute the centroid and make the
  // vector from the Bc-Cc line point upwards (+Y).
  var cent = findCentroid(geo);
  if ((cent == null) || vnear(cent,new THREE.Vector3(0,0,0))) {
    switch(solid) {
    case "CUBE":
      cent = new THREE.Vector3(0,1,0);
      break;
    case "DODECAHEDRON":
      cent = findFaceCentroid(geo,geo.faces[0]);
      break;
    case "ICOSAHEDRON":
      cent = findFaceCentroid(geo,geo.faces[0]);
      break;
    default:
      console.log("INTERNAL ERROR",solid);
      break;
  }

  }

//  console.log("centroid ",cent,O);

  var oxy = new THREE.Vector2(O.x-B.x,O.y-B.y);

  // angle is measured against
  var z_angle = oxy.angle() -Math.PI/2;



  let rot_z = new THREE.Matrix4().makeRotationZ(-z_angle);

  // WARNING :: I am traslating back to the center.
  // I really need to translate the BC vector back to the
  // the Z axis. This is really worth a subroutine here!

//  console.log("rot_z",z_angle * 180 / Math.PI);

  // Now, sadly, we have to pute the object back at the axis to
  // rotate...
  let trans_down = new THREE.Matrix4().makeTranslation(0,B.y,0);
  let trans_down_i = new THREE.Matrix4().getInverse(trans_down);
  geo.applyMatrix(trans_down_i);
  Bc.applyMatrix4(trans_down_i);
  Cc.applyMatrix4(trans_down_i);
  O.applyMatrix4(trans_down_i);

  geo.applyMatrix(rot_z);
  Bc.applyMatrix4(rot_z);
  Cc.applyMatrix4(rot_z);
  Bnl.applyMatrix4(rot_z);
  Cnl.applyMatrix4(rot_z);
  O.applyMatrix4(rot_z);

  geo.applyMatrix(trans_down);
  Bc.applyMatrix4(trans_down);
  Cc.applyMatrix4(trans_down);
  O.applyMatrix4(trans_down);

  let Fx = Cc.clone();
  let Tx = Bc.clone();
  // I hate non-functional math....
  let Bnx = Bnl.clone();
  let Cnx = Cnl.clone();
//  let dir = Cc.clone().sub(Bc.clone());
//  dir.normalize();


  var lgeometry = new THREE.Geometry();
  lgeometry.vertices.push(
    new THREE.Vector3(Bc.x,Bc.y,Bc.z-0.1),
    new THREE.Vector3(Cc.x,Cc.y,Cc.z+0.3)
  );

  var lmaterial = new THREE.LineBasicMaterial({
	color: 0xff00
  });

  var line = new THREE.Line( lgeometry, lmaterial );

  var group = new THREE.Group();

  var obj = new THREE.Mesh( geo, new THREE.MeshNormalMaterial(
    { transparent: true,
      opacity: 0.5 }
  ));
  group.add(obj);
  group.name = "SUPERSTRUCTURE";

  var osphere = createSphere(1/20, new THREE.Vector3(0,0,0), "yellow");
  group.add(line);
  group.add(osphere);

  return [group,Bnx,Cnx,SCALE_TO_UNIT_EDGE];
}

$( document ).ready(function() {
  runUnitTests();
//  $("#construct_via_norms").prop('checked', true);

  $( "#angle_omega_slider" ).slider( "value",INIT_RHO );
  $( "#angle_omega" ).val( INIT_RHO );
  ANGLE_OMEGA_d = INIT_RHO;
  $( "#angle_rho_slider" ).slider( "value",INIT_OMEGA );
  $( "#angle_rho" ).val( INIT_OMEGA );
  ANGLE_RHO_d = INIT_OMEGA;

  setup_input_molecule("#tau_slider","#tau_d","#tau_txt",TAU_d,(v => TAU_d = v));

  setup_input_molecule("#normal_b_x_slider","#normal_b_x",
                       "#b_x",NORMAL_B_X,(v => NORMAL_B_X = v));
  setup_input_molecule("#normal_b_y_slider","#normal_b_y",
                       "#b_y",NORMAL_B_Y,(v => NORMAL_B_Y = v));
  setup_input_molecule("#normal_b_z_slider","#normal_b_z",
                       "#b_z",NORMAL_B_Z,(v => NORMAL_B_Z = v));
  setup_input_molecule("#normal_c_x_slider","#normal_c_x",
                       "#c_x",NORMAL_C_X,(v => NORMAL_C_X = v));
  setup_input_molecule("#normal_c_y_slider","#normal_c_y",
                       "#c_y",NORMAL_C_Y,(v => NORMAL_C_Y = v));
  setup_input_molecule("#normal_c_z_slider","#normal_c_z",
                       "#c_z",NORMAL_C_Z,(v => NORMAL_C_Z = v));

  function faceValueChanged( event, ui ) {
    var face = $( "#face-spinner" ).spinner( "value" );
    var solid = getPlatonicSolidInput();
    if (solid) {
      face = Math.max(1,face);
      face = Math.min(numFaces(solid)-1,face);
      $( "#face-spinner" ).spinner( "value",face );
      var taus = updateLegalTauValues(solid,face);
      tau_v = getSelectedTaus(taus);
      var modex = $("input:radio[name=radio-2]:checked")
      modex.removeAttr('checked');
      modex.button( "refresh" )
    }
    onComputeDelix();
  }

  function linkValueChanged( event, ui ) {
    var links = $( "#link-spinner" ).spinner( "value" );
    links = Math.min(links,100);
    links = Math.max(links,1);
    NUM_PRISMS = links;
    NUM_SEGMENTS = (2 * NUM_PRISMS) + 1;
    onComputeDelix();
  }

   $( function() {
     var face_spinner = $( "#face-spinner" ).spinner();
     $( "#face-spinner" ).on( "spinchange",faceValueChanged);
     $( "#face-spinner" ).on( "spin",faceValueChanged);
     $( "#face-spinner" ).spinner( "value", 1 );
     $("#face-spinner-span").hide();

     var link_spinner = $( "#link-spinner" ).spinner();
     $( "#link-spinner" ).on( "spinchange",linkValueChanged);
     $( "#link-spinner" ).on( "spin",linkValueChanged);
     $( "#link-spinner" ).spinner( "value", NUM_PRISMS );

      $("#platonic_parameters").hide();
} );

  $(function () { main(); });

  onComputeDelix();

  populatePlatonicHelixTable();

  $('.collapse').trigger('click');


});
