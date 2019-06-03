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

// for testing, we need to know when somethigns is "closeto a target"
// to deal with roundoff error

function near(x, y, e = 1e-5) {
  return Math.abs(x - y) <= e;
}

function vnear(a, b, e = 1e-5) {
  return near(a.x,b.x,e) && near(a.y,b.y,e) && near(a.z,b.z,e);
}


function QFromRhoOmega(L,rho,omega) {
  return L/(Math.sqrt(1 + Math.tan(rho) ** 2 + Math.tan(omega) ** 2));
}
function ChordFromLD(L,di) {
  return Math.sqrt(L ** 2 - di**2);
}
function ThetaFromRC(r,c) {
  return Math.abs(2 * Math.asin(c/(2*r)));
}

// The "Sup" here is a "superstructure"
function AbstractPrism(L,Nb,Nc,Sup = null) {
  this.L = L;
  this.Nb = Nb.clone().normalize();
  this.Nc = Nc.clone().normalize();
  this.superstructure_prototype = Sup;
}

function ChordFromLDaxis(L,Da) {
  return Math.sqrt((L**2) - (Da**2));
}
// Note this has a "special" case...
function RotationFromRadiusChord(R,C) {
  let v = C / (2 * R);
  if (near(v,1,1e-4)) {
    return Math.PI;
  } else {
    let theta = 2 * Math.asin(v);
    if (isNaN(theta)) {
      debugger;
    }
    return 2 * Math.asin(v);
  }
}

function testAxis(B,Ba,H,da) {
    // Now test that H*da + Ba = Ca;
  let Ca = new THREE.Vector3(-Ba.x,Ba.y,-Ba.z);
  let X = Ba.clone();
  let Hn = H.clone().normalize();
  let Hc = Hn.clone();
  // I don't understand this!!!
  // The only way this makes sense is if I have da reversed EVERYWHERE.
  Hc.multiplyScalar(da);
  X.add(Hc);
  if (!vnear(X,Ca)) {
//    debugger;
  }
  console.assert(vnear(X,Ca));

  var Ba_m_B = Ba.clone();
  // Ba_minus_B
  Ba_m_B.sub(B);
  console.assert(near(0,H.dot(Ba_m_B)),"dot product");
}



// D is the third physical point on Z-axis centered segmented helix, L is the length.
// The return values are:
// [r,theta,da,c,phi,H,Ba];
// r -- the radius
// theta -- the signed rotation of one point to the next around the axis.
// da -- absolute value along the axis when going from one axis perpencicular to the next
// (that is, Ba to Ca)
// c -- the chord length along the axis
// phi -- the signed angle of the Helix axis against a segment
// H -- the vector of the axis of the helix
// Ba -- the point on the axis closest to point B

function KahnAxis(L,D) {
  let x = D.x;
  let y = D.y;
  let z = D.z;
  let A = new THREE.Vector3(-x,y,-z);
  let B = new THREE.Vector3(0,0,-L/2);
  let C = new THREE.Vector3(0,0,L/2);

  let Bb = A.clone();
  Bb.add(C);
  Bb.multiplyScalar(1/2);
  Bb.sub(B);


  // If Bb and Cb have zero length, the cross product
  // is undefined, and the axis of the helix may be
  // taken to be anywhere (is not uniquely defined.)
  // This is the case of a "straight line".
  if (near(Bb.length(),0)) {
    // We will take H to be the vector point from B to C,
    // use a zero radius.
    let H = new THREE.Vector3(0,0,L);
    let da = L;
    let r = 0;
    let c = 0;
    // I need to define this more clearly....
    //    let phi = Math.PI;
    let phi = 0;
    // We don't have enough information to define theta!!!
    let theta = null;
    return [r,null,da,c,phi,H,null];
  } else {
    let CmB = C.clone();
    CmB.sub(B);
    CmB = new THREE.Vector3(0,0,L);
    if (near(y,0,1e-4)) { // flat case

      // Note: I think at a minimum the compuation of Ba below is wrong.
      // This case needs to be tested more extensively.
      // This also should be integrated with the code below if possible.

      // If the flat case, the axis is parallel to the vector
      // AC or BD.
      let CmA = C.clone();
      CmA.sub(A);
      let H = CmA.clone();
      H.normalize();
      da = CmB.dot(H);
      let r = Bb.length() / 2;
      let c = 2*r;
      // How do I know phi should be positive?
      let phix = Math.acos(da/L);
      // This seems a better alternative....
      let phi = Math.atan2(H.z,H.x) - Math.PI/2;
      let theta = Math.PI;
      if (A.x < 0) {
        theta = -theta;
      }
      let Bax = Math.sqrt(1 - da**2/L**2) * Math.abs(da) /2;
      if (A.x < 0) { // this is becasue of da
        Bax = -Bax;
      }
      let Ba = new THREE.Vector3(Bax,0, -(da**2)/(2*L));
      testAxis(B,Ba,H,da);
      return [r,theta,da,c,phi,H,Ba];
    } else {
      let Cb = new THREE.Vector3(-Bb.x,Bb.y,-Bb.z);
//      let H = new THREE.Vector3(0,0,0);
      // H.crossVectors(Bb,Cb);
      // Hd is H computed direction....
      // This is optimization of that...
      let H = new THREE.Vector3(-2 * Bb.y * Bb.z,0, 2 * Bb.y * Bb.x);
//      console.assert(Hd.x == H.x && Hd.y == H.y && Hd.z == H.z);

//      H = Hd.clone();
      H.normalize();

      // da is the length of the projection of BC onto H
      //      let dax = CmB.dot(H);
      let da = L * Math.abs(Bb.x) / (Math.sqrt(Bb.x**2 + Bb.z**2));
      // If the sense is CW, the travel is negative...
      if (A.x > 0) {
        da = -da;
      }
      // phi is now the angle of H with the z axis
      // phi loses information, is sometimes off by PI
      // let phix = Math.acos(da/L);

      // phi is actually correct as measured against the z axis....
      // TODO: This value of phi requires us to change the rendering
      let phi = -(Math.atan2(H.z,H.x) - Math.PI/2);
      if (near(phi,-Math.PI)) {
        console.log("PHI is 180! Correcting!");
        phi = 0;
      }

      if (phi < 0) phi = -phi;

      // This sin could be computed with a perp dot....
      // Note: Sin(acos(x)) = sqrt(1 - x^2) https://socratic.org/questions/how-do-you-simplify-sin-arccos-x-1
      //          let Baxold = Math.sin(phi) * da /2;

      // How do we determine the sign here?
      // If phi is > 180 (and less than 360), then Bax is negative.
      let Bax = Math.sqrt(1 - da**2/L**2) * Math.abs(da) /2;

      if (A.x < 0) { // this is becasue of da
        Bax = -Bax;
      }
      if (!(Math.sign(Bax) == Math.sign(A.x))) {
        debugger;
      }
      console.assert(Math.sign(Bax) == Math.sign(A.x));
      // NOTE: Try to compute the point Ba
      // from the H vector and da!
      // This be undefinied if phi = PI/2.
      // if phi = PI/2, then B on the axis is
      // the intersection of the Bb and Cb,
      // which will also have x and z = 0.
      var Ba;
      if (near(Math.abs(phi),Math.PI/2,1e-4)) {
        // let psi = Math.atan2(Bb.y,Bb.z);
        Ba = new THREE.Vector3(0,
                               // Math.tan(psi)* L/2,
                               (Bb.y * L) / (Bb.z * 2),
                               0);
        console.assert(near(Math.tan(Math.atan2(Bb.y,Bb.z)),Bb.y/Bb.z));
      } else {
        Ba = new THREE.Vector3(Bax,
                               Bb.y * ( Bax / Bb.x),
                               // -Math.cos(phi) * da /2
                               // Note: substitution here
                               // of da will allows us to express better.
                               // is this really always negative?
                               // if the face angles are steep,
                               // then Ba could be on the ohter side of Ca.
                               // TODO: Investigate with steep angles and
                               // create a test case later.
                               -(da**2)/(2*L));
      }

      // here we assert that Ba is on the axis...
      // To make sure we have the sign right, in our model
      // the y value of Ba must be below the y = 0 plane...
      console.assert(Ba.y <= 0);
      var Ba_m_B = Ba.clone();
      // Ba_minus_B
      Ba_m_B.sub(B);
      let r = Ba_m_B.length();
      let c = ChordFromLDaxis(L,da);
      let theta = RotationFromRadiusChord(r,c);
      testAxis(B,Ba,H,da);
      return [r,theta,da,c,phi,H,Ba];
    }
  }
}

// This code taken from:
// https://stackoverflow.com/questions/22360936/will-three-js-object3d-clone-create-a-deep-copy-of-the-geometry

    /** Gives the aptitude for an object3D to clone recursively with its material cloned (normal clone does not clone material)*/

    THREE.Object3D.prototype.GdeepCloneMaterials = function() {
        var object = this.clone( new THREE.Object3D(), false );

        for ( var i = 0; i < this.children.length; i++ ) {

            var child = this.children[ i ];
            if ( child.GdeepCloneMaterials ) {
                object.add( child.GdeepCloneMaterials() );
            } else {
                object.add( child.clone() );
            }

        }
        return object;
    };

    THREE.Mesh.prototype.GdeepCloneMaterials = function( object, recursive ) {
        if ( object === undefined ) {
            object = new THREE.Mesh( this.geometry, this.material.clone() );
        }

        THREE.Object3D.prototype.GdeepCloneMaterials.call( this, object, recursive );

        return object;
    };
// Add the novel prism to the old prism by placing the B face
// aginst the C face of the old with a twist of tau and return.
// note that old is a prism instance, not a mesh instance,
// and it has a link to the abstract prism inside it.
// NOTE: with faces reversed, this still produces
// positive z.
// NOTE: A posssible goal here is to return the
// transforms needed for this transformation, to apply to other
// Geometry objects in the same frame.
// Note: I believe the up vector can be computed from
// tb and tc.
function adjoinPrism(old,tau,joinToC,debug = false) {
  // First, we copy the old prism in exactly the same position

  // I was using clone here, but it seems to be unreliable!!!
  var nu = Object.assign({},old);
  nu.p = Object.assign({}, old.p);
  nu.b = new THREE.Vector3(old.b.x,old.b.y,old.b.z);
  nu.c = new THREE.Vector3(old.c.x,old.c.y,old.c.z);

  nu.tb = new THREE.Vector3(old.tb.x,old.tb.y,old.tb.z);
  nu.lb = new THREE.Vector3(old.lb.x,old.lb.y,old.lb.z);
  nu.rb = new THREE.Vector3(old.rb.x,old.rb.y,old.rb.z);

  nu.tc = new THREE.Vector3(old.tc.x,old.tc.y,old.tc.z);
  nu.lc = new THREE.Vector3(old.lc.x,old.lc.y,old.lc.z);
  nu.rc = new THREE.Vector3(old.rc.x,old.rc.y,old.rc.z);

  // Then we translate along the axis of the old prism
  var av;
  if (joinToC)
  {
    var temp = old.c.clone();
    temp.sub(nu.b);
    av = temp;
  }
  else {
    var temp = old.b.clone();
    temp.sub(nu.c);
    av = temp;
  }
  var trans = new THREE.Matrix4();
  trans.makeTranslation(av.x,av.y,av.z);

  applyMatrix4ToPrism(nu,trans);

  var p_trans = new THREE.Matrix4();
  (joinToC) ?
    p_trans.makeTranslation(-nu.b.x,-nu.b.y,-nu.b.z):
    p_trans.makeTranslation(-nu.c.x,-nu.c.y,-nu.c.z);
  var p_trans_r = new THREE.Matrix4();
  (joinToC) ?
    p_trans_r.makeTranslation(nu.b.x,nu.b.y,nu.b.z) :
    p_trans_r.makeTranslation(nu.c.x,nu.c.y,nu.c.z);

  applyMatrix4ToPrism(nu,p_trans);

  var av =
      (joinToC) ?
      new THREE.Vector3(0,0,0).subVectors(nu.c,nu.b):
      new THREE.Vector3(0,0,0).subVectors(nu.b,nu.c);

  // At this, point b had better be at the origin...
  if (joinToC) {
    console.assert(near(nu.b.length(),0,1e-4));
  } else {
    console.assert(near(nu.c.length(),0,1e-4));
  }

  // At this point I may be making an incorrect assumptiong that
  // the up vectors are the same. It is not clear to me what
  // the quaternion is doing with the up vectors. The basic goal
  // here is to rotate the new prism so that it is face-to-face
  // with the old one.

  // Then we rotate about the joint (we may actually do this first)
  var Q_to_Nb = new THREE.Quaternion();
  var Q_to_Nc = new THREE.Quaternion();
  if (joinToC) {
    Q_to_Nb.setFromUnitVectors(nu.p.Nb,
                               new THREE.Vector3(0,0,-1));
    Q_to_Nc.setFromUnitVectors(new THREE.Vector3(0,0,1),
                               nu.p.Nc);
  } else {
    Q_to_Nc.setFromUnitVectors(nu.p.Nc,
                               new THREE.Vector3(0,0,1));
    Q_to_Nb.setFromUnitVectors(new THREE.Vector3(0,0,-1),
                               nu.p.Nb);
  }

  var Rqb = new THREE.Matrix4().makeRotationFromQuaternion(Q_to_Nb);
  var Rqc = new THREE.Matrix4().makeRotationFromQuaternion(Q_to_Nc);
  var Unknown = new THREE.Matrix4().identity();
  if (joinToC) {
    applyQuaternionToPrism(nu,Q_to_Nb);
    applyQuaternionToPrism(nu,Q_to_Nc);
    Unknown.multiplyMatrices(Rqc,Rqb);
  } else {
    applyQuaternionToPrism(nu,Q_to_Nc);
    applyQuaternionToPrism(nu,Q_to_Nb);
    Unknown.multiplyMatrices(Rqb,Rqc);
  }

  // Apparently, Unknown fundamentally makes sure we are z-aligned
  // By taking into about both normals. This exposes
  // the possiblity that we could construct this matrix
  // analytically (though I don't know how) which would
  // let us analytically apply the latter functions.
  // In this way we could maybe analytically produce answers
  // for all platonic solids.  What would that be worth? a week,
  // for sure!

  var rt = new THREE.Matrix4();
  rt.makeRotationAxis(joinToC ? nu.p.Nc : nu.p.Nb,tau);

  applyMatrix4ToPrism(nu,rt);

  applyMatrix4ToPrism(nu,p_trans_r);

  var finalMatrix = new THREE.Matrix4().identity();
  var rotations = new THREE.Matrix4().identity();
  // Question: Does this work if we take out the translational
  // elements? YES IT DOES.
  finalMatrix.premultiply(trans);
  finalMatrix.premultiply(p_trans);
  rotations.premultiply(Unknown);
  rotations.premultiply(rt);
  finalMatrix.premultiply(rotations);
  finalMatrix.premultiply(p_trans_r);

  // now we need to check that this matrix applied
  // to the original prism puts it int the same place...
  // if we apply this matrix to the point B, we expect
  // and we are joining to B, we expect it in end up
  // in the same place as point nu.b.
  if (true) {
    // Now, let me try to extract the data and check it...

    // Note: it is not clear that this is computing
    // the axis in the same direction that I do. However,
    // it seems to be correct modulo that.
    if (debug){
      var btest = old.b.clone();
      btest.applyMatrix4(finalMatrix);
      console.assert(vnear(btest,nu.b));
      console.log("btest,nu.b",btest,nu.b);
      if (!vnear(btest,nu.b)) {
        console.log("SHOULD BE NEAR",btest,nu.b);
        debugger;
      }
      console.log(finalMatrix);
      const L = nu.b.distanceTo(nu.c);
      [radius,thetaP,da,chord,phi,u] = computeThetaAxisFromMatrix4(L,finalMatrix,nu.b);

      console.log("L,RADIUS",L,radius);
      console.log("THETA, phi",
                thetaP * 180 / Math.PI,
                phi * 180 / Math.PI);
      console.log("BBB da, chord",
                  da,chord);
    }
  }

  var diff = new THREE.Vector3();
  if (joinToC) {
    diff.subVectors(old.c,nu.b);
  } else {
    diff.subVectors(old.b,nu.c);
  }
  console.assert(near(diff.length(),0,1e-4));
  if (!near(diff.length(),0,1e-4)) {
    console.log("DIFF",diff);
  }
  // So that the normal of the C face is the opposite of the B face.
  return [nu,finalMatrix];
}

// TODO--I think this function needs to apply
// the transformation to the normals of the abstract prism
// as well, if it is a rotation about z.
function applyMatrix4ToPrism(nu,trans) {
  nu.tb.applyMatrix4(trans);
  nu.lb.applyMatrix4(trans);
  nu.rb.applyMatrix4(trans);

  nu.tc.applyMatrix4(trans);
  nu.lc.applyMatrix4(trans);
  nu.rc.applyMatrix4(trans);

  nu.b.applyMatrix4(trans);
  nu.c.applyMatrix4(trans);
  // if (nu.sup) {
  //    console.log("TRANS before:",nu.sup);
  //    nu.sup.applyMatrix(trans);
  //    console.log("TRANS after:",nu.sup);
  // }
}

function applyQuaternionToPrism(nu,q) {
  nu.tb.applyQuaternion(q);
  nu.lb.applyQuaternion(q);
  nu.rb.applyQuaternion(q);

  nu.tc.applyQuaternion(q);
  nu.lc.applyQuaternion(q);
  nu.rc.applyQuaternion(q);

  nu.b.applyQuaternion(q);
  nu.c.applyQuaternion(q);
  // if (nu.sup) {
  //    console.log("QUAT before:",nu.sup);
  //    nu.sup.applyQuaternion(q);
  //    console.log("QUAT after:",nu.sup);
  // }
}

// absp is an abstract prism.
// posZ => c is in positive z direction.
function CreatePrism(absp,PRISM_FACE_RATIO_LENGTH) {
  // This object is a triangular prism. We will use
  // parallelpiped and spheres to create it.

  // First we have to compute the three
  // points at the ends be rotating by the normals.
  // We begin with Equilateral triangles at the ends.
  // The axis is Z-aligned, and will have one edge
  // pointing up.

  const W = TRIANGLE_WIDTH = 1;
  const H = TRIANGLE_HEIGHT = Math.sqrt(3)/2;
  const BASE = -(1/3);

  var TB = new THREE.Vector3(0, H + BASE, 0);
  var LB = new THREE.Vector3(-W/2, BASE, 0);
  var RB = new THREE.Vector3(W/2, BASE, 0);

  let PRISM_FACE_LENGTH = absp.L * PRISM_FACE_RATIO_LENGTH;
  TB.multiplyScalar(PRISM_FACE_LENGTH);
  LB.multiplyScalar(PRISM_FACE_LENGTH);
  RB.multiplyScalar(PRISM_FACE_LENGTH);

  var TC = TB.clone();
  var LC = LB.clone();
  var RC = RB.clone();

  // now we have to rotate these points.
  var Q_to_Nb = new THREE.Quaternion();
  Q_to_Nb.setFromUnitVectors(new THREE.Vector3(0,0,-1),
                             absp.Nb);
  var Q_to_Nc = new THREE.Quaternion();
  Q_to_Nc.setFromUnitVectors(new THREE.Vector3(0,0,1),
                             absp.Nc);
  TB.applyQuaternion(Q_to_Nb);
  LB.applyQuaternion(Q_to_Nb);
  RB.applyQuaternion(Q_to_Nb);

  TC.applyQuaternion(Q_to_Nc);
  LC.applyQuaternion(Q_to_Nc);
  RC.applyQuaternion(Q_to_Nc);

  // Now translate them...
  var transB = new THREE.Matrix4();
  var transC = new THREE.Matrix4();

  let b = new THREE.Vector3(0,0,0);
  let c = new THREE.Vector3(0,0,0);

  transB.makeTranslation(0,0,-absp.L/2);
  transC.makeTranslation(0,0,absp.L/2);
  b.applyMatrix4(transB);
  c.applyMatrix4(transC);

  TB.applyMatrix4(transB);
  LB.applyMatrix4(transB);
  RB.applyMatrix4(transB);

  TC.applyMatrix4(transC);
  LC.applyMatrix4(transC);
  RC.applyMatrix4(transC);

  // Now the superstruction prototype should be cloned...
  var sup = null;
  // If the suprestructure exists, we presumably
  // computed the face normal from it, so we don't
  // presumbaly TB,LB,RB and TC,LC,RC lie on those
  // faces already
  if (absp.superstructre_prototype) {
    sup = absp.superstructre_prototype.clone();
  }

  // for convenience, should I make the instance vectors
  // part of this object? Or should I add the matrices,
  // effective make it a "screw" transformation
  return {p: absp,
          b: b,
          c: c,
          nb: absp.Nb,
          nc: absp.Nc,
          tb: TB, lb : LB, rb : RB,
          tc : TC, lc : LC, rc : RC,
          sup: sup};
}

function condition_angle(angle) {
  if (angle < (-2 * Math.PI)) {
    return condition_angle(angle + (2 * Math.PI));
  } else if (angle > (2 * Math.PI)) {
    return condition_angle(angle - (2 * Math.PI));
  } else {
    return angle;
  }
}

function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}


// These functions are taken from the Mathematica code.
function V0fromLNB(L,Nb,Nc,delta) {
  var rt = new THREE.Matrix4();
  rt.makeRotationAxis(new THREE.Vector3(1,0,0),-delta);
  var v0 = Nb.clone();
  v0.applyMatrix4(rt);
  v0.normalize();
  return v0.multiplyScalar(L);
}

function ADirFromParam(L,v0,tau,Nb) {
  var k = Nb.clone();
  k.normalize();
  var rt = new THREE.Matrix4();
  rt.makeRotationAxis(k,tau);
  vc = v0.clone();
  vc.applyMatrix4(rt);
  vc.normalize();
  vc.multiplyScalar(L);
  return vc;
}


// Compute the A point of an ideally centered prism
// with normal faces NBu, Ncu.
// Maybe this should normalize the input for robustness?
function AfromLtauNbNc(L,tau,NBu,NCu,debug = false) {
  var abs0 = new AbstractPrism(L,NBu,NCu);
  if (!abs0) {
    debugger;
  }
  var p_i = CreatePrism(abs0,PRISM_FACE_RATIO_LENGTH);

  var rt = new THREE.Matrix4();
  var p_b = adjoinPrism(p_i,tau,false,debug)[0];
  var p_c = adjoinPrism(p_i,tau,true,debug)[0];

  let fb = new THREE.Vector2(p_b.b.x,p_b.b.y);
  let fc = new THREE.Vector2(p_c.c.x,p_c.c.y);

  function compute_angle_midpoint(fb,fc) {
    let fbc = new THREE.Vector2().addVectors(fb,fc);
    // if fbc is zero, that means the vectors perfectly
    // oppose each other. We should produce a psi to
    // rotate so that b is in the X - Z plane

    if (near(fbc.length(),0,1e-3)) { // The angles are in opposition
      let fba = fb.angle();
      let fca = fc.angle();
      let psi = -( (fba > fca) ? fba - Math.PI : fca - Math.PI );
      if (tau == -Math.PI) {
        return psi+Math.PI;
      }
      return psi;
    }
    // The ".angle()" function in THREE measures against the X-axis
    // I assume this is "reverse on the clock" operation!
    let psi = -(fbc.angle() - ( 3 * Math.PI / 2));
    return psi;
  }


  // what is the meaning of psi? It is the angle midpoint
  // between the projects of fb and fc
  let psi = compute_angle_midpoint(fb,fc);

  // That is why we rotate around the z axis to get down
  // to find the A point.
  rt.makeRotationAxis(new THREE.Vector3(0,0,1),psi);

  applyMatrix4ToPrism(p_i,rt);
  p_i.p.Nb.applyMatrix4(rt);
  p_i.p.Nc.applyMatrix4(rt);

  var p_b = adjoinPrism(p_i,tau,false,debug);
  var p_c = adjoinPrism(p_i,tau,true,debug);

  // We are assuming the rotation goes from C to D, even though we are returning A
  return [p_b[0].b,psi,p_c[1]];
}


// The fundamental problem here is that theta
// is NOT being computed as the best theta that can be
// returned. It is off by a few degrees.  This is a
// problem in ComputingBalance!!!
function testAfromLfailureCase1() {
  let Nb = new THREE.Vector3(-0.36496094853172817,0,-0.9310228278870616);
  let Nc = new THREE.Vector3(-0.39460613482779394,-0.40696332067835117,0.8238123900371481);
  let L = 1;
  let tau = 0;
  let res = AfromLtauNbNc(L,tau,Nb,Nc);
}

function testAfromLtauMultiple() {
  // my goal here is to test with a variety of
  // vaules, and in each case to make sure that
  // we have found a rotation angle that "balances"
  // This means that A and D should have the same
  // x and y value.  We have two ways of
  // testing this, which provides us an independent calculation.
  const NUM_TEST = 13;

  for(var i = 0; i < NUM_TEST; i++) {
    var Nbx = -1 + i*(2/NUM_TEST);
    for(var j = 0; j < NUM_TEST; j++) {
      var Nby = -1 + j*(2/NUM_TEST);
      for(var k = 0; k < NUM_TEST; k++) {
        var Ncx = -1 + k*(2/NUM_TEST);
        for(var l = 0; l < NUM_TEST; l++) {
          var Ncy = -1 + l*(2/NUM_TEST);
          let Nb = new THREE.Vector3(Nbx,Nby,-1);
          let Nc = new THREE.Vector3(Ncx,Ncy,1);
          Nb.normalize();
          Nc.normalize();
          let L = 1;
          let tau = 0;
          let res = AfromLtauNbNc(L,tau,Nb,Nc);
          let A = res[0];
          //                    let vs = res[1];
          let theta = res[1];
          // This is code based on building the physical prism.
          // in a sense it is more accurate. The fact
          // that these values don't match is a serious problem.

          var abs = new AbstractPrism(L,Nb,Nc);
          var p_i = CreatePrism(abs,PRISM_FACE_RATIO_LENGTH);

          var rt = new THREE.Matrix4();

          rt.makeRotationAxis(new THREE.Vector3(0,0,1),theta);

          applyMatrix4ToPrism(p_i,rt);
          p_i.p.Nb.applyMatrix4(rt);
          p_i.p.Nc.applyMatrix4(rt);


          var p_b = adjoinPrism(p_i,tau,false,false)[0];
          var p_c = adjoinPrism(p_i,tau,true,false)[0];
          return;
        }
      }
    }
  }
}


function testAfromLtauNbNn(tau) {
  let L = 1;
  let NB0 = new THREE.Vector3(0,-Math.sin(Math.PI/8),-Math.cos(Math.PI/8));
  let NC0 = new THREE.Vector3(0,-Math.sin(Math.PI/8),Math.cos(Math.PI/8));
  let tetAngle = Math.atan(Math.sqrt(2)/2);
  let NB1 = new THREE.Vector3(0,-Math.sin(tetAngle),-Math.cos(tetAngle));
  let NC1 = new THREE.Vector3(0,-Math.sin(tetAngle),Math.cos(tetAngle));
  let AA = AfromLtauNbNc(L,tau,NB1,NC1);
  return AA;
}

function testAfromLtauNbNnContinuityOfTau() {
  var a = testAfromLtauNbNn(134 * Math.PI / 180);
  var b = testAfromLtauNbNn(135 * Math.PI / 180);
  console.assert(near(a[0].x,b[0].x,0.2));
  console.assert(near(a[0].y,b[0].y,0.2));
}

function testAfromLtauNbNnKeepsYPure() {
  let L = 1;
  let NB0 = new THREE.Vector3(0, -0.25254006068835666,-0.967586439418991);
  let NC0 = new THREE.Vector3(0,0.31337745420390184,0.9496286491027328);
  let tau = 0;
  let A = AfromLtauNbNc(L,tau,NB0,NC0)[0];
  console.assert(near(A.x,0,1e-6));
}




// This tests the well known situation
// of the Boerdijk-Coxeter tetrahelix.
// This models a tetrahedron of side length 1.
function testRegularTetsKahnAxis()
{
  let L0 = 1/3;
  let tetAngle = Math.atan(Math.sqrt(2)/2);
  let NB1 = new THREE.Vector3(0,-Math.sin(tetAngle),-Math.cos(tetAngle));
  let NC1 = new THREE.Vector3(0,-Math.sin(tetAngle),Math.cos(tetAngle));
  let tau = 2 * Math.PI /3;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let B = new THREE.Vector3(0,0,-L0/2);
  // let D = new THREE.Vector3(-A.x,A.y,-A.z);
  // let res = KahnAxis(L0,D);
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  // I am not sure the negative sign here is justified,
  // I need to check on this.
  let theta_exp = Math.acos(-2/3);
  // We know (thanks to Coxeter) the expected theta!
  console.assert(near(theta,theta_exp,0.00001));

  console.assert(H.x > 0);
  // We assert in a situation like this that da should be positive...
  console.assert(da > 0);

  testAxis(B,Ba,H,da);
}

// This tests a common special case: when the helix
// is degenerate, forming a polygon (or a torus, if you
// consider the object to have physical dimensions.)
function testKahnAxisYTorus()
{
  let N = 17;
  let angle = Math.PI/17;
  for(let a = angle; a < Math.PI*2; a += angle) {
    testKahnAxisYTorusAux(angle);
  }
}
function testKahnAxisYTorusAux(angle)
{
  let L0 = 2;
  let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
  let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
  let tau = 0;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let B = new THREE.Vector3(0,0,-L0/2);
  //  let D = new THREE.Vector3(-A.x,A.y,-A.z);
  //  let res = KahnAxis(L0,D);
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];

    // Now we want to check that H and Da are correct...
  // H should have only only an x component...
  console.assert(near(0,H.y) && near(0,H.z));
  console.assert(near(1,H.x));
  console.assert(near(da,0));

  testAxis(B,Ba,H,da);
}
// That that tau = 180 is not degenerate..
function testKahnAxisTau180()
{
  let L0 = 2;
  let angle = Math.PI/7;
//  let angle = 0;
  let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
  // We add in a little here so we are not a torus
  let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
  let tau = Math.PI;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
//  let B = new THREE.Vector3(0,0,-L0/2);
//  let D = new THREE.Vector3(-A.x,A.y,-A.z);
//  let res = KahnAxis(L0,D);

  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  // What did I mean here?
  console.assert(!near(da,0,0.1));
  console.assert(!isNaN(theta));
}

function testKahnAxisFull()
{
  const N = 10;
  const J = 3;
  const K = 3;
  let angle = Math.PI/7;
  let NB1 = new THREE.Vector3(Math.sin(angle),0,-Math.cos(angle)).normalize();
  for(let k = -(K-1); k < K; k++) {
    for(let j = -(J-1); j < J; j++) {
      for(let i = -(N-1); i < N; i++) {
        // tau is limited to within +- 180.
        let tau = -((Math.PI) * i) / (N-1)/10
        let NC1 = new THREE.Vector3(j,k,1).normalize();
        testKahnAxisFullAux(tau, NB1, NC1);
      }
    }
  }
}
function testKahnAxisFullAux(tau,NB1,NC1)
{
  let L0 = 2;
//  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
//  let B = new THREE.Vector3(0,0,-L0/2);
//  let D = new THREE.Vector3(-A.x,A.y,-A.z);
//  let res = KahnAxis(L0,D);
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  let B = res[7];
  console.assert(!isNaN(theta));

  // Now test that H*da + Ba = Ca;
  let Ca = new THREE.Vector3(-Ba.x,Ba.y,-Ba.z);
  let X = Ba.clone();

  //  WHY DOES THIS WORK?????
  // This suggests Ba is correct but H points the wrong
  // way.
  //  let Hc = (tau < 0) ? H.clone().negate() : H.clone();
  testAxis(B,Ba,H,da);
}

// This is to test the Flat, Zig-Zag Case
function testKahnAxisFlat()
{
  let base = Math.PI/7;
  const N = 10;
  for(let i = -(N-1); i < N; i++) {
    let angle = -((base) * i) / ((N-1)/N);
    testKahnAxisFlatAux(angle);
  }
}
function compareMethods(L,tau,NB1,NC1) {
  let AR = AfromLtauNbNc(L,tau,NB1,NC1);
  let A = AR[0];
  let B = new THREE.Vector3(0,0,-L/2);
  let C = new THREE.Vector3(0,0,L/2);
  let D = new THREE.Vector3(-A.x,A.y,-A.z);
  let R = AR[2];
  let res = KahnAxis(L,D);
  [r,thetaP,da,chord,phi,u,Ba] = computeThetaAxisFromMatrix4(L,R,B);
  const BC = new THREE.Vector3().subVectors(C,B);
  let phi_p = Math.acos(da/BC.length());

  if (res[1] == null) {
    // In this case, KahnAxis is undefined and we cannot compare!!!
    return [...res,B];
  }
  console.assert(near(r,res[0]));
  console.assert(near(res[1],0) || near(thetaP,res[1]));
  if (!((near(res[1],0) || near(thetaP,res[1])))) {
    debugger;
  }
  console.assert(near(da,res[2]));
  console.assert(near(chord,res[3]));

  // TODO: This is producing a sign error....
  console.assert(near(phi,res[4]));
  if (!(near(phi,res[4]))) {
    console.log("phi,Res[4]", phi * 180 / Math.PI, res[4] * 180 / Math.PI);
    debugger;
  }

  // Although possibly I should make the axis the length of da,
  // but instead

  console.assert(vnear(u.clone().normalize(),res[5].clone().normalize()));
  if (!(vnear(u.clone().normalize(),res[5].clone().normalize()))) {
    debugger;
  }
  // TODO: Apparently computeThetaAxis is not computeing Ba correctly!!!
 // console.assert(vnear(Ba,res[6]));

  return [...res,B];
}
function testKahnAxisFlatAux(angle)
{
  let L0 = 2;

  let NB1 = new THREE.Vector3(-Math.sin(angle),0,-1);
  // We add in a little here so we are not a torus
  let NC1 = new THREE.Vector3(Math.sin(angle),0,1);
  let tau = Math.PI;
//  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
//  let B = new THREE.Vector3(0,0,-L0/2);
//  let D = new THREE.Vector3(-A.x,A.y,-A.z);
  //  let res = KahnAxis(L0,D);
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  let B = res[7];
  console.assert(!isNaN(theta));

  // Now test that H*da + Ba = Ca;
  if (Ba) {
    let Ca = new THREE.Vector3(-Ba.x,Ba.y,-Ba.z);

    // in the flat case H will be in the y = 0 plane
    console.assert(near(H.y,0));

    let ntheta = theta/Math.PI;
    console.assert(near(ntheta,1) || near(ntheta,-1));

    testAxis(B,Ba,H,da);
  } else {
      console.assert(!Ba);
  }
}

function testThreeIsRightHanded() {
  let A = new THREE.Vector3(1,0,0);
  let B = new THREE.Vector3(0,1,0);
  let Z = A.cross(B);
  // in a right-handed system, A cross B provides positve Z.
  console.assert(near(Z.x,0));
  console.assert(near(Z.y,0));
  console.assert(near(Z.z,1));
}


function TraceMatrix4(A) {
  let t = 0;
  let xAxis = new THREE.Vector3();
  let yAxis = new THREE.Vector3();
  let zAxis = new THREE.Vector3();
  A.extractBasis(xAxis,yAxis,zAxis);
  t = xAxis.x + yAxis.y + zAxis.z;
  // for(var i = 0; i < 4; i++) {
  //   console.log(i + i * 4);
  //   t += A.elements[i + i*4];
  // }
  return t;
}

function testTrace() {
  const theta = Math.PI/11;
  // now create a rotation matrix based on this angle...
  const axis = new THREE.Vector3(0,0,1);
  const R = new THREE.Matrix4().makeRotationAxis(axis,theta);
  const t0 = TraceMatrix4(R);
  console.log("TRACE0",t0);
}

function SubtractMatrices(A,B) {
  const Ae = A.elements();
  const Be = B.elements();
  const Ce = Ae.clone();
  for(var i = 0; i < Be.elements.length; i++) {
    Ce[i] -= Be[i];
  }
  const C = A.clone();
  C.fromArray(Ce);
  return C;
}


function computeThetaUDa(R) {
  if (near(R.determinant(),-1)) {
    console.log("R.det",R.determinant());
  }
  // now attempt to recover parameters from the rotation matrix....
  const val = (1/2) * (TraceMatrix4(R) - 1);
  const thetaP = (near(val,-1)) ? Math.PI :  Math.acos(val);

  // We need an arbitray point.. technically we should check
  // this is not on the axis, but that is unlikely....
  const B = new THREE.Vector3(-4,4,-4/2);
  const C = B.clone().applyMatrix4(R);
  let u = null;
  if (thetaP == Math.PI) {
    // In this case our normal method is numerically unstable
    // (it defines an axis vector of zero length.). However,
    // in this "zig-zag" case, u is easily computed as R^2(P) - P,
    // where P is an arbitrary point.
    // Sine we need B and C in other cases, we will compute D
    // as a further application of R to C.
    const D = C.clone().applyMatrix4(R);

    // This appears to somehow be conceptually wrong.
    u = new THREE.Vector3().subVectors(D,B);
    u.normalize();
  } else {

    // Now we will attempt to extract the screw axis from this matrix...
    // unfortunately THREE doesn't seem to implment matrix addition and substraction...

    let xAxis = new THREE.Vector3();
    let yAxis = new THREE.Vector3();
    let zAxis = new THREE.Vector3();
    R.extractBasis(xAxis,yAxis,zAxis);

    // Note: if there is no rotation, we have the special "straight line" case....

    // const Ro = R.clone();
    // const Rt = R.transpose();
    // const d = SubtractMatrices(Ro,Rt);
    const scale = 1/ (Math.sin(thetaP));
    // from https://en.wikipedia.org/wiki/Rotation_matrix
    // note element is column-major
    const b = yAxis.x;
    const c = zAxis.x;
    const d = xAxis.y;
    const f = zAxis.y;
    const g = xAxis.z;
    const h = yAxis.z;
    u = new THREE.Vector3(h-f,c-g,d-b);
    // This may not be the best way to detect it, but I believe
    // that if u is "vnear" zero we have an undefined "straight line" case.
    if (vnear(u, new THREE.Vector3(0,0,0))) {
      console.log("WARNING!!! NUMERICAL Instability!",thetaP);
      // This is occuring in the "flat" or "zig-zag" case!!!
      // I need to figure out a different way to compute this,
      // even though it seems to work, which is weird.
      //     debugger;
    }

    u.multiplyScalar(1/(2 * Math.sin(thetaP)));
  }

  //  const SIGN_ADJUST = -1;
  // this just accomplishes normalization, so it is unclear where we need it!
  //  u.multiplyScalar(SIGN_ADJUST);
  //  console.log('AXIS',u,scale);

  // Now I believe da is indepednet of L and any point.
  //  const B = new THREE.Vector3(-4,4,-4/2);
  //  const C = B.clone().applyMatrix4(R);
  const BC = new THREE.Vector3().subVectors(C,B);

  // Now BC.u sould give us da...
  const unorm = u.clone().normalize();
  console.log("UNORM",unorm);
  const da_scalar_project = BC.dot(unorm);

  return [thetaP,u,da_scalar_project];
}


// https://www.seas.upenn.edu/~meam520/notes02/EulerChasles4.pdf
// Note: This may solve my problem: https://en.wikipedia.org/wiki/Screw_axis#Computing_a_point_on_the_screw_axis
function computeRChord(L,thetaP,da) {
  let chord = ChordFromLDaxis(L,da);
  const r = chord / (2 * Math.sin(thetaP/2));
  return [r,chord];
}

function computePointIndependentParameters(L,R) {
  [thetaP,u,da] = computeThetaUDa(R);
  [r,chord] = computeRChord(L,thetaP,da);
  // Now this has in fact given us everyting we need until such time
  // as we want to compute Ba, which requires a point B to be specified....
  // Can all of this be true?
  return [r,thetaP,da,chord,u];
}

function computeThetaAxisFromMatrix4(L,R,B) {

  [r,thetaP,da,chord,u] = computePointIndependentParameters(L,R);

  // Now I believe da is indepednet of L and any point.
  // const B = new THREE.Vector3(-4,4,-L/2);
  const C = B.clone().applyMatrix4(R);
  const BC = new THREE.Vector3().subVectors(C,B);

  // NOTE: phi is NOT DEFINED until a point B is given!!!!
//  let phi = (near(da,BC.length())) ? 0 : Math.acos(da / BC.length());

  let phi = (near(da,L)) ? 0 : Math.acos(da / L);


  // I really have to understanad this---why is this failing!
  console.assert(near(L,BC.length()));
  if (!near(L,BC.length())) {
    debugger;
  }

  console.log("da,bc",da,BC.length());
  console.log("phi",phi * 180/Math.PI);
  if (isNaN(phi)) debugger;

  {
    // Note: This math was more or less my own invention,
    // by applying the matrix to an arbitrary point!
    // I'm not sure this is docuemnted anywhere else,
    // and I'm not sure how it relates to other methods.
    // To be more general, B should be an input
    // to this routine.
    // IMPLICATIONS: We can finally find a point
    // on the axis from an ARBITRARY point.
    const uu = u.clone().multiplyScalar(1).setLength(da);

    // I think this point has to be the correct radius
    // away from the axis---but how do we know that?
    //    const B = new THREE.Vector3(0,2,-L/2);
//    const B = new THREE.Vector3(-4,4,-L/2);
//    const C = B.clone().applyMatrix4(R);
//    const BC = new THREE.Vector3().subVectors(C,B);
    const phi_p = u.angleTo(BC);
    const L_p = B.distanceTo(C);
    let da_p = L_p * Math.cos(phi_p);
    let chord_p = ChordFromLDaxis(L_p,da_p);
    const r_p = chord_p / (2 * Math.sin(thetaP/2));
    // console.log("chord,chord_p",chord,chord_p,chord_p/chord);
    // console.log("L,L_p", L, L_p,L_p/L);
    // console.log("da,da_p", da, da_p,da_p/da);
    // console.log("r,r_p", r, r_p,r_p/r);
    // console.log("u",u);
    // console.log("phi, cos(phi)",phi_p * 180/Math.PI,Math.cos(phi_p));

    const Bp = new THREE.Vector3().subVectors(B,uu);
    // Now we want to create a vector of length r
    // on the x axis, and rotate by
   const Cp = new THREE.Vector3().subVectors(C,uu);

    const Bm = new THREE.Vector3().addVectors(B,Bp);
    Bm.multiplyScalar(1/2);

    const Cm = new THREE.Vector3().addVectors(C,Cp);
    Cm.multiplyScalar(1/2);

    const Mv = new THREE.Vector3().subVectors(B,C);
    const Mp = new THREE.Vector3().addVectors(B,C);
    Mp.multiplyScalar(1/2);

  //  const BBp = new THREE.Vector3().subVectors(B,Bp);
    const Q = Mv.clone()
    Q.cross(u);
    const ql = Math.sqrt(r_p**2 - (chord_p/2)**2);
//    console.log("ql:",ql);
    Q.setLength(-ql);
//    console.log("length:",Q.length());
  //  console.log("Bm,Q",Bm,Q);
    const Ba = Mp.clone().add(Q);

    const uuhalf = uu.clone();
    uuhalf.multiplyScalar(-1/2);
    const nBa = new THREE.Vector3().addVectors(Ba,uuhalf);
    P = nBa;
 //   console.log("B,C,Bp,Cp",B,C,Bp,Cp);
//    console.log("Bm,Cm,Q,Ba",Bm,Cm,Q,Ba);
    // Now in theory C is on the axis is
    return [r,thetaP,da,chord,phi,u,P,[B,C,Mp,Q]];
  }


//  return [r,thetaP,da,chord,phi,u,P];
}

function testAbilityToGetScrewAxisFromRotationMatrix() {
  const theta = Math.PI/11;
  // now create a rotation matrix based on this angle...
  const axis = new THREE.Vector3(0,0,1);
  //  const R = new THREE.Matrix4().makeRotationAxis(axis,theta);
  const q = new THREE.Quaternion().setFromAxisAngle(axis,theta);
  const R = new THREE.Matrix4().makeRotationFromQuaternion(q);
  // Now we will attempt to extract the screw axis from this matrix...
  // unfortunately THREE doesn't seem to implment matrix addition and substraction...
  const L = 1;
  const B = new THREE.Vector3(4,4,4);
  [r,thetaP,da,chord,phi,u] = computeThetaAxisFromMatrix4(L,R,B);

  // now check that the parameters match...
  console.assert(near(theta,thetaP));
  if (!vnear(axis,u)) {
    console.log("axis,U",axis,u);
  }
  console.assert(vnear(axis,u));
}

function runChaslesTests() {
  testTrace();
  testAbilityToGetScrewAxisFromRotationMatrix();
}

function runUnitTests() {
  runChaslesTests();
  testThreeIsRightHanded();
  testAfromLtauNbNnKeepsYPure();
  testAfromLtauNbNnContinuityOfTau();
  testRegularTetsKahnAxis();
  testKahnAxisYTorus();
  testKahnAxisFull();
  testKahnAxisFlat()
  testKahnAxisTau180();
  testAfromLfailureCase1();
  testAfromLtauMultiple();
}
