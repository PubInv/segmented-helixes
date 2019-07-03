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
  try {
    return near(a.x,b.x,e) && near(a.y,b.y,e) && near(a.z,b.z,e);
  } catch (e) {
    debugger;
    return new THREE.Vector3();
  }
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
  // this is a numerical nicety, for when Da is erroneously greater than Da
  if (Da > L) return 0;
  return Math.sqrt((L**2) - (Da**2));
}
// Note this has a "special" case...
function RotationFromRadiusChord(R,C) {
  let v = C / (2 * R);
  if (v > 1) {
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
  if (!vnear(X,Ca,1e-3)) {
    debugger;
  }
  console.assert(vnear(X,Ca,1e-3));

  var Ba_m_B = Ba.clone();
  // Ba_minus_B
//  console.log(Ba_m_B,B);
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

function PointAxis(L,A) {
  let x = A.x;
  let y = A.y;
  let z = A.z;
//  let D = new THREE.Vector3(-x,y,-z);
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
    let phi = 0;
    // We don't have enough information to define theta!!!
    let theta = null;
    return [r,null,da,c,phi,H,B];
  } else {
    let da = -L * Bb.x / (Math.sqrt(Bb.x**2 + Bb.z**2));
    // I wish this were tighter, but there is a bit of trouble...
    if (near(y,0,1e-3)) { // flat case
      let CmB = new THREE.Vector3(0,0,L);
      // If the flat case, the axis is parallel to the vector
      // AC or BD.
      let CmA = C.clone();
      CmA.sub(A);
      let H = CmA.clone();
      H.normalize();
//      da = CmB.dot(H);
      if (A.x > 0) {
        H.multiplyScalar(-1);
//        da = -da;
      }
      let r = Bb.length() / 2;
      let c = 2*r;
      // How do I know phi should be positive?
      let phix = Math.acos(da/L);
      // This seems a better alternative....
      let phi = (Math.atan2(-H.z,-H.x) - Math.PI/2);
      if (A.x < 0) {
        phi = -phi;
      }
      // We choose to measure phi so that it is positive
      // against the z axis...
      if (A.x > 0) {
        phi = Math.PI - phi;
      }
      let theta = Math.PI;
      // This is a similar choice when theta = PI made computeThetaUDa!!!
      let Bax = -Math.sqrt(1 - da**2/L**2) * da /2;

      let Ba = new THREE.Vector3(Bax,0, -(da**2)/(2*L));
      testAxis(B,Ba,H,da);
      return [r,theta,da,c,phi,H,Ba];
    } else {
//      let Cb = new THREE.Vector3(-Bb.x,Bb.y,-Bb.z);
//      let H = new THREE.Vector3(0,0,0);
      // H.crossVectors(Bb,Cb);
      // Hd is H computed direction....
      // This is optimization of that...
      let H = new THREE.Vector3(-2 * Bb.y * Bb.z,0, 2 * Bb.y * Bb.x);
      H.normalize();
      // da is the length of the projection of BC onto H
      //      let dax = CmB.dot(H);
      //      let da = L * Math.abs(Bb.x) / (Math.sqrt(Bb.x**2 + Bb.z**2));
      // This is negative due to our sign convention that
      // counterclockwise motion is represented by negative da
      //      let da = -L * Bb.x / (Math.sqrt(Bb.x**2 + Bb.z**2));
//      da = daq;

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
      let Bax = - Math.sqrt(1 - da**2/L**2) * da /2;

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

      let Bay = (near(Bb.x,0)) ?
      // When Bb.x is near zero, we have a torioidal
      // case. Then Bb.y / Bb.z is proportional
      // to the height from the axis.
          (Bb.y * L) / (Bb.z * 2) :
          Bb.y * ( Bax / Bb.x);

      const  Ba = new THREE.Vector3(
        Bax,
        Bay,
        -(da**2)/(2*L));
      // here we assert that Ba is on the axis...
      // To make sure we have the sign right, in our model
      // the y value of Ba must be below the y = 0 plane...
      console.assert(Ba.y <= 0);
      if (Ba.y > 0) {
        debugger;
      }
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

function checkNonScaling(M) {
  let t = 0;
  let xAxis = new THREE.Vector3();
  let yAxis = new THREE.Vector3();
  let zAxis = new THREE.Vector3();
  M.extractBasis(xAxis,yAxis,zAxis);
  let sx = new THREE.Vector3(xAxis.x,yAxis.x,zAxis.x);
  let sy = new THREE.Vector3(xAxis.y,yAxis.y,zAxis.y);
  let sz = new THREE.Vector3(xAxis.z,yAxis.z,zAxis.z);
  return (near(sx.length(),1) && near(sy.length(),1) && near(sz.length(),1));
}
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

  console.assert(!vnear(new THREE.Vector3(0,0,0),
                        (joinToC ? nu.p.Nc : nu.p.Nb)));

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
  console.assert(!vnear(new THREE.Vector3(0,0,0),(joinToC ? nu.p.Nc : nu.p.Nb)));

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
  console.assert(!vnear(new THREE.Vector3(0,0,0),(joinToC ? nu.p.Nc : nu.p.Nb)));
  var rt = new THREE.Matrix4();

  // I believe this is the "interception point" where we could compute
  // the tau need to make A.x = D.x in our standard model. Instead of
  // calling makeRotationAxis(v,a) with a known angle, we would set up
  // the linear algebra equations with "a" treated as the unknown,
  // and equation A.x and D.x.  That might be clearer in the matrix
  // operations below. We will have to "analyze" the transformation
  // matrix of course by taking it apart. However, we have nicely
  // functioning code up to this point, so it should be possible
  // to do a copypasta and produce a unit test around it; then
  // we ought to be able to work out the linear algebra with assurety.
  // Note that in principle, we could create a formula for tau as a function of psi.
  // This would allow us to enter psi = PI for example, in order for us
  // to learn the ways of wisdom.
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

  console.assert(checkNonScaling(finalMatrix));
  if (!checkNonScaling(finalMatrix)) {
    debugger;
  }
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

// Find the point ont the line l0l1 closest to point p
function closestPoint(P,A,B) {
  var D = B.clone().sub( A ).normalize();
  var d = P.clone().sub( A ).dot( D );
  var X = A.clone().add( D.clone().multiplyScalar( d ) );
  return X;
}



// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
// x0 is the point, x1 and x2 are points on the line
function pointToLineDist(x0,x1,x2) {
  const denom = new THREE.Vector3().subVectors(x2,x1);
  const x0mx1 = new THREE.Vector3().subVectors(x0,x1);
  const x0mx2 = new THREE.Vector3().subVectors(x0,x2);
  const numercross = new THREE.Vector3().crossVectors(x0mx1,x0mx2);
  return numercross.length()/denom.length();
}

function testClosestPoint()
{
  var A = new THREE.Vector3( 1, 1, 1 );
  var B = new THREE.Vector3( 2, 2, 2 );
  for(var i = 0; i < 10; i++) {
    var P = new THREE.Vector3( i, 0, 0 );
    var X = closestPoint(P,A,B);
    console.log(P,X);
    console.assert(near(pointToLineDist(P,A,B),X.distanceTo(P)));
  }
}


// compute the rejection of a onto b
function scalarProjectionOnto(a,b) {
  const nb = b.length();
  const dp = a.dot(b);
  return dp/nb;
}
function vectorProjectionOnto(a,b) {
  const bu = b.clone().normalize();
  const sp = scalarProjectionOnto(a,b);
  return bu.multiplyScalar(sp);
}
function vectorRejectionOnto(a,b) {
  var diff = new THREE.Vector3();
  const ap = vectorProjectionOnto(a,b);
  diff.subVectors(a,ap);
  console.assert(near(diff.dot(b),0));
  return diff;
}

// compute the signed distance between a point x
// and a plane with norm pn containing the origin.
// Note this is only correct for planes through the
// origin.
function distancePointToPlaneThroughOrigin(pn,x) {
  const numer = pn.dot(x);
  const denom = pn.length();
  return numer/denom;
}
function testDistancePointToPlane() {
  const pn = new THREE.Vector3(0,0,1);
  const p0 = new THREE.Vector3(0,0,0);
  const x = new THREE.Vector3(10,100,37);
  const d = distancePointToPlaneThroughOrigin(pn,x);
  console.assert(near(d,37));
}

// Compute the points where a circle
// centered at c with normal cn of radius r
// intersects a plane of norma pn containing point p
// http://mathforum.org/library/drmath/view/69136.html
// (This appears to be missing an oppportunity!
function computeRotationIntoPlane(A,c,cn,r,pn,p) {
  console.assert(vnear(p,new THREE.Vector3(0,0,0)));
  // First we need to construct a vector perpendicular to cn
  // of length r...we will take the cross product with a unit vector...
//  const A = new THREE.Vector3(0,3/5 ,2);
  U = new THREE.Vector3().crossVectors(cn,A);
  if (vnear(U,new THREE.Vector3(0,0,0))) {
    conosole.log("failed to produce a good vector!",cn,A,U);
    debugger;
  }
  // Check that we got something...
  console.assert(!vnear(U,new THREE.Vector3(0,0,0)));
  V = new THREE.Vector3().crossVectors(cn,U);
  // Now U and V are perpendicular in the plane of the ring...
  // Set their length to r....
  U.normalize().multiplyScalar(r);
  V.normalize().multiplyScalar(r);
  const UdotN = U.dot(pn);
  const VdotN = V.dot(pn);
  const Q = Math.sqrt(UdotN**2 + VdotN**2);
  if (near(Q,0)) {
    console.log("Q is near 0, a numerical instability!");
    // This means that U and V are both perpendicular to the plane,
    // so any angle keeps us in the plane I guess
    debugger;
  }

  // It is very unclear what the signs of this should be...
  const pndc_q = Math.max(-1,Math.min(1,pn.dot(c)/Q));
  const t = Math.asin(-pndc_q) - Math.asin(UdotN/Q);
  if (isNaN(t)) {
    debugger;
  }

  var C = c.clone();
  const Cprime = C.clone();

//  console.log("t ",t * 180 / Math.PI);

  C.add(U.clone().multiplyScalar(Math.cos(t)));
  Cprime.add(U.clone().multiplyScalar(Math.cos(-t)));

  // I am worried about the sign of this coefficient;
  // I do not know how to assure it.

  C.add(V.clone().multiplyScalar(Math.sin(-t)));
  Cprime.add(V.clone().multiplyScalar(Math.sin(t)));

  var Cf;
  var tf;
  if (near(pn.dot(C),0)) {
    Cf = C;
    tf = t;
  } else if (near(pn.dot(Cprime),0)) {
    Cf = Cprime;
    tf = -t;
  } else {
    console.log("Catastrophe---we couldn't get to the plane!");
    console.log("C,Cprime",C,Cprime);
    console.log("t",t * 180/Math.PI);
//    debugger;
    // Possibly do to the "balance" problem we can't
    // actually reach the plane. We should return
    // the value that is closest?
    const d = distancePointToPlaneThroughOrigin(pn,C);
    const dprime = distancePointToPlaneThroughOrigin(pn,Cprime);
    if (Math.abs(d) < Math.abs(dprime)) {
      Cf = C;
      tf = t;
    } else {
      Cf = Cprime;
      tf = -t;
    }
  }

  // What we really want to do here is to assert that
  // C truly is in the plane <pn,p>

//  console.assert(near(pn.dot(Cf),0));
//   if (!near(pn.dot(Cf),0)) {
//     console.log(" xxx ", pn.dot(c)/Q);
//     console.log("angle ",tf * 180/Math.PI);
//     console.log("pn,Cf",pn,Cf);
// //    debugger;
//     return null;

//   }
  return [Cf,tf];
}


function computeRotationAngleIntoPlane(A,c,cn,r,pn,p) {
  if (near(r,0)) {
    return 0;
  }
  const [P,t] = computeRotationIntoPlane(A,c,cn,r,pn,p);
  if (!P) return null;
  // now we want to copute the angle between cP and cA.
  const Ac = new THREE.Vector3().subVectors(A,c);
  console.assert(near(r,Ac.length()));
  if (!near(r,Ac.length())) {
    debugger;
  }
  const Pc = new THREE.Vector3().subVectors(P,c);
  const angle = Ac.angleTo(Pc);
  // Now would should be able to rotate about cn move A into P.
  // If not, we are doing something wrong!
  // now we will check that the rotation about the axis really
  // places A in the plane
  var trans = new THREE.Matrix4();
  trans.makeTranslation(-c.x,-c.y,-c.z);
  var trans_i = new THREE.Matrix4();
  trans_i.makeTranslation(c.x,c.y,c.z);
  var rt = new THREE.Matrix4();
  // Now here is a problem....what should the sign of this angle
  // be?  I am sure it is not computed as a singed angle!

  rt.makeRotationAxis(cn,angle);

  var finalMatrix = new THREE.Matrix4().identity();
  finalMatrix.premultiply(trans);
  finalMatrix.premultiply(rt);
  finalMatrix.premultiply(trans_i);
  const Aprime = A.clone().applyMatrix4(finalMatrix);

  const nangle = -angle;

  var nrt = new THREE.Matrix4();
  nrt.makeRotationAxis(cn,nangle);

  var nfinalMatrix = new THREE.Matrix4().identity();
  nfinalMatrix.premultiply(trans);
  nfinalMatrix.premultiply(nrt);
  nfinalMatrix.premultiply(trans_i);
  const nAprime = A.clone().applyMatrix4(nfinalMatrix);

//  console.assert(near(nAprime.x,0));
//  console.log("A,Aprime,nAprime",A,Aprime,nAprime);
  if (near(pn.dot(nAprime),0)) {
    return nangle;
  } else {
    return angle;
  }

}

// At present this seems to be rotating
// completely around the z axis (not taking the first
// opportunity reach the first YZ plane
function testComputeRotationAngleIntoPlane() {
  const NUMTEST = 5;
  // We will places the axis in the Z direction,
  // at a point a little to the left and above
  const axis = new THREE.Vector3(0,0,1);
  const c = new THREE.Vector3(1,1,1);
  const c2 = new THREE.Vector3().addVectors(c,axis);

  // We will use the Y-Z plane
  const pn = new THREE.Vector3(1,0,0);
  const p = new THREE.Vector3(0,0,0);
  // Now we will var y the points in the Y direction
  for(var i = 0; i < NUMTEST; i++) {
    const A = new THREE.Vector3(1,i+2,1);
    const r = pointToLineDist(A,c,c2);
    console.log("i,r",i,r);
    console.log("A",A);
    console.log("c,c2",c,c2);
    var angle = computeRotationAngleIntoPlane(A,c,axis,r,pn,p);
    if (angle > Math.PI/2) {
      angle = Math.PI - angle;
    }
    console.log("angle ", angle * 180 / Math.PI);
    // now we will check that the rotation about the axis really
    // places A in the plane
    var trans = new THREE.Matrix4();
    trans.makeTranslation(-c.x,-c.y,-c.z);
    var trans_i = new THREE.Matrix4();
    trans_i.makeTranslation(c.x,c.y,c.z);
    var rt = new THREE.Matrix4();
    rt.makeRotationAxis(axis,angle);

    var finalMatrix = new THREE.Matrix4().identity();
    finalMatrix.premultiply(trans);
    finalMatrix.premultiply(rt);
    finalMatrix.premultiply(trans_i);
    const Aprime = A.clone().applyMatrix4(finalMatrix);

    console.assert(near(Aprime.x,0));
    console.log("A,Aprime",A,Aprime);
  }
}

function testComputeRotationAngleIntoPlaneSpecific() {
  const axis = new THREE.Vector3(-0.2817571131513068,-0.7663416965396759,-0.5773502691896258);

  const c = new THREE.Vector3(-0.19923236535683433,-0.5418854103292079,-0.4082482904638629);

//  const c = new THREE.Vector3(1,1,1);
  const c2 = new THREE.Vector3().addVectors(c,axis);

  // We will use the Y-Z plane
  const pn = new THREE.Vector3(1,0,0);
  const p = new THREE.Vector3(0,0,0);
  // Now we will var y the points in the Y direction

  const A = new THREE.Vector3(0.1704379832499971,-0.9853684051488933,2.220446049250313e-16);

  const r = pointToLineDist(A,c,c2);
  console.assert(near(r,0.7071067811865477));
  console.log("A",A);
  console.log("c,c2",c,c2);
  var angle = computeRotationAngleIntoPlane(A,c,axis,r,pn,p);
  // This is likely to be wrong half the time, I guess.
//  if (angle > Math.PI/2) {
//    angle = Math.PI - angle;
  //  }
  console.log("angle ", angle * 180 / Math.PI);
  // now we will check that the rotation about the axis really
  // places A in the plane
  var trans = new THREE.Matrix4();
  trans.makeTranslation(-c.x,-c.y,-c.z);
  var trans_i = new THREE.Matrix4();
  trans_i.makeTranslation(c.x,c.y,c.z);
  var rt = new THREE.Matrix4();
  rt.makeRotationAxis(axis,angle);

  var finalMatrix = new THREE.Matrix4().identity();
  finalMatrix.premultiply(trans);
  finalMatrix.premultiply(rt);
  finalMatrix.premultiply(trans_i);
  const Aprime = A.clone().applyMatrix4(finalMatrix);

//  console.assert(near(Aprime.x,0));
//  console.log("A,Aprime",A,Aprime);
}
function testComputeRotationIntoPlane() {
  const c = new THREE.Vector3(1,1,1);
  const cn = new THREE.Vector3(0,0,1).normalize();
  const r = 3;
  const p = new THREE.Vector3(0,0,0);
  const pn = new THREE.Vector3(0,1,0).normalize();
  const A = new THREE.Vector3(1,2,3);

  var [pnt,t] = computeRotationIntoPlane(A,c,cn,r,pn,p);
  console.log("Point: ",pnt);
  // The should be two points...
//  console.assert(pnts.length == 2);
  // The points are in the plane...
  console.assert(near(pnt.z,1));
//  console.assert(near(pnts[1].y,0));
  // the points have x values of 1/2
//  console.assert(near(Math.abs(pnt.x),1));
  console.assert(near(Math.abs(pnt.y),0));

}

function testComputeRotationIntoPlane2() {
  const c = new THREE.Vector3(0.2763168314822806,
                              -0.7646590249937488,
                              -0.5749149571305298);
  const cn = new THREE.Vector3(0.2774872962675397,
                               -0.7678980837824754,
                               -0.5773502691896258).normalize();
  const r = 3;
  const p = new THREE.Vector3(0,0,0);
  const pn = new THREE.Vector3(0,1,0).normalize();
  const A = new THREE.Vector3(1,2,3);
  var [pnt,t] = computeRotationIntoPlane(A,c,cn,r,pn,p);

  console.assert(near(pnt.dot(pn),0));
}

// compute the angle of rotation needed to move the from vector
// to the to vector by rotating about the axis.
// Assume the axis is specified by a vector from the origin,
// as are from and to.
function computeRotation(axis,from,to) {
  // My plan is to use "angleBetween" on the rejections...
  const an = axis.clone().normalize();
  const rfrom = vectorRejectionOnto(from,axis);
  const rto = vectorRejectionOnto(to,axis);
  // hear I want to make sure the cross product of from and to is parallel to axis.
  const c = new THREE.Vector3().crossVectors(rto,rfrom).normalize();
  console.log("axis, cross",an,c);
  return rfrom.angleTo(rto);
}

function testComputeRotation() {
  const axis = new THREE.Vector3(1,0,1).normalize();
  const from = new THREE.Vector3(2,0,1).normalize();
  const to = new THREE.Vector3(1,0,2).normalize();
  const theta = computeRotation(axis,from,to);
  console.log("theta = ", theta * 180/Math.PI);
  const rt = new THREE.Matrix4();
  rt.makeRotationAxis(axis,theta);
  const fc = from.normalize().clone().applyMatrix4(rt);
  console.assert(vnear(fc,to.clone().normalize()));
}
function testComputeRotationFull() {
  const NUMTEST = 5;
  const DELTA = Math.PI / 13;
  const axis = new THREE.Vector3(1,2,3).normalize();
  for(var i = 0; i < NUMTEST; i++) {
    const from = new THREE.Vector3(2,1,i).normalize();
    const theta = DELTA * i;
    const rt = new THREE.Matrix4();
    rt.makeRotationAxis(axis,theta);
    const to = from.clone().applyMatrix4(rt);
    const theta_test = computeRotation(axis,from,to);
    console.assert(near(theta,theta_test));
  }
}

// Return tau to achieve this psi
//
function computeTauToAchievePsi(NBu,NCu,psi,initTau = 0,L = 1,debug = false) {
  // First, we copy the old prism in exactly the same position

  var abs = new AbstractPrism(L,NBu,NCu);
  if (!abs) {
    debugger;
  }
  // In fact this algorithm more or less assumes these are in balance...

  const old = CreatePrism(abs,PRISM_FACE_RATIO_LENGTH);

  // By setting tau = 0 here, I am doing something weird....
  // I'm not sure this is legitimate. My purpose is
  // really to put the thing in balance...I might need a special
  // routine for that...
  const [A,psi_balance,D] = AfromLtauNbNc(L,initTau,NBu,NCu);
  if (debug) console.log("INIT_TAU",initTau);

  const rbal = new THREE.Matrix4();
  rbal.makeRotationAxis(new THREE.Vector3(0,0,1),psi_balance);

  applyMatrix4ToPrism(old,rbal);
  old.p.Nb.applyMatrix4(rbal);
  old.p.Nc.applyMatrix4(rbal);

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

  let joinToC = false;

  console.assert(!vnear(new THREE.Vector3(0,0,0),
                        (joinToC ? nu.p.Nc : nu.p.Nb)));

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
  console.assert(!vnear(new THREE.Vector3(0,0,0),(joinToC ? nu.p.Nc : nu.p.Nb)));

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

  console.assert(!vnear(new THREE.Vector3(0,0,0),(joinToC ? nu.p.Nc : nu.p.Nb)));


  // http://mathforum.org/library/drmath/view/69136.html

  // Now to use this we need parameters.
  // The axis is the norm of circle, cn.
  const axis = nu.p.Nb;
//  axis = NBu; // ???
  const cn = axis;
  // The plane we want to rotate into the is X-Z plane (for psi = 0)
  const pn = new THREE.Vector3(1,0,0);
  const p = new THREE.Vector3(0,0,0);
  // What is the center of the circle?
  // We can pick any point on the axis...
  const c = axis.clone().normalize();
  // Now how do we find the radius?
  const p1 = p;
  const p2 = axis;
  const P = nu.b.clone(); // or A?
  const r = pointToLineDist(P,p1,p2);
  const M = closestPoint(P,p1,p2);

  var compensationAngle = computeRotationAngleIntoPlane(P,M,cn,r,pn,p);

  if (debug) {
    console.log("axis",axis);
    console.log("P ",P);
    console.log("M ",M);
    console.log("r ",r);
    console.log("cn,r,pn,p ",cn,r,pn,p);
    console.log("compensationAngle", compensationAngle * 180/Math.PI);
  }

  if (compensationAngle === null) {
    return [null,null];
  }

  const rt = new THREE.Matrix4();
  rt.makeRotationAxis(nu.p.Nb,compensationAngle);

  const [computedTau,u] = computeThetaU(rt);
  if (debug) {
    console.assert(near(computedTau,compensationAngle));
    console.log("computedTau, compensationAngle",computedTau * 180/Math.PI,compensationAngle * 180/Math.PI);
  }

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
  return [computedTau,finalMatrix]
}

function iterativeComputeTauToAchievePsi(NBu,NCu,psi,initTau = 0,L = 1,MAX_ITERS = 10,debug = false) {

  var tau = initTau;
  for(var i = 0; i < MAX_ITERS; i++) {
    var [tau,rt] = computeTauToAchievePsi(NBu,NCu,psi,tau,L,debug);
  }
  return tau;
}

function bruteForceTauComputation(Nb,Nc) {
  const NUM = 360;
  const DELTA = (2 * Math.PI) / 360;
  const L0 = 1;
  const P0 = new AbstractPrism(
    L0,
    Nb,
    Nc,
    null);

  const PRISM_FACE_RATIO_LENGTH = 1/2;
  var p_temp = CreatePrism(P0,PRISM_FACE_RATIO_LENGTH);

  let B = new THREE.Vector3(0,0,-L0/2);
  let C = new THREE.Vector3(0,0,L0/2);

  var min_tightness_tau = -1;
  var max_tightness_tau = -1;
  var min_tightness = 99999;
  var max_tightness = -1;

  for(var tau = 0; tau < 2 * Math.PI; tau += DELTA) {
//    console.log("TAU",tau * 180 / Math.PI);

    var [resK,resM,local_p_i,rt] =
      computeInternal(L0,B,C,tau,p_temp,Nb,Nc);

    const r = resM[0];
    const d = resM[2];
    const tightness = Math.abs(d / r);

    if (tightness > max_tightness) {
      max_tightness = tightness;
      max_tightness_tau = tau;
    }
    if (tightness < min_tightness) {
      min_tightness = tightness;
      min_tightness_tau = tau;
    }
  }
  return [[max_tightness,max_tightness_tau],
          [min_tightness,min_tightness_tau]];
}

function testBruteForceTauComputation() {
  const NUM_TEST = 4;
  for(var i = 0; i < NUM_TEST; i++) {
    var Nbx = -1 + i*(2/NUM_TEST);
    for(var j = 0; j < NUM_TEST; j++) {
      var Nby = -1 + j*(2/NUM_TEST);
      for(var k = 0; k < NUM_TEST; k++) {
        var Ncx = -1 + k*(2/NUM_TEST);
        for(var l = 0; l < NUM_TEST; l++) {
          var Ncy = -1 + l*(2/NUM_TEST);
          if ((Nbx != Ncx) || (Nby != Ncy)) {
            let Nb = new THREE.Vector3(Nbx,Nby,-1);
            let Nc = new THREE.Vector3(Ncx,Ncy,1);
            Nb.normalize();
            Nc.normalize();
            const [[max_tightness,max_tightness_tau],
                   [min_tightness,min_tightness_tau]]
                  = bruteForceTauComputation(Nb,Nc);
            console.assert(max_tightness > 1);
            // because we are using only degree resolution, this is a bit course
            console.assert(near(min_tightness,0,0.10));
          }
        }
      }
    }
  }
}

function testIterativeComputeTauToAchievePsi() {

  {
    const NBu = new THREE.Vector3(-1/2,-1/2,-1).normalize();
    const NCu = new THREE.Vector3(0,0,1).normalize();

    const tau = iterativeComputeTauToAchievePsi(NBu,NCu,0,0,1,10,false);
    console.assert(near(tau,0));
  }

  {
    const NBu = new THREE.Vector3(-1/2,-1/2,-1).normalize();
    const NCu = new THREE.Vector3(1,0,1).normalize();

    const tau = iterativeComputeTauToAchievePsi(NBu,NCu,0,0,1,10,false);
    const tau_d = tau * 180 / Math.PI;
    console.assert(tau_d > 11.7);
    console.assert(tau_d < 11.8);
  }

}

function testComputeTauToAchievePsi() {

  const NUM_TEST = 2;
  const L = 1;
  B = new THREE.Vector3(0,0,-1/2);
  C = new THREE.Vector3(0,0,1/2);
  var n = 0;
  var e = 0;

  for(var i = 0; i < NUM_TEST; i++) {
    var Nbx = -1 + i*(2/NUM_TEST);
    for(var j = 0; j < NUM_TEST; j++) {
      var Nby = -1 + j*(2/NUM_TEST);
      for(var k = 0; k < NUM_TEST; k++) {
        var Ncx = -1 + k*(2/NUM_TEST);
        for(var l = 0; l < NUM_TEST; l++) {
          var Ncy = -1 + l*(2/NUM_TEST);
          if (!((Nbx == 0) && (Nby == 0) && (Ncx == 0) && (Ncy ==0))) {
            let Nb = new THREE.Vector3(Nbx,Nby,-1);
            let Nc = new THREE.Vector3(Ncx,Ncy,1);
            Nb.normalize();
            Nc.normalize();
            if (!vnear(Nb,Nc.clone().multiplyScalar(-1))) {
              const initTau = Math.PI/36;
              const [tau,rt] = computeTauToAchievePsi(Nb,Nc,0,initTau);
              if (near(tau, Math.PI/2)) {
//                                debugger;
              }
              if (tau === null) {
                console.log("FAILED: ",Nb,Nc);
                e++;
              } else {
                n++;

                console.log(" computed tau for 0:",Nb,Nc,tau*180/Math.PI);
                // now test that rt (being a rotation of Tau)
                // produces a toriod: rt(C) = D, rt(D.x,D.y,-D.z) = B
                let D = C.clone().applyMatrix4(rt);
                let A = new THREE.Vector3(-D.x,D.y,-D.z);
                let Ap = A.clone().applyMatrix4(rt);
                console.assert(vnear(Ap,B));
              }
            }
          }
        }
      }
    }
  }
  console.log("ERRORS : ",e,n);
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

    // TODO: explacing this with a dot-product is probably
    // nicer.
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
function testRegularTetsPointAxis()
{
  let L0 = 1/3;
  let tetAngle = Math.atan(Math.sqrt(2)/2);
  let NB1 = new THREE.Vector3(0,-Math.sin(tetAngle),-Math.cos(tetAngle));
  let NC1 = new THREE.Vector3(0,-Math.sin(tetAngle),Math.cos(tetAngle));
  let tau = 2 * Math.PI /3;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let B = new THREE.Vector3(0,0,-L0/2);
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
function testPointAxisYTorus()
{
  let N = 17;
  let angle = Math.PI/17;
  for(let a = angle; a < Math.PI*2; a += angle) {
    testPointAxisYTorusAux(angle);
  }
}
function testPointAxisYTorusAux(angle)
{
  let L0 = 2;
  let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
  let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
  let tau = 0;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let B = new THREE.Vector3(0,0,-L0/2);
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
function testPointAxisTau180()
{
  let L0 = 2;
  let angle = Math.PI/7;
//  let angle = 0;
  let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
  // We add in a little here so we are not a torus
  let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
  let tau = Math.PI;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  console.assert(near(r,0));
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  let zero = new THREE.Vector3(0,0,0);
  console.assert(!near(da,0,0.1));
  console.assert(!isNaN(theta));
  console.assert(near(Ba.x,0));
  console.assert(near(Ba.y,0));
  let zaxis = new THREE.Vector3(0,0,1);
  let nzaxis = new THREE.Vector3(0,0,-1);
  console.assert(vnear(H,zaxis) || vnear(H,nzaxis));
  console.assert(phi >= 0);
}

function testPointAxisTauNeg180()
{
  let L0 = 2;
  let angle = Math.PI/7;
//  let angle = 0;
  let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
  // We add in a little here so we are not a torus
  let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
  let tau = -Math.PI;
  let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  console.assert(near(r,0));
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  let zero = new THREE.Vector3(0,0,0);
  console.assert(!near(da,0,0.1));
  console.assert(!isNaN(theta));
  console.assert(near(Ba.x,0));
  console.assert(near(Ba.y,0));
  let zaxis = new THREE.Vector3(0,0,1);
  console.assert(vnear(H,zaxis));
  console.assert(phi >= 0);
}

function testPointAxisFull()
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
        testPointAxisFullAux(tau, NB1, NC1);
      }
    }
  }
}
function testPointAxisFullAux(tau,NB1,NC1)
{
  let L0 = 2;
  let res = compareMethods(L0,tau,NB1,NC1);
  let r = res[0];
  let theta = res[1];
  let da = res[2];
  let c = res[3];
  let phi = res[4];
  let H = res[5];
  let Ba = res[6];
  let B = res[7];
  console.assert(" B = ", B);
  console.assert(!isNaN(theta));
  try {
//    console.log(B.x);
  } catch (e) {
    debugger;
  }
  testAxis(B,Ba,H,da);
}

// This is to test the Flat, Zig-Zag Case
function testPointAxisFlat()
{
  let base = Math.PI/7;
  const N = 10;
  for(let i = -(N-1); i < N; i++) {
    let angle = -((base) * i) / ((N-1)/N);
    testPointAxisFlatAux(angle);
  }
}
function compareMethods(L,tau,NB1,NC1) {
  let AR = AfromLtauNbNc(L,tau,NB1,NC1);
  let A = AR[0];
  let B = new THREE.Vector3(0,0,-L/2);
  let C = new THREE.Vector3(0,0,L/2);
  let D = new THREE.Vector3(-A.x,A.y,-A.z);
  let R = AR[2];
  let res = PointAxis(L,A);
  [r,thetaP,da,chord,phi,u,Ba] = computeThetaAxisFromMatrix4(L,R,B);
  const BC = new THREE.Vector3().subVectors(C,B);

  if (res[1] == null) {
    // In this case, PointAxis is undefined and we cannot compare!!!
    return [r,thetaP,da,chord,phi,u,Ba,B];
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

  console.assert(vnear(Ba,res[6]));
  if (!(vnear(Ba,res[6]))) {
    debugger;
  }


  // TODO: Apparently computeThetaAxis is not computeing Ba correctly!!!
 // console.assert(vnear(Ba,res[6]));

  return [r,thetaP,da,chord,phi,u,Ba,B];
}
function testPointAxisFlatAux(angle)
{
  let L0 = 2;

  let NB1 = new THREE.Vector3(-Math.sin(angle),0,-1);
  // We add in a little here so we are not a torus
  let NC1 = new THREE.Vector3(Math.sin(angle),0,1);
  let tau = Math.PI;
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

function compute_da(u,B,C) {
  // Now I believe da is indepednet of L and any point.
  const BC = new THREE.Vector3().subVectors(C,B);

  // The computation of this should be removed to a different function.
  // That will allow us to remove the arbitrary points out, I think.
  // Now BC.u sould give us da...
  const unorm = u.clone().normalize();
  const da_scalar_project = BC.dot(unorm);
  return da_scalar_project;
}

function compute_da_X(R,u,B,C) {
  // Now I believe da is indepednet of L and any point.
  const BC = new THREE.Vector3().subVectors(C,B);

  // The computation of this should be removed to a different function.
  // That will allow us to remove the arbitrary points out, I think.
  // Now BC.u sould give us da...
  const unorm = u.clone().normalize();
  const da_scalar_project = BC.dot(unorm);

  // This is all test code.
  {
  const ele = R.elements;
  const a = ele[0];
  const e = ele[1];
  const i = ele[2];

  const b = ele[4];
  const f = ele[5];
  const j = ele[6];

  const c = ele[8];
  const g = ele[9];
  const k = ele[10];

  const d = ele[12];
  const h = ele[13];
  const l = ele[14];

  const x = (a+b+c+d);
  const y = (e+f+g+h);
  const z = (i+j+k+l);
  console.assert(a+b+c+x,1);
  const BCp = new THREE.Vector3(x-1,y-1,z-1);
  console.assert(vnear(BC,BCp));
  const da = BCp.dot(unorm);
  console.assert(near(da,da_scalar_project));
  }
  return da_scalar_project;
}

// https://math.stackexchange.com/questions/1053105/know-if-a-4x4-matrix-is-a-composition-of-rotations-and-translations-quaternions
function isRigidTransformation(R) {
  const det = R.determinant()
  if (near(det,-1)) {
    console.log("R.det :",R.determinant());
    return false;
  }
  if (det < 0) {
    console.log("R is a reflection! :",det);
    return false;
  }
  console.assert(checkNonScaling(R));
  if (!checkNonScaling(R)) {
    debugger;
    return false;
  }
  const s = R.elements[15];
  console.assert(near(s,1));

  let xA = new THREE.Vector3();
  let yA = new THREE.Vector3();
  let zA = new THREE.Vector3();
  R.extractBasis(xA,yA,zA);
  const A = new THREE.Matrix3();
  A.setFromMatrix4(R);
  const At = A.clone().transpose();
  const T = new THREE.Matrix3().multiplyMatrices(A,At);
  console.assert(near(T.elements[0],1));
  console.assert(near(T.elements[1],0));
  console.assert(near(T.elements[2],0));

  console.assert(near(T.elements[3],0));
  console.assert(near(T.elements[4],1));
  console.assert(near(T.elements[5],0));

  console.assert(near(T.elements[6],0));
  console.assert(near(T.elements[7],0));
  console.assert(near(T.elements[8],1));
  if (!near(T.elements[0],1)) {
    return false;
  }
  if (!near(T.elements[1],0)) {
    return false;
  }
  if (!near(T.elements[2],0)) {
    return false;
  }
  if (!near(T.elements[3],0)) {
    return false;
  }
  if (!near(T.elements[4],1)) {
    return false;
  }
  if (!near(T.elements[5],0)) {
    return false;
  }
  if (!near(T.elements[6],0)) {
    return false;
  }
  if (!near(T.elements[7],0)) {
    return false;
  }
  if (!near(T.elements[8],1)) {
    return false;
  }

  return true;
}

function computeThetaU(R) {
  console.assert(isRigidTransformation(R));
  // now attempt to recover parameters from the rotation matrix....
  const tr = TraceMatrix4(R);
  const val = (1/2) * (tr - 1);

  // I don't know what to do if this condition holds
  console.assert(Math.abs(val) < 1.01);

  // Possibly this needs to be signed based on something
  // else we compute (A.x) when val = -1.
  const thetaP = (near(val,-1)) ? Math.PI :  Math.acos(val);

  // var PROBLEM = false;
  // if (near(val,-1)) {
  //   console.log("WARNING!!! TRACE DEGENERATE");
  //   PROBLEM = true;
  //   debugger;
  // }

  console.assert(!isNaN(thetaP));

  // We need an arbitray point.. technically we should check
  // this is not on the axis, but that is unlikely....
  const B = new THREE.Vector3(0,4,-4/2);
  const C = B.clone().applyMatrix4(R);
  let u = null;
  if (thetaP == Math.PI || thetaP == 0) {
    // In this case our normal method is numerically unstable
    // (it defines an axis vector of zero length.). However,
    // in this "zig-zag" case, u is easily computed as R^2(P) - P,
    // where P is an arbitrary point.
    // Sine we need B and C in other cases, we will compute D
    // as a further application of R to C.
    const D = C.clone().applyMatrix4(R);
    let Ri = new THREE.Matrix4().getInverse(R);
    const A = B.clone().applyMatrix4(Ri);

    u = new THREE.Vector3().subVectors(D,B);
    u.normalize();
    if (A.x > 0) {
      u.multiplyScalar(-1);
    }
    return [thetaP,u];
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
    // Note: I am defining my axis with a right hand rule from B to C,
    // and with a negative travel to represent the other direction.

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
    return [thetaP,u];
  }
}


// https://www.seas.upenn.edu/~meam520/notes02/EulerChasles4.pdf
// Note: This may solve my problem: https://en.wikipedia.org/wiki/Screw_axis#Computing_a_point_on_the_screw_axis
function computeRChord(L,thetaP,da) {
  let chord = ChordFromLDaxis(L,da);
  if (near(0,chord,1e-3)) {
    return [0,0];
  } else {
    if (near(0,Math.sin(thetaP/2))) {
      console.log("WARNING!!!");
      debugger;
    }
    const r = chord / (2 * Math.sin(thetaP/2));
    return [r,chord];
  }
}

function computePointIndependentParameters(L,R) {
  [thetaP,u] = computeThetaU(R);

  // We need an arbitray point.. technically we should check
  // this is not on the axis, but that is unlikely....
  // TODO: See if this can be determined without using an aribitrary point..
  const B = new THREE.Vector3(1,1,1);
  const C = B.clone().applyMatrix4(R);
//  const da = compute_da(u,B,C);
  const da = compute_da_X(R,u,B,C);

  if (near(da,L)) {
    [r,chord] = [0,0];
  } else {
    [r,chord] = computeRChord(L,thetaP,da);
  }
  return [r,thetaP,da,chord,u];
}

function computePhi(L,R,B) {
  const C = B.clone().applyMatrix4(R);
  const BC = new THREE.Vector3().subVectors(C,B);

  console.assert(near(L,BC.length()));
  // This assertion most likely indicates B does not match the rotation matrix
  if (!near(L,BC.length())) {
    debugger;
  }
  //  let phi = (near(da,L)) ? 0 : Math.acos(da / L);
  let phi;
  if (da/L < -1) {
    phi = -Math.PI;
  } else if (da/L > 1) {
    phi = 0;
  } else {
    phi = (near(Math.abs(da),Math.abs(L))) ? 0 : Math.acos(da / L);
  }
  if (isNaN(phi)) debugger;
  return phi;
}


// Possibly this is the same math, or a different approach to
// the formula given: https://en.wikipedia.org/wiki/Screw_axis
// C = b x d - b b x (b x d) / (2b . b).
// where C is a point on the axis, d is the displacement,
// and b is the Rodrigues' vector.
// Questions remain:
// Is one better than the other?
// Is does that computed point correspond to a perpendicular to B, as this
// produces (which is a major advantage) (though I suppose that approach doesn't
// require us to namea point at all.
function computeAxisPoint(r,chord,u,R,B) {
  const C = B.clone().applyMatrix4(R);
  const BC = new THREE.Vector3().subVectors(C,B);
  const Mp = new THREE.Vector3().addVectors(B,C);

  Mp.multiplyScalar(1/2);
  const Q = BC.clone()
  Q.cross(u);

  const ql = Math.sqrt(r**2 - (chord/2)**2);
  // ql is effectively the radius of the helix measure from the midpoint
  // of a member, which is slightly less that the joint radius
  Q.setLength(-ql);
  const Ba = Mp.clone().add(Q);

  // At this point Ba is on the axis, but does not match B yet...
  const uuhalf = u.clone().multiplyScalar(-da/2);
  const nBa = new THREE.Vector3().addVectors(Ba,uuhalf);
  return nBa;
}

function computeThetaAxisFromMatrix4(L,R,B) {
  [r,thetaP,da,chord,u] = computePointIndependentParameters(L,R);

  // Can we note compute \phi now that we know $r$?
  let phi = computePhi(L,R,B);

  let nBa = computeAxisPoint(r,chord,u,R,B);

  return [r,thetaP,da,chord,phi,u,nBa,B];
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

  const B = new THREE.Vector3(4,4,4);
  const C = B.clone().applyMatrix4(R);
  const L = B.distanceTo(C);
  [r,thetaP,da,chord,phi,u] = computeThetaAxisFromMatrix4(L,R,B);

  // now check that the parameters match...
  console.assert(near(theta,thetaP));
  if (!vnear(axis,u)) {
    console.log("axis,U",axis,u);
  }
  console.assert(vnear(axis,u));
}

// The purpose of this test is to the computation of the rotation matrix
//
function testPointOnAxisRotationMatrix() {
  const theta = Math.PI/11;
  // now create a rotation matrix based on this angle...
  const zaxis = new THREE.Vector3(0,0,1);
  const B = new THREE.Vector3(0,1,0);
  const O = new THREE.Vector3(0,0,0);
  const re = B.length();

  function compute_from_arbitrary_axis(B,R) {
    const C = B.clone().applyMatrix4(R);
    const L = B.distanceTo(C);
    return  computeThetaAxisFromMatrix4(L,R,B);
  }

  const N = 4;
  // now let's rotate though various Euler angles and test
  for(var a = 1; a < N+1; a++) {
    for(var b = 1; b < N+1; b++) {
      for(var c = 1; c < N+1; c++) {
        const axis_pre = new THREE.Vector3(a,b,c).normalize();
        const pre = new THREE.Matrix4().makeRotationAxis(axis_pre,theta)
//        const T = new THREE.Matrix4().makeTranslation(a,b,c);
//        const RT = R.multiply(T);
        const axis = zaxis.clone().applyMatrix4(pre);
        const Bpre = B.clone().applyMatrix4(pre);
        const R = new THREE.Matrix4().makeRotationAxis(axis,theta)
        const RT = R;
        console.assert(near(Bpre.length(),1));
        [r,thetaP,da,chord,phi,u,Ba,Bx] = compute_from_arbitrary_axis(Bpre,RT);
        console.assert(near(re,r));
        console.assert(near(thetaP,theta));
        console.assert(vnear(u,axis));
        console.assert(vnear(Ba,O));
      }
    }
  }

  console.log("POINT ON AXIS:",Ba);

}



function runChaslesTests() {
  testTrace();
  testAbilityToGetScrewAxisFromRotationMatrix();
}

function runUnitTests() {
  testClosestPoint();
  testComputeRotation();
  testComputeRotationFull();

  testComputeRotationIntoPlane();
  testComputeRotationIntoPlane2();

  runChaslesTests();
  testThreeIsRightHanded();
  testPointOnAxisRotationMatrix();

  testAfromLtauNbNnKeepsYPure();
  testAfromLtauNbNnContinuityOfTau();
  testRegularTetsPointAxis();
  testPointAxisYTorus();
  testPointAxisFull();
  testPointAxisFlat()
  testPointAxisTau180();
  testPointAxisTauNeg180();
  testAfromLfailureCase1();
  testAfromLtauMultiple();
  testComputeTauToAchievePsi();

  testComputeRotationAngleIntoPlane();
  testComputeRotationAngleIntoPlaneSpecific();
  testIterativeComputeTauToAchievePsi();

  testDistancePointToPlane();

  testBruteForceTauComputation();
}
