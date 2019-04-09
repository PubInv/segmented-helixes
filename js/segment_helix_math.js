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

function QFromRhoOmega(L,rho,omega) {
    return L/(Math.sqrt(1 + Math.tan(rho) ** 2 + Math.tan(omega) ** 2));
}
function ChordFromLD(L,di) {
    return Math.sqrt(L ** 2 - di**2);
}
function ThetaFromRC(r,c) {
    return Math.abs(2 * Math.asin(c/(2*r)));
}

function AbstractPrism(L,Nb,Nc) {
    this.L = L;
    this.Nb = Nb.clone().normalize();
    this.Nc = Nc.clone().normalize();
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
        let theta = 2 * Math.asin(C / (2 * R));
        if (isNaN(theta)) {
            debugger;
        }
        return 2 * Math.asin(C / (2 * R));
    }
}

// D is a the point, L is the length.
// The return values are:
// r -- the radius
// theta -- the signed rotation of one point to the next around the axis.
// da -- signed travel along the axis when going from one segment to the next
// c -- the chord length along the axis
// phi -- the angle of the Helix axis against a segment
// H -- the vector of the axis of the helix
function KahnAxis(L,D) {
    let x = D.x;
    let y = D.y;
    let z = D.z;
    let A = new THREE.Vector3(-x,y,-z);

    let B = new THREE.Vector3(0,0,-L/2);
    let C = new THREE.Vector3(0,0,L/2);

    let Cb = B.clone();
    Cb.add(D);
    Cb.multiplyScalar(1/2);
    Cb.sub(C);

    let Bb = A.clone();
    Bb.add(C);
    Bb.multiplyScalar(1/2);
    Bb.sub(B);


    // If Bb and Cb have zero length, the cross product
    // is undefined, and the axis of the helix may be
  // taken to be anywhere (is not uniquely defined.)
  // This is the case of a "straight line".
    if (Bb.length() == 0) {
        // We will take H to be the vector point from B to C,
        // use a zero radius.
        let H = new THREE.Vector3(0,0,L);
        let da = L;
        let r = 0;
        let c = 0;
        // I need to define this more clearly....
        let phi = Math.PI;
        let theta = 0;
        return [r,theta,da,c,phi,H];
    } else {
        let CmB = C.clone();
        CmB.sub(B);
        if (near(y,0,1e-4)) { // This is the Zig-Zag case
            let CmA = C.clone();
            CmA.sub(A);
            let H = CmA.clone();
            H.normalize();
//            H.multiplyScalar(-1);
            // I think I need to understand why this is negative
            let da = (CmA.x >= 0) ? CmB.dot(H) : -CmB.dot(H);
            let r = Bb.length() / 2;
            let c = 2*r;
            let phi = Math.acos(da/L);
            let theta = Math.PI;
            return [r,theta,da,c,phi,H];
        } else {
            let H = new THREE.Vector3(0,0,0);
            H.crossVectors(Bb,Cb);
            H.normalize();
            // This is "rejection of CmB along H".
            // Negative indicates a negative z travel
            let da = CmB.dot(H);
            let phi = Math.acos(da/L);
            let Bax = Math.sin(phi) * da /2;

            // This be undefinied if phi = PI/2.
            // if phi = PI/2, then B on the axis is
            // the intersection of the Bb and Cb,
            // which will also have x and z = 0.
            var Ba;
            if (near(phi,Math.PI/2,1e-4)) {
                let psi = Math.atan2(Bb.y,Bb.z);
                Ba = new THREE.Vector3(0,
                                       Math.tan(psi)* L/2,
                                       0);
            } else {
                Ba = new THREE.Vector3(Bax,
                                       Bb.y * ( Bax / Bb.x),
                                       -Math.cos(phi) * da /2);
            }
            var Ba_m_B = Ba.clone();
            Ba_m_B.sub(B);
            let r = Ba_m_B.length();
            let c = ChordFromLDaxis(L,da);
            let theta = RotationFromRadiusChord(r,c);
            console.assert(!isNaN(theta));
            return [r,theta,da,c,phi,H];
        }
    }
}

// Add the novel prism to the old prism by placing the B face
// aginst the C face of the old with a twist of tau and return.
// note that old is a prism instance, not a mesh instance,
// and it has a link to the abstract prism inside it.
// NOTE: with faces reversed, this still produces
// positive z.
function adjoinPrism(old,tau,joinToC) {
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

    if (joinToC) {
        applyQuaternionToPrism(nu,Q_to_Nb);
        applyQuaternionToPrism(nu,Q_to_Nc);
    } else {
        applyQuaternionToPrism(nu,Q_to_Nc);
        applyQuaternionToPrism(nu,Q_to_Nb);
    }

    var rt = new THREE.Matrix4();
    rt.makeRotationAxis(joinToC ? nu.p.Nc : nu.p.Nb,tau);

    applyMatrix4ToPrism(nu,rt);

    applyMatrix4ToPrism(nu,p_trans_r);

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
    return nu;
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

    // for convenience, should I make the instance vectors
    // part of this object? Or should I add the matrices,
    // effective make it a "screw" transformation
    return {p: absp,
            b: b,
            c: c,
            nb: absp.Nb,
            nc: absp.Nc,
            tb: TB, lb : LB, rb : RB,
            tc : TC, lc : LC, rc : RC};

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



function AfromLtauNbNc(L,tau,NBu,NCu) {
    var abs0 = new AbstractPrism(L,NBu,NCu);
    if (!abs0) {
        debugger;
    }
    var p_i = CreatePrism(abs0,PRISM_FACE_RATIO_LENGTH);

    var rt = new THREE.Matrix4();
    var p_b = adjoinPrism(p_i,tau,false);
    var p_c = adjoinPrism(p_i,tau,true);

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
            console.log(psi * 180 / Math.PI);
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

    let psi = compute_angle_midpoint(fb,fc);

    rt.makeRotationAxis(new THREE.Vector3(0,0,1),psi);

    applyMatrix4ToPrism(p_i,rt);
    p_i.p.Nb.applyMatrix4(rt);
    p_i.p.Nc.applyMatrix4(rt);

    var p_b = adjoinPrism(p_i,tau,false);
    var p_c = adjoinPrism(p_i,tau,true);
    // The x values of these two should be the oppoosite

    console.assert(near(p_b.b.x,-p_c.c.x,1e-4));
    console.assert(near(p_b.b.y,p_c.c.y,1e-4));

    if (!near(p_b.b.x,-p_c.c.x,1e-4)) {
        debugger;
    }
    if (!near(p_b.b.y,p_c.c.y,1e-4)) {
        debugger;
    }

    return [p_b.b,psi];
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


                    var p_b = adjoinPrism(p_i,tau,false);
                    var p_c = adjoinPrism(p_i,tau,true);
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


// This deosn't work---Unified need to be separated
function testRegularTetsUnified() {
    let L0 = 1;
    let tetAngle = Math.atan(Math.sqrt(2)/2);
    /* Rotation by an aribitrary amount should give us the same number. */
    testRot = Math.PI / 5;
    var rt = new THREE.Matrix4();
    rt.makeRotationAxis(new THREE.Vector3(0,0,1),testRot);
    let plainB = new THREE.Vector3(0,-Math.sin(tetAngle),-Math.cos(tetAngle));
    let plainC = new THREE.Vector3(0,-Math.sin(tetAngle),Math.cos(tetAngle));
    let tau = 2 * Math.PI /3;
    var NB1 = plainB.clone();
    NB1.applyMatrix4(rt);
    var NC1 = plainC.clone();
    NC1.applyMatrix4(rt);
    var A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
    let B = new THREE.Vector3(0,0,-L0/2);
    let D = new THREE.Vector3(-A.x,A.y,-A.z);
//    let res = UnifiedComp(L0,B,D);
    let res = KahnAxis(L0,D);
    let r = res[0];
    let theta = res[1];
    let da = res[2];
    let c = res[3];
    let phi = res[4];
    let theta_exp = Math.acos(-2/3);
    console.assert(near(theta,theta_exp,0.00001));
}

function testRegularTetsKahnAxis()
{
    let L0 = 1/3;
    let tetAngle = Math.atan(Math.sqrt(2)/2);
    let NB1 = new THREE.Vector3(0,-Math.sin(tetAngle),-Math.cos(tetAngle));
    let NC1 = new THREE.Vector3(0,-Math.sin(tetAngle),Math.cos(tetAngle));
    let tau = 2 * Math.PI /3;
    let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
    let B = new THREE.Vector3(0,0,-L0/2);
    let D = new THREE.Vector3(-A.x,A.y,-A.z);
    let res = KahnAxis(L0,D);
    let r = res[0];
    let theta = res[1];
    let da = res[2];
    let c = res[3];
    let phi = res[4];
    let theta_exp = Math.acos(-2/3);
    console.assert(near(theta,theta_exp,0.00001));
}

function testKahnAxisYTorus()
{
    let L0 = 2;
    let angle = Math.PI/7;
    let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
    let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
    let tau = 0;
    let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
    let B = new THREE.Vector3(0,0,-L0/2);
    let D = new THREE.Vector3(-A.x,A.y,-A.z);
    let res = KahnAxis(L0,D);
    let r = res[0];
    let theta = res[1];
    let da = res[2];
    let c = res[3];
    let phi = res[4];
    console.assert((theta * 180/ Math.PI) != 180);
    console.assert(near(da,0,0.01));
}
// That that tau = 180 is not degenerate..
function testKahnAxisTau180()
{
    let L0 = 2;
    //    let angle = Math.PI/7;
    let angle = 0;
    let NB1 = new THREE.Vector3(0,-Math.sin(angle),-Math.cos(angle));
    // We add in a little here so we are not a torus
    let NC1 = new THREE.Vector3(0,-Math.sin(angle),Math.cos(angle));
    let tau = Math.PI;
    let A = AfromLtauNbNc(L0,tau,NB1,NC1)[0];
    let B = new THREE.Vector3(0,0,-L0/2);
    let D = new THREE.Vector3(-A.x,A.y,-A.z);
    let res = KahnAxis(L0,D);
    let r = res[0];
    let theta = res[1];
    let da = res[2];
    let c = res[3];
    let phi = res[4];
//    console.assert((theta * 180/ Math.PI) != 180);
    console.assert(!near(da,0,0.1));
    console.assert(!isNaN(theta));
}


function runUnitTests() {
    testAfromLtauNbNnKeepsYPure();
    testAfromLtauNbNnContinuityOfTau();
    testRegularTetsKahnAxis();
    testKahnAxisYTorus();
    testKahnAxisTau180();
    testAfromLfailureCase1();
    testAfromLtauMultiple();
}
