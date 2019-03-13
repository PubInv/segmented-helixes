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

// return r,theta,d,c,phi    
function UnifiedComp(L,rho,omega) {
    const B = new THREE.Vector3(0,0,-L/2);
    const C = new THREE.Vector3(0,0,L/2);            
    const q = QFromRhoOmega(L,rho,omega);
    const x = q * Math.tan(omega);
    const y = q * Math.tan(rho);
    const u = q - L;
    const L2q = (L/2) + q; 
    const A = new THREE.Vector3(-x,y,-L2q);
    const D = new THREE.Vector3(x,y,L2q);
    console.log("ABCD",A,B,C,D);
    const Bb = new THREE.Vector3(x,-y,u);
    const Cb = new THREE.Vector3(-x,-y,-u);
    const BbXCb = new THREE.Vector3().crossVectors(Bb,Cb);
    const n2 = new THREE.Vector3().crossVectors(Cb,BbXCb);
    const BbN2inner = Bb.dot(n2);
    const n2z = n2.z;            
    const Bn2 = n2z / (BbN2inner);
    const LBn2 = L * Bn2;
    const C1 = new THREE.Vector3(Bb.x,Bb.y,Bb.z).multiplyScalar(LBn2);
    C1.add(B);
    // now, if we haven't messed up, phi can be computed from C1...
    const phi = Math.atan2(Math.abs(C1.x),Math.abs(C1.z));
    const di = 2 * Math.sqrt(C1.x ** 2 + C1.z ** 2);
    // Have to be careful, multiplyScalar changes its value.
    const r = Bb.multiplyScalar(LBn2).length();
    const c = ChordFromLD(L,di);
    const theta2 = ThetaFromRC(r,c);
    return [r,theta2,di,c,phi];
}
function ChordFromLDaxis(L,Da) {
        return Math.sqrt((L**2) - (Da**2));
}
// Note this has a "special" case...
function RotationFromRadiusChord(R,C) {
    let v = C / (2 * R);
    if (v == 1) {
        return Math.PI;
    } else {
        return 2 * Math.asin(C / (2 * R));
    }
}

// D is a the point, L is the 
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
    // Conjecture--the lengths are always the same?
    console.assert(near(Bb.length(),Cb.length(),1e-3));
    // Now if these are zero what do we do?
    if (Bb.length() == 0) {
        console.log("STRAIGHT LINE!");
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
    let H = new THREE.Vector3(0,0,0);
    H.crossVectors(Bb,Cb);
    H.normalize();
    let CmB = C.clone();
    CmB.sub(B);
    console.log("CmB,H :",CmB,H);
    let da = CmB.dot(H);
    let phi = Math.acos(da/L);
    let Bax = Math.sin(phi) * da /2;
    // This be undefinied if phi = PI/2.
    // if phi = PI/2, then B on the axis is
    // the intersection of the Bb and Cb,
    // which will also have x and z = 0.
    var Ba;
    if (near(phi,Math.PI/2,1e-4)) {
        console.log("XXXX",phi);
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
        return [r,theta,da,c,phi,H];
    }
}

// Add the novel prism to the old prism by placing the B face
// aginst the C face of the old with a twist of tau and return.
// note that old is a prism instance, not a mesh instance,
// and it has a link to the abstract prism inside it.
// NOTE: This is not currently applying tau!
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
        console.log("CONDITIONED ANGLE");
        return condition_angle(angle + (2 * Math.PI));
    } else if (angle > (2 * Math.PI)) {
        console.log("CONDITIONED ANGLE");        
        return condition_angle(angle - (2 * Math.PI));
    } else {
        return angle;
    }
}

function isNumeric(n) {
    return !isNaN(parseFloat(n)) && isFinite(n);
}


// I am starting to think this is conceptually flawed.
// We don't actually care about balancing the normals, we care about
// balancing the joint vectors. I have been confused about this.
// This needs to take tau, and either compute or balance the face vectors!
// That will mean this has to be a function of tau.
function ComputeBalancingRotation(Nb,Nc) {
    var b = Nb.clone();
    b.normalize();
    var c = Nc.clone();
    c.normalize();
    // Note: if b and c are collinear, we want to
    // use the angle of the one with longest projection in XY.
    // If they are opposite, we have a zig-zig,
    // if they point the same way we have a 
    
    // Note: Mathemtaica uses a different order
    var pb = new THREE.Vector2(b.x,b.y);
    var pc = new THREE.Vector2(c.x,c.y);    
    var ba = pb.angle();
    var ca = pc.angle();
    var theta;
    
    if (near(ba,ca,1e-3)) {
        theta = -ba;
    } else if (near(ba-ca, Math.PI, 1e-3) || near(ba-ca, -Math.PI, 1e-3)) {
        if (pb.length() == pc.length()) {
            console.log("WARING---collinear opposing normals!");
            theta = -ba;
        } else if (pb.length() < pc.length()) {
            theta = -ca;
        } else {
            theta = -ba;
        }
    } else {
    // This should perhaps be a subtraction between ba and ca!
    // the -90 degrees is to make the average point straight down.
    // Do we want the average or the difference here?
        // theta = -(ba+ca)/2;
        // Is it better to find the midpoint between
        // the normalized vectors, find it's angle,
        // and rotate by the amount that puts it straight down?
        let ban = pb.clone().normalize();
        let can = pc.clone().normalize();
        let m = ban.clone().add(can);
        theta = -m.angle();
    }
    // now, instead of pointing at the x-axis, we point at the negative y axis.
    theta += -Math.PI/2;
    theta = condition_angle(theta);
    var rt = new THREE.Matrix4();
    rt.makeRotationAxis(new THREE.Vector3(0,0,1),theta);
    b.applyMatrix4(rt);
    c.applyMatrix4(rt);
    console.log("BALANCE",theta);
    return [theta,b,c];
}

function testComputeBalancingRotation0() {
    let L = 1;
    // This should be 45 degrees
    let Nb = new THREE.Vector3(1/2,1/2,-1);
    // This should be 0
    let Nc = new THREE.Vector3(0,0,1);
    // our rotation should be -22.5 + 90
    var res = ComputeBalancingRotation(Nb,Nc);
    var theta = res[0];
    var Br = res[1];
    var Cr = res[2];
    var Bp = new THREE.Vector2(Br.x,Br.y);
    var Cp = new THREE.Vector2(Cr.x,Cr.y);
    Bp.normalize();
    Cp.normalize();
    //        console.log(Bp,Cp);
    console.assert(Bp.length() == 0 || Cp.length() == 0 || near(Bp.x,-Cp.x,1e-4));
    console.assert(Bp.y <= 0);
    console.assert(Cp.y <= 0);
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



// This is a serious question here...
// if we are really attempting to compute "A"
// in space, then we cannot throw away the
// rotational angle from Compute Balancing Rotation....
// We need to rotate back by that amount, I guess.
// We have code that renders the adjoining prisms...
// The ultimate test is really that this must
// return the same thing as that code! We know
// the Adjoining Prism method is right, becase
// we an see it lining up in its rendering.
// So tomorrow: Write a test to make sure this
// matches the adjoinging prism. Then
// figure out what is wrong with this code.
function AfromLtauNbNc(L,tau,NBu,NCu) {
    
    // let res = ComputeBalancingRotation(NBu,NCu);
    // let theta = res[0];
    // let Nb = res[1];
    // let Nc = res[2];


    // I now believe this approach is completely wrong.
    // We actually need to balance the vectors of the adjoined prisms.
    // Although it may be too much, I need to somehow start with the
    // math from adjoining and simplify it!
    // This is the alternative means of computation.
    // let Z = new THREE.Vector3(0,0,1);
    // let delta = Z.angleTo(Nc);
    // let v0 = V0fromLNB(L,Nb,Nc,delta);
    // console.assert(near(L,v0.length(),1e-5));
    // let Ad = ADirFromParam(L,v0,tau,Nb);
    // let B = new THREE.Vector3(0,0,-L/2);
    // var result = Ad.clone();
    // result.add(B);



    // This is code based on building the physical prism.
    // in a sense it is more accurate. The fact
    // that these values don't match is a serious problem.

    // NOTE: This is all pretty horrible.
    // There are several problems:
    // Why am I not computing the angle correctly on the first
    // try?
    // Why does the code above, which should compute it with
    // having to do all this physical simulation, not work?
    // Why I am very slighly off (the helix doesn't quite
    // pierce the joint in a centered way)?
    // Why do some values clearly get far out of whack?

    var abs0 = new AbstractPrism(L,NBu,NCu);
    if (!abs0) {
        debugger;
    }
    var p_i = CreatePrism(abs0,PRISM_FACE_RATIO_LENGTH);

    var rt = new THREE.Matrix4();
//    rt.makeRotationAxis(new THREE.Vector3(0,0,1),theta);
    
//    applyMatrix4ToPrism(p_i,rt);
//    p_i.p.Nb.applyMatrix4(rt);
//    p_i.p.Nc.applyMatrix4(rt);
    
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
            console.log(fba * 180 / Math.PI);
            console.log(fca * 180 / Math.PI);                        
            let psi = -( (fba > fca) ? fba - Math.PI : fca - Math.PI );
            console.log(psi * 180 / Math.PI);
            return psi;
        }
        // The ".angle()" function in THREE measures against the X-axis
        // I assume this is "reverse on the clock" operation!
        let psi = -(fbc.angle() - ( 3 * Math.PI / 2));
        return psi;
    }

    // This is likely the problem....
    let psi = compute_angle_midpoint(fb,fc);

    console.log("PSI = ",psi * 180 / Math.PI);

    // The fact that psi is not Beta somehow means my balancing computation
    // above is not correct (by definition.)
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

    //   return result;
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

function testComputeBalancingRotation3() {
    var B = new THREE.Vector3(1,1/2,1);
    var C = new THREE.Vector3(-1/2,-1,-1);
    var res = ComputeBalancingRotation(B,C);
    
    var theta = res[0];
    var Br = res[1];
    var Cr = res[2];
    var Bp = new THREE.Vector2(Br.x,Br.y);
    var Cp = new THREE.Vector2(Cr.x,Cr.y);
    Bp.normalize();
    Cp.normalize();
    //        console.log(Bp,Cp);
    console.assert(Bp.length() == 0 || Cp.length() == 0 || near(Bp.x,-Cp.x,1e-4));
    console.assert(Bp.y <= 1e-6);
    console.assert(Cp.y <= 1e-6);
}

// This needs to work. We also must test other
// examples of colinearity between the projection
// of Nb and Nc
function testComputeBalancingRotation1() {
    let Nb = new THREE.Vector3(0, -0.25254006068835666,-0.967586439418991);
    let Nc = new THREE.Vector3(0,0.31337745420390184,0.9496286491027328);
    var res = ComputeBalancingRotation(Nb,Nc);

    var theta = res[0];
    var Br = res[1];
    var Cr = res[2];
    var Bp = new THREE.Vector2(Br.x,Br.y);
    var Cp = new THREE.Vector2(Cr.x,Cr.y);
    Bp.normalize();
    Cp.normalize();
    //        console.log(Bp,Cp);
    console.assert(Bp.length() == 0 || Cp.length() == 0 || near(Bp.x,-Cp.x,1e-4));
    console.assert((Bp.y + Cp.y) <= 1e-6);
}

function testComputeBalancingRotation2() {
    
    let Nb = new THREE.Vector3(-1/3,-1 ,-1);
    let Nc = new THREE.Vector3(-1/6,-1,1);
    Nb.normalize();
    Nc.normalize();
    for(var i = 0; i < 10; i++) {
        // we'll test through 20 degrees
        let phi = 2 * i * Math.PI / 180;
        var nb = Nb.clone();
        var nc = Nc.clone();
        var rt = new THREE.Matrix4();
        rt.makeRotationAxis(new THREE.Vector3(0,0,1),phi);
        nb.applyMatrix4(rt);
        nc.applyMatrix4(rt);
        var res = ComputeBalancingRotation(nb,nc);
        var theta = res[0];
        var Br = res[1];
        var Cr = res[2];
//        console.log(theta * 180 / Math.PI, Br, Cr);
        var Bp = new THREE.Vector2(Br.x,Br.y);
        var Cp = new THREE.Vector2(Cr.x,Cr.y);
        Bp.normalize();
        Cp.normalize();
//        console.log(Bp,Cp);
        console.assert(near(Bp.x,-Cp.x,1e-4));
        console.assert(Bp.y <= 1e-6);
        console.assert(Cp.y <= 1e-6);
    }
}

function test_UnifiedComp() {
    const res = UnifiedComp(1,Math.PI/10,Math.PI/30);
    console.log(res);
    console.log(res[4] * 180 / Math.PI);
    console.assert(near(res[0],0.688145,1e-4));
    console.assert(near(res[1],0.704444,1e-4));        
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
    console.assert((theta * 180/ Math.PI) != 180);
    console.assert(!near(da,0,0.1));
    console.assert(!isNaN(theta));
}


function runUnitTests() {
    test_UnifiedComp();
    testComputeBalancingRotation0();
    testComputeBalancingRotation1();
    testComputeBalancingRotation2();
    testComputeBalancingRotation3();    
    testAfromLtauNbNnKeepsYPure();
    testAfromLtauNbNnContinuityOfTau();    
    testRegularTetsKahnAxis();
    testKahnAxisYTorus();
    testKahnAxisTau180();
    testAfromLfailureCase1();
    testAfromLtauMultiple();
}
