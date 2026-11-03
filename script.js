
(function(){
"use strict";
/* ---------- fluids & state ---------- */
const FLUIDS={air:{nu:1.516e-5,rho:1.204},water:{nu:1.004e-6,rho:998},glycerin:{nu:1.12e-3,rho:1260}};
const DEF={fluid:"air",U:5,L:1,Recr:5e5,xs:0.6};
const S=Object.assign({tracers:true,profiles:true,env:true,cubic:true,paused:false},DEF);
const nu=()=>FLUIDS[S.fluid].nu, rho=()=>FLUIDS[S.fluid].rho;

/* ---------- Blasius solution by RK4 ---------- */
const BH=0.005, BN=2400; // eta 0..12
const BF=(function(){
  const t=new Float64Array(BN+1); let y=[0,0,0.332057];
  const d=v=>[v[1],v[2],-0.5*v[0]*v[2]];
  for(let i=1;i<=BN;i++){
    const k1=d(y), k2=d(y.map((v,j)=>v+BH/2*k1[j])), k3=d(y.map((v,j)=>v+BH/2*k2[j])), k4=d(y.map((v,j)=>v+BH*k3[j]));
    y=y.map((v,j)=>v+BH/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j])); t[i]=y[1];
  }
  return t;
})();
function blasius(eta){ if(eta<=0) return 0; if(eta>=BN*BH) return 1; const f=eta/BH, i=Math.floor(f), r=f-i; return Math.min(1,BF[i]*(1-r)+BF[i+1]*r); }

/* ---------- physics: ONE shared module, used by the 3D simulation, the readout and the equations page ---------- */
// P = {U [m/s], L [m], nu [m²/s], rho [kg/m³], Recr [-]}. Lengths always in metres.
const PHYS={
  Re:(P,x)=>P.U*x/P.nu,
  xcr:P=>P.Recr*P.nu/P.U,
  isTurb:(P,x)=>PHYS.Re(P,x)>=P.Recr,
  lam:{ // Blasius, laminar
    d:(P,x)=>4.91*x/Math.sqrt(PHYS.Re(P,x)),
    ds:(P,x)=>1.7208*x/Math.sqrt(PHYS.Re(P,x)),
    th:(P,x)=>0.664*x/Math.sqrt(PHYS.Re(P,x)),
    cf:(P,x)=>0.664/Math.sqrt(PHYS.Re(P,x)),
    CD:P=>1.328/Math.sqrt(PHYS.Re(P,P.L)),
    eta:(P,x,y)=>y*Math.sqrt(P.U/(P.nu*x)),
    u:(P,x,y)=>blasius(PHYS.lam.eta(P,x,y))
  },
  tur:{ // 1/7-power law, turbulent from the leading edge
    d:(P,x)=>0.37*x/Math.pow(PHYS.Re(P,x),0.2),
    ds:(P,x)=>PHYS.tur.d(P,x)/8,
    th:(P,x)=>7*PHYS.tur.d(P,x)/72,
    cf:(P,x)=>0.059/Math.pow(PHYS.Re(P,x),0.2),
    CDfull:P=>0.074/Math.pow(PHYS.Re(P,P.L),0.2),
    A:P=>P.Recr*(0.074/Math.pow(P.Recr,0.2)-1.328/Math.sqrt(P.Recr)),
    CD:P=>PHYS.tur.CDfull(P)-PHYS.tur.A(P)/PHYS.Re(P,P.L),
    u:(P,x,y)=>{ const d=PHYS.tur.d(P,x); return y>=d?1:Math.pow(y/d,1/7); }
  },
  cub:{ // cubic polynomial (laminar approximation)
    d:(P,x)=>4.64*x/Math.sqrt(PHYS.Re(P,x)),
    u:(P,x,y)=>{ const r=y/PHYS.cub.d(P,x); return r>=1?1:1.5*r-0.5*r*r*r; }
  },
  // applied values: regime chosen by local Re_x against Re_cr, as the simulation does
  delta:(P,x)=>x<=0?0:(PHYS.isTurb(P,x)?PHYS.tur.d(P,x):PHYS.lam.d(P,x)),
  ds:(P,x)=>x<=0?0:(PHYS.isTurb(P,x)?PHYS.tur.ds(P,x):PHYS.lam.ds(P,x)),
  th:(P,x)=>x<=0?0:(PHYS.isTurb(P,x)?PHYS.tur.th(P,x):PHYS.lam.th(P,x)),
  cf:(P,x)=>PHYS.isTurb(P,x)?PHYS.tur.cf(P,x):PHYS.lam.cf(P,x),
  CD:P=>PHYS.Re(P,P.L)<=P.Recr?PHYS.lam.CD(P):PHYS.tur.CD(P),
  u:(P,x,y)=>{ if(x<=0) return 1; if(y<=0) return 0; return PHYS.isTurb(P,x)?PHYS.tur.u(P,x,y):PHYS.lam.u(P,x,y); }
};
// simulation binding: SP mirrors the slider state; refreshed at the start of update()
const SP={U:S.U,L:S.L,nu:nu(),rho:rho(),Recr:S.Recr};
function syncSP(){ SP.U=S.U; SP.L=S.L; SP.nu=nu(); SP.rho=rho(); SP.Recr=S.Recr; }
const Rex=x=>PHYS.Re(SP,x), xTr=()=>PHYS.xcr(SP);
const delta=x=>PHYS.delta(SP,x), dispThick=x=>PHYS.ds(SP,x), momThick=x=>PHYS.th(SP,x), cfx=x=>PHYS.cf(SP,x), CD=()=>PHYS.CD(SP);
const uRatio=(x,y)=>PHYS.u(SP,x,y), cubic=(y,x)=>PHYS.cub.u(SP,x,y);

/* ---------- formatting ---------- */
const SUP={"-":"⁻","0":"⁰","1":"¹","2":"²","3":"³","4":"⁴","5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹"};
function sci(v,d=2){ if(!isFinite(v)||v===0) return "0"; const e=Math.floor(Math.log10(Math.abs(v))); if(e>=-1&&e<4) return (+v.toPrecision(3)).toString(); return (v/Math.pow(10,e)).toFixed(d)+"×10"+String(e).split("").map(c=>SUP[c]).join(""); }
function len(m){ if(m>=1) return m.toFixed(3)+" m"; if(m>=1e-3) return (m*1e3).toPrecision(3)+" mm"; return (m*1e6).toPrecision(3)+" µm"; }

/* ---------- sliders ---------- */
const $=id=>document.getElementById(id);
const UMIN=0.05,UMAX=30, RCMIN=1e5,RCMAX=3e6;
const logTo=(t,a,b)=>a*Math.pow(b/a,t/1000), logFrom=(v,a,b)=>Math.round(1000*Math.log(v/a)/Math.log(b/a));

function syncInputs(){
  $("sU").value=logFrom(S.U,UMIN,UMAX); $("sL").value=S.L; $("sRc").value=logFrom(S.Recr,RCMIN,RCMAX); $("sX").value=S.xs;
  document.querySelectorAll("#fluidSeg button").forEach(b=>b.setAttribute("aria-pressed",b.dataset.f===S.fluid));
  $("tTr").checked=S.tracers;$("tPr").checked=S.profiles;$("tEn").checked=S.env;$("tCu").checked=S.cubic;
}
$("sU").addEventListener("input",e=>{S.U=logTo(+e.target.value,UMIN,UMAX);update();});
$("sL").addEventListener("input",e=>{S.L=+e.target.value;update();});
$("sRc").addEventListener("input",e=>{S.Recr=logTo(+e.target.value,RCMIN,RCMAX);update();});
$("sX").addEventListener("input",e=>{S.xs=+e.target.value;update();});
document.querySelectorAll("#fluidSeg button").forEach(b=>b.addEventListener("click",()=>{S.fluid=b.dataset.f;syncInputs();update();}));
[["tTr","tracers"],["tPr","profiles"],["tEn","env"],["tCu","cubic"]].forEach(([id,k])=>$(id).addEventListener("change",e=>{S[k]=e.target.checked;update();}));
$("bPlay").addEventListener("click",()=>{S.paused=!S.paused;$("bPlay").textContent=S.paused?"Play":"Pause";});
$("bReset").addEventListener("click",()=>{Object.assign(S,DEF);syncInputs();update();resetView();});
const PRESETS={lam:{fluid:"air",U:2,L:1,Recr:5e5,xs:0.6},tr:{fluid:"air",U:15,L:1,Recr:5e5,xs:0.8},water:{fluid:"water",U:0.3,L:1,Recr:5e5,xs:0.6},low:{fluid:"glycerin",U:1,L:0.5,Recr:5e5,xs:0.5}};
document.querySelectorAll(".preset").forEach(b=>b.addEventListener("click",()=>{Object.assign(S,PRESETS[b.dataset.p]);syncInputs();update();}));

/* ---------- three.js scene ---------- */
const host=$("three");
const renderer=new THREE.WebGLRenderer({antialias:true});
renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));
renderer.setClearColor(0x0a0e11,1);
host.appendChild(renderer.domElement);
const scene=new THREE.Scene();
scene.fog=new THREE.Fog(0x0a0e11,26,48);
const camera=new THREE.PerspectiveCamera(42,1,0.1,200);
const controls=new THREE.OrbitControls(camera,renderer.domElement);
controls.enableDamping=true; controls.minDistance=4; controls.maxDistance=34; controls.maxPolarAngle=Math.PI/2-0.03;
function resetView(){camera.position.set(3.2,6.2,13.5);controls.target.set(5,1.0,0);controls.update();}
resetView();
scene.add(new THREE.AmbientLight(0xffffff,0.55));
const dl=new THREE.DirectionalLight(0xffffff,0.7); dl.position.set(4,10,6); scene.add(dl);

const PX=10, PZ=4; // plate display size
const plate=new THREE.Mesh(new THREE.BoxGeometry(PX,0.08,PZ),new THREE.MeshStandardMaterial({color:0x5b6872,metalness:0.6,roughness:0.45}));
plate.position.set(PX/2,-0.04,0); scene.add(plate);
// plate grid (L/10)
(function(){
  const p=[]; for(let i=0;i<=10;i++){p.push(i,0.002,-PZ/2,i,0.002,PZ/2);} p.push(0,0.002,-PZ/2,PX,0.002,-PZ/2,0,0.002,PZ/2,PX,0.002,PZ/2);
  const g=new THREE.BufferGeometry(); g.setAttribute("position",new THREE.Float32BufferAttribute(p,3));
  scene.add(new THREE.LineSegments(g,new THREE.LineBasicMaterial({color:0x34424c})));
})();
// leading edge highlight
(function(){const g=new THREE.BufferGeometry();g.setAttribute("position",new THREE.Float32BufferAttribute([0,0.004,-PZ/2,0,0.004,PZ/2],3));scene.add(new THREE.Line(g,new THREE.LineBasicMaterial({color:0xdbe3e8})));})();
const pickPlane=new THREE.Mesh(new THREE.PlaneGeometry(PX,PZ),new THREE.MeshBasicMaterial({visible:false}));
pickPlane.rotation.x=-Math.PI/2; pickPlane.position.set(PX/2,0.01,0); scene.add(pickPlane);

// display mapping
let SX=PX/S.L, SY=1, dMax=1;
const YTOP=2.4;
function remap(){ SX=PX/S.L; dMax=0; for(let i=1;i<=400;i++){dMax=Math.max(dMax,delta(S.L*i/400));} SY=YTOP/dMax; }

/* colour map u/U -> rgb */
const CM=[[0,[.165,.184,.561]],[.25,[.184,.498,.902]],[.5,[.122,.827,.608]],[.75,[.941,.761,.173]],[1,[.886,.255,.165]]];
function cmap(t,out){ t=Math.max(0,Math.min(1,t)); for(let i=1;i<CM.length;i++){ if(t<=CM[i][0]){const a=CM[i-1],b=CM[i],r=(t-a[0])/(b[0]-a[0]); out[0]=a[1][0]+(b[1][0]-a[1][0])*r; out[1]=a[1][1]+(b[1][1]-a[1][1])*r; out[2]=a[1][2]+(b[1][2]-a[1][2])*r; return out;} } out[0]=CM[4][1][0];out[1]=CM[4][1][1];out[2]=CM[4][1][2]; return out; }

/* tracers */
const N=2600, XMIN=-2.2, XMAX=PX+0.3;
const pos=new Float32Array(N*3), col=new Float32Array(N*3), py0=new Float32Array(N), ph=new Float32Array(N);
function spawnY(){ const r=Math.random(); return 0.015+Math.pow(r,1.6)*3.5; } // denser near wall
for(let i=0;i<N;i++){pos[i*3]=XMIN+Math.random()*(XMAX-XMIN); py0[i]=spawnY(); pos[i*3+2]=(Math.random()-0.5)*PZ*0.96; ph[i]=Math.random()*Math.PI*2;}
const tg=new THREE.BufferGeometry(); tg.setAttribute("position",new THREE.BufferAttribute(pos,3)); tg.setAttribute("color",new THREE.BufferAttribute(col,3));
const dot=(function(){const c=document.createElement("canvas");c.width=c.height=32;const x=c.getContext("2d");const g=x.createRadialGradient(16,16,0,16,16,16);g.addColorStop(0,"#fff");g.addColorStop(.55,"#fff");g.addColorStop(1,"rgba(255,255,255,0)");x.fillStyle=g;x.fillRect(0,0,32,32);return new THREE.CanvasTexture(c);})();
const tracers=new THREE.Points(tg,new THREE.PointsMaterial({size:0.085,vertexColors:true,map:dot,transparent:true,alphaTest:0.2,depthWrite:false}));
scene.add(tracers);

/* envelope */
const EN=200;
const envGeo=new THREE.BufferGeometry(), envPos=new Float32Array((EN+1)*4*3); const idx=[];
for(let i=0;i<EN;i++){ const a=i*4,b=(i+1)*4; idx.push(a,b,a+1, a+1,b,b+1); idx.push(a+2,b+2,a+3, a+3,b+2,b+3); }
envGeo.setAttribute("position",new THREE.BufferAttribute(envPos,3)); envGeo.setIndex(idx);
const envMesh=new THREE.Mesh(envGeo,new THREE.MeshBasicMaterial({color:0x58cfe0,transparent:true,opacity:0.12,side:THREE.DoubleSide,depthWrite:false}));
scene.add(envMesh);
const edgeGeo=new THREE.BufferGeometry(), edgePos=new Float32Array((EN+1)*3);
edgeGeo.setAttribute("position",new THREE.BufferAttribute(edgePos,3));
const edge=new THREE.Line(edgeGeo,new THREE.LineBasicMaterial({color:0x58cfe0})); scene.add(edge);
const edgeGeo2=new THREE.BufferGeometry(), edgePos2=new Float32Array((EN+1)*3);
edgeGeo2.setAttribute("position",new THREE.BufferAttribute(edgePos2,3));
const edge2=new THREE.Line(edgeGeo2,new THREE.LineBasicMaterial({color:0x2c6c76})); scene.add(edge2);

function buildEnvelope(){
  for(let i=0;i<=EN;i++){
    const xd=PX*i/EN, xp=S.L*i/EN, yd=Math.min(delta(xp)*SY,YTOP*1.6), o=i*12;
    // top sheet (z=-2,z=+2), front curtain (z=+2: y=0, y=yd)
    envPos.set([xd,yd,-PZ/2, xd,yd,PZ/2, xd,0.003,PZ/2, xd,yd,PZ/2],o);
    edgePos.set([xd,yd,PZ/2],i*3); edgePos2.set([xd,yd,-PZ/2],i*3);
  }
  envGeo.attributes.position.needsUpdate=true; envGeo.computeBoundingSphere();
  edgeGeo.attributes.position.needsUpdate=true; edgeGeo2.attributes.position.needsUpdate=true;
  edgeGeo.computeBoundingSphere(); edgeGeo2.computeBoundingSphere();
}

/* profiles */
let profLines=null, selLines=null;
const STATIONS=[0.12,0.32,0.52,0.72,0.92], K=14, ARL=1.35;
function profileSegs(xFrac,seg,colArr,c){
  const xp=S.L*xFrac, xd=PX*xFrac, d=delta(xp), top=1.25*d; let prev=null;
  const push=(a,b)=>{seg.push(...a,...b);colArr.push(...c,...c);};
  push([xd,0.003,0],[xd,top*SY,0]);
  for(let j=0;j<=K;j++){
    const yp=top*j/K, yd=yp*SY, u=uRatio(xp,yp), tip=[xd+u*ARL,yd,0];
    if(u>0.02){ push([xd,yd,0],tip); const h=Math.min(0.09,u*ARL*0.3); push(tip,[tip[0]-h,yd+h*0.5,0]); push(tip,[tip[0]-h,yd-h*0.5,0]); }
    if(prev) push(prev,tip); prev=tip;
  }
  // fine curve
  let pv=null; for(let j=0;j<=60;j++){const yp=top*j/60,tip=[xd+uRatio(xp,yp)*ARL,yp*SY,0]; if(pv) push(pv,tip); pv=tip;}
}
function lineObj(seg,colArr){ const g=new THREE.BufferGeometry(); g.setAttribute("position",new THREE.Float32BufferAttribute(seg,3)); g.setAttribute("color",new THREE.Float32BufferAttribute(colArr,3)); return new THREE.LineSegments(g,new THREE.LineBasicMaterial({vertexColors:true})); }
function buildProfiles(){
  [profLines,selLines].forEach(o=>{ if(o){scene.remove(o);o.geometry.dispose();o.material.dispose();} });
  let seg=[],c=[]; STATIONS.forEach(f=>{ if(Math.abs(f-S.xs)>0.04) profileSegs(f,seg,c,[0.86,0.89,0.91]); });
  profLines=lineObj(seg,c); profLines.visible=S.profiles; scene.add(profLines);
  seg=[];c=[]; profileSegs(S.xs,seg,c,[0.95,0.66,0.23]);
  const xd=PX*S.xs; seg.push(xd,0.006,-PZ/2,xd,0.006,PZ/2); c.push(0.95,0.66,0.23,0.95,0.66,0.23);
  selLines=lineObj(seg,c); scene.add(selLines);
}

/* transition marker */
const trGeo=new THREE.BufferGeometry(); trGeo.setAttribute("position",new THREE.Float32BufferAttribute([0,0.008,-PZ/2,0,0.008,PZ/2],3));
const trLine=new THREE.Line(trGeo,new THREE.LineBasicMaterial({color:0xef5a4c})); scene.add(trLine);

/* picking */
const ray=new THREE.Raycaster(), mv=new THREE.Vector2(); let down=null;
renderer.domElement.addEventListener("pointerdown",e=>{down=[e.clientX,e.clientY];});
renderer.domElement.addEventListener("pointerup",e=>{
  if(!down) return; const moved=Math.hypot(e.clientX-down[0],e.clientY-down[1]); down=null; if(moved>5) return;
  const r=renderer.domElement.getBoundingClientRect(); mv.set((e.clientX-r.left)/r.width*2-1,-(e.clientY-r.top)/r.height*2+1);
  ray.setFromCamera(mv,camera); const hit=ray.intersectObject(pickPlane)[0];
  if(hit){ S.xs=Math.max(0.01,Math.min(1,Math.round(hit.point.x/PX*100)/100)); $("sX").value=S.xs; update(); }
});

/* resize */
function resize(){ const w=host.clientWidth,h=host.clientHeight; if(!w||!h) return; renderer.setSize(w,h,false); camera.aspect=w/h; camera.updateProjectionMatrix(); drawGraphs(); }
new ResizeObserver(resize).observe(host);
const chartRO=new ResizeObserver(()=>drawGraphs()); chartRO.observe($("cProf")); chartRO.observe($("cGrow"));

/* ---------- collapsible sidebar ---------- */
const bench=document.querySelector(".bench"), side=$("sidebar"), sbBtn=$("bSide"), backdrop=$("sbBackdrop");
const mobile=window.matchMedia("(max-width:880px)");
function savedOpen(){ try{ return localStorage.getItem("bl-sidebar")!=="0"; }catch(e){ return true; } }
function setSidebar(open,save){
  bench.classList.toggle("sb-closed",!open);
  sbBtn.setAttribute("aria-expanded",String(open)); sbBtn.setAttribute("aria-label",open?"Hide controls":"Show controls");
  side.inert=!open; backdrop.hidden=!(open&&mobile.matches);
  if(save&&!mobile.matches){ try{ localStorage.setItem("bl-sidebar",open?"1":"0"); }catch(e){} }
}
sbBtn.addEventListener("click",()=>setSidebar(bench.classList.contains("sb-closed"),true));
$("bSideClose").addEventListener("click",()=>setSidebar(false,false));
backdrop.addEventListener("click",()=>setSidebar(false,false));
document.addEventListener("keydown",e=>{ if(e.key==="Escape"&&mobile.matches&&!bench.classList.contains("sb-closed")){ setSidebar(false,false); sbBtn.focus(); } });
(mobile.addEventListener?mobile.addEventListener.bind(mobile,"change"):mobile.addListener.bind(mobile))(()=>setSidebar(mobile.matches?false:savedOpen(),false));
bench.addEventListener("transitionend",e=>{ if(e.target===bench) resize(); });
setSidebar(mobile.matches?false:savedOpen(),false);

/* HUD tags */
const v3=new THREE.Vector3();
function placeTag(el,x,y,z){ v3.set(x,y,z).project(camera); if(v3.z>1){el.style.display="none";return;} el.style.display=""; el.style.left=((v3.x+1)/2*host.clientWidth)+"px"; el.style.top=((1-v3.y)/2*host.clientHeight)+"px"; }

/* ---------- 2D graphs ---------- */
function prepCanvas(c){ const r=c.getBoundingClientRect(), d=Math.min(window.devicePixelRatio||1,2); c.width=Math.round(r.width*d); c.height=Math.round(r.height*d); const x=c.getContext("2d"); x.setTransform(d,0,0,d,0,0); return [x,r.width,r.height]; }
const C={fg:"#dbe3e8",mu:"#82929d",ln:"#253139",am:"#f2a93b",cy:"#58cfe0",rd:"#ef5a4c"};
function axes(x,W,H,m,xmax,ymax,xt,yt,xl,yl,yfmt){
  x.font="11px 'Chivo Mono', monospace"; x.lineWidth=1;
  for(const v of xt){const px=m.l+(W-m.l-m.r)*v/xmax; x.strokeStyle=C.ln; x.beginPath(); x.moveTo(px,m.t); x.lineTo(px,H-m.b); x.stroke(); x.fillStyle=C.mu; x.textAlign="center"; x.fillText(+v.toPrecision(3),px,H-m.b+14);}
  for(const v of yt){const py=H-m.b-(H-m.t-m.b)*v/ymax; x.strokeStyle=C.ln; x.beginPath(); x.moveTo(m.l,py); x.lineTo(W-m.r,py); x.stroke(); x.fillStyle=C.mu; x.textAlign="right"; x.fillText(yfmt(v),m.l-6,py+4);}
  x.fillStyle=C.mu; x.textAlign="center"; x.fillText(xl,m.l+(W-m.l-m.r)/2,H-4);
  x.save(); x.translate(11,m.t+(H-m.t-m.b)/2); x.rotate(-Math.PI/2); x.fillText(yl,0,0); x.restore();
}
function niceTicks(max,n){ const raw=max/n, p=Math.pow(10,Math.floor(Math.log10(raw))), s=[1,2,2.5,5,10].find(k=>k*p>=raw)*p; const t=[]; for(let v=0;v<=max*1.0001;v+=s) t.push(+v.toPrecision(6)); return t; }
function unitFor(m){ return m>=0.1?[1,"m"]:m>=1e-4?[1e3,"mm"]:[1e6,"µm"]; }

function drawProfile(){
  const c=$("cProf"); const [x,W,H]=prepCanvas(c); x.clearRect(0,0,W,H);
  const xp=S.L*S.xs, d=delta(xp), ymax=1.5*d, [uf,un]=unitFor(ymax), Y=ymax*uf;
  const m={l:52,r:12,t:10,b:34}, xw=W-m.l-m.r, yh=H-m.t-m.b;
  axes(x,W,H,m,1.1,Y,[0,0.25,0.5,0.75,1],niceTicks(Y,4),"u / U∞","y ("+un+")",v=>+v.toPrecision(3));
  const X=u=>m.l+xw*u/1.1, YY=y=>H-m.b-yh*(y*uf)/Y;
  x.setLineDash([2,3]); x.strokeStyle=C.mu; x.beginPath(); x.moveTo(m.l,YY(d)); x.lineTo(W-m.r,YY(d)); x.stroke();
  x.fillStyle=C.mu; x.textAlign="left"; x.fillText("δ₉₉ = "+len(d),m.l+6,YY(d)-5); x.setLineDash([]);
  const turb=Rex(xp)>=S.Recr;
  if(S.cubic&&!turb){ x.setLineDash([6,4]); x.strokeStyle=C.cy; x.lineWidth=1.6; x.beginPath(); for(let i=0;i<=120;i++){const y=ymax*i/120; const px=X(cubic(y,xp)), py=YY(y); i?x.lineTo(px,py):x.moveTo(px,py);} x.stroke(); x.setLineDash([]); }
  x.strokeStyle=C.am; x.lineWidth=2.2; x.beginPath(); for(let i=0;i<=200;i++){const y=ymax*i/200; const px=X(uRatio(xp,y)), py=YY(y); i?x.lineTo(px,py):x.moveTo(px,py);} x.stroke();
  $("lgMain").textContent=turb?"1/7-power law (turbulent)":"Blasius (RK4 solution)";
  $("lgCub").hidden=turb||!S.cubic;
  $("profX").textContent="x = "+len(xp)+" · "+(turb?"turbulent":"laminar");
}
function drawGrowth(){
  const c=$("cGrow"); const [x,W,H]=prepCanvas(c); x.clearRect(0,0,W,H);
  const ymax=dMax*1.12, [uf,un]=unitFor(ymax), Y=ymax*uf;
  const m={l:52,r:12,t:10,b:34}, xw=W-m.l-m.r, yh=H-m.t-m.b;
  axes(x,W,H,m,S.L,Y,niceTicks(S.L,4),niceTicks(Y,4),"x (m)","δ ("+un+")",v=>+v.toPrecision(3));
  const X=v=>m.l+xw*v/S.L, YY=y=>H-m.b-yh*(y*uf)/Y;
  const xt=xTr();
  if(xt<S.L){ x.setLineDash([5,4]); x.strokeStyle=C.rd; x.lineWidth=1.2; x.beginPath(); x.moveTo(X(xt),m.t); x.lineTo(X(xt),H-m.b); x.stroke(); x.setLineDash([]); }
  const curve=(f,col,dash,w)=>{ x.setLineDash(dash); x.strokeStyle=col; x.lineWidth=w; x.beginPath(); for(let i=0;i<=300;i++){const xv=S.L*i/300; const px=X(xv), py=YY(f(xv)); i?x.lineTo(px,py):x.moveTo(px,py);} x.stroke(); x.setLineDash([]); };
  curve(dispThick,C.mu,[5,4],1.4); curve(delta,C.cy,[],2.2);
  const xs=S.L*S.xs; x.strokeStyle=C.am; x.lineWidth=1.2; x.beginPath(); x.moveTo(X(xs),m.t); x.lineTo(X(xs),H-m.b); x.stroke();
  x.fillStyle=C.am; x.beginPath(); x.arc(X(xs),YY(delta(xs)),4,0,Math.PI*2); x.fill();
}
function drawGraphs(){ drawProfile(); drawGrowth(); }

/* ---------- readout ---------- */
function readout(){
  const xp=S.L*S.xs, R=Rex(xp), RL=Rex(S.L), q=0.5*rho()*S.U*S.U, cf=cfx(xp), cd=CD(), xt=xTr();
  const rows=[
    ["sec","Station x = "+len(xp)],
    ["Local Reynolds number <small>Re_x</small>",sci(R)],
    ["Regime",R>=S.Recr?"turbulent":"laminar"],
    ["Thickness <small>δ₉₉</small>",len(delta(xp))],
    ["Displacement thickness <small>δ*</small>",len(dispThick(xp))],
    ["Momentum thickness <small>θ</small>",len(momThick(xp))],
    ["Shape factor <small>H = δ*/θ</small>",(dispThick(xp)/momThick(xp)).toFixed(2)],
    ["Skin friction <small>C_f,x</small>",sci(cf,3)],
    ["Wall shear stress <small>τ_w</small>",sci(cf*q,3)+" Pa"],
    ["sec","Whole plate, one side"],
    ["Plate Reynolds number <small>Re_L</small>",sci(RL)],
    ["Transition point <small>x_cr</small>",xt<S.L?len(xt):"beyond plate"],
    ["Drag coefficient <small>C_D</small>",sci(cd,3)],
    ["Drag per metre width <small>F/b</small>",sci(cd*q*S.L,3)+" N/m"],
  ];
  const html=a=>a.map(r=>r[0]==="sec"?`<tr><td class="sec" colspan="2">${r[1]}</td></tr>`:`<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join("");
  $("tbl").innerHTML=html(rows.slice(0,9)); $("tbl2").innerHTML=html(rows.slice(9));
  const ch=$("chipRegime"); const trOn=xt<S.L;
  ch.textContent=trOn?"TRANSITION AT x = "+len(xt):"LAMINAR OVER WHOLE PLATE"; ch.className="chip "+(trOn?"tr":"lam");
  $("chipRe").textContent="Re_L = "+sci(RL); $("chipWarn").hidden=RL>=1000;
}

/* ---------- update ---------- */
function update(){
  syncSP();
  $("oU").textContent=S.U.toFixed(S.U<1?3:2)+" m/s"; $("oL").textContent=S.L.toFixed(2)+" m";
  $("oRc").textContent=sci(S.Recr,1); $("oX").textContent=S.xs.toFixed(2);
  $("fluidProps").textContent="ν = "+sci(nu(),3)+" m²/s · ρ = "+rho()+" kg/m³";
  remap(); buildEnvelope(); buildProfiles();
  envMesh.visible=edge.visible=edge2.visible=S.env; tracers.visible=S.tracers;
  const xt=xTr(); trLine.visible=xt<S.L; trLine.position.x=xt*SX;
  const ex=SY/SX; $("hExag").textContent="×"+(ex>=10?Math.round(ex):ex.toFixed(1));
  readout(); drawGraphs(); if(typeof renderEq==="function") renderEq();
}

/* ---------- loop ---------- */
const clock=new THREE.Clock(), rgb=[0,0,0]; let tAcc=0;
function frame(){
  requestAnimationFrame(frame);
  if(viewSim.hidden) return; // equations tab open: keep state, skip drawing
  const dt=Math.min(clock.getDelta(),0.05);
  if(!S.paused&&S.tracers){
    tAcc+=dt; const vb=1.0+2.6*Math.log(S.U/UMIN)/Math.log(UMAX/UMIN), xt=xTr();
    for(let i=0;i<N;i++){
      const k=i*3; let xd=pos[k]; const xp=xd/SX; let yd=py0[i];
      const u=uRatio(Math.min(xp,S.L),yd/SY);
      if(xp>xt&&xp<=S.L&&u<0.995){ yd+=Math.sin(tAcc*5+ph[i])*0.06*(1-u)*YTOP; if(yd<0.01) yd=0.01; }
      xd+=Math.max(u,0.012)*vb*dt;
      if(xd>XMAX){ xd=XMIN+Math.random()*0.4; py0[i]=spawnY(); pos[k+2]=(Math.random()-0.5)*PZ*0.96; }
      pos[k]=xd; pos[k+1]=yd; cmap(u,rgb); col[k]=rgb[0]; col[k+1]=rgb[1]; col[k+2]=rgb[2];
    }
    tg.attributes.position.needsUpdate=true; tg.attributes.color.needsUpdate=true;
  }
  controls.update(); renderer.render(scene,camera);
  placeTag($("tagLE"),0,0,PZ/2); placeTag($("tagFlow"),-1.4,2.2,0);
  const xt=xTr(); const tr=$("tagTR"); tr.hidden=!(xt<S.L);
  if(!tr.hidden){ tr.textContent="Transition Re_x = "+sci(S.Recr,1); placeTag(tr,xt*SX,0,-PZ/2); }
}

/* ---------- equations & calculations view (uses the same PHYS module) ---------- */
const viewSim=$("viewSim"), viewEq=$("viewEq");
const EF={fluid:$("eFluid"),U:$("eU"),L:$("eL"),x:$("eX"),Recr:$("eRc"),y:$("eY")};
let eqMode="live", manualSeeded=false; const MAN={};
function nf(v,sig=4){
  if(!isFinite(v)) return "—"; if(v===0) return "0";
  const a=Math.abs(v), sg=v<0?"−":""; let e=Math.floor(Math.log10(a));
  if(e>=-3&&e<5) return sg+(+a.toPrecision(sig)).toString();
  let m=+(a/Math.pow(10,e)).toPrecision(sig); if(m>=10){m=+(m/10).toPrecision(sig);e++;}
  return sg+m+"×10<sup>"+(e<0?"−":"")+Math.abs(e)+"</sup>";
}
const lenH=m=>nf(m)+" m = "+nf(m*1e3)+" mm";
const fr=(n,d)=>`<span class="fr"><span>${n}</span><span>${d}</span></span>`;
const rt=x=>`<span class="rt">√<span>${x}</span></span>`;
const Ui="U<sub>∞</sub>", Rx="Re<sub>x</sub>", RLs="Re<sub>L</sub>", Rc="Re<sub>cr</sub>";
const ON=t=>["on",t], OFF=t=>["off",t];
function calc(o){
  return `<article class="calc"><div class="ch"><h3>${o.name}</h3>${o.tag?`<span class="st ${o.tag[0]}">${o.tag[1]}</span>`:""}</div>`+
    `<div class="eqn">${o.eq}</div><ol class="steps">${o.steps.map(s=>`<li>${s}</li>`).join("")}</ol>`+
    `<div class="res">${o.sym} = <b>${o.res}</b></div>${o.note?`<p class="note">${o.note}</p>`:""}</article>`;
}
const sec=(t,lede,items,warn)=>`<section class="card eqsec"><h2>${t}</h2>${lede?`<p class="lede${warn?" warn":""}">${lede}</p>`:""}<div class="calcs">${items.join("")}</div></section>`;

function fillFromSim(){ EF.fluid.value=S.fluid; EF.U.value=+S.U.toPrecision(4); EF.L.value=+S.L.toPrecision(4); EF.x.value=+(S.xs*S.L).toPrecision(4); EF.Recr.value=+S.Recr.toPrecision(3); }
function setMode(m){
  if(eqMode==="manual"&&m==="live") ["fluid","U","L","x","Recr"].forEach(k=>MAN[k]=EF[k].value);
  eqMode=m; const live=m==="live";
  document.querySelectorAll("#modeSeg button").forEach(b=>b.setAttribute("aria-pressed",String(b.dataset.m===m)));
  ["U","L","x","Recr"].forEach(k=>EF[k].readOnly=live); EF.fluid.disabled=live; $("bCopy").hidden=live;
  $("modeHint").textContent=live?"Values follow the Simulation tab sliders and update live. Height y is set here only.":"Your own values. The simulation does not change.";
  if(live) fillFromSim(); else if(!manualSeeded){ fillFromSim(); manualSeeded=true; } else Object.keys(MAN).forEach(k=>EF[k].value=MAN[k]);
  renderEq();
}
document.querySelectorAll("#modeSeg button").forEach(b=>b.addEventListener("click",()=>setMode(b.dataset.m)));
$("bCopy").addEventListener("click",()=>{fillFromSim();renderEq();});
Object.values(EF).forEach(el=>el.addEventListener(el.tagName==="SELECT"?"change":"input",renderEq));

function readNum(el,label,allowZero){
  const raw=el.value.trim(); let v=NaN, msg="";
  if(raw==="") msg=label+": enter a value.";
  else { v=Number(raw); if(!isFinite(v)) msg=label+": enter a number."; else if(allowZero?v<0:v<=0) msg=label+(allowZero?" cannot be negative.":" must be greater than 0."); }
  el.setAttribute("aria-invalid",msg?"true":"false"); return [v,msg];
}

function renderEq(){
  if(viewEq.hidden) return;
  const live=eqMode==="live"; if(live) fillFromSim();
  const errs=[], rd=(k,lab,z)=>{const r=readNum(EF[k],lab,z); if(r[1]) errs.push(r[1]); return r[0];};
  let fluid, U, L, x, Recr;
  if(live){ fluid=S.fluid; U=S.U; L=S.L; x=S.xs*S.L; Recr=S.Recr; ["U","L","x","Recr"].forEach(k=>EF[k].setAttribute("aria-invalid","false")); }
  else { fluid=EF.fluid.value; U=rd("U","U∞"); L=rd("L","Plate length L"); x=rd("x","Distance x"); Recr=rd("Recr","Critical Re"); }
  const yMM=rd("y","Height y",true), F=FLUIDS[fluid];
  $("eProps").innerHTML=`ν = <b>${nf(F.nu)} m²/s</b> (kinematic)<br>ρ = <b>${nf(F.rho)} kg/m³</b><br>μ = ρν = <b>${nf(F.rho*F.nu)} Pa·s</b> (dynamic)`;
  if(errs.length){ $("eErr").innerHTML=errs.join("<br>"); $("eErr").hidden=false; $("eqOut").innerHTML=`<section class="card"><p class="lede warn">Fix the highlighted input to see the calculations.</p></section>`; return; }
  $("eErr").hidden=true;
  const P={U,L,nu:F.nu,rho:F.rho,Recr}, y=yMM/1000;
  const Re=PHYS.Re(P,x), ReL=PHYS.Re(P,L), sRe=Math.sqrt(Re), sReL=Math.sqrt(ReL), Re5=Math.pow(Re,0.2), ReL5=Math.pow(ReL,0.2);
  const turbX=PHYS.isTurb(P,x), lamPlate=ReL<=Recr, xcr=PHYS.xcr(P), q=0.5*F.rho*U*U;
  const tagX=t=>t?ON("Used at this x"):OFF("Not used at this x");
  const out=[];

  // 2 · Reynolds numbers
  out.push(sec("2 · Reynolds numbers",
    (x>L?"x is past the end of the plate (x > L). Local values are still calculated. ":"")+`Regime at x: <b>${turbX?"turbulent":"laminar"}</b> (${Rx} ${turbX?"≥":"<"} ${Rc} = ${nf(Recr,3)}). ${ReL<1000?"Re<sub>L</sub> < 1000: boundary-layer theory is not reliable here.":""}`,[
    calc({name:"Local Reynolds number",sym:Rx,eq:`${Rx} = ${fr(Ui+" · x","ν")}`,
      steps:[`= ${fr(nf(U)+" × "+nf(x),nf(F.nu))}`,`= ${fr(nf(U*x)+" m²/s",nf(F.nu)+" m²/s")}`],res:nf(Re),tag:turbX?["off","turbulent"]:["on","laminar"]}),
    calc({name:"Plate Reynolds number",sym:RLs,eq:`${RLs} = ${fr(Ui+" · L","ν")}`,
      steps:[`= ${fr(nf(U)+" × "+nf(L),nf(F.nu))}`,`= ${fr(nf(U*L)+" m²/s",nf(F.nu)+" m²/s")}`],res:nf(ReL),tag:lamPlate?["on","laminar plate"]:["off","transition on plate"]}),
    calc({name:"Transition point",sym:"x<sub>cr</sub>",eq:`x<sub>cr</sub> = ${fr(Rc+" · ν",Ui)}`,
      steps:[`= ${fr(nf(Recr)+" × "+nf(F.nu),nf(U))}`],res:lenH(xcr),note:xcr>=L?"x<sub>cr</sub> ≥ L: the whole plate stays laminar.":"Flow is laminar for x < x<sub>cr</sub> and turbulent after it."})
  ],ReL<1000||x>L));

  // 3 · Laminar
  const lamC=(name,sym,c,f,unitLen)=>calc({name,sym,tag:tagX(!turbX),
    eq:unitLen?`${sym} = ${fr(c+" x",rt(Rx))}`:`${sym} = ${fr(c,rt(Rx))}`,
    steps:[unitLen?`= ${fr(c+" × "+nf(x),rt(nf(Re)))}`:`= ${fr(c,rt(nf(Re)))}`,`√${Rx} = ${nf(sRe)}`,unitLen?`= ${fr(nf(c*x)+" m",nf(sRe))}`:`= ${fr(c,nf(sRe))}`],
    res:unitLen?lenH(f):nf(f)+" (–)"});
  const eta=PHYS.lam.eta(P,x,y), ub=PHYS.lam.u(P,x,y);
  out.push(sec("3 · Laminar boundary layer (Blasius)",`Valid for ${Rx} < ${Rc}. The simulation uses these values for every station before the transition point.`,[
    lamC("Thickness δ₉₉","δ<sub>99</sub>",4.91,PHYS.lam.d(P,x),true),
    lamC("Displacement thickness","δ*",1.7208,PHYS.lam.ds(P,x),true),
    lamC("Momentum thickness","θ",0.664,PHYS.lam.th(P,x),true),
    calc({name:"Shape factor",sym:"H",tag:tagX(!turbX),eq:`H = ${fr("δ*","θ")} = ${fr("1.7208","0.664")}`,
      steps:[`= ${fr(nf(PHYS.lam.ds(P,x))+" m",nf(PHYS.lam.th(P,x))+" m")}`],res:nf(PHYS.lam.ds(P,x)/PHYS.lam.th(P,x))+" (–)",note:"Constant for Blasius flow."}),
    lamC("Local skin-friction coefficient","C<sub>f,x</sub>",0.664,PHYS.lam.cf(P,x),false),
    calc({name:"Average drag coefficient (whole plate)",sym:"C<sub>D</sub>",tag:lamPlate?ON("Used for this plate"):OFF("Not used: Re_L > Re_cr"),
      eq:`C<sub>D</sub> = ${fr("1.328",rt(RLs))}`,steps:[`= ${fr("1.328",rt(nf(ReL)))}`,`√${RLs} = ${nf(sReL)}`],res:nf(PHYS.lam.CD(P))+" (–)",
      note:"Average over the plate length, one side. Not the same as the local C<sub>f,x</sub>."}),
    calc({name:"Velocity at height y (Blasius)",sym:"u/U<sub>∞</sub>",tag:tagX(!turbX),
      eq:`η = y ${rt(fr(Ui,"ν x"))}, &nbsp; u/U<sub>∞</sub> = f′(η)`,
      steps:[`${fr(Ui,"ν x")} = ${fr(nf(U),nf(F.nu)+" × "+nf(x))} = ${nf(U/(F.nu*x))} m⁻²`,`η = ${nf(y)} m × ${nf(Math.sqrt(U/(F.nu*x)))} m⁻¹ = ${nf(eta)}`,`f′(${nf(eta,3)}) from the RK4 solution of f‴ + ½ f f″ = 0`],
      res:nf(ub)+" (–)",note:"η ≈ 4.91 is the edge of the layer (u = 0.99 U∞)."})
  ]));

  // 4 · Turbulent
  const dT=PHYS.tur.d(P,x), uT=PHYS.tur.u(P,x,y), A=PHYS.tur.A(P);
  out.push(sec("4 · Turbulent boundary layer (1/7-power law)",`Valid for ${Rx} ≥ ${Rc}. These formulas assume turbulence from the leading edge, which is why δ jumps at the transition point. They are usually quoted for ${Rx} up to about 10<sup>7</sup>; check the limit in your textbook.`,[
    calc({name:"Thickness δ",sym:"δ",tag:tagX(turbX),eq:`δ = ${fr("0.37 x",Rx+"<sup>1/5</sup>")}`,
      steps:[`= ${fr("0.37 × "+nf(x),"("+nf(Re)+")<sup>1/5</sup>")}`,`${Rx}<sup>1/5</sup> = ${nf(Re5)}`,`= ${fr(nf(0.37*x)+" m",nf(Re5))}`],res:lenH(dT)}),
    calc({name:"Velocity at height y",sym:"u/U<sub>∞</sub>",tag:tagX(turbX),eq:`u/U<sub>∞</sub> = (y/δ)<sup>1/7</sup> &nbsp; for y < δ`,
      steps:[`y/δ = ${fr(nf(y)+" m",nf(dT)+" m")} = ${nf(y/dT)}`, y>=dT?"y ≥ δ, so the point is outside the layer: u/U∞ = 1":`(${nf(y/dT)})<sup>1/7</sup>`],res:nf(uT)+" (–)"}),
    calc({name:"Displacement and momentum thickness",sym:"δ*, θ",tag:tagX(turbX),eq:`δ* = ${fr("δ","8")}, &nbsp; θ = ${fr("7 δ","72")}`,
      steps:[`δ* = ${fr(nf(dT),"8")} = ${nf(PHYS.tur.ds(P,x)*1e3)} mm`,`θ = ${fr("7 × "+nf(dT),"72")} = ${nf(PHYS.tur.th(P,x)*1e3)} mm`,`H = δ*/θ = 72/56 = ${nf(PHYS.tur.ds(P,x)/PHYS.tur.th(P,x))}`],
      res:nf(PHYS.tur.ds(P,x)*1e3)+" mm, "+nf(PHYS.tur.th(P,x)*1e3)+" mm",note:"Both follow from integrating the 1/7-power profile."}),
    calc({name:"Local skin-friction coefficient",sym:"C<sub>f,x</sub>",tag:tagX(turbX),eq:`C<sub>f,x</sub> = ${fr("0.059",Rx+"<sup>1/5</sup>")}`,
      steps:[`= ${fr("0.059",nf(Re5))}`],res:nf(PHYS.tur.cf(P,x))+" (–)"}),
    calc({name:"Drag coefficient, turbulent from leading edge",sym:"C<sub>D,turb</sub>",tag:OFF("Reference only"),eq:`C<sub>D,turb</sub> = ${fr("0.074",RLs+"<sup>1/5</sup>")}`,
      steps:[`${RLs}<sup>1/5</sup> = (${nf(ReL)})<sup>1/5</sup> = ${nf(ReL5)}`,`= ${fr("0.074",nf(ReL5))}`],res:nf(PHYS.tur.CDfull(P))+" (–)",note:"Ignores the laminar part at the front of the plate. The simulation uses the mixed value instead."}),
    calc({name:"Laminar-part correction constant",sym:"A",tag:OFF("Depends on Re_cr only"),eq:`A = ${Rc} ( ${fr("0.074",Rc+"<sup>1/5</sup>")} − ${fr("1.328",rt(Rc))} )`,
      steps:[`= ${nf(Recr)} × ( ${fr("0.074",nf(Math.pow(Recr,0.2)))} − ${fr("1.328",nf(Math.sqrt(Recr)))} )`,`= ${nf(Recr)} × ( ${nf(0.074/Math.pow(Recr,0.2))} − ${nf(1.328/Math.sqrt(Recr))} )`],res:nf(A),note:"Gives the textbook value A ≈ 1742 when Re_cr = 5×10⁵."}),
    calc({name:"Average drag coefficient, mixed laminar + turbulent",sym:"C<sub>D</sub>",tag:lamPlate?OFF("Not used: Re_L ≤ Re_cr"):ON("Used for this plate"),
      eq:`C<sub>D</sub> = ${fr("0.074",RLs+"<sup>1/5</sup>")} − ${fr("A",RLs)}`,
      steps:[`= ${fr("0.074",nf(ReL5))} − ${fr(nf(A),nf(ReL))}`,`= ${nf(PHYS.tur.CDfull(P))} − ${nf(A/ReL)}`],res:nf(PHYS.tur.CD(P))+" (–)",note:"Average over the plate, one side, for a laminar front part and turbulent rear part."})
  ]));

  // 5 · Cubic
  const dC=PHYS.cub.d(P,x), uC=PHYS.cub.u(P,x,y), r=y/dC;
  out.push(sec("5 · Cubic velocity approximation",`Polynomial fit for laminar flow, shown as the dashed line on the profile graph. Compare it with the exact Blasius value.${turbX?" This x is turbulent, so the cubic profile does not apply here.":""}`,[
    calc({name:"Thickness (cubic)",sym:"δ<sub>c</sub>",tag:tagX(false),eq:`δ<sub>c</sub> = ${fr("4.64 x",rt(Rx))}`,
      steps:[`= ${fr("4.64 × "+nf(x),nf(sRe))}`],res:lenH(dC),note:"Reference only: the simulation uses the Blasius δ₉₉."}),
    calc({name:"Velocity at height y (cubic)",sym:"u/U<sub>∞</sub>",tag:OFF("Graph overlay"),eq:`u/U<sub>∞</sub> = 1.5 ${fr("y","δ<sub>c</sub>")} − 0.5 (${fr("y","δ<sub>c</sub>")})<sup>3</sup>`,
      steps:[`y/δ<sub>c</sub> = ${fr(nf(y),nf(dC))} = ${nf(r)}`, r>=1?"y ≥ δc, so u/U∞ = 1":`= 1.5 × ${nf(r)} − 0.5 × ${nf(r)}<sup>3</sup> = ${nf(1.5*r)} − ${nf(0.5*r*r*r)}`],res:nf(uC)+" (–)"}),
    calc({name:"Cubic versus Blasius at y",sym:"Δ",tag:OFF("Comparison"),eq:`Δ = (u/U<sub>∞</sub>)<sub>cubic</sub> − (u/U<sub>∞</sub>)<sub>Blasius</sub>`,
      steps:[`= ${nf(uC)} − ${nf(ub)}`],res:nf(uC-ub,3)+(ub>0?` (${nf(100*(uC-ub)/ub,3)} % of Blasius)`:""),note:"Small values show the cubic profile is a good approximation."})
  ]));

  // 6 · Derived & final (same functions as the simulation readout)
  const cfA=PHYS.cf(P,x), cdA=PHYS.CD(P), uA=PHYS.u(P,x,y);
  out.push(sec("6 · Derived quantities and final results",`These use the regime the simulation applies at this x (${turbX?"turbulent":"laminar"}) and for the plate (${lamPlate?"laminar":"mixed"}).${live?" They match the Readout on the Simulation tab.":""}`,[
    calc({name:"Dynamic pressure",sym:"q",eq:`q = ½ ρ ${Ui}<sup>2</sup>`,steps:[`= 0.5 × ${nf(F.rho)} × ${nf(U)}<sup>2</sup>`],res:nf(q)+" Pa"}),
    calc({name:"Wall shear stress",sym:"τ<sub>w</sub>",eq:`τ<sub>w</sub> = C<sub>f,x</sub> · q`,steps:[`C<sub>f,x</sub> = ${nf(cfA)} (${turbX?"turbulent":"laminar"})`,`= ${nf(cfA)} × ${nf(q)} Pa`],res:nf(cfA*q)+" Pa"}),
    calc({name:"Drag force per metre width, one side",sym:"F/b",eq:`F/b = C<sub>D</sub> · q · L`,steps:[`C<sub>D</sub> = ${nf(cdA)} (${lamPlate?"laminar":"mixed"})`,`= ${nf(cdA)} × ${nf(q)} Pa × ${nf(L)} m`],res:nf(cdA*q*L)+" N/m"}),
    calc({name:"Local speed at height y",sym:"u",eq:`u = (u/U<sub>∞</sub>) · ${Ui}`,steps:[`u/U<sub>∞</sub> = ${nf(uA)} (${turbX?"1/7-power law":"Blasius"})`,`= ${nf(uA)} × ${nf(U)} m/s`],res:nf(uA*U)+" m/s"})
  ])+"");
  const rows=[[Rx,nf(Re)],["Regime at x",turbX?"turbulent":"laminar"],["δ₉₉ (applied)",nf(PHYS.delta(P,x)*1e3)+" mm"],["δ*",nf(PHYS.ds(P,x)*1e3)+" mm"],["θ",nf(PHYS.th(P,x)*1e3)+" mm"],["H",nf(PHYS.ds(P,x)/PHYS.th(P,x),3)],
    ["C<sub>f,x</sub> (local)",nf(cfA)],["τ<sub>w</sub>",nf(cfA*q)+" Pa"],[RLs,nf(ReL)],["x<sub>cr</sub>",xcr<L?nf(xcr*1e3)+" mm":"beyond plate"],["C<sub>D</sub> (plate average)",nf(cdA)],["F/b",nf(cdA*q*L)+" N/m"],["u at y = "+nf(yMM)+" mm",nf(uA*U)+" m/s"]];
  out[out.length-1]=out[out.length-1].replace(/<\/section>$/,`<table class="ftable"><tbody>${rows.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join("")}</tbody></table></section>`);
  $("eqOut").innerHTML=out.join("");
}

/* ---------- view routing (#simulation / #equations) ---------- */
function showView(){
  const v=location.hash==="#equations"?"eq":location.hash==="#theory"?"th":"sim", sim=v==="sim", eq=v==="eq";
  viewSim.hidden=!sim; viewEq.hidden=!eq; $("viewTh").hidden=v!=="th";
  [["tabSim","sim"],["tabEq","eq"],["tabTh","th"]].forEach(([id,k])=>{ if(k===v) $(id).setAttribute("aria-current","page"); else $(id).removeAttribute("aria-current"); });
  sbBtn.hidden=!sim; if(!sim&&mobile.matches) setSidebar(false,false);
  if(eq){ if(eqMode==="live") fillFromSim(); renderEq(); } else if(sim){ resize(); }
}
window.addEventListener("hashchange",showView);
setMode("live"); showView();

if(window.matchMedia&&matchMedia("(prefers-reduced-motion: reduce)").matches){S.paused=true;$("bPlay").textContent="Play";}
syncInputs(); update(); resize(); frame();
// first paint of tracer colours when paused
(function(){const xt=0; for(let i=0;i<N;i++){const k=i*3; pos[k+1]=py0[i]; cmap(uRatio(Math.min(pos[k]/SX,S.L),py0[i]/SY),rgb); col[k]=rgb[0];col[k+1]=rgb[1];col[k+2]=rgb[2];} tg.attributes.position.needsUpdate=true; tg.attributes.color.needsUpdate=true;})();
document.fonts&&document.fonts.ready.then(drawGraphs);
})();
