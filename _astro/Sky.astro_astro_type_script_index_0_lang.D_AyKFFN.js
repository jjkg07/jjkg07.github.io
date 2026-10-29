var e=class{w;h;stride;u;v;dye;u0;v0;dye0;p;div;curl;iterations=14;vorticity=10;velocityDecay=.5;dyeDecay=.12;windU=0;windV=0;constructor(e,t){this.w=e,this.h=t,this.stride=e+2;let n=(e+2)*(t+2);this.u=new Float32Array(n),this.v=new Float32Array(n),this.dye=new Float32Array(n),this.u0=new Float32Array(n),this.v0=new Float32Array(n),this.dye0=new Float32Array(n),this.p=new Float32Array(n),this.div=new Float32Array(n),this.curl=new Float32Array(n)}setWind(e,t){this.u.fill(e),this.v.fill(t),this.windU=e,this.windV=t}splat(e,t,r,i,a,o){let{w:s,h:c,stride:l,u,v:d,dye:f}=this,p=o*o,m=Math.min(Math.ceil(o*2.5),(Math.min(s,c)>>1)-1),h=Math.floor(e),g=Math.floor(t);for(let o=g-m;o<=g+m;o++){let g=o+.5-t,_=l*n(o,c);for(let t=h-m;t<=h+m;t++){let o=t+.5-e,c=Math.exp(-(o*o+g*g)/p);if(c<.01)continue;let l=_+n(t,s);u[l]+=r*c,d[l]+=i*c,f[l]=Math.min(f[l]+a*c,4)}}}stir(e,t,r,i,a,o){let{w:s,h:c,stride:l,u,v:d}=this,f=o*o,p=Math.min(Math.ceil(o*2.5),(Math.min(s,c)>>1)-1),m=Math.floor(e),h=Math.floor(t);for(let o=h-p;o<=h+p;o++){let h=o+.5-t,g=l*n(o,c);for(let t=m-p;t<=m+p;t++){let o=t+.5-e,c=Math.exp(-(o*o+h*h)/f)*a;if(c<.01)continue;let l=g+n(t,s);u[l]+=(r-u[l])*c,d[l]+=(i-d[l])*c}}}sample(e,t,n){n[0]=this.bilinear(this.u,e,t),n[1]=this.bilinear(this.v,e,t)}sampleDye(e,t){return this.bilinear(this.dye,e,t)}step(e){this.setBoundary(this.u),this.setBoundary(this.v),this.setBoundary(this.dye),this.confineVorticity(e),this.project(),this.u0.set(this.u),this.v0.set(this.v),this.advect(this.u,this.u0,this.u0,this.v0,e),this.advect(this.v,this.v0,this.u0,this.v0,e),this.project();let t=Math.exp(-e*this.velocityDecay),n=Math.exp(-e*this.dyeDecay),{u:r,v:i,windU:a,windV:o}=this;for(let e=0;e<r.length;e++)r[e]=a+(r[e]-a)*t,i[e]=o+(i[e]-o)*t;this.dye0.set(this.dye),this.advect(this.dye,this.dye0,this.u,this.v,e);let s=this.dye;for(let e=0;e<s.length;e++)s[e]*=n}bilinear(e,n,r){let i=t(n+.5,this.w),a=t(r+.5,this.h),o=i|0,s=a|0,c=i-o,l=a-s,u=this.stride,d=o+u*s;return(1-c)*((1-l)*e[d]+l*e[d+u])+c*((1-l)*e[d+1]+l*e[d+1+u])}advect(e,n,r,i,a){let{w:o,h:s,stride:c}=this;for(let l=1;l<=s;l++)for(let u=1;u<=o;u++){let d=u+c*l,f=t(u-a*r[d],o),p=t(l-a*i[d],s),m=f|0,h=p|0,g=f-m,_=p-h,v=1-g,y=1-_,b=m+c*h;e[d]=v*(y*n[b]+_*n[b+c])+g*(y*n[b+1]+_*n[b+1+c])}this.setBoundary(e)}project(){let{w:e,h:t,stride:n,u:r,v:i,p:a,div:o}=this;for(let s=1;s<=t;s++)for(let t=1;t<=e;t++){let e=t+n*s;o[e]=-.5*(r[e+1]-r[e-1]+i[e+n]-i[e-n]),a[e]=0}this.setBoundary(o),this.setBoundary(a);for(let r=0;r<this.iterations;r++){for(let r=1;r<=t;r++)for(let t=1;t<=e;t++){let e=t+n*r;a[e]=(o[e]+a[e-1]+a[e+1]+a[e-n]+a[e+n])*.25}this.setBoundary(a)}for(let o=1;o<=t;o++)for(let t=1;t<=e;t++){let e=t+n*o;r[e]-=.5*(a[e+1]-a[e-1]),i[e]-=.5*(a[e+n]-a[e-n])}this.setBoundary(r),this.setBoundary(i)}confineVorticity(e){let{w:t,h:n,stride:r,u:i,v:a,curl:o}=this;for(let e=1;e<=n;e++)for(let n=1;n<=t;n++){let t=n+r*e;o[t]=.5*(a[t+1]-a[t-1]-(i[t+r]-i[t-r]))}this.setBoundary(o);let s=this.vorticity*e;for(let e=1;e<=n;e++)for(let n=1;n<=t;n++){let t=n+r*e,c=.5*(Math.abs(o[t+1])-Math.abs(o[t-1])),l=.5*(Math.abs(o[t+r])-Math.abs(o[t-r])),u=Math.hypot(c,l)+1e-5,d=o[t];i[t]+=s*l*d/u,a[t]-=s*c*d/u}}setBoundary(e){let{w:t,h:n,stride:r}=this,i=r,a=r*n,o=r*(n+1);for(let n=1;n<=t;n++)e[n]=e[a+n],e[o+n]=e[i+n];for(let i=0;i<=n+1;i++){let n=r*i;e[n]=e[n+t],e[n+t+1]=e[n+1]}}};function t(e,t){let n=(e-.5)%t;return(n<0?n+t:n)+.5}function n(e,t){let n=e%t;return(n<0?n+t:n)+1}var r={perception:90,separation:20,minSpeed:70,maxSpeed:150,maxForce:65,alignWeight:1.5,cohesionWeight:.6,separationWeight:3.7,wrapMargin:80,beatRate:1.7,flapEffort:1,glide:4.2,beats:2.5},i=220,a=.3,o=class{count=0;x=new Float32Array;y=new Float32Array;vx=new Float32Array;vy=new Float32Array;ax=new Float32Array;ay=new Float32Array;boost=new Float32Array;wing=new Float32Array;depth=new Float32Array;phase=new Float32Array;amp=new Float32Array;effort=new Float32Array;burst=new Float32Array;tempo=new Float32Array;params;width=0;height=0;worldW=0;worldH=0;cellW=0;cellH=0;cellStart=new Int32Array;cellOf=new Int32Array;order=new Int32Array;cols=0;rows=0;gridPerception=0;constructor(e=r){this.params=e}resize(e,t,n){let r=this.count,i=Math.min(r,n),a=e=>{let t=new Float32Array(n);return t.set(e.subarray(0,i)),t};this.x=a(this.x),this.y=a(this.y),this.vx=a(this.vx),this.vy=a(this.vy),this.wing=a(this.wing),this.depth=a(this.depth),this.phase=a(this.phase),this.amp=a(this.amp),this.effort=a(this.effort),this.burst=a(this.burst),this.tempo=a(this.tempo),this.ax=new Float32Array(n),this.ay=new Float32Array(n),this.boost=new Float32Array(n),this.cellOf=new Int32Array(n),this.order=new Int32Array(n);let{minSpeed:o,maxSpeed:c,wrapMargin:l}=this.params;this.width=e,this.height=t,this.worldW=e+l*2,this.worldH=t+l*2;let u=Array.from({length:Math.max(3,Math.round(n/35))},()=>({x:Math.random()*e,y:Math.random()*t,heading:Math.random()*Math.PI*2}));for(let e=0;e<n;e++)if(e>=i){let t=u[e%u.length],n=Math.sqrt(Math.random())*45,r=Math.random()*Math.PI*2;this.x[e]=s(t.x+Math.cos(r)*n,l,this.worldW),this.y[e]=s(t.y+Math.sin(r)*n,l,this.worldH);let i=t.heading+(Math.random()-.5)*.8,a=o+Math.random()*(c-o);this.vx[e]=Math.cos(i)*a,this.vy[e]=Math.sin(i)*a,this.depth[e]=Math.random(),this.tempo[e]=.85+Math.random()*.3,this.phase[e]=Math.random()*Math.PI*2,this.burst[e]=Math.random()<.4?Math.random()*this.burstLength():-Math.random()*this.glideLength()}else this.x[e]=s(this.x[e],l,this.worldW),this.y[e]=s(this.y[e],l,this.worldH);this.count=n,this.layoutGrid()}layoutGrid(){let{perception:e}=this.params;this.gridPerception=e,this.cols=Math.max(1,Math.floor(this.worldW/e)),this.rows=Math.max(1,Math.floor(this.worldH/e)),this.cellW=this.worldW/this.cols,this.cellH=this.worldH/this.rows,this.cellStart=new Int32Array(this.cols*this.rows+1)}buildGrid(){this.params.perception!==this.gridPerception&&this.layoutGrid();let{count:e,x:t,y:n,cols:r,rows:i,cellW:a,cellH:o,cellStart:s,cellOf:c,order:l}=this,u=this.params.wrapMargin;s.fill(0);for(let l=0;l<e;l++){let e=Math.min(r-1,Math.max(0,(t[l]+u)/a|0))+Math.min(i-1,Math.max(0,(n[l]+u)/o|0))*r;c[l]=e,s[e+1]++}for(let e=0;e<r*i;e++)s[e+1]+=s[e];let d=s.slice(0,r*i);for(let t=0;t<e;t++)l[d[c[t]]++]=t}step(e){this.buildGrid();let{count:t,x:n,y:r,vx:i,vy:a,ax:o,ay:d,boost:f,cols:p,rows:m,cellStart:h,cellOf:g,order:_,worldW:v,worldH:y}=this,b=this.params,x=v/2,S=y/2,C=new Int32Array(Math.min(3,p)),w=new Int32Array(Math.min(3,m)),T=b.perception*b.perception,E=b.separation*b.separation;for(let s=0;s<t;s++){let t=n[s],D=r[s],O=0,k=0,A=0,j=0,M=0,N=0,P=0,F=g[s];c(F%p,p,C),c(F/p|0,m,w);for(let e of w)for(let o of C){let c=o+e*p;for(let e=h[c];e<h[c+1];e++){let o=_[e];if(o===s)continue;let c=n[o]-t,l=r[o]-D;c>x?c-=v:c<-x&&(c+=v),l>S?l-=y:l<-S&&(l+=y);let u=c*c+l*l;u>T||(O++,k+=i[o],A+=a[o],j+=c,M+=l,u<E&&u>1e-4&&(N-=c/u,P-=l/u))}}let I=o[s],L=d[s];O>0&&(u(k,A,i[s],a[s],b,b.alignWeight),u(j,M,i[s],a[s],b,b.cohesionWeight),(N!==0||P!==0)&&u(N,P,i[s],a[s],b,b.separationWeight),I+=l[0],L+=l[1],l[0]=l[1]=0);let R=i[s]+I*e,z=a[s]+L*e,B=Math.hypot(R,z)||1,V=b.maxSpeed+f[s];B>V?(R*=V/B,z*=V/B):B<b.minSpeed&&(R*=b.minSpeed/B,z*=b.minSpeed/B),this.flap(s,R,z,e),i[s]=R,a[s]=z,o[s]=0,d[s]=0,f[s]*=Math.exp(-e*1.5)}let D=b.wrapMargin;for(let o=0;o<t;o++)n[o]=s(n[o]+i[o]*e,D,v),r[o]=s(r[o]+a[o]*e,D,y)}burstLength(){return this.params.beats*(.6+Math.random()*.8)/Math.max(this.params.beatRate,.1)}glideLength(){return this.params.glide*(.5+Math.random())}flap(e,t,n,r){let o=this.params,s=this.vx[e],c=this.vy[e],l=Math.hypot(s,c)||1,u=(t-s)/r,d=(n-c)/r,f=Math.max(0,(u*s+d*c)/l),p=Math.abs(s*d-c*u)/l,m=Math.max(0,-c/l)*.15,h=this.boost[e]/Math.max(o.maxSpeed,1),g=Math.min(1,((f+.6*p)/i+h)*o.flapEffort+m),_=this.effort[e]+=(g-this.effort[e])*(1-Math.exp(-r/.25)),v=this.burst[e];v>0?(v-=r,v<=0&&(v=-this.glideLength())):(v+=r,v>=0&&(v=this.burstLength())),this.burst[e]=v;let y=Math.max(+(v>0),Math.min(1,Math.max(0,(_-a)*2.5))),b=y>this.amp[e]?8:2.5,x=this.amp[e]+=(y-this.amp[e])*(1-Math.exp(-r*b)),S=(this.phase[e]+Math.PI*2*o.beatRate*this.tempo[e]*(1+.6*_)*r)%(Math.PI*2);this.phase[e]=S,this.wing[e]=x*Math.sin(S)}};function s(e,t,n){let r=(e+t)%n;return(r<0?r+n:r)-t}function c(e,t,n){if(t<3)for(let e=0;e<t;e++)n[e]=e;else n[0]=(e+t-1)%t,n[1]=e,n[2]=(e+1)%t}var l=new Float32Array(2);function u(e,t,n,r,i,a){let o=Math.hypot(e,t)||1,s=e/o*i.maxSpeed-n,c=t/o*i.maxSpeed-r,u=Math.hypot(s,c),d=(u>i.maxForce?i.maxForce/u:1)*a;l[0]+=s*d,l[1]+=c*d}var d={amount:.15,density:.7,size:6.5,stretch:8,puffs:12,life:32},f=class{fluid;params;clouds=[];constructor(e,t){this.fluid=e,this.params=t;for(let e=0;e<this.target;e++){let e=this.spawn();e.age=Math.random()*e.life*.7;let t=e.age/e.life,n=(t/2-Math.sin(2*Math.PI*t)/(4*Math.PI))*e.life;this.emit(e,n),this.clouds.push(e)}}get target(){let{w:e,h:t}=this.fluid;return Math.max(1,Math.round(e*t*this.params.density/1e3))}step(e){let t=[0,0];for(let n=this.clouds.length-1;n>=0;n--){let r=this.clouds[n];if(r.age+=e,r.age>=r.life){this.clouds.splice(n,1);continue}this.fluid.sample(r.x,r.y,t),r.x+=t[0]*e,r.y+=t[1]*e;let i=r.age/r.life;this.emit(r,(.5-.5*Math.cos(2*Math.PI*i))*e)}this.clouds.length<this.target&&this.clouds.push(this.spawn())}emit(e,t){let n=this.params.amount*t;for(let t of e.puffs)this.fluid.splat(e.x+t.dx,e.y+t.dy,0,0,n*t.weight,t.radius)}spawn(){let{w:e,h:t}=this.fluid,{size:n,stretch:r,puffs:i,life:a}=this.params,o=n*(.65+Math.random()*.7),s=Math.sqrt(r)*1.4,c=1.4/Math.sqrt(r),l=Math.max(1,Math.round(i*(.65+Math.random()*.7))),u=Array.from({length:l},()=>({dx:(Math.random()-.5)*o*s,dy:(Math.random()-.5)*o*c,radius:o*(.35+Math.random()*.35),weight:.5+Math.random()*.5}));return{x:Math.random()*e,y:Math.random()*t,age:0,life:a*(.7+Math.random()*.6),puffs:u}}};function p(e){let t=[0];for(let n=1;n<e;n*=2){let e=Array(n*2*n*2);for(let r=0;r<n;r++)for(let i=0;i<n;i++){let a=t[r*n+i]*4;e[r*2*n+i]=a,e[r*2*n+i+n]=a+2,e[(r+n)*2*n+i]=a+3,e[(r+n)*2*n+i+n]=a+1}t=e}let n=e*e;return{size:e,values:Float32Array.from(t,e=>(e+.5)/n)}}var m=p(8);function h(e=64,t=1.5,n=7){let r=`${e}:${t}:${n}`,i=g.get(r);if(i)return i;let a=e*e,o=_(e,t),s=new Uint8Array(a),c=new Float32Array(a),l=new Int32Array(a),u=(t,n)=>{let r=t%e,i=t/e|0;for(let t=0;t<e;t++){let a=(t-i+e)%e*e,s=t*e;for(let t=0;t<e;t++)c[s+t]+=n*o[a+(t-r+e)%e]}},d=()=>{let e=-1;for(let t=0;t<a;t++)s[t]&&(e<0||c[t]>c[e])&&(e=t);return e},f=()=>{let e=-1;for(let t=0;t<a;t++)!s[t]&&(e<0||c[t]<c[e])&&(e=t);return e},p=(e,t)=>{s[e]=+!!t,u(e,t?1:-1)},m=v(n),h=Math.max(1,Math.round(a*.1));for(let e=0;e<h;){let t=Math.floor(m()*a);s[t]||(p(t,!0),e++)}for(;;){let e=d();p(e,!1);let t=f();if(p(t,!0),t===e)break}let y=s.slice(),b=c.slice();for(let e=h;e>0;e--){let t=d();p(t,!1),l[t]=e-1}s.set(y),c.set(b);for(let e=h;e<a;e++){let t=f();p(t,!0),l[t]=e}let x={size:e,values:Float32Array.from(l,e=>(e+.5)/a)};return g.set(r,x),x}var g=new Map;function _(e,t){let n=new Float32Array(e*e),r=2*t*t;for(let t=0;t<e;t++){let i=Math.min(t,e-t);for(let a=0;a<e;a++){let o=Math.min(a,e-a);n[t*e+a]=Math.exp(-(o*o+i*i)/r)}}return n}function v(e){let t=e>>>0;return()=>{t=t+1831565813>>>0;let e=t;return e=Math.imul(e^e>>>15,e|1),e^=e+Math.imul(e^e>>>7,e|61),((e^e>>>14)>>>0)/4294967296}}var y={opacity:.6,thickness:1.4,birdSize:5.6},b=.05,x=(e,t)=>Math.max(0,(1-Math.exp(-e*t)-b)/.95),S=m.values,C=([e,t,n],r=255)=>(r<<24|n<<16|t<<8|e)>>>0;function w(e,t,n,r,i,a,o,s){let{count:c,x:l,y:u,vx:d,vy:f}=e,p=o,m=1.5/3.6*p;for(let e=0;e<c;e++){let o=Math.hypot(d[e],f[e])||1,c=d[e]/o,h=f[e]/o,g=l[e]/i,_=u[e]/a,v=Math.max(0,Math.floor(g-p-1)),y=Math.min(n-1,Math.ceil(g+p+1)),b=Math.max(0,Math.floor(_-p-1)),x=Math.min(r-1,Math.ceil(_+p+1));for(let e=b;e<=x;e++)for(let r=v;r<=y;r++){let i=r+.5-g,a=e+.5-_,o=-(i*c+a*h)+.6,l=Math.abs(-i*h+a*c);o>=0&&o<=p&&l<=o/p*m&&(t[r+e*n]=s)}}}function T(e,t,n,r=3){let i=e.getContext(`2d`),{size:a,values:o}=t,s=a-1,c=Math.log2(a),l=document.createElement(`canvas`),u=l.getContext(`2d`),d,f,p=0,m=0,h=r,g=r,_=r,v=0,y=0;return{resize(e,t,n,i){h=Math.max(1,Math.round(n/e*r)),g=h*e/n,_=h*t/i,p=Math.ceil(n/h),m=Math.ceil(i/h),l.width=p,l.height=m,d=u.createImageData(p,m),f=new Uint32Array(d.data.buffer)},setColors(e,t){v=C(e),y=C(t)},draw(t,r,a){let{opacity:b,thickness:S}=n,C=g/r,T=_/r,E=0;for(let e=0;e<m;e++){let n=(e+.5)*T,r=(e&s)<<c;for(let e=0;e<p;e++,E++){let i=x(t.sampleDye((e+.5)*C,n),S)*b;f[E]=i>o[r|e&s]?y:0}}w(a,f,p,m,g,_,n.birdSize,v),u.putImageData(d,0,0),i.setTransform(1,0,0,1,0,0),i.clearRect(0,0,e.width,e.height),i.imageSmoothingEnabled=!1,i.drawImage(l,0,0,p*h,m*h)}}}function E(e,t,n=8,r=14,i=12){let a=e.getContext(`2d`),o=document.createElement(`canvas`),s=o.getContext(`2d`),c=0,l=0,u=new Uint8Array,d=1,f=[0,0,0],p=[0,0,0];function m(){let e=Math.round(n*d),t=Math.round(r*d);o.width=e*18,o.height=t,s.setTransform(d,0,0,d,0,0),s.font=`${i}px ui-monospace, "SF Mono", Menlo, monospace`,s.textAlign=`center`,s.textBaseline=`middle`;for(let e=0;e<18;e++){let[t,i,a]=e<10?p:f;s.fillStyle=`rgb(${t} ${i} ${a})`,s.fillText(` .:-=+*#%@→↘↓↙←↖↑↗`[e],e*n+n/2,r/2)}}return{resize(e,t,i){c=Math.ceil(e/n),l=Math.ceil(t/r),u=new Uint8Array(c*l),d=i/e,m()},setColors(e,t){f=e,p=t,m()},draw(i,s,f){a.setTransform(1,0,0,1,0,0),a.clearRect(0,0,e.width,e.height),a.setTransform(d,0,0,d,0,0);let{opacity:p,thickness:m}=t;for(let e=0,t=0;e<l;e++){let a=(e+.5)*r/s,o=(e&7)<<3;for(let e=0;e<c;e++,t++){let r=x(i.sampleDye((e+.5)*n/s,a),m)*Math.min(1,p/.6);u[t]=Math.min(9,Math.floor(r*9+S[o|e&7]))}}let{count:h,x:g,y:_,vx:v,vy:y}=f;for(let e=0;e<h;e++){let t=Math.min(c-1,g[e]/n|0),i=Math.min(l-1,_[e]/r|0),a=Math.round(Math.atan2(y[e],v[e])/(Math.PI/4))&7;u[t+i*c]=10+a}let b=o.width/18,C=o.height;for(let e=0,t=0;e<l;e++)for(let i=0;i<c;i++,t++){let s=u[t];s!==0&&a.drawImage(o,s*b,0,b,C,i*n,e*r,n,r)}}}}var D=`
// ---- PCG hash ------------------------------------------------------------------
uvec3 pcg3d(uvec3 v) {
  v = v * 1664525u + 1013904223u;
  v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
  v ^= v >> 16u;
  v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
  return v;
}

/** Three independent uniform randoms in [0, 1) for an integer coordinate. */
vec3 hash33(uvec3 v) {
  return vec3(pcg3d(v)) * (1.0 / 4294967296.0);
}

// ---- simplex noise (Ashima / Gustavson) --------------------------------------------
vec2 mod289(vec2 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec3 mod289(vec3 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec4 mod289(vec4 x) { return x - floor(x * (1.0 / 289.0)) * 289.0; }
vec3 permute(vec3 x) { return mod289(((x * 34.0) + 1.0) * x); }
vec4 permute(vec4 x) { return mod289(((x * 34.0) + 1.0) * x); }
vec4 taylorInvSqrt(vec4 r) { return 1.79284291400159 - 0.85373472095314 * r; }

/** 2D simplex noise, roughly [-1, 1]. */
float snoise(vec2 v) {
  const vec4 C = vec4(0.211324865405187, 0.366025403784439, -0.577350269189626, 0.024390243902439);
  vec2 i = floor(v + dot(v, C.yy));
  vec2 x0 = v - i + dot(i, C.xx);
  vec2 i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;
  i = mod289(i);
  vec3 p = permute(permute(i.y + vec3(0.0, i1.y, 1.0)) + i.x + vec3(0.0, i1.x, 1.0));
  vec3 m = max(0.5 - vec3(dot(x0, x0), dot(x12.xy, x12.xy), dot(x12.zw, x12.zw)), 0.0);
  m = m * m;
  m = m * m;
  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;
  m *= 1.79284291400159 - 0.85373472095314 * (a0 * a0 + h * h);
  vec3 g;
  g.x = a0.x * x0.x + h.x * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}

/** 3D simplex noise, roughly [-1, 1]. Use z as time to evolve a 2D field smoothly. */
float snoise(vec3 v) {
  const vec2 C = vec2(1.0 / 6.0, 1.0 / 3.0);
  const vec4 D = vec4(0.0, 0.5, 1.0, 2.0);
  vec3 i = floor(v + dot(v, C.yyy));
  vec3 x0 = v - i + dot(i, C.xxx);
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min(g.xyz, l.zxy);
  vec3 i2 = max(g.xyz, l.zxy);
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy;
  vec3 x3 = x0 - D.yyy;
  i = mod289(i);
  vec4 p = permute(permute(permute(i.z + vec4(0.0, i1.z, i2.z, 1.0)) + i.y + vec4(0.0, i1.y, i2.y, 1.0)) + i.x + vec4(0.0, i1.x, i2.x, 1.0));
  float n_ = 0.142857142857;
  vec3 ns = n_ * D.wyz - D.xzx;
  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);
  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_);
  vec4 x = x_ * ns.x + ns.yyyy;
  vec4 y = y_ * ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);
  vec4 b0 = vec4(x.xy, y.xy);
  vec4 b1 = vec4(x.zw, y.zw);
  vec4 s0 = floor(b0) * 2.0 + 1.0;
  vec4 s1 = floor(b1) * 2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));
  vec4 a0 = b0.xzyw + s0.xzyw * sh.xxyy;
  vec4 a1 = b1.xzyw + s1.xzyw * sh.zzww;
  vec3 p0 = vec3(a0.xy, h.x);
  vec3 p1 = vec3(a0.zw, h.y);
  vec3 p2 = vec3(a1.xy, h.z);
  vec3 p3 = vec3(a1.zw, h.w);
  vec4 norm = taylorInvSqrt(vec4(dot(p0, p0), dot(p1, p1), dot(p2, p2), dot(p3, p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;
  vec4 m = max(0.6 - vec4(dot(x0, x0), dot(x1, x1), dot(x2, x2), dot(x3, x3)), 0.0);
  m = m * m;
  return 42.0 * dot(m * m, vec4(dot(p0, x0), dot(p1, x1), dot(p2, x2), dot(p3, x3)));
}

// ---- fBm ---------------------------------------------------------------------------
// each octave is rotated ~37° and shifted, so lattice directions and origins never line up
const mat2 OCTAVE_ROTATE = mat2(0.8, -0.6, 0.6, 0.8);
const mat3 OCTAVE_ROTATE3 = mat3(0.00, 0.80, 0.60, -0.80, 0.36, -0.48, -0.60, -0.48, 0.64);

/** Fractal Brownian motion over 2D simplex noise, roughly [-1, 1]. octaves ≤ 8. */
float fbm(vec2 p, int octaves) {
  float v = 0.0, a = 0.5;
  for (int i = 0; i < 8; i++) {
    if (i >= octaves) break;
    v += a * snoise(p);
    p = OCTAVE_ROTATE * p * 2.03 + vec2(13.7, 5.1);
    a *= 0.5;
  }
  return v;
}

/** Fractal Brownian motion over 3D simplex noise, roughly [-1, 1]. octaves ≤ 8. */
float fbm(vec3 p, int octaves) {
  float v = 0.0, a = 0.5;
  for (int i = 0; i < 8; i++) {
    if (i >= octaves) break;
    v += a * snoise(p);
    p = OCTAVE_ROTATE3 * p * 2.03 + vec3(13.7, 5.1, 3.3);
    a *= 0.5;
  }
  return v;
}

// ---- film grain --------------------------------------------------------------------
/**
 * Per-pixel, per-frame grain in [-1, 1] with a triangular distribution (sum of two
 * uniforms), which reads as film rather than TV static. Integer-hashed, so no lattice.
 */
float grain(ivec2 px, uint frame) {
  vec3 h = hash33(uvec3(uvec2(px), frame));
  return h.x + h.y - 1.0;
}
`,O=[0,.2,.52,1],k={name:`night`,stops:[`#1f2641`,`#363e5c`,`#4c5371`,`#555c76`],glow:`#3c4668`,glowAmount:0,horizon:-.01,mist:`#495170`,cloudSun:`#000000`,cloudShade:[1.3,1.32,1.42],bird:`#151a2b`},A=[{hour:3.4,...k},{hour:4.6,name:`blue hour`,stops:[`#5c5d86`,`#4f5682`,`#465078`,`#424a6e`],glow:`#a07aa4`,glowAmount:.5,horizon:.02,mist:`#6c6e96`,cloudSun:`#4a3a58`,cloudShade:[1.28,1.28,1.4],bird:`#1d1f33`},{hour:5.45,name:`dawn`,stops:[`#e8acaa`,`#bc9db6`,`#8e8fb2`,`#6f7ca3`],glow:`#f7b29e`,glowAmount:.75,horizon:.03,mist:`#d3bccb`,cloudSun:`#e9a2ab`,cloudShade:[1.16,1.08,1.16],bird:`#2c2840`},{hour:6.4,name:`sunrise`,stops:[`#ffc6a2`,`#efbfb2`,`#b4b8c4`,`#8ba3b9`],glow:`#ffb88f`,glowAmount:1,horizon:.04,mist:`#f6dccf`,cloudSun:`#ffc2a6`,cloudShade:[.98,.95,1],bird:`#3a3448`},{hour:7.9,name:`morning`,stops:[`#f4dfca`,`#d2d8d6`,`#a6c1cf`,`#80a7c2`],glow:`#ffe0bc`,glowAmount:.55,horizon:.04,mist:`#e8ebe6`,cloudSun:`#fff0e0`,cloudShade:[.86,.89,.95],bird:`#2f3646`},{hour:12.5,name:`midday`,stops:[`#e0eaee`,`#bad5e2`,`#8dbad3`,`#6a9cc2`],glow:`#fff3d8`,glowAmount:.45,horizon:.04,mist:`#e4edf2`,cloudSun:`#fff8ec`,cloudShade:[.86,.9,.96],bird:`#2a3442`},{hour:16.6,name:`afternoon`,stops:[`#efe4d2`,`#cfd6e0`,`#97bbd2`,`#6898b8`],glow:`#ffe4b4`,glowAmount:.6,horizon:.04,mist:`#eeeadc`,cloudSun:`#fff0d8`,cloudShade:[.87,.89,.93],bird:`#2c3440`},{hour:19.6,name:`golden hour`,stops:[`#ffc47c`,`#f3cfa6`,`#bcbcc4`,`#5f90aa`],glow:`#fdac2a`,glowAmount:.7,horizon:.03,mist:`#f5e6cc`,cloudSun:`#ffcf9c`,cloudShade:[1.02,.98,.95],bird:`#352f34`},{hour:20.35,name:`sunset`,stops:[`#fd8a74`,`#f4a096`,`#cf8e8c`,`#a07a7e`],glow:`#f3624b`,glowAmount:1,horizon:.02,mist:`#f2c6bd`,cloudSun:`#ffb4a0`,cloudShade:[1.14,1.04,1.04],bird:`#3a2630`},{hour:21.05,name:`twilight`,stops:[`#58678f`,`#6c6e9b`,`#8672a3`,`#a3689b`],glow:`#d07c98`,glowAmount:.55,horizon:.01,mist:`#9e8db6`,cloudSun:`#e6a2c8`,cloudShade:[1.24,1.14,1.26],bird:`#272238`},{hour:21.8,name:`dusk`,stops:[`#2f3659`,`#3f4571`,`#554f80`,`#5d4f7e`],glow:`#7a5c8e`,glowAmount:.3,horizon:0,mist:`#4c4c74`,cloudSun:`#5a4468`,cloudShade:[1.3,1.26,1.42],bird:`#1a1a2e`},{hour:23,...k}],j=e=>{let t=parseInt(e.slice(1),16);return[(t>>16&255)/255,(t>>8&255)/255,(t&255)/255]},M=e=>e<=.04045?e/12.92:((e+.055)/1.055)**2.4,N=e=>e<=.0031308?e*12.92:1.055*e**(1/2.4)-.055;function P(e){let[t,n,r]=e.map(M),i=Math.cbrt(.4122214708*t+.5363325363*n+.0514459929*r),a=Math.cbrt(.2119034982*t+.6806995451*n+.1073969566*r),o=Math.cbrt(.0883024619*t+.2817188376*n+.6299787005*r);return[.2104542553*i+.793617785*a-.0040720468*o,1.9779984951*i-2.428592205*a+.4505937099*o,.0259040371*i+.7827717662*a-.808675766*o]}function F([e,t,n]){let r=(e+.3963377774*t+.2158037573*n)**3,i=(e-.1055613458*t-.0638541728*n)**3,a=(e-.0894841775*t-1.291485548*n)**3;return[4.0767416621*r-3.3077115913*i+.2309699292*a,-1.2684380046*r+2.6097574011*i-.3413193965*a,-.0041960863*r-.7034186147*i+1.707614701*a].map(e=>Math.min(1,Math.max(0,N(Math.max(0,e)))))}function I(e,t,n){let r=Math.hypot(e[1],e[2]),i=Math.hypot(t[1],t[2]),a=Math.atan2(e[2],e[1]),o=Math.atan2(t[2],t[1]);r<.012&&(a=o),i<.012&&(o=a);let s=o-a;if(s>Math.PI&&(s-=2*Math.PI),s<-Math.PI&&(s+=2*Math.PI),B(a,s,L)&&!z(a,L)&&!z(o,L)){if(Math.min(r,i)>.035&&Math.abs(s)>135*Math.PI/180)s-=Math.sign(s)*2*Math.PI;else return V(e,t,n)}let c=e[0]+(t[0]-e[0])*n,l=r+(i-r)*n,u=a+s*n;return[c,l*Math.cos(u),l*Math.sin(u)]}var L=130*Math.PI/180,R=e=>Math.atan2(Math.sin(e),Math.cos(e)),z=(e,t)=>Math.abs(R(e-t))<35*Math.PI/180;function B(e,t,n){let r=R(n-e);return t>0?r>0&&r<t:r<0&&r>t}var V=(e,t,n)=>[e[0]+(t[0]-e[0])*n,e[1]+(t[1]-e[1])*n,e[2]+(t[2]-e[2])*n],H=new Map,U=e=>{let t=H.get(e);return t||H.set(e,t=P(j(e))),t},W=NaN,G;function K(e){return e!==W&&(W=e,G=ee(e)),G}function ee(e){let t=(e%24+24)%24,n=A.length,r=n-1;for(let e=0;e<n;e++)A[e].hour<=t&&(r=e);let i=A[r],a=A[(r+1)%n],o=(a.hour-i.hour+24)%24||24,s=((t-i.hour)%24+24)%24/o,c=s*s*(3-2*s),l=(e,t)=>I(U(e),U(t),c);return{stops:i.stops.map((e,t)=>l(e,a.stops[t])),glow:l(i.glow,a.glow),glowAmount:i.glowAmount+(a.glowAmount-i.glowAmount)*c,horizon:i.horizon+(a.horizon-i.horizon)*c,mist:l(i.mist,a.mist),cloudSun:V(j(i.cloudSun),j(a.cloudSun),c),cloudShade:V(i.cloudShade,a.cloudShade,c),bird:((e,t)=>F(l(e,t)))(i.bird,a.bird)}}var q=.5,te=1.15,ne=.55;function re(e){let t=14.5,n=((e-6)%24+24)%24,r,i;if(n<=t){let e=n/t;r=te*Math.sin(Math.PI*e),i=ne*(2*e-1)}else{let e=(n-t)/9.5;r=-.45*Math.sin(Math.PI*e),i=ne*(1-2*e)}let a=Math.cos(r);return[Math.sin(i)*a,Math.sin(r),Math.cos(i)*a]}function ie(e){let t=e[0],n=-(e[1]-q/2),r=Math.hypot(t,n)||1;return[t/r,n/r]}var ae=e=>Math.min(1,Math.max(0,(e[1]+.02)/.08)),oe=()=>{let e=new Date;return e.getHours()+e.getMinutes()/60},se=`#version 300 es
void main() {
  // one triangle that covers the viewport
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;function J(e,t,n=se){let r=e.createProgram();for(let[i,a]of[[e.VERTEX_SHADER,n],[e.FRAGMENT_SHADER,t]]){let t=e.createShader(i);if(e.shaderSource(t,a),e.compileShader(t),!e.getShaderParameter(t,e.COMPILE_STATUS)){let n=e.getShaderInfoLog(t)??`shader error`;throw console.error(n,`
`+a.split(`
`).map((e,t)=>`${t+1}: ${e}`).join(`
`)),Error(n)}e.attachShader(r,t)}if(e.linkProgram(r),!e.getProgramParameter(r,e.LINK_STATUS))throw Error(e.getProgramInfoLog(r)??`link error`);let i=new Map;return{program:r,u:t=>(i.has(t)||i.set(t,e.getUniformLocation(r,t)),i.get(t))}}function Y(e,t,n=e.CLAMP_TO_EDGE){let r=e.createTexture();return e.bindTexture(e.TEXTURE_2D,r),e.texParameteri(e.TEXTURE_2D,e.TEXTURE_MIN_FILTER,t),e.texParameteri(e.TEXTURE_2D,e.TEXTURE_MAG_FILTER,t),e.texParameteri(e.TEXTURE_2D,e.TEXTURE_WRAP_S,n),e.texParameteri(e.TEXTURE_2D,e.TEXTURE_WRAP_T,n),r}var ce=()=>({time:oe(),glow:.45,glowSize:.2,saturation:.8,soft:.2,brightness:1,sunlight:.48}),le=2.5/6,ue=1;function de(e){let t=e.time,n=((t-(t>=23||t<3?20:3))%24+24)%24,r=matchMedia(`(prefers-reduced-motion: reduce)`).matches,i=NaN;return a=>{if(!r){if(!Number.isNaN(i)&&e.time!==i)return void(r=!0);n=Math.max(0,n-a*le*Math.max(.1,Math.min(1,n/ue))),n<.01&&(n=0),i=e.time=((t-n)%24+24)%24,n===0&&(r=!0)}}}function fe(e){return[{label:`time`,min:0,max:24,step:.05,hint:`Hour of day; picks the palette and moves the sun. Defaults to your local time. ~5–7 sunrise, ~19–21 sunset.`,get:()=>e.time,set:t=>e.time=t},{label:`glow`,min:0,max:2,step:.01,hint:`Strength of the coloured glow around the sun. 0 = flat palette gradient.`,get:()=>e.glow,set:t=>e.glow=t},{label:`glow size`,min:.05,max:1,step:.01,hint:`How far the sun's glow spreads. Small = a tight halo in the corner, large = it washes the whole sky.`,get:()=>e.glowSize,set:t=>e.glowSize=t},{label:`saturation`,min:0,max:2,step:.01,hint:`Colour intensity of the palette. Lower = dustier, more muted.`,get:()=>e.saturation,set:t=>e.saturation=t},{label:`soft`,min:0,max:1,step:.01,hint:`Pastel haze: pulls toward the hour's own pale tint (peach at dawn, aqua at noon), not toward grey.`,get:()=>e.soft,set:t=>e.soft=t},{label:`brightness`,min:.5,max:1.5,step:.01,hint:`Overall lightness of the sky; hue is kept.`,get:()=>e.brightness,set:t=>e.brightness=t},{label:`sunlight`,min:0,max:2,step:.01,hint:`How strongly the sun tints the cloud tops facing it.`,get:()=>e.sunlight,set:t=>e.sunlight=t}]}var pe=e=>e.toFixed(4),me=`#version 300 es
precision highp float;
const float SPAN = ${pe(q)};
const float STOP_POS[4] = float[4](${O.map(pe).join(`, `)});

uniform vec2 uSize;
uniform float uAspect;
uniform vec3 uSun;          // unit vector toward the sun, y up, z into the screen
uniform vec3 uStops[4];     // OKLab gradient, horizon → zenith
uniform vec3 uGlowLab;      // OKLab glow colour
uniform float uGlowStrength;
uniform float uGlowSize;
uniform float uHorizon;
uniform vec3 uMist;         // OKLab
uniform float uSoft;
uniform float uSaturation;
uniform float uBrightness;
out vec4 outColor;

vec3 oklabToSrgb(vec3 c) {
  float l = c.x + 0.3963377774 * c.y + 0.2158037573 * c.z;
  float m = c.x - 0.1055613458 * c.y - 0.0638541728 * c.z;
  float s = c.x - 0.0894841775 * c.y - 1.2914855480 * c.z;
  l = l * l * l; m = m * m * m; s = s * s * s;
  vec3 lin = vec3(
    4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
   -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
   -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s);
  lin = max(lin, 0.0);
  vec3 lo = lin * 12.92;
  vec3 hi = 1.055 * pow(lin, vec3(1.0 / 2.4)) - 0.055;
  return clamp(mix(lo, hi, step(0.0031308, lin)), 0.0, 1.0);
}

// OKLab blend that doesn't go grey: where the two hues oppose (warm glow over cool sky) the
// midpoint is nudged toward rose, as real skies are between a peach glow and blue. Unlike a
// hue lerp this is continuous everywhere, so no seams where the hue difference crosses 180°.
vec3 mixWarm(vec3 p, vec3 q, float t) {
  vec3 c = mix(p, q, t);
  float cp = length(p.yz), cq = length(q.yz);
  float oppose = max(-dot(p.yz, q.yz), 0.0) / max(cp * cq, 1e-4);
  c.yz += vec2(0.9, 0.25) * min(cp, cq) * oppose * 4.0 * t * (1.0 - t);
  // then restore the chroma the straight line lost, capped so it stays continuous
  float lost = length(c.yz);
  c.yz *= min(mix(cp, cq, t) / max(lost, 1e-4), 2.0);
  return c;
}

vec3 gradient(float y) {
  vec3 c = uStops[0];
  for (int i = 0; i < 3; i++) {
    float t = clamp((y - STOP_POS[i]) / (STOP_POS[i + 1] - STOP_POS[i]), 0.0, 1.0);
    if (y >= STOP_POS[i]) c = mix(uStops[i], uStops[i + 1], t);
  }
  return c;
}

void main() {
  vec2 uv = gl_FragCoord.xy / uSize; // y up: 0 = horizon
  float elevation = uv.y * SPAN;
  float azimuth = (uv.x - 0.5) * SPAN * uAspect;
  vec3 dir = vec3(sin(azimuth) * cos(elevation), sin(elevation), cos(azimuth) * cos(elevation));

  vec3 col = gradient(uv.y);

  // horizon brightening: the long, hazy path near the ground
  col.x += uHorizon * exp(-uv.y / 0.12);

  // sun glow: a Mie-like forward lobe (tight core + broad skirt) around the sun's true
  // position, plus a band along the horizon under it, which keeps a corner lit while the
  // sun is just below the horizon
  float theta = acos(clamp(dot(dir, uSun), -1.0, 1.0));
  float core = exp(-theta / uGlowSize);
  float skirt = exp(-theta / (uGlowSize * 2.5));
  float sunAz = atan(uSun.x, uSun.z);
  float band = exp(-elevation / (0.35 * uGlowSize)) * exp(-abs(azimuth - sunAz) / (uGlowSize * 3.0));
  float g = uGlowStrength * (0.8 * core + 0.15 * skirt + 0.4 * band);
  float w = min(1.0 - exp(-1.4 * g), 0.95);
  col = mixWarm(col, uGlowLab, w);
  col.x += 0.05 * w * w;

  // pastel: lightness toward the palette's own mist, chroma relaxed but tinted by the mist
  col.x = mix(col.x, uMist.x, 0.6 * uSoft);
  col.yz = col.yz * (1.0 - 0.4 * uSoft) + uMist.yz * (0.3 * uSoft);

  col.yz *= uSaturation;
  col.x *= uBrightness;
  outColor = vec4(oklabToSrgb(col), 1.0);
}`,he=e=>[e.time,e.glow,e.glowSize,e.saturation,e.soft,e.brightness].join();function ge(e,t,n,r){let i=K(n.time);e.uniform1f(t(`uAspect`),r),e.uniform3f(t(`uSun`),...re(n.time)),e.uniform3fv(t(`uStops`),i.stops.flat()),e.uniform3f(t(`uGlowLab`),...i.glow),e.uniform1f(t(`uGlowStrength`),n.glow*i.glowAmount),e.uniform1f(t(`uGlowSize`),n.glowSize),e.uniform1f(t(`uHorizon`),i.horizon),e.uniform3f(t(`uMist`),...i.mist),e.uniform1f(t(`uSoft`),n.soft),e.uniform1f(t(`uSaturation`),n.saturation),e.uniform1f(t(`uBrightness`),n.brightness)}var _e=`
uniform vec3 uSunColour;
uniform vec3 uShadeColour;
uniform vec2 uSunScreen;
uniform vec3 uSunDir;
uniform float uSunUp;
uniform vec3 uBird;
`;function ve(e,t,n){let r=re(n.time),i=K(n.time),a=Math.min(1,Math.max(0,(r[1]+.15)/.2)),o=i.cloudSun.map(e=>e*n.sunlight*a),[s,c]=ie(r);e.uniform3f(t(`uSunColour`),...o),e.uniform3f(t(`uShadeColour`),...i.cloudShade),e.uniform2f(t(`uSunScreen`),s,c),e.uniform3f(t(`uSunDir`),...r),e.uniform1f(t(`uSunUp`),ae(r)),e.uniform3f(t(`uBird`),...i.bird)}var ye=()=>({shape:2,span:1,sweep:.5,beat:.8,stroke:1.1,spread:.35});function be(e){return[{label:`shape`,min:0,max:3,step:1,hint:`0 = silhouette from below with swept wings, 1 = gull: two thin painted strokes, 2 = side profile with a beating wing, 3 = half triangle flipping up and down.`,get:()=>e.shape,set:t=>e.shape=t},{label:`span`,min:.4,max:2.5,step:.05,hint:`Wing length relative to the body. 1 = as designed; the whole bird's size is birds → size.`,get:()=>e.span,set:t=>e.span=t},{label:`sweep`,min:0,max:1.5,step:.05,hint:`Silhouette: how far back the wings sweep (high = swift-like crescent). Gull: how deep the resting V is. Profile: how far back the wing reaches.`,get:()=>e.sweep,set:t=>e.sweep=t},{label:`beat depth`,min:0,max:1,step:.05,hint:`How far the wings swing on each beat. 0 = no visible flap.`,get:()=>e.beat,set:t=>e.beat=t},{label:`stroke`,min:.4,max:3,step:.05,hint:`Line weight of the gull shape, CSS px.`,get:()=>e.stroke,set:t=>e.stroke=t},{label:`depth spread`,min:0,max:.8,step:.05,hint:`How much smaller the farthest birds are drawn. 0 = all the same size.`,get:()=>e.spread,set:t=>e.spread=t}]}var xe=`#version 300 es
layout(location = 0) in vec4 aBird; // x, y (CSS px), unit "forward" of the shape
layout(location = 1) in vec3 aPose; // wing -1..1, bird scale (CSS px), x scale (side-on facing)
uniform vec2 uCss;
uniform float uPadCss;
uniform float uSpan;
out vec2 vLocal; // CSS px: x across the wings, y forward (side-on shapes: x forward, y up)
out float vWing;
out float vHalf;

void main() {
  // triangle strip over the quad's corners
  vec2 corner = vec2(float(gl_VertexID & 1), float(gl_VertexID >> 1)) * 2.0 - 1.0;
  vLocal = corner * (aPose.y * 1.5 * max(uSpan, 1.0) + uPadCss);
  vec2 h = aBird.zw;
  vec2 css = aBird.xy + vec2(-h.y, h.x) * (vLocal.x * aPose.z) + h * vLocal.y;
  gl_Position = vec4(css.x / uCss.x * 2.0 - 1.0, 1.0 - css.y / uCss.y * 2.0, 0.0, 1.0);
  vWing = aPose.x;
  vHalf = aPose.y;
}`,Se=`#version 300 es
precision highp float;
in vec2 vLocal;
in float vWing;
in float vHalf;
uniform float uCssPerPx; // CSS px per device px
uniform float uShape;
uniform float uSweep;
uniform float uSpan;    // wing length multiplier
uniform float uBeat;
uniform float uStroke;  // CSS px
uniform float uGrain;
uniform uint uFrame;
out vec4 outColor;
${D}
${_e}

/** Distance to a segment whose radius tapers from ra at a to rb at b. */
float taper(vec2 p, vec2 a, vec2 b, float ra, float rb) {
  vec2 ab = b - a;
  float t = clamp(dot(p - a, ab) / dot(ab, ab), 0.0, 1.0);
  return length(p - a - ab * t) - mix(ra, rb, t);
}

/** Side-on wing outline, convex, root front → tip → root back; x forward, y = fraction raised. */
const vec2 WING[9] = vec2[](
  vec2(0.49, 0.0), vec2(0.47, 0.27), vec2(0.38, 0.62), vec2(0.2, 1.0), // leading edge
  vec2(0.0, 0.92), vec2(-0.2, 0.77),                                   // tip, falling away
  vec2(-0.28, 0.63), vec2(-0.3, 0.35), vec2(-0.26, 0.0));              // trailing edge

/** Union of two distances, rounded over a band of width k so the join has no crease. */
float smin(float a, float b, float k) {
  float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return mix(b, a, h) - k * h * (1.0 - h);
}

/** Approximate distance to an axis-aligned ellipse with semi-axes r, centred on the origin. */
float ellipse(vec2 p, vec2 r) {
  float k0 = length(p / r);
  float k1 = length(p / (r * r));
  return k0 * (k0 - 1.0) / k1;
}

/** Signed distance to the triangle p0 p1 p2 (either winding). */
float triangle(vec2 p, vec2 p0, vec2 p1, vec2 p2) {
  vec2 e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
  vec2 v0 = p - p0, v1 = p - p1, v2 = p - p2;
  vec2 pq0 = v0 - e0 * clamp(dot(v0, e0) / dot(e0, e0), 0.0, 1.0);
  vec2 pq1 = v1 - e1 * clamp(dot(v1, e1) / dot(e1, e1), 0.0, 1.0);
  vec2 pq2 = v2 - e2 * clamp(dot(v2, e2) / dot(e2, e2), 0.0, 1.0);
  float s = sign(e0.x * e2.y - e0.y * e2.x);
  vec2 d = min(min(vec2(dot(pq0, pq0), s * (v0.x * e0.y - v0.y * e0.x)),
                   vec2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x))),
                   vec2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x)));
  return -sqrt(d.x) * sign(d.y);
}

void main() {
  // shape space: half the designed wingspan = 1 (wings themselves stretch by uSpan), y forward (up, for the gull); the shape is mirror-symmetric
  vec2 p = vLocal / vHalf;
  vec2 q = vec2(abs(p.x), p.y);
  float px = uCssPerPx / vHalf;
  float lift = vWing * uBeat * 1.15; // wing elevation, rad; + = up
  float reach = cos(lift);           // a tilted wing looks shorter from below
  float d;
  float fade = 1.0;

  if (uShape < 0.5) {
    // upstroke folds the wings back, downstroke reaches forward
    float fold = uSweep + max(lift, 0.0) * 0.4 - max(-lift, 0.0) * 0.2;
    vec2 shoulder = vec2(0.05, 0.08);
    vec2 wrist = vec2(0.04 + 0.4 * reach * uSpan, 0.12 - 0.16 * fold);
    vec2 tip = vec2(reach, -0.6 * fold) * uSpan;
    d = min(taper(q, shoulder, wrist, 0.12, 0.08), taper(q, wrist, tip, 0.08, 0.012));
    d = min(d, taper(q, vec2(0.0, 0.3), vec2(0.0, -0.32), 0.08, 0.035)); // head to tail
    d = min(d, taper(q, vec2(0.0, -0.28), vec2(0.07, -0.5), 0.035, 0.008)); // forked tail
  } else if (uShape > 2.5) {
    // x forward, y up: half a block bird, cut along its spine, as a wing that flips from above
    // the spine to below it on each beat (edge-on as it crosses); held raised in a glide
    float flip = clamp(0.3 + vWing * uBeat * 1.3, -1.0, 1.0) * 0.9 * uSpan;
    vec2 nose = vec2(0.9, 0.0), tail = vec2(-0.9, 0.0);
    d = abs(flip) < 0.01 ? taper(p, nose, tail, 0.0, 0.0) : triangle(p, nose, tail, vec2(-0.9, flip));
  } else if (uShape > 1.5) {
    // x forward, y up; proportions after a magpie in flight. The near wing is its outline
    // traced from the photo (WING: x forward, v = fraction of the way up the raised wing),
    // hinged along the back. Side-on its chord always shows; the beat only stretches it up or
    // down, and a glide holds it a little raised. Sweep leans it further back with height
    float rise = sin(lift + 0.3) * 1.3 * uSpan;
    float lean = (uSweep - 0.5) * 0.25;
    vec2 w[9];
    for (int i = 0; i < 9; i++) w[i] = vec2(WING[i].x - lean * WING[i].y, 0.14 + WING[i].y * rise);
    d = 1e9;
    for (int i = 1; i < 8; i++) d = min(d, triangle(p, w[0], w[i], w[i + 1]));
    d -= 0.03;
    // body from ellipses fitted to circles drawn over the photo: the back, plus a head whose
    // crown arcs all the way down to the beak tip with the chin cut out under it, plus a tail
    // squeezed between two cut-outs. Each cut-out shapes only its own part. Shape units here are
    // photo px × 0.00512, origin at (340, 300) in the 612×408 reference, facing right
    float back = ellipse(p - vec2(0.076, 0.04), vec2(0.378, 0.166));
    float head = ellipse(p - vec2(0.478, -0.158), vec2(0.455, 0.341));
    head = max(head, -ellipse(p - vec2(0.761, -0.275), vec2(0.497, 0.332))); // chin and underside of the beak
    head = max(head, -ellipse(p - vec2(0.102, -0.435), vec2(0.461, 0.471)));  // throat: trims the crown circle back to the belly
    float tail = ellipse(p - vec2(-0.461, 0.097), vec2(0.435, 0.184)); // tucked inside the back's underside
    tail = max(tail, -ellipse(p - vec2(-0.846, 0.261), vec2(0.797, 0.161))); // top of the tail
    tail = max(tail, -ellipse(p - vec2(-0.893, -0.274), vec2(0.789, 0.343))); // underside of the tail
    // a hair of blending where the parts meet, so the joins don't crease
    float body = smin(smin(back, head, 0.03), tail, 0.03);
    d = min(d, body);
  } else {
    // drawn upright (y = up on screen): wings up = a V, down = a ^; each wing sags between body and tip, the painted gull's two bowls
    float tipY = sin(lift) * 0.55 + uSweep * 0.3;
    vec2 centre = vec2(0.0, -0.04);
    vec2 tip = vec2(reach, tipY) * uSpan;
    vec2 bow = mix(centre, tip, 0.45) - vec2(0.0, 0.18 * uSpan);
    // never thinner than a device pixel; thinner strokes fade instead so they don't shimmer
    float r = uStroke * 0.5 / vHalf;
    float rr = max(r, px * 0.5);
    fade = r / rr;
    // the quadratic curve centre → bow → tip, as three segments
    vec2 b1 = mix(mix(centre, bow, 0.33), mix(bow, tip, 0.33), 0.33);
    vec2 b2 = mix(mix(centre, bow, 0.67), mix(bow, tip, 0.67), 0.67);
    d = min(min(taper(q, centre, b1, rr * 1.15, rr), taper(q, b1, b2, rr, rr * 0.85)), taper(q, b2, tip, rr * 0.85, rr * 0.55));
  }

  float a = clamp(0.5 - d / px, 0.0, 1.0) * fade;
  if (a <= 0.0) discard;
  // the same film grain the blit puts on the sky underneath
  float l = dot(uBird, vec3(0.2126, 0.7152, 0.0722));
  float g = grain(ivec2(gl_FragCoord.xy), uFrame) * uGrain * (0.4 + 2.4 * l * (1.0 - l));
  outColor = vec4(uBird + g * vec3(1.06, 1.0, 0.88), a);
}`,X=7,Ce=1.5,we=.5,Te=6,Ee=class{gl;look;program;vao;buffer;data=new Float32Array;facing=new Float32Array;constructor(e,t){this.gl=e,this.look=t,this.program=J(e,Se,xe),this.vao=e.createVertexArray(),this.buffer=e.createBuffer(),e.bindVertexArray(this.vao),e.bindBuffer(e.ARRAY_BUFFER,this.buffer),e.enableVertexAttribArray(0),e.vertexAttribPointer(0,4,e.FLOAT,!1,28,0),e.vertexAttribDivisor(0,1),e.enableVertexAttribArray(1),e.vertexAttribPointer(1,3,e.FLOAT,!1,28,16),e.vertexAttribDivisor(1,1),e.bindVertexArray(null)}draw(e,t,n,r,i,a,o,s,c){let{gl:l,look:u}=this,{count:d,x:f,y:p,vx:m,vy:h,wing:g,depth:_}=e;if(d===0)return;if(this.data.length<d*X&&(this.data=new Float32Array(d*X)),this.facing.length<d){let e=new Float32Array(d);e.set(this.facing),this.facing=e}let{data:v,facing:y}=this,b=n*Ce,x=Math.round(u.shape),S=1-Math.exp(-c*Te);for(let e=0,t=0;e<d;e++,t+=X){let n=Math.hypot(m[e],h[e])||1;v[t]=f[e],v[t+1]=p[e];let r=m[e]/n,i=h[e]/n,a=1;if(x===1){let e=.7*r*i;r=Math.sin(e),i=-Math.cos(e)}else if(x>=2){let t=y[e]===0?Math.sign(r)||1:r>.15?1:r<-.15?-1:Math.sign(y[e]);y[e]=y[e]===0?t:y[e]+(t-y[e])*S,a=y[e];let n=Math.max(-.5,Math.min(we,Math.atan2(-i,Math.abs(r))))*y[e];r=-Math.sin(n),i=-Math.cos(n)}v[t+2]=r,v[t+3]=i,v[t+4]=g[e],v[t+5]=b*(1-u.spread*_[e]),v[t+6]=a}let{program:C,u:w}=this.program,T=r/a;l.useProgram(C),l.bindVertexArray(this.vao),l.bindBuffer(l.ARRAY_BUFFER,this.buffer),l.bufferData(l.ARRAY_BUFFER,v.subarray(0,d*X),l.STREAM_DRAW),l.uniform2f(w(`uCss`),r,i),l.uniform1f(w(`uPadCss`),2*T),l.uniform1f(w(`uCssPerPx`),T),l.uniform1f(w(`uShape`),u.shape),l.uniform1f(w(`uSweep`),u.sweep),l.uniform1f(w(`uSpan`),u.span),l.uniform1f(w(`uBeat`),u.beat),l.uniform1f(w(`uStroke`),u.stroke),l.uniform1f(w(`uGrain`),o),l.uniform1ui(w(`uFrame`),s),ve(l,w,t),l.enable(l.BLEND),l.blendFunc(l.SRC_ALPHA,l.ONE_MINUS_SRC_ALPHA),l.drawArraysInstanced(l.TRIANGLE_STRIP,0,4,d),l.disable(l.BLEND),l.bindVertexArray(null)}},De=14,Oe=1.2,ke=.25,Ae=class{gl;w=0;h=0;dx=new Float32Array;dy=new Float32Array;tx=new Float32Array;ty=new Float32Array;gu=new Float32Array;gv=new Float32Array;rgba=new Float32Array;source=null;tex;maxGust=0;constructor(e){this.gl=e,this.tex=Y(e,e.LINEAR)}update(e,t,n){(e!==this.source||e.w!==this.w||e.h!==this.h)&&this.reset(e);let{w:r,h:i,gu:a,gv:o}=this,s=r+2,{u:c,v:l,windU:u,windV:d}=e,f=0;for(let e=0;e<c.length;e++){let t=c[e]-u,n=l[e]-d,r=Math.sqrt(t*t+n*n),i=r>Oe?1-Oe/r:0;r>f&&(f=r),a[e]=t*i,o[e]=n*i}if(this.maxGust=f,t>0){let{dx:e,dy:u,tx:d,ty:f}=this,p=Math.exp(-t/Math.max(n,.05));for(let n=1;n<=i;n++)for(let i=1;i<=r;i++){let r=i+s*n,m=i-.5-c[r]*t,h=n-.5-l[r]*t;d[r]=(this.sample(e,m,h)+a[r]*t)*p,f[r]=(this.sample(u,m,h)+o[r]*t)*p}Z(d,r,i),Z(f,r,i);for(let t=1;t<=i;t++)for(let n=1;n<=r;n++){let r=n+s*t,i=d[r]+ke*((d[r-1]+d[r+1]+d[r-s]+d[r+s])*.25-d[r]),a=f[r]+ke*((f[r-1]+f[r+1]+f[r-s]+f[r+s])*.25-f[r]),o=i*i+a*a;if(o>196){let e=De/Math.sqrt(o);i*=e,a*=e}e[r]=i,u[r]=a}Z(e,r,i),Z(u,r,i)}this.upload()}upload(){let{gl:e,dx:t,dy:n,gu:r,gv:i,rgba:a}=this;for(let e=0;e<t.length;e++)a[4*e]=t[e],a[4*e+1]=n[e],a[4*e+2]=r[e],a[4*e+3]=i[e];e.bindTexture(e.TEXTURE_2D,this.tex),e.texSubImage2D(e.TEXTURE_2D,0,0,0,this.w+2,this.h+2,e.RGBA,e.FLOAT,a)}reset(e){this.source=e,this.w=e.w,this.h=e.h;let t=(e.w+2)*(e.h+2);this.dx=new Float32Array(t),this.dy=new Float32Array(t),this.tx=new Float32Array(t),this.ty=new Float32Array(t),this.gu=new Float32Array(t),this.gv=new Float32Array(t),this.rgba=new Float32Array(t*4);let{gl:n}=this;n.bindTexture(n.TEXTURE_2D,this.tex),n.texImage2D(n.TEXTURE_2D,0,n.RGBA16F,e.w+2,e.h+2,0,n.RGBA,n.FLOAT,this.rgba)}sample(e,t,n){let{w:r,h:i}=this,a=je(t+.5,r),o=je(n+.5,i),s=a|0,c=o|0,l=a-s,u=o-c,d=r+2,f=e[s+d*c],p=e[s+1+d*c],m=e[s+d*(c+1)],h=e[s+1+d*(c+1)];return(f+(p-f)*l)*(1-u)+(m+(h-m)*l)*u}};function je(e,t){return e=(e-1)%t,e<0&&(e+=t),e+1}function Z(e,t,n){let r=t+2;for(let i=1;i<=n;i++)e[r*i]=e[t+r*i],e[t+1+r*i]=e[1+r*i];for(let i=0;i<=t+1;i++)e[i]=e[i+r*n],e[i+r*(n+1)]=e[i+r]}var Me=3.2,Ne=`
uniform vec2 uCloudScreen;      // canvas size, CSS px
uniform float uCloudTime;
uniform vec2 uCloudDrift;       // world drift: x = banks, y = cirrus
uniform sampler2D uCloudBend;   // fluid grid incl. wrapped border: rg = displacement, ba = gust (cells, cells/s)
uniform float uCloudCell;       // CSS px per fluid cell
uniform float uCloudBendGain;
uniform float uCloudCirrus;
uniform float uCloudFine;
uniform float uCloudCover;
uniform float uCloudVolume;
uniform float uCloudReach;
uniform float uCloudDensity;

const float CLOUD_VIEW = ${q.toFixed(4)};
const float BANK_BASE = ${1 .toFixed(3)};
const float CIRRUS_HEIGHT = ${Me.toFixed(3)};

/** Bend texture at a screen point: rg = displacement, ba = gust, in CSS px and CSS px/s. */
vec4 cloudFlow(vec2 css) {
  vec2 t = (css / uCloudCell + 1.0) / vec2(textureSize(uCloudBend, 0));
  return texture(uCloudBend, t) * uCloudCell;
}

/**
 * Fluid-driven displacement at a screen point, CSS px. Vertical pushes are damped: clouds lie
 * flat on their planes, so their noise is squashed vertically on screen and the same push
 * up or down tears through far more of a cloud's shape than one sideways.
 */
vec2 cloudBend(vec2 css) {
  return cloudFlow(css).rg * uCloudBendGain * vec2(1.0, 0.3);
}

/** View angles through a screen point, rad: x = azimuth (0 at centre), y = elevation (0 at the bottom edge). */
vec2 cloudAngles(vec2 css) {
  return vec2(css.x - 0.5 * uCloudScreen.x, uCloudScreen.y - css.y) * (CLOUD_VIEW / uCloudScreen.y);
}

/** View ray through a screen point: y up, z into the screen. Elevation is kept just above the horizon. */
vec3 cloudRay(vec2 css) {
  vec2 v = cloudAngles(css);
  float a = v.x, e = max(v.y, 0.012);
  return vec3(sin(a) * cos(e), sin(e), cos(a) * cos(e));
}

// ---- banks -------------------------------------------------------------------------
/** fBm over 2D simplex whose finer octaves churn faster, so the banks evolve as they drift. */
float bankNoise(vec2 p, int octaves) {
  float v = 0.0, a = 0.5;
  vec2 churn = vec2(0.021, -0.013) * uCloudTime;
  for (int i = 0; i < 6; i++) {
    if (i >= octaves) break;
    v += a * snoise(p + churn);
    p = OCTAVE_ROTATE * p * 2.07 + vec2(19.3, 7.9);
    churn = OCTAVE_ROTATE * churn * 1.6;
    a *= 0.5;
  }
  return v;
}

/** Bank density 0..1 at a point in bank noise space. */
float bankDensity(vec2 p, int octaves, float cover) {
  float n = bankNoise(p, octaves) * 0.5 + 0.5;
  float threshold = 1.0 - cover;
  return smoothstep(threshold, threshold + 0.22, n);
}

/** Point on the bank plane seen through a screen point; z eases off near the horizon so it never blows up. */
vec2 bankPlane(vec2 css) {
  vec2 v = cloudAngles(css);
  float z = BANK_BASE / (tan(max(v.y, 0.0)) + 0.18);
  return vec2(v.x * z, z);
}

/** Bank coverage at screen height e (0 = horizon, 1 = top): banks crowd toward the horizon. */
float bankCover(float e) {
  return uCloudCover * (1.0 - 0.75 * smoothstep(uCloudReach * 0.55, uCloudReach + 0.3, e));
}

/** Bank noise-space position for a bank-plane point: rolls stretched along x, with a slow warp so outlines wander. */
vec2 bankSpace(vec2 plane) {
  vec2 p = (plane + vec2(uCloudDrift.x, 0.0)) * vec2(1.25, 2.0);
  return p + vec2(snoise(p * 0.3 + vec2(3.1, 8.7)), snoise(p * 0.3 + vec2(11.4, 1.3))) * 0.35;
}

/** Distance fade of the banks into the sky (aerial perspective). */
float bankHaze(vec2 plane) {
  return exp(-plane.y * 0.12);
}

// ---- cirrus ------------------------------------------------------------------------
/** Thin line along the zero-contour of n, antialiased: sub-pixel lines fade rather than break up. */
float hairline(float n, float width) {
  float fw = fwidth(n);
  float w = max(width, fw * 1.1);
  return (1.0 - smoothstep(0.0, w, abs(n))) * (width / w);
}

/** Hairline density on the cirrus sheet. */
float cirrusThreads(vec2 p) {
  // long threads: zero-contours of strongly x-stretched simplex, two scales, gently warped
  vec2 q = p * vec2(0.2, 4.2);
  q.y += 0.7 * snoise(p * vec2(0.16, 0.8) + vec2(5.2, 1.7));
  float w = mix(0.2, 0.05, uCloudFine);
  float t1 = hairline(snoise(q), w);
  float t2 = hairline(snoise(OCTAVE_ROTATE * q * vec2(1.7, 2.1) + vec2(41.0, 13.0)), w * 0.7);
  return max(t1, t2 * 0.45);
}

/** Point on the cirrus sheet seen through a screen point, drifted. */
vec2 cirrusPlane(vec2 css) {
  vec3 dir = cloudRay(css);
  return dir.xz * (CIRRUS_HEIGHT / dir.y) + vec2(uCloudDrift.y, 0.0);
}

/** Where cirrus threads gather: a coarse fBm groups them into a few drifting fields, fading toward the horizon. */
float cirrusMask(vec2 p, float e) {
  float patchy = fbm(p * vec2(0.14, 0.6) + vec2(0.0, uCloudTime * 0.004), 3);
  return smoothstep(0.05, 0.5, patchy) * uCloudCirrus * smoothstep(0.1, 0.55, e);
}

/** Visible cloud amount at a screen point (both layers, bent by b = cloudBend(css)), for the vapour source. */
float cloudAmount(vec2 css, vec2 b) {
  float e = 1.0 - css.y / uCloudScreen.y;
  vec2 plane = bankPlane(css - b);
  float banks = bankDensity(bankSpace(plane), 4, bankCover(e)) * bankHaze(plane);
  vec2 cp = cirrusPlane(css - b * 0.6);
  float threads = cirrusMask(cp, e) * cirrusThreads(cp);
  return max(banks, threads * 0.6);
}
`;function Pe(e,t,n){let{look:r}=n;e.uniform2f(t(`uCloudScreen`),n.width,n.height),e.uniform1f(t(`uCloudTime`),n.time),e.uniform2f(t(`uCloudDrift`),n.drift[0],n.drift[1]),e.uniform1f(t(`uCloudCell`),n.cell),e.uniform1f(t(`uCloudBendGain`),r.bend),e.uniform1f(t(`uCloudCirrus`),r.cirrus),e.uniform1f(t(`uCloudFine`),r.fine),e.uniform1f(t(`uCloudCover`),r.cover),e.uniform1f(t(`uCloudVolume`),r.volume),e.uniform1f(t(`uCloudReach`),r.reach),e.uniform1f(t(`uCloudDensity`),r.density)}var Fe=3,Ie=`
const float VAPOUR_CSS = ${Fe.toFixed(1)};

float vapourAt(sampler2D vap, vec2 css) {
  return texture(vap, css / (VAPOUR_CSS * vec2(textureSize(vap, 0)))).r;
}
`,Le=`#version 300 es
precision highp float;
uniform sampler2D uVapPrev;
uniform vec2 uVapWind;      // steady wind, CSS px/s
uniform float uVapDt;
uniform float uVapDecay;    // per-frame factor
uniform float uVapGain;
out vec4 outColor;
${D}
${Ne}
${Ie}

void main() {
  vec2 css = gl_FragCoord.xy * VAPOUR_CSS;
  vec4 flow = cloudFlow(css);
  // carried by the whole flow (wind + gust), back-traced
  vec2 vel = flow.ba + uVapWind;
  float v = vapourAt(uVapPrev, css - vel * uVapDt) * uVapDecay;

  // shed where a real gust meets visible cloud
  float gust = length(flow.ba);
  float tear = smoothstep(${40 .toFixed(1)}, 220.0, gust);
  if (tear > 0.0 && uVapGain > 0.0) {
    float wisp = smoothstep(-0.15, 0.55, snoise(vec3(css * vec2(0.011, 0.035), uCloudTime * 1.3)));
    v += cloudAmount(css, flow.rg * uCloudBendGain) * tear * wisp * uVapGain * uVapDt;
  }
  outColor = vec4(min(v, 1.5), 0.0, 0.0, 1.0);
}`,Re=class{gl;program=null;tex=[];fbo=[];read=0;w=0;h=0;supported;liveFor=0;constructor(e){this.gl=e,this.supported=!!e.getExtension(`EXT_color_buffer_float`),this.supported&&(this.program=J(e,Le));for(let t=0;t<2;t++)this.tex.push(Y(e,e.LINEAR)),this.fbo.push(e.createFramebuffer())}get texture(){return this.tex[this.read]}resize(e,t){let{gl:n}=this,r=Math.max(1,Math.ceil(e/Fe)),i=Math.max(1,Math.ceil(t/Fe));if(r===this.w&&i===this.h)return;this.w=r,this.h=i;let a=new Float32Array(r*i);for(let e=0;e<2;e++)n.bindTexture(n.TEXTURE_2D,this.tex[e]),this.supported?(n.texImage2D(n.TEXTURE_2D,0,n.R16F,r,i,0,n.RED,n.FLOAT,a),n.bindFramebuffer(n.FRAMEBUFFER,this.fbo[e]),n.framebufferTexture2D(n.FRAMEBUFFER,n.COLOR_ATTACHMENT0,n.TEXTURE_2D,this.tex[e],0)):n.texImage2D(n.TEXTURE_2D,0,n.R8,r,i,0,n.RED,n.UNSIGNED_BYTE,new Uint8Array(r*i));n.bindFramebuffer(n.FRAMEBUFFER,null)}step(e,t,n,r,i,a,o,s){let{gl:c,program:l}=this;if(this.resize(e.width,e.height),!l||i<=0)return;if(s)this.liveFor=o*2.5;else if((this.liveFor-=i)<=0)return;let u=1-this.read;c.bindFramebuffer(c.FRAMEBUFFER,this.fbo[u]),c.viewport(0,0,this.w,this.h),c.useProgram(l.program),c.activeTexture(c.TEXTURE0),c.bindTexture(c.TEXTURE_2D,this.tex[this.read]),c.activeTexture(c.TEXTURE1),c.bindTexture(c.TEXTURE_2D,t);let d=l.u;c.uniform1i(d(`uVapPrev`),0),c.uniform1i(d(`uCloudBend`),1),Pe(c,d,e),c.uniform2f(d(`uVapWind`),n,r),c.uniform1f(d(`uVapDt`),i),c.uniform1f(d(`uVapDecay`),Math.exp(-3*i/Math.max(o,.05))),c.uniform1f(d(`uVapGain`),a),c.drawArrays(c.TRIANGLES,0,3),c.bindFramebuffer(c.FRAMEBUFFER,null),c.activeTexture(c.TEXTURE0),this.read=u}},ze=()=>({cirrus:.8,fine:.5,cover:.48,volume:.38,reach:.65,density:2.6,drift:1,bend:.35,heal:2.5,puff:.25,puffLife:.9});function Be(e){return[{label:`cirrus`,min:0,max:2,step:.01,hint:`Thin high threads in the upper sky. 0 = none.`,get:()=>e.cirrus,set:t=>e.cirrus=t},{label:`fine`,min:0,max:1,step:.01,hint:`Cirrus thread width. Higher = finer hairlines.`,get:()=>e.fine,set:t=>e.fine=t},{label:`cover`,min:0,max:1,step:.01,hint:`How much of the sky the lower banks cover. 0 = clear.`,get:()=>e.cover,set:t=>e.cover=t},{label:`volume`,min:.1,max:1,step:.01,hint:`How deep the banks look. Higher = brighter tops, darker undersides and cores.`,get:()=>e.volume,set:t=>e.volume=t},{label:`reach`,min:0,max:1,step:.01,hint:`How far up the screen the banks climb from the horizon.`,get:()=>e.reach,set:t=>e.reach=t},{label:`density`,min:.2,max:3,step:.01,hint:`How solid the banks are. Higher = opaque with darker cores, lower = misty.`,get:()=>e.density,set:t=>e.density=t},{label:`drift`,min:0,max:4,step:.05,hint:`How fast the wind carries the clouds. Lower layers move faster (parallax).`,get:()=>e.drift,set:t=>e.drift=t},{label:`bend`,min:0,max:3,step:.05,hint:`How much the cursor's gusts push and swirl the clouds. 0 = untouchable.`,get:()=>e.bend,set:t=>e.bend=t},{label:`heal`,min:.5,max:10,step:.1,hint:`Seconds for bent clouds to drift back into shape.`,get:()=>e.heal,set:t=>e.heal=t},{label:`puff`,min:0,max:3,step:.05,hint:`Vapour torn off a cloud when the cursor sweeps through it. 0 = none.`,get:()=>e.puff,set:t=>e.puff=t},{label:`puff life`,min:.2,max:3,step:.05,hint:`Seconds for torn-off vapour to dissolve.`,get:()=>e.puffLife,set:t=>e.puffLife=t}]}var Ve=`
${Ne}
${Ie}
uniform sampler2D uCloudVapour;  // VAPOUR_CSS CSS px per texel, row 0 at the top

float cloudLuma(vec3 c) {
  return dot(c, vec3(0.2126, 0.7152, 0.0722));
}

/** The local sky, desaturated a little and brightened by k: the ambient light on a cloud top. */
vec3 cloudLift(vec3 sky, float k) {
  return mix(sky, vec3(cloudLuma(sky)), 0.25) * k;
}

/**
 * Clouds may shade, but never read as holes in the sky: luminance is held to within a few
 * percent of the sky behind (a little more room by day, when bases really are darker).
 */
vec3 cloudFloor(vec3 col, vec3 sky) {
  float floorL = cloudLuma(sky) * mix(0.95, 0.88, uSunUp);
  float l = cloudLuma(col);
  return l < floorL ? col * (floorL / max(l, 1e-4)) : col;
}

/**
 * Cloud body colour, shared by banks and vapour so torn wisps match the cloud they came from.
 * top: 0 at the shaded base, 1 at the lit top. sun: how much direct sunlight reaches here.
 * Ambient is the local sky, lifted toward the top; brightening beyond the sky rolls off per
 * channel into the headroom the sky leaves, so bright tops never clip into flat white.
 */
vec3 cloudColour(vec3 sky, float top, float sun) {
  vec3 ambient = mix(sky * uShadeColour * 0.9, cloudLift(sky, 1.08), mix(0.2, 1.0, top));
  vec3 lit = ambient + uSunColour * sun;
  vec3 room = max(1.0 - sky, vec3(0.1));
  vec3 over = max(lit - sky, 0.0);
  lit = min(lit, sky) + room * (1.0 - exp(-over / room));
  return cloudFloor(lit, sky);
}

/**
 * 2.5D banks: density on the bank plane, light marched a few steps toward the sun through the
 * same field. Shading is built so bodies stay bright: ambient sky light lifted toward the top,
 * a multi-scatter term that grows with thickness (cores glow rather than die), darkening kept
 * to the sun-away base, and a floor so nothing drops far below the sky. Forward scattering
 * gives thin edges between us and the sun a warm rim. Returns premultiplied colour and coverage.
 */
vec4 banks(vec3 sky, vec2 css, float e) {
  float cover = bankCover(e);
  if (cover < 0.01) return vec4(0.0);

  vec2 plane = bankPlane(css);
  vec2 p = bankSpace(plane);
  float d = bankDensity(p, 5, cover);
  if (d < 0.004) return vec4(0.0);

  // sun on the bank plane (x right, +y = farther = up the screen). Height pulls the light
  // toward "from above", so by day the near (lower) side of each bank is its base.
  vec3 sun = uSunDir;
  vec2 L = normalize(vec2(sun.x, sun.z * 0.35 + max(sun.y, 0.0) * 1.6) + 1e-4);
  float sigma = 4.0 * uCloudDensity * (0.5 + uCloudVolume);
  float od = 0.0;
  const float STEP = 0.07;
  for (int i = 1; i <= 4; i++) {
    od += bankDensity(p + L * (STEP * 0.6 * float(i * i)), 3, cover) * float(2 * i - 1);
  }
  od *= sigma * STEP * 0.6;
  float direct = exp(-od);
  // multiple scattering: thick bodies fill with soft light that fades slowly with depth
  float scatter = (1.0 - exp(-2.5 * d)) * exp(-od * 0.18) * 0.55;
  // 0 at the base (deep toward the sun), 1 at the lit top
  float top = clamp(exp(-od * 0.35) + 0.3 * (1.0 - d), 0.0, 1.0);

  // forward scattering (Henyey-Greenstein) for edges seen against the sun
  vec3 ray = cloudRay(css);
  float mu = dot(ray, sun);
  const float G = 0.6;
  float hg = (1.0 - G * G) / pow(1.0 + G * G - 2.0 * G * mu, 1.5) / (4.0 * 3.14159);
  float rim = hg * (1.0 - d) * 1.6;

  vec3 lit = cloudColour(sky, top, direct * 0.85 + scatter + top * 0.3 + rim);

  float a = d * mix(0.55, 0.95, smoothstep(0.0, 1.0, uCloudDensity)) * bankHaze(plane);
  return vec4(lit * a, a);
}

vec3 cirrus(vec3 sky, vec2 css, vec2 bend, float e) {
  vec2 bent = css - bend * 0.6;
  vec2 p = cirrusPlane(bent);
  float mask = cirrusMask(p, e);
  if (mask < 0.003) return sky;

  float d0 = cirrusThreads(p);
  if (d0 < 0.003) return sky;
  // the same threads a step toward the sun: lit on the sun side, shadowed on the far side
  vec2 toSun = uSunScreen * max(fwidth(css.y) * 1.5, 2.5);
  float d1 = cirrusThreads(cirrusPlane(bent + toSun));
  float light = clamp(0.6 + (d0 - d1) * 1.4, 0.0, 1.0);

  vec3 shade = sky * uShadeColour * 0.9;
  vec3 lit = sky * 0.85 + uSunColour * 0.7;
  return mix(sky, cloudFloor(mix(shade, lit, light), sky), clamp(d0 * mask, 0.0, 1.0) * 0.7);
}

/**
 * Torn-off wisps, shaded like thin bank cloud: the vapour a step toward the sun stands in for
 * the banks' light march, so the sun side of each wisp is lit and the far side shaded.
 */
vec3 vapour(vec3 col, vec3 sky, vec2 css) {
  float v = vapourAt(uCloudVapour, css);
  if (v < 0.004) return col;
  float ahead = vapourAt(uCloudVapour, css + uSunScreen * 7.0);
  float od = 4.0 * uCloudDensity * (0.5 + uCloudVolume) * (ahead + v * 0.5) * 0.15;
  float top = clamp(exp(-od * 0.35) + 0.3 * (1.0 - v), 0.0, 1.0);
  vec3 body = cloudColour(sky, top, exp(-od) * 0.85 * 0.6 + top * 0.3);
  return mix(col, body, (1.0 - exp(-v * 2.4)) * 0.8);
}

vec3 composeClouds(vec3 sky, vec2 css, vec2 uv) {
  vec2 bend = cloudBend(css);
  vec3 col = cirrus(sky, css, bend, uv.y);
  vec4 b = banks(sky, css - bend, uv.y);
  col = col * (1.0 - b.a) + b.rgb;
  return vapour(col, sky, css);
}
`,He=class{look;bend;vapour;shape;constructor(e,t){this.look=t,this.bend=new Ae(e),this.vapour=new Re(e),this.shape={width:1,height:1,time:0,drift:[0,0],cell:1,look:t}}update(e,t,n,r,i){let{shape:a,look:o}=this;a.time+=n,a.width=r,a.height=i,a.cell=t;let s=e.windU*t*n*o.drift*(q/Math.max(i,1));a.drift[0]+=s*1*4.5,a.drift[1]+=s*Me*2.2,this.bend.update(e,n,o.heal);let c=o.puff>0&&this.bend.maxGust*t>40;this.vapour.step(a,this.bend.tex,e.windU*t,e.windV*t,n,o.puff*7,o.puffLife,c)}setUniforms(e,t,n){Pe(e,t,this.shape),e.activeTexture(e.TEXTURE0+n),e.bindTexture(e.TEXTURE_2D,this.bend.tex),e.uniform1i(t(`uCloudBend`),n),e.activeTexture(e.TEXTURE0+n+1),e.bindTexture(e.TEXTURE_2D,this.vapour.texture),e.uniform1i(t(`uCloudVapour`),n+1),e.activeTexture(e.TEXTURE0)}},Ue=()=>({dither:0,levels:16,grain:.075});function We(e){return[{label:`dither`,min:0,max:1,step:1,hint:`0 = smooth gradients with grain, 1 = Bayer 8×8 ordered dither in 3px dots (brightness only, hue kept).`,get:()=>e.dither,set:t=>e.dither=t},{label:`levels`,min:2,max:16,step:1,hint:`Brightness steps when dithering; hue is kept. 2 = stark, 8+ = nearly smooth.`,get:()=>e.levels,set:t=>e.levels=t},{label:`grain`,min:0,max:.2,step:.005,hint:`Warm film grain over everything, as it looks on a 2× screen; scaled so it reads the same on sharper and coarser ones. Smooth mode only.`,get:()=>e.grain,set:t=>e.grain=t}]}var Ge=()=>({palette:ce(),clouds:ze(),birds:ye(),finish:Ue()}),Ke=e=>[{title:`palette`,knobs:fe(e.palette)},{title:`sky clouds`,knobs:Be(e.clouds)},{title:`bird shape`,knobs:be(e.birds)},{title:`finish`,knobs:We(e.finish)}],qe=`#version 300 es
precision highp float;
uniform sampler2D uSky;
uniform sampler2D uBirds;
uniform vec2 uTarget;   // size of this pass in pixels
uniform vec2 uCssPerPx; // CSS px covered by one pixel of this pass
uniform vec2 uCss;      // canvas size in CSS px
uniform vec2 uDotCss;   // CSS px per bird-mask dot
uniform float uDither;
uniform float uLevels;
out vec4 outColor;
${D}
${_e}
${Ve}

float bayer8(ivec2 p) {
  int x = p.x & 7, y = p.y & 7, xy = x ^ y;
  int v = ((xy & 1) << 5) | ((y & 1) << 4) | ((xy & 2) << 2) | ((y & 2) << 1) | ((xy & 4) >> 1) | ((y & 4) >> 2);
  return (float(v) + 0.5) / 64.0;
}

void main() {
  vec2 px = vec2(gl_FragCoord.x, uTarget.y - gl_FragCoord.y); // y down
  vec2 css = px * uCssPerPx;
  vec2 uv = vec2(css.x / uCss.x, 1.0 - css.y / uCss.y);

  vec3 col = composeClouds(texture(uSky, uv).rgb, css, uv);

  ivec2 birdDot = ivec2(floor(css / uDotCss));
  if (texelFetch(uBirds, birdDot, 0).r > 0.5) col = uBird;

  if (uDither > 0.5) {
    // dither brightness only and keep the local hue: per-channel rounding snaps to a fixed
    // RGB cube and sprays unrelated hues (red, yellow) into a soft gradient
    float steps = max(uLevels - 1.0, 1.0);
    float l = max(dot(col, vec3(0.2126, 0.7152, 0.0722)), 1e-3);
    float lq = floor(l * steps + bayer8(ivec2(px))) / steps;
    col = clamp(col * (lq / l), 0.0, 1.0);
  }
  outColor = vec4(col, 1.0);
}`,Je=`#version 300 es
precision highp float;
uniform sampler2D uScene;
uniform vec2 uDevice;
uniform vec2 uTarget;
uniform float uDotPx;
uniform float uDither;
uniform float uGrain;
uniform uint uFrame;
out vec4 outColor;
${D}

void main() {
  if (uDither > 0.5) {
    // whole device pixels per dot, so the grid stays exact at any pixel ratio
    ivec2 d = ivec2(floor(vec2(gl_FragCoord.x, uDevice.y - gl_FragCoord.y) / uDotPx));
    outColor = vec4(texelFetch(uScene, ivec2(d.x, int(uTarget.y) - 1 - d.y), 0).rgb, 1.0);
    return;
  }
  vec3 col = texture(uScene, gl_FragCoord.xy / uDevice).rgb;
  // grain is strongest in the midtones, like film, and slightly warm
  float l = dot(col, vec3(0.2126, 0.7152, 0.0722));
  float g = grain(ivec2(gl_FragCoord.xy), uFrame) * uGrain * (0.4 + 2.4 * l * (1.0 - l));
  outColor = vec4(col + g * vec3(1.06, 1.0, 0.88), 1.0);
}`,Q=320,Ye=24,Xe=2;function Ze(e,t,n,r=3){let i=e.getContext(`webgl2`,{antialias:!1,alpha:!1});if(!i)return null;let{palette:a,finish:o}=n,s=de(a),c=J(i,me),l=J(i,qe),u=J(i,Je);i.bindVertexArray(i.createVertexArray()),i.pixelStorei(i.UNPACK_ALIGNMENT,1);let d=new He(i,n.clouds),f=new Ee(i,n.birds),p=Y(i,i.LINEAR),m=i.createFramebuffer(),h=Y(i,i.LINEAR),g=i.createFramebuffer(),_=Y(i,i.NEAREST),v=0,y=0,b=0,x=0,S=r,C=r,T=r,E=0,D=0,O=new Uint8Array,k=!1,A=0,j=``,M=``,N=0,P=0,F=0;function I(e,t,n,r){i.bindTexture(i.TEXTURE_2D,t),i.texImage2D(i.TEXTURE_2D,0,i.RGBA8,n,r,0,i.RGBA,i.UNSIGNED_BYTE,null),i.bindFramebuffer(i.FRAMEBUFFER,e),i.framebufferTexture2D(i.FRAMEBUFFER,i.COLOR_ATTACHMENT0,i.TEXTURE_2D,t,0),i.bindFramebuffer(i.FRAMEBUFFER,null)}function L(){let e=o.dither>.5,t=`${e}:${v}x${y}:${E}x${D}`;t!==M&&(M=t,N=e?E:Math.ceil(v),P=e?D:Math.ceil(y),I(g,h,N,P))}function R(){let e=`${he(a)}:${A}:${v/y}`;if(e===j)return;j=e;let{program:t,u:n}=c;i.useProgram(t),i.bindFramebuffer(i.FRAMEBUFFER,m),i.viewport(0,0,Q,A),i.uniform2f(n(`uSize`),Q,A),ge(i,n,a,v/y),i.drawArrays(i.TRIANGLES,0,3),z(),i.bindFramebuffer(i.FRAMEBUFFER,null)}function z(){let e=new Uint8Array(4);i.readPixels(160,A-1,1,1,i.RGBA,i.UNSIGNED_BYTE,e),document.documentElement.style.setProperty(`--sky-top`,`rgb(${e[0]} ${e[1]} ${e[2]})`)}return{resize(e,t,n,a){v=e,y=t,b=n,x=a,S=Math.max(1,Math.round(n/e*r)),C=S*e/n,T=S*t/a,E=Math.ceil(n/S),D=Math.ceil(a/S),O=new Uint8Array(E*D),i.bindTexture(i.TEXTURE_2D,_),i.texImage2D(i.TEXTURE_2D,0,i.R8,E,D,0,i.RED,i.UNSIGNED_BYTE,O),A=Math.max(2,Math.round(Q*t/e)),I(m,p,Q,A),j=``,M=``},setColors(){},draw(e,n,r,c){F+=c,s(c),R(),L(),d.update(e,n,c,v,y);let m=o.dither>.5;(m||k)&&(O.fill(0),m&&w(r,O,E,D,C,T,t.birdSize,255),i.bindTexture(i.TEXTURE_2D,_),i.texSubImage2D(i.TEXTURE_2D,0,0,0,E,D,i.RED,i.UNSIGNED_BYTE,O),k=m);let A=l;i.useProgram(A.program),i.bindFramebuffer(i.FRAMEBUFFER,g),i.viewport(0,0,N,P),i.activeTexture(i.TEXTURE0),i.bindTexture(i.TEXTURE_2D,p),i.activeTexture(i.TEXTURE1),i.bindTexture(i.TEXTURE_2D,_),i.uniform1i(A.u(`uSky`),0),i.uniform1i(A.u(`uBirds`),1),i.uniform2f(A.u(`uTarget`),N,P),o.dither>.5?i.uniform2f(A.u(`uCssPerPx`),C,T):i.uniform2f(A.u(`uCssPerPx`),v/N,y/P),i.uniform2f(A.u(`uCss`),v,y),i.uniform2f(A.u(`uDotCss`),C,T),i.uniform1f(A.u(`uDither`),o.dither),i.uniform1f(A.u(`uLevels`),o.levels),ve(i,A.u,a),d.setUniforms(i,A.u,2),i.drawArrays(i.TRIANGLES,0,3);let j=u;i.useProgram(j.program),i.bindFramebuffer(i.FRAMEBUFFER,null),i.viewport(0,0,b,x),i.activeTexture(i.TEXTURE0),i.bindTexture(i.TEXTURE_2D,h),i.uniform1i(j.u(`uScene`),0),i.uniform2f(j.u(`uDevice`),b,x),i.uniform2f(j.u(`uTarget`),N,P),i.uniform1f(j.u(`uDotPx`),S),i.uniform1f(j.u(`uDither`),o.dither);let M=o.grain*(b/v)/Xe,I=Math.floor(F*Ye);i.uniform1f(j.u(`uGrain`),M),i.uniform1ui(j.u(`uFrame`),I),i.drawArrays(i.TRIANGLES,0,3),m||f.draw(r,a,t.birdSize,v,y,b,M,I,c)}}}function Qe(e){let t=document.createElement(`aside`);t.className=`knobs`,t.setAttribute(`aria-label`,`Tuning`);let n=document.createElement(`div`);n.className=`knobs-tip`,n.id=`knobs-tip`,n.setAttribute(`role`,`tooltip`),n.hidden=!0;let r=new $e(n);for(let n of e)t.append(et(n,r));return document.body.append(t,n),()=>{t.remove(),n.remove()}}var $e=class{tip;constructor(e){this.tip=e}attach(e,t,n){if(!n)return;let r=()=>this.show(e,n),i=()=>this.hide();e.addEventListener(`pointerenter`,r),e.addEventListener(`pointerleave`,i),t.addEventListener(`focus`,r),t.addEventListener(`blur`,i),t.setAttribute(`aria-describedby`,this.tip.id)}show(e,t){let{tip:n}=this;n.textContent=t,n.hidden=!1;let r=e.getBoundingClientRect(),i=e.closest(`.knobs`).getBoundingClientRect(),a=i.left-n.offsetWidth-8;a<8&&(a=Math.min(i.right+8,innerWidth-n.offsetWidth-8));let o=Math.max(8,Math.min(r.top+r.height/2-n.offsetHeight/2,innerHeight-n.offsetHeight-8));n.style.left=`${a}px`,n.style.top=`${o}px`}hide(){this.tip.hidden=!0}};function et({title:e,knobs:t,actions:n=[]},r){let i=document.createElement(`section`);i.setAttribute(`aria-label`,e);let a=document.createElement(`p`);a.className=`knobs-title`,a.textContent=e,i.append(a);for(let e of t){let t=document.createElement(`label`),n=document.createElement(`span`);n.textContent=e.label;let a=document.createElement(`input`);a.type=`range`,a.min=String(e.min),a.max=String(e.max),a.step=String(e.step),a.value=String(e.get());let o=document.createElement(`output`),s=Math.max(0,-Math.floor(Math.log10(e.step))),c=()=>o.textContent=e.get().toFixed(s);a.addEventListener(`input`,()=>{e.set(Number(a.value)),c()}),c(),t.append(n,a,o),r.attach(t,a,e.hint),i.append(t)}let o=document.createElement(`div`);o.className=`knobs-actions`;let s=document.createElement(`button`);s.type=`button`,s.textContent=`copy`,s.addEventListener(`click`,async()=>{let e=Object.fromEntries(t.map(e=>[e.label,e.get()]));await navigator.clipboard.writeText(JSON.stringify(e,null,2)),s.textContent=`copied`,setTimeout(()=>s.textContent=`copy`,1200)}),r.attach(s,s,`Copy the ${e} values as JSON, to paste back as new defaults.`),o.append(s);for(let e of n){let t=document.createElement(`button`);t.type=`button`,t.textContent=e.label,t.addEventListener(`click`,e.run),r.attach(t,t,e.hint),o.append(t)}return i.append(o),i}var tt={"void-and-cluster":()=>h(64),bayer2:()=>p(2),bayer4:()=>p(4),bayer8:()=>p(8)},nt=9,rt=26e3,it=.8,at=.07,ot=.12,st=3.5,ct=.35,lt=320,ut=1.3,dt={wake:.01,gustGain:1,gustMaxAccel:260,gustBoost:.9},ft=(e,t)=>Math.round(Math.min(150,Math.max(20,e*t*76/1405320)));function pt(t,{mode:n=`dither`,matrix:i=`bayer8`,knobs:a=!1}={}){let s={...d},c={...y},l={wind:it,dissolve:ot},u={...r},p={...dt},m=null,h=Ge(),g=n===`sky`?Ze(t,c,h):null;g&&(document.documentElement.dataset.mode=`sky`);let _=g??(n===`ascii`?E(t,c):T(t,tt[i](),c)),v=new o(u),b,x,S=nt,C=0,w=0,D=matchMedia(`(prefers-reduced-motion: reduce)`),O=matchMedia(`(prefers-color-scheme: dark)`);function k(e,t){let n=e.getPropertyValue(t).trim().split(/[\s,]+/).map(Number);return n.length===3&&n.every(Number.isFinite)?n:null}function A(){let e=getComputedStyle(document.documentElement),t=k(e,`--fg`),n=k(e,`--cloud`);t&&n&&_.setColors(t,n)}function j(){b.windU=l.wind,b.windV=l.wind*at}function M(t=b.w,n=b.h){b=new e(t,n),b.setWind(l.wind,l.wind*at),b.dyeDecay=l.dissolve,x=new f(b,s)}function N(e,n,r,i){t.width=r,t.height=i,_.resize(e,n,r,i),(e!==C||n!==w)&&(C=e,w=n,S=Math.max(nt,Math.sqrt(C*w/rt)),M(Math.ceil(C/S),Math.ceil(w/S)),v.resize(C,w,m??ft(C,w)))}let P={x:0,y:0,lastX:0,lastY:0,active:!1,moved:!1};function F(e){if(e.target instanceof Element&&e.target.closest(`.knobs`)){P.active=!1;return}let n=t.getBoundingClientRect(),r=e.clientX-n.left,i=e.clientY-n.top;P.active||=(P.lastX=r,P.lastY=i,!0),P.x=r,P.y=i,P.moved=!0}function I(){P.active=!1}function L(e){e.pointerType!==`mouse`&&(P.active=!1)}function R(e){if(!P.moved)return;P.moved=!1;let t=P.x-P.lastX,n=P.y-P.lastY;P.lastX=P.x,P.lastY=P.y;let r=Math.hypot(t,n);if(r<.5)return;let i=Math.min(r/e,lt),a=t/r*(i/S),o=n/r*(i/S),s=Math.max(1,Math.ceil(r/S/1.5));for(let e=1;e<=s;e++){let r=e/s,i=P.x-t*(1-r),c=P.y-n*(1-r);b.stir(i/S,c/S,a,o,ct,st)}}let z=new Float32Array(2);function B(){let{count:e,x:t,y:n,vx:r,vy:i,ax:a,ay:o,boost:s}=v;for(let c=0;c<e;c++){if(t[c]<0||t[c]>C||n[c]<0||n[c]>w)continue;let e=t[c]/S,l=n[c]/S;b.splat(e,l,r[c]/S*p.wake,i[c]/S*p.wake,0,ut),b.sample(e,l,z);let u=(z[0]-b.windU)*S,d=(z[1]-b.windV)*S,f=Math.hypot(u,d);if(f<.001)continue;let m=Math.min(f*p.gustGain,p.gustMaxAccel);a[c]+=u/f*m,o[c]+=d/f*m;let h=f*p.gustBoost;h>s[c]&&(s[c]=h)}}let V=performance.now(),H=0;function U(e){let t=D.matches?.35:1,n=Math.min((e-V)/1e3,1/30)*t;V=e,n>0&&(R(n/t),B(),g||x.step(n),b.step(n),v.step(n)),_.draw(b,S,v,n),H=requestAnimationFrame(U)}let W=new ResizeObserver(([e])=>{let{width:t,height:n}=e.contentRect,r=e.devicePixelContentBoxSize?.[0],i=window.devicePixelRatio||1;N(t,n,r?r.inlineSize:Math.round(t*i),r?r.blockSize:Math.round(n*i))});try{W.observe(t,{box:`device-pixel-content-box`})}catch{W.observe(t)}A();let G=window.devicePixelRatio||1;N(t.clientWidth,t.clientHeight,Math.round(t.clientWidth*G),Math.round(t.clientHeight*G)),window.addEventListener(`pointermove`,F,{passive:!0}),window.addEventListener(`pointerdown`,F,{passive:!0}),window.addEventListener(`pointerup`,L,{passive:!0}),document.documentElement.addEventListener(`pointerleave`,I),window.addEventListener(`blur`,I),O.addEventListener(`change`,A),H=requestAnimationFrame(U);let K=a?Qe([...g?[...Ke(h),{title:`air`,knobs:[{label:`wind`,min:0,max:6,step:.1,hint:`Steady breeze that carries the whole sky sideways.`,get:()=>l.wind,set:e=>(l.wind=e,j())}]}]:[{title:`clouds`,knobs:[{label:`amount`,min:0,max:.4,step:.005,hint:`How much vapour each cloud puts out. Higher = thicker clouds. Builds up over ~8s.`,get:()=>s.amount,set:e=>s.amount=e},{label:`density`,min:0,max:3,step:.05,hint:`How many clouds are in the sky at once. New ones fade in; extras wait to die off.`,get:()=>s.density,set:e=>s.density=e},{label:`size`,min:1,max:20,step:.5,hint:`Typical cloud size. Applies to newly formed clouds — hit 'new sky' to see it now.`,get:()=>s.size,set:e=>s.size=e},{label:`stretch`,min:.5,max:8,step:.1,hint:`Width ÷ height of a cloud. Higher = long flat streaks. New clouds only.`,get:()=>s.stretch,set:e=>s.stretch=e},{label:`puffs`,min:1,max:16,step:1,hint:`Blobs per cloud. More = lumpier, fuller shapes. New clouds only.`,get:()=>s.puffs,set:e=>s.puffs=e},{label:`life`,min:4,max:90,step:1,hint:`Seconds a cloud lives: forms, holds, then thins out. New clouds only.`,get:()=>s.life,set:e=>s.life=e},{label:`dissolve`,min:0,max:1,step:.01,hint:`How fast vapour fades once released. Higher = crisper clouds, shorter trails.`,get:()=>l.dissolve,set:e=>l.dissolve=b.dyeDecay=e},{label:`wind`,min:0,max:6,step:.1,hint:`Steady breeze that carries the whole sky sideways.`,get:()=>l.wind,set:e=>(l.wind=e,j())},{label:`opacity`,min:0,max:1,step:.01,hint:`How filled-in the thickest part of a cloud is drawn. Instant.`,get:()=>c.opacity,set:e=>c.opacity=e},{label:`thickness`,min:.2,max:5,step:.1,hint:`How quickly thin vapour reads as solid. Higher = puffier, denser look. Instant.`,get:()=>c.thickness,set:e=>c.thickness=e}],actions:[{label:`new sky`,hint:`Start a fresh sky right away using the current settings.`,run:()=>M()}]}],{title:`birds`,knobs:[{label:`count`,min:0,max:1e3,step:10,hint:`Number of birds. New ones arrive in small flocks.`,get:()=>v.count,set:e=>(m=e,v.resize(C,w,e))},{label:`size`,min:1.5,max:10,step:.1,hint:`How big each bird is drawn.`,get:()=>c.birdSize,set:e=>c.birdSize=e},{label:`min speed`,min:0,max:300,step:1,hint:`Slowest a bird will cruise. Can't go above max speed.`,get:()=>u.minSpeed,set:e=>u.minSpeed=Math.min(e,u.maxSpeed)},{label:`max speed`,min:5,max:500,step:1,hint:`Fastest a bird cruises on its own (gusts can push past this).`,get:()=>u.maxSpeed,set:e=>u.maxSpeed=Math.max(e,u.minSpeed)},{label:`agility`,min:5,max:400,step:5,hint:`How sharply birds can turn. Low = lazy arcs, high = twitchy.`,get:()=>u.maxForce,set:e=>u.maxForce=e},{label:`sight`,min:20,max:400,step:5,hint:`How far a bird can see its neighbours. Higher = bigger, more connected flocks.`,get:()=>u.perception,set:e=>u.perception=e},{label:`spacing`,min:2,max:100,step:1,hint:`Personal space: how close a neighbour can get before the bird moves away.`,get:()=>u.separation,set:e=>u.separation=e},{label:`align`,min:0,max:8,step:.1,hint:`How strongly a bird matches its neighbours' heading.`,get:()=>u.alignWeight,set:e=>u.alignWeight=e},{label:`cohesion`,min:0,max:8,step:.1,hint:`How strongly a bird is drawn toward the middle of its group.`,get:()=>u.cohesionWeight,set:e=>u.cohesionWeight=e},{label:`avoid`,min:0,max:15,step:.1,hint:`How strongly a bird pushes off neighbours inside its personal space.`,get:()=>u.separationWeight,set:e=>u.separationWeight=e},{label:`beat rate`,min:.5,max:12,step:.1,hint:`Wing beats per second while flapping. Hard work (turning, speeding up, gusts) beats up to 60% faster.`,get:()=>u.beatRate,set:e=>u.beatRate=e},{label:`flap effort`,min:0,max:4,step:.05,hint:`How readily turning, speeding up and gusts make a bird flap. 0 = only the cruising rhythm.`,get:()=>u.flapEffort,set:e=>u.flapEffort=e},{label:`glide`,min:0,max:10,step:.1,hint:`Average seconds a cruising bird glides between bursts of flapping. 0 = never stops flapping.`,get:()=>u.glide,set:e=>u.glide=e},{label:`beats`,min:0,max:12,step:.5,hint:`Average wing beats in a cruising burst before the next glide. 0 = only flaps when working.`,get:()=>u.beats,set:e=>u.beats=e},{label:`gust push`,min:0,max:10,step:.1,hint:`How hard moving air (your cursor) shoves birds.`,get:()=>p.gustGain,set:e=>p.gustGain=e},{label:`push cap`,min:0,max:1500,step:10,hint:`Ceiling on that shove, so a big swipe nudges rather than flings.`,get:()=>p.gustMaxAccel,set:e=>p.gustMaxAccel=e},{label:`gust speed`,min:0,max:3,step:.05,hint:`Extra speed a gust lends a bird, as a fraction of the gust's speed. Fades over ~1s.`,get:()=>p.gustBoost,set:e=>p.gustBoost=e},{label:`wake`,min:0,max:.5,step:.01,hint:`How much each bird stirs the air (and clouds) as it flies.`,get:()=>p.wake,set:e=>p.wake=e}],actions:[{label:`new flock`,hint:`Scatter a fresh flock using the current settings.`,run:()=>(v.resize(C,w,0),v.resize(C,w,m??ft(C,w)))}]}]):()=>{};return()=>{K(),cancelAnimationFrame(H),W.disconnect(),window.removeEventListener(`pointermove`,F),window.removeEventListener(`pointerdown`,F),window.removeEventListener(`pointerup`,L),document.documentElement.removeEventListener(`pointerleave`,I),window.removeEventListener(`blur`,I),O.removeEventListener(`change`,A)}}var mt=document.getElementById(`scene`),ht=new URLSearchParams(location.search),gt=ht.get(`mode`),_t=gt===`ascii`||gt===`dither`?gt:`sky`,$=ht.get(`dither`);pt(mt,{mode:_t,matrix:$===`bayer2`||$===`bayer4`||$===`void-and-cluster`?$:`bayer8`,knobs:ht.has(`knobs`)});var vt=()=>scrollTo(0,parseFloat(getComputedStyle(document.documentElement).getPropertyValue(`--tuck`))||0);vt(),addEventListener(`resize`,vt),document.addEventListener(`astro:after-swap`,vt),document.addEventListener(`astro:before-swap`,e=>{let t=document.documentElement.style.getPropertyValue(`--sky-top`);t&&e.newDocument.documentElement.style.setProperty(`--sky-top`,t)});