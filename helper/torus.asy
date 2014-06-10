//asy -f pdf -noprc -render=4 torus.asy 

size(200);
import solids;
import graph3;
import settings;
dotfactor=6;

pen surfPen=lightgrey;
pen arcPen=black+0.7bp;

currentprojection=perspective(5,4,4);

real R=2;
real a=0.8;

triple fs(pair t) {
  return ((R+a*Cos(t.y))*Cos(t.x),(R+a*Cos(t.y))*Sin(t.x), a*Sin(t.y));
}

surface s=surface(fs,(0,0),(360,360),8,8,Spline);
draw(s,surfPen,render(compression=Low,merge=true));

int m=20;
int n=10;
real arcFactor=1.0;

pair p,q,v;

for(int i=1;i<=n;++i){
  for(int j=0;j<m;++j){
    p=(j*360/m,(i%n)*360/n);
    q=(((j+arcFactor)%m)*360/m,i*360/n);
    v=(((j+arcFactor/2)%m)*360/m,i*360/n);
    draw(fs(p)..fs(v)..fs(q),arcPen);
    q=(j*360/m,((i%n)-arcFactor)*360/n);
    draw(fs(p)..fs((p+q)/2)..fs(q),arcPen);
    dot(fs(p));
    revolution b=sphere(fs(p),0.06);
    if (j<m/2+4 && j>1)
        draw(surface(b),.9blue+opacity(.9));
    else
        draw(surface(b),.9red+opacity(.9));
  }
}
