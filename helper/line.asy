//asy -f pdf -noprc -render=4 line.asy 

size(80);
import solids;
import graph3;
import settings;
dotfactor=4;

currentprojection=orthographic(4,-4,0);

pen arcPen=black+0.7bp;
draw((0,0,0)--(7,0,0),arcPen);

for(int i=0;i<8;++i){
    revolution b=sphere((i,0,0),0.15);
    if (i<4)
        draw(surface(b),.9red+opacity(.9));
    else
        draw(surface(b),.9blue+opacity(.9));
}
