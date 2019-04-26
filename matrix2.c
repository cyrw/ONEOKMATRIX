#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct{
	int r;
	int c;
	float *d;
}matrix, *pma;

struct objs;
typedef struct objs obj;
typedef struct objs *pobj;
struct objs{
    pobj parent;
    pma pm;
    pma ro;
};

struct cam{
    float fov;
    float wfov, wwfov;
    float hfov, hhfov;
    float n;
    float n2;
    pobj camo;
};
typedef struct cam cam;
typedef cam *pcam;

cam mcam;

pma initm(int r, int c){
	pma pm = malloc(sizeof(matrix));
	pm->d = malloc(sizeof(float) * r*c);
	pm->r = r;
	pm->c = c;
	return pm;
}

void freem(pma pma){
    if(pma == 0)
        return;
    if(pma->d != 0)
	    free(pma->d);
	free(pma);
}

pma multip(pma pma1, pma pma2){
	pma res;
	int r=0;
	int c=0;
	int tmp=0;
	
	if(pma1->c != pma2->r){
		printf("error diff matrix can't multip\n");
		return 0;
	}
	res = initm(pma1->r, pma2->c); // = M1 * M2

	for(r=0;r<pma1->r;r++){
		for(c = 0 ;c<pma2->c;c++){
		    //printf("\n");
			for(tmp = 0;tmp<pma1->c;tmp++){
				res->d[r*res->c+c] += pma1->d[r*pma1->c+tmp] * pma2->d[tmp*pma2->c+c];
				printf("rct %d,%d,%d  [%d][%d]%f = up + [%d][%d]=[%d]%f * [%d][%d]=[%d]%f\n",r,c,tmp   ,r,c,res->d[r*res->c+c],  r,tmp,r*pma1->c+tmp,pma1->d[r*pma1->c+tmp],   tmp,c,tmp*pma2->c+c,pma2->d[tmp*pma2->c+c]);
			}
			//printf("\n");
		}
	}
	return res;
/*
	res->d[0][0] += pma1->d[0][0] * pma2->d[0][0];
	res->d[0][0] += pma1->d[0][1] * pma2->d[1][0];

	res->d[0][1] += pma1->d[0][0] * pma2->d[1][0];
	res->d[0][1] += pma1->d[0][1] * pma2->d[1][1];

	res->d[1][0] += pma1->d[1][0] * pma2->d[0][0];
	res->d[1][0] += pma1->d[1][1] * pma2->d[1][0];

	res->d[1][1] += pma1->d[1][0] * pma2->d[1][0];
	res->d[1][1] += pma1->d[1][1] * pma2->d[1][1];
*/
}

void setdata(pma pma,float *d){
	int r = 0;
	int c = 0;
	for(;r < pma->r; r++){
		for(c = 0; c < pma->c; c++){
			pma->d[r*pma->c+c] = d[r*pma->c+c];
			//printf("%d,%d = %f ",r,c,d[r*pma->c+c]);
		}
	}
}

void setunit(pma pm){
    memset(pm->d, 0, pm->r*pm->c*sizeof(pm->d[0]));
    int ind = pm->r < pm->c ? pm->r : pm->c;
    int i = 0;
    for(;i < ind; i++){
        pm->d[i*pm->c+i] = 1;
        //printf("%d,%d = %f ",r,c,d[r*pma->c+c]);
    }
}

void showpma(pma pm){
	int r = 0;
	int c = 0;
	if(pm == 0)
	    return;
	printf("\n");
	for(;r < pm->r; r++){
		for(c = 0 ;c < pm->c; c++){
			printf("%0.2f\t",pm->d[r*pm->c+c]);
		}
		printf("\n");
	}
}

pma tmatrix(pma pm){
    pma res = initm(pm->c, pm->r);
    int r = 0;
    int c = 0;
    for(;r<pm->r; r++){
        for(c=0;c<pm->c; c++){
            //printf("%d %d = %d %d \n",c*pm->r,r,r*pm->c,c);
            res->d[c*pm->r+r] = pm->d[r*pm->c+c];
        }
    }
    return res;
}

//对称matrix
int issymmetric(pma pm){
    pma tm = tmatrix(pm);
    //showpma(tm);
    return issame(pm, tm);
}

int issame(pma pm1, pma pm2){
    if(pm1->r != pm2->r)
        return 0;
    if(pm1->c != pm2->c)
        return 0;

    //printf("sizeof(pm1->d[0]) : %d, r*c*~ = %d\n",sizeof(pm1->d[0]), pm1->r * pm1->c * sizeof(pm1->d[0]));
    //int i = 0;
    //for(;i< pm1->r * pm1->c * sizeof(pm1->d[0]); i++){
    //    printf("%x %x  ",((char *)pm1->d)[i], ((char *)pm2->d)[i]);
    //}
    return memcmp(pm1->d, pm2->d, pm1->r * pm1->c * sizeof(pm1->d[0])) == 0;
}

pma mmove(pma pm, float x, float y){
    pma dm = initm(3,2);
    float del[][2]={
    {1,0},
    {0,1},
    {x,y},
    };
    setdata(dm, del);
    printf("move : %2f %2f + %2f %2f = ",pm->d[0],pm->d[1],x,y);
    pma res = multip(pm, dm);
    printf("%2f %2f\n",res->d[0],res->d[1]);
    return res;
}

pma mscale(pma pm, float scalex,float scaley){
    pma dm = initm(3,2);
    float del[][2]={
    {scalex,0},
    {0,scaley},
    {0,0},
    };
    setdata(dm, del);
    pma res = multip(pm, dm);
    return res;
}

float rad2angle(float rad){
    return rad * (360.0f / (2 * M_PI));
}

float angle2rad(float angle){
    return angle * (2 * M_PI / 360.0f);
}


pma mrotate(pma pm, float rad){//c = sqr(x**2 + y**2)    // a = arctan(y/x)       xn = cos(a+b) * c    yn = sin(a+b)  * c
    pma dm = initm(3,2);   //xn = cos(b) * c      yn = sin(b) * c
    float del[][2]={  // xn = cos(arctan(y/x) + b) * sqr(x**2 + y**2)
    {cosf(-rad), -sinf(-rad)},
    {sinf(-rad), cosf(-rad)},
    {0,0},
    };
    setdata(dm, del);
    pma res = multip(pm, dm);
    return res;
}

pma mrotatea(pma pm, float angle){
    return mrotate(pm, angle2rad(angle));
}


pma init2dm(int x,int y){
    pma pma1 = initm(2,3);
    float d[][3]={
    {x,y,1},
    {0,0,0},
    };
    setdata(pma1, d);
    return pma1;
}

void show2dm(pma pm){

}

void pp(float x){
    int r = 40;
    int i = r + x;
    //printf("i %d\n",i);
    char *s = malloc(i+2);
    memset(s, '0', i);
    s[i]='\n';
    s[i+1]='\0';
    printf("%s",s);
    free(s);
}

void _2dtest(){
    pma p1 = init2dm(3,4);
    showpma(p1);
    pma res = mmove(p1,6,7);
    showpma(res);
    show2dm(res); // = 4,5
    freem(res);

    printf("scale 2,2 \n");
    res = mscale(p1, 0.5, 0.5);
    showpma(res);
    freem(res);

    res = mrotatea(p1,180);
    showpma(res);
    freem(res);


    pma p2 = init2dm(40,0);

    struct timespec req = {0}, rem = {0};
    req.tv_nsec = 20 * 1000000L;
    req.tv_sec = 0;
    int i = 0;
    for(;i < 360*8+1; i++){
        res = mrotatea(p2,i);
        //showpma(res);
        pp(res->d[0]);
        //printf("%2f\n",res->d[0]);
        freem(res);
        nanosleep(&req, &rem);
    }
}

pma init3dm(float x, float y, float z){
    pma res = initm(4,1);
    float d[][1]={
    {x},
    {y},
    {z},
    {1}
    };
    setdata(res, d);
    return res;
}

pma move3d(pma in, float dx, float dy, float dz){
    pma pre = initm(4,4);
    float d[][4] = {
    {1,0,0,dx},
    {0,1,0,dy},
    {0,0,1,dz},
    {0,0,0,1},
    };
    setdata(pre, d);
    pma res = multip(pre ,in);
    return res;
}

pma scale3d(pma in, float sx,float sy,float sz){
    pma pre = initm(4,4);
    float d[][4] = {
    {sx,0,0,0},
    {0,sy,0,0},
    {0,0,sz,0},
    {0,0,0,1},
    };
    setdata(pre, d);
    pma res = multip(pre ,in);
    return res;
}


float Q_rsqrt( float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                       // evil floating point bit level hacking（对浮点数的邪恶位元hack）
	i  = 0x5f3759df - ( i >> 1 );               // what the fuck?（这他妈的是怎么回事？）
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration （第一次迭代）
//      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed（第二次迭代，可以删除）

	return y;
}

float len3d(pma in){
    return sqrtf(in->d[0]*in->d[0] + in->d[1]*in->d[1] + in->d[2]*in->d[2]);
    //return Q_rsqrt(pma->d[0]*pma->d[0] + pma->d[1]*pma->d[1] + pma->d[2]*pma->d[2]);
}



pma normaliz3d(pma in){
    pma pre = initm(4,4);
    float len = len3d(in);
    float dlen = 1/len;
    printf("len %2f   dlen %2f \n",len, dlen);
    //dlen = Q_rsqrt(in->d[0]*in->d[0] + in->d[1]*in->d[1] + in->d[2]*in->d[2]);

    float d[][4] = {
    {dlen,0,0,0},
    {0,dlen,0,0},
    {0,0,dlen,0},
    {0,0,0,1},
    };
    setdata(pre, d);
    pma res = multip(pre ,in);
    return res;
}

//left hand cordi
pma rotate3dr(pma in, float x ,float y, float z){
        //z
        pma dm = initm(4,4);   //xn = cos(b) * c      yn = sin(b) * c

        float del[][4]={  
        {cosf(z), -sinf(z),0,0},
        {sinf(z), cosf(z),0,0},
        {0,0,1,0},
        {0,0,0,1},
        };
        setdata(dm, del);
        showpma(dm);
        printf("* ");
        showpma(in);
        pma res = multip(dm, in);
        printf("= ");
        showpma(res);

        //y
        //dm = initm(4,4);
        pma in2 = res;
        float dely[][4]={  
        {cosf(y),0,sinf(y),0},
        {0,1,0,0},
        {-sinf(y),0,cosf(y),0},
        {0,0,0,1},
        };
        setdata(dm, dely);
        showpma(dm);
        printf("* ");
        showpma(in2);
        res = multip(dm, in2);
        printf("= ");
        showpma(res);
        //freem(dm);
        //freem(in2);

        //x
        pma in3 = res;
        float delx[][4]={
        {1,0,0,0},
        {0,cosf(-x), sinf(-x),0},
        {0,-sinf(-x), cosf(-x),0},
        {0,0,0,1},
        };
        setdata(dm, delx);
        showpma(dm);
        printf("* ");
        showpma(in2);
        res = multip(dm, in3);
        printf("= ");
        showpma(res);

        freem(dm);

        return res;
}

//angle
pma rotate3d(pma in, float x ,float y, float z){
        return rotate3dr(in,angle2rad(x), angle2rad(y), angle2rad(z));
}

//pma rero3d(pma in, float x float y, float z){
   
//}

// ----> z
pma came3d(float n, pma in){
    //point
    //rotate  z->
    //fov = (0,180)
    //n, deep?
    float sf = (n / in->d[2]);
    float xn = in->d[0] * sf;
    float yn = in->d[1] * sf;
    printf("o: %2f,%2f,%2f n: %2f,%2f,%2f \n", in->d[0],in->d[1],in->d[2], xn,yn,n);
    pma res = init3dm(xn,yn,in->d[2]);
    return res;
}

pma came3d2(pma in){
    //n n2
    //x,y compose ,enqueue
    //dequeue set pixel color ,draw
    pma res;
    int i = 0;
    int w = 80;
    int h = 24;
    
     if(mcam.n > in->d[2] || mcam.n2 < in->d[2])
        return;
     
     float mtanx = tanf(angle2rad(mcam.wwfov));
     float mtany = tanf(angle2rad(mcam.hhfov));
     float ptanx = in->d[0] / in->d[2];
     float ptany = in->d[1] / in->d[2];
     if(ptanx > mtanx || ptanx < -mtanx)
        return;
     if(ptany > mtany || ptany < -mtany)
        return;
        
     //x,y compose
     float npx = ptanx / mtanx;
     float npy = in->d[2] * mtany * in->d[1];
     
     //mapping to n
     float sf = mcam.n / in->d[2];
     float npx2 = npx * sf;
     float npy2 = npy * sf;

    res = init3dm(npx2, npy2, in->d[2]);
    return res;
}

pma came3d3(pma in){
    //1. point's z is    n < z < n2 ?
    //2. x,y,z compose to 0~1f ,0~1f ,0~1f 
    //3. return
    //4. (out)dequeue set pixel color ,draw
    pma res;
    int i = 0;
    if(in == 0)
        return 0;
        
    float distn = mcam.n + mcam.camo->pm->d[2];
    float distn2 = mcam.n2 + mcam.camo->pm->d[2];
    if(distn > in->d[2] || distn2 < in->d[2])
       return 0;
     
     float nd = distn2 - distn;
     float mtanx = tanf(angle2rad(mcam.wwfov));
     float mtany = tanf(angle2rad(mcam.hhfov));
     float p2cz = in->d[2] - mcam.camo->pm->d[2];
     float p2cx = in->d[0] - mcam.camo->pm->d[0];
     float p2cy = in->d[1] - mcam.camo->pm->d[1];
     float ptanx = p2cx / p2cz; //相对于cam作为0，0，0 算tan
     float ptany = p2cy / p2cz;
     if(ptanx > mtanx || ptanx < -mtanx)
        return 0;
     if(ptany > mtany || ptany < -mtany)
        return 0;
    
    float tmpz = p2cz * 2;
    float camd[][4]={
    {1/(tmpz*mtanx),0,0,0.5},
    {0,-1/(tmpz*mtany),0,0.5},
    {0,0,(1-distn/p2cz)/nd,0},
    {0,0,0,1},
    };
    /*
    printf("mcam.wwfov %f\n",mcam.wwfov);
    printf("mtanx %f\n",mtanx);
    printf("tmpz %f\n",tmpz);
    printf("tmpz*mtanx %f\n",tmpz*mtanx);
    printf("1/(tmpz*mtanx) %f\n",1/(tmpz*mtanx));
    */
    pma camp = initm(4,4);
    setdata(camp, camd);
    res = multip(camp, in);
    freem(camp);
    
    return res;
}

void show3ds(pma *in, int len){
    //w 80 h 24 // c : 40,12
    int r = 0;
    int c = 0;
    int w = 80;
    int h = 24;

    //int x = (int)in->d[0];
    //int y = (int)in->d[1];
    int i = 0;
    int hasp = 0;
    //printf("%d,%d\n",x,y);
    
    printf("len %d ",len);
    printf("(int)in[i]->d[0] %d",(int)in[i]->d[0]);
    
    char b = ' ';
    char p = '@';
    char p2 ='0';
    char p3 ='*';
    char p4 ='.';
    char z = '+';
    for(r=0;r<h;r++){
        for(c=0;c<w;c++){
            hasp = 0;
            //printf("%d,%d ",c-(80/2),r-(24/2));
            for(i=0;i<len;i++){
                if((int)in[i]->d[0] == c-(w/2) && (int)in[i]->d[1] == (h/2)-r){
                    if((int)in[i]->d[2] > 0 && (int)in[i]->d[2] < 4)
                        write(1,&p,1);
                    else if((int)in[i]->d[2] < 8){
                        write(1,&p2,1);
                    }else if((int)in[i]->d[2] < 12){
                        write(1,&p3,1);
                    }else{
                        write(1,&p4,1);
                    }
                    hasp = 1;
                    //printf("%d ",i);
                    break;
                }
            }
            if(hasp == 0){
                if(c == w/2 && r == h/2)
                    write(1,&z,1);
                else
                    write(1,&b,1);
            }
        }
        printf("\n");
    }
}

void show3ds2(pma *in, int len){
    //w 80 h 24 // c : 40,12
    int r = 0;
    int c = 0;
    int w = 80;
    int h = 24;

    //int x = (int)in->d[0];
    //int y = (int)in->d[1];
    int i = 0;
    int hasp = 0;
    //printf("%d,%d\n",x,y);
    
    printf("len %d ",len);
    if(in != 0 && in[0] != 0 && in[0]->d != 0)
        printf("x,y,z %2f,%2f,%2f\n",in[i]->d[0],in[i]->d[1],in[i]->d[2]);
    
    char b = ' ';
    char p = '@';
    char p2 ='0';
    char p3 ='*';
    char p4 ='.';
    char z = '+';
    for(r=0;r<h;r++){
        for(c=0;c<w;c++){
            hasp = 0;
            //printf("%d,%d ",c-(80/2),r-(24/2));
            for(i=0;i<len;i++){
                if(in != 0 && in[i] != 0 && in[i]->d != 0 && (int)in[i]->d[0] == c && (int)in[i]->d[1] == r){
                    if(in[i]->d[2] > 0 && in[i]->d[2] < 1.0f/4)
                        write(1,&p,1);
                    else if(in[i]->d[2] < 1.0f/4*2){
                        write(1,&p2,1);
                    }else if(in[i]->d[2] < 1.0f/4*3){
                        write(1,&p3,1);
                    }else{
                        write(1,&p4,1);
                    }
                    hasp = 1;
                    //printf("%d ",i);
                    break;
                }
            }
            if(hasp == 0){
                if(c == w/2 && r == h/2)
                    write(1,&z,1);
                else
                    write(1,&b,1);
            }
        }
        printf("\n");
    }
}

pobj initobj(pobj par){
    pobj po = malloc(sizeof(obj));
    po->parent = par;
    return po;
}

pma getwordp(pobj obj){
    pma pm = init3dm(0,0,0);
    do{
        printf("===\n");
        showpma(pm);
        showpma(obj->pm);
        pma res = move3d(pm, obj->pm->d[0],obj->pm->d[1],obj->pm->d[2]);
        freem(pm);
        showpma(res);
        pm = res;
        obj = obj->parent;
    }while(obj != 0);
    return pm;
}

pma fitpix(pma in,int w,int h){
    if(in == 0)
        return 0;
        
    pma fitm = initm(3,4);
    float fitd[][4]={
    {w,0,0,0},
    {0,h,0,0},
    {0,0,1,0}
    };
    setdata(fitm, fitd);
    return multip(fitm, in);
}

void _3dtest(){
//world
    pma p1;
    pobj world = initobj(0);
    world->pm = init3dm(0,0,0);

    pobj o1 = initobj(world);
    o1->pm = init3dm(0,8,22);
    
    pobj obj1 = initobj(o1);
    obj1->pm = init3dm(14,0,0);

    pobj obj2 = initobj(obj1);
    obj2->pm = init3dm(4,0,0);
    
//cam
    mcam.fov = 90;
    mcam.wfov = mcam.fov;
    mcam.hfov = mcam.fov;//mcam.wfov / 80 * 24;//
    mcam.wwfov = mcam.wfov / 2;
    mcam.hhfov = mcam.hfov / 2;
    mcam.n = 2;
    mcam.n2 = 42;
    
    mcam.camo = initobj(world);
    mcam.camo->pm = init3dm(7,0,0);
    
    
    //camo->ro = ;

//frame rate
struct timespec ts;
ts.tv_sec = 0;
ts.tv_nsec = 100000000L * 1;


pma spma1[3];
pma spma2[3];
pma spma3[3];
do{
    p1 = rotate3d(obj1->pm, 0,5,0);
    freem(obj1->pm);
    //showpma(p1);
    obj1->pm = p1;

    //p1 = rotate3d(obj2->pm, 0,0,30);
    //freem(obj2->pm);
    //obj2->pm = p1;

    printf("world point ======\n");
    spma1[0] = getwordp(obj1);
    //spma1[1] = getwordp(obj2);
    //printf("world point spma1 ======\n");
    showpma(spma1[0]);
    //TODO CAM +XYZ 
    //FOV splash
        
    spma2[0] = came3d3(spma1[0]);
    //spma2[1] = came3d3(spma1[1]);
    //printf("normalize point =======\n");
    //showpma(spma2[0]);
    //spma2[1] = came3d2(spma1[1]);
    freem(spma1[0]);
    //freem(spma1[1]);
    
    //printf("fit pix =======\n");
    spma3[0] = fitpix(spma2[0], 80, 24);
    //spma3[1] = fitpix(spma2[1], 80, 24);
    //showpma(spma3[0]);
    show3ds2(&spma3,1);
    freem(spma2[0]);
    //freem(spma2[1]);
    freem(spma3[0]);
    //freem(spma3[1]);

    //freem(spma2[1]);
    nanosleep(&ts, 0);
}while(1);
return;



    showpma(p1);
    pma p2 = move3d(p1,0,-2,-3);
    showpma(p2);
    //原点对称
    pma p3 = scale3d(p2, -1,-1,-1);
    showpma(p3);

    pma p4 = normaliz3d(p3);
    showpma(p4);

    printf("====== ro \n");
    pma p5 = init3dm(2,2,2);
    showpma(p5);
    pma p51 = rotate3d(p5, 0,90,90); //90 2,-2,-2  //2,0,0 > 0,0,-2
    showpma(p51);
    freem(p51);

    p51 = rotate3d(p5, 180,180,0); //90 2,-2,-2  //2,0,0 > 0,0,-2
    showpma(p51);
    freem(p51);

    p51 = rotate3d(p5, 270,270,0); //90 2,-2,-2  //2,0,0 > 0,0,-2
    showpma(p51);
    freem(p51);

    pma poa[3];
    pma poao[3];
    pma poat[3];
    pma poar[3];

    poa[0] = init3dm(10,3,9);
    poa[1] = init3dm(10,3,2);
    poa[2] = init3dm(10,3,5);
    poat[0] = poa[0];
    poat[1] = poa[1];
    poat[2] = poa[2];


    
    float rx = 0;
    float ry = 0;
    float rz = 0;
    //cam in 000 ,rotate 000
    char in = '1';
    while(in != 'q'){
        printf("c %c \n",in);

        poat[0] = rotate3d(poa[0], rx, ry, rz);
        poat[1] = rotate3d(poa[1], rx, ry, rz);
        poat[2] = rotate3d(poa[2], rx, ry, rz);
        poar[0] = came3d(2,poat[0]);
        poar[1] = came3d(2,poat[1]);
        poar[2] = came3d(2,poat[2]);
        freem(poa[0]);
        freem(poa[1]);
        freem(poa[2]);
        poa[0]=poat[0];
        poa[1]=poat[1];
        poa[2]=poat[2];
        //freem(poat[0]);
        //freem(poat[1]);
        //freem(poat[2]);
        //showpma(po2);
        show3ds(poar,3);
        freem(poar[0]);
        freem(poar[1]);
        freem(poar[2]);
        ry=0,rx=0,rz=0;
        in = getchar();
        //rz+=10;
        switch (in){
            case 'a':
                ry+=45;
               // poa[0]->d[0]--;
               // poa[1]->d[0]--;
               // poa[2]->d[0]--;
                break;
            case 's':
                rx+=45;
                //poa[0]->d[1]++;
              //  poa[1]->d[1]++;
              //  poa[2]->d[1]++;
                break;
            case 'd':
                ry-=45;
             //   poa[0]->d[0]++;
            //    poa[1]->d[0]++;
              //  poa[2]->d[0]++;
                break;
            case 'w':
                rx-=45;
             //   poa[0]->d[1]--;
             //   poa[1]->d[1]--;
              //  poa[2]->d[1]--;
                break;
        }
    }
}

void main(){
    int i = 0;
    setvbuf(stdout,0,_IONBF,0);
/*
	pma pma1 = initm(2,2);
	pma pma2 = initm(2,2);
	float d[][2]={{1,2},{3,4}};
	setdata(pma1, d);
	float d2[][2]={{2,3},{4,5}};
	setdata(pma2, d2);
	showpma(pma1);
	showpma(pma2);
	pma res = multip(pma1,pma2);
	showpma(res);
	freem(res);

	res = initm(3,7);
	setunit(res);
	showpma(res);
	pma res2 = tmatrix(res);
	showpma(res2);
	freem(res2);
	freem(res);

	res = initm(3,3);
	setunit(res);
	showpma(res);
    i = issymmetric(res);
    printf("issymmetric %d\n",i);
    */

    //_2dtest();
    _3dtest();
}


