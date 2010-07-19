#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;


enum {DIM=2};

double square(double x){  return x*x;}
typedef struct P P;
struct P {
  double x[DIM];
  P(const P&p){
    for(int i=0;i<DIM;i++)
      x[i]=p.x[i];
  }
  P&operator=(const P&p){
    for(int i=0;i<DIM;i++)
      x[i]=p.x[i];
    return *this;
  }
  bool operator==(const P&p)const{
    for(int i=0;i<DIM;i++)
      if(x[i]!=p.x[i])
	return false;
    return true;
  }
  double&operator[](int i){
    return x[i];
  }
  P(double a=.0,double b=.0,double c=.0){
    x[0]=a;
    if(DIM>1) x[1]=b;
    if(DIM>2) x[2]=c;
  }
  P operator-(const P&y){
    P q;
    for(int i=0;i<DIM;i++)
      q.x[i]=x[i]-y.x[i];
    return q;
  }
  P operator+(const P&y){
    P q;
    for(int i=0;i<DIM;i++)
      q.x[i]=x[i]+y.x[i];
    return q;
  }
  double operator*(const P&y){
    double sum=0.;
    for(int i=0;i<DIM;i++)
      sum+=x[i]*y.x[i];
    return sum;
  }
  double norm(){
    double sum=0.;
    for(int i=0;i<DIM;i++)
      sum+=x[i]*x[i];
    return sqrt(sum);
  }
};
  
double dist(P&p,const P&q){
  return (p-q).norm();
}

/*
int
main()
{
  P p(1.,1.),q(2.,1.),d=p-q;
  cout << d.norm() << endl;
}
*/


typedef struct B B;
struct B{
  P min,max;
  B(){}
  B(const P&min_,const P&max_):min(min_),max(max_){}
};

double dist(B&b,P&p)
{
  double sum=0;
  for(int i=0;i<DIM;i++){
    if(p[i]<b.min[i])sum+=square(p[i]-b.min[i]);
    if(b.max[i]<p[i])sum+=square(p[i]-b.max[i]);
  }
  return sum;
}
typedef struct N N;
struct N:B{
  int mom,left,right,ptmin,ptmax;
  N(){}
  N(P pi,P pa,int m,int l,int r,int mi,int ma):
    B(pi,pa),mom(m),left(l),right(r),ptmin(mi),ptmax(ma){}
};

struct T{
  static const double BIG;
  vector<P > &points;
  int nb,n;
  N*boxes;
  vector<int>ptindex,rptindex;
  double*coords;
  T(vector<P>&points_);
  int locate(P p);
  int nearest(P p);
};

const double T::BIG(1e20);


int selecti(const int k,int*in,int n,double*ar)
{
#define swap(a,b) do{int h=in[a];in[a]=in[b];in[b]=h;}while(0)
#define r(l) ar[in[l]]
  int l=0,ir=n-1;
  for(;;){
    if(ir<=l+1){
      if(ir==l+1 && r(ir)<r(l))
	swap(l,ir);
      return in[k];
    }else{
      int mid=(1+ir)/2;
      swap(mid,l+1);
      if(r(l)>r(ir)) swap(l,ir);
      if(r(l+1)>r(ir)) swap(l+1,ir);
      if(r(l)>r(l+1)) swap(l,l+1);
      int i=l+1,j=ir,ia=in[i];
      double a=ar[ia];
      for(;;){
	do i++; while(r(i)<a);
	do j--; while(r(j)>a);
	if(j<i) break;
	swap(i,j);
      }
      in[l+1]=in[j];
      in[j]=ia;
      if(j>=k) ir=j-1;
      if(j<=k) l=i;
    }   
  }
#undef swap
#undef r
}
/*
int main(){ 
  const int n=7;
  int in[n]={0,1,2,3,4,5,6};
  double ar[n]={.0,.8,.1,.3,.3,.9,.6};
  cout << "selecti=" << selecti(n/2,in,n,ar) <<endl;
  for(int i=0;i<n;i++)
    cout << in[i] << " " << ar[in[i]] << endl;
  return 0;
}
*/

int next_axis(int axis){
  return (1+axis)%DIM;
}

T::T(vector<P>&points_):points(points_),n(points_.size()),ptindex(n),rptindex(n){
  for(int k=0;k<n;k++)
    ptindex[k]=k;
  int m=1<<(int)(ceil(log2(n)));
  nb=min(m,2*n-m/2)-1;
  boxes=new N[nb];
  coords=new double[DIM*n];
  for(int j=0;j<DIM;j++)
    for(int i=0;i<n;i++)
      coords[i+DIM*j]=points[i][j];
  P min(-BIG,-BIG,-BIG),max(BIG,BIG,BIG);
  boxes[0]=N(min,max,0,0,0,0,n-1);
  int jbox=0,taskmom[50],taskaxis[50],nowtask=1;
  taskaxis[nowtask]=0;
  taskmom[nowtask]=0;
  while(nowtask){
    int tmom=taskmom[nowtask],
      axis=taskaxis[nowtask--],
      ptmin=boxes[tmom].ptmin,
      ptmax=boxes[tmom].ptmax,
      *hp=&ptindex[ptmin],
      np=ptmax-ptmin+1,
      kk=(np-1)/2;
    double *cp=&coords[axis*n];
    selecti(kk,hp,np,cp);
    max=boxes[tmom].max;
    min=boxes[tmom].min;
    max[axis]=min[axis]=coords[axis*n+hp[kk]];
    boxes[jbox++]=N(boxes[tmom].min,max,tmom,0,0,ptmin,ptmin+kk);
    boxes[jbox++]=N(min,boxes[tmom].max,tmom,0,0,ptmin+kk+1,ptmax);
    boxes[tmom].left=jbox-1;
    boxes[tmom].right=jbox;
    if(kk>1){
      nowtask++;
      taskmom[nowtask]=jbox-1;
      taskaxis[nowtask]=next_axis(axis);
    }
    if(np-kk>2){
      nowtask++;
      taskmom[nowtask]=jbox;
      taskaxis[nowtask]=next_axis(axis);
    }
  }
  for(int j=0;j<n;j++)
    rptindex[ptindex[j]]=j;
  delete [] coords;
}

int T::locate(P p)
{
  int b,axis;
  b=axis=0;
  while(boxes[b].left){
    int l=boxes[b].left;
    if(p[axis]<=boxes[l].max[axis])
      b=l;
    else
      b=boxes[b].right;
    axis=next_axis(axis);
  }
  return b;
}

int T::nearest(P p)
{
  int k=locate(p),nearest=0;
  double nearest_dist=(points[nearest]-p).norm();
  for(int i=boxes[k].ptmin+1;i<=boxes[k].ptmax;i++){
    double d=(points[ptindex[i]]-p).norm();
    if(d<nearest_dist){
      nearest_dist=d;
      nearest=ptindex[i];
    }
  }
  int ntask=1,task[50];
  task[ntask]=0;
  while(ntask){
    k=task[ntask--];
    if(dist(boxes[k],p)<nearest_dist){
      if(boxes[k].left){
	task[++ntask]=boxes[k].left;
	task[++ntask]=boxes[k].right;
      } else {
	for(int i=boxes[k].ptmin;i<=boxes[k].ptmax;i++){
	  double d=(points[ptindex[i]]-p).norm();
	  if(d<nearest_dist){
	    nearest_dist=d;
	    nearest=ptindex[i];
	  }
	}
      }
    }
  }
    
  return nearest;
}

int main()
{
  int n=10000000;
  vector<P> points(n);
  for(int i=0;i<n;i++)
    points[i]=P(drand48(),drand48());
  cout << points.size() << endl;
  T t(points);
  
  cout << t.nearest(P(.5,.5)) << endl;
  return 0;
}
