function[res]=apply1(NPAD,P,leafsize,inc)
PL=P{1}{1};pttnL=P{1}{2};
PR=P{2}{1};pttnR=P{2}{2};
PM=P{3}{1};pttnM=P{3}{2};
sz=P{4};shinv=P{5};
N1=sz(1);N2=sz(2);N3=sz(3);
inc=reshape(inc,sz);
NP=numel(pttnL);
res=inc;
%left solve
if(1)
    b=1;
    fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
    res(:,:,fm:to)=db_solve(N1,cnt,leafsize,PL{b},res(:,:,fm:to));
    transL=res(:,:,to);
end
for b=2:NP
    fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
    res(:,:,fm)=res(:,:,fm)+shinv.*transL;
    tmp=zeros([N1,N2,cnt+NPAD-1]);
    tmp(:,:,NPAD:end)=res(:,:,fm:to);
    tmp=db_solve(N1,cnt+NPAD-1,leafsize,PL{b},tmp);
    res(:,:,fm:to)=tmp(:,:,NPAD:end);
    transL=res(:,:,to);
end
%right solve
if(1)
    b=NP;
    fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
    res(:,:,fm:to)=db_solve(N1,cnt,leafsize,PR{b},res(:,:,fm:to));
    transR=res(:,:,fm);
end
for b=NP-1:-1:1
    fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
    res(:,:,to)=res(:,:,to)+shinv.*transR;
    tmp=zeros([N1,N2,cnt+NPAD-1]);
    tmp(:,:,1:cnt)=res(:,:,fm:to);
    tmp=db_solve(N1,cnt+NPAD-1,leafsize,PR{b},tmp);
    res(:,:,fm:to)=tmp(:,:,1:cnt);
    transR=res(:,:,fm);
end
%middle solve
if(1)
    fm=pttnM(1);to=pttnM(2);cnt=to-fm+1;
    res(:,:,fm)=res(:,:,fm)+shinv.*transL;
    res(:,:,to)=res(:,:,to)+shinv.*transR;
    tmp=zeros([N1,N2,cnt+2*(NPAD-1)]);
    tmp(:,:,NPAD:cnt+(NPAD-1))=res(:,:,fm:to);
    tmp=db_solve(N1,cnt+2*(NPAD-1),leafsize,PM,tmp);
    res(:,:,fm:to)=tmp(:,:,NPAD:cnt+(NPAD-1));
    transL=res(:,:,fm);
    transR=res(:,:,to);
end
%left back solve
for b=NP:-1:2
    fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
    tmp=zeros([N1,N2,cnt+NPAD-1]);
    tmp(:,:,end)=shinv.*transL;
    tmp=db_solve(N1,cnt+NPAD-1,leafsize,PL{b},tmp);
    res(:,:,fm:to)=res(:,:,fm:to)+tmp(:,:,NPAD:end);
    transL=res(:,:,fm);
end
if(1)
    b=1;
    fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
    tmp=zeros([N1,N2,cnt]);
    tmp(:,:,end)=shinv.*transL;
    tmp=db_solve(N1,cnt,leafsize,PL{b},tmp);
    res(:,:,fm:to)=res(:,:,fm:to)+tmp;
end
%right back solve
for b=1:NP-1
    fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
    tmp=zeros([N1,N2,cnt+NPAD-1]);
    tmp(:,:,1)=shinv.*transR;
    tmp=db_solve(N1,cnt+NPAD-1,leafsize,PR{b},tmp);
    res(:,:,fm:to)=res(:,:,fm:to)+tmp(:,:,1:cnt);
    transR=res(:,:,to);
end
if(1)
    b=NP;
    fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
    tmp=zeros([N1,N2,cnt]);
    tmp(:,:,1)=shinv.*transR;
    tmp=db_solve(N1,cnt,leafsize,PR{b},tmp);
    res(:,:,fm:to)=res(:,:,fm:to)+tmp;
end

    res=res(:);
end

