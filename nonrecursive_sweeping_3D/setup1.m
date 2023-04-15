function[P]=setup1(NPML,NLPD,NPAD,pL,pR,h,ksq,s1,s2,s3,leafsize)
[N1,N2,N3]=size(ksq);

%generate the partitions
M=(N3+1)/2;
ordL=NLPD+NPML-1:NLPD:M-1;
NP=length(ordL);
pttnL=cell(1,NP);
pttnL{1}=[1,NLPD+NPML-1];
for g=2:NP
    pttnL{g}=[ordL(g-1)+1,ordL(g)];
end
ordR=N3+1-ordL(end:-1:1);
pttnR=cell(1,NP);
for g=1:NP-1
    pttnR{g}=[ordR(g),ordR(g+1)-1];
end
pttnR{NP}=[N3+1-(NLPD+NPML-1),N3];
pttnM=[ordL(end)+1,ordR(1)-1];

    %construct left sweeping matrices
    PL=cell(1,NP);
    if(1)
        b=1;
        fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
        ksqnew=ksq(:,:,fm:to);
        s3new=s3(2*fm-1:2*to+1);
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PL{b}=db_setup(N1,cnt,leafsize,tmp);
        fprintf('complete 1: PL fm=%d to=%d\n',fm,to);
    end
    for b=2:NP
        fm=pttnL{b}(1);to=pttnL{b}(2);cnt=to-fm+1;
        ksqnew=ksq(:,:,fm-(NPAD-1):to);
        s3new=[pL,s3(2*fm:2*to+1)];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PL{b}=db_setup(N1,cnt+NPAD-1,leafsize,tmp);
        fprintf('complete 1: PL fm=%d to=%d\n',fm,to);
    end
    %construct right sweeping matrices
    PR=cell(1,NP);
    for b=1:NP-1
        fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
        ksqnew=ksq(:,:,fm:to+(NPAD-1));
        s3new=[s3(2*fm-1:2*to),pR];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PR{b}=db_setup(N1,cnt+NPAD-1,leafsize,tmp);
        fprintf('complete 1: PR fm=%d to=%d\n',fm,to);
    end
    if(1)
        b=NP;
        fm=pttnR{b}(1);to=pttnR{b}(2);cnt=to-fm+1;
        ksqnew=ksq(:,:,fm:to);
        s3new=s3(2*fm-1:2*to+1);
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PR{b}=db_setup(N1,cnt,leafsize,tmp);
        fprintf('complete 1: PR fm=%d to=%d\n',fm,to);
    end
    %construct middle layers
    if(1)
        fm=pttnM(1);to=pttnM(2);cnt=to-fm+1;
        ksqnew=ksq(:,:,fm-(NPAD-1):to+(NPAD-1));
        s3new=[pL,s3(2*fm:2*to),pR];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PM=db_setup(N1,cnt+2*(NPAD-1),leafsize,tmp);
        fprintf('complete 1: PM fm=%d to=%d\n',fm,to);
    end

    shinv=reshape(-1/(h*h)*kron(1./s2(2:2:2*N2),1./s1(2:2:2*N1)),[N1,N2,1]);
    P={{PL,pttnL},{PR,pttnR},{PM,pttnM},size(ksq),shinv};
end