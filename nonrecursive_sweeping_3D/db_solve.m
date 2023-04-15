function u = db_solve(N,NLS,P,Tree,f)
  
  NG = N+1;
  L = round(log(NG/P)/log(2)) + 1;
  
  u = f;
  for ell=1:L
    W = 2^(ell-1)*P;
    nck = NG/W;
    %fprintf(1,'%d\n',ell);
    %tic;
    for i=1:nck
      for j=1:nck
        TT = Tree{ell}{i,j};        %LINS = TT.LINS;        LBDS = TT.LBDS;
        GINS = TT.GINS;        GBDS = TT.GBDS;
        invA = TT.invA;        B = TT.B;        S = TT.S;
        %load
        in = u(GINS);
        bd = u(GBDS);
        tt = invA*in;
        %L
        bd = bd - B*tt;
        %solve
        in = tt;
        if(ell==L)
          bd = S*bd;
        end
        %save
        u(GINS) = in;
        u(GBDS) = bd;
      end
    end
    %toc;
  end
  %-----------------
  for ell=L:-1:1
    W = 2^(ell-1)*P;
    nck = NG/W;
    %fprintf(1,'%d\n',ell);
    %tic;
    for i=1:nck
      for j=1:nck
        TT = Tree{ell}{i,j};        %LINS = TT.LINS;        LBDS = TT.LBDS;
        GINS = TT.GINS;        GBDS = TT.GBDS;
        invA = TT.invA;        B = TT.B;
        %load
        in = u(GINS);
        bd = u(GBDS);
        %Lt
        in = in - invA*(B.'*bd);
        %save
        u(GINS) = in;
        u(GBDS) = bd;
      end
    end
    %toc;
  end
  %-----------------
  
  