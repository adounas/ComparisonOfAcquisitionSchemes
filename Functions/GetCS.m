function CS = GetCS(opts,fm,TE,ts)

if opts.flag == 1
    Ns = length(ts);
    Ne = length(TE);
    Nm = length(fm);
    
    t = reshape(ts+TE,[Ns*Ne,1]);
    
    %% define chemical shift operator
    CS = exp(-1i*2*pi*t*fm);
    CS = reshape(CS,[Ns,Ne,1,1,Nm]);
else
    CS = 1;
end

end