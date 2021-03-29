function T2 = GetT2(opts,TE,ts)

if opts.flag == 1 
    T2s = opts.val;
    
    Nm = length(opts.val);
    Ns = length(ts);
    Ne = length(TE);
    
    t = reshape(ts+TE,[Ns*Ne,1]);
    
    %% define T2* operator
    T2 = exp(-t./T2s); 
    T2 = reshape(T2,[Ns,Ne,1,1,Nm]);
else
    T2 = 1;
end

end
