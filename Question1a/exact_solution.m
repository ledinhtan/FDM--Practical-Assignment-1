function u_ex=exact_solution(x,k)
if (k==1)
    u_ex=2*x^4-x^3+1;
end

if(k==2)
    u_ex=-x^4+x^3+x+1;
end