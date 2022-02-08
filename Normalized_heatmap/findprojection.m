function p=findprojection(p1,p2,q)
%% p1 p2 is two point one the line, q is the point to be projected, p is the projection point

syms x y
eqn1= (p1(2)-p2(2))*(p1(1)-x)==(p1(1)-p2(1))*(p1(2)-y);
eqn2= (q(1)-x)*(p1(1)-p2(1))+(q(2)-y)*(p1(2)-p2(2))==0;
sol=solve([eqn1, eqn2],[x,y],'Real',true);
p(1)=double(sol.x);
p(2)=double(sol.y);

end
