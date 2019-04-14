function [ res ] = gradF( v )
   res = J(v)'*F(v);
end

