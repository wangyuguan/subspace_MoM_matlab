function [real_t_complex, complex_t_real] = get_view_real_t_complex(P)

vec_idx = @(p, u) 2*p^2+p+u+1;
n = 2*P^2+3*P+1;
real_t_complex = zeros([n,n]);

for p=0:P
    for u=-2*p:2*p

        if u==0 
            real_t_complex(vec_idx(p,u), vec_idx(p,u)) = 1;
        elseif u>0 
            real_t_complex(vec_idx(p,u), vec_idx(p,u)) = 1;
            real_t_complex(vec_idx(p,u), vec_idx(p,-u)) = (-1)^(u)*1j;
        else
            real_t_complex(vec_idx(p,u), vec_idx(p,u)) = -1j;
            real_t_complex(vec_idx(p,u), vec_idx(p,-u)) = (-1)^(u);
        end

    end

end

complex_t_real = inv(real_t_complex);


end