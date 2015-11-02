function p = modulo(m,n)

p = mod(m,n);
if any(p == 0)
    p(p == 0) = n;
end