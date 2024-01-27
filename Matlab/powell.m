function f = powell(x)
	A = x(1) + 10.0*x(2);
	B = x(3) - x(4);
	C = x(2) - 2.0*x(3);
	D = x(1) - x(4);
    
	f = A^2 + 5.0*B^2 + C^4 + 10.0*D^4;
    return
end