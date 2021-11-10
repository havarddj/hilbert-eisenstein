/*

October 31, 2018.

Code to compute signs of a ray class character

*/

D:=12;
p:=23;
K:=QuadraticField(D);
O:=Integers(K);

pp:=ideal<O|p>;
Fs:=Factorisation(pp);
P:=Fs[1][1];


// Computes the kernels of the quotient maps
// from the narrow class group to the
// two quite-narrow class groups.

// I believe the signs of a character on
// the narrow group are given by its orders
// when restricted to the two kernels:
// order 1 = +1; order 2 = -1.
// DO YOU AGREE JAN? 

function from_narrow_kernels(P)

	G12,m12:=RayClassGroup(P,[1,2]); // narrow group
	G1,m1:=RayClassGroup(P,[1]); // quite-narrow group
	G2,m2:=RayClassGroup(P,[2]); // quite-narrow group

	m1inv:=Inverse(m1);
	m2inv:=Inverse(m2);
	m12inv:=Inverse(m12);

	// Map from narrow to narrow1 and narrow2

	ker1:=[];
	for g in G12 do
		if m1inv(m12(g)) eq G1!0 then
    		Append(~ker1,g);
    	end if;
	end for;

	ker2:=[];
	for g in G12 do
		if m2inv(m12(g)) eq G2!0 then
    		Append(~ker2,g);
    	end if;
	end for;

	return ker1,ker2,G12,m12,G1,m1,G2,m2;

end function;

//





